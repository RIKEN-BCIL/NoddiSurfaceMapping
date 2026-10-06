#!/bin/bash
# matlab_noddi_bin.sh <bin dir> <out dir> <rician 0|1> [rician dir]   2026-09-19
# env: NODDIDIR, NIFTIMATLIB, MATLABBIN, DPAR (parallel diffusivity, m^2/s; 1.1e-9 default = cortex, 1.7e-9 for white matter)
# One bin (or chunk) of a NODDI fit with MATLAB: <bin dir> holds data.nii.gz,
# nodif_brain_mask.nii.gz, bvals, bvecs (the bin's own b-table when it is a
# gradient-nonlinearity bin). Writes cuDIMOT-named maps into <out dir>
# (mean_fintra, mean_fiso, mean_kappa, OD, dyads1, fmin) so that the merge and
# the mapping stages are shared with the cuDIMOT fitters, then DONE (or FAIL).
#   rician 0 : the NODDI toolbox (NODDIDIR), WatsonSHStickTortIsoV_B0, batch_fitting_single
#   rician 1 : the Rician NODDI at <rician dir>, which must provide
#              noddi_rician_bin(datafile, maskfile, bvalsfile, bvecsfile, outprefix)
#              writing <outprefix>_ficvf.nii, _fiso.nii, _kappa.nii, _odi.nii,
#              _fibredirs_{x,y,z}vec.nii and _fmin.nii like SaveParamsAsNIfTI does
set -u
B="$1"; O="$2"; RIC="${3:-0}"; RDIR="${4:-}"
NODDIDIR=${NODDIDIR:?export NODDIDIR}; NIFTIMATLIB=${NIFTIMATLIB:?export NIFTIMATLIB}; MATLABBIN=${MATLABBIN:-matlab}; DPAR=${DPAR:-1.1e-9}
export FSLDIR=${FSLDIR:-/usr/local/fsl}; export PATH=$FSLDIR/bin:$PATH
mkdir -p $O; rm -f $O/DONE $O/FAIL
for f in data nodif_brain_mask; do [ -e $O/$f.nii ] || gunzip -c $B/$f.nii.gz > $O/$f.nii; done   # the toolbox reads plain .nii
# the toolbox's FSL2Protocol takes every distinct b as a shell and b<=threshold as b0: HCP-style
# tables (b0 = 5, b = 990/995/1000/1005 ...) must be rounded to their nominal shells first
python3 - "$B/bvals" "$O/bvals_shells" <<'PYB'
import sys, numpy as np
b = np.loadtxt(sys.argv[1]); r = np.where(b < 100, 0, np.round(b / 100.0) * 100); np.savetxt(sys.argv[2], r[None], fmt="%d")
print("bvals rounded to shells:", sorted(set(r.astype(int).tolist())))
PYB
cat > $O/fit_bin.m <<M
try
  addpath(genpath('$NODDIDIR')); addpath(genpath('$NIFTIMATLIB'));
  if $RIC
    addpath(genpath('$RDIR'));
    noddi_rician_bin('$O/data.nii', '$O/nodif_brain_mask.nii', '$O/bvals_shells', '$B/bvecs', '$O/noddi');
  else
    CreateROI('$O/data.nii', '$O/nodif_brain_mask.nii', '$O/roi.mat');
    protocol = FSL2Protocol('$O/bvals_shells', '$B/bvecs', 50);
    model = MakeModel('WatsonSHStickTortIsoV_B0');
    % the intrinsic (parallel) diffusivity: the toolbox default 1.7e-9 is the white-matter value;
    % the cortex and subcortex are fitted at 1.1e-9 as cuDIMOT's CorticalNODDI_Watson is (env DPAR, m^2/s)
    diIdx = GetParameterIndex('WatsonSHStickTortIsoV_B0', 'di');
    model.GS.fixedvals(diIdx) = $DPAR; model.GD.fixedvals(diIdx) = $DPAR;
    batch_fitting_single('$O/roi.mat', protocol, model, '$O/params.mat');
    SaveParamsAsNIfTI('$O/params.mat', '$O/roi.mat', '$O/nodif_brain_mask.nii', '$O/noddi');
  end
  fid = fopen('$O/MATLAB_OK', 'w'); fclose(fid);
catch err
  disp(getReport(err)); fid = fopen('$O/FAIL', 'w'); fprintf(fid, '%s\n', getReport(err)); fclose(fid);
end
exit;
M
$MATLABBIN -nodesktop -nosplash -nodisplay -r "run('$O/fit_bin.m')" > $O/matlab.log 2>&1
[ -e $O/MATLAB_OK ] || { [ -e $O/FAIL ] || echo "matlab did not finish" > $O/FAIL; exit 1; }
# cuDIMOT names
for pair in ficvf:mean_fintra fiso:mean_fiso kappa:mean_kappa odi:OD fmin:fmin; do
  src=$O/noddi_${pair%%:*}.nii; [ -e $src ] || continue; fslmaths $src -mas $B/nodif_brain_mask.nii.gz $O/${pair##*:}
done
# the toolbox stores kappa scaled by 1/10 (its odi = 2/pi * atan(1/(10 kappa))): bring it to cuDIMOT's unit
[ -e $O/mean_kappa.nii.gz ] && fslmaths $O/mean_kappa -mul 10 $O/mean_kappa
[ -e $O/noddi_fibredirs_xvec.nii ] && fslmerge -t $O/dyads1 $O/noddi_fibredirs_xvec.nii $O/noddi_fibredirs_yvec.nii $O/noddi_fibredirs_zvec.nii && fslmaths $O/dyads1 -mas $B/nodif_brain_mask.nii.gz $O/dyads1
rm -f $O/data.nii $O/nodif_brain_mask.nii $O/MATLAB_OK
[ -e $O/mean_fintra.nii.gz ] && touch $O/DONE || { echo "no mean_fintra written" > $O/FAIL; exit 1; }
