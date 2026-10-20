#!/bin/bash
# matlab_noddi_bin.sh <rician 0|1> <pool size> <dpar|-> <fit tag> <bin dir> [<bin dir> ...]   2026-09-19
# env: NODDIDIR, NIFTIMATLIB, MATLABBIN, NODDI_TISSUE (GM [default] = cortex and subcortex, di 1.1e-9; WM = white matter, di 1.7e-9)
# NODDI fits with MATLAB of one or more bins (or chunks), in ONE MATLAB session: the data
# is read once and the parpool started once for all of them. Each <bin dir> holds
# data.nii.gz, nodif_brain_mask.nii.gz, bvals, bvecs (the bin's own b-table when it is a
# gradient-nonlinearity bin); the bins share one data file, unpacked by
# NoddiSurfaceMapping.sh as <bin dir>/../data.nii. Writes cuDIMOT-named maps into
# <bin dir>.<fit tag> (mean_fintra, mean_fiso, mean_kappa, OD, dyads1, fmin, AIC, BIC) so
# that the merge and the mapping stages are shared with the cuDIMOT fitters, then DONE (or
# FAIL) there; a bin that fails does not stop the others, and the exit status is 1 if any did.
# The NODDI toolbox (NODDIDIR), WatsonSHStickTortIsoV_B0, batch_fitting on a parpool of
# <pool size> workers (NoddiSurfaceMapping.sh passes its --nthr cap, or 1 with --matlab-fslsub):
#   rician 0 (--fitter=matlab)        : sigmaScaling fixed at 100, sigma = (b0 SD)/100 ('fixed')
#   rician 1 (--fitter=matlab-rician) : sigmaScaling fitted by maximum likelihood ('MLE')
#   dpar : fixed intrinsic diffusivity, m^2/s ('-': the tissue's, 1.1e-9 GM / 1.7e-9 WM)
set -u
RIC="$1"; POOLSIZE="$2"; DPAR="$3"; TAG="$4"; shift 4
[ $# -gt 0 ] || { echo "usage: $0 <rician 0|1> <pool size> <dpar|-> <fit tag> <bin dir> [<bin dir> ...]"; exit 1; }
SIGMASCALING=fixed; [ "$RIC" = 1 ] && SIGMASCALING=MLE
[ "$DPAR" = - ] && DPAR="[]"
NODDIDIR=${NODDIDIR:?export NODDIDIR}; NIFTIMATLIB=${NIFTIMATLIB:?export NIFTIMATLIB}; MATLABBIN=${MATLABBIN:-matlab}; NODDI_TISSUE=${NODDI_TISSUE:-GM}
export FSLDIR=${FSLDIR:-/usr/local/fsl}; export PATH=$FSLDIR/bin:$PATH
# the toolbox reads plain .nii: the data is unpacked once for all the bins, by NoddiSurfaceMapping.sh
# (<bin dir>/../data.nii), and read in place
G=$(dirname $1); DATA=$G/data.nii
[ -e $DATA ] || { for B in "$@" ; do mkdir -p $B.$TAG; echo "ERROR: $DATA not found (NoddiSurfaceMapping.sh unpacks it before the bins)" > $B.$TAG/FAIL; done; exit 1; }

# per bin: the mask as plain .nii. The b-table goes to FSL2Protocol as it is, not rounded to
# shells: FSL2Protocol gives every distinct b its own gradient strength, so each volume keeps its
# exact b (HCP's 990/995/1000/1005, and the per-volume b of a gradient-nonlinearity bin,
# b |(I+L)g|^2, which rounding would undo); b <= 8 is b0. The many "shells" this makes matter
# only to the toolbox's plotting (PlotFittedModel, VoxelDataViewer), not to the fit.
BINS=""; OUTS=""
for B in "$@" ; do
	# a bin is given only when it is not DONE: whatever an earlier run left in its output (a
	# params.mat, finished or not, roi.mat, maps) is removed and the bin is fitted from the start
	O=$B.$TAG; rm -rf $O; mkdir -p $O
	[ -e $O/nodif_brain_mask.nii ] || gunzip -c $B/nodif_brain_mask.nii.gz > $O/nodif_brain_mask.nii
	BINS="$BINS '$B'"; OUTS="$OUTS '$O'"
done

M=$1.$TAG/fit_bins.m; L=$1.$TAG/matlab_session.log	# in the first bin's output: unique to this call (several run at once with --matlab-fslsub)
cat > $M <<EOF
B = {$BINS};
O = {$OUTS};
try
  addpath(genpath('$NODDIDIR')); addpath(genpath('$NIFTIMATLIB'));
  % the data, read once for all the bins, voxels x volumes
  dwi = nifti('$DATA'); D = dwi.dat(:,:,:,:); sz = size(D); D = reshape(D, prod(sz(1:3)), sz(4)); clear dwi
catch err
  for n = 1:numel(O), fid = fopen(fullfile(O{n}, 'FAIL'), 'w'); fprintf(fid, '%s\n', getReport(err)); fclose(fid); end
  disp(getReport(err)); exit;
end
for n = 1:numel(B)
  diary(fullfile(O{n}, 'matlab.log')); diary on
  try
    fprintf('bin %s (%d of %d)\n', B{n}, n, numel(B));
    % the ROI as CreateROI makes it (roi, mask, idx; voxels in its i-j-k loop order), cut from
    % the data in memory instead of reading the data file again for every bin
    msk = nifti(fullfile(O{n}, 'nodif_brain_mask.nii')); msk = msk.dat(:,:,:);
    [I, J, K] = ind2sub(size(msk), find(msk > 0)); idx = sortrows([I J K]);
    lin = sub2ind(size(msk), idx(:,1), idx(:,2), idx(:,3));
    roi = double(D(lin, :)); mask = zeros(size(msk)); mask(lin) = 1:numel(lin);
    roifile = fullfile(O{n}, 'roi.mat'); save(roifile, 'roi', 'mask', 'idx'); clear roi mask
    protocol = FSL2Protocol(fullfile(B{n}, 'bvals'), fullfile(B{n}, 'bvecs'), 8);
    % the modified toolbox (matlab/NODDI-R): the tissue fixes the intrinsic (parallel)
    % diffusivity, GM 1.1e-9 as cuDIMOT's CorticalNODDI_Watson, WM 1.7e-9 as its NODDI_Watson (env NODDI_TISSUE),
    % unless dpar sets it; sigmaScaling (sigma = b0 SD / sigmaScaling) is 'fixed' at 100 or fitted by maximum likelihood ('MLE')
    model = MakeModel('WatsonSHStickTortIsoV_B0', '$NODDI_TISSUE', '$SIGMASCALING', $DPAR);
    % batch_fitting starts the parpool on the first bin and reuses it for the others
    batch_fitting(roifile, protocol, model, fullfile(O{n}, 'params.mat'), $POOLSIZE, false);
    SaveParamsAsNIfTI(fullfile(O{n}, 'params.mat'), roifile, fullfile(O{n}, 'nodif_brain_mask.nii'), fullfile(O{n}, 'noddi'));
    fid = fopen(fullfile(O{n}, 'MATLAB_OK'), 'w'); fclose(fid);
  catch err
    disp(getReport(err)); fid = fopen(fullfile(O{n}, 'FAIL'), 'w'); fprintf(fid, '%s\n', getReport(err)); fclose(fid);
  end
  diary off
end
delete(gcp('nocreate'));
exit;
EOF
$MATLABBIN -nodesktop -nosplash -nodisplay -r "run('$M')" > $L 2>&1

# per bin: cuDIMOT names, then DONE
st=0
for B in "$@" ; do
	O=$B.$TAG
	[ -e $O/MATLAB_OK ] || { [ -e $O/FAIL ] || echo "matlab did not finish (see $L)" > $O/FAIL; st=1; continue; }
	# (the modified toolbox's SaveParamsAsNIfTI already writes kappa in cuDIMOT's unit, odi = 2/pi * atan(1/kappa))
	for pair in ficvf:mean_fintra fiso:mean_fiso kappa:mean_kappa odi:OD fmin:fmin AIC:AIC BIC:BIC; do
		src=$O/noddi_${pair%%:*}.nii; [ -e $src ] || continue; fslmaths $src -mas $B/nodif_brain_mask.nii.gz $O/${pair##*:}
	done
	[ -e $O/noddi_fibredirs_xvec.nii ] && fslmerge -t $O/dyads1 $O/noddi_fibredirs_xvec.nii $O/noddi_fibredirs_yvec.nii $O/noddi_fibredirs_zvec.nii && fslmaths $O/dyads1 -mas $B/nodif_brain_mask.nii.gz $O/dyads1
	rm -f $O/nodif_brain_mask.nii $O/MATLAB_OK
	if [ -e $O/mean_fintra.nii.gz ] ; then touch $O/DONE ; else echo "no mean_fintra written" > $O/FAIL; st=1; fi
done
exit $st
