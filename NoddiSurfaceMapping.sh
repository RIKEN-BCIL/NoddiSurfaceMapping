#! /bin/bash
# =============================================================================
# NoddiSurfaceMapping.sh  v2.0.0-beta
#   NODDI of the cerebral cortex, fitted per vertex on the native cortical
#   surface (or the hippocampal surfaces), and mapped to CIFTI grayordinates.
#
# Author:   Takuya Hayashi, RIKEN BDR Laboratory for Brain Connectomics Imaging
# Citation: Fukutomi H, Glasser MF, Zhang H, et al. (2018) Neurite imaging reveals
#           microstructural variations in human cerebral cortical gray matter.
#           NeuroImage 182:488-499. doi:10.1016/j.neuroimage.2018.02.017
#           Surface fit and CIFTI mapping: Hayashi et al., in prep.;
#           Yoshida et al., in prep.
# License:  MIT (see LICENSE)
# =============================================================================

set -e
NSM_VERSION=2.0.0-beta
CMD=`echo $0 | sed -e 's/^\(.*\)\/\([^\/]*\)/\2/'`

#########################################################
# Usage & Exit
#########################################################

UsageExit () {

 echo ""
 echo " NoddiSurfaceMapping v$NSM_VERSION -- NODDI fitted on the cortical surface, mapped to CIFTI"
 echo " Cite: Fukutomi et al., NeuroImage 2018; Hayashi et al., in prep.; Yoshida et al., in prep."
 echo " Usage: $CMD [options] <StudyFolder> <SubjectID>"
 echo ""
 echo "  What it does by default: one cubic sample of the diffusion signal per vertex on a native-mesh"
 echo "  layer surface a quarter of the way from white to pial (--equi-frac-distance=0.25), NODDI"
 echo "  fitted per vertex, the fitted maps smoothed on that layer surface by its own mean vertex"
 echo "  spacing, and nothing excluded. The maps a judgement would rest on are written out instead"
 echo "  (noddi_snr, noddi_gmpv) so it can be made downstream. See the comments in the script for why"
 echo "  each of these is what it is."
 echo ""
 echo "  Options:"
 echo "    -a <num>  : species atlas (0 Human [default], 1 Macaque hybrid (Mac30BS), 1.1 MacaqueCyno,"
 echo "                1.2 MacaqueRhesus, 1.3 MacaqueSnow, 2 Marmoset, 3 NightMonkey)"
 echo "    -d <str>  : diffusion folder under T1w/ and MNINonLinear/Results/ (default: Diffusion)"
 echo "    -t <a,b,c>: DTI b-value thresholds: upper, lower, b=0 upper (default 1050,100,50)"
 echo "    -M        : RegName=MSMAll (default MSMSulc)"
 echo "    -n <num>  : NJOBS passed to cuDIMOT (default 4). Unrelated to --nthr"
 echo "    -S <num>  : flag vertices whose b0 SNR is below <num> x the cortical median as dropout"
 echo "                (orbitofrontal, inferior temporal), where the fit inflates NDI. DEFAULT 0 =" 
 echo "                flag nothing; 0.4 flags 3-7 % of the cortex in HCP-YA."
 echo "                noddi_snr and noddi_gmpv are written whatever this is set to."
 echo "    -N        : leave flagged vertices as NaN instead of filling them from their neighbours"
 echo ""
 echo "    --mapping=<what> : fit route and (for the surface route) the ordinate target"
 echo "                cortex                     : cortical surface only; outputs named *_cortex"
 echo "                grayordinate [default]     : cortical surface + subcortical grey of wmparc"
 echo "                grayordinate+whiteordinate : + the white matter as volume voxels (\$CUDIMOTDIR_WM)"
 echo "                hippunfold                 : hippocampal (HippUnfold) surfaces, not the cortex;"
 echo "                     hipp+dentate, finest density fitted then downsampled to the rest on the"
 echo "                     surface; --equi-frac-distance defaults to 0.5; needs wb_command >= 2.1.0."
 echo "                     Folder <Study>/<Subject>/T1w/HippUnfold/T1w; outputs under"
 echo "                     MNINonLinear/Results/<dwi>/HippUnfold/<den>/<S>.hippocampus_<map>_<tag>*.dscalar.nii"
 echo "                volume                     : obsolete whole-brain volume fit + myelin-style mapping"
 echo "                                             (Fukutomi et al., 2018)"
 echo "    --equi-frac-distance=<f> : sampling depth = fraction f of white->pial (0 white, 1 pial);"
 echo "                DEFAULT 0.25, 0.5 for --mapping=hippunfold (inner->outer).         tag fracdist<f>"
 echo "    --equi-vol-distance=<f>  : equivolume fraction f from white.                    tag voldist<f>"
 echo "    --equi-abs-distance=<mm> : <mm> inside white along its normal, capped at --abs-distance-cap"
 echo "                (default 0.25) x local thickness.  tag absdist<mm>   (at most one of the three)"
 echo "    --tmpdir=<dir>  : where data.nii.gz is gunzipped once for the sampling (default the fit folder)"
 echo "    --runmode=<prep,fit,map> : which stages to run (default all three, any comma-separated subset)"
 echo "    --fitter=<name> : cudimot-mcmc [default] MCMC with Rician noise; cudimot-mle Levenberg-"
 echo "                Marquardt only; matlab the NODDI toolbox; matlab-rician a Rician NODDI in MATLAB (--rnoddidir)."
 echo "                Each writes its own folder and CIFTI suffix"
 echo "    --gradnonlin=<file> : gradient nonlinearity tensor (grad_dev.nii.gz, 9 volumes, on the"
 echo "                data grid). Every DTIFit gets it, and the NODDI fits are done in bins of the"
 echo "                effective b-value scale |(I+L)g|^2, each bin with its own b-table. Every"
 echo "                fitter honours it, MATLAB included, and the bins are also the unit of"
 echo "                parallelism"
 echo "    --gnlbin=<num>  : width of those bins (default 0.01; bins under 200 vertices are merged)"
 echo "    --nchunk=<num>  : how many equal chunks to split the data into for the MATLAB fitters"
 echo "                when --gradnonlin is NOT used (default 20). With it, the bins are the chunks"
 echo "    --matlab=<bin>  : MATLAB executable (default matlab)"
 echo "    --noddidir=<dir>: NODDI toolbox (default <this dir>/NODDI/NODDI_toolbox)"
 echo "    --rnoddidir=<dir>: the Rician NODDI implementation, for --fitter=matlab-rician: a folder with"
 echo "                noddi_rician_bin(datafile, maskfile, bvalsfile, bvecsfile, outprefix), writing"
 echo "                <outprefix>_ficvf/_fiso/_kappa/_odi/_fibredirs_{x,y,z}vec/_fmin.nii as the toolbox does"
 echo "    --nthr=<num>    : cap the CPU threads this run may use (default: the queue NSLOTS, or 1"
 echo "                outside a queue; 0 lifts the cap). wb_command takes every core it finds"
 echo "    --no-dti  : skip DTIFit and DiffusionStats (volume measures the surface fit does not read)"
 echo ""

 exit 1

}

CheckPath () {
# HCP PIPELINE
if [ -z $HCPPIPEDIR ] ; then
        echo "Please export HCPPIPEDIR"
        exit 1
fi
if [ -z $EnvironmentScript ] ; then
        echo "Please export EnvironmentScript"
        exit 1
fi
}

if [ "$2" = "" ] ; then CheckPath; UsageExit; fi

#########################################################
# Setup
#########################################################
RunModes="prep,fit,map"
Fitter="cudimot-mcmc"
NChunk=20
MATLABBIN="matlab"
NODDIDIR="${NODDIDIR:-}"
RNODDIDIR="${RNODDIDIR:-}"
REGNAME="MSMSulc"
Species="0"
THR="1050,100,50" #b-value upper and lower threshold, b=0 upper threshold for DTI
MappingMethod="surface"		# surface (fit per vertex) | volume (obsolete: volume fit + myelin-style mapping)
Coord="grayordinate"	# cortex | grayordinate | grayordinate+whiteordinate | hippunfold
MappingArg=""			# raw --mapping value (resolved to MappingMethod + Coord below)
CortexOnly="NO"		# both derived from --coord below
WMFit="NO"
NJOB=4
GradNonlin=""
GnlBinWidth=0.01
DoDTI="YES"
NThr=${NSLOTS:-1}
DiffusionName="Diffusion"
SNRfrac=0
KeepNaN="NO"
# sampling depth of the surface fit (2026-09-30): mode frac|abs|vol, value as given
DepthMode=""; DepthValue=""; DepthCap=""; DepthCount=0
TmpDir=""
# --coord=hippunfold (2026-10): fit NODDI per vertex on the hippocampal
# (HippUnfold) surfaces at the finest density, then downsample on the surface to
# the coarser ones. The HippUnfold space folder is <Study>/<Subject>/T1w/HippUnfold/T1w;
# hipp and dentate are both done, and every density present there is written.
HippUnfoldDir=""; HippDens=""; HippStructs="hipp dentate"
export LC_ALL=C  # suppress warning of locale

# long options first (--gradnonlin=<file>, --gnlbin=<num>), then the single-letter ones
ARGS=()
for arg in "$@"; do
	case "$arg" in
		--mapping=*)    MappingArg="${arg#*=}";;
		--gradnonlin=*) GradNonlin="${arg#*=}";;
		--no-dti) DoDTI="NO";;
		--nthr=*) NThr="${arg#*=}";;
		--gnlbin=*)     GnlBinWidth="${arg#*=}";;
		--runmode=*)    RunModes="${arg#*=}";;
		--fitter=*)     Fitter="${arg#*=}";;
		--nchunk=*)     NChunk="${arg#*=}";;
		--matlab=*)     MATLABBIN="${arg#*=}";;
		--noddidir=*)   NODDIDIR="${arg#*=}";;
		--rnoddidir=*)  RNODDIDIR="${arg#*=}";;
		--equi-frac-distance=*) DepthMode=frac; DepthValue="${arg#*=}"; DepthCount=$((DepthCount+1));;
		--equi-abs-distance=*)  DepthMode=abs;  DepthValue="${arg#*=}"; DepthCount=$((DepthCount+1));;
		--equi-vol-distance=*)  DepthMode=vol;  DepthValue="${arg#*=}"; DepthCount=$((DepthCount+1));;
		--abs-distance-cap=*)   DepthCap="${arg#*=}";;
		--tmpdir=*)     TmpDir="${arg#*=}";;
		--mapping|--gradnonlin|--gnlbin|--runmode|--fitter|--nchunk|--matlab|--noddidir|--rnoddidir|--equi-frac-distance|--equi-abs-distance|--equi-vol-distance|--abs-distance-cap|--tmpdir) echo "ERROR: $arg needs =<value>"; exit 1;;
		--*) echo "ERROR: unknown option $arg"; exit 1;;
		*) ARGS+=("$arg");;
	esac
done
set -- "${ARGS[@]}"

# The thread cap. wb_command -volume-to-surface-mapping is the stage that costs:
# it maps every diffusion volume onto the native midthickness and, left alone,
# runs on every core of the node. A batch of these submitted one slot each then
# oversubscribes the machine by whatever factor OpenMP chose -- measured at
# 12x on 2026-09-23, which put the GPU queue into alarm and stopped it
# scheduling anything more. numpy's BLAS does the same in the packing step.
#
# The default follows NSLOTS so that the run takes what it asked for and no
# more. One slot per task is usually the better trade: the stage scales poorly
# (four threads measured about three times one), so ten single-threaded tasks
# on ten slots finish more work per hour than two four-threaded ones.
if [ "${NThr:-0}" -gt 0 ] 2>/dev/null ; then
	export OMP_NUM_THREADS=$NThr OPENBLAS_NUM_THREADS=$NThr MKL_NUM_THREADS=$NThr
	export NUMEXPR_NUM_THREADS=$NThr
fi
while getopts "Ma:t:n:d:S:N" OPT
do
	case "$OPT" in
		"a" ) Species="$OPTARG";;
		"M" ) REGNAME="MSMAll";;
		"t" ) THR="$OPTARG";;
		"n" ) NJOB="$OPTARG";;
		"d" ) DiffusionName="$OPTARG";;
		"S" ) SNRfrac="$OPTARG";;
		"N" ) KeepNaN="YES";;
		* )  UsageExit;;
	esac
done
shift `expr $OPTIND - 1`

# --mapping is the single route/target selector: it carries both the fit route
# (surface|volume) and, for the surface route, the ordinate target that used to be
# --coord (cortex|grayordinate|grayordinate+whiteordinate|hippunfold). --coord is
# kept as a deprecated alias for those targets.
case "$MappingArg" in
	"" )                                                   ;;   # not given: keep defaults / --coord
	surface | volume )                MappingMethod="$MappingArg" ;;
	cortex | grayordinate | grayordinate+whiteordinate | hippunfold ) MappingMethod=surface; Coord="$MappingArg" ;;
	* ) echo "ERROR: --mapping must be surface, volume, cortex, grayordinate, grayordinate+whiteordinate or hippunfold"; exit 1 ;;
esac
case "$MappingMethod" in surface|volume) ;; *) echo "ERROR: --mapping must be surface or volume"; exit 1;; esac
# The sampling depth (2026-09-30). One of the three options; none = equidistant
# 0.25 from white, the depth a 5-subject HCP-YA pilot chose (depth_pilot2). The
# tag carries the number as given, so that it can be matched by eye to the
# command line, and names every output of the surface route.
if [ "$DepthCount" -gt 1 ] ; then echo "ERROR: --equi-frac-distance, --equi-abs-distance and --equi-vol-distance are mutually exclusive"; exit 1; fi
# default sampling depth: 0.25 of white->pial for cortex; the midthickness (0.5
# between HippUnfold inner and outer) for the hippocampal surfaces
[ -n "$DepthMode" ] || { DepthMode=frac; DepthValue=0.25; [ "$Coord" = hippunfold ] && DepthValue=0.5; }
_isnum () { echo "$1" | grep -Eq '^([0-9]+\.?[0-9]*|\.[0-9]+)$'; }
_isnum "$DepthValue" || { echo "ERROR: the sampling depth must be a non-negative number (got '$DepthValue')"; exit 1; }
if [ -n "$DepthCap" ] ; then
	[ "$DepthMode" = abs ] || { echo "ERROR: --abs-distance-cap goes with --equi-abs-distance only"; exit 1; }
	_isnum "$DepthCap" && awk -v c="$DepthCap" 'BEGIN{exit !(c>0 && c<=1)}' || { echo "ERROR: --abs-distance-cap must be in (0,1] (got '$DepthCap')"; exit 1; }
fi
case "$DepthMode" in
	frac|vol) awk -v f="$DepthValue" 'BEGIN{exit !(f>=0 && f<=1)}' || { echo "ERROR: --equi-$DepthMode-distance must be in [0,1] (got $DepthValue)"; exit 1; }
		DepthTag=${DepthMode}dist${DepthValue} ;;
	abs)	awk -v f="$DepthValue" 'BEGIN{exit !(f>0)}' || { echo "ERROR: --equi-abs-distance must be > 0 mm (got $DepthValue)"; exit 1; }
		DepthTag=absdist${DepthValue}
		# a cap other than the default changes the surface, so it changes the name
		if [ -n "$DepthCap" ] && awk -v c="$DepthCap" 'BEGIN{exit !(c!=0.25)}' ; then DepthTag=${DepthTag}cap${DepthCap}; fi
		[ -n "$DepthCap" ] || DepthCap=0.25 ;;
esac
if [ "$MappingMethod" = volume ] ; then
	[ "$DepthCount" -gt 0 ] && echo "WARNING: the sampling-depth options apply to --mapping=surface only; ignored"
	DepthTag=""; NTag=""
else
	NTag="_${DepthTag}"
fi
# which ordinates the fit covers, in the CIFTI sense
case "$Coord" in
	cortex)                     CortexOnly=YES; WMFit=NO  ;;
	grayordinate)               CortexOnly=NO;  WMFit=NO  ;;
	grayordinate+whiteordinate) CortexOnly=NO;  WMFit=YES ;;
	hippunfold)                 CortexOnly=NO;  WMFit=NO  ;;   # hippocampal surfaces (HippUnfold)
	*) echo "ERROR: --mapping must be cortex, grayordinate, grayordinate+whiteordinate, hippunfold or volume"; exit 1;;
esac
case "$Fitter" in cudimot-mcmc|cudimot-mle|matlab|matlab-rician) ;; *) echo "ERROR: --fitter must be cudimot-mcmc, cudimot-mle, matlab or matlab-rician"; exit 1;; esac
for st in $(echo $RunModes | tr ',' ' ') ; do case "$st" in prep|fit|map) ;; *) echo "ERROR: --runmode stages are prep, fit, map (got $st)"; exit 1;; esac; done
DoPrep=NO; DoFit=NO; DoMap=NO
case ",$RunModes," in *,prep,*) DoPrep=YES;; esac; case ",$RunModes," in *,fit,*) DoFit=YES;; esac; case ",$RunModes," in *,map,*) DoMap=YES;; esac

SetUp () {

CheckPath
source $EnvironmentScript
# Workbench: CARET7DIR as the HCP environment script sets it. The layer
# construction needs -surface-average -weight and -surface-cortex-layer (>= 1.5).
"${CARET7DIR:-/nonexistent}"/wb_command -version > /dev/null 2>&1 || { echo "ERROR: no working wb_command (CARET7DIR=${CARET7DIR:-unset})"; exit 1; }
# --coord=hippunfold needs the native HIPPOCAMPUS[_DENTATE] CIFTI structures, which
# Connectome Workbench added in v2.1.0. Prefer a >=2.1.0 build here; fail if none.
if [ "$Coord" = hippunfold ] ; then
	_wbver () { "$1"/wb_command -version 2>/dev/null | awk '/^Version/{print $2; exit}'; }
	_wbge210 () { awk -v v="$1" 'BEGIN{n=split(v,a,"."); exit !((a[1]>2)||(a[1]==2&&a[2]>=1))}'; }
	if ! _wbge210 "$(_wbver "$CARET7DIR")" ; then
		for d in "${WORKBENCHDIR_HIPP:-}" ; do
			[ -n "$d" ] && "$d"/wb_command -version >/dev/null 2>&1 && _wbge210 "$(_wbver "$d")" && { export CARET7DIR=$d; break; }
		done
	fi
	_wbge210 "$(_wbver "$CARET7DIR")" || { echo "ERROR: --mapping=hippunfold needs wb_command >= 2.1.0 (HIPPOCAMPUS_DENTATE CIFTI structures); found $(_wbver "$CARET7DIR") at $CARET7DIR. Set WORKBENCHDIR_HIPP=<dir with a >=2.1.0 wb_command>."; exit 1; }
	echo "HippUnfold: using wb_command $(_wbver "$CARET7DIR") at $CARET7DIR"
fi
# the scripts it calls (dwistats, the MATLAB
# wrapper) are used in place: NODDIHCP may point at the shared install while
# this file lives elsewhere (NODDIHCP_SCRIPTS_DIR overrides)
NODDIHCP=$(dirname $(readlink -f $0))				# path to NoddiSurfaceMapping
[ -d $NODDIHCP/scripts ] || NODDIHCP=$(dirname $NODDIHCP)	# (a copy under scripts/, e.g. the .next being tested)
[ -n "${NODDIHCP_SCRIPTS_DIR:-}" ] && NODDIHCP=$NODDIHCP_SCRIPTS_DIR

# The fitters. cuDIMOT (MCMC or Levenberg-Marquardt) on the GPU, or MATLAB on the
# CPU: the NODDI toolbox, or a Rician NODDI supplied with --rnoddidir. FitTag names the fit folder
# (<fit dir>.<FitTag>) and FitSuffix the CIFTI (noddi_ficvf<OutSuffix><FitSuffix>...),
# so that the fitters' results coexist.
case "$Fitter" in
	cudimot-mcmc)  FitTag=CorticalNODDI_Watson;       FitSuffix="";;
	cudimot-mle)   FitTag=CorticalNODDI_Watson_MLE;   FitSuffix="_mle";;
	matlab)        FitTag=MatlabNODDI;                FitSuffix="_matlab";;
	matlab-rician) FitTag=MatlabRicianNODDI;          FitSuffix="_matlabrician";;
esac
# --gradnonlin gives a different answer from the nominal b-table, so it names
# its outputs differently (2026-09-23). Without this the two share every CIFTI
# name and a run with the option silently overwrites the run without it, which
# is exactly the pair a comparison needs to keep. The suffix matches the name
# the separate surffit_gnl.sh wrote before the option was folded in here, so
# the CIFTIs of the two routes are interchangeable.
[ -n "$GradNonlin" ] && { FitTag="${FitTag}_gnl"; FitSuffix="${FitSuffix}_gnl"; }
[ -n "$NODDIDIR" ] || NODDIDIR=${NODDIHCP}/NODDI/NODDI_toolbox
NIFTIMATLIB=${NIFTIMATLIB:-${NODDIHCP}/NODDI/niftimatlib/matlab}
if [ "$Fitter" = matlab ] && [ ! -d $NODDIDIR/fitting ] ; then echo "ERROR: NODDI toolbox not found at $NODDIDIR (--noddidir)"; exit 1; fi
if [ "$Fitter" = matlab-rician ] ; then
	[ -f "$RNODDIDIR/noddi_rician_bin.m" ] || { echo "ERROR: --fitter=matlab-rician needs --rnoddidir=<dir> containing noddi_rician_bin.m (got '${RNODDIDIR}')"; exit 1; }
fi

# CUDIMOT-CorticalNODDI
# requires compilation of CUDIMOT Watson NODDI (w/ modelname=CorticalNODDI) with default value of parallel diffusivity set to 1.1
# also need to set compute mode to Exclusive_Process (e.g. 'nvidia-smi -c 3') and persistent mode to ON (e.g. 'nvidia-smi -pm 1') to efficiently run GPU in parallel  

CUDIMOTMODELNAME=CorticalNODDI_Watson
CUDIMOTLD_LIBRARY_PATH=""
CUDIMOTDIR="${CUDIMOTDIR:-/usr/local/CUDIMOT}"
# the white-matter build: parallel diffusivity 1.7 (model name NODDI_Watson in that tree)
CUDIMOTDIR_WM="${CUDIMOTDIR_WM:-/usr/local/CUDIMOT_dpar1p7}"
CUDIMOTMODELNAME_WM="${CUDIMOTMODELNAME_WM:-NODDI_Watson}"
echo "========================================="
echo "START NoddiSurfaceMapping.sh v$NSM_VERSION"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Host:$(hostname) Job:${JOB_ID:-N/A} GPU:${CUDA_VISIBLE_DEVICES:-not set}"
echo "========================================="
# --- Site-specific GPU assignment (example, disabled) ---------------------------
# For an SGE whose GPU prolog writes "export CUDA_VISIBLE_DEVICES=<id>" to
# /tmp/sge_gpu_env_<host>_<job>_<task> instead of setting it for the job.
# Uncomment and adapt if your scheduler works that way; otherwise every fit runs
# on the devices CUDA_VISIBLE_DEVICES already names (or all of them).
#export SGE_ROOT=/usr/local/sge
#
#export FSLGECUDAQ="cuda.q"
#export PATH=$SGE_ROOT/bin/lx-amd64:$PATH
#
## --- GPU Configuration (SGE Prolog Integration) ---
## Load the GPU device assignment determined by the SGE Prolog script.
## In a shared /tmp environment, we use hostname-specific files to prevent conflicts between BCIL-G11 and BCIL-G12.
## added 2026/01/19
## Path must match the ENV_FILE defined in the Prolog/Epilog scripts
#
#MY_HOST="$(hostname)"
#ts="$(date '+%Y-%m-%d %H:%M:%S')"
#queue_name="${QUEUE%@*}"
#log_prefix="[$ts] Host:$MY_HOST Job:$JOB_ID Queue:${QUEUE:-N/A}"
#
#echo "========================================="
#echo "START NoddiSurfaceMapping.sh"
#echo "$log_prefix"
#echo "========================================="
#
## ---- queue-specific GPU handling ----
#if [[ "$queue_name" == "cuda.q" ]]; then
#  # the prolog names the file with the array task id (2026-09: array-job aware);
#  # the task-less name is what it wrote before that
#  TASK_ID="${SGE_TASK_ID:-0}"; [ "$TASK_ID" = "undefined" ] && TASK_ID=0
#  ENV_FILE="/tmp/sge_gpu_env_${MY_HOST}_${JOB_ID}_${TASK_ID}"
#  [[ -f "$ENV_FILE" ]] || ENV_FILE="/tmp/sge_gpu_env_${MY_HOST}_${JOB_ID}"
#
#  if [[ -f "$ENV_FILE" ]]; then
#    source "$ENV_FILE"
#    if [[ -z "${CUDA_VISIBLE_DEVICES:-}" ]]; then
#    echo "FATAL: ENV_FILE loaded but CUDA_VISIBLE_DEVICES is empty. Aborting."
#    exit 101
#    fi
#    echo "Running on cuda.q | Assigned GPU: ${CUDA_VISIBLE_DEVICES:-UNSET}"
#  else
#    echo "Running on cuda.q | WARNING: GPU environment file not found ($ENV_FILE)"
#    echo "Check if the SGE Prolog script is properly configured for this queue."
#    echo "FATAL: cuda.q job started without GPU assignment. Aborting to avoid mis-execution."
#    exit 100
#  fi
#
#else
#  echo "NOT on cuda.q | GPU will NOT be assigned by SGE prolog in this queue."
#fi
## --------------------------------------------------

# NJOB specifies the number of data splits in $CUDIMOT/bin/jobs_wrapper.sh.
# cuDIMOT's jobs_wrapper.sh runs the parts one after another on the same GPU, so
# increasing this value does not speed up a local run.
# However, caution is required when N=1: if the data is not split, the entire 
# dataset is loaded at once, which may lead to a GPU Out of Memory (OOM) error.

# References:
# NODDI: http://mig.cs.ucl.ac.uk/index.php?n=Tutorial.NODDImatlab
# CUDIMOT:https://users.fmrib.ox.ac.uk/~moisesf/cudimot
# Surface fit (2026-09-17): the diffusion signal is mapped onto the cortical
# surface first and NODDI is fitted per vertex, so that no voxel straddling
# white matter or CSF enters the fit (the idea of Fukutomi et al. 2018, with
# the mapping HCP's fMRISurface uses today instead of the myelin-style one).

# Folder Name
T1wFOLDER="T1w"
DWIT1wFOLDER="$DiffusionName"
DWINativeFOLDER="$DiffusionName"
AtlasSpaceNativeFOLDER="Native"
AtlasSpaceFOLDER="MNINonLinear"
AtlasSpaceResultsDWIFOLDER="$DWIT1wFOLDER"
AtlasSpaceFOLDER="MNINonLinear"

FreeSurferSubjectFolder=$StudyFolder/$Subject/$T1wFOLDER
FreeSurferSubjectID=$Subject
DWINativeFolder=$StudyFolder/$Subject/$DWINativeFOLDER
DtiRegDir=$DWINativeFolder/reg
T1wFolder=$StudyFolder/$Subject/$T1wFOLDER
DWIT1wFolder=$T1wFolder/$DWIT1wFOLDER
T1wNativeFolder=$T1wFolder/$AtlasSpaceNativeFOLDER
# --coord=hippunfold: the HippUnfold space folder (same space as data.nii.gz) and
# the densities to write, auto-detected from the <den>/ subfolders present there
if [ "$Coord" = hippunfold ] ; then
	HippUnfoldDir=$T1wFolder/HippUnfold/$T1wFOLDER
	HippDens=$(for d in "$HippUnfoldDir"/*/ ; do den=$(basename "$d")
		[ -e "$HippUnfoldDir/$den/$Subject.L.hipp_midthickness.$den.surf.gii" ] && echo "$den"
	done | awk '{n=$1; sub(/k$/,"",n); if($1 ~ /k$/) n=n*1000; print n" "$1}' | sort -rn | awk '{printf "%s ",$2}')
fi

AtlasSpaceFolder=$StudyFolder/$Subject/$AtlasSpaceFOLDER
AtlasSpaceNativeFolder=$AtlasSpaceFolder/$AtlasSpaceNativeFOLDER
AtlasSpaceResultsDWIFolder=$AtlasSpaceFolder/Results/$AtlasSpaceResultsDWIFOLDER
AtlasSpaceFolder=$StudyFolder/$Subject/$AtlasSpaceFOLDER

# Surface mapping
ribbonLlabel=3
ribbonRlabel=42
ROIFolder=$AtlasSpaceFolder/ROIs
# where the per-vertex maps live before resampling: one folder per method,
# so that old (myelin-style) and new (surface fit) results never overwrite each other
if [ "$MappingMethod" = "surface" ] ; then
	# every surface-route folder carries the sampling depth (2026-09-30); the
	# subcortical fit does not depend on it and stays shared
	MapFolder=$AtlasSpaceResultsDWIFolder/SurfaceFitMapping${FitSuffix}${NTag}
	SurfFitFolder=${DWIT1wFolder}.SurfaceNODDI${NTag}
	SubFitFolder=${DWIT1wFolder}.SubcorticalNODDI
else
	MapFolder=$AtlasSpaceResultsDWIFolder/RibbonVolumeToSurfaceMapping${FitSuffix}
	SurfFitFolder=""
fi
WMFitFolder=${DWIT1wFolder}.WhiteMatterNODDI
# output suffix: the same names as before when the subcortical volume is in the
# CIFTI, *_cortex when it is not
if [ "$CortexOnly" = "YES" ] ; then OutSuffix="_cortex" ; else OutSuffix="" ; fi
# and the sampling depth (surface route only): noddi_ficvf_fracdist0.25_gnl.32k_fs_LR.dscalar.nii
OutSuffix="${OutSuffix}${NTag}"
# the layer surface the surface route samples on, per hemisphere (see _make_layer)
LayerSurf () { echo "$T1wNativeFolder/$Subject.$1.${DepthTag}.native.surf.gii"; }
TmpWork=${TmpDir:-$SurfFitFolder}

# Which species
case $Species in 
 0)	export SPECIES=Human
	;;
 1)	export SPECIES=MacaqueMac30BS
	;;
 1.1)	export SPECIES=MacaqueCyno
	;;
 1.2)	export SPECIES=MacaqueRhesus
	;;
 1.3)	export SPECIES=MacaqueSnow
	;;
 2)	export SPECIES=Marmoset
	;;
 3)	export SPECIES=NightMonkey
	;;
 *)	echo "Not yet supportted atlas species: $Species"; exit 1
esac

#source $EnvironmentScript
source $HCPPIPEDIR/global/scripts/log.shlib  # Logging related functions
if [ ! "$SPECIES" = Human ] ; then
	source $HCPPIPEDIR/Examples/Scripts/SetUpSPECIES.sh $SPECIES
else
	LowResMeshes="32"
	GrayordinatesResolutions="1.60"
	GrayordinatesSpaceDIR=$HCPPIPEDIR/global/templates/standard_mesh_atlases
fi
DiffRes="`fslval $DWIT1wFolder/data.nii.gz pixdim1 | awk '{printf "%0.2f",$1}'`"
NODDIMappingFWHM="`echo "$DiffRes * 2.5" | bc -l`"
NODDIMappingSigma=`echo "$NODDIMappingFWHM / ( 2 * ( sqrt ( 2 * l ( 2 ) ) ) )" | bc -l`
SmoothingFWHM="$DiffRes"
SmoothingSigma=`echo "$SmoothingFWHM / ( 2 * ( sqrt ( 2 * l ( 2 ) ) ) )" | bc -l`
LowResMeshes=(`echo $LowResMeshes | sed -e 's/@/ /g'`)
GrayordinatesResolutions=(`echo $GrayordinatesResolutions | sed -e 's/@/ /g'`)
SNRthr=17     # Fukutomi et al., NeuroImage 2018
SNRthr=0

if [ ! -e "$AtlasSpaceFolder" ] ; then
 echo "Error: Cannot find $AtlasSpaceFolder"; exit 1;
fi

if [ ! "`imtest $AtlasSpaceFolder/ribbon.nii.gz`" = 1 ] ; then
 echo "ERROR: cannot find ribbon.nii.gz in $AtlasSpaceFolder"; exit 1;
fi

if [ "$REGNAME" = "MSMAll" ] ; then
	Reg="_MSMAll"
else 
	Reg=""
fi
log_Msg "HCPPIPDIR=$HCPPIPEDIR"
log_Msg "Diffusion Resolution: $DiffRes"
if [ -n "$GradNonlin" ] && [ `imtest $GradNonlin` -eq 0 ] ; then echo "ERROR: --gradnonlin file not found: $GradNonlin"; exit 1; fi
GNLOPT=""; [ -n "$GradNonlin" ] && GNLOPT="--gradnonlin=$GradNonlin"
log_Msg "MappingMethod: $MappingMethod   Coord: $Coord   OutSuffix: '$OutSuffix'   SNRfrac: $SNRfrac   KeepNaN: $KeepNaN"
[ "$MappingMethod" = surface ] && log_Msg "Sampling depth: $DepthMode $DepthValue${DepthCap:+ (cap $DepthCap)}   tag: $DepthTag   layer: $(LayerSurf L|sed "s/\.L\./.<H>./")   wb_command: $CARET7DIR"
log_Msg "GradNonlin: ${GradNonlin:-none}   bin width: $GnlBinWidth"
log_Msg "RunModes: $RunModes   Fitter: $Fitter (FitTag $FitTag, CIFTI suffix '$FitSuffix')   NChunk: $NChunk"
log_Msg "NODDIMappingFWHM (myelin-style only): $NODDIMappingFWHM"
log_Msg "LowResMeshes: $LowResMeshes"
log_Msg "GrayordinatesResolutsions: $GrayordinatesResolutions"
}

#########################################################
# Calculate DTI and NODDI and do surface mapping
#########################################################
DTIFit () {

log_Msg "Start: DTIFit"

thr=(`echo $THR | sed -e 's/,/ /g'`)
buthresh="${thr[0]}" # b-value upper threshold for DTI
blthresh="${thr[1]}" # b-value lower threshold for DTI
b0thresh="${thr[2]}" # b0 volume threshold for DTI b0

log_Msg "b-value upper threshhold: $buthresh"
log_Msg "b-value lower threshhold: $blthresh"
log_Msg "b=0 volume threshhold: $b0thresh"

if [ -e $DWIT1wFolder/dti_dwi.txt ] ; then rm $DWIT1wFolder/dti_dwi.txt; fi
j=0;for i in `cat $DWIT1wFolder/bvals` ; do if [ `echo $i | awk '{printf "%d",$1}'` -le $buthresh ] && [ `echo $i | awk '{printf "%d",$1}'` -ge $blthresh ] ; then j=`zeropad $j 4`; echo -n "vol${j} ">> $DWIT1wFolder/dti_dwi.txt;fi;j=`expr $j + 1`;done


if [ ! -e $DWIT1wFolder/dti_dwi.txt ] ; then
  echo "ERROR: cannnot create dti_dwi.txt. Consider increasing b-value upper threshold (default=1050)"
  exit 1;
else
  log_Msg "Number of dwi volumes=`cat $DWIT1wFolder/dti_dwi.txt | wc -w`"
fi 

if [ -e $DWIT1wFolder/dti_b0.txt ] ; then rm $DWIT1wFolder/dti_b0.txt; fi
j=0;for i in `cat $DWIT1wFolder/bvals` ; do if [ `echo $i | awk '{printf "%d",$1}'` -le $b0thresh ] ; then j=`zeropad $j 4`; echo -n "vol${j} ">> $DWIT1wFolder/dti_b0.txt;fi;j=`expr $j + 1`;done
if [ ! -e $DWIT1wFolder/dti_b0.txt ] ; then
  echo "ERROR: cannnot create dti_dwi.txt. Consider increasing b=0 upper threshold (default=50)"
  exit 1;
else
  log_Msg "Number of b0 volumes=`cat $DWIT1wFolder/dti_b0.txt | wc -w`"
fi

if [ -e $DWIT1wFolder/dti_vol.txt ] ; then rm $DWIT1wFolder/dti_vol.txt; fi
for i in $(cat $DWIT1wFolder/dti_b0.txt $DWIT1wFolder/dti_dwi.txt) ; do echo $i; done | sort | awk '{printf "%s ",$1}' > $DWIT1wFolder/dti_vol.txt

if [ -e $DWIT1wFolder/dti_bvecs ] ; then rm $DWIT1wFolder/dti_bvecs;fi
touch $DWIT1wFolder/dti_bvecs;
if [ -e $DWIT1wFolder/dti_bvals ] ; then rm $DWIT1wFolder/dti_bvals;fi
touch $DWIT1wFolder/dti_bvals;

for i in `cat $DWIT1wFolder/dti_vol.txt`; do
  j=`echo $i | sed -e 's/vol//g'`; j=`expr $j + 1`
  cat $DWIT1wFolder/bvals | awk '{printf "%f ", '$`echo $j`'}' >> $DWIT1wFolder/dti_bvals
  cat $DWIT1wFolder/bvecs | awk '{printf "%f \n", '$`echo $j`'}' > $DWIT1wFolder/dti_bvecstmp
  paste $DWIT1wFolder/dti_bvecs $DWIT1wFolder/dti_bvecstmp >  $DWIT1wFolder/dti_bvecstmp2
  mv $DWIT1wFolder/dti_bvecstmp2 $DWIT1wFolder/dti_bvecs
done

fslsplit $DWIT1wFolder/data $DWIT1wFolder/vol
if [ -e $DWIT1wFolder/dti_vollist.txt ] ; then rm $DWIT1wFolder/dti_vollist.txt;fi
for i in `cat $DWIT1wFolder/dti_vol.txt`; do echo $DWIT1wFolder/$i >> $DWIT1wFolder/dti_vollist.txt; done
fslmerge -t $DWIT1wFolder/dti_data `cat $DWIT1wFolder/dti_vollist.txt`
dtifit  -k $DWIT1wFolder/dti_data.nii.gz -o $DWIT1wFolder/dti -m $DWIT1wFolder/nodif_brain_mask -r $DWIT1wFolder/dti_bvecs -b $DWIT1wFolder/dti_bvals --sse --save_tensor $GNLOPT
imrm $DWIT1wFolder/vol????.nii.gz
rm $DWIT1wFolder/dti_bvecstmp $DWIT1wFolder/dti_vollist.txt

}

#########################################################
# cuDIMOT: one fit of one directory (data, bvals, bvecs, nodif_brain_mask)
#   _cudimot <subjectdir> <cudimot install dir> <model name>
#########################################################
# _cudimot_fit <dir> <cuDIMOT dir> <model> <njobs> <grad_dev or ""> [final-step options]
# cuDIMOT's three-step NODDI fit, run in place, one stage after another. The
# steps, their options and their initialisation are those of cuDIMOT's own
# Pipeline_<model>.sh; they are written out here because that script hands every
# stage to fsl_sub when SGE is present -- which, called from inside a GPU job,
# queues further GPU jobs behind the one holding the device -- and because conda
# builds of cuDIMOT do not install it. Of cuDIMOT only its programs are called:
# split_parts_<model>, <model>, merge_parts_<model> and cart2spherical. The DTI
# that initialises the orientation is run here too (what cuDIMOT's Run_dtifit.sh
# does), so that with --gradnonlin it gets the tensor per voxel and the nominal
# b-table, exactly as the production fits did; with the bin's b-table instead
# the start moves by a median 1.4 deg and ODI by 0.0025 on average (HCP-YA).
_cudimot_step () {	# <bin> <model> <njobs> <outdir> <options...>
	local bin="$1" model="$2" nj="$3" W="$4"; shift 4
	local io="--outputdir=$W --partsdir=$W/diff_parts --nParts=$nj"
	mkdir -p $W/diff_parts $W/logs
	$bin/split_parts_$model "$@" $io --idPart=0 --logdir=$W/logs/preProcess
	local p; for ((p = 0; p < nj; p++)) ; do
		$bin/$model "$@" $io --idPart=$p --logdir=$W/logs/${model}_$(printf %04d $p)
	done
	$bin/merge_parts_$model "$@" $io --idPart=0 --logdir=$W/logs/postProcess
}
_cudimot_fit () {
	local S="$1" bin="$2/bin" model="$3" nj="$4" gd="$5"; shift 5
	local O=$S.$model G1=$S.$model/GridSearch G2=$S.$model/FitFractions
	local opts="--bi=1000 --nj=1250 --se=25 --data=$S/data --maskfile=$O/nodif_brain_mask --forcedir --CFP=$O/CFP --FixP=$O/FixP"
	local th=$O/Dtifit/dtifit_V1_th.nii.gz ph=$O/Dtifit/dtifit_V1_ph.nii.gz
	rm -rf $O; mkdir -p $O/Dtifit $O/logs $G1 $G2
	cp $S/nodif_brain_mask.nii.gz $S/bvecs $S/bvals $O/
	printf '%s\n%s\n' $O/bvecs $O/bvals > $O/CFP
	# S0, fixed during the fit: the mean of the volumes with b <= 50
	local vols=$(awk '{for (i = 1; i <= NF; i++) if ($i + 0 <= 50) printf "%s%d", (n++ ? "," : ""), i - 1}' $S/bvals)
	[ -n "$vols" ] || { echo "ERROR: no b <= 50 volume in $S/bvals"; return 1; }
	${FSLDIR}/bin/fslselectvols -i $S/data -o $O/S0 --vols=$vols -m
	echo $O/S0 > $O/FixP
	# initial orientation: DTI of the same data, with the gradient nonlinearity
	# per voxel and the nominal b-table when there is one
	local bv=$S/bvals bc=$S/bvecs gnl=""
	if [ -n "$gd" ] && [ -e "$gd" ] ; then
		gnl="--gradnonlin=$gd"
		[ -e $S/bvals_nominal ] && bv=$S/bvals_nominal
		[ -e $S/bvecs_nominal ] && bc=$S/bvecs_nominal
	fi
	${FSLDIR}/bin/dtifit -k $S/data -m $S/nodif_brain_mask -r $bc -b $bv -o $O/Dtifit/dtifit --save_tensor $gnl
	# cart2spherical exits 1 even when it has written its output: judge it by the files
	$bin/cart2spherical $O/Dtifit/dtifit_V1 $O/Dtifit/dtifit_V1 || true
	[ -s $th ] && [ -s $ph ] || { echo "ERROR: cart2spherical wrote no $th / $ph"; return 1; }
	# 1 grid search over fiso, fintra, kappa with the orientation held at DTI's
	printf '\n\n\n%s\n%s\n' $th $ph > $G1/InitializationParameters
	printf '%s\n' "search[0]=(0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0)" \
		"search[1]=(0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0)" "search[2]=(1,2,3,4,5,6,7,8)" > $G1/GridSearch
	_cudimot_step $bin $model $nj $G1 $opts --gridSearch=$G1/GridSearch --no_LevMar --init_params=$G1/InitializationParameters --fixed=3,4
	# 2 the two fractions, kappa and orientation held
	printf '%s\n' $G1/Param_0_samples $G1/Param_1_samples $G1/Param_2_samples $th $ph > $G2/InitializationParameters
	_cudimot_step $bin $model $nj $G2 $opts --init_params=$G2/InitializationParameters --fixed=2,3,4
	# 3 every parameter (MCMC with Rician noise unless cudimot-mle)
	printf '%s\n' $G2/Param_0_samples $G2/Param_1_samples $G1/Param_2_samples $th $ph > $O/InitializationParameters
	_cudimot_step $bin $model $nj $O $opts --init_params=$O/InitializationParameters "$@"
	# what cuDIMOT's <model>_finish.sh does: samples renamed, means, dyads, OD
	local i n=(fiso fintra kappa th ph)
	for i in 0 1 2 3 4 ; do mv $O/Param_${i}_samples.nii.gz $O/${n[$i]}_samples.nii.gz ; done
	for i in fiso fintra kappa ; do ${FSLDIR}/bin/fslmaths $O/${i}_samples -Tmean $O/mean_$i ; done
	${FSLDIR}/bin/make_dyadic_vectors $O/th_samples $O/ph_samples $O/nodif_brain_mask.nii.gz $O/dyads1
	${FSLDIR}/bin/fslmaths $O/mean_kappa -recip -atan -mul 0.636619772367581 $O/OD
}

_cudimot () {
	local dir="$1" cudir="$2" model="$3"

	# cuDIMOT may have been built against another FSL than the one in FSLDIR:
	# CUDIMOT_FSLDIR selects it for the fit only
	local DEFAULT_FSLDIR=$FSLDIR
	[ -n "${CUDIMOT_FSLDIR:-}" ] && export FSLDIR=$CUDIMOT_FSLDIR
	# cuDIMOT's binaries resolve FSL's libraries by RPATH; anything else on
	# LD_LIBRARY_PATH ahead of them (another FSL, CUDA) shadows the wrong ones and
	# cart2spherical dies with "undefined symbol" (seen 2026-09-16)
	local DEFAULT_LD_LIBRARY_PATH=$LD_LIBRARY_PATH
	if [ ! -z "${CUDIMOTLD_LIBRARY_PATH}" ] ; then
		export LD_LIBRARY_PATH=${CUDIMOTLD_LIBRARY_PATH}
	else
		export LD_LIBRARY_PATH=
	fi
	export CUDIMOT=${cudir}
	log_Msg "CUDIMOT: ${cudir}  model: ${model}  dir: ${dir}"
	if [ `imtest ${dir}.${model}/OD.nii.gz` -eq 1 ] ; then
		imrm ${dir}.${model}/OD.nii.gz
	fi
	local last=""
	[ "$Fitter" = cudimot-mcmc ] && last="--runMCMC --rician"
	local need; for need in cart2spherical split_parts_${model} ${model} merge_parts_${model} ; do
		[ -x ${cudir}/bin/$need ] || { echo "ERROR: ${cudir}/bin/$need not found (see README, cuDIMOT model)"; exit 1; }
	done
	log_Msg "CUDIMOT fit: ${model} in ${cudir}, ${NJOB} parts, final step: ${last:-Levenberg-Marquardt} --BIC_AIC"
	LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:${FSLDIR}/lib" _cudimot_fit ${dir} ${cudir} ${model} $NJOB "${4:-}" $last --BIC_AIC > ${dir}.pipeline.log 2>&1 \
		|| { echo "ERROR: cuDIMOT failed in ${dir}; see ${dir}.pipeline.log"; exit 1; }
	[ `imtest ${dir}.${model}/OD.nii.gz` -eq 1 ] || { echo "ERROR: cuDIMOT left no OD in ${dir}.${model}; see ${dir}.pipeline.log"; exit 1; }
	export LD_LIBRARY_PATH=$DEFAULT_LD_LIBRARY_PATH
	export FSLDIR=$DEFAULT_FSLDIR
}

#########################################################
# The fit, in three stages (2026-09-19):
#   prep : the fit directories (<dir>: data, nodif_brain_mask, bvals, bvecs) and,
#          with --gradnonlin, the bins of the effective b-scale under
#          <dir>.gnl_bins/bin_XX (or --nchunk equal chunks for the MATLAB fitters)
#   fit  : every bin with the chosen fitter, into <bin>.<FitTag>
#   map  : the bins merged into <dir>.<FitTag> (fslmaths -add over their
#          disjoint masks), then unpacked / copied for the surface mapping
# cuDIMOT takes ONE b-table per dataset (bvals/bvecs are the model's common
# fixed parameters), hence the bins: each holds the data (a link), a mask of its
# own voxels, bvals scaled by |M g|^2 and bvecs rotated to M g with M = I + the
# bin's mean L, and the nominal b-table for dtifit --gradnonlin (exact per
# voxel). Bins under 200 voxels are merged into their nearest neighbour: cuDIMOT
# splits a dataset into NJOBS parts and a part with no voxel leaves no output.
#########################################################

# _prep_bins <dir> <grad_dev on this grid, or "">
_prep_bins () {
	local dir="$1" gdev="$2" G=${1}.gnl_bins
	rm -rf $G; mkdir -p $G
	python3 - "$dir" "$gdev" "$G" "$GnlBinWidth" "$NChunk" "$Fitter" <<'PY'
import sys, os, numpy as np, nibabel as nib
D, GD, G, BW, NCH, FIT = sys.argv[1], sys.argv[2], sys.argv[3], float(sys.argv[4]), int(sys.argv[5]), sys.argv[6]
mask = nib.load(D + "/nodif_brain_mask.nii.gz"); m = np.asarray(mask.dataobj) > 0; idx = np.nonzero(m); nvox = int(m.sum())
bval = np.loadtxt(D + "/bvals"); bvec = np.loadtxt(D + "/bvecs"); g = bvec.T; dw = bval > 100
def link(src, dst):
    if not os.path.lexists(dst): os.symlink(os.path.abspath(src), dst)
if GD:
    L = np.asanyarray(nib.load(GD).dataobj, dtype=np.float32)
    if L.shape[:3] != m.shape: raise SystemExit("grad_dev grid %s does not match the data grid %s" % (L.shape[:3], m.shape))
    # The nine components are written by calc_grad_perc_dev as grad_dev_x, _y, _z
    # merged in that order, so reading them row-major gives the Jacobian
    # d(w_i)/d(r_j). FSL reads the same nine as L = [[d1,d4,d7],[d2,d5,d8],
    # [d3,d6,d9]] -- the transpose -- and left-multiplies the gradient direction.
    # Verified against dtifit --gradnonlin on synthetic data built from a known
    # tensor: the column-major reading recovers the tensor exactly, the row-major
    # one leaves the first eigenvector 10.8 deg off, worse than not correcting at
    # all (4.4 deg). On HCP-YA the two readings differ by up to 6.8 deg in the
    # effective directions, so the swapaxes is not cosmetic.
    Lm = L[idx].reshape(-1, 3, 3).swapaxes(-1, -2); Mm = np.eye(3, dtype=np.float32)[None] + Lm
    scale = ((np.einsum("vij,nj->vni", Mm, g[dw]) ** 2).sum(2)).mean(1)
    edges = np.arange(np.floor(scale.min() / BW) * BW, scale.max() + BW, BW); nb = len(edges) - 1
    b = np.clip(((scale - edges[0]) / BW).astype(int), 0, nb - 1); counts = np.bincount(b, minlength=nb); MINB = 200
    for k in np.argsort(counts):
        if counts[k] == 0 or counts[k] >= MINB: continue
        others = [j for j in range(nb) if j != k and counts[j] >= MINB]
        if not others: break
        j = min(others, key=lambda j: abs(j - k)); b[b == k] = j; counts[j] += counts[k]; counts[k] = 0
    sc = np.zeros(m.shape, np.float32); sc[idx] = scale; nib.save(nib.Nifti1Image(sc, mask.affine, mask.header), G + "/bscale.nii.gz")
    kind = "gradient-nonlinearity bins (b-scale %.3f-%.3f, width %.3f)" % (scale.min(), scale.max(), BW)
elif FIT.startswith("matlab") and NCH > 1:
    nb = min(NCH, nvox); b = (np.arange(nvox) * nb) // nvox; Mm = None; kind = "%d chunks for MATLAB" % nb
else:
    nb = 1; b = np.zeros(nvox, int); Mm = None; kind = "one bin (no gradient nonlinearity)"
open(G + "/bins.txt", "w").write("")
for k in range(nb):
    sel = b == k; n = int(sel.sum())
    if n == 0: continue
    d = "%s/bin_%02d" % (G, k); os.makedirs(d, exist_ok=True)
    mk = np.zeros(m.shape, np.uint8); mk[tuple(a[sel] for a in idx)] = 1
    nib.save(nib.Nifti1Image(mk, mask.affine, mask.header), d + "/nodif_brain_mask.nii.gz")
    link(D + "/data.nii.gz", d + "/data.nii.gz")
    if Mm is not None:
        Mb = Mm[sel].mean(0); gb = (Mb @ g.T).T; bb = bval * (gb ** 2).sum(1); gb = gb / np.maximum(np.linalg.norm(gb, axis=1, keepdims=True), 1e-9); gb[~dw] = 0
        np.savetxt(d + "/bvals", np.round(bb)[None].astype(int), fmt="%d"); np.savetxt(d + "/bvecs", gb.T, fmt="%.6f")
        link(D + "/bvals", d + "/bvals_nominal"); link(D + "/bvecs", d + "/bvecs_nominal")
        open(G + "/bins.txt", "a").write("bin_%02d %d %.4f\n" % (k, n, scale[sel].mean()))
    else:
        link(D + "/bvals", d + "/bvals"); link(D + "/bvecs", d + "/bvecs")
        open(G + "/bins.txt", "a").write("bin_%02d %d 1.0000\n" % (k, n))
print("prep %s: %d voxels/vertices in %s" % (os.path.basename(D), nvox, kind))
PY
	[ $? -eq 0 ] || { echo "ERROR: binning failed for ${dir}"; exit 1; }
}

# _fit_bins <dir> <cudimot dir> <cudimot model> <grad_dev on this grid, or "">
#   every bin of <dir> with the chosen fitter; the bins' outputs are <bin>.<FitTag>
_fit_bins () {
	local dir="$1" cudir="$2" model="$3" gdev="$4" G=${1}.gnl_bins
	[ -s $G/bins.txt ] || { echo "ERROR: no bins under $G; run --runmode=prep first"; exit 1; }
	local d
	case "$Fitter" in
	cudimot-*)
		for d in $G/bin_* ; do
			[ `imtest $d.${FitTag}/OD.nii.gz` -eq 1 ] && continue
			log_Msg "fit $(basename $d): $(awk -v b=$(basename $d) '$1==b{print $2" voxels, b-scale "$3}' $G/bins.txt)"
			_cudimot $d ${cudir} ${model} "${gdev}"
			[ "$FitTag" = "$model" ] || { rm -rf $d.${FitTag}; mv $d.${model} $d.${FitTag}; }
		done ;;
	matlab*)
		local ric=0; [ "$Fitter" = matlab-rician ] && ric=1
		local tasks=$G/matlab_tasks.txt; : > $tasks; mkdir -p $G/logs
		for d in $G/bin_* ; do
			[ -e $d.${FitTag}/DONE ] && continue
			echo "${NODDIHCP}/scripts/matlab_noddi_bin.sh $d $d.${FitTag} $ric $RNODDIDIR" >> $tasks
		done
		local n=$(wc -l < $tasks)
		if [ $n -gt 0 ] ; then
			export NODDIDIR NIFTIMATLIB MATLABBIN
			local jid=$(${FSLDIR}/bin/fsl_sub ${NSM_MATLAB_QUEUE:+-q $NSM_MATLAB_QUEUE} -l $G/logs -N noddi_${Subject} -t $tasks)
			log_Msg "MATLAB NODDI: $n bins submitted with fsl_sub (job $jid); waiting"
			while : ; do
				local done=0 fail=0
				for d in $G/bin_* ; do [ -e $d.${FitTag}/DONE ] && done=$((done+1)); [ -e $d.${FitTag}/FAIL ] && fail=$((fail+1)); done
				[ $fail -gt 0 ] && { echo "ERROR: $fail MATLAB bins failed under $G (see <bin>.${FitTag}/FAIL and matlab.log)"; exit 1; }
				[ $done -ge $(ls -d $G/bin_* | wc -l) ] && break
				sleep 60
			done
		fi ;;
	esac
}

# _merge_bins <dir>: <dir>.<FitTag> from the bins (their masks are disjoint, so the maps add)
_merge_bins () {
	local dir="$1" G=${1}.gnl_bins out=${1}.${FitTag}
	[ -s $G/bins.txt ] || { echo "ERROR: no bins under $G"; exit 1; }
	rm -rf $out; mkdir -p $out
	local f d
	for f in mean_fintra mean_fiso mean_kappa OD dyads1 dyads1_dispersion AIC BIC S0 fmin ; do
		local parts=""
		for d in $G/bin_* ; do [ `imtest $d.${FitTag}/$f` -eq 1 ] && parts="$parts $d.${FitTag}/$f" ; done
		[ -n "$parts" ] || continue
		set -- $parts; ${FSLDIR}/bin/imcp $1 $out/$f; shift
		local f2; for f2 in "$@" ; do ${FSLDIR}/bin/fslmaths $out/$f -add $f2 $out/$f ; done
	done
	[ -e $G/bscale.nii.gz ] && cp $G/bscale.nii.gz $out/gnl_bscale.nii.gz
	echo "$Fitter" > $out/fitter.txt
	log_Msg "merged $(ls -d $G/bin_* | wc -l) bins into $out"
}

#########################################################
# The whole-brain volume fit of the obsolete volume route (--mapping=volume)
#########################################################
Prep_WholeBrain () {
	log_Msg "prep: whole-brain volume fit directory"
	_prep_bins $DWIT1wFolder "$GradNonlin"
}
Fit_WholeBrain () {
	log_Msg "fit: whole brain ($Fitter)"
	_fit_bins $DWIT1wFolder ${CUDIMOTDIR} ${CUDIMOTMODELNAME} "$GradNonlin"
}
Map_WholeBrain () {
	_merge_bins $DWIT1wFolder
	local F=${DWIT1wFolder}.${FitTag}
	fslmaths $F/mean_fintra $DWIT1wFolder/noddi_ficvf${FitSuffix}
	fslmaths $F/mean_fiso $DWIT1wFolder/noddi_fiso${FitSuffix}
	fslmaths $F/OD $DWIT1wFolder/noddi_odi${FitSuffix}
	fslmaths $F/dyads1 $DWIT1wFolder/noddi_dir${FitSuffix}
	fslmaths $F/mean_kappa $DWIT1wFolder/noddi_kappa${FitSuffix}
	[ `imtest $F/AIC` = 1 ] && fslmaths $F/AIC $DWIT1wFolder/noddi_AIC${FitSuffix}
	[ `imtest $F/BIC` = 1 ] && fslmaths $F/BIC $DWIT1wFolder/noddi_BIC${FitSuffix}
}


#########################################################
# The cortex: the diffusion signal sampled onto a layer surface between white
# and pial (the sampling depth, 2026-09-30), kept as CIFTI dense scalars on the
# native mesh, handed to the fitter as the NIFTI wb_command -cifti-convert makes
# of them, fitted per vertex, and converted back with the same CIFTI as template
#########################################################

# _make_layer: the layer surface of each hemisphere, $(LayerSurf <H>), on the
# native mesh in T1w/Native (the geometry the shared script samples on: T1w
# space, the same space as data.nii.gz, so no registration). Built once and
# reused by every run at the same depth; the other stages only read it.
#   frac f : w + f (p - w), as wb_command -surface-average with weights 1-f, f
#            (f = 0.5 is the native midthickness, the shared script's surface)
#   abs  d : a = min(d / t_perp, cap), t_perp = (p - w).n_white, and w + a (p - w):
#            d mm inside white measured along the white-surface normal, the vertex
#            moved along p - w, never beyond cap x the local thickness (as
#            deep37/make_layer.py; Workbench's vertex normal is the normalised sum
#            of the unit triangle normals, deep37's the area-weighted sum, which
#            moves the layer by a median 0.7 um, p99 0.06 mm on HCP-YA 103818)
#   vol  f : wb_command -surface-cortex-layer at volume fraction f from white
_make_layer () {
	local H w p out tmp T
	for H in L R ; do
		out=$(LayerSurf $H)
		[ -s $out ] && continue
		w=$T1wNativeFolder/$Subject.$H.white.native.surf.gii
		p=$T1wNativeFolder/$Subject.$H.pial.native.surf.gii
		[ -e $w ] && [ -e $p ] || { echo "ERROR: $w or $p not found"; exit 1; }
		tmp=${out%.surf.gii}.tmp$$.surf.gii
		case $DepthMode in
		frac)
			${CARET7DIR}/wb_command -surface-average $tmp -surf $w -weight `awk -v f=$DepthValue 'BEGIN{printf "%.12g", 1-f}'` -surf $p -weight $DepthValue ;;
		vol)
			${CARET7DIR}/wb_command -surface-cortex-layer $w $p $DepthValue $tmp ;;
		abs)
			T=$SurfFitFolder/_layer.$H.$$; mkdir -p $T
			${CARET7DIR}/wb_command -surface-coordinates-to-metric $w $T/w.func.gii
			${CARET7DIR}/wb_command -surface-coordinates-to-metric $p $T/p.func.gii
			${CARET7DIR}/wb_command -surface-normals $w $T/n.func.gii
			${CARET7DIR}/wb_command -metric-math "min($DepthValue / max((px - wx) * nx + (py - wy) * ny + (pz - wz) * nz, 0.000001), $DepthCap)" $T/a.func.gii \
				-var px $T/p.func.gii -column 1 -var py $T/p.func.gii -column 2 -var pz $T/p.func.gii -column 3 \
				-var wx $T/w.func.gii -column 1 -var wy $T/w.func.gii -column 2 -var wz $T/w.func.gii -column 3 \
				-var nx $T/n.func.gii -column 1 -var ny $T/n.func.gii -column 2 -var nz $T/n.func.gii -column 3 > /dev/null
			local c
			for c in 1 2 3 ; do
				${CARET7DIR}/wb_command -metric-math "w + a * (p - w)" $T/c$c.func.gii -var w $T/w.func.gii -column $c -var p $T/p.func.gii -column $c -var a $T/a.func.gii > /dev/null
			done
			${CARET7DIR}/wb_command -metric-merge $T/xyz.func.gii -metric $T/c1.func.gii -metric $T/c2.func.gii -metric $T/c3.func.gii
			${CARET7DIR}/wb_command -surface-set-coordinates $w $T/xyz.func.gii $tmp
			log_Msg "layer $H: `${CARET7DIR}/wb_command -metric-stats $T/a.func.gii -reduce MEDIAN` median fraction of the thickness; `${CARET7DIR}/wb_command -metric-math "a >= $DepthCap - 0.000001" $T/cap.func.gii -var a $T/a.func.gii > /dev/null; ${CARET7DIR}/wb_command -metric-stats $T/cap.func.gii -reduce MEAN` of the vertices at the cap"
			rm -rf $T ;;
		esac
		mv -f $tmp $out
		log_Msg "layer surface: $out ($DepthMode $DepthValue${DepthCap:+ cap $DepthCap})"
	done
}

# _volfile <path with or without .nii.gz>: the file itself
_volfile () { local f; for f in "$1" "$1.nii.gz" "$1.nii" ; do [ -f "$f" ] && { echo "$f"; return 0; }; done; echo "ERROR: $1 not found" >&2; return 1; }

# _sample_dscalar <volume> <out .dscalar.nii>: one cubic sample per vertex of
# every volume, on the layer surface of both hemispheres, as a CIFTI dense
# scalar on the native mesh (one map per volume; the vertices of the native
# roi.native.shape.gii, i.e. without the medial wall). A .nii.gz is gunzipped
# once to $TmpWork first: wb_command reading HCP-YA's 1.2 GB data.nii.gz directly
# took ~40 min per hemisphere (deep37, 2026-09-28).
_sample_dscalar () {
	local vol out raw H
	vol=`_volfile "$1"`; out="$2"
	[ -s "$out" ] && { log_Msg "exists: $out"; return 0; }
	mkdir -p $TmpWork
	case "$vol" in
	*.nii.gz)
		raw=$TmpWork/_nsm_${Subject}_${DepthTag}_$$_`basename $vol .nii.gz`.nii
		trap "rm -f '$raw'" EXIT
		gunzip -c "$vol" > "$raw" ;;
	*)	raw="$vol" ;;
	esac
	for H in L R ; do
		${CARET7DIR}/wb_command -volume-to-surface-mapping "$raw" $(LayerSurf $H) ${out%.dscalar.nii}._$H.func.gii -cubic
	done
	[ "$raw" = "$vol" ] || { rm -f "$raw"; trap - EXIT; }
	${CARET7DIR}/wb_command -cifti-create-dense-scalar ${out%.dscalar.nii}._tmp.dscalar.nii \
		-left-metric ${out%.dscalar.nii}._L.func.gii -roi-left "$AtlasSpaceNativeFolder"/"$Subject".L.roi.native.shape.gii \
		-right-metric ${out%.dscalar.nii}._R.func.gii -roi-right "$AtlasSpaceNativeFolder"/"$Subject".R.roi.native.shape.gii
	mv -f ${out%.dscalar.nii}._tmp.dscalar.nii "$out"
	rm -f ${out%.dscalar.nii}._?.func.gii
	log_Msg "sampled `basename $vol` onto the $DepthTag layer: $out"
}

Prep_Surface () {
log_Msg "prep: surface (per-vertex) fit directory, sampling depth $DepthTag"
mkdir -p $SurfFitFolder $MapFolder
local W=$SurfFitFolder T=$DepthTag
_make_layer
# ONE cubic sample on the layer surface. Averaging the ribbon first, as this
# did until 2026-09-27, mixes signal from orientations that no single NODDI
# model fits, and the misfit returns as bias and as variance both: on HCP-YA
# 103818 the ribbon average sat +0.027 in ODI and -0.007 in NDI against the
# volume fit mapped the same way, where the mid sample sits +0.004 / +0.003,
# and it was 16 % noisier as well. Nothing gates the sampling: the b0 CoV
# mask it used to pass as -volume-roi is an option single-point mapping does
# not take, and it was measured to flag 97 % of what -S already flags anyway.
# noddi_gmpv records what each sample is made of, for weighting downstream.
# The samples are kept as CIFTI (2026-09-30), so they can be read, masked and
# resampled with Workbench like any other per-vertex data.
_sample_dscalar $DWIT1wFolder/data.nii.gz $W/data_${T}.native.dscalar.nii
[ -n "$GradNonlin" ] && _sample_dscalar "$GradNonlin" $W/grad_dev_${T}.native.dscalar.nii
_confmaps
# The fit mask, as the shared script's packing made it: the native roi (the
# CIFTI's vertices) where every sample is finite and the first volume positive
# (the sum over volumes is NaN if any sample is NaN or infinite)
${CARET7DIR}/wb_command -cifti-reduce $W/data_${T}.native.dscalar.nii SUM $W/_sum.dscalar.nii
${CARET7DIR}/wb_command -cifti-math '(x0 > 0) * ((s - s) == 0)' $W/mask_${T}.native.dscalar.nii -fixnan 0 \
	-var x0 $W/data_${T}.native.dscalar.nii -select 1 1 -var s $W/_sum.dscalar.nii > /dev/null
rm -f $W/_sum.dscalar.nii
# The NIFTI the fitter reads: wb_command -cifti-convert -to-nifti lays the CIFTI
# rows (the grayordinates: left cortex, then right) along the spatial dimensions
# in order (i fastest), zero-padded. -smaller-dims makes the grid near-cubic
# (HCP-YA 103818: 232838 vertices -> 62 x 62 x 61), far inside NIFTI-1's 32767,
# and, unlike the default or -smaller-file (32767 x 8 x 1, 29105 x 8 x 1), leaves
# no singleton spatial dimension: FSL drops a trailing one when it writes a 3-D
# image, and the mask (rewritten by fslmaths below as the 3-D image the binning
# and cuDIMOT expect; wb writes a one-map CIFTI as X x Y x Z x 1) would then no
# longer match the grid of grad_dev in _prep_bins. The mask stays float32: the
# binning saves bscale.nii.gz with the mask's header, and a uint8 mask (as the
# shared script's packing wrote) truncates that map to 0/1. The fitters see the
# vertices only through the mask (cuDIMOT's split_parts, dtifit, the bins'
# masks), so the layout is theirs to ignore; the shared script's own packing
# (512 x N x 1) was the same idea by hand. Outside the mask the data and
# grad_dev are zero, as that packing made them.
${CARET7DIR}/wb_command -cifti-convert -to-nifti $W/mask_${T}.native.dscalar.nii $W/mask_${T}.nii.gz -smaller-dims
${FSLDIR}/bin/fslmaths $W/mask_${T}.nii.gz -bin $W/mask_${T}.nii.gz
local src
for src in data ${GradNonlin:+grad_dev} ; do
	${CARET7DIR}/wb_command -cifti-math 'x * m' $W/_${src}_masked.dscalar.nii -fixnan 0 \
		-var x $W/${src}_${T}.native.dscalar.nii -var m $W/mask_${T}.native.dscalar.nii -select 1 1 -repeat > /dev/null
	${CARET7DIR}/wb_command -cifti-convert -to-nifti $W/_${src}_masked.dscalar.nii $W/${src}_${T}.nii.gz -smaller-dims
	rm -f $W/_${src}_masked.dscalar.nii
done
# the names cuDIMOT (data, nodif_brain_mask) and the binning (grad_dev) look for
ln -sfn data_${T}.nii.gz $W/data.nii.gz
ln -sfn mask_${T}.nii.gz $W/nodif_brain_mask.nii.gz
if [ -n "$GradNonlin" ] ; then ln -sfn grad_dev_${T}.nii.gz $W/grad_dev.nii.gz ; else rm -f $W/grad_dev.nii.gz ; fi
cp $DWIT1wFolder/bvals $W/bvals
cp $DWIT1wFolder/bvecs $W/bvecs
{ echo "depth $DepthMode $DepthValue${DepthCap:+ cap $DepthCap}   tag $T"
  echo "layer $(LayerSurf L) $(LayerSurf R)"
  echo "CIFTI rows (native roi):" `${CARET7DIR}/wb_command -file-information $W/mask_${T}.native.dscalar.nii -no-map-info | grep -E 'Cortex(Left|Right): .* out of'`
  echo "fit mask: `${CARET7DIR}/wb_command -cifti-stats $W/mask_${T}.native.dscalar.nii -reduce COUNT_NONZERO` vertices"
  echo "nifti `${FSLDIR}/bin/fslval $W/data_${T}.nii.gz dim1` x `${FSLDIR}/bin/fslval $W/data_${T}.nii.gz dim2` x `${FSLDIR}/bin/fslval $W/data_${T}.nii.gz dim3` x `${FSLDIR}/bin/fslval $W/data_${T}.nii.gz dim4` (wb_command -cifti-convert -to-nifti -smaller-dims)"
} | tr -s ' ' > $W/depth_${T}.txt
log_Msg "`cat $W/depth_${T}.txt | tr '\n' ';'`"
local gd=""; [ -n "$GradNonlin" ] && gd=$W/grad_dev.nii.gz
_prep_bins $W "$gd"
}

Fit_Surface () {
log_Msg "fit: surface ($Fitter), sampling depth $DepthTag"
START_TIME=`date +%s`
local gd=""; [ -n "$GradNonlin" ] && gd=$SurfFitFolder/grad_dev.nii.gz
_fit_bins $SurfFitFolder ${CUDIMOTDIR} ${CUDIMOTMODELNAME} "$gd"
END_TIME=`date +%s`; SS=`echo ${END_TIME} - ${START_TIME} | bc`
log_Msg "Calculation time for the surface NODDI fit: `printf "%02d:%02d:%02d" $((SS/3600)) $(((SS%3600)/60)) $((SS%60))`"
}

Map_Surface () {
_merge_bins $SurfFitFolder
mkdir -p $MapFolder
_make_layer
# back to CIFTI: the fitted NIFTIs have the mask's grid, so the mask dscalar is
# their template. The native dscalar keeps the unsmoothed fit; the func.gii of
# each hemisphere, under the names the volume route used plus the depth tag,
# goes on to everything downstream (0 outside the roi, as the fit left it).
local F=${SurfFitFolder}.${FitTag} tmpl=$SurfFitFolder/mask_${DepthTag}.native.dscalar.nii pair src dst out
[ -e $tmpl ] || { echo "ERROR: $tmpl not found; run --runmode=prep at this depth first"; exit 1; }
for pair in mean_fintra:noddi_ficvf mean_fiso:noddi_fiso mean_kappa:noddi_kappa OD:noddi_odi AIC:noddi_AIC BIC:noddi_BIC fmin:noddi_fmin gnl_bscale:gnl_bscale ; do
	src=${pair%%:*}; dst=${pair##*:}
	[ `imtest $F/$src` = 1 ] || continue
	out=$MapFolder/"$Subject".${dst}${NTag}.native.dscalar.nii
	${CARET7DIR}/wb_command -cifti-convert -from-nifti `_volfile $F/$src` $tmpl $out
	${CARET7DIR}/wb_command -set-map-names $out -map 1 "${Subject}_${dst}${NTag}${FitSuffix}"
	${CARET7DIR}/wb_command -cifti-separate $out COLUMN \
		-metric CORTEX_LEFT $MapFolder/"$Subject".L.${dst}${NTag}.native.func.gii \
		-metric CORTEX_RIGHT $MapFolder/"$Subject".R.${dst}${NTag}.native.func.gii
done
log_Msg "fitted maps back on the native mesh: $MapFolder/$Subject.<map>${NTag}.native.dscalar.nii and .<H>.<map>${NTag}.native.func.gii"
# Smoothing AFTER the fit is free; smoothing the signal before it is not. The
# ribbon average is exactly that mistake made radially, and it costs bias and
# variance both. Here the fitted parameters are smoothed geodesically on the
# layer surface that was sampled, which lowers the noise and leaves the value alone.
for Hemisphere in L R ; do
	surf=$(LayerSurf $Hemisphere)
	# sigma = this mesh's own mean vertex spacing, so it follows the species
	# rather than a number chosen on human data
	sigma=`${CARET7DIR}/wb_command -surface-information $surf | awk '/^Spacing:/{f=1} f&&/Mean:/{print $2; exit}'`
	[ -n "$sigma" ] || { log_Msg "cannot read the vertex spacing of $surf: skipping the post-fit smoothing"; continue; }
	log_Msg "post-fit smoothing, $Hemisphere: sigma $sigma mm (mean vertex spacing of the $DepthTag layer)"
	for v in noddi_ficvf noddi_fiso noddi_kappa noddi_odi ; do
		f=$MapFolder/"$Subject"."$Hemisphere".${v}${NTag}.native.func.gii
		[ -e $f ] || continue
		${CARET7DIR}/wb_command -metric-smoothing $surf $f $sigma ${f%.func.gii}._sm.func.gii
		mv ${f%.func.gii}._sm.func.gii $f
	done
done
}

#########################################################
# noddi_gmpv: what fraction of the diffusion voxel each vertex's sample is
# cortical grey matter. The ribbon is binarised at the structural resolution and
# blurred by the box filter of one diffusion voxel (SD = voxel/sqrt(12)), then
# read at the vertex the same way the data was (on the same layer surface). It
# is a confidence map, not a mask: carry it and weight or threshold on it
# downstream rather than dropping vertices here. It matters most where the voxel
# is large against the ribbon -- on 2026-09-27 the midthickness sample was 0.999
# grey in human (1.25 mm voxel, 2.79 mm cortex) and 0.998 in macaque
# (0.90 / 1.72), but in marmoset (0.80 / 1.29) 12.5 % of vertices fall below 0.9
# and 6 % below 0.8. Deeper layers sit nearer the white boundary and lower it.
#########################################################
_confmaps () {
	local vox sd
	[ -e $T1wFolder/ribbon.nii.gz ] || { log_Msg "no ribbon.nii.gz: skipping noddi_gmpv"; return 0; }
	for Hemisphere in L R ; do
		[ -e $MapFolder/"$Subject"."$Hemisphere".noddi_gmpv${NTag}.native.func.gii ] && continue
		if [ ! -e $SurfFitFolder/gmpv.nii.gz ] ; then
			vox=`${FSLDIR}/bin/fslval $DWIT1wFolder/data.nii.gz pixdim1 | tr -d ' '`
			sd=`echo "$vox / 3.4641016" | bc -l`
			log_Msg "noddi_gmpv: diffusion voxel ${vox} mm, box-filter SD ${sd} mm"
			${FSLDIR}/bin/fslmaths $T1wFolder/ribbon.nii.gz -thr $ribbonLlabel -uthr $ribbonLlabel -bin $SurfFitFolder/_gml
			${FSLDIR}/bin/fslmaths $T1wFolder/ribbon.nii.gz -thr $ribbonRlabel -uthr $ribbonRlabel -bin $SurfFitFolder/_gmr
			${FSLDIR}/bin/fslmaths $SurfFitFolder/_gml -add $SurfFitFolder/_gmr -bin -s $sd $SurfFitFolder/gmpv
			${FSLDIR}/bin/imrm $SurfFitFolder/_gml $SurfFitFolder/_gmr
		fi
		mkdir -p $MapFolder
		${CARET7DIR}/wb_command -volume-to-surface-mapping $SurfFitFolder/gmpv.nii.gz \
			$(LayerSurf $Hemisphere) \
			$MapFolder/"$Subject"."$Hemisphere".noddi_gmpv${NTag}.native.func.gii -cubic
		${CARET7DIR}/wb_command -set-map-name $MapFolder/"$Subject"."$Hemisphere".noddi_gmpv${NTag}.native.func.gii 1 "$Subject"_"$Hemisphere"_noddi_gmpv${NTag}
	done
}

#########################################################
# The subcortical grey matter: the volume fit restricted to the FreeSurfer
# subcortical labels of wmparc (dilated one voxel), for the volume part of the CIFTI
#########################################################
Prep_Subcortical () {
log_Msg "prep: subcortical fit directory"
mkdir -p $SubFitFolder
if [ ! -e $SubFitFolder/nodif_brain_mask.nii.gz ] ; then
	flirt -in $T1wFolder/wmparc.nii.gz -ref $DWIT1wFolder/nodif_brain_mask.nii.gz -applyxfm -usesqform -interp nearestneighbour -out $SubFitFolder/wmparc_diffgrid.nii.gz
	${CARET7DIR}/wb_command -volume-label-import $SubFitFolder/wmparc_diffgrid.nii.gz $HCPPIPEDIR/global/config/FreeSurferSubcorticalLabelTableLut.txt $SubFitFolder/subcortical_labels.nii.gz -discard-others -drop-unused-labels
	fslmaths $SubFitFolder/subcortical_labels.nii.gz -bin -dilM -mas $DWIT1wFolder/nodif_brain_mask.nii.gz $SubFitFolder/nodif_brain_mask.nii.gz
fi
ln -sf $DWIT1wFolder/data.nii.gz $SubFitFolder/data.nii.gz
cp $DWIT1wFolder/bvals $SubFitFolder/bvals; cp $DWIT1wFolder/bvecs $SubFitFolder/bvecs
_prep_bins $SubFitFolder "$GradNonlin"
}
Fit_Subcortical () {
log_Msg "fit: subcortical ($Fitter)"
_fit_bins $SubFitFolder ${CUDIMOTDIR} ${CUDIMOTMODELNAME} "$GradNonlin"
}
Map_Subcortical () {
_merge_bins $SubFitFolder
local F=${SubFitFolder}.${FitTag}
fslmaths $F/mean_fintra $DWIT1wFolder/noddi_ficvf_subcortical${FitSuffix}
fslmaths $F/mean_fiso $DWIT1wFolder/noddi_fiso_subcortical${FitSuffix}
fslmaths $F/mean_kappa $DWIT1wFolder/noddi_kappa_subcortical${FitSuffix}
[ `imtest $F/AIC` = 1 ] && fslmaths $F/AIC $DWIT1wFolder/noddi_AIC_subcortical${FitSuffix}
[ `imtest $F/BIC` = 1 ] && fslmaths $F/BIC $DWIT1wFolder/noddi_BIC_subcortical${FitSuffix}
}

#########################################################
# The white matter (-w): the cuDIMOT built for parallel diffusivity 1.7
# (or the MATLAB fitters); white-matter voxels only, volume outputs
#########################################################
Prep_WM () {
log_Msg "prep: white-matter fit directory"
mkdir -p $WMFitFolder
if [ ! -e $WMFitFolder/nodif_brain_mask.nii.gz ] ; then
	flirt -in $T1wFolder/wmparc.nii.gz -ref $DWIT1wFolder/nodif_brain_mask.nii.gz -applyxfm -usesqform -interp nearestneighbour -out $WMFitFolder/wmparc_diffgrid.nii.gz
	# cerebral white matter (2, 41), the corpus callosum (251-255) and the
	# white-matter parcels of wmparc (3000-4999)
	${CARET7DIR}/wb_command -volume-math '((x==2)+(x==41)+(x>=251)*(x<=255)+(x>=3000)*(x<5000))>0' $WMFitFolder/nodif_brain_mask.nii.gz -var x $WMFitFolder/wmparc_diffgrid.nii.gz
	fslmaths $WMFitFolder/nodif_brain_mask.nii.gz -mas $DWIT1wFolder/nodif_brain_mask.nii.gz $WMFitFolder/nodif_brain_mask.nii.gz
fi
ln -sf $DWIT1wFolder/data.nii.gz $WMFitFolder/data.nii.gz
cp $DWIT1wFolder/bvals $WMFitFolder/bvals; cp $DWIT1wFolder/bvecs $WMFitFolder/bvecs
_prep_bins $WMFitFolder "$GradNonlin"
}
Fit_WM () {
log_Msg "fit: white matter ($Fitter, dPar 1.7 build ${CUDIMOTDIR_WM} for cuDIMOT)"
case "$Fitter" in cudimot-*) [ -x ${CUDIMOTDIR_WM}/bin/${CUDIMOTMODELNAME_WM} ] || { echo "ERROR: ${CUDIMOTDIR_WM}/bin/${CUDIMOTMODELNAME_WM} not found"; exit 1; } ;; esac
_fit_bins $WMFitFolder ${CUDIMOTDIR_WM} ${CUDIMOTMODELNAME_WM} "$GradNonlin"
}
Map_WM () {
_merge_bins $WMFitFolder
local F=${WMFitFolder}.${FitTag}
for pair in mean_fintra:noddi_wm_ficvf mean_fiso:noddi_wm_fiso mean_kappa:noddi_wm_kappa OD:noddi_wm_odi dyads1:noddi_wm_dir AIC:noddi_wm_AIC BIC:noddi_wm_BIC ; do
	src=${pair%%:*}; dst=${pair##*:}${FitSuffix}
	[ `imtest $F/$src` = 1 ] && fslmaths $F/$src $DWIT1wFolder/$dst
done
# and in MNI space, for the volume users
for dst in noddi_wm_ficvf noddi_wm_fiso noddi_wm_kappa noddi_wm_odi ; do
	[ `imtest $DWIT1wFolder/${dst}${FitSuffix}` = 1 ] && ${CARET7DIR}/wb_command -volume-warpfield-resample $DWIT1wFolder/${dst}${FitSuffix}.nii.gz $AtlasSpaceFolder/xfms/acpc_dc2standard.nii.gz $AtlasSpaceFolder/T1w_restore.nii.gz CUBIC $AtlasSpaceResultsDWIFolder/${dst}${FitSuffix}.nii.gz -fnirt $T1wFolder/T1w_acpc_dc_restore.nii.gz
done
}

DiffusionStats () {

log_Msg "Start: DiffusionStats"

${NODDIHCP}/scripts/dwistats $DWIT1wFolder/data.nii.gz $DWIT1wFolder/bvals $DWIT1wFolder/data $DWIT1wFolder/nodif_brain_mask.nii.gz $b0thresh

}

#########################################################
# map one volume metric onto the native surfaces, by the chosen method
#   _mapvolume <vol name in DWIT1wFolder>
#########################################################
_mapvolume () {
	local vol="$1" outname="${2:-$1}"
	if [ "$MappingMethod" = "volume" ] ; then
		${CARET7DIR}/wb_command -volume-warpfield-resample $DWIT1wFolder/${vol}.nii.gz $AtlasSpaceFolder/xfms/acpc_dc2standard.nii.gz $AtlasSpaceFolder/T1w_restore.nii.gz CUBIC $AtlasSpaceResultsDWIFolder/${vol}.nii.gz -fnirt $T1wFolder/T1w_acpc_dc_restore.nii.gz
		for Hemisphere in L R ; do
		${CARET7DIR}/wb_command -volume-to-surface-mapping $AtlasSpaceResultsDWIFolder/${vol}.nii.gz "$AtlasSpaceNativeFolder"/"$Subject"."$Hemisphere".midthickness.native.surf.gii $MapFolder/"$Subject"."$Hemisphere".${outname}${NTag}.native.func.gii -myelin-style $AtlasSpaceFolder/ribbon_"$Hemisphere".nii.gz "$AtlasSpaceNativeFolder"/"$Subject"."$Hemisphere".thickness.native.shape.gii "$NODDIMappingSigma"
		done
	else
		# the same sampling the fit used (the layer surface of the sampling
		# depth): a volume metric has to sit on the same support as the fitted
		# ones or the two cannot be read against each other on the same surface
		for Hemisphere in L R ; do
		${CARET7DIR}/wb_command -volume-to-surface-mapping $DWIT1wFolder/${vol}.nii.gz $(LayerSurf $Hemisphere) $MapFolder/"$Subject"."$Hemisphere".${outname}${NTag}.native.func.gii -cubic
		done
	fi
}

DiffusionSurfaceMapping() {

	log_Msg "Start: DiffusionSurfaceMapping ($MappingMethod)"

	# --- For BALSA preprocessed data ---
	# Split the native-space cortical thickness CIFTI into per-hemisphere metric files
	# if either the left or right hemisphere file is missing.
	LThickness="$AtlasSpaceNativeFolder"/"$Subject".L.thickness.native.shape.gii
	RThickness="$AtlasSpaceNativeFolder"/"$Subject".R.thickness.native.shape.gii
	if [ ! -e "$LThickness" ] || [ ! -e "$RThickness" ]; then
	${CARET7DIR}/wb_command -cifti-separate "$AtlasSpaceNativeFolder"/"$Subject".thickness.native.dscalar.nii COLUMN -metric CORTEX_LEFT  "$LThickness" -metric CORTEX_RIGHT "$RThickness"
	fi
	# --------------------------------------------------

	fslmaths $AtlasSpaceFolder/ribbon.nii.gz -thr $ribbonLlabel -uthr $ribbonLlabel -bin $AtlasSpaceFolder/ribbon_L.nii.gz
	fslmaths $AtlasSpaceFolder/ribbon.nii.gz -thr $ribbonRlabel -uthr $ribbonRlabel -bin $AtlasSpaceFolder/ribbon_R.nii.gz
	for GrayordinatesResolution in ${GrayordinatesResolutions[@]} ; do
	if [ ! -e "$ROIFolder"/Atlas_ROIs."$GrayordinatesResolution".nii.gz ] ; then
	cp ${GrayordinatesSpaceDIR}/Atlas_ROIs."$GrayordinatesResolution".nii.gz "$ROIFolder"/Atlas_ROIs."$GrayordinatesResolution".nii.gz
	fi
	if [ ! -e $AtlasSpaceFolder/T1w_restore."$GrayordinatesResolution".nii.gz ] ; then
	flirt -in $AtlasSpaceFolder/T1w_restore.nii.gz -applyisoxfm "$GrayordinatesResolution" -ref $AtlasSpaceFolder/ROIs/Atlas_ROIs."$GrayordinatesResolution".nii.gz -o $AtlasSpaceFolder/T1w_restore."$GrayordinatesResolution".nii.gz -interp sinc
	fi
	done

	mkdir -p $MapFolder

	# SNR surface mapping (a volume metric: mapped by the chosen method)
	if [ `imtest $DWIT1wFolder/data_snr.nii.gz` = 1 ] ; then
	_mapvolume data_snr
	for Hemisphere in L R ; do
	${CARET7DIR}/wb_command -metric-mask $MapFolder/"$Subject"."$Hemisphere".data_snr${NTag}.native.func.gii "$AtlasSpaceNativeFolder"/"$Subject"."$Hemisphere".roi.native.shape.gii $MapFolder/"$Subject"."$Hemisphere".data_snr${NTag}.native.func.gii
	${CARET7DIR}/wb_command -metric-math 'x>'$SNRthr'' $MapFolder/"$Subject"."$Hemisphere".goodvertex${NTag}.native.func.gii -var x  $MapFolder/"$Subject"."$Hemisphere".data_snr${NTag}.native.func.gii
	done
	else
	for Hemisphere in L R ; do
	cp "$AtlasSpaceNativeFolder"/"$Subject"."$Hemisphere".roi.native.shape.gii $MapFolder/"$Subject"."$Hemisphere".goodvertex${NTag}.native.func.gii
	done
	fi

	# NODDI validity (-S). Where the signal has dropped out (orbitofrontal,
	# inferior temporal: the sinuses and the mastoid air) the fit puts most of the
	# signal into the isotropic compartment and NDI, defined within the remainder,
	# inflates. A b0 coefficient-of-variation mask cannot catch this: dropout is
	# uniformly dark, not variable. The flag is the vertex's b0 SNR (mean / sd of the b0 volumes
	# mapped onto the vertex; the volume SNR map when the fit was not on the
	# surface) below SNRfrac x the cortical median of the subject -- relative, so
	# that it follows each dataset's own noise level. It is the only criterion:
	# an fiso threshold and the b0 CoV mask were both measured to be the same
	# thing seen from further away (97 % and 60 % overlap), and noddi_gmpv, which
	# is nearly independent of it, is left as a map to weight on downstream. On 2026-09-17 (HCP-YA 122317 and Rlife H17013003) 0.4 x median
	# flagged 3-7 % of the cortex, two thirds of it on the ventral surfaces, and
	# removed most of the abnormal NDI there; fiso > 0.3 flagged 6-19 % scattered
	# over the whole cortex and was made optional.
	if [ "$MappingMethod" = "volume" ] && [ `imtest $DWIT1wFolder/noddi_fiso.nii.gz` = 1 ] && [ ! -e $MapFolder/"$Subject".L.noddi_fiso${NTag}.native.func.gii ] ; then
		_mapvolume noddi_fiso
	fi
	python3 - $MapFolder/"$Subject" "$AtlasSpaceNativeFolder"/"$Subject" "$SurfFitFolder" "$SNRfrac" "$NTag" <<'PYV'
import sys, os, numpy as np, nibabel as nib
pre, npre, sff, sfrac, tag = sys.argv[1], sys.argv[2], sys.argv[3], float(sys.argv[4]), sys.argv[5]
g = lambda p: nib.load(p).darrays[0].data.astype(np.float32)
roi, good, snr = {}, {}, {}
# the b0 samples of the fit, per hemisphere, from the CIFTI the prep stage kept
# (data_<depth tag>.native.dscalar.nii; NaN off the native roi)
b0 = {}
dp = "%s/data%s.native.dscalar.nii" % (sff, tag)
if sff and tag and os.path.exists(dp) and os.path.exists(sff + "/bvals"):
    ib0 = np.nonzero(np.loadtxt(sff + "/bvals") < 50)[0]
    img = nib.load(dp); ax = img.header.get_axis(1)
    X = np.asarray(img.dataobj, dtype=np.float32)[ib0]
    for name, slc, bm in ax.iter_structures():
        H = {"CIFTI_STRUCTURE_CORTEX_LEFT": "L", "CIFTI_STRUCTURE_CORTEX_RIGHT": "R"}.get(name)
        if H is None: continue
        A = np.full((ax.nvertices[name], len(ib0)), np.nan, np.float32); A[bm.vertex] = X[:, slc].T; b0[H] = A
for H in "LR":
    roi[H] = g("%s.%s.roi.native.shape.gii" % (npre, H)) > 0
    try: good[H] = g("%s.%s.goodvertex%s.native.func.gii" % (pre, H, tag)) > 0
    except Exception: good[H] = roi[H].copy()
    if H in b0:
        snr[H] = b0[H].mean(1) / np.maximum(b0[H].std(1, ddof=1), 1e-6); src = "b0 volumes mapped onto the vertex"
    else:
        try: snr[H] = g("%s.%s.data_snr%s.native.func.gii" % (pre, H, tag)); src = "volume SNR map"
        except Exception: snr[H] = None; src = "none"
med = None
if sfrac > 0 and all(snr[H] is not None for H in "LR"):
    allsnr = np.concatenate([snr[H][roi[H]] for H in "LR"]); med = float(np.median(allsnr[np.isfinite(allsnr) & (allsnr > 0)]))
for H in "LR":
    valid = roi[H] & good[H]
    if med is not None: valid &= np.isfinite(snr[H]) & (snr[H] >= sfrac * med)
    da = nib.gifti.GiftiDataArray(valid.astype(np.float32), intent="NIFTI_INTENT_NORMAL")
    nib.save(nib.gifti.GiftiImage(darrays=[da]), "%s.%s.noddi_valid%s.native.func.gii" % (pre, H, tag))
    # the SNR -S thresholds, written out whatever -S is set to, so the same
    # judgement can be made downstream instead
    if snr[H] is not None:
        sa = nib.gifti.GiftiDataArray(np.nan_to_num(snr[H]).astype(np.float32), intent="NIFTI_INTENT_NORMAL")
        nib.save(nib.gifti.GiftiImage(darrays=[sa]), "%s.%s.noddi_snr%s.native.func.gii" % (pre, H, tag))
    print("%s %s: %d of %d cortex vertices valid (%.1f %% flagged; SNR < %.2f x median %s [%s])" % (os.path.basename(pre), H, valid.sum(), roi[H].sum(), 100.0 * (roi[H] & ~valid).sum() / max(roi[H].sum(), 1), sfrac, ("%.1f" % med) if med is not None else "n/a", src))
PYV
	for Hemisphere in L R ; do
	${CARET7DIR}/wb_command -set-map-name $MapFolder/"$Subject"."$Hemisphere".noddi_valid${NTag}.native.func.gii 1 "$Subject"_"$Hemisphere"_noddi_valid${NTag}
	[ -e $MapFolder/"$Subject"."$Hemisphere".noddi_snr${NTag}.native.func.gii ] && \
		${CARET7DIR}/wb_command -set-map-name $MapFolder/"$Subject"."$Hemisphere".noddi_snr${NTag}.native.func.gii 1 "$Subject"_"$Hemisphere"_noddi_snr${NTag}
	done

	# Diffusion metrics on the native surface. Volume metrics (DTI; and every
	# NODDI metric on the volume route) are mapped; on the surface route the
	# NODDI metrics are already per-vertex, from the fit. The NODDI metrics are
	# masked by the validity flag and the flagged vertices filled from their
	# neighbours (nearest value, no distance limit: the dropout patches are
	# 1-2 cm across), as HCP's fMRISurface fills its goodvoxels holes; with -N
	# they are set to NaN on the 32k mesh instead. noddi_valid says which.
	for vol in dti_FA dti_MD noddi_kappa noddi_ficvf noddi_fiso noddi_odi noddi_AIC noddi_BIC; do
	case $vol in
		dti_*) [ `imtest $DWIT1wFolder/${vol}.nii.gz` = 1 ] || continue; _mapvolume $vol ;;
		*) if [ "$MappingMethod" = "volume" ] ; then
			[ `imtest $DWIT1wFolder/${vol}${FitSuffix}.nii.gz` = 1 ] || continue
			[ -e $MapFolder/"$Subject".L.${vol}${NTag}.native.func.gii ] || _mapvolume ${vol}${FitSuffix} ${vol}
		   else
			[ -e $MapFolder/"$Subject".L.${vol}${NTag}.native.func.gii ] || continue
		   fi ;;
	esac
	case $vol in noddi_*) maskname=noddi_valid; dist=1000 ;; *) maskname=goodvertex; dist=20 ;; esac
	for Hemisphere in L R ; do
	${CARET7DIR}/wb_command -metric-mask $MapFolder/"$Subject"."$Hemisphere".${vol}${NTag}.native.func.gii $MapFolder/"$Subject"."$Hemisphere".${maskname}${NTag}.native.func.gii $MapFolder/"$Subject"."$Hemisphere".${vol}${NTag}.native.func.gii
	${CARET7DIR}/wb_command -metric-dilate $MapFolder/"$Subject"."$Hemisphere".${vol}${NTag}.native.func.gii "$AtlasSpaceNativeFolder"/"$Subject"."$Hemisphere".midthickness.native.surf.gii $dist $MapFolder/"$Subject"."$Hemisphere".${vol}${NTag}.native.func.gii -nearest
	${CARET7DIR}/wb_command -metric-mask $MapFolder/"$Subject"."$Hemisphere".${vol}${NTag}.native.func.gii "$AtlasSpaceNativeFolder"/"$Subject"."$Hemisphere".roi.native.shape.gii $MapFolder/"$Subject"."$Hemisphere".${vol}${NTag}.native.func.gii
	${CARET7DIR}/wb_command -set-map-name $MapFolder/"$Subject"."$Hemisphere".${vol}${NTag}.native.func.gii 1 "$Subject"_"$Hemisphere"_"$vol"${NTag}
	${CARET7DIR}/wb_command -metric-palette $MapFolder/"$Subject"."$Hemisphere".${vol}${NTag}.native.func.gii MODE_AUTO_SCALE_PERCENTAGE -pos-percent 4 96 -interpolate true -palette-name videen_style -disp-pos true -disp-neg false -disp-zero false
	done
	done

	_surfaceMappingbyReg() {
		# noddi_valid first: the NODDI metrics need its 32k version to place their NaNs
		for vol in noddi_valid dti_FA dti_MD noddi_kappa noddi_ficvf noddi_fiso noddi_odi data_snr noddi_AIC noddi_BIC noddi_gmpv noddi_snr; do
		[ -e $MapFolder/"$Subject".L.${vol}${NTag}.native.func.gii ] || continue

		for Hemisphere in L R ; do

		#LowResMesh
		for LowResMesh in ${LowResMeshes[@]}; do
			DownsampleFolder=$AtlasSpaceFolder/fsaverage_LR${LowResMesh}k
			${CARET7DIR}/wb_command -metric-resample $MapFolder/"$Subject"."$Hemisphere".${vol}${NTag}.native.func.gii "$AtlasSpaceNativeFolder"/"$Subject"."$Hemisphere".sphere.${REGNAME}.native.surf.gii "$DownsampleFolder"/"$Subject"."$Hemisphere".sphere."$LowResMesh"k_fs_LR.surf.gii ADAP_BARY_AREA $AtlasSpaceResultsDWIFolder/"$Subject"."$Hemisphere".${vol}${NTag}${Reg}."$LowResMesh"k_fs_LR.func.gii -area-surfs "$AtlasSpaceNativeFolder"/"$Subject"."$Hemisphere".midthickness.native.surf.gii "$DownsampleFolder"/"$Subject"."$Hemisphere".midthickness${Reg}."$LowResMesh"k_fs_LR.surf.gii -current-roi "$AtlasSpaceNativeFolder"/"$Subject"."$Hemisphere".roi.native.shape.gii
			${CARET7DIR}/wb_command -metric-mask $AtlasSpaceResultsDWIFolder/"$Subject"."$Hemisphere".${vol}${NTag}${Reg}."$LowResMesh"k_fs_LR.func.gii "$DownsampleFolder"/"$Subject"."$Hemisphere".atlasroi."$LowResMesh"k_fs_LR.shape.gii $AtlasSpaceResultsDWIFolder/"$Subject"."$Hemisphere".${vol}${NTag}${Reg}."$LowResMesh"k_fs_LR.func.gii 
			${CARET7DIR}/wb_command -metric-smoothing "$DownsampleFolder"/"$Subject"."$Hemisphere".midthickness${Reg}."$LowResMesh"k_fs_LR.surf.gii $AtlasSpaceResultsDWIFolder/"$Subject"."$Hemisphere".${vol}${NTag}${Reg}."$LowResMesh"k_fs_LR.func.gii "$SmoothingSigma" $AtlasSpaceResultsDWIFolder/"$Subject"."$Hemisphere".${vol}${NTag}${Reg}_s"$SmoothingFWHM"."$LowResMesh"k_fs_LR.func.gii -roi "$DownsampleFolder"/"$Subject"."$Hemisphere".atlasroi."$LowResMesh"k_fs_LR.shape.gii
			# the validity mask travels unsmoothed (a resampled mask, thresholded at 0.5);
			# the NODDI metrics are set to NaN where it is 0
			case $vol in
			noddi_gmpv|noddi_snr)
				# geometry, not a fit: it must survive where the fit is invalid,
				# because explaining those vertices is the whole point of it
				;;
			noddi_valid)
				${CARET7DIR}/wb_command -metric-math '(x>0.5)' $AtlasSpaceResultsDWIFolder/"$Subject"."$Hemisphere".${vol}${NTag}${Reg}_s"$SmoothingFWHM"."$LowResMesh"k_fs_LR.func.gii -var x $AtlasSpaceResultsDWIFolder/"$Subject"."$Hemisphere".${vol}${NTag}${Reg}."$LowResMesh"k_fs_LR.func.gii > /dev/null ;;
			noddi_*)
				[ "$KeepNaN" = "YES" ] || continue
				python3 - $AtlasSpaceResultsDWIFolder/"$Subject"."$Hemisphere".${vol}${NTag}${Reg}_s"$SmoothingFWHM"."$LowResMesh"k_fs_LR.func.gii $AtlasSpaceResultsDWIFolder/"$Subject"."$Hemisphere".noddi_valid${NTag}${Reg}_s"$SmoothingFWHM"."$LowResMesh"k_fs_LR.func.gii <<'PYN'
import sys, numpy as np, nibabel as nib
mp, vp = sys.argv[1:3]
try:
    v = nib.load(vp).darrays[0].data > 0.5
except Exception:
    sys.exit(0)
im = nib.load(mp); d = im.darrays[0].data.astype(np.float32); d[~v] = np.nan; im.darrays[0].data = d; nib.save(im, mp)
PYN
				;;
			esac
		done
		done

		# the volume part: the subcortical grey matter, unless cortex only
		if [ "$CortexOnly" != "YES" ] ; then
		# which volume carries this metric: on the surface route the subcortical
		# fit (noddi_*_subcortical), on the volume route the whole-brain fit
		subvol=$vol
		if [ "$MappingMethod" = "surface" ] ; then
			case $vol in noddi_*) subvol=${vol}_subcortical${FitSuffix} ;; esac
		else
			case $vol in noddi_*) subvol=${vol}${FitSuffix} ;; esac
		fi
		if [ `imtest $DWIT1wFolder/${subvol}.nii.gz` = 1 ] ; then
		for GrayordinatesResolution in ${GrayordinatesResolutions[@]} ; do
		${CARET7DIR}/wb_command -volume-warpfield-resample $DWIT1wFolder/${subvol}.nii.gz $AtlasSpaceFolder/xfms/acpc_dc2standard.nii.gz $AtlasSpaceFolder/T1w_restore."$GrayordinatesResolution".nii.gz CUBIC $AtlasSpaceResultsDWIFolder/${vol}${NTag}."$GrayordinatesResolution".nii.gz -fnirt $T1wFolder/T1w_acpc_dc_restore.nii.gz
		if [ ! -e "$ROIFolder"/ROIs."$GrayordinatesResolution".nii.gz ] ; then
			applywarp --interp=nn -i "$AtlasSpaceFolder"/wmparc.nii.gz -r "$ROIFolder"/Atlas_ROIs."$GrayordinatesResolution".nii.gz -o "$ROIFolder"/ROIs.${GrayordinatesResolution}.nii.gz
			${CARET7DIR}/wb_command -volume-label-import "$ROIFolder"/ROIs.${GrayordinatesResolution}.nii.gz $HCPPIPEDIR/global/config/FreeSurferSubcorticalLabelTableLut.txt "$ROIFolder"/ROIs.${GrayordinatesResolution}.nii.gz -discard-others -drop-unused-labels
		fi
		${CARET7DIR}/wb_command -volume-parcel-resampling $AtlasSpaceResultsDWIFolder/${vol}${NTag}."$GrayordinatesResolution".nii.gz "$ROIFolder"/ROIs."$GrayordinatesResolution".nii.gz "$ROIFolder"/Atlas_ROIs."$GrayordinatesResolution".nii.gz $SmoothingSigma $AtlasSpaceResultsDWIFolder/${vol}${NTag}_AtlasSubcortical."$GrayordinatesResolution"_s"$SmoothingFWHM".nii.gz -fix-zeros
		done
		havevol=YES
		else
		havevol=NO
		fi
		else
		havevol=NO
		fi

		# Merge surface (and subcortical volume) to create cifti
		i=0
		for LowResMesh in ${LowResMeshes[@]}; do
		GrayordinatesResolution="${GrayordinatesResolutions[$i]}"
		DownsampleFolder=$AtlasSpaceFolder/fsaverage_LR${LowResMesh}k
		case $vol in noddi_*) fs=$FitSuffix ;; *) fs="" ;; esac
		out=$AtlasSpaceResultsDWIFolder/${vol}${OutSuffix}${fs}${Reg}."$LowResMesh"k_fs_LR.dscalar.nii
		if [ "$havevol" = "YES" ] ; then
		${CARET7DIR}/wb_command -cifti-create-dense-scalar $out -volume $AtlasSpaceResultsDWIFolder/${vol}${NTag}_AtlasSubcortical.${GrayordinatesResolution}_s"$SmoothingFWHM".nii.gz "$ROIFolder"/Atlas_ROIs.${GrayordinatesResolution}.nii.gz -left-metric $AtlasSpaceResultsDWIFolder/"$Subject".L.${vol}${NTag}${Reg}_s"$SmoothingFWHM"."$LowResMesh"k_fs_LR.func.gii -roi-left "$DownsampleFolder"/"$Subject".L.atlasroi."$LowResMesh"k_fs_LR.shape.gii -right-metric $AtlasSpaceResultsDWIFolder/"$Subject".R.${vol}${NTag}${Reg}_s"$SmoothingFWHM"."$LowResMesh"k_fs_LR.func.gii -roi-right "$DownsampleFolder"/"$Subject".R.atlasroi."$LowResMesh"k_fs_LR.shape.gii
		else
		${CARET7DIR}/wb_command -cifti-create-dense-scalar $out -left-metric $AtlasSpaceResultsDWIFolder/"$Subject".L.${vol}${NTag}${Reg}_s"$SmoothingFWHM"."$LowResMesh"k_fs_LR.func.gii -roi-left "$DownsampleFolder"/"$Subject".L.atlasroi."$LowResMesh"k_fs_LR.shape.gii -right-metric $AtlasSpaceResultsDWIFolder/"$Subject".R.${vol}${NTag}${Reg}_s"$SmoothingFWHM"."$LowResMesh"k_fs_LR.func.gii -roi-right "$DownsampleFolder"/"$Subject".R.atlasroi."$LowResMesh"k_fs_LR.shape.gii
		fi
		${CARET7DIR}/wb_command -set-map-names $out -map 1 "${Subject}_${vol}${OutSuffix}${fs}"
		${CARET7DIR}/wb_command -cifti-palette $out MODE_AUTO_SCALE_PERCENTAGE $out -pos-percent 4 96 -interpolate true -palette-name videen_style -disp-pos true -disp-neg false -disp-zero false
		i=`expr $i + 1`
		done
		done

		# Convert kappa to odi (the volume route has no odi surface of its own)
		for LowResMesh in ${LowResMeshes[@]}; do
		if [ ! -e $AtlasSpaceResultsDWIFolder/noddi_odi${OutSuffix}${FitSuffix}${Reg}."$LowResMesh"k_fs_LR.dscalar.nii ] && [ -e $AtlasSpaceResultsDWIFolder/noddi_kappa${OutSuffix}${FitSuffix}${Reg}."$LowResMesh"k_fs_LR.dscalar.nii ] ; then
		${CARET7DIR}/wb_command -cifti-math 'max(2*atan(1/kappa)/PI,0)' $AtlasSpaceResultsDWIFolder/noddi_odi${OutSuffix}${FitSuffix}${Reg}."$LowResMesh"k_fs_LR.dscalar.nii -var kappa $AtlasSpaceResultsDWIFolder/noddi_kappa${OutSuffix}${FitSuffix}${Reg}."$LowResMesh"k_fs_LR.dscalar.nii
		fi
		done

		# Remove files
		for vol in dti_FA dti_MD noddi_kappa noddi_ficvf noddi_fiso noddi_odi data_snr noddi_AIC noddi_BIC noddi_valid noddi_gmpv noddi_snr; do
		for Hemisphere in R L; do
		for LowResMesh in ${LowResMeshes[@]} ; do
			\rm -rf $AtlasSpaceResultsDWIFolder/"$Subject".${Hemisphere}.${vol}${NTag}${Reg}."$LowResMesh"k_fs_LR.func.gii
			\rm -rf $AtlasSpaceResultsDWIFolder/"$Subject".${Hemisphere}.${vol}${NTag}${Reg}_s"$SmoothingFWHM"."$LowResMesh"k_fs_LR.func.gii
		done
		done
		for GrayordinatesResolution in ${GrayordinatesResolutions[@]} ; do
		\rm -rf $AtlasSpaceResultsDWIFolder/${vol}${NTag}."$GrayordinatesResolution".nii.gz
		\rm -rf $AtlasSpaceResultsDWIFolder/${vol}${NTag}_AtlasSubcortical."$GrayordinatesResolution"_s"$SmoothingFWHM".nii.gz
		done
		done
	}

	_surfaceMappingbyReg
	
	#If REGNAME is MSMAll, perform surface mapping for both MSMAll and MSMSulc registrations
	if [ "$REGNAME" = "MSMAll" ]; then
		log_Msg "REGNAME=MSMAll: additionally performing surface mapping with MSMSulc registration"
		REGNAME="MSMSulc"
		Reg=""
		_surfaceMappingbyReg

		#restore
		REGNAME="MSMAll"
		Reg="_MSMAll"
	fi
}

#########################################################
# HippUnfold route: per-vertex NODDI on the hippocampal surfaces.
# The diffusion signal is sampled (cubic) on the midthickness-depth layer of the
# FINEST density, packed as a CIFTI over HIPPOCAMPUS_LEFT/RIGHT, fitted per vertex
# with the same cuDIMOT machinery as the cortex (via the -cifti-convert -to-nifti
# fake NIFTI), and the fitted maps are then downsampled on the surface to the
# coarser densities by interpolating in the shared unfold space. grad_dev is
# honoured exactly as for the cortex (sampled on the same layer, binned by the
# effective b-scale). hipp and dentate are fitted TOGETHER in one CIFTI that uses
# the native HIPPOCAMPUS[_DENTATE]_LEFT/RIGHT structures (requires wb_command
# >= 2.1.0, which SetUp selects for this mode).
#########################################################
_hipp_finest () { echo $HippDens | awk '{print $1}'; }
# a surface of this structure in the (T1w) HippUnfold tree: <struct> <inner|outer|midthickness> <den> <H>
_hipp_surf   () { echo "$HippUnfoldDir/$3/$Subject.$4.$1_$2.$3.surf.gii"; }
# the unfold-space midthickness (BIDS layout under hippunfold/): <struct> <den> <H>
_hipp_unfold () { echo "$HippUnfoldDir/hippunfold/sub-$Subject/surf/sub-${Subject}_hemi-$3_space-unfold_den-$2_label-${1}_midthickness.surf.gii"; }
# the CIFTI structure name of a HippUnfold structure/hemisphere: <struct> <H>
_hipp_cstruct () {
	case "$1.$2" in
	hipp.L) echo HIPPOCAMPUS_LEFT;; hipp.R) echo HIPPOCAMPUS_RIGHT;;
	dentate.L) echo HIPPOCAMPUS_DENTATE_LEFT;; dentate.R) echo HIPPOCAMPUS_DENTATE_RIGHT;;
	*) echo "ERROR: unknown HippUnfold structure $1" >&2; return 1;;
	esac
}

# the sampling layer of one structure/hemisphere at the finest density, from its
# inner and outer surfaces (frac: inner + f(outer-inner); vol: equivolume)
_hipp_layer () {
	local st=$1 den=$2 H=$3 out=$4 inn outs
	inn=$(_hipp_surf $st inner $den $H); outs=$(_hipp_surf $st outer $den $H)
	[ -e "$inn" ] && [ -e "$outs" ] || { echo "ERROR: HippUnfold $st inner/outer not found at den $den: $inn"; exit 1; }
	case $DepthMode in
	frac) ${CARET7DIR}/wb_command -surface-average "$out" -surf "$inn" -weight `awk -v f=$DepthValue 'BEGIN{printf "%.12g",1-f}'` -surf "$outs" -weight $DepthValue ;;
	vol)  ${CARET7DIR}/wb_command -surface-cortex-layer "$inn" "$outs" $DepthValue "$out" ;;
	abs)  echo "ERROR: --equi-abs-distance is not supported for HippUnfold; use --equi-frac-distance / --equi-vol-distance"; exit 1 ;;
	esac
}

# interpolate a per-vertex metric from one HippUnfold density onto another, on the
# surface, in the shared unfold space (cubic, nearest fallback; multi-map):
#   <in.func.gii> <out.func.gii> <src unfold surf> <tgt unfold surf>
_hipp_interp () {
	python3 - "$1" "$2" "$3" "$4" <<'PY'
import sys, numpy as np, nibabel as nib
from scipy.interpolate import griddata
inf, outf, su, tu = sys.argv[1:5]
V = np.asarray(nib.load(inf).agg_data(), dtype=np.float64)
if V.ndim == 1: V = V[:, None]
us = np.asarray(nib.load(su).darrays[0].data, dtype=np.float64)[:, :2]
ut = np.asarray(nib.load(tu).darrays[0].data, dtype=np.float64)[:, :2]
out = np.empty((ut.shape[0], V.shape[1]), np.float32)
for k in range(V.shape[1]):
    o = griddata(us, V[:, k], ut, method='cubic'); bad = ~np.isfinite(o)
    if bad.any(): o[bad] = griddata(us, V[:, k], ut[bad], method='nearest')
    out[:, k] = o
da = [nib.gifti.GiftiDataArray(out[:, k].astype(np.float32), intent='NIFTI_INTENT_NORMAL') for k in range(out.shape[1])]
nib.GiftiImage(darrays=da).to_filename(outf)
PY
}

# downsample a finest-density combined hippocampal dscalar onto a coarser density:
# separate each structure, interpolate it in unfold space, reassemble. <src> <out> <tgt den>
_hipp_downsample () {
	local src=$1 out=$2 tden=$3 finest=$(_hipp_finest) st H cs args="" tmp
	tmp=`mktemp -d "${TmpWork:-/tmp}/_hippds.XXXXXX"`
	for st in $HippStructs ; do for H in L R ; do
		cs=$(_hipp_cstruct $st $H)
		${CARET7DIR}/wb_command -cifti-separate "$src" COLUMN -metric $cs $tmp/s.$st.$H.func.gii
		_hipp_interp $tmp/s.$st.$H.func.gii $tmp/t.$st.$H.func.gii $(_hipp_unfold $st $finest $H) $(_hipp_unfold $st $tden $H)
		args="$args -metric $cs $tmp/t.$st.$H.func.gii"
	done; done
	${CARET7DIR}/wb_command -cifti-create-dense-scalar "$out" $args
	rm -rf $tmp
}

# the -metric arguments that pack the per-(structure,hemisphere) files <pfx>.<st>.<H>.func.gii
_hipp_metricargs () {
	local pfx=$1 st H a=""
	for st in $HippStructs ; do for H in L R ; do a="$a -metric $(_hipp_cstruct $st $H) ${pfx}.$st.$H.func.gii"; done; done
	echo "$a"
}

HippUnfold_Prep () {
	local den=$(_hipp_finest) st H W layer raw rawgd
	[ -d "$HippUnfoldDir" ] || { echo "ERROR: HippUnfold folder not found: $HippUnfoldDir"; exit 1; }
	[ -n "$HippDens" ] || { echo "ERROR: no HippUnfold densities found under $HippUnfoldDir (expected <den>/$Subject.L.hipp_midthickness.<den>.surf.gii)"; exit 1; }
	W=$T1wFolder/${DWIT1wFOLDER}.HippUnfoldNODDI${NTag}
	mkdir -p $W $TmpWork
	log_Msg "HippUnfold prep: structures [$HippStructs], den $den, depth $DepthTag -> $W"
	# gunzip the 4-D data once (reading data.nii.gz directly is slow)
	raw=$TmpWork/_hipp_${Subject}_$$_data.nii; gunzip -c $DWIT1wFolder/data.nii.gz > "$raw"
	trap "rm -f '$raw'" EXIT
	[ -n "$GradNonlin" ] && rawgd=`_volfile "$GradNonlin"`
	# sample the diffusion signal (and grad_dev) on every structure's layer
	for st in $HippStructs ; do for H in L R ; do
		layer=$W/$Subject.$H.${st}_layer.${DepthTag}.${den}.surf.gii
		[ -s "$layer" ] || _hipp_layer $st $den $H "$layer"
		${CARET7DIR}/wb_command -volume-to-surface-mapping "$raw" "$layer" $W/_data.$st.$H.func.gii -cubic
		[ -n "$GradNonlin" ] && ${CARET7DIR}/wb_command -volume-to-surface-mapping "$rawgd" "$layer" $W/_grad.$st.$H.func.gii -cubic
	done; done
	# one combined CIFTI over all structures (native HIPPOCAMPUS[_DENTATE] labels)
	${CARET7DIR}/wb_command -cifti-create-dense-scalar $W/data_${DepthTag}.native.dscalar.nii $(_hipp_metricargs $W/_data)
	[ -n "$GradNonlin" ] && ${CARET7DIR}/wb_command -cifti-create-dense-scalar $W/grad_dev_${DepthTag}.native.dscalar.nii $(_hipp_metricargs $W/_grad)
	rm -f $W/_data.*.func.gii $W/_grad.*.func.gii
	# mask: vertices where every sample is finite and the first volume is > 0
	${CARET7DIR}/wb_command -cifti-reduce $W/data_${DepthTag}.native.dscalar.nii SUM $W/_sum.dscalar.nii
	${CARET7DIR}/wb_command -cifti-math '(x0 > 0) * ((s - s) == 0)' $W/mask_${DepthTag}.native.dscalar.nii -fixnan 0 \
		-var x0 $W/data_${DepthTag}.native.dscalar.nii -select 1 1 -var s $W/_sum.dscalar.nii > /dev/null
	rm -f $W/_sum.dscalar.nii
	${CARET7DIR}/wb_command -cifti-convert -to-nifti $W/mask_${DepthTag}.native.dscalar.nii $W/mask_${DepthTag}.nii.gz -smaller-dims
	${FSLDIR}/bin/fslmaths $W/mask_${DepthTag}.nii.gz -bin $W/mask_${DepthTag}.nii.gz
	local s
	for s in data ${GradNonlin:+grad_dev} ; do
		${CARET7DIR}/wb_command -cifti-math 'x * m' $W/_${s}_masked.dscalar.nii -fixnan 0 \
			-var x $W/${s}_${DepthTag}.native.dscalar.nii -var m $W/mask_${DepthTag}.native.dscalar.nii -select 1 1 -repeat > /dev/null
		${CARET7DIR}/wb_command -cifti-convert -to-nifti $W/_${s}_masked.dscalar.nii $W/${s}_${DepthTag}.nii.gz -smaller-dims
		rm -f $W/_${s}_masked.dscalar.nii
	done
	ln -sfn data_${DepthTag}.nii.gz $W/data.nii.gz
	ln -sfn mask_${DepthTag}.nii.gz $W/nodif_brain_mask.nii.gz
	if [ -n "$GradNonlin" ] ; then ln -sfn grad_dev_${DepthTag}.nii.gz $W/grad_dev.nii.gz ; else rm -f $W/grad_dev.nii.gz ; fi
	cp $DWIT1wFolder/bvals $W/bvals; cp $DWIT1wFolder/bvecs $W/bvecs
	log_Msg "HippUnfold: `${CARET7DIR}/wb_command -cifti-stats $W/mask_${DepthTag}.native.dscalar.nii -reduce COUNT_NONZERO` fit vertices; nifti `${FSLDIR}/bin/fslval $W/data_${DepthTag}.nii.gz dim1`x`${FSLDIR}/bin/fslval $W/data_${DepthTag}.nii.gz dim2`x`${FSLDIR}/bin/fslval $W/data_${DepthTag}.nii.gz dim3`x`${FSLDIR}/bin/fslval $W/data_${DepthTag}.nii.gz dim4`"
	local gd=""; [ -n "$GradNonlin" ] && gd=$W/grad_dev.nii.gz
	rm -f "$raw"; trap - EXIT
	_prep_bins $W "$gd"
}

HippUnfold_Fit () {
	local W=$T1wFolder/${DWIT1wFOLDER}.HippUnfoldNODDI${NTag} gd=""
	[ -n "$GradNonlin" ] && gd=$W/grad_dev.nii.gz
	log_Msg "HippUnfold fit: [$HippStructs] ($Fitter)"
	_fit_bins $W ${CUDIMOTDIR} ${CUDIMOTMODELNAME} "$gd"
}

HippUnfold_Map () {
	local finest=$(_hipp_finest) W F tmpl pair srcn dst den st H cs
	W=$T1wFolder/${DWIT1wFOLDER}.HippUnfoldNODDI${NTag}
	_merge_bins $W
	F=$W.${FitTag}; tmpl=$W/mask_${DepthTag}.native.dscalar.nii
	[ -e "$tmpl" ] || { echo "ERROR: $tmpl not found; run prep first"; exit 1; }
	# fitted maps back to a combined CIFTI on the finest mesh (same grid as the mask)
	local mapset=""
	for pair in mean_fintra:noddi_ficvf mean_fiso:noddi_fiso mean_kappa:noddi_kappa OD:noddi_odi AIC:noddi_AIC BIC:noddi_BIC fmin:noddi_fmin ; do
		srcn=${pair%%:*}; dst=${pair##*:}
		[ `imtest $F/$srcn` = 1 ] || continue
		${CARET7DIR}/wb_command -cifti-convert -from-nifti `_volfile $F/$srcn` $tmpl $W/_fit_${dst}.dscalar.nii
		mapset="$mapset $dst"
	done
	# every density: the finest directly, the coarser by unfold-space resample
	for den in $HippDens ; do
		local RD=$AtlasSpaceFolder/Results/$AtlasSpaceResultsDWIFOLDER/HippUnfold/$den
		local TD=$T1wFolder/$DWIT1wFOLDER/HippUnfold/$den
		mkdir -p $RD $TD
		for dst in $mapset ; do
			local nm=$Subject.hippocampus_${dst}${NTag}${FitSuffix}.${den}.dscalar.nii
			if [ "$den" = "$finest" ] ; then cp $W/_fit_${dst}.dscalar.nii $RD/$nm
			else _hipp_downsample $W/_fit_${dst}.dscalar.nii $RD/$nm $den ; fi
			${CARET7DIR}/wb_command -set-map-names $RD/$nm -map 1 "${Subject}_hippocampus_${dst}${NTag}${FitSuffix}.${den}"
			# per-structure/hemisphere shape.gii (compat with existing HippUnfold outputs)
			for st in $HippStructs ; do for H in L R ; do
				cs=$(_hipp_cstruct $st $H)
				${CARET7DIR}/wb_command -cifti-separate $RD/$nm COLUMN -metric $cs $RD/$Subject.$H.${st}_${dst}${NTag}.${den}.shape.gii
				ln -sfn $RD/$Subject.$H.${st}_${dst}${NTag}.${den}.shape.gii $TD/$Subject.$H.${st}_${dst}${NTag}.${den}.shape.gii
			done; done
			ln -sfn $RD/$nm $TD/$nm
		done
		log_Msg "HippUnfold map: den $den ->$(for d in $mapset; do echo -n " $d"; done)  ($RD)"
	done
	rm -f $W/_fit_*.dscalar.nii
}

#########################################################
# Run
#########################################################

Run () {
	StudyFolder=$1
	if [ ! "`echo ${StudyFolder} | head -c 1`" = "/" ]; then StudyFolder=${PWD}/${StudyFolder}; fi
	Subject=$2
	if [ "`echo ${Subject: -1}`" = "/" ] ; then Subject=`echo ${Subject/%?/}` ; fi
	{
		SetUp;
		log_Msg "Start $CMD for subject: $Subject at `date -R`"
		log_Msg "SPECIES=$SPECIES"
		log_Msg "HostName: $(hostname)"
		log_Msg "PATH: $PATH"
		case "$Fitter" in cudimot-*)
			log_Msg "CUDIMOT: $CUDIMOTDIR"
			log_Msg "LD_LIBRARY_PATH: $LD_LIBRARY_PATH"
			log_Msg "CUDA_VISIBLE_DEVICES: $CUDA_VISIBLE_DEVICES"
			log_Msg "NJOB: $NJOB" ;;
		matlab*)
			log_Msg "MATLAB: $MATLABBIN   NODDI toolbox: $NODDIDIR   Rician: ${RNODDIDIR:-n/a}" ;;
		esac

		# ---- prep
		if [ "$DoPrep" = YES ] ; then
			log_Msg "==== stage prep"
			if [ "$Coord" = hippunfold ] ; then
				# the hippocampal (HippUnfold) route is independent of the cortex:
				# it fits the finest density and downsamples on the surface. The
				# cortical/subcortical/WM fits are not run in this mode.
				[ "$DoDTI" = YES ] && { DTIFit; DiffusionStats; }
				HippUnfold_Prep
			else
			# DTIFit and DiffusionStats are whole-brain volume measures and no
			# part of the surface fit reads them: the per-vertex b0 SNR the
			# validity flag uses is computed from the mapped data itself. The
			# -o of the pre-09-20 script skipped them and the batches made with
			# it carry no DTI; --no-dti restores that, so a rerun can match a
			# cohort processed then (2026-09-23).
			if [ "$DoDTI" = YES ] ; then DTIFit; DiffusionStats
			else log_Msg "skipping DTIFit and DiffusionStats (--no-dti)"; fi
			if [ "$MappingMethod" = "surface" ] ; then
				Prep_Surface
				[ "$CortexOnly" = YES ] || Prep_Subcortical
			else
				Prep_WholeBrain
			fi
			[ "$WMFit" = YES ] && Prep_WM
			fi
		fi
		# ---- fit
		if [ "$DoFit" = YES ] ; then
			log_Msg "==== stage fit ($Fitter)"
			if [ "$Coord" = hippunfold ] ; then
				HippUnfold_Fit
			elif [ "$MappingMethod" = "surface" ] ; then
				Fit_Surface
				[ "$CortexOnly" = YES ] || Fit_Subcortical
			else
				Fit_WholeBrain
			fi
			[ "$Coord" != hippunfold ] && [ "$WMFit" = YES ] && Fit_WM
		fi
		# ---- map
		if [ "$DoMap" = YES ] ; then
			log_Msg "==== stage map"
			if [ "$Coord" = hippunfold ] ; then
				HippUnfold_Map
			else
			if [ "$MappingMethod" = "surface" ] ; then
				Map_Surface
				[ "$CortexOnly" = YES ] || Map_Subcortical
			else
				Map_WholeBrain
			fi
			[ "$WMFit" = YES ] && Map_WM
			log_Msg "RegName: $REGNAME"
			DiffusionSurfaceMapping
			fi
		fi

		log_Msg "Finish $CMD for subject: $Subject at `date -R`"
	}
	exit 0;
}

Run $@;
