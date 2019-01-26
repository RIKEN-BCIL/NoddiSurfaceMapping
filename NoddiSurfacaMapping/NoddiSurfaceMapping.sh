#! /bin/bash

set -e
CMD=`echo $0 | sed -e 's/^\(.*\)\/\([^\/]*\)/\2/'`

#########################################################
# Usage & Exit
#########################################################

UsageExit () {

 echo "Calculate NODDI and do surface-mapping for HCP data"
 echo "Usage: $CMD <StudyFolder> <SubjectID 1> <SubjectID 2> ..."
 echo ""
 echo "  Options:"
 echo "    -a <num> : species atlas (0. Human [default], 1. Macaque, 2. Marmoset)"
 echo "    -M       : RegName=MSMAll (default: MSMSulc)"
 echo "    -t <num>,<num>,<num> : b-value upper and lower threshold, and b=0 upper threshold"
 echo ""
 exit 1;

}

#########################################################
# Setup
#########################################################

SetUp () {

# HCP PIPELINE
HCPPIPEDIR=/mnt/pub/devel/git/Pipelines
source $HCPPIPEDIR/Examples/Scripts/SetUpHCPPipeline_RIKEN.sh
source $HCPPIPEDIR/global/scripts/log.shlib  # Logging related functions

# NODDI, requires  'AMICO', 'Camino' and 'matlab' or 'python and pythonspams'
NODDIHCP="/mnt/pub/devel/HCP-RIKEN/NODDI"				# path to NoddiSurfaceMapping
AMICODIR="/mnt/pub/devel/HCP-RIKEN/NODDI/AMICO-master"		# path to AMICO dir
AMICODATADIR="/tmp"                  	# path to AMICO data dir
AMICOPROTOCOL="PROTOCOL_$$"
RunMode="0"  							# 0: Matlab, 1: Python

# References:
# NODDI: http://mig.cs.ucl.ac.uk/index.php?n=Tutorial.NODDImatlab
# AMICO: https://github.com/daducci/AMICO

RegName="MSMSulc"
Species="0"
thr="3100,100,50" #b-value upper and lower threshold, b=0 upper threshold for HCP
while getopts Ma:t: OPT
 do
 case "$OPT" in
   "a" ) export Species="$OPTARG";;
   "M" ) export RegName="MSMAll";;
   "t" ) thr="$OPTARG";;
    * )  Usage_exit;;
 esac
done;
shift `expr $OPTIND - 1`

# Folder Name
T1wFOLDER="T1w"
DWIT1wFOLDER="Diffusion"
DWINativeFOLDER="Diffusion"
AtlasSpaceNativeFOLDER="Native"
AtlasSpaceFOLDER="MNINonLinear"
AtlasSpaceResultsDWIFOLDER="$DWIT1wFOLDER"
AtlasSpaceFOLDER="MNINonLinear"

FreeSurferSubjectFolder=$StudyFolder/$Subject/$T1wFOLDER
FreeSurferSubjectID=$Subject
DWIT1wFolder=$StudyFolder/$Subject/$T1wFOLDER/$DWIT1wFOLDER
DWINativeFolder=$StudyFolder/$Subject/$DWINativeFOLDER
DtiRegDir=$DWINativeFolder/reg
DWIT1wFolder=$StudyFolder/$Subject/$T1wFOLDER/$DWIT1wFOLDER
T1wFolder=$StudyFolder/$Subject/$T1wFOLDER
AtlasSpaceFolder=$StudyFolder/$Subject/$AtlasSpaceFOLDER
AtlasSpaceNativeFolder=$AtlasSpaceFolder/$AtlasSpaceNativeFOLDER
AtlasSpaceResultsDWIFolder=$AtlasSpaceFolder/Results/$AtlasSpaceResultsDWIFOLDER
AtlasSpaceFolder=$StudyFolder/$Subject/$AtlasSpaceFOLDER

# Surface mapping
ribbonLlabel=3
ribbonRlabel=42
ROIFolder=$AtlasSpaceFolder/ROIs

# Which SPECIES
case $Species in
 0)	export SPECIES=Human
	HighResMesh=164
	LowResMeshes=32  # Separate with "@" if needed multiple meshes (e.g. 32@10) with the grayordinate mesh at the last
	BrainOrdinatesResolution=2
	;;
 1) 	export SPECIES=Macaque
	HighResMesh=164
	LowResMeshes=10  # Separate with "@" if needed multiple meshes (e.g. 32@10) with the grayordinate mesh at the last
	BrainOrdinatesResolution=1.25
	;;
 2)	export SPECIES=Marmoset
	HighResMesh=164
	LowResMeshes=2  # Separate with "@" if needed multiple meshes (e.g. 32@10) with the grayordinate mesh at the last
	BrainOrdinatesResolution=1.0
	;;
 *) echo "Not yet supportted atlas species: $Species"; exit 1
esac

NODDIMappingFWHM="`echo "$BrainOrdinatesResolution * 2.5" | bc -l`"
NODDIMappingSigma=`echo "$NODDIMappingFWHM / ( 2 * ( sqrt ( 2 * l ( 2 ) ) ) )" | bc -l`
SmoothingFWHM="$BrainOrdinatesResolution"
SmoothingSigma=`echo "$SmoothingFWHM / ( 2 * ( sqrt ( 2 * l ( 2 ) ) ) )" | bc -l`

if [ ! -e "$AtlasSpaceFolder" ] ; then
 echo "Error: Cannot find $AtlasSpaceFolder"; exit 1;
fi

if [ ! "`imtest $AtlasSpaceFolder/ribbon.nii.gz`" = 1 ] ; then
 echo "ERROR: cannot find ribbon.nii.gz in $AtlasSpaceFolder"; exit 1;
fi

}


#########################################################
# Calculate DTI and NODDI and do surface mapping
#########################################################
DTIFit () {

log_Msg "Start: DTIFit"

thr=(`echo $thr | sed -e 's/,/ /g'`)
buthresh="${thr[0]}" # b-value threshold for DTI
blthresh="${thr[1]}" # b-value threshold for DTI
b0thresh="${thr[2]}"   # b0 volume threshold for DTI b0

log_Msg "b-value upper threshhold: $buthresh"
log_Msg "b-value lower threshhold: $blthresh"
log_Msg "b=0 volume threshhold: $b0thresh"


if [ -e $DWIT1wFolder/dti_dwi.txt ] ; then rm $DWIT1wFolder/dti_dwi.txt; fi
j=0;for i in `cat $DWIT1wFolder/bvals` ; do if [ `echo $i | awk '{printf "%d",$1}'` -le $buthresh ] && [ `echo $i | awk '{printf "%d",$1}'` -ge $blthresh ] ; then j=`zeropad $j 4`; echo -n "vol${j} ">> $DWIT1wFolder/dti_dwi.txt;fi;j=`expr $j + 1`;done
if [ -e $DWIT1wFolder/dti_b0.txt ] ; then rm $DWIT1wFolder/dti_b0.txt; fi
j=0;for i in `cat $DWIT1wFolder/bvals` ; do if [ `echo $i | awk '{printf "%d",$1}'` -le $b0thresh ] ; then j=`zeropad $j 4`; echo -n "vol${j} ">> $DWIT1wFolder/dti_b0.txt;fi;j=`expr $j + 1`;done
if [ -e $DWIT1wFolder/dti_vol.txt ] ; then rm $DWIT1wFolder/dti_vol.txt; fi
cat $DWIT1wFolder/dti_b0.txt $DWIT1wFolder/dti_dwi.txt | sort > $DWIT1wFolder/dti_vol.txt

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
dtifit  -k $DWIT1wFolder/dti_data.nii.gz -o $DWIT1wFolder/dti -m $DWIT1wFolder/nodif_brain_mask -r $DWIT1wFolder/dti_bvecs -b $DWIT1wFolder/dti_bvals --sse
imrm $DWIT1wFolder/vol????.nii.gz $DWIT1wFolder/dti_data.nii.gz
rm $DWIT1wFolder/dti_bvecstmp $DWIT1wFolder/dti_bvecs $DWIT1wFolder/dti_bvals $DWIT1wFolder/dti_vollist.txt

}


NODDIFit () {

log_Msg "Start: NODDIFit"

protocol=$AMICOPROTOCOL
subjdir=subject_$$

if [ -e $AMICODATADIR/$protocol ] ; then
 rm -rf $AMICODATADIR/$protocol;
 if [ "$?" = "1" ] ; then echo "ERROR: canot remove $AMICODATADIR/$protocol. Exit."; exit 1; fi
 echo "Re-newing protocol directory: $AMICODATADIR/$protocol"
fi
mkdir $AMICODATADIR/$protocol
mkdir -p $AMICODATADIR/$protocol/$subjdir
fslchfiletype NIFTI $DWIT1wFolder/data.nii.gz $AMICODATADIR/$protocol/$subjdir/data.nii
fslchfiletype NIFTI $DWIT1wFolder/nodif_brain_mask.nii.gz $AMICODATADIR/$protocol/$subjdir/nodif_brain_mask.nii

${NODDIHCP}/scripts/fsl2scheme.sh $DWIT1wFolder/bvals $DWIT1wFolder/bvecs 1 $AMICODATADIR/$protocol/$subjdir/dwi.scheme -r

if [ "$RunMode" = "0" ] ; then
	command=$AMICODATADIR/$protocol/$subjdir/noddifit.m;
	if [ -e $command ] ; then rm $command;fi

cat <<EOF >> $command
addpath('${AMICODIR}/matlab');
AMICO_Setup;
AMICO_PrecomputeRotationMatrices();
AMICO_SetSubject('$protocol','$subjdir');
CONFIG.dwiFilename    = fullfile( CONFIG.DATA_path, 'data.nii' );
CONFIG.maskFilename   = fullfile( CONFIG.DATA_path, 'nodif_brain_mask.nii' );
CONFIG.schemeFilename = fullfile( CONFIG.DATA_path, 'dwi.scheme' );
AMICO_LoadData;
AMICO_SetModel('NODDI');
AMICO_GenerateKernels(true);
AMICO_ResampleKernels();
AMICO_Fit();
EOF

	matlab -nodesktop -nosplash -r "run $command;quit;"

elif [ "$RunMode" = "1" ] ; then
	command=$AMICODATADIR/$protocol/$subjdir/noddifit.py;
	if [ -e $command ] ; then rm $command;fi

cat <<EOF >> $command;
import sys;
sys.path.append('$AMICODIR/python');
sys.path.append('$AMICODATADIR');
import amico;
amico.core.setup();
ae = amico.Evaluation("$protocoldir","$subjdir");
ae.load_data(dwi_filename = "data.nii" , scheme_filename = "dwi.scheme", mask_filename = "nodif_brain_mask.nii" , b0_thr = 0);
ae.set_model("NODDI");
ae.generate_kernels();
ae.load_kernels();
ae.fit();
ae.save_results();
EOF

	python $command

else
	echo "ERROR: cannot find which RunMode (matlab or python) should be used!"; exit 1;
fi

if [ ! -e $AMICODATADIR/$protocol/$subjdir/AMICO/NODDI/FIT_ICVF.nii ] ; then
 echo "ERROR: NODDI FIT_ICVF is not calculated. Exit"; exit 1;
fi
fslchfiletype NIFTI_GZ $AMICODATADIR/$protocol/$subjdir/AMICO/NODDI/FIT_ICVF.nii $DWIT1wFolder/noddi_ficvf.nii.gz
fslchfiletype NIFTI_GZ $AMICODATADIR/$protocol/$subjdir/AMICO/NODDI/FIT_ISOVF.nii $DWIT1wFolder/noddi_fiso.nii.gz
fslchfiletype NIFTI_GZ $AMICODATADIR/$protocol/$subjdir/AMICO/NODDI/FIT_OD.nii $DWIT1wFolder/noddi_odi.nii.gz
fslchfiletype NIFTI_GZ $AMICODATADIR/$protocol/$subjdir/AMICO/NODDI/FIT_dir.nii $DWIT1wFolder/noddi_dir.nii.gz
\rm -rf $AMICODATADIR/$protocol/$subjdir
\rm -rf $AMICODATADIR/$protocol

}

DiffusionStats () {

log_Msg "DiffusionStats"

${NODDIHCP}/scripts/dwistats $DWIT1wFolder/data.nii.gz $DWIT1wFolder/bvals $DWIT1wFolder/data $DWIT1wFolder/nodif_brain_mask.nii.gz

}

DiffusionSurfaceMapping () {

log_Msg "Start: DiffusionSurfaceMapping"

fslmaths $AtlasSpaceFolder/ribbon.nii.gz -thr $ribbonLlabel -uthr $ribbonLlabel -bin $AtlasSpaceFolder/ribbon_L.nii.gz
fslmaths $AtlasSpaceFolder/ribbon.nii.gz -thr $ribbonRlabel -uthr $ribbonRlabel -bin $AtlasSpaceFolder/ribbon_R.nii.gz
if [ ! -e $AtlasSpaceFolder/T1w_restore."$BrainOrdinatesResolution".nii.gz ] ; then
 flirt -in $AtlasSpaceFolder/T1w_restore.nii.gz -applyisoxfm "$BrainOrdinatesResolution" -ref $AtlasSpaceFolder/ROIs/Atlas_ROIs."$BrainOrdinatesResolution".nii.gz -o $AtlasSpaceFolder/T1w_restore."$BrainOrdinatesResolution".nii.gz -interp sinc
fi

#fslmaths $DWIT1wFolder/noddi_odi.nii.gz  -mul 3.1415926535897 -div 2 -tan -recip -thr 0 $DWIT1wFolder/noddi_kappa.nii.gz
wb_command -volume-math 'max(1/tan((odi*PI)/2),0)' $DWIT1wFolder/noddi_kappa.nii.gz -var odi $DWIT1wFolder/noddi_odi.nii.gz 1>/dev/null

mkdir -p $AtlasSpaceResultsDWIFolder/RibbonVolumeToSurfaceMapping

${CARET7DIR}/wb_command -volume-warpfield-resample $DWIT1wFolder/data_snr.nii.gz $AtlasSpaceFolder/xfms/acpc_dc2standard.nii.gz $AtlasSpaceFolder/T1w_restore.nii.gz CUBIC $AtlasSpaceResultsDWIFolder/data_snr.nii.gz -fnirt $T1wFolder/T1w_acpc_dc_restore.nii.gz

for Hemisphere in L R ; do
   ${CARET7DIR}/wb_command -volume-to-surface-mapping $AtlasSpaceResultsDWIFolder/data_snr.nii.gz "$AtlasSpaceNativeFolder"/"$Subject"."$Hemisphere".midthickness.native.surf.gii  $AtlasSpaceResultsDWIFolder/RibbonVolumeToSurfaceMapping/"$Subject"."$Hemisphere".data_snr.native.func.gii -myelin-style $AtlasSpaceFolder/ribbon_"$Hemisphere".nii.gz "$AtlasSpaceNativeFolder"/"$Subject"."$Hemisphere".thickness.native.shape.gii "$NODDIMappingSigma"
   ${CARET7DIR}/wb_command -metric-mask $AtlasSpaceResultsDWIFolder/RibbonVolumeToSurfaceMapping/"$Subject"."$Hemisphere".data_snr.native.func.gii "$AtlasSpaceNativeFolder"/"$Subject"."$Hemisphere".roi.native.shape.gii $AtlasSpaceResultsDWIFolder/RibbonVolumeToSurfaceMapping/"$Subject"."$Hemisphere".data_snr.native.func.gii
   ${CARET7DIR}/wb_command -metric-math 'x>10' $AtlasSpaceResultsDWIFolder/RibbonVolumeToSurfaceMapping/"$Subject"."$Hemisphere".goodvertex.native.func.gii -var x  $AtlasSpaceResultsDWIFolder/RibbonVolumeToSurfaceMapping/"$Subject"."$Hemisphere".data_snr.native.func.gii
done

for vol in dti_FA dti_MD noddi_kappa noddi_ficvf ; do
 if [ `imtest $DWIT1wFolder/${vol}.nii.gz` = 1 ] ; then
  # Volume-to-surface-mapping
  ${CARET7DIR}/wb_command -volume-warpfield-resample $DWIT1wFolder/${vol}.nii.gz $AtlasSpaceFolder/xfms/acpc_dc2standard.nii.gz $AtlasSpaceFolder/T1w_restore.nii.gz CUBIC $AtlasSpaceResultsDWIFolder/${vol}.nii.gz -fnirt $T1wFolder/T1w_acpc_dc_restore.nii.gz &>/dev/null
  for Hemisphere in L R ; do
   ${CARET7DIR}/wb_command -volume-to-surface-mapping $AtlasSpaceResultsDWIFolder/${vol}.nii.gz "$AtlasSpaceNativeFolder"/"$Subject"."$Hemisphere".midthickness.native.surf.gii  $AtlasSpaceResultsDWIFolder/RibbonVolumeToSurfaceMapping/"$Subject"."$Hemisphere".${vol}.native.func.gii -myelin-style $AtlasSpaceFolder/ribbon_"$Hemisphere".nii.gz "$AtlasSpaceNativeFolder"/"$Subject"."$Hemisphere".thickness.native.shape.gii "$NODDIMappingSigma"
   ${CARET7DIR}/wb_command -metric-mask $AtlasSpaceResultsDWIFolder/RibbonVolumeToSurfaceMapping/"$Subject"."$Hemisphere".${vol}.native.func.gii $AtlasSpaceResultsDWIFolder/RibbonVolumeToSurfaceMapping/"$Subject"."$Hemisphere".goodvertex.native.func.gii $AtlasSpaceResultsDWIFolder/RibbonVolumeToSurfaceMapping/"$Subject"."$Hemisphere".${vol}.native.func.gii
   ${CARET7DIR}/wb_command -metric-dilate $AtlasSpaceResultsDWIFolder/RibbonVolumeToSurfaceMapping/"$Subject"."$Hemisphere".${vol}.native.func.gii "$AtlasSpaceNativeFolder"/"$Subject"."$Hemisphere".midthickness.native.surf.gii 20 $AtlasSpaceResultsDWIFolder/RibbonVolumeToSurfaceMapping/"$Subject"."$Hemisphere".${vol}.native.func.gii -nearest
   ${CARET7DIR}/wb_command -metric-mask $AtlasSpaceResultsDWIFolder/RibbonVolumeToSurfaceMapping/"$Subject"."$Hemisphere".${vol}.native.func.gii "$AtlasSpaceNativeFolder"/"$Subject"."$Hemisphere".roi.native.shape.gii $AtlasSpaceResultsDWIFolder/RibbonVolumeToSurfaceMapping/"$Subject"."$Hemisphere".${vol}.native.func.gii
   ${CARET7DIR}/wb_command -set-map-name $AtlasSpaceResultsDWIFolder/RibbonVolumeToSurfaceMapping/"$Subject"."$Hemisphere".${vol}.native.func.gii 1 "$Subject"_"$Hemisphere"_"$vol"
   ${CARET7DIR}/wb_command -metric-palette $AtlasSpaceResultsDWIFolder/RibbonVolumeToSurfaceMapping/"$Subject"."$Hemisphere".${vol}.native.func.gii MODE_AUTO_SCALE_PERCENTAGE -pos-percent 4 96 -interpolate true -palette-name videen_style -disp-pos true -disp-neg false -disp-zero false
   #LowResMesh
   for LowResMesh in `echo $LowResMeshes | sed -e 's/@/ /g'`; do
    DownsampleFolder=$AtlasSpaceFolder/fsaverage_LR${LowResMesh}k
  ${CARET7DIR}/wb_command -metric-resample $AtlasSpaceResultsDWIFolder/RibbonVolumeToSurfaceMapping/"$Subject"."$Hemisphere".${vol}.native.func.gii "$AtlasSpaceNativeFolder"/"$Subject"."$Hemisphere".sphere.${RegName}.native.surf.gii "$DownsampleFolder"/"$Subject"."$Hemisphere".sphere."$LowResMesh"k_fs_LR.surf.gii ADAP_BARY_AREA $AtlasSpaceResultsDWIFolder/"$Subject"."$Hemisphere".${vol}."$LowResMesh"k_fs_LR.func.gii -area-surfs "$AtlasSpaceNativeFolder"/"$Subject"."$Hemisphere".midthickness.native.surf.gii "$DownsampleFolder"/"$Subject"."$Hemisphere".midthickness."$LowResMesh"k_fs_LR.surf.gii -current-roi "$AtlasSpaceNativeFolder"/"$Subject"."$Hemisphere".roi.native.shape.gii
  ${CARET7DIR}/wb_command -metric-mask $AtlasSpaceResultsDWIFolder/"$Subject"."$Hemisphere".${vol}."$LowResMesh"k_fs_LR.func.gii "$DownsampleFolder"/"$Subject"."$Hemisphere".atlasroi."$LowResMesh"k_fs_LR.shape.gii $AtlasSpaceResultsDWIFolder/"$Subject"."$Hemisphere".${vol}."$LowResMesh"k_fs_LR.func.gii
    ${CARET7DIR}/wb_command -metric-smoothing "$DownsampleFolder"/"$Subject"."$Hemisphere".midthickness."$LowResMesh"k_fs_LR.surf.gii $AtlasSpaceResultsDWIFolder/"$Subject"."$Hemisphere".${vol}."$LowResMesh"k_fs_LR.func.gii "$SmoothingSigma" $AtlasSpaceResultsDWIFolder/"$Subject"."$Hemisphere".${vol}_s"$SmoothingFWHM"."$LowResMesh"k_fs_LR.func.gii -roi "$DownsampleFolder"/"$Subject"."$Hemisphere".atlasroi."$LowResMesh"k_fs_LR.shape.gii
   done
   #HighResMesh
   ${CARET7DIR}/wb_command -metric-resample $AtlasSpaceResultsDWIFolder/RibbonVolumeToSurfaceMapping/"$Subject"."$Hemisphere".${vol}.native.func.gii "$AtlasSpaceNativeFolder"/"$Subject"."$Hemisphere".sphere.${RegName}.native.surf.gii "$AtlasSpaceFolder"/"$Subject"."$Hemisphere".sphere."$HighResMesh"k_fs_LR.surf.gii ADAP_BARY_AREA $AtlasSpaceResultsDWIFolder/"$Subject"."$Hemisphere".${vol}."$HighResMesh"k_fs_LR.func.gii -area-surfs "$AtlasSpaceNativeFolder"/"$Subject"."$Hemisphere".midthickness.native.surf.gii "$AtlasSpaceFolder"/"$Subject"."$Hemisphere".midthickness."$HighResMesh"k_fs_LR.surf.gii -current-roi "$AtlasSpaceNativeFolder"/"$Subject"."$Hemisphere".roi.native.shape.gii
   ${CARET7DIR}/wb_command -metric-mask $AtlasSpaceResultsDWIFolder/"$Subject"."$Hemisphere".${vol}."$HighResMesh"k_fs_LR.func.gii "$AtlasSpaceFolder"/"$Subject"."$Hemisphere".atlasroi."$HighResMesh"k_fs_LR.shape.gii $AtlasSpaceResultsDWIFolder/"$Subject"."$Hemisphere".${vol}."$HighResMesh"k_fs_LR.func.gii
   ${CARET7DIR}/wb_command -metric-smoothing "$AtlasSpaceFolder"/"$Subject"."$Hemisphere".midthickness."$HighResMesh"k_fs_LR.surf.gii $AtlasSpaceResultsDWIFolder/"$Subject"."$Hemisphere".${vol}."$HighResMesh"k_fs_LR.func.gii "$SmoothingSigma" $AtlasSpaceResultsDWIFolder/"$Subject"."$Hemisphere".${vol}_s"$SmoothingFWHM"."$HighResMesh"k_fs_LR.func.gii -roi "$AtlasSpaceFolder"/"$Subject"."$Hemisphere".atlasroi."$HighResMesh"k_fs_LR.shape.gii
  done

  # Do volume parcel resampling for subcortical gray matter
  ${CARET7DIR}/wb_command -volume-warpfield-resample $DWIT1wFolder/${vol}.nii.gz $AtlasSpaceFolder/xfms/acpc_dc2standard.nii.gz $AtlasSpaceFolder/T1w_restore."$BrainOrdinatesResolution".nii.gz CUBIC $AtlasSpaceResultsDWIFolder/${vol}."$BrainOrdinatesResolution".nii.gz -fnirt $T1wFolder/T1w_acpc_dc_restore.nii.gz &> /dev/null
  ${CARET7DIR}/wb_command -volume-parcel-resampling $AtlasSpaceResultsDWIFolder/${vol}."$BrainOrdinatesResolution".nii.gz "$ROIFolder"/ROIs."$BrainOrdinatesResolution".nii.gz "$ROIFolder"/Atlas_ROIs."$BrainOrdinatesResolution".nii.gz $SmoothingSigma $AtlasSpaceResultsDWIFolder/${vol}_AtlasSubcortical_s"$SmoothingFWHM".nii.gz -fix-zeros
  # Merge surface and subcortical volume to create cifti
  LowReMesh=`echo $LowResMeshes | sed -e 's/@/ /g' | awk '{print $NF}'`
  ${CARET7DIR}/wb_command -cifti-create-dense-scalar $AtlasSpaceResultsDWIFolder/${vol}.dscalar.nii -volume $AtlasSpaceResultsDWIFolder/${vol}_AtlasSubcortical_s"$SmoothingFWHM".nii.gz "$ROIFolder"/Atlas_ROIs."$BrainOrdinatesResolution".nii.gz -left-metric $AtlasSpaceResultsDWIFolder/"$Subject".L.${vol}_s"$SmoothingFWHM"."$LowResMesh"k_fs_LR.func.gii -roi-left "$DownsampleFolder"/"$Subject".L.atlasroi."$LowResMesh"k_fs_LR.shape.gii -right-metric $AtlasSpaceResultsDWIFolder/"$Subject".R.${vol}_s"$SmoothingFWHM"."$LowResMesh"k_fs_LR.func.gii -roi-right "$DownsampleFolder"/"$Subject".R.atlasroi."$LowResMesh"k_fs_LR.shape.gii
  ${CARET7DIR}/wb_command -set-map-names $AtlasSpaceResultsDWIFolder/${vol}.dscalar.nii -map 1 "${Subject}_${vol}"
  ${CARET7DIR}/wb_command -cifti-palette $AtlasSpaceResultsDWIFolder/${vol}.dscalar.nii MODE_AUTO_SCALE_PERCENTAGE $AtlasSpaceResultsDWIFolder/${vol}.dscalar.nii -pos-percent 4 96 -interpolate true -palette-name videen_style -disp-pos true -disp-neg false -disp-zero false
 fi
done

# Convert kappa to odi
for Hemisphere in R L; do
  ${CARET7DIR}/wb_command -metric-math 'max(2*atan(1/kappa)/PI,0)' $AtlasSpaceResultsDWIFolder/"$Subject"."$Hemisphere".noddi_odi${Reg}."$HighResMesh"k_fs_LR.func.gii -var kappa $AtlasSpaceResultsDWIFolder/"$Subject"."$Hemisphere".noddi_kappa${Reg}."$HighResMesh"k_fs_LR.func.gii
done
${CARET7DIR}/wb_command -cifti-math 'max(2*atan(1/kappa)/PI,0)' $AtlasSpaceResultsDWIFolder/noddi_odi${Reg}.dscalar.nii -var kappa $AtlasSpaceResultsDWIFolder/noddi_kappa${Reg}.dscalar.nii

}

#########################################################
# Run
#########################################################

Run () {

if [ "$2" = "" ] ; then UsageExit; fi

StudyFolder=$1;
if [ ! "`echo ${StudyFolder} | head -c 1`" = "/" ]; then StudyFolder=${PWD}/${StudyFolder}; fi
shift;
Subjects="$@";

for Subject in $Subjects ; do
 if [ "`echo ${Subject: -1}`" = "/" ] ; then Subject=`echo ${Subject/%?/}` ;fi
 log_Msg "Start $CMD for subject: $Subject at `date -R`"
 SetUp;
 log_Msg "SPECIES=$SPECIES"
 DTIFit;
 NODDIFit;
 DiffusionStats
 DiffusionSurfaceMapping;
 log_Msg "Finished subject: $Subject at `date -R`"

done

exit 0;

}

Run $@;
