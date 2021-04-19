#!/bin/bash

echo "REMEMBER, MUST RUN AS LOCAL USER"

if [ ! -d "nrg_improc" ]; then 
        echo git clone https://github.com/mmilch01/nrg_improc
        git clone https://github.com/mmilch01/nrg_improc
fi        

pushd nrg_improc &> /dev/null
        echo git pull
        git pull
popd &> /dev/null
mkdir -p nrg-improc; rm -rf nrg-improc/*; T=`pwd`/nrg-improc
mkdir -p nrg-improc/ATLAS nrg-improc/matlab/routines nrg-improc/matlab/surf


nrg_improc_tools="run_facemasking2_xnat mask_face mask_face_nomatlab cmpanalyze analyze2dcm dcm2analyze dcm_nii nii_img hdbet_wrapper dcm2nii2013 dcm2nii.ini xnat2loc dcm_sort dcminfo"
atlas_resources="ATLAS/*.mat ATLAS/*.lst ATLAS/CAPIIO.hdr ATLAS/CAPIIO.img"
matlab_surf_resources="RectangularMesh.m blur3_thin.m maskVol.m mask_surf.m plane_eq.m projectVol.m showmesh.m test_rect_mesh.m"
matlab_routines_resources="avw_hdr_make.m avw_hdr_read.m avw_hdr_write.m avw_img_read.m avw_img_write.m blur3.m dispvol.m dispvol3D.m get_seg_stats.m saveVol.m save_vol.m select_threshold.m"


pushd nrg_improc &> /dev/null
        echo cp -f $nrg_improc_tools $T/
        cp -f $nrg_improc_tools $T/
        echo cp -f $atlas_resources $T/ATLAS/
        cp -f $atlas_resources $T/ATLAS/
        
        pushd matlab/routines &> /dev/null
                echo cp -rf $matlab_routines_resources $T/matlab/routines/
                cp -rf $matlab_routines_resources $T/matlab/routines/
        popd &> /dev/null
        pushd matlab/surf &> /dev/null
                echo cp -rf $matlab_surf_resources $T/matlab/surf/
                cp -rf $matlab_surf_resources $T/matlab/surf/
        popd &> /dev/null
popd &> /dev/null

rm -rf fsl

mkdir -p fsl/bin; T=`pwd`/fsl/bin

fslbin=`which flirt`
fslbin=${fslbin%/flirt} 
fsl_home=${fslbin%/bin}
mkdir -p fsl/data/standard fsl/lib
#fslpython/envs/fslpython/bin/python3.7
fsl_resources=(data/standard/MNI152lin_T1_1mm.nii.gz LICENCE etc/flirtsch/ident.mat)
for res in ${fsl_resources[*]}; do
        mkdir -p fsl/`dirname $res`
        echo cp $fsl_home/$res fsl/$res
        cp $fsl_home/$res fsl/$res
done


if [ ! -d "$fslbin" -o ! -f "$fsl_home/data/standard/MNI152lin_T1_1mm.nii.gz" ]; then
        echo "prepare_context.sh ERROR: FSL tools/resources not found on path, exiting!"
        exit -1
fi

fsl_tools="flirt fslmaths img2imgcoord fslcpgeom convert_xfm fslreorient2std bet remove_ext tmpnam imtest fslval fslhd avscale fslswapdim imcp fslchfiletype fslchfiletype_exe"
pushd $fslbin &> /dev/null
        echo cp -f $fsl_tools $T/
        cp -f $fsl_tools $T/
popd &> /dev/null

sed -i "s|.*fslpython.*|#!/usr/local/miniconda3/bin/python|g" fsl/bin/imcp

fsl_libs=(libgfortran.so.3 libopenblas.so.0)
libs64=(libquadmath.so.0.0.0)
libs64_targ=(libquadmath.so.0)


for lib in ${fsl_libs[*]}; do
        cp $fsl_home/lib/$lib fsl/lib/
done

for i in ${!libs64[*]}; do
        cp /lib64/${libs64[i]} fsl/lib/${libs64_targ[i]}
done

echo chmod -R o+rX `pwd`
chmod -R o+rX `pwd`
