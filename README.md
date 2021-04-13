## Summary
Facemasking2 algorithm to de-identify high resolution MRI of the head. This is a major update of the original <a href="https://download.nrg.wustl.edu/pub/FaceMasking/">facemasking</a> tool. 

## Installing as a local tool
Environment requirements:
1. 64-bit Linux
2. <a href="https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FslInstallation">FSL 5 or 6</a> on path
3. MATLAB R2016 or later on path, or <a href="https://www.mathworks.com/products/compiler/matlab-runtime.html">MATLAB Runtime (MCR) v.9.4</a> installed locally.
4. <a href="https://imagemagick.org/index.php">ImageMagick convert</a>, bc, which, zip and unzip on path
5. (Optional) HD BET tool: https://github.com/NeuroAI-HD/HD-BET

**Installation:** copy the contents of nrg-improc to a local directory and make sure it's on path. For an MCR installation, also copy the contents of 'mcr/for_redistribution_files_only' directory (that would be your deployed application directory). 

## Installing as a Docker container
Make sure the up to date version of <a href="https://docs.docker.com/get-docker/">Docker</a> is installed and Docker daemon is running.<br>
**sudo docker pull** xnat/facemasking2:0.9

## Running
For a full list of options, run **mask_face** without parameters. If you wish to use HDBET to mask out the brain, GPU is recommended. 

1. As a standalone application with MATLAB on path:<br>
**mask_face** <dicom_series_dir> [options]
2. As a deployed MCR application:<br>
**mask_face_nomatlab** <dicom_series_dir> -mcr_home <Matlab runtime environment v 9.4 installation dir> -deployed_home <deployed app dir copied from mcr/for_redistribution_files_only> <dicom_series_dir> [options]
3. As a Docker container: <br>
**sudo docker run** -u $(id -u ${USER}):$(id -g ${USER}) -v \`pwd\`:/docker_mount --rm xnat/facemasking2:0.9 **mask_face_nomatlab** <dicom_series_dir> [options] <br>
To add GPU support to run HDBET, use --gpus <gpu_id|all> option with **docker run** command.

## Reviewing output
1. A mosaic of the facemasked image is saved under maskface/<original_dicom_dir_name>\_normfilter.png
2. A 3D surface rendering of facemasked image is saved under maskface/<original_dicom_dir_name>\_normfilter_surf.png
3. Facemasked image in Analyze 7.5 format is saved under maskface/<original_dicom_dir_name>\_full\_normfilter.img
4. MATLAB output logs are saved under maskface/*.out
5. Defaced DICOM's are saved under DICOM_DEFACED dir.
