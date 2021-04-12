#!/bin/bash

if [ -z "$1" ]; then 
	echo "usage: debug_image.sh <container args>"
	exit -1
fi

echo docker run -u $(id -u ${USER}):$(id -g ${USER}) -v `pwd`:/docker_mount -it --rm --gpus all registry.nrg.wustl.edu/docker/nrg-repo/facemasking2:latest mask_face_nomatlab $@
docker run -u $(id -u ${USER}):$(id -g ${USER}) -v `pwd`:/docker_mount -it --rm --gpus all registry.nrg.wustl.edu/docker/nrg-repo/facemasking2:latest mask_face_nomatlab $@

