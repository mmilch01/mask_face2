#!/bin/bash

echo "building Docker image for Facemasking2"
echo docker build . -t registry.nrg.wustl.edu/docker/nrg-repo/facemasking2 
docker build . -t registry.nrg.wustl.edu/docker/nrg-repo/facemasking2 
