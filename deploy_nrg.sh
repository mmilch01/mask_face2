sudo docker tag registry.nrg.wustl.edu/docker/nrg-repo/facemasking2:latest \
    registry.nrg.wustl.edu/docker/nrg-repo/facemasking2:0.9
sudo docker login registry.nrg.wustl.edu
sudo docker push registry.nrg.wustl.edu/docker/nrg-repo/facemasking2:latest:0.9

