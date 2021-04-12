sudo docker tag egistry.nrg.wustl.edu/docker/nrg-repo/facemasking2:latest \
    registry.nrg.wustl.edu/docker/nrg-repo/facemasking2:latest:1.0
sudo docker login registry.nrg.wustl.edu
sudo docker push registry.nrg.wustl.edu/docker/nrg-repo/facemasking2:latest:1.0

