FROM centos:7.5.1804

ENV MCR_HOME=/usr/local/MATLAB_Runtime \	
	MASKFACE_HOME=/usr/local/maskface2 \
	MCR_CACHE_ROOT=/tmp \
	FSLDIR=/nrgpackages/packages/fsl \
	PATH=/usr/local/miniconda3/bin:/nrgpackages/packages/fsl/bin:/nrgpackages/tools/nrg-improc:$PATH

RUN curl https://nvidia.github.io/nvidia-docker/centos7/nvidia-docker.repo > /etc/yum.repos.d/nvidia-docker.repo && \ 
	yum -y install nvidia-container-toolkit ImageMagick wget bzip2 git bc which libX11 libXt java 

RUN mkdir -p /usr/local/MATLAB_Runtime && \
	mkdir -p /docker_mount && \
	mkdir -p /usr/local/maskface2 && \
	mkdir -p /input && \
	cd /tmp && wget https://repo.continuum.io/miniconda/Miniconda3-4.5.1-Linux-x86_64.sh; chmod +x Miniconda3-4.5.1-Linux-x86_64.sh && \
    ./Miniconda3-4.5.1-Linux-x86_64.sh -u -b -p /usr/local/miniconda3 && \
    rm -rf /tmp/* && \
    pip install --upgrade pip
    
RUN cd /var/local; git clone https://github.com/MIC-DKFZ/HD-BET; cd HD-BET; pip install -e . ; cd HD_BET; \
    wget -O 0.model https://zenodo.org/record/2540695/files/0.model?download=1; \    
    sed -i "s#folder_with_parameter_files.*#folder_with_parameter_files = \'/var/local/HD-BET/HD_BET\'#g" paths.py; \
    pip install fslpy

#RUN yum -y remove  wget bzip2 git

COPY  mcr/for_redistribution_files_only/* /usr/local/maskface2/
COPY  mcr/MATLAB_MCR /usr/local/MATLAB_MCR/
COPY nrg-improc /nrgpackages/tools/nrg-improc/
COPY fsl	/nrgpackages/packages/fsl/

WORKDIR /docker_mount
