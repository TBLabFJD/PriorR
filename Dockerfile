## Define base image where to build our image
FROM condaforge/mambaforge:4.12.0-2

## Documentation of the image
LABEL authors="Raquel Romero" \
      description="PriorR"

# Install mamba to manage conda (recommended when installing a high number of dependencies)
#RUN conda install mamba -c conda-forge

# Install the conda environment
COPY environment.yml /
RUN mamba env create -f /environment.yml && conda clean -a

# Add conda installation dir to PATH (instead of doing 'conda activate')
ENV PATH /opt/conda/envs/Priorr/bin:$PATH


# Dump the details of the installed packages to a file for posterity
RUN conda env export --name minimal_docker > minimal_docker.yml


RUN apt-get update && apt-get install -y \
gcc \
g++ \
make \
perl \
cpanminus \
libmysqlclient-dev \
unzip \
libbz2-dev \
liblzma-dev \
bcftools \
tabix \
bzip2 \
bedtools

RUN cpanm Archive::Zip DBD::mysql LWP::Simple


RUN wget https://github.com/Ensembl/ensembl-vep/archive/release/105.zip && unzip 105.zip
WORKDIR ensembl-vep-release-105
RUN perl INSTALL.pl -a a --NO_UPDATE --CACHE_VERSION 105


ADD AutoMap/ /home/docker/AutoMap/
RUN chmod 777 /home/docker


# copy the app directory into the image
COPY app/* /home/app/
COPY app/* /session/app/



WORKDIR /

EXPOSE 8888

CMD R -e "shiny::runApp('/home/app/PRIORR3.1_docker.R', host='0.0.0.0', port=8888)"




# docker build -t priorr_gonzalo .
