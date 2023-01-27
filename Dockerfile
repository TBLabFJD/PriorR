## Define base image where to build our image
FROM nfcore/base:1.13.3 
## Documentation of the image
LABEL authors="Raquel Romero" \
      description="PriorR"

# Install mamba to manage conda (recommended when installing a high number of dependencies)
RUN conda install mamba -c conda-forge

# Install the conda environment
COPY environment.yml /
RUN mamba env create -f /environment.yml && conda clean -a


# Add conda installation dir to PATH (instead of doing 'conda activate')
ENV PATH /opt/conda/envs/nextflow_r/bin:$PATH

# Dump the details of the installed packages to a file for posterity
RUN conda env export --name minimal_docker > minimal_docker.yml


CMD