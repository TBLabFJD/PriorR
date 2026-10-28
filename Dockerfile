FROM condaforge/miniforge3:24.9.2-0

LABEL authors="Raquel Romero, Pablo Minguez" \
      description="PriorR v4.0: prioritization of disease-linked genetic variants" \
      org.opencontainers.image.source="https://github.com/TBLabFJD/PriorR" \
      org.opencontainers.image.vendor="Bioinformatics Unit - Genetics and Genomics Department, IIS-Fundación Jiménez Díaz" \
      org.opencontainers.image.version="4.0"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a -y
ENV PATH=/opt/conda/envs/Priorr/bin:$PATH

COPY app/ /srv/priorr/app/
COPY tests/ /srv/priorr/tests/
WORKDIR /srv/priorr/app

# PriorR usually runs on each user's computer (see README); this image runs it as a server.
# Mount persistent volumes for:
#   /data       variant files that can be opened from the app (docker run -v /path/to/files:/data ...)
#   /var/priorr users, comments, profiles, saved sessions and the WGS cache (docker run -v priorr:/var/priorr ...)
# Example:  docker run -d -p 8888:8888 -v priorr:/var/priorr -v /path/to/variant/files:/data tblabfjd/priorr
# AlphaMissense (optional, once): docker exec <container> Rscript Program_files/scripts/get_alphamissense.R
ENV PRIORR_DIR=/srv/priorr/app/Program_files \
    PRIORR_SHARED_DIR=/var/priorr/shared \
    PRIORR_CACHE_DIR=/var/priorr/cache \
    PRIORR_DATA_DIRS=Data=/data
RUN mkdir -p /var/priorr /data

EXPOSE 8888

CMD ["Rscript", "-e", "shiny::runApp('/srv/priorr/app', host = '0.0.0.0', port = 8888)"]
