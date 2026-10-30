# PriorR as a server (one computer, the users connect with the browser). Usually PriorR runs on each
# user's computer instead (see README). Base: R on a minimal Ubuntu (rocker), no conda.
FROM rocker/r-ver:4.6.1

LABEL authors="Raquel Romero, Pablo Minguez" \
      description="PriorR v4.2.1: prioritization of disease-linked genetic variants" \
      org.opencontainers.image.source="https://github.com/TBLabFJD/PriorR" \
      org.opencontainers.image.vendor="Bioinformatics Unit - Genetics and Genomics Department, IIS-Fundación Jiménez Díaz" \
      org.opencontainers.image.version="4.2.1"

# In one layer: the system libraries up to date (security fixes), libsodium for the passwords and libuv
# for the file dialog (fs); the R packages of PriorR (binaries of Posit Package Manager, as configured
# by rocker) and testthat for the tests of the image; then the compilers and kernel headers of rocker
# are removed (binary packages don't need them; they were most of the vulnerabilities reported by the
# scanners), keeping the runtime libraries that R and the packages use.
RUN apt-get update && apt-get upgrade -y && apt-get install -y --no-install-recommends libsodium23 libuv1 \
    && install2.r --error --skipinstalled shiny DT dplyr data.table jsonlite ggplot2 shinydashboard \
         shinydashboardPlus shinyWidgets shinyalert shinyjs shinyFiles sodium DBI duckdb commonmark testthat \
    && apt-mark manual libgomp1 libquadmath0 libatomic1 \
    && apt-get purge -y --auto-remove linux-libc-dev binutils gcc g++ gfortran cpp make libc6-dev \
    && rm -rf /var/lib/apt/lists/* /tmp/downloaded_packages /tmp/Rtmp*

COPY app/ /srv/priorr/app/
COPY tests/ /srv/priorr/tests/
COPY docker/priorr-start.sh /usr/local/bin/priorr-start
WORKDIR /srv/priorr/app

# PriorR runs as the user priorr, not root. Mount persistent volumes for:
#   /data       variant files that can be opened from the app (docker run -v /path/to/files:/data ...)
#   /var/priorr users, comments, profiles, projects, saved sessions and the caches (docker run -v priorr:/var/priorr ...)
# Example:  docker run -d -p 8888:8888 -v priorr:/var/priorr -v /path/to/variant/files:/data tblabfjd/priorr
# AlphaMissense (optional, once): docker exec -u priorr <container> Rscript Program_files/scripts/get_alphamissense.R
ENV PRIORR_DIR=/srv/priorr/app/Program_files \
    PRIORR_SHARED_DIR=/var/priorr/shared \
    PRIORR_CACHE_DIR=/var/priorr/cache \
    PRIORR_DATA_DIRS=Data=/data
RUN useradd --system --create-home --shell /usr/sbin/nologin priorr \
    && mkdir -p /var/priorr /data && chown priorr:priorr /var/priorr

EXPOSE 8888

ENTRYPOINT ["/usr/local/bin/priorr-start"]
CMD ["Rscript", "-e", "shiny::runApp('/srv/priorr/app', host = '0.0.0.0', port = 8888)"]
