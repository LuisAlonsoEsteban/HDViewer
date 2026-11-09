# Build from hd_viewer/analysis/deconvolution/api/.
# docker build -t ghcr.io/luisalonsoesteban/visium-deconv-spacexr:0.2.0 -f plugins/spacexr/Dockerfile .
FROM rocker/r-ver:4.4.3
ENV PYTHONNOUSERSITE=1
RUN apt-get update && apt-get install -y --no-install-recommends python3 python3-venv libcurl4-openssl-dev libssl-dev libxml2-dev libgsl-dev liblapack-dev libblas-dev gfortran cmake && rm -rf /var/lib/apt/lists/*
RUN python3 -m venv /opt/python
ENV PATH="/opt/python/bin:$PATH"
# Use the GitHub Reference/SpatialRNA API, not the Bioconductor spacexr API.
# spacexr with the Visium HD speed-ups of dmcable/spacexr pull request 206 (doublet mode:
# foreach/doParallel workers, vectorized result gathering), pinned at the PR's head commit.
# Return to dmcable/spacexr once it is merged there.
ARG SPACEXR_REPO=jpromeror/spacexr
ARG SPACEXR_REF=33d375cc5d7b7b5db97ee5cc1d3dd32b682afc9e
RUN R -e 'options(timeout=1200); install.packages(c("remotes", "jsonlite"), repos="https://cloud.r-project.org"); remotes::install_github(paste0("'"$SPACEXR_REPO"'", "@", "'"$SPACEXR_REF"'"), upgrade="never", build_vignettes=FALSE); stopifnot(requireNamespace("spacexr", quietly=TRUE), requireNamespace("Matrix", quietly=TRUE), requireNamespace("jsonlite", quietly=TRUE)); stopifnot(all(c("Reference", "SpatialRNA", "create.RCTD", "run.RCTD") %in% getNamespaceExports("spacexr")))'
WORKDIR /app
COPY plugins/spacexr/requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt
COPY sdk/ /app/
COPY plugins/spacexr/ /app/
ENTRYPOINT ["python", "/app/run.py"]
