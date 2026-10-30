# Build from hd_viewer/board/operations/layer/deconvolution/api/.
# docker build -t ghcr.io/luisalonsoesteban/visium-deconv-music:0.1.0 -f plugins/music/Dockerfile .
FROM rocker/r-ver:4.4.3
ENV PYTHONNOUSERSITE=1
RUN apt-get update && apt-get install -y --no-install-recommends python3 python3-venv libcurl4-openssl-dev libssl-dev libxml2-dev libhdf5-dev libgsl-dev && rm -rf /var/lib/apt/lists/*
RUN python3 -m venv /opt/python
ENV PATH="/opt/python/bin:$PATH"
RUN R -e 'install.packages(c("BiocManager", "remotes", "jsonlite"), repos="https://cloud.r-project.org"); BiocManager::install(c("SingleCellExperiment", "TOAST"), ask=FALSE, update=FALSE); remotes::install_github("xuranw/MuSiC", upgrade="never")'
WORKDIR /app
COPY plugins/music/requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt
COPY sdk/ /app/
COPY plugins/music/ /app/
ENTRYPOINT ["python", "/app/run.py"]
