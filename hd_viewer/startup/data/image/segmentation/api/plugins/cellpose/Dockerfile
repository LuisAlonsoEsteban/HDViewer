# Build from hd_viewer/startup/data/image/segmentation/api/:
#   docker build -t ghcr.io/luisalonsoesteban/visium-seg-cellpose:0.1.0 -f plugins/cellpose/Dockerfile .
# CPU-only torch keeps the image small; swap the index URL in requirements.txt for CUDA.
FROM python:3.11-slim
# Never pick up packages from a mounted home directory (~/.local)
ENV PYTHONNOUSERSITE=1

WORKDIR /app
COPY plugins/cellpose/requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Bake the weights in so runs work offline and as any user
ENV CELLPOSE_LOCAL_MODELS_PATH=/opt/models/cellpose
RUN python -c "from cellpose import models; \
models.Cellpose(gpu=False, model_type='cyto3'); \
models.Cellpose(gpu=False, model_type='nuclei')" \
 && chmod -R a+rX /opt/models

COPY sdk/vhd_seg.py plugins/cellpose/run.py /app/
ENTRYPOINT ["python", "/app/run.py"]
