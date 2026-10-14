# Build from the segmentation/ directory:
#   docker build -t ghcr.io/luisalonsoesteban/visium-seg-stardist:0.1.0 -f plugins/stardist/Dockerfile .
FROM python:3.11-slim

WORKDIR /app
COPY plugins/stardist/requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Bake the pretrained weights in so runs work offline and as any user
ENV STARDIST_MODELS=/opt/models/stardist
RUN python -c "import shutil; from stardist.models import StarDist2D; \
m = StarDist2D.from_pretrained('2D_versatile_he'); \
shutil.copytree(m.logdir, '/opt/models/stardist/2D_versatile_he')" \
 && chmod -R a+rX /opt/models

COPY sdk/vhd_seg.py plugins/stardist/run.py /app/
ENTRYPOINT ["python", "/app/run.py"]
