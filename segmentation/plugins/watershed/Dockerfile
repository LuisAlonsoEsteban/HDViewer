# Build from the segmentation/ directory:
#   docker build -t ghcr.io/luisalonsoesteban/visium-seg-watershed:0.1.0 -f plugins/watershed/Dockerfile .
FROM python:3.11-slim

WORKDIR /app
COPY plugins/watershed/requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

COPY sdk/vhd_seg.py plugins/watershed/run.py /app/
ENTRYPOINT ["python", "/app/run.py"]
