# Build from hd_viewer/board/operations/layer/deconvolution/api/.
# docker build -t ghcr.io/luisalonsoesteban/visium-deconv-cell2location:0.1.0 -f plugins/cell2location/Dockerfile .
FROM python:3.11-slim
ENV PYTHONNOUSERSITE=1
WORKDIR /app
COPY plugins/cell2location/requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt
COPY sdk/ /app/
COPY plugins/cell2location/ /app/
ENTRYPOINT ["python", "/app/run.py"]
