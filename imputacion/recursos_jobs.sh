#!/bin/bash
# recursos_jobs.sh <jobid> [<jobid> ...]
#
# Contabilidad de SLURM de jobs YA TERMINADOS. Tiene que ejecutarse fuera del
# job: mientras el job corre, su step sigue en RUNNING y sacct devuelve MaxRSS
# y TotalCPU vacios (por eso las tablas de los .out salian en blanco).
#
# Sin argumentos usa los jobs t_clust_* del usuario en las ultimas 24 h.

set -u
FMT=JobID,JobName%20,Elapsed,TotalCPU,MaxRSS,MaxVMSize,ReqMem,State,ExitCode

if [ $# -gt 0 ]; then
    JOBS=$(IFS=,; echo "$*")
    sacct -j "$JOBS" --format=$FMT --units=M
else
    sacct -u "$USER" -S "$(date -d '24 hours ago' +%Y-%m-%dT%H:%M)" \
          --name=t_clust_knn,t_clust_spamean,t_clust_stagate,t_clust_magic \
          --format=$FMT --units=M
fi

echo
echo "--- pico de RSS por job (linea .batch: es la que lo trae) ---"
for j in "$@"; do
    seff "$j" 2>/dev/null | grep -E "Job ID|Memory Utilized|CPU Utilized|Job Wall" || true
done
