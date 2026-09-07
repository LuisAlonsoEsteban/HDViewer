#!/bin/bash
#SBATCH --job-name=t_load
#SBATCH --output=logs/t_load_%j.out
#SBATCH --error=logs/t_load_%j.err
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:15:00

PYTHON_BIN=/scratch/tcondegerma/conda-envs/tfm-impute/bin/python

echo "Python en uso: $PYTHON_BIN"
$PYTHON_BIN -c "import anndata, numpy, scipy; print('imports OK')" || { echo "FALLO IMPORTS"; exit 1; }

cd /scratch/tcondegerma/VisiumVisualizer/imputacion/methods

echo "JOB: $SLURM_JOB_ID   NODO: $(hostname)   FECHA: $(date)"
/usr/bin/time -v python time_load.py
seff $SLURM_JOB_ID 2>/dev/null || echo "[seff no disponible]"