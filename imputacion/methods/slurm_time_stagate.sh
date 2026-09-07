#!/bin/bash
#SBATCH --job-name=t_stagate
#SBATCH --output=logs/t_stagate_%j.out
#SBATCH --error=logs/t_stagate_%j.err
#SBATCH --gres=gpu:p40:1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=01:00:00

PYTHON_BIN=/scratch/tcondegerma/conda-envs/tfm-impute/bin/python

echo "Python en uso: $PYTHON_BIN"
$PYTHON_BIN -c "import anndata, numpy, scipy; print('imports OK')" || { echo "FALLO IMPORTS"; exit 1; }


cd /scratch/tcondegerma/VisiumVisualizer/imputacion/methods

echo "==================================================================="
echo "JOB: $SLURM_JOB_ID   NODO: $(hostname)   FECHA: $(date)"
echo "METODO: STAGATE_pyG   MODO: solo imputacion (sin evaluacion)"
echo "GPU asignada:"
nvidia-smi --query-gpu=name,memory.total,memory.used --format=csv
echo "==================================================================="

/usr/bin/time -v python time_methods.py --metodo stagate --repeticiones 5

echo "==================================================================="
echo "Fin: $(date)"
echo "==================================================================="
echo "==================================================================="
echo "SACCT (resumen del job):"
sacct -j $SLURM_JOB_ID --format=JobID,JobName,Partition,Elapsed,TotalCPU,MaxRSS,MaxVMSize,ReqMem,State,ExitCode --units=M
echo "==================================================================="
seff $SLURM_JOB_ID 2>/dev/null || echo "[seff no disponible]"