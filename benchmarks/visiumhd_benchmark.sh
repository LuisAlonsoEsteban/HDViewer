#!/bin/bash
#SBATCH --job-name=visiumhd_benchmark
#SBATCH --output=/scratch/tcondegerma/VisiumVisualizer/benchmarks/results/slurm_%j.out
#SBATCH --error=/scratch/tcondegerma/VisiumVisualizer/benchmarks/results/slurm_%j.err
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH --qos=regular

echo "========================================"
echo "  VisiumHD Benchmark"
echo "  Job ID   : $SLURM_JOB_ID"
echo "  Node     : $SLURM_NODELIST"
echo "  Start    : $(date '+%Y-%m-%d %H:%M:%S')"
echo "========================================"

source /dipc/tcondegerma/.conda/envs/tfm/etc/conda/activate.d/env_vars.sh 2>/dev/null || true
export PATH="/dipc/tcondegerma/.conda/envs/tfm/bin:$PATH"
export CONDA_PREFIX="/dipc/tcondegerma/.conda/envs/tfm"

cd /scratch/tcondegerma/VisiumVisualizer

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK
export TF_CPP_MIN_LOG_LEVEL=2

echo ""
echo "--- Python ---"
python --version
echo ""

# --- Benchmark principal ---
/usr/bin/time -v python benchmarks/benchmark_visiumhd.py
EXIT_CODE=$?

echo ""
echo "========================================"
echo "  Fin: $(date '+%Y-%m-%d %H:%M:%S')"
echo "  Exit code: $EXIT_CODE"
echo "========================================"