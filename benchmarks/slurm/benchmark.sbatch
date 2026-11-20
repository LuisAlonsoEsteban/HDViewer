#!/bin/bash
# Run the benchmark on a whole compute node (exclusive, so other jobs do not distort timings).
#
#   sbatch benchmarks/slurm/benchmark.sbatch                              # every dataset, one job
#   sbatch --export=ALL,DATASETS=mouse_brain benchmarks/slurm/benchmark.sbatch
#   sbatch --export=ALL,DATASETS=syn_100k,REPEATS=1,NAME=smoke_syn100k benchmarks/slurm/benchmark.sbatch
#
# Results go to benchmarks/results/$NAME.jsonl. Resubmitting the same NAME resumes: runs already
# recorded are skipped, so a job that hits the time limit loses at most the run in progress.
#SBATCH --job-name=hdv-bench
#SBATCH --qos=xlong
#SBATCH --exclusive
#SBATCH --nodes=1
#SBATCH --mem=240G
#SBATCH --time=3-00:00:00
#SBATCH --output=benchmarks/slurm/logs/%x-%j.out
set -euo pipefail

REPO="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")/../.." && pwd)}"
cd "$REPO"
DATASETS="${DATASETS:-mouse_brain,syn_100k,syn_500k,syn_2M,syn_6M}"
NAME="${NAME:-run_$(echo "$DATASETS" | tr ',' '_')}"
EXTRA=()
[ -n "${REPEATS:-}" ] && EXTRA+=(--repeats "$REPEATS")
[ -n "${TASKS:-}" ] && EXTRA+=(--tasks "$TASKS")
[ -n "${TOOLS:-}" ] && EXTRA+=(--tools "$TOOLS")

echo "node $(hostname) | $(nproc) cpus | $(free -g | awk '/Mem:/{print $2}') GB | datasets $DATASETS | results $NAME"
# Chromium and Xvfb need a writable /tmp; keep it on the node's local disk
export TMPDIR="${TMPDIR:-/tmp}"
.venv/bin/python benchmarks/run.py --datasets "$DATASETS" --name "$NAME" --resume "${EXTRA[@]}"
.venv/bin/python benchmarks/report.py
