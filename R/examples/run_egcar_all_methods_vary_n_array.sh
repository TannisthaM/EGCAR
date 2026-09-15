#!/bin/bash
#SBATCH --job-name=egcar_all_n
#SBATCH --output=logs/egcar_all_n_%A_%a.out
#SBATCH --error=logs/egcar_all_n_%A_%a.err
#SBATCH --array=1-150%40
#SBATCH --time=12:00:00
#SBATCH --partition=caslake
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=20G
#SBATCH --account=pi-cdonnat

# EXTENDED GRID: rank now sweeps 1..10 (was 1,2,5) and every (panel, n, rank)
# combination is additionally crossed with 3 signal strengths (0.3, 0.5, 0.8;
# was one fixed signal=0.8). n/p_per_block/panel values are UNCHANGED.
#   2 panels x 10 ranks x 3 signals x 6 n-values x 10 reps = 3,600 tasks
# (10x the original 360). task_info now returns a 6th field (signal).
#
# --array=1-150%40 covers one 150-task chunk at up to 40 concurrent tasks
# ("parallelize over as many arrays as possible", raised from the original
# script's 60-task/%5 default). Submit 24 chunks to cover all 3,600 -- see
# the offset loop below. The real ceiling on concurrency is whichever is
# smaller of your account's QOS limits and RCC's site-wide MaxArraySize:
#   sacctmgr show qos format=Name,MaxSubmitJobsPerUser,MaxTRESPerUser
#   scontrol show config | grep -i maxarraysize
# If those allow more, raise both the --array upper bound and the %N
# concurrency together.
#
# 20G/12h is UNVERIFIED at this task count and grid (rank now goes up to 10,
# which increases the number of CV candidates evaluated per task compared to
# the rank<=5 grid this script was originally sized for). Run a one-task
# smoke submission first (--array=1-1) and check
# `sacct -j JOBID --format=MaxRSS,Elapsed` before trusting the full run.
#
# All seven CV methods use the SAME fold-worker pool in a task:
# min(EGCAR_CV_WORKERS, 5 folds, SLURM_CPUS_PER_TASK - 1). Default: 5 workers.
# Methods run sequentially; BLAS is single-threaded.
# Packages must be installed once BEFORE submission, not by concurrent jobs.
# egcar (not SGCAR) is now required -- run `Rscript ... install_packages`
# with EGCAR_LOCAL_SOURCE set to your egcar zip/tarball first if needed.

set -euo pipefail

export EGCAR_R_MODULE="${EGCAR_R_MODULE:-R/4.2.0}"
module load "$EGCAR_R_MODULE"
export R_LIBS_USER="${R_LIBS_USER:-$HOME/Rlibs}"
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1 BLIS_NUM_THREADS=1

if [[ -n "${EGCAR_ROOT:-}" ]]; then
  ROOT="$EGCAR_ROOT"
else
  : "${SCRATCH:?SCRATCH is unset; set EGCAR_ROOT to the experiment directory}"
  ROOT="$SCRATCH/$USER/CCA-experiments"
fi

RSCRIPT="$ROOT/r/experiments/cluster/run_egcar_all_methods_vary_n.R"
[[ -f "$RSCRIPT" ]] || { echo "Missing R script: $RSCRIPT" >&2; exit 2; }

RUN_ID="${EGCAR_RUN_ID:-job_${SLURM_ARRAY_JOB_ID:-${SLURM_JOB_ID:-manual}}}"
[[ "$RUN_ID" =~ ^[A-Za-z0-9][A-Za-z0-9_.-]*$ ]] || { echo "Invalid EGCAR_RUN_ID" >&2; exit 2; }

OUTDIR="$ROOT/egcar_all_methods_vary_n_outputs/$RUN_ID"
mkdir -p "$ROOT/logs" "$OUTDIR"/{metrics,fits,cv,plots,logs,metadata,completed,progress,loading_visualizations}

cd "$ROOT"

ACTION="${EGCAR_ACTION:-worker}"
if [[ "$ACTION" == "aggregate" ]]; then
  Rscript "$RSCRIPT" aggregate 0 0 "$OUTDIR"
  exit 0
fi
[[ "$ACTION" == "worker" ]] || { echo "Unknown EGCAR_ACTION=$ACTION" >&2; exit 2; }

TASK_OFFSET="${EGCAR_TASK_OFFSET:-0}"
local_task="${SLURM_ARRAY_TASK_ID:?Submit this script as a SLURM job array}"
[[ "$TASK_OFFSET" =~ ^[0-9]+$ && "$local_task" =~ ^[0-9]+$ ]] || exit 2
global_task=$(( 10#$TASK_OFFSET + 10#$local_task ))

expected=$(Rscript "$RSCRIPT" expected_tasks)
[[ "$expected" =~ ^[0-9]+$ ]] || { echo "Cannot read expected task count" >&2; exit 2; }

if (( global_task < 1 || global_task > expected )); then
  echo "Invalid global task $global_task; expected 1..$expected" >&2; exit 2
fi

# task_info now returns 6 tab-separated fields (signal added at the end).
read -r config_id rep_id n p_per_block rank signal < <(Rscript "$RSCRIPT" task_info "$global_task")

echo "JOB=${SLURM_JOB_ID:-NA} ARRAY_JOB=${SLURM_ARRAY_JOB_ID:-NA} LOCAL_TASK=$local_task"
echo "offset=$TASK_OFFSET global_task=$global_task config_id=$config_id rep_id=$rep_id"
echo "n=$n p_per_block=$p_per_block rank=$rank signal=$signal cpus=${SLURM_CPUS_PER_TASK:-NA} output=$OUTDIR"

if (( p_per_block >= 1000 )); then
  echo "WARNING: dense, full-dimensional reference CV. Memory/time at this p are unprofiled." >&2
fi
if (( p_per_block >= 5000 )); then
  echo "WARNING: p_total=15000. One dense double p_total-by-p_total matrix alone is 1.8 GB." >&2
fi

Rscript "$RSCRIPT" check_packages
Rscript "$RSCRIPT" worker "$config_id" "$rep_id" "$OUTDIR"

completed=$(find "$OUTDIR/completed" -maxdepth 1 -type f -name 'task_*.done' | wc -l | tr -d ' ')
echo "Committed tasks: $completed/$expected"

if [[ "$completed" -eq "$expected" && ! -f "$OUTDIR/AGGREGATION_COMPLETE.txt" ]]; then
  LOCK="$OUTDIR/.aggregation_lock"
  if mkdir "$LOCK" 2>/dev/null; then
    trap 'rmdir "$LOCK" 2>/dev/null || true' EXIT
    if [[ ! -f "$OUTDIR/AGGREGATION_COMPLETE.txt" ]]; then
      Rscript "$RSCRIPT" aggregate 0 0 "$OUTDIR"
    fi
  else
    echo "Another task is aggregating."
  fi
fi
