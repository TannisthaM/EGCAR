#!/bin/bash
#SBATCH --job-name=egcar_vp_t2
#SBATCH --output=logs/egcar_vp_t2_%A_%a.out
#SBATCH --error=logs/egcar_vp_t2_%A_%a.err
#SBATCH --array=1-50%10
#SBATCH --time=24:00:00
#SBATCH --partition=caslake
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=40G
#SBATCH --account=pi-cdonnat

# 50 tasks total (n=150 fixed, rank=1:5, 10 replications). Still one
# array in one sbatch call -- no offset-chunking needed at this size.
#
# p1=p2=p3=500, p_total=1500 -- this OOM'd at 20G in the original combined
# tier1 even at low rank; 40G is a real increase, watch MaxRSS on task 1.
# All CV methods share one fold-worker pool per task:
# min(SLURM_CPUS_PER_TASK - 1, 5 folds). Default here: 5 workers.
# Packages must be installed once BEFORE submission -- run
# `Rscript ... install_packages` with EGCAR_LOCAL_SOURCE set to your fixed
# egcar_0_2_14_l21_calibrated.zip first if not already installed
# at >=0.2.14 (this script's check_packages step enforces that version).

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

RSCRIPT="$ROOT/r/experiments/cluster/run_egcar_vary_p_tier2.R"
[[ -f "$RSCRIPT" ]] || { echo "Missing R script: $RSCRIPT" >&2; exit 2; }

RUN_ID="${EGCAR_RUN_ID:-job_${SLURM_ARRAY_JOB_ID:-${SLURM_JOB_ID:-manual}}}"
[[ "$RUN_ID" =~ ^[A-Za-z0-9][A-Za-z0-9_.-]*$ ]] || { echo "Invalid EGCAR_RUN_ID" >&2; exit 2; }

OUTDIR="$ROOT/egcar_vary_p_tier2_outputs/$RUN_ID"
mkdir -p "$ROOT/logs" "$OUTDIR"/{metrics,raw,completed,plots_all_methods,plots_without_oracle1}

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

read -r config_id rep_id n p_per_block rank signal < <(Rscript "$RSCRIPT" task_info "$global_task")

echo "JOB=${SLURM_JOB_ID:-NA} ARRAY_JOB=${SLURM_ARRAY_JOB_ID:-NA} LOCAL_TASK=$local_task"
echo "offset=$TASK_OFFSET global_task=$global_task config_id=$config_id rep_id=$rep_id"
echo "n=$n p_per_block=$p_per_block rank=$rank signal=$signal cpus=${SLURM_CPUS_PER_TASK:-NA} output=$OUTDIR"

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
