#!/usr/bin/env bash
# run_08c_remaining.sh — fire all remaining 8c runs sequentially.
#
# Covers (in order):
#   1. brain subclusters: immune, opc_oligodendrocytes, astrocytes_ependymal
#   2. placenta main
#   3. placenta subclusters: dsc, endothelium, myeloid, nk
#
# (Brain main + brain-main per-cell already complete; skipped here.)
#
# Pre-flight checks ALL inputs before any run; continue past per-run failures;
# per-run logs under logs/08c_remaining_<timestamp>/; status TSV with rc per
# run; final summary.
#
# Usage from tmux (recommended):
#   tmux new -d -s 08c_remaining 'bash run_08c_remaining.sh \
#     2>&1 | tee logs/08c_remaining_batch.log'
#   tail -f logs/08c_remaining_batch.log

set -u
trap 'echo "Interrupted"; exit 130' INT TERM

DRY_RUN=0
if [[ "${1:-}" == "--dry-run" ]]; then
  DRY_RUN=1
fi

REPO=/home/poller/BAP-BrainPlacenta
cd "$REPO" || { echo "ERROR: cannot cd to $REPO"; exit 1; }

N_JOBS=16
MIN_FREE_GB=30

# Job spec: "<tissue>:<subcluster_slug>". Empty subcluster = main run.
# Edit this array to skip/reorder runs.
JOBS=(
  "brain:immune"
  "brain:opc_oligodendrocytes"
  "brain:astrocytes_ependymal"
  "placenta:"                       # placenta main
  "placenta:dsc"
  "placenta:endothelium"
  "placenta:myeloid"
  "placenta:nk"
)

LOGDIR="logs/08c_remaining_$(date +%Y%m%d_%H%M%S)"
mkdir -p "$LOGDIR"
STATUS="$LOGDIR/status.tsv"
printf "run\tstart\tend\telapsed_s\trc\n" > "$STATUS"

echo "==============================================="
echo "08c batch: brain subclusters + placenta (main + subclusters)"
if [[ "$DRY_RUN" == "1" ]]; then
  echo "  *** DRY-RUN MODE: pre-flight only, no actual 8c runs ***"
fi
echo "  Jobs (${#JOBS[@]}):"
for J in "${JOBS[@]}"; do echo "    - $J"; done
echo "  --n-jobs:   $N_JOBS"
echo "  Log dir:    $LOGDIR"
echo "  Status TSV: $STATUS"
echo "==============================================="

# ---- Pre-flight ---------------------------------------------------------
echo
echo "[pre-flight]"

[[ -f scripts/08c_pathways.py ]] || { echo "  ERROR: scripts/08c_pathways.py missing"; exit 1; }
echo "  ✓ scripts/08c_pathways.py present"

FREE_GB=$(df --output=avail -BG "$REPO" | tail -1 | tr -dc '0-9')
if [[ "$FREE_GB" -lt "$MIN_FREE_GB" ]]; then
  echo "  ERROR: only ${FREE_GB}G free in $REPO (need ${MIN_FREE_GB}G)"
  exit 1
fi
echo "  ✓ ${FREE_GB}G free disk"

check_job_inputs() {
  local tissue=$1 sub=$2
  local config="config/${tissue}.yaml"
  [[ -f "$config" ]] || { echo "  ERROR: $config missing"; return 1; }

  local de h5 obs_needed
  if [[ -n "$sub" ]]; then
    de="results/${tissue}/tables/08b_de/08b_de_results_subcluster_${sub}.csv"
    h5="results/${tissue}/h5ad/08c_subclustered/${sub}.h5ad"
    obs_needed="donor_id sample_id sex age group pool subcluster_name"
  else
    de="results/${tissue}/tables/08b_de/08b_de_results.csv"
    h5="results/${tissue}/h5ad/08_annotated/all_samples.h5ad"
    case "$tissue" in
      brain)    obs_needed="donor_id sample_id sex age group pool celltypist_broad" ;;
      placenta) obs_needed="donor_id sample_id sex age group pool celltype_majority" ;;
    esac
  fi
  [[ -f "$de" ]] || { echo "  ERROR: $de missing"; return 1; }
  [[ -f "$h5" ]] || { echo "  ERROR: $h5 missing"; return 1; }

  uv run python -c "
import anndata as ad, sys
a = ad.read_h5ad('$h5', backed='r')
need = '$obs_needed'.split()
miss = [c for c in need if c not in a.obs.columns]
a.file.close()
if miss:
    sys.exit(f'missing obs cols: {miss}')
" || { echo "  ERROR: $h5 obs sanity failed"; return 1; }
  return 0
}

PRE_FAIL=0
for J in "${JOBS[@]}"; do
  TISSUE=${J%%:*}
  SUB=${J#*:}
  if check_job_inputs "$TISSUE" "$SUB"; then
    if [[ -n "$SUB" ]]; then
      echo "  ✓ $TISSUE / $SUB: inputs OK"
    else
      echo "  ✓ $TISSUE / main: inputs OK"
    fi
  else
    PRE_FAIL=$((PRE_FAIL + 1))
  fi
done

if [[ "$PRE_FAIL" -gt 0 ]]; then
  echo
  echo "ABORT: $PRE_FAIL pre-flight check(s) failed. Fix inputs before rerunning."
  exit 1
fi

# ---- Sequential runs ----------------------------------------------------
echo
echo "[runs]"

run_one() {
  local tissue=$1 sub=$2
  local label
  if [[ -n "$sub" ]]; then
    label="${tissue}_sub_${sub}"
  else
    label="${tissue}_main"
  fi
  local config="config/${tissue}.yaml"
  local log="$LOGDIR/${label}.log"
  local start=$(date +%s)

  echo
  echo "=== [$(date '+%H:%M:%S')] starting: $label"
  echo "    log: $log"

  local cmd_args=(--config "$config" --n-jobs "$N_JOBS")
  if [[ -n "$sub" ]]; then
    cmd_args+=(--subcluster "$sub")
  fi

  if [[ "$DRY_RUN" == "1" ]]; then
    echo "    [dry-run] would invoke: uv run python -u scripts/08c_pathways.py ${cmd_args[*]}"
    printf "%s\t%s\t%s\t%d\t%d\n" "$label" "$start" "$start" 0 0 >> "$STATUS"
    echo "+++ [$(date '+%H:%M:%S')] DRY-RUN OK: $label"
    return
  fi

  uv run python -u scripts/08c_pathways.py "${cmd_args[@]}" > "$log" 2>&1
  local rc=$?
  local end=$(date +%s)
  local elapsed=$((end - start))
  printf "%s\t%s\t%s\t%d\t%d\n" "$label" "$start" "$end" "$elapsed" "$rc" >> "$STATUS"

  if [[ $rc -ne 0 ]]; then
    echo "!!! [$(date '+%H:%M:%S')] FAILED: $label (rc=$rc, ${elapsed}s) — tail of log:"
    tail -n 30 "$log" | sed 's/^/    | /'
  else
    echo "+++ [$(date '+%H:%M:%S')] OK: $label (${elapsed}s) — tail of log:"
    tail -n 8 "$log" | sed 's/^/    | /'
  fi
}

for J in "${JOBS[@]}"; do
  TISSUE=${J%%:*}
  SUB=${J#*:}
  run_one "$TISSUE" "$SUB"
done

# ---- Final summary ------------------------------------------------------
echo
echo "==============================================="
echo "[final status]"
column -t -s $'\t' "$STATUS" || cat "$STATUS"
NFAIL=$(awk 'NR>1 && $5 != 0' "$STATUS" | wc -l)
NRUN=$(awk 'NR>1' "$STATUS" | wc -l)
echo
if [[ "$NFAIL" -eq 0 ]]; then
  echo "✓ All ${NRUN} runs completed successfully."
  echo "  Per-run logs: $LOGDIR/*.log"
else
  echo "✗ ${NFAIL}/${NRUN} runs failed. See logs in $LOGDIR."
  echo "  To resume a failed subcluster (skips AUCell if h5ad was written):"
  echo "    uv run python scripts/08c_pathways.py --config config/<tissue>.yaml \\"
  echo "      --subcluster <slug> --per-cell-only --n-jobs $N_JOBS"
fi
echo "==============================================="
exit "$NFAIL"
