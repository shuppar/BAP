#!/usr/bin/env bash
# queue_human_downloads.sh — wait for the ABC atlas tmux session to finish,
# then trigger the human-validation download script.
#
# Idempotent: if abc_ref is already gone, runs downloads immediately.
# If the reference outputs don't exist after abc_ref ends, ABC failed —
# we still run downloads but log a warning.
#
# Usage (workstation):
#   tmux new -s human_dl
#   ./scripts/queue_human_downloads.sh 2>&1 | tee logs/human_downloads.log
#   # Ctrl-b d to detach

set -euo pipefail

PROJ=/home/poller/BAP-BrainPlacenta
ABC_SESSION="abc_ref"
SLEEP_SEC=120

log() { printf "\033[1;34m[queue]\033[0m %s\n" "$*"; }

# -- Wait for abc_ref tmux session to end --
if tmux ls 2>/dev/null | grep -q "^${ABC_SESSION}:"; then
  log "Waiting for tmux session '${ABC_SESSION}' to end (polling every ${SLEEP_SEC}s)..."
  while tmux ls 2>/dev/null | grep -q "^${ABC_SESSION}:"; do
    sleep "$SLEEP_SEC"
  done
  log "tmux session '${ABC_SESSION}' has ended."
else
  log "No '${ABC_SESSION}' tmux session running — proceeding immediately."
fi

# -- Sanity-check ABC outputs (informational, not fatal) --
if [[ -f "$PROJ/refs/abc_brain_ref.h5ad" && -f "$PROJ/refs/celltypist_brain_adult.pkl" ]]; then
  log "ABC reference build outputs present — looks like it succeeded."
else
  log "WARN: ABC reference build outputs MISSING — it may have errored."
  log "      refs/abc_brain_ref.h5ad        : $(test -f $PROJ/refs/abc_brain_ref.h5ad && echo present || echo missing)"
  log "      refs/celltypist_brain_adult.pkl: $(test -f $PROJ/refs/celltypist_brain_adult.pkl && echo present || echo missing)"
  log "      Downloads will still proceed (they're independent)."
fi

# -- Kick off the downloads --
cd "$PROJ"
log "Starting human-validation downloads..."
./scripts/download_human_validation.sh
log "Queue runner finished."
