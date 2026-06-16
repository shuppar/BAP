#!/usr/bin/env bash
# run_07b_all.sh — Phase 7b subcluster queue (all GPU, sequential, one tmux wrap).
#
# Subcluster ONLY lumped types where data-driven STATE DISCOVERY pays off.
# Neuron subtypes (deep-layer Ex, PV/SST/VIP, DA) are already resolved at the
# subclass tier → handled by focal DE in 8b, NOT here. Trophoblast subtypes are
# likewise already resolved (SpT/SynTI/II/S-TGC/GC/JZP/LaTP) → focal DE in 8b.
#
#   Brain    (celltypist_broad):  Immune, OPC/Oligodendrocytes, Astrocytes/Ependymal
#   Placenta (celltype_majority): DSC, Endothelium, Myeloid, NK
#
# All run GPU, SEQUENTIAL, in one tmux session. The GPU is serialized for VRAM,
# so one queue is safe; and since placenta scVI on CPU is slow, a single GPU
# queue beats the old GPU/CPU split. ~1.5–2 h total on the RTX 4500 Ada.
# set -e → stops on first failure (don't burn GPU hours on a broken run).
#
# Launch (from WS project root), wrapped in tmux per project convention:
#   tmux new -d -s 07b 'bash scripts/run_07b_all.sh 2>&1 | tee logs/07b_queue.log'
#   tmux attach -t 07b           # watch live
#   tail -f logs/07b_queue.log   # or follow the master log
#   tail -f logs/07b_immune.log  # per-type log
#
# Smoke-first (run ONE type, eyeball outputs, then launch the rest):
#   bash scripts/run_07b_all.sh Immune

set -euo pipefail
cd /home/poller/BAP-BrainPlacenta
mkdir -p logs

PYTHON="uv run python -u scripts/07b_subcluster.py"

# entry = "celltype:label_key:config:resolution"   (no celltype contains a space or ':')
# Resolution notes:
#   Immune 0.4              — microglia states + BAM/perivascular mac
#   OPC/Oligodendrocytes 0.5 — OPC / iOL / mOL axis + stressed OL
#   Astrocytes/Ependymal 0.4 — astro states + ependymal/choroid split
#   DSC 0.4                  — homeostatic vs stressed decidual stroma
#   Endothelium 0.3         — fetal vs maternal EC
#   Myeloid 0.4             — Hofbauer M1/M2/transitional
#   NK 0.4                  — dNK1/2/3
JOBS=(
    "Immune:celltypist_broad:config/brain.yaml:0.4"
    "OPC/Oligodendrocytes:celltypist_broad:config/brain.yaml:0.5"
    "Astrocytes/Ependymal:celltypist_broad:config/brain.yaml:0.3"
    "DSC:celltype_majority:config/placenta.yaml:0.4"
    "Endothelium:celltype_majority:config/placenta.yaml:0.3"
    "Myeloid:celltype_majority:config/placenta.yaml:0.4"
    "NK:celltype_majority:config/placenta.yaml:0.4"
)

# Optional CLI args = subset of celltypes to run (smoke-test one before the full queue)
SELECT=("$@")
run_this() {
    [ ${#SELECT[@]} -eq 0 ] && return 0
    local ct="$1"; local s
    for s in "${SELECT[@]}"; do [ "$s" = "$ct" ] && return 0; done
    return 1
}

for entry in "${JOBS[@]}"; do
    IFS=':' read -r CT LK CFG RES <<< "$entry"
    run_this "$CT" || { echo "--- skip ${CT}"; continue; }
    SLUG=$(echo "$CT" | tr '[:upper:]' '[:lower:]' \
        | sed 's/[^0-9a-zA-Z]/_/g; s/__*/_/g; s/^_//; s/_$//')
    LOG="logs/07b_${SLUG}.log"
    echo ">>> START ${CT}  (label=${LK}, res=${RES}, log=${LOG})"
    ${PYTHON} --config "${CFG}" --celltype "${CT}" \
        --label-key "${LK}" --resolution "${RES}" 2>&1 | tee "${LOG}"
    echo ">>> DONE ${CT}"
done

echo "=== ALL 7b JOBS DONE ==="
