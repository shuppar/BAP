#!/usr/bin/env bash
# download_human_validation.sh — pull all open-access human snRNA-seq datasets
# for Phase 9 cross-species RRHO2 validation of the mouse prenatal-stress study.
#
# Idempotent: per-dataset directories skipped if non-empty.
# Run after the mouse pipeline reference build (ABC atlas) finishes.
#
# Datasets (all fully open, no dbGaP/Synapse application needed):
#   1. Nagy 2020 (GSE144136)    — adult dlPFC MDD M, 17/17, ~80K nuclei
#   2. Maitra 2023 (GSE213982)  — adult dlPFC MDD F+M, 71 donors, >160K nuclei
#   3. Velmeshev 2019           — ASD PFC+ACC, ages 4-22, 31 donors (UCSC)
#   4. Herring 2022 (GSE168408) — developmental PFC, gestation→adult, 17 donors
#   5. Marsh 2022 (GSE198373)   — mid-gestation placenta, 8 samples (matches E12)
#
# Total expected size: ~30-50 GB. Lands on NVMe (/home/poller/BAP-BrainPlacenta/).
#
# Stage 2 (controlled access) NOT covered here. User must apply via dbGaP /
# Synapse to their institutional signing official:
#   - ECHO-PATHWAYS placental bulk (phs003619 CANDLE + phs003620 GAPPS)
#   - PsychENCODE brainSCOPE (Synapse, SAGE)
#   - Hwang/Girgenti PTSD 2025 (PsychENCODE/Synapse)
#   - Pique-Regi/Garcia-Flores term placenta (phs001886)
#
# Usage:
#   ./scripts/download_human_validation.sh

set -euo pipefail

BASE=/home/poller/BAP-BrainPlacenta/data/human_validation

log()  { printf "\033[1;34m[hv-dl]\033[0m %s\n" "$*"; }
warn() { printf "\033[1;33m[warn ]\033[0m %s\n" "$*" >&2; }

# Idempotency helper: skip if dir has any files, else download into it.
maybe_download() {
  local dest="$1"; shift
  local label="$1"; shift
  # remaining args = URLs to download
  mkdir -p "$dest"
  if [ -n "$(ls -A "$dest" 2>/dev/null)" ]; then
    log "[$label] already populated ($dest) — skipping."
    return 0
  fi
  log "[$label] downloading to $dest ..."
  for url in "$@"; do
    wget -nv --show-progress -c -P "$dest" "$url" || {
      warn "[$label] FAILED on $url — partial files left in place; can resume with -c"
      return 1
    }
  done
  log "[$label] done."
}

# Mirror an entire GEO suppl/ directory (handles arbitrary file lists per series).
mirror_geo_suppl() {
  local gse="$1"          # e.g. GSE144136
  local dest="$2"
  mkdir -p "$dest"
  if [ -n "$(ls -A "$dest" 2>/dev/null)" ]; then
    log "[$gse] already populated ($dest) — skipping."
    return 0
  fi
  local stub="${gse:0:$(( ${#gse} - 3 ))}nnn"   # e.g. GSE144nnn
  local url="https://ftp.ncbi.nlm.nih.gov/geo/series/${stub}/${gse}/suppl/"
  log "[$gse] mirroring $url -> $dest"
  # -r recursive, -l1 one level, -np no parent, -nd no dir hierarchy,
  # -A '*' accept all, -c continue on retry.
  wget -nv --show-progress -r -l1 -np -nd -c -A '*' -P "$dest" "$url" || {
    warn "[$gse] mirror FAILED — partial files left for resume"
    return 1
  }
  # wget leaves an index.html behind from the directory listing — drop it
  rm -f "$dest/index.html" "$dest/robots.txt"
  log "[$gse] done."
}

# ============================================================================
# 1. Nagy 2020 — MDD dlPFC male (GSE144136)
# ============================================================================
mirror_geo_suppl GSE144136 "$BASE/brain/nagy_2020_GSE144136"

# ============================================================================
# 2. Maitra 2023 — MDD dlPFC female + reanalysed male (GSE213982)
# ============================================================================
mirror_geo_suppl GSE213982 "$BASE/brain/maitra_2023_GSE213982"

# ============================================================================
# 3. Velmeshev 2019 — ASD PFC + ACC (UCSC Cell Browser autism dataset)
# ============================================================================
# Raw 10X MTX format (matrix.mtx + barcodes + genes + meta) + log2 expr matrix
maybe_download "$BASE/brain/velmeshev_2019_autism" "velmeshev_2019" \
  "https://cells.ucsc.edu/autism/rawMatrix.zip" \
  "https://cells.ucsc.edu/autism/exprMatrix.tsv.gz" \
  "https://cells.ucsc.edu/autism/meta.tsv"

# ============================================================================
# 4. Herring 2022 — developmental PFC (GSE168408)
# ============================================================================
# This is the largest in the set (snRNA + snATAC). Expect 10-20 GB.
mirror_geo_suppl GSE168408 "$BASE/brain/herring_2022_GSE168408"

# ============================================================================
# 5. Marsh 2022 — mid-gestation placenta (GSE198373)
# ============================================================================
mirror_geo_suppl GSE198373 "$BASE/placenta/marsh_2022_GSE198373"

# ============================================================================
# Summary
# ============================================================================
echo
log "=========================================="
log "Human validation downloads complete."
log "=========================================="
echo
log "Disk usage:"
du -sh "$BASE"/brain/* "$BASE"/placenta/* 2>/dev/null
echo
log "Next:"
log "  1. Inspect file lists per dataset; some come as MTX/barcodes/features,"
log "     others as .h5 or .csv. Phase 9 RRHO2 script will handle each format."
log "  2. Apply for Stage-2 controlled-access data via your institutional"
log "     signing official (dbGaP eRA Commons / Synapse SAGE)."
log "     See: refs/dbgap_application_checklist.md"
