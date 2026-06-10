#!/usr/bin/env bash
# download_human_validation.sh — pull all open-access human snRNA-seq datasets
# for Phase 9 cross-species RRHO2 validation.
#
# GEO FTP blocks wget via robots.txt, so we use curl + User-Agent to scrape the
# suppl/ directory listing, then download each file individually.
#
# Datasets (all fully open, no dbGaP/Synapse application needed):
#   1. Nagy 2020 (GSE144136)    — adult dlPFC MDD M
#   2. Maitra 2023 (GSE213982)  — adult dlPFC MDD F+M
#   3. Velmeshev 2019           — ASD PFC+ACC (UCSC)
#   4. Herring 2022 (GSE168408) — developmental PFC
#   5. Marsh 2022 (GSE198373)   — mid-gestation placenta
#
# Idempotent: per-dataset directories skipped if non-empty (has files, not just
# dir entries).
#
# Usage:
#   bash scripts/download_human_validation.sh

set -euo pipefail

BASE=/home/poller/BAP-BrainPlacenta/data/human_validation
UA="Mozilla/5.0"

log()  { printf "\033[1;34m[hv-dl]\033[0m %s\n" "$*"; }
warn() { printf "\033[1;33m[warn ]\033[0m %s\n" "$*" >&2; }

# Check if directory has any real files (not just metadata)
dir_has_files() {
    local d="$1"
    [[ -d "$d" ]] && [[ "$(find "$d" -maxdepth 1 -type f 2>/dev/null | wc -l)" -gt 0 ]]
}

# Download a single file via curl (resumable, with User-Agent)
fetch_file() {
    local url="$1"
    local dest="$2"
    if [[ -f "$dest" ]] && [[ -s "$dest" ]]; then
        log "  already have: $(basename "$dest")"
        return 0
    fi
    log "  fetching: $(basename "$dest")"
    curl -sS -A "$UA" -L -C - -o "$dest" "$url" || {
        warn "  FAILED: $url"
        return 1
    }
}

# Scrape GEO suppl directory and download all listed files
fetch_geo_suppl() {
    local gse="$1"
    local dest="$2"
    mkdir -p "$dest"
    if dir_has_files "$dest"; then
        log "[$gse] already populated — skipping."
        return 0
    fi
    local stub="${gse:0:$(( ${#gse} - 3 ))}nnn"
    local index_url="https://ftp.ncbi.nlm.nih.gov/geo/series/${stub}/${gse}/suppl/"
    log "[$gse] scraping index: $index_url"
    # Get the list of files (anything ending in a known data extension)
    local files
    files=$(curl -sS -A "$UA" "$index_url" | \
        grep -oE 'href="[^"]+\.(gz|tar|h5|csv|mtx|tsv|zip|xlsx|rds|rda|h5ad)"' | \
        sed 's/href="//;s/"$//' | \
        grep -viE 'premrna|reference|genome|cellranger|GRCh|GRCm|annotation' | \
        sort -u)
    if [[ -z "$files" ]]; then
        warn "[$gse] no files found at $index_url"
        return 1
    fi
    log "[$gse] $(echo "$files" | wc -l) files to fetch"
    local f
    while IFS= read -r f; do
        [[ -z "$f" ]] && continue
        # f is a relative filename like "GSE144136_X.csv.gz"
        local file_url="${index_url}${f}"
        local file_dest="${dest}/${f}"
        fetch_file "$file_url" "$file_dest" || warn "  continuing despite failure on $f"
    done <<< "$files"
    log "[$gse] done."
}

# ============================================================================
# 1. Nagy 2020 (GSE144136) — MDD dlPFC male
# ============================================================================
fetch_geo_suppl GSE144136 "$BASE/brain/nagy_2020_GSE144136"

# ============================================================================
# 2. Maitra 2023 (GSE213982) — MDD dlPFC female + reanalysed male
# ============================================================================
fetch_geo_suppl GSE213982 "$BASE/brain/maitra_2023_GSE213982"

# ============================================================================
# 3. Velmeshev 2019 — ASD PFC + ACC (UCSC Cell Browser, not GEO)
# ============================================================================
velmeshev_dest="$BASE/brain/velmeshev_2019_autism"
mkdir -p "$velmeshev_dest"
if dir_has_files "$velmeshev_dest"; then
    log "[velmeshev_2019] already populated — skipping."
else
    log "[velmeshev_2019] downloading from UCSC..."
    fetch_file "https://cells.ucsc.edu/autism/rawMatrix.zip"    "$velmeshev_dest/rawMatrix.zip"
    fetch_file "https://cells.ucsc.edu/autism/exprMatrix.tsv.gz" "$velmeshev_dest/exprMatrix.tsv.gz"
    fetch_file "https://cells.ucsc.edu/autism/meta.tsv"          "$velmeshev_dest/meta.tsv"
fi

# ============================================================================
# 4. Herring 2022 (GSE168408) — developmental PFC
# ============================================================================
fetch_geo_suppl GSE168408 "$BASE/brain/herring_2022_GSE168408"

# ============================================================================
# 5. Marsh 2022 (GSE198373) — mid-gestation placenta
# ============================================================================
fetch_geo_suppl GSE198373 "$BASE/placenta/marsh_2022_GSE198373"

# ============================================================================
# Summary
# ============================================================================
echo
log "=========================================="
log "Human validation downloads complete."
log "=========================================="
echo
log "Disk usage per dataset:"
du -sh "$BASE"/brain/* "$BASE"/placenta/* 2>/dev/null
echo
log "Next steps:"
log "  1. Inspect file lists per dataset; some come as MTX/barcodes/features,"
log "     others as .h5 or .csv. Phase 9 RRHO2 script will handle each format."
log "  2. Apply for Stage-2 controlled-access data via your institutional"
log "     signing official (dbGaP eRA Commons / Synapse SAGE)."
log "     See: refs/dbgap_application_checklist.md"
