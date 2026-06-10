#!/usr/bin/env bash
# download_ms_validation.sh — pull open-access human MS snRNA-seq datasets for
# Phase 9 cross-species validation (MS arm).
#
# Scientific framing (per Phase 9 strategy): MS is used as a CELL-TYPE
# SIGNATURE REFERENCE for stressed microglia and oligodendrocyte states. NOT
# as etiologic causation. Convergent disease-associated microglia
# (MIMS-foamy/iron ↔ DAM/stressed microglia) and convergent OL stress
# responses (Jäkel Oligo5/Oligo6 ↔ mouse prenatal-stress OPC/COP signatures)
# are the comparisons of interest. See STATUS_06-05.md §8 for full framing.
#
# Datasets (3 open Stage-1; Schirmer 2019 deferred — raw FASTQ only on SRA):
#   1. Macnair 2025 (Zenodo 10.5281/zenodo.8338963) — 632K nuclei, 156 samples,
#      54 MS + 28 controls, cleaned annotated counts matrices
#      Neuron 113(3):396-410, https://doi.org/10.1016/j.neuron.2024.11.016
#   2. Absinta 2021 (GSE180759) — MIMS-iron/MIMS-foamy microglia
#      Nature 597:709-714, https://doi.org/10.1038/s41586-021-03892-7
#   3. Jäkel 2019 (GSE118257 processed) — Oligo5/Oligo6 disease-associated OL
#      Nature 566:543-547, https://doi.org/10.1038/s41586-019-0903-2
#
# Idempotent: per-dataset directories skipped if non-empty.
#
# Usage:
#   bash scripts/download_ms_validation.sh

set -euo pipefail

BASE=/home/poller/BAP-BrainPlacenta/data/human_validation_ms
UA="Mozilla/5.0"

log()  { printf "\033[1;34m[ms-dl]\033[0m %s\n" "$*"; }
warn() { printf "\033[1;33m[warn ]\033[0m %s\n" "$*" >&2; }

dir_has_files() {
    local d="$1"
    [[ -d "$d" ]] && [[ "$(find "$d" -maxdepth 1 -type f 2>/dev/null | wc -l)" -gt 0 ]]
}

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

# GEO suppl directory mirror (same pattern as download_human_validation.sh).
# Excludes Cell Ranger reference bundles (premrna/genome) that some series ship.
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
        local file_url="${index_url}${f}"
        local file_dest="${dest}/${f}"
        fetch_file "$file_url" "$file_dest" || warn "  continuing despite failure on $f"
    done <<< "$files"
    log "[$gse] done."
}

# Zenodo record fetcher — uses the public REST API to list files, then
# downloads each via the API's direct download links.
fetch_zenodo_record() {
    local record_id="$1"
    local dest="$2"
    mkdir -p "$dest"
    if dir_has_files "$dest"; then
        log "[zenodo:$record_id] already populated — skipping."
        return 0
    fi
    local api="https://zenodo.org/api/records/${record_id}"
    log "[zenodo:$record_id] querying API: $api"
    local meta
    meta=$(curl -sS -A "$UA" "$api") || {
        warn "[zenodo:$record_id] API call failed"
        return 1
    }
    # Parse the JSON to extract file URLs + names. python3 is available system-wide.
    local file_list
    file_list=$(echo "$meta" | python3 -c "
import json, sys
data = json.load(sys.stdin)
for f in data.get('files', []):
    name = f.get('key', '')
    url = f.get('links', {}).get('self', '')
    if name and url:
        print(f'{name}\t{url}')
")
    if [[ -z "$file_list" ]]; then
        warn "[zenodo:$record_id] no files in record"
        return 1
    fi
    local n=$(echo "$file_list" | wc -l)
    log "[zenodo:$record_id] $n files to fetch"
    while IFS=$'\t' read -r name url; do
        [[ -z "$name" ]] && continue
        fetch_file "$url" "${dest}/${name}" || warn "  continuing despite failure on $name"
    done <<< "$file_list"
    log "[zenodo:$record_id] done."
}

# ============================================================================
# 1. Macnair 2025 (Zenodo 10.5281/zenodo.8338963)
# ============================================================================
# Contains: cleaned annotated counts matrices + sample metadata + cell type
# annotations for both cohorts. Likely Seurat .rds or h5ad files; Phase 9
# loader handles either.
fetch_zenodo_record 8338963 "$BASE/macnair_2025_zenodo8338963"

# ============================================================================
# 2. Absinta 2021 (GSE180759) — MIMS-iron/foamy microglia
# ============================================================================
fetch_geo_suppl GSE180759 "$BASE/absinta_2021_GSE180759"

# ============================================================================
# 3. Jäkel 2019 (GSE118257) — Oligo5/Oligo6 disease-associated OL
# ============================================================================
fetch_geo_suppl GSE118257 "$BASE/jakel_2019_GSE118257"

# ============================================================================
# Summary
# ============================================================================
echo
log "=========================================="
log "MS validation downloads complete."
log "=========================================="
echo
log "Disk usage per dataset:"
du -sh "$BASE"/* 2>/dev/null
echo
log "NOT included in this download (intentionally):"
log "  - Schirmer 2019 (SRA PRJNA544731): raw FASTQ only, requires Cell Ranger"
log "    reprocessing. Defer unless reviewers specifically request it."
log "  - Mouse OL trajectory anchors (Marques 2016 GSE75330, Falcão 2018"
log "    GSE113973): would be downloaded separately if added to Phase 8d."
echo
log "Next: update config/cross_species_celltype_map.yaml with the MS subtype"
log "labels (MIMS_iron, MIMS_foamy, Oligo5, Oligo6, Micro_A..Micro_E)."
log "See STATUS_06-05.md §8 for full scientific framing of the MS arm."
