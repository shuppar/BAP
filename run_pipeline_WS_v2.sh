# =============================================================================
# run_pipeline_WS_v2.sh — WORKSTATION runbook (brain + placenta, real data)
# =============================================================================
# Supersedes run_pipeline_WS.sh. UPDATED Friday June 5 2026 with:
#   - CellBender Phase 1 REMOVED (weakref bug, see STATUS_06-05_v2.md §3)
#   - Phase 5 UMAP seed sweep step added (5b)
#   - Brain Phase 7 uses per-age CellTypist models (P1+4W+3mo all mapped)
#   - Placenta Phase 7 uses curated Marsh+Simmons markers (placenta_markers.yaml
#     already merged into placenta.yaml)
#   - Phase 9 split into ARM A (psychiatric/neurodev) and ARM B (MS as
#     stressed-cell signature reference, NOT etiology)
#   - Smoke-test policy enforced for every >10 min phase
#
# Not executable; a runbook to walk through commands one-by-one in tmux.
#
# Conventions:
#   - All scripts subcluster-aware via `--subcluster <slug>`
#   - 8c REQUIRES --tf (gates 8f view 5 TF concordance and 8g view 3)
#   - 8e/8f/8g write under per-phase subfolders
#   - Anything >1 min → tmux. `<phase>_<tissue>` naming.
#   - python -u for live output in tee logs


# =============================================================================
# Pre-flight (refs already in place from earlier setup-remote.sh)
# =============================================================================
nvidia-smi
which uv
which Rscript
test -f refs/msigdb_mouse.tsv          && echo "msigdb_mouse.tsv present"          || echo "MISSING"
test -f refs/abc_brain_ref.h5ad        && echo "abc_brain_ref.h5ad present"        || echo "MISSING"
test -f refs/celltypist_brain_adult.pkl && echo "celltypist_brain_adult.pkl present" || echo "MISSING"

# NO LONGER NEEDED (CellBender skipped): .venv-cellbender/, ambient-RNA step.


# =============================================================================
# Focal cell types per tissue
# =============================================================================
# Brain — edit between brain and placenta runs
CELL_TYPES_BRAIN=(
  "Excitatory neurons"
  "Inhibitory neurons"
  "Microglia"
  "Oligodendrocytes"
  "Astrocytes"
  "OPC"
  # P1 brain: add "Radial glia / NPCs" after Phase 7 review
)
CELL_TYPES_PLACENTA=(
  "Syncytiotrophoblast-II (SynT-II)"
  "Spongiotrophoblast (SpT)"
  "Hofbauer cell (fetal macrophage)"
  "Fetal endothelium (labyrinth)"
  "Decidual stromal"
  # Fill after Phase 7 placenta review
)

slugify() { echo "$1" | tr '[:upper:] ' '[:lower:]_' | tr -d ',/'; }

# Set CONFIG once, re-run script for the other tissue
CONFIG=config/brain.yaml          # or config/placenta.yaml
grep -c "^- id:" "$CONFIG"        # 34 brain or 23 placenta


# =============================================================================
# Phase 0 — Validation (5 min, no compute)
# =============================================================================
tmux new -s phase0_brain -d \
  "uv run python scripts/01_validate.py --config config/brain.yaml \
   2>&1 | tee logs/phase0_brain.log"
tmux new -s phase0_placenta -d \
  "uv run python scripts/01_validate.py --config config/placenta.yaml \
   2>&1 | tee logs/phase0_placenta.log"
# Inspect: results/<tissue>/validation/validation_report.txt


# =============================================================================
# Phase 0b — Add assigned_sex column to sex_check.csv (one-time, ~10 sec)
# =============================================================================
# 10x Flex under-detects Xist, so most declared-F samples come back as
# 'ambiguous' in inferred_sex. assigned_sex resolves this: ambiguous → F
# (absence of Y = female by biology). C5 brain stays M (real Y mismatch).
# ALL downstream sex-stratified analyses (8a propeller, 8b DE design) must
# use assigned_sex, NOT declared_sex (has unknowns) or inferred_sex
# (has ambiguous nulls).
uv run python <<'PYEOF'
import pandas as pd
from pathlib import Path
for tissue in ['brain', 'placenta']:
    p = Path(f'results/{tissue}/validation/sex_check.csv')
    if not p.exists():
        print(f'[{tissue}] sex_check.csv missing — run Phase 0 first')
        continue
    df = pd.read_csv(p)
    df['assigned_sex'] = df['inferred_sex'].where(df['inferred_sex'] != 'ambiguous', 'F')
    df.to_csv(p, index=False)
    n_amb = (df['inferred_sex'] == 'ambiguous').sum()
    n_swap = ((df['declared_sex'] != df['assigned_sex']) &
              (df['declared_sex'] != 'unknown')).sum()
    print(f'[{tissue}] {len(df)} samples | {n_amb} ambiguous → F | '
          f'{n_swap} declared/assigned mismatches (review)')
PYEOF


# =============================================================================
# Phase 1 — Ambient RNA (CellBender) — SKIPPED 2026-06-05
# =============================================================================
# weakref.ReferenceType pickle bug in CellBender 0.3.0/0.3.2 across all
# torch+pyro+numpy combos. To revisit: official Docker image or SoupX.
# Phase 2 reads sample_filtered_feature_bc_matrix.h5 directly.


# =============================================================================
# Phase 2 — QC (per-sample MAD + hard caps)
# =============================================================================
# Smoke test on E6 (smallest brain sample) before launching full:
uv run python -c "
import yaml
cfg = yaml.safe_load(open('config/brain.yaml'))
cfg['samples'] = [s for s in cfg['samples'] if s['id']=='E6']
cfg['results_dir'] = 'results/brain_smoketest'
yaml.safe_dump(cfg, open('config/brain_smoketest.yaml','w'), sort_keys=False)
"
uv run python scripts/02_qc.py --config config/brain_smoketest.yaml
# If clean:
rm -rf results/brain_smoketest config/brain_smoketest.yaml

tmux new -s qc_brain -d \
  "uv run python -u scripts/02_qc.py --config config/brain.yaml \
   2>&1 | tee logs/qc_brain.log"
tmux new -s qc_placenta -d \
  "uv run python -u scripts/02_qc.py --config config/placenta.yaml \
   2>&1 | tee logs/qc_placenta.log"


# =============================================================================
# Phase 3 — Doublet detection (scDblFinder per pool, R subprocess)
# =============================================================================
# IMPORTANT: scDblFinder runs PER POOL. Pool3 has both brain (2 P1 Late) and
# placenta (14 E12.5) samples — running brain.yaml and placenta.yaml separately
# means scDblFinder doesn't see cross-tissue doublets. Documented limitation;
# in practice cross-tissue doublets have implausible expression and are caught
# by the artificial-doublet simulation anyway.
tmux new -s doublets_brain -d \
  "uv run python -u scripts/03_doublets.py --config config/brain.yaml \
   2>&1 | tee logs/doublets_brain.log"
tmux new -s doublets_placenta -d \
  "uv run python -u scripts/03_doublets.py --config config/placenta.yaml \
   2>&1 | tee logs/doublets_placenta.log"


# =============================================================================
# Phase 4 — Concat + HVG + cell cycle
# =============================================================================
# Hemo genes excluded from HVGs (critical for placenta where Hbb dominates PCs)
tmux new -s prep_brain -d \
  "uv run python -u scripts/04_integration_prep.py --config config/brain.yaml \
   2>&1 | tee logs/prep_brain.log"
tmux new -s prep_placenta -d \
  "uv run python -u scripts/04_integration_prep.py --config config/placenta.yaml \
   2>&1 | tee logs/prep_placenta.log"


# =============================================================================
# Phase 5 — scVI integration (GPU, BF16, ~30-45 min per tissue at 5 sec/epoch)
# =============================================================================
# Smoke test with max_epochs=5 on a 1-sample subset FIRST.
# Sequential, NOT parallel (single GPU). Brain first, then placenta.
# Cell cycle: did NOT drive global clustering in either tissue (verified
# from umap_post_phase.png). condition_cell_cycle stays OFF.

# Brain smoke test (3-5 min):
uv run python -c "
import yaml
cfg = yaml.safe_load(open('config/brain.yaml'))
cfg.setdefault('scvi', {})['max_epochs'] = 5
yaml.safe_dump(cfg, open('config/brain_scvi_smoketest.yaml','w'), sort_keys=False)
"
CUDA_VISIBLE_DEVICES=0 uv run python scripts/05_integration.py \
  --config config/brain_scvi_smoketest.yaml
# Verify outputs in results/brain/h5ad/06_integrated/, then:
rm -rf results/brain/h5ad/06_integrated/ \
       results/brain/plots/05_integration/ \
       results/brain/tables/05_integration/ \
       config/brain_scvi_smoketest.yaml

# Full brain (sequential — single GPU):
tmux new -s scvi_brain -d \
  "CUDA_VISIBLE_DEVICES=0 uv run python -u scripts/05_integration.py \
   --config config/brain.yaml 2>&1 | tee logs/scvi_brain.log"
# Wait for completion BEFORE launching placenta:
while ! grep -q 'Phase 5 complete' logs/scvi_brain.log 2>/dev/null; do sleep 60; done

tmux new -s scvi_placenta -d \
  "CUDA_VISIBLE_DEVICES=0 uv run python -u scripts/05_integration.py \
   --config config/placenta.yaml 2>&1 | tee logs/scvi_placenta.log"


# =============================================================================
# Phase 5b — UMAP seed sweep (NEW; CPU only, ~30 min per tissue)
# =============================================================================
# UMAP is non-deterministic without random_state. We don't retrain scVI.
# Loads integrated h5ad, recomputes neighbors + UMAP for 5 seeds, saves a
# 5-panel comparison plot per covariate. Pick the cleanest UMAP and lock it.
# (Script `05b_umap_sweep.py` to be written — TODO next session.)
for tissue in brain placenta; do
  tmux new -s umap_sweep_${tissue} -d \
    "uv run python -u scripts/05b_umap_sweep.py \
     --config config/${tissue}.yaml \
     --seeds 42,0,7,123,2024 \
     2>&1 | tee logs/umap_sweep_${tissue}.log"
done
# Review plots in results/<tissue>/plots/05b_umap_sweep/
# Choose a seed and edit config/<tissue>.yaml: scvi.umap_seed: <chosen>
# Then re-run 05b with --apply to overwrite .obsm['X_umap'] in the h5ad.


# =============================================================================
# Phase 6 — Clustering (Leiden, igraph backend)
# =============================================================================
# Multi-resolution sweep [0.2, 0.4, 0.6, 0.8, 1.0, 1.5, 2.0]; auto-pick by
# geometric knee. Override with --resolution if the auto-pick looks off.
# Defaults gave: brain 51 clusters @ res 0.4, placenta 29 clusters @ res 0.6.
tmux new -s cluster_brain -d \
  "uv run python -u scripts/06_clustering.py --config config/brain.yaml \
   2>&1 | tee logs/cluster_brain.log"
tmux new -s cluster_placenta -d \
  "uv run python -u scripts/06_clustering.py --config config/placenta.yaml \
   2>&1 | tee logs/cluster_placenta.log"


# =============================================================================
# Phase 7 — Annotation
# =============================================================================
# BRAIN — per-age CellTypist:
#   P1 → Developing_Mouse_Brain.pkl (built-in)
#   4W → refs/celltypist_brain_adult.pkl (ABC-trained, this session)
#   3mo → refs/celltypist_brain_adult.pkl
# Verify config/brain.yaml has all three mapped:
grep -A4 "celltypist_models:" config/brain.yaml
# Must show P1, 4W, 3mo. If only P1 (the original bug), fix:
uv run python -c "
import yaml
cfg = yaml.safe_load(open('config/brain.yaml'))
cfg['annotation'] = {'celltypist_models': {
    'P1': 'Developing_Mouse_Brain.pkl',
    '4W': 'refs/celltypist_brain_adult.pkl',
    '3mo': 'refs/celltypist_brain_adult.pkl'}}
yaml.safe_dump(cfg, open('config/brain.yaml','w'), sort_keys=False)
"

# PLACENTA — no CellTypist model. Uses curated markers (Marsh+Simmons+others)
# from placenta_markers.yaml, already merged into placenta.yaml under
# annotation.markers. Verify:
grep -c "annotation:" config/placenta.yaml      # ≥1
grep -c "markers:" config/placenta.yaml         # ≥1

tmux new -s annot_brain -d \
  "uv run python -u scripts/07_annotation.py --config config/brain.yaml \
   2>&1 | tee logs/annot_brain.log"
tmux new -s annot_placenta -d \
  "uv run python -u scripts/07_annotation.py --config config/placenta.yaml \
   2>&1 | tee logs/annot_placenta.log"
# REVIEW: cluster_purity.csv (any cluster <60% majority needs manual check)


# =============================================================================
# Phase 7b — Subclustering (loop over focal cell types)
# =============================================================================
for ct in "${CELL_TYPES_BRAIN[@]}"; do      # or PLACENTA
  uv run python scripts/07b_subcluster.py --config "$CONFIG" --celltype "$ct"
done


# =============================================================================
# Phase 7d — Subcluster annotation
# =============================================================================
for ct in "${CELL_TYPES_BRAIN[@]}"; do
  uv run python scripts/07d_subcluster_annotate.py --config "$CONFIG" \
      --celltype "$ct" --markers config/subcluster_markers.yaml
done


# =============================================================================
# Phase 7c — scANVI reference label transfer (WORKSTATION ONLY)
# =============================================================================
# SMOKE TEST FIRST with reduced epochs. Brain uses refs/abc_brain_ref.h5ad.
# Placenta has no ref — script exits cleanly.
uv run python scripts/07c_label_transfer.py --config "$CONFIG"


# =============================================================================
# Phase 8a — Composition (propeller via R)
# =============================================================================
uv run python scripts/08a_composition.py --config "$CONFIG"


# =============================================================================
# Phase 8b — Pseudobulk DE (PyDESeq2). Main + subcluster.
# =============================================================================
# NOTE: brain — also run a sensitivity DE excluding C5 (sex swap flag)
uv run python scripts/08b_de.py --config "$CONFIG"
for ct in "${CELL_TYPES_BRAIN[@]}"; do
  slug=$(slugify "$ct")
  uv run python scripts/08b_de.py --config "$CONFIG" --subcluster "$slug"
done


# =============================================================================
# Phase 8c — Pathway/GSEA + TF activity. MAIN + SUBCLUSTER. --tf REQUIRED.
# =============================================================================
# SMOKE TEST on one cell type first. --tf gates 8f view 5 and 8g view 3.
uv run python scripts/08c_pathways.py --config "$CONFIG" --tf
for ct in "${CELL_TYPES_BRAIN[@]}"; do
  slug=$(slugify "$ct")
  uv run python scripts/08c_pathways.py --config "$CONFIG" --tf --subcluster "$slug"
done


# =============================================================================
# Phase 8d — Trajectory (PAGA + DPT)
# =============================================================================
uv run python scripts/08d_trajectory.py --config "$CONFIG"


# =============================================================================
# Phase 8e — Cell-cell communication. MAIN + SUBCLUSTER.
# =============================================================================
uv run python scripts/08e_communication.py --config "$CONFIG" --zscore-rows
for ct in "${CELL_TYPES_BRAIN[@]}"; do
  slug=$(slugify "$ct")
  uv run python scripts/08e_communication.py --config "$CONFIG" \
      --subcluster "$slug" --zscore-rows
done


# =============================================================================
# CHECKPOINT — run 0-8e for BOTH tissues before continuing
# =============================================================================


# =============================================================================
# Phase 8f — Cross-tissue (placenta → brain). Six views.
# =============================================================================
uv run python scripts/08f_cross_tissue.py \
    --brain-config config/brain.yaml \
    --placenta-config config/placenta.yaml


# =============================================================================
# Phase 8g — Cross-age persistence (brain only). Six views.
# =============================================================================
uv run python scripts/08g_cross_age.py --config config/brain.yaml


# =============================================================================
# Phase 9 — Cross-species (TWO ARMS)
# =============================================================================
# ARM A — psychiatric/neurodev (Nagy + Maitra + Velmeshev + Herring + Marsh)
# ARM B — MS as stressed-cell signature ref (Macnair + Absinta + Jäkel)
#         "MS as REFERENCE for stressed microglia/OL states, NOT etiology"
#
# ARM A downloads complete. ARM B downloads NOT YET (run below):
tmux new -s ms_downloads -d \
  "bash scripts/download_ms_validation.sh \
   2>&1 | tee logs/ms_downloads.log"
# Wait for completion (~1-2 hours, network-bound), then:

# SMOKE TEST Phase 9 on Velmeshev first (most-complete ARM A dataset, simplest format)
uv run python scripts/09_cross_species_validation.py \
    --config config/brain.yaml \
    --celltype-map config/cross_species_celltype_map.yaml \
    --datasets velmeshev_2019_autism \
    --arm A

# Then full ARM A:
uv run python scripts/09_cross_species_validation.py \
    --config config/brain.yaml \
    --celltype-map config/cross_species_celltype_map.yaml \
    --arm A

# Then ARM B:
uv run python scripts/09_cross_species_validation.py \
    --config config/brain.yaml \
    --celltype-map config/cross_species_celltype_map.yaml \
    --arm B
# Outputs framed as SEPARATE scientific questions per arm.


# =============================================================================
# Headline outputs for paper
# =============================================================================
#   results/brain/tables/08f_cross_tissue/08f_lr_cross_tissue.csv
#   results/brain/tables/08g_cross_age/08g_core_signature_genes.csv
#   results/brain/tables/08b_de/08b_de_results.csv
#   results/brain/tables/08c_pathways/08c_pathway_results.csv
#   results/brain/tables/08c_pathways/08c_tf_activity.csv
#   results/brain/tables/09_cross_species_armA/*_rrho2.csv
#   results/brain/tables/09_cross_species_armB/*_rrho2.csv   ← MS arm
#   results/brain/tables/09_cross_species_armB/*_subset_rrho2.csv  ← MIMS-iron etc.
#
# Plot refactor (PNG+PDF hybrid) targeted at the 5-10 paper figures only,
# after Phase 8 results land. Don't refactor all 13 plot scripts.
