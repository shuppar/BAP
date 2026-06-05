# =============================================================================
# run_pipeline_WS.sh — WORKSTATION runbook (brain + placenta, real data)
# =============================================================================
# Not executable; a runbook to walk through commands one-by-one in tmux.
#
# Conventions match dev runbook (run_pipeline.sh):
#   - All scripts subcluster-aware via `--subcluster <slug>`
#   - 8c REQUIRES --tf (gates 8f view 5 TF concordance and 8g view 3)
#   - 8e/8f/8g write under per-phase subfolders for tables and plots
#
# Pre-flight: confirm GPU is available and dependencies are installed.
# (./setup-remote.sh, run once, has already installed: uv venv, R packages,
#  CellBender sidecar, refs/msigdb_mouse.tsv, refs/abc_brain_ref.h5ad,
#  refs/celltypist_brain_adult.pkl. If any of these are missing, re-run
#  ./setup-remote.sh — it's idempotent.)
nvidia-smi
which uv
which Rscript
test -f refs/msigdb_mouse.tsv          && echo "msigdb_mouse.tsv present"          || echo "MISSING — re-run ./setup-remote.sh"
test -f refs/abc_brain_ref.h5ad        && echo "abc_brain_ref.h5ad present"        || echo "MISSING — re-run ./setup-remote.sh"
test -f refs/celltypist_brain_adult.pkl && echo "celltypist_brain_adult.pkl present" || echo "MISSING — re-run ./setup-remote.sh"

# One-time per-machine bootstrap (if not already done):
#   ./setup-remote.sh
# This installs uv, runs `uv sync`, installs R + CRAN/Bioc packages, builds
# the CellBender sidecar venv (.venv-cellbender).


# =============================================================================
# Focal cell types — used by every loop below. Edit for placenta run.
# =============================================================================
# Brain default. Replace this list before the placenta run.
CELL_TYPES=(
  "Excitatory neurons"
  "Inhibitory neurons"
  "Microglia"
  "Oligodendrocytes"
  "Astrocytes"
  "OPC"
  # add "Radial glia / NPCs" for P1 brain
  # placenta: "Trophoblast" "Decidual" "EVT" "Endothelial" "Immune"
)

# Helper: slugify a cell type label the same way scripts do internally
slugify() { echo "$1" | tr '[:upper:] ' '[:lower:]_' | tr -d ',/'; }


# =============================================================================
# Pick the config — set this ONCE, then re-run the script for the other tissue.
# =============================================================================
CONFIG=config/brain.yaml          # change to config/placenta.yaml for placenta
# Sanity check the config picks up the right samples
grep -c "^- id:" "$CONFIG"        # should be 34 (brain) or 23 (placenta)


# =============================================================================
# Phase 0 — Validation (mandatory; 5 min, no compute)
# =============================================================================
# Catches sample swaps, metadata gaps, confounds before any GPU time.
uv run python scripts/01_validate.py --config "$CONFIG"
# Inspect: results/<tissue>/00_validation_report.html in browser


# =============================================================================
# Phase 1 — Ambient RNA (CellBender, GPU, ~1-2h/sample × parallel 2)
# =============================================================================
# Runs from the CellBender sidecar venv (separate from main env).
# Re-does cell-vs-empty calling and removes ambient contamination.
# Serialize GPU with scVI later: don't run them simultaneously.
CUDA_VISIBLE_DEVICES=0 .venv-cellbender/bin/python scripts/02_ambient.py \
    --config "$CONFIG"
# Inspect: results/<tissue>/plots/01_ambient/summary_ambient_fraction.png


# =============================================================================
# Phase 2 — QC (per-sample MAD + hard floors/caps + cohort flag)
# =============================================================================
uv run python scripts/02_qc.py --config "$CONFIG"
cat results/*/tables/02_qc/summary_qc.csv  # check no failed-prep flags


# =============================================================================
# Phase 3 — Doublet detection (scDblFinder per pool, R subprocess)
# =============================================================================
uv run python scripts/03_doublets.py --config "$CONFIG"


# =============================================================================
# Phase 4 — Concat + HVG + cell cycle
# =============================================================================
uv run python scripts/04_integration_prep.py --config "$CONFIG"


# =============================================================================
# Phase 5 — scVI integration (GPU, BF16, ~2-3h)
# =============================================================================
# Wrap in tmux so SSH drops don't kill it.
tmux new -s scvi -d "CUDA_VISIBLE_DEVICES=0 uv run python scripts/05_integration.py --config $CONFIG"
# Monitor:   tmux attach -t scvi    /    watch nvidia-smi
# After done: inspect plots/05_integration/post_integration_umap_by_pool.png


# =============================================================================
# Phase 6 — Clustering (Leiden, igraph backend)
# =============================================================================
uv run python scripts/06_clustering.py --config "$CONFIG"


# =============================================================================
# Phase 7 — Annotation (CellTypist + markers, per-cluster majority)
# =============================================================================
# Adult brain (4W, 3mo): no built-in CellTypist model — uses the ABC-trained
# .pkl pointed to by annotation.celltypist_models.{4W,3mo} in YAML.
uv run python scripts/07_annotation.py --config "$CONFIG"


# =============================================================================
# Phase 7b — Subclustering (loop over ALL focal cell types)
# =============================================================================
for ct in "${CELL_TYPES[@]}"; do
  uv run python scripts/07b_subcluster.py --config "$CONFIG" --celltype "$ct"
done


# =============================================================================
# Phase 7d — Subcluster annotation (CellTypist + marker scoring; cluster-level)
# =============================================================================
for ct in "${CELL_TYPES[@]}"; do
  uv run python scripts/07d_subcluster_annotate.py --config "$CONFIG" \
      --celltype "$ct" --markers config/subcluster_markers.yaml
done


# =============================================================================
# Phase 7c — scANVI reference label transfer (WORKSTATION ONLY)
# =============================================================================
# Requires reference: block in YAML + labeled reference h5ad. For brain:
# Allen Brain Cell Atlas build via prepare_reference.py. For placenta:
# Marsh & Blelloch 2020 staged similarly.
uv run python scripts/07c_label_transfer.py --config "$CONFIG"


# =============================================================================
# Phase 8a — Composition (propeller via R subprocess)
# =============================================================================
uv run python scripts/08a_composition.py --config "$CONFIG"


# =============================================================================
# Phase 8b — Pseudobulk DE (PyDESeq2). Main + subcluster per cell type.
# =============================================================================
uv run python scripts/08b_de.py --config "$CONFIG"
for ct in "${CELL_TYPES[@]}"; do
  slug=$(slugify "$ct")
  uv run python scripts/08b_de.py --config "$CONFIG" --subcluster "$slug"
done


# =============================================================================
# Phase 8c — Pathway/GSEA + TF activity. MAIN + SUBCLUSTER, --tf REQUIRED.
# =============================================================================
# --tf gates 8f view 5 (TF concordance) and 8g view 3 (TF persistence).
# Cannot be recovered without re-running 8c. ALWAYS pass --tf.
uv run python scripts/08c_pathways.py --config "$CONFIG" --tf
for ct in "${CELL_TYPES[@]}"; do
  slug=$(slugify "$ct")
  uv run python scripts/08c_pathways.py --config "$CONFIG" --tf --subcluster "$slug"
done


# =============================================================================
# Phase 8d — Trajectory (PAGA + DPT). NO velocity, NO CellRank (10x Flex).
# =============================================================================
uv run python scripts/08d_trajectory.py --config "$CONFIG"
# Optional: --root-celltype "Radial glia / NPCs" to set DPT root explicitly


# =============================================================================
# Phase 8e — Cell-cell communication. MAIN + SUBCLUSTER loop.
# =============================================================================
# Three arms in one script: baseline + differential + per-donor.
# --n-perms 1000 is the default (use lower only on dev).
uv run python scripts/08e_communication.py --config "$CONFIG" --zscore-rows
for ct in "${CELL_TYPES[@]}"; do
  slug=$(slugify "$ct")
  uv run python scripts/08e_communication.py --config "$CONFIG" \
      --subcluster "$slug" --zscore-rows
done


# =============================================================================
# CHECKPOINT — Run the above through Phase 8e for BOTH tissues before continuing.
# =============================================================================
# Phases 8f and 8g operate across tissues / across ages and need both 8b and
# 8c (with --tf) finished for the corresponding tissue(s).
#   - 8f cross-tissue: needs BOTH brain and placenta finished through 8c
#   - 8g cross-age:    brain only (placenta has incomplete cross-age factorial;
#                      script exits cleanly with a warning if tissue: placenta)


# =============================================================================
# Phase 8f — Cross-tissue (placenta → brain cascades). Six views.
# =============================================================================
# REAL run on the workstation: NO --dev-test, pass real configs.
uv run python scripts/08f_cross_tissue.py \
    --brain-config config/brain.yaml \
    --placenta-config config/placenta.yaml
# Key output:  results/brain/tables/08f_cross_tissue/08f_lr_cross_tissue.csv


# =============================================================================
# Phase 8g — Cross-age persistence (brain only). Six views.
# =============================================================================
uv run python scripts/08g_cross_age.py --config config/brain.yaml
# Key output:  results/brain/tables/08g_cross_age/08g_core_signature_genes.csv


# =============================================================================
# Done. Headline files for paper assembly:
# =============================================================================
#   results/brain/tables/08f_cross_tissue/08f_lr_cross_tissue.csv
#       Placental ligand × brain receptor mechanistic hypotheses, with
#       stress_axis column flagging GR/CRH/cytokine axis genes.
#
#   results/brain/tables/08g_cross_age/08g_core_signature_genes.csv
#       Genes persistent in BOTH Early and Late arms, same direction
#       throughout development. The most robust core stress signature.
#
#   results/brain/tables/08b_de/08b_de_results.csv
#       Master DE table; powers everything downstream.
#
#   results/brain/tables/08c_pathways/08c_pathway_results.csv
#   results/brain/tables/08c_pathways/08c_tf_activity.csv
#       Pathway and TF activity, joined to DE on (celltype, gene).
#
# Per-cell-type subcluster equivalents of 8b/8c/8e live in suffixed files:
#   08b_de_results_subcluster_<slug>.csv
#   08c_pathway_results_subcluster_<slug>.csv
#   tables/08e_communication_subcluster_<slug>/08e_*.csv
#
# Phase 9 figure assembly: deferred to a post-workstation notebook step.
# Phase 10 provenance: optional pre/parallel — manifest.json + provenance.py.
