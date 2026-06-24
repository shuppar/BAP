# =============================================================================
# run_pipeline_WS.sh — WORKSTATION runbook (brain + placenta, real data)
# =============================================================================
# UPDATED 2026-06-05: CellBender step removed. See STATUS_06-05.md §3 for
# rationale. Pipeline reads Cell Ranger filtered h5 directly from each
# sample's `h5` field in the YAML.
#
# Not executable; a runbook to walk through commands one-by-one in tmux.
#
# Conventions:
#   - All scripts subcluster-aware via `--subcluster <slug>`
#   - 8c REQUIRES --tf (gates 8f view 5 TF concordance and 8g view 3)
#   - 8e/8f/8g write under per-phase subfolders for tables and plots
#   - Anything >1 min → tmux. Session naming: <phase>_<tissue>.
#
# Pre-flight: confirm GPU is available and dependencies are installed.
# (./setup-remote.sh, run once, has already installed: uv venv, R packages,
#  refs/msigdb_mouse.tsv, refs/abc_brain_ref.h5ad,
#  refs/celltypist_brain_adult.pkl. If any of these are missing, re-run
#  ./setup-remote.sh — it's idempotent.
#  NOTE: the CellBender sidecar venv .venv-cellbender/ is no longer required.
#  See STATUS_06-05.md for details.)
nvidia-smi
which uv
which Rscript
test -f refs/msigdb_mouse.tsv          && echo "msigdb_mouse.tsv present"          || echo "MISSING — re-run ./setup-remote.sh"
test -f refs/abc_brain_ref.h5ad        && echo "abc_brain_ref.h5ad present"        || echo "MISSING — re-run ./setup-remote.sh"
test -f refs/celltypist_brain_adult.pkl && echo "celltypist_brain_adult.pkl present" || echo "MISSING — re-run ./setup-remote.sh"


# =============================================================================
# Focal cell types — used by every loop below. Edit for placenta run.
# =============================================================================
CELL_TYPES=(
  "Excitatory neurons"
  "Inhibitory neurons"
  "Microglia"
  "Oligodendrocytes"
  "Astrocytes"
  "OPC"
  # placenta: "Trophoblast" "Decidual" "EVT" "Endothelial" "Immune"
)

slugify() { echo "$1" | tr '[:upper:] ' '[:lower:]_' | tr -d ',/'; }


# =============================================================================
# Pick the config — set this ONCE, then re-run the script for the other tissue.
# =============================================================================
CONFIG=config/brain.yaml          # change to config/placenta.yaml for placenta
grep -c "^- id:" "$CONFIG"        # should be 34 (brain) or 23 (placenta)


# =============================================================================
# Phase 0 — Validation (mandatory; ~5 min)
# =============================================================================
tmux new -s phase0_brain -d \
  "uv run python scripts/01_validate.py --config config/brain.yaml \
   2>&1 | tee logs/phase0_brain.log"
tmux new -s phase0_placenta -d \
  "uv run python scripts/01_validate.py --config config/placenta.yaml \
   2>&1 | tee logs/phase0_placenta.log"
# Inspect: results/<tissue>/validation/validation_report.txt


# =============================================================================
# Phase 1 — Ambient RNA (CellBender) — SKIPPED 2026-06-05
# =============================================================================
# CellBender failed across every torch/pyro/numpy combo we tried due to a
# known weakref pickling bug in checkpoint save (GitHub #371, #386, #395).
# snRNA-seq has low ambient anyway; Phase 2 QC catches contaminated cells.
# To revisit: use the official Docker image us.gcr.io/broad-dsde-methods/cellbender:0.3.0
# or switch to SoupX in R. See STATUS_06-05.md §3.


# =============================================================================
# Phase 2 — QC (per-sample MAD + hard floors/caps)
# =============================================================================
# Smoke test first:
#   uv run python -c "import yaml; cfg=yaml.safe_load(open('config/brain.yaml')); \
#     cfg['samples'] = [s for s in cfg['samples'] if s['id']=='E6']; \
#     cfg['results_dir']='results/brain_smoketest'; \
#     yaml.safe_dump(cfg, open('config/brain_smoketest.yaml','w'))"
#   uv run python scripts/02_qc.py --config config/brain_smoketest.yaml
# If clean, run for real:
tmux new -s qc_brain -d \
  "uv run python scripts/02_qc.py --config config/brain.yaml \
   2>&1 | tee logs/qc_brain.log"
tmux new -s qc_placenta -d \
  "uv run python scripts/02_qc.py --config config/placenta.yaml \
   2>&1 | tee logs/qc_placenta.log"
# Inspect: results/<tissue>/tables/02_qc/summary_qc.csv


# =============================================================================
# Phase 3 — Doublet detection (scDblFinder per pool, R subprocess)
# =============================================================================
tmux new -s doublets_brain -d \
  "uv run python scripts/03_doublets.py --config config/brain.yaml \
   2>&1 | tee logs/doublets_brain.log"
tmux new -s doublets_placenta -d \
  "uv run python scripts/03_doublets.py --config config/placenta.yaml \
   2>&1 | tee logs/doublets_placenta.log"


# =============================================================================
# Phase 4 — Concat + HVG + cell cycle
# =============================================================================
# NOTE for placenta: ensure hemo genes (Hbb, Hba) are excluded from HVGs.
# Placenta has biological hemoglobin signal that would dominate PCs.
uv run python scripts/04_integration_prep.py --config "$CONFIG"


# =============================================================================
# Phase 5 — scVI integration (GPU, BF16, ~2-3h)
# =============================================================================
# SMOKE TEST FIRST with --max-epochs 5 on a small subset before full run.
tmux new -s scvi_brain -d \
  "CUDA_VISIBLE_DEVICES=0 uv run python scripts/05_integration.py \
   --config config/brain.yaml 2>&1 | tee logs/scvi_brain.log"
# Wait for brain scVI to finish before starting placenta (single GPU):
tmux new -s scvi_placenta -d \
  "CUDA_VISIBLE_DEVICES=0 uv run python scripts/05_integration.py \
   --config config/placenta.yaml 2>&1 | tee logs/scvi_placenta.log"


# =============================================================================
# Phase 6 — Clustering (Leiden, igraph backend)
# =============================================================================
uv run python scripts/06_clustering.py --config "$CONFIG"


# =============================================================================
# Phase 7 — Annotation (CellTypist + markers, per-cluster majority)
# =============================================================================
uv run python scripts/07_annotation.py --config "$CONFIG"


# =============================================================================
# Phase 7b — Subclustering (loop over ALL focal cell types)
# =============================================================================
for ct in "${CELL_TYPES[@]}"; do
  uv run python scripts/07b_subcluster.py --config "$CONFIG" --celltype "$ct"
done


# =============================================================================
# Phase 7d — Subcluster annotation
# =============================================================================
for ct in "${CELL_TYPES[@]}"; do
  uv run python scripts/07d_subcluster_annotate.py --config "$CONFIG" \
      --celltype "$ct" --markers config/subcluster_markers.yaml
done


# =============================================================================
# Phase 7c — scANVI reference label transfer (WORKSTATION ONLY)
# =============================================================================
# SMOKE TEST FIRST. Brain uses refs/abc_brain_ref.h5ad. Placenta has no ref
# yet — script will exit cleanly.
uv run python scripts/07c_label_transfer.py --config "$CONFIG"


# =============================================================================
# Phase 8a — Composition (propeller via R subprocess)
# =============================================================================
uv run python scripts/08a_composition.py --config "$CONFIG"


# =============================================================================
# Phase 8b — Pseudobulk DE (PyDESeq2). Main + subcluster.
# =============================================================================
# NOTE: For brain, run a sensitivity analysis excluding C5 (sex swap flag).
uv run python scripts/08b_de.py --config "$CONFIG"
for ct in "${CELL_TYPES[@]}"; do
  slug=$(slugify "$ct")
  uv run python scripts/08b_de.py --config "$CONFIG" --subcluster "$slug"
done


# =============================================================================
# Phase 7e — Cell-type counts diagnostic (brain + placenta; ~1 min)
# =============================================================================
# Per-donor × cell-type count CSV used for sanity-checking 8a propeller
# inputs and for paper Table S?. Brain: 3 granularities × (whole + regions).
# Placenta: whole only.
uv run python scripts/07e_celltype_counts.py --config "$CONFIG"
# Output: results/<tissue>/tables/07_annotation/07e_celltype_counts.csv


# =============================================================================
# Phase 8b FOLLOW-UPS — brain only (placenta lacks within_group_across_age)
# =============================================================================
# Three scripts that operate on the master 08b_de_results.csv (NO re-running
# of DE needed). Skip these entirely for placenta.

# --- 1. Developmental disruption analysis ---
# Classifies genes in the within_group_across_age contrast into 5 direction
# classes: universal / relaxed_only (=LOST) / stress_shared (=GAINED) /
# early_only / late_only. Writes summary + long-form gene tables.
uv run python scripts/08b_developmental_disruption.py --config config/brain.yaml
for ct in immune opc_oligodendrocytes astrocytes_ependymal; do
  uv run python scripts/08b_developmental_disruption.py --config config/brain.yaml --subcluster "$ct"
done

# --- 2. Disruption mirror plot + stress-consistency stacked bars ---
# Per (sex × level): mirror bar of LOST (red, left) vs GAINED (blue, right)
# with paired |LFC| boxplots showing effect-size collapse in stress; PLUS
# stress-consistency stacked bars (Early-only / Both-sig / Late-only) per age.
uv run python scripts/08b_followup_plots.py --config config/brain.yaml
for ct in immune opc_oligodendrocytes astrocytes_ependymal; do
  uv run python scripts/08b_followup_plots.py --config config/brain.yaml --subcluster "$ct"
done

# --- 3. k-preserving null shuffle test + within-stratum binomial breakdown ---
# Per (sex × level × celltype): tests whether observed LOST/GAINED counts
# exceed (or fall below) the k-preserving null. Reports binomial enrichment +
# depletion p-values for all 6 disjoint sig-pattern categories (R-only /
# E-only / L-only / R∩E / R∩L / E∩L). Headline figure = 2-panel: delta mirror
# bar + within-stratum 6-bar breakdown per cell type.
# CPU-bound, use processes (built-in via parallel_map use_threads=False).
uv run python scripts/08b_disruption_shuffle_test.py --config config/brain.yaml --n-perm 1000 --n-jobs 16
for ct in immune opc_oligodendrocytes astrocytes_ependymal; do
  uv run python scripts/08b_disruption_shuffle_test.py --config config/brain.yaml --subcluster "$ct" --n-perm 1000 --n-jobs 16
done
# Outputs:
#   results/brain/tables/08b_de/08b_developmental_disruption_summary.csv
#   results/brain/tables/08b_de/08b_developmental_disruption_genes.csv
#   results/brain/tables/08b_de/08b_disruption_shuffle_test.csv (+ subcluster variants)
#   results/brain/plots/08b_de/summary/{disruption,consistency,shuffle_test}/{sex}/{level}.png


# =============================================================================
# Phase 8c — Pathway/GSEA + TF activity. MAIN + SUBCLUSTER, --tf REQUIRED.
# =============================================================================
# SMOKE TEST with one cell type before running full loop.
uv run python scripts/08c_pathways.py --config "$CONFIG" --tf
for ct in "${CELL_TYPES[@]}"; do
  slug=$(slugify "$ct")
  uv run python scripts/08c_pathways.py --config "$CONFIG" --tf --subcluster "$slug"
done


# =============================================================================
# Phase 8d — Trajectory (PAGA + diffusion pseudotime). Supplementary/mechanistic.
# =============================================================================
# Two modes per the trajectory: block in each YAML:
#   - whole-tissue PAGA (no --subcluster): celltypist_broad / celltype_majority
#     connectivity graph + edge diagnostics. Does NOT modify X_umap.
#   - focal-lineage DPT (--subcluster <lineage>): per-donor median pseudotime +
#     mature fraction, MW-U pairwise + KW omnibus (animal as unit; sex strata
#     combined/M/F; all low_n at n≈4). DPT is WITHIN-AGE both tissues (a cross-age
#     brain axis tracks the P1-vs-adult gap, not OPC→MOL maturation — smoke test
#     2026-06-22). No RNA velocity/CellRank (10x Flex exon-only).
# Lineage names = 08c subcluster h5ad basenames (brain) / trophoblast (placenta),
# NOT the CELL_TYPES display names. immune runs two roots (PAM_ATM + Homeostatic;
# state axis, not a lineage). Lineage objects are small (26-66K cells); whole-tissue
# PAGA loads the full annotated object (~661K brain / ~397K placenta) — the long pole.
BRAIN_LINEAGES=(opc_oligodendrocytes astrocytes_ependymal immune)
PLACENTA_LINEAGES=(trophoblast)

# --- focal-lineage DPT (picks the array matching your $CONFIG) ---
if [[ "$CONFIG" == *placenta* ]]; then
  LINEAGES=("${PLACENTA_LINEAGES[@]}")
else
  LINEAGES=("${BRAIN_LINEAGES[@]}")
fi
for lin in "${LINEAGES[@]}"; do
  uv run python scripts/08d_trajectory.py --config "$CONFIG" --subcluster "$lin"
done

# --- whole-tissue PAGA (tmux — neighbors on the full object is ~15-30 min) ---
tmux new -s traj_whole -d \
  "uv run python scripts/08d_trajectory.py --config $CONFIG \
   2>&1 | tee logs/08d_wholePAGA.log"
# Outputs:
#   results/<tissue>/tables/08d_trajectory{,_subcluster_<lin>}/*_dpt_*.csv
#   results/<tissue>/tables/08d_trajectory{,_subcluster_<lin>}/*_paga_edge_diagnostics.csv
#   results/<tissue>/plots/08d_trajectory_subcluster_<lin>/{paga,pseudotime/<age>/<root>,per_donor/<sex>}/


# =============================================================================
# Phase 8e — Cell-cell communication. MAIN + SUBCLUSTER.
# =============================================================================
uv run python scripts/08e_communication.py --config "$CONFIG" --zscore-rows
for ct in "${CELL_TYPES[@]}"; do
  slug=$(slugify "$ct")
  uv run python scripts/08e_communication.py --config "$CONFIG" \
      --subcluster "$slug" --zscore-rows
done


# =============================================================================
# CHECKPOINT — Run the above through 8e for BOTH tissues before 8f/8g.
# =============================================================================


# =============================================================================
# Phase 8f — Cross-tissue (placenta → brain cascades). Six views.
# =============================================================================
uv run python scripts/08f_cross_tissue.py \
    --brain-config config/brain.yaml \
    --placenta-config config/placenta.yaml


# =============================================================================
# Phase 8g — Cross-age persistence (brain only). Six views.
# =============================================================================
uv run python scripts/08g_cross_age.py --config config/brain.yaml


# =============================================================================
# Phase 9 — Cross-species RRHO2 validation
# =============================================================================
# Requires Stage-1 human downloads in data/human_validation/.
# SMOKE TEST FIRST on Velmeshev only (the most-complete dataset).
# Per-dataset loaders in 09_cross_species_validation.py were stubs at last
# check — verify before running full thing.
uv run python scripts/09_cross_species_validation.py \
    --config config/brain.yaml \
    --celltype-map config/cross_species_celltype_map.yaml


# =============================================================================
# Headline files for paper assembly:
# =============================================================================
#   results/brain/tables/08f_cross_tissue/08f_lr_cross_tissue.csv
#   results/brain/tables/08g_cross_age/08g_core_signature_genes.csv
#   results/brain/tables/08b_de/08b_de_results.csv
#   results/brain/tables/08b_de/08b_developmental_disruption_summary.csv
#   results/brain/tables/08b_de/08b_disruption_shuffle_test.csv
#   results/brain/tables/08c_pathways/08c_pathway_results.csv
#   results/brain/tables/08c_pathways/08c_tf_activity.csv
#   results/{tissue}/tables/09_cross_species/*_rrho2.csv  (Phase 9)
#
# Paper-quality figures:
#   results/brain/plots/08b_de/summary/disruption/combined/whole.png
#       — mirror bar of LOST vs GAINED + effect-size collapse boxes
#   results/brain/plots/08b_de/summary/shuffle_test/combined/whole.png
#       — k-preserving null test + within-stratum 6-bar breakdown
