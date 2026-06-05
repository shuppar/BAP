# run_pipeline.sh — manual, not history. Tight start-to-finish.
# Not executable. Edit configs, run commands by hand in the order below.
# Git history is the changelog; this file is "what to run today, in what order".
#
# Two workflows below: DEV (laptop smoke test) and WORKSTATION (real run).
# All phase scripts are identical in both — only the config differs.


# ============================================================================
# DEV setup (laptop only) — one-time pre-step before Phase 0
# ============================================================================
# dev has 1 sample per group; pseudobulk needs n>=3. The splitter reads the
# 3 dev .h5 files, writes 9 split .h5 files (random barcode partition), and
# emits config/dev_split.yaml listing the 9 pseudo-samples (donor_id suffixed
# _ps1/2/3, group/age/sex/pool inherited from parent). NOT part of workstation.
# Pseudo-donors are random partitions of one animal -> numbers MEANINGLESS,
# smoke test of code paths only.
uv run python dev_split_h5.py --config config/dev.yaml --n 3
ls data/dev_split/                          # 9 .h5 files
cat config/dev_split.yaml | head -40
# From here, use --config config/dev_split.yaml for every dev phase.


# ============================================================================
# Phase 0: validate (manifest, balance, sex-check, fingerprints) — no compute
# ============================================================================
# Hard-fails on sample-swaps, missing metadata, declared-vs-inferred sex mismatch.
# Run THIS before any compute-heavy phase, per project doc §5.
uv run python scripts/01_validate.py --config config/dev_split.yaml
open results/dev/validation/manifest_balance_matrix.png
open results/dev/validation/sex_check_scatter.png
cat results/dev/validation/sex_check.csv


# ============================================================================
# Phase 1: ambient RNA (CellBender) — WORKSTATION ONLY, GPU
# ============================================================================
# Skipped on laptop dev. On workstation:
#   .venv-cellbender (separate venv — CellBender pins incompatible torch).
#   GPU-parallel 2 samples, ~1-2h each, 150 epochs.
#   uv run --python .venv-cellbender python scripts/02_ambient.py --config config/brain.yaml
# (script still TBD)


# ============================================================================
# Phase 2: per-sample QC
# ============================================================================
# Per-sample MAD bounds + hard floors (min_counts=500, min_genes=200) + hard
# caps (mt<=1%, hemo<=5%). Cohort-outlier flag (sample median UMI/genes >3
# cohort-MADs below median) catches failed-prep samples whose own MAD made
# their per-sample bounds permissive.
uv run python scripts/02_qc.py --config config/dev_split.yaml
open results/dev/plots/02_qc/E1-4WkM1_ps1_violin_prepost.png
cat results/dev/tables/02_qc/02_qc_summary.csv
git add scripts/02_qc.py && git commit -m "phase 2: per-sample QC"


# ============================================================================
# Phase 3: doublets (scDblFinder per pool, via R subprocess)
# ============================================================================
# Doublets form within a capture -> run per pool, not per sample. R subprocess
# (not rpy2) for crash isolation. scDblFinder gets samples= arg so simulated
# doublets respect within-sample boundaries. Classified doublets removed.
uv run python scripts/03_doublets.py --config config/dev_split.yaml
open results/dev/plots/03_doublets/Pool1_doublet_score_dist.png
cat results/dev/tables/03_doublets/summary_doublets.csv
git add scripts/03_doublets.py scripts/run_scdblfinder.R && git commit -m "phase 3: doublets"


# ============================================================================
# Phase 4: concat + HVG + cell-cycle scoring
# ============================================================================
# Builds layers['lognorm'] (Phase 5 drops it after pre-int UMAP; downstream
# recomputes via _utils.add_lognorm). HVGs via seurat_v3 on raw counts,
# batch_key=pool. Generic exclusions: mito/ribo/hemo/sex-linked. Tirosh
# cell-cycle scores added to .obs (cc conditioning optional in scVI).
uv run python scripts/04_integration_prep.py --config config/dev_split.yaml
open results/dev/plots/04_integration_prep/hvg_dispersion.png
open results/dev/plots/04_integration_prep/pre_integration_umap_by_pool.png
cat results/dev/tables/04_integration_prep/04_integration_prep_summary.csv
git add scripts/04_integration_prep.py && git commit -m "phase 4: HVG + cell cycle"


# ============================================================================
# Phase 5: scVI integration
# ============================================================================
# batch_key=pool. categorical covariates = none (biology stays in latent).
# Continuous: pct_counts_mt. Accelerator+precision auto-picked (BF16 GPU /
# FP32 CPU). Dev caps max_epochs at 50 (tiny data). Workstation: 400 epochs,
# bs=1024, early stop patience 30.
uv run python scripts/05_integration.py --config config/dev_split.yaml
open results/dev/plots/05_integration/umap_post_integration_by_pool.png
open results/dev/plots/05_integration/scvi_loss_curve.png
cat results/dev/tables/05_integration/05_integration_scvi_history.csv
git add scripts/05_integration.py && git commit -m "phase 5: scVI integration"


# ============================================================================
# Phase 6: clustering (multi-res Leiden, igraph backend)
# ============================================================================
# Sweep resolutions, auto-pick via geometric knee. Override:  --resolution 0.6
# Uses igraph backend (faster, scanpy's future default).
uv run python scripts/06_clustering.py --config config/dev_split.yaml
open results/dev/plots/06_clustering/resolution_selection.png
open results/dev/plots/06_clustering/clustree.png
open results/dev/plots/06_clustering/cluster_composition_by_sample.png
git add scripts/06_clustering.py && git commit -m "phase 6: leiden (igraph)"


# ============================================================================
# Phase 7: annotation (CellTypist + markers, per-cluster majority)
# ============================================================================
# Per-cluster majority voting (CellTypist convention) — not per-cell argmax.
# Low-purity (<60% majority) clusters announced in stdout for manual review.
# Dev (4W brain) has no built-in adult mouse CellTypist model -> marker track
# only. Workstation: train .pkl on Allen Brain Cell Atlas, point at it.
uv run python scripts/07_annotation.py --config config/dev_split.yaml
open results/dev/plots/07_annotation/umap_celltype_annotation.png
open results/dev/plots/07_annotation/marker_dotplot.png
cat results/dev/tables/07_annotation/07_annotation_summary.csv
git add scripts/07_annotation.py && git commit -m "phase 7: annotation"


# ============================================================================
# Phase 7b: subclustering — LOOP THROUGH ALL FOCAL CELL TYPES
# ============================================================================
# Re-runs HVG+scVI on the subset to resolve subtypes. Output 'subcluster' is
# INTEGER ids; 7d names them. Needs >=50 cells of the type.
#
# RUN FOR EVERY FOCAL CELL TYPE — don't stop at one. The full pipeline
# expects subclustered objects for all of these. The exact list is dataset-
# specific; below is the brain default. Edit before running placenta.
CELL_TYPES=(
  "Excitatory neurons"
  "Inhibitory neurons"
  "Microglia"
  "Oligodendrocytes"
  "Astrocytes"
  "OPC"
  # add "Radial glia / NPCs" for P1 brain; "trophoblast" / "decidual" for placenta
)
for ct in "${CELL_TYPES[@]}"; do
  uv run python scripts/07b_subcluster.py --config config/dev_split.yaml \
      --celltype "$ct"
done
# Skipped types print "ERROR: only N cells — too few to subcluster reliably"
# and exit cleanly; that's expected for rare types.
open results/dev/plots/07b_subcluster/excitatory_neurons/umap_subclusters.png
cat results/dev/tables/07b_subcluster/07b_subcluster_excitatory_neurons_markers.csv
git add scripts/07b_subcluster.py && git commit -m "phase 7b: subclustering"


# ============================================================================
# Phase 7d: subcluster annotation — LOOP THROUGH ALL FOCAL CELL TYPES
# ============================================================================
# Track A: CellTypist majority_voting per cluster (if a model is configured).
# Track B: literature marker scoring from config/subcluster_markers.yaml,
# aggregated per cluster (mean score -> argmax). Already cluster-level by
# construction. Writes obs['subcluster_name'] back to the 7b h5ad.
for ct in "${CELL_TYPES[@]}"; do
  uv run python scripts/07d_subcluster_annotate.py --config config/dev_split.yaml \
      --celltype "$ct" --markers config/subcluster_markers.yaml
done
# placenta marker-only (no CellTypist model):  --no-celltypist
open results/dev/plots/07b_subcluster/excitatory_neurons/subcluster_names_umap.png
cat results/dev/tables/07d_subcluster_annotate/07d_subcluster_excitatory_neurons_annotation.csv
git add scripts/07d_subcluster_annotate.py config/subcluster_markers.yaml \
    && git commit -m "phase 7d: subcluster names"


# ============================================================================
# Phase 7c: scANVI reference label transfer — WORKSTATION ONLY
# ============================================================================
# Needs a reference: block in YAML + a labeled reference h5ad (multi-GB).
# Exits cleanly on dev with "no 'reference:' block" — by design.
# On workstation:
#   uv pip install "git+https://github.com/alleninstitute/abc_atlas_access.git"
#   # build labeled ref h5ad (cell-type + region on .obs), then:
#   uv run python scripts/prepare_reference.py --config config/brain.yaml \
#       --source /path/to/raw_reference.h5ad --max-cells-per-label 2000
#   uv run python scripts/07c_label_transfer.py --config config/brain.yaml
# NOTE: the same Allen BCA build also trains CellTypist .pkl for adult 4W/3mo
# brain (no built-in adult model exists). Point annotation.celltypist_models
# at the .pkl and re-run Phase 7.
git add scripts/07c_label_transfer.py scripts/prepare_reference.py \
        scripts/_utils.py scripts/build_yaml.py config/ \
    && git commit -m "phase 7c: scANVI ref transfer (workstation)"


# ============================================================================
# Phase 8a: composition (propeller / speckle via R subprocess)
# ============================================================================
# Animal (donor_id) is the unit. propeller = arcsin/logit + limma moderation
# (good for small n). 2-group=t, 3-group=ANOVA F. sex/pool enter the design.
# Driven by the declarative `contrasts:` block in config.
# R deps once: R -e 'BiocManager::install(c("speckle","limma"))'
uv run python scripts/08a_composition.py --config config/dev_split.yaml
open results/dev/plots/08a_composition/early_vs_relaxed_per_age/age-4W/stacked_bar.png
open results/dev/plots/08a_composition/early_vs_relaxed_per_age/age-4W/propeller_effects.png
cat results/dev/tables/08a_composition/08a_composition_results.csv
git add scripts/08a_composition.py scripts/run_propeller.R \
    && git commit -m "phase 8a: composition (propeller)"


# ============================================================================
# Phase 8b: pseudobulk DE (PyDESeq2)
# ============================================================================
# Pseudobulk ONLY — sum raw counts per donor per cell type. Design from
# contrast spec (~ sex + pool + group). Volcanoes label top genes via
# var['symbol']. flag/reliability/note carried per row. Also writes
# per-sample expression matrix of DE genes for OFFLINE audit
# (08b_de_gene_expression_per_sample.csv) — join to 8c leading-edge on
# (celltype, gene) for per-sample levels.
#
# Main run (coarse cell types):
uv run python scripts/08b_de.py --config config/dev_split.yaml
# REQUIRED — also run subcluster DE for every focal cell type that finished 7b.
# Each writes results suffixed: 08b_de_results_subcluster_{slug}.csv
# (don't overwrite the main 08b_de_results.csv). 8e/8f/8g read whichever is
# applicable; subcluster DE is what powers subcluster-level pathway/CCC.
for ct in "${CELL_TYPES[@]}"; do
  slug=$(echo "$ct" | tr '[:upper:] ' '[:lower:]_' | tr -d ',/')
  uv run python scripts/08b_de.py --config config/dev_split.yaml --subcluster "$slug"
done
# --no-expr-matrix to skip offline-audit CSV; --expr-sig-only to trim it.
open results/dev/plots/08b_de/early_vs_relaxed_per_age/age-4W/Excitatory_neurons/volcano.png
cat results/dev/tables/08b_de/08b_de_results.csv
cat results/dev/tables/08b_de/08b_de_results_subcluster_excitatory_neurons.csv
git add scripts/08b_de.py && git commit -m "phase 8b: DE + per-sample expr matrix"


# ============================================================================
# Phase 8c: pathway/GSEA (decoupler) + TF activity (CollecTRI)
# ============================================================================
# GSEA on DE Wald stats. Mouse MSigDB via msigdbr (MH+M2+M5+M8), exported once
# to refs/msigdb_mouse.tsv (decoupler's MSigDB-mouse fetch is broken).
# FDR within each collection (kosher); FDR_pooled kept as ref. Side-by-side
# per-collection panels (dotplot, volcano, celltype heatmap) + running-enrich
# per top hit. M8 (cell-type sets) esp. useful for subcluster mode.
# OFFLINE audit: pathway_leading_edge.csv — genes driving each significant
# pathway with log2FC + direction (join to 8b expr matrix for per-sample).
#
# --tf IS REQUIRED for the full downstream chain. 8f (cross-tissue) and 8g
# (cross-age) both read 08c_tf_activity.csv when available; without --tf
# the TF concordance view in 8f is silently skipped, and the TF arm of 8g
# can't run. ALWAYS pass --tf in production unless deliberately deferring.
# (YAML alternative: pathways.run_tf_activity: true)
Rscript scripts/fetch_genesets.R --out refs/msigdb_mouse.tsv      # one-time
uv run python scripts/08c_pathways.py --config config/dev_split.yaml --tf
# Subcluster pathways + TF — run for every focal cell type:
for ct in "${CELL_TYPES[@]}"; do
  slug=$(echo "$ct" | tr '[:upper:] ' '[:lower:]_' | tr -d ',/')
  uv run python scripts/08c_pathways.py --config config/dev_split.yaml --tf \
      --subcluster "$slug"
done
open results/dev/plots/08c_pathways/early_vs_relaxed_per_age/age-4W/celltype_pathway_heatmap_panels.png
open results/dev/plots/08c_pathways/early_vs_relaxed_per_age/age-4W/Excitatory_neurons/gsea_volcano_panels.png
open results/dev/plots/08c_pathways/early_vs_relaxed_per_age/age-4W/Excitatory_neurons/tf_activity_barplot.png
open results/dev/plots/08c_pathways/early_vs_relaxed_per_age/age-4W/tf_activity_heatmap.png
cat results/dev/tables/08c_pathways/08c_pathway_results.csv
cat results/dev/tables/08c_pathways/08c_pathway_leading_edge.csv
cat results/dev/tables/08c_pathways/08c_tf_activity.csv
git add scripts/08c_pathways.py scripts/fetch_genesets.R \
    && git commit -m "phase 8c: GSEA + leading-edge + TF activity"


# ============================================================================
# Phase 8d: trajectory (PAGA + diffusion pseudotime)
# ============================================================================
# NO RNA velocity: 10x Flex is probe-based (exon-only), can't resolve
# spliced/unspliced. No CellRank either (without velocity it duplicates PAGA).
# Works on CELLS (not donors) — runs on the same annotated object.
# All ages treated identically; DPT group comparison runs pooled + per-age,
# age-split rows carry pool_age_confound caveat. OFFLINE audit:
# trajectory_paga_edge_diagnostics.csv — per cell-type-pair edge: connectivity,
# shared top genes, ambient_driven flag, doublet rate, %mt. Use to audit
# surprising PAGA edges offline (no workstation needed).
uv run python scripts/08d_trajectory.py --config config/dev_split.yaml
# --root-celltype "Radial glia / NPCs"  to set DPT root explicitly
open results/dev/plots/08d_trajectory/paga/paga_by_celltype.png
open results/dev/plots/08d_trajectory/pseudotime/dpt_umap.png
cat results/dev/tables/08d_trajectory/08d_trajectory_paga_edge_diagnostics.csv
cat results/dev/tables/08d_trajectory/08d_trajectory_dpt_group_comparison.csv
git add scripts/08d_trajectory.py && git commit -m "phase 8d: PAGA + DPT"


# ============================================================================
# Phase 8e: cell-cell communication (LIANA+ consensus)
# ============================================================================
# Three analytical arms in one script:
#   1. Baseline   — rank_aggregate per group×age (pooled cells)
#   2. Differential — df_to_lr on 8b Wald stats (ES-v-Rel, LS-v-Rel, ES-v-LS)
#   3. Per-donor  — rank_aggregate per donor → group-level statistics
#
# Output tree (5 subfolders):
#   01_overview/          chord comparisons, pathway heatmap, LR trajectory
#   02_baseline_per_group/{group}_{age}/  chord, network, dotplots
#   03_differential/{contrast}_{age}/     volcano, dotplot
#   03_differential/delta_heatmaps/       CLUSTERED Δ heatmaps (focused + full)
#   03_differential/rank_rank/            signature concordance scatters
#   04_sender_receiver/   bubble + Δ heatmaps for all 3 group pairs
#   05_per_donor/{age}/   stripplots + Δ bar with FDR annotations
#
# FLAGS THAT MUST FIRE for full output:
#   default: baseline + differential + per-donor all run
#   --n-perms 1000        permutation specificity for CellPhoneDB (workstation;
#                         use 10 on dev — too few cells for stable permutations)
#   --zscore-rows         adds Z-scored variants of clustered Δ heatmaps
#                         (pattern view alongside absolute Δ)
# Optional skip:
#   --skip-per-donor      skip arm 3 (faster smoke test; loses statistics)
#
# Run main (all coarse cell types):
uv run python scripts/08e_communication.py --config config/dev_split.yaml \
    --n-perms 10 --zscore-rows
# Subcluster runs — LOOP through every focal cell type (separate output tree
# per subcluster so they don't collide; tables/plots auto-suffixed with slug):
for ct in "${CELL_TYPES[@]}"; do
  slug=$(echo "$ct" | tr '[:upper:] ' '[:lower:]_' | tr -d ',/')
  uv run python scripts/08e_communication.py --config config/dev_split.yaml \
      --subcluster "$slug" --n-perms 10 --zscore-rows
done
open results/dev/plots/08e_communication/01_overview/chord_comparison_4W.png
open results/dev/plots/08e_communication/03_differential/delta_heatmaps/delta_lr_heatmap_Early_Stress_vs_Relaxed_4W_focused.png
open results/dev/plots/08e_communication/03_differential/rank_rank/rank_rank_Early_Stress-Relaxed_vs_Late_Stress-Relaxed_4W.png
cat results/dev/tables/08e_communication/08e_lr_baseline.csv
cat results/dev/tables/08e_communication/08e_lr_quantified.csv
git add scripts/08e_communication.py scripts/_08e_plots_*.py \
    && git commit -m "phase 8e: CCC (baseline + differential + per-donor)"


# ============================================================================
# Phase 8f: cross-tissue (placenta → brain cascades)
# ============================================================================
# Two biologically aligned arms: E12.5 placenta (Early) → P1/4W/3mo brain Early;
# E18.5 placenta (Late) → P1/4W/3mo brain Late. P1 Late carries pool-confound
# flag (propagated to every output row).
#
# Six analytical views in one script:
#   1. DEG overlap         hypergeom per ct_pair × direction
#   2. RRHO                rank-rank hypergeometric (custom NumPy)
#   3. Pathway concordance NES-sign-based, from 8c pathway tables
#   4. LR cross-tissue     placental ligand × brain receptor (KEY FILE)
#                          — stress_axis column flags GR/MR/CRH/cytokine genes
#   5. TF concordance      mirror of view 3 on 8c TF activity table
#                          — REQUIRES 8c run with --tf
#   6. Overlap enrichment  ORA of cross-tissue overlap genes vs MSigDB
#                          — REQUIRES refs/msigdb_mouse.tsv (built by 8c step)
#
# UPSTREAM REQUIREMENTS — confirm BEFORE running:
#   ✓ Both tissues completed Phase 8b (08b_de_results.csv)
#   ✓ Both tissues completed Phase 8c WITH --tf (else view 5 silently skips)
#   ✓ refs/msigdb_mouse.tsv exists (else view 6 silently skips)
#
# Dev smoke-test only — duplicates brain dir, runs cross-tissue against itself.
# DO NOT use --dev-test on real data; on workstation drop the flag.
cp -r results/dev results/dev_placenta
sed 's|results_dir: results/dev|results_dir: results/dev_placenta|' \
    config/dev_split.yaml > config/dev_placenta.yaml
uv run python scripts/08f_cross_tissue.py \
    --brain-config config/dev_split.yaml \
    --placenta-config config/dev_placenta.yaml \
    --dev-test
open results/dev/plots/08f_cross_tissue/02_deg_overlap/Early_devtest_4W/effect_size_scatter.png
open results/dev/plots/08f_cross_tissue/05_lr_cross_tissue/Early_devtest_4W/lr_cross_tissue_scatter.png
cat results/dev/tables/08f_cross_tissue/08f_lr_cross_tissue.csv         # KEY
cat results/dev/tables/08f_cross_tissue/08f_tf_concordance.csv
cat results/dev/tables/08f_cross_tissue/08f_overlap_enrichment.csv
git add scripts/08f_cross_tissue.py && git commit -m "phase 8f: cross-tissue cascades"


# ============================================================================
# Phase 8g: cross-age persistence (derived from 8b/8c tables; lightest)
# ============================================================================
# Operates entirely on existing 8b/8c CSVs — no re-running of DE/GSEA.
# Brain-only by design (placenta has incomplete cross-age factorial; script
# exits cleanly with a warning if tissue: placenta).
#
# Six analytical views in one script:
#   1. Gene-level persistence       — classify each (celltype, gene) per arm
#   2. Pathway-level persistence    — same classification on 8c GSEA
#   3. TF-level persistence         — same on 8c TF activity (needs 8c --tf)
#   4. Effect-size trajectories     — top persistent features, log2FC/NES vs age
#   5. Early vs Late at each age    — hypergeometric overlap + Spearman ρ
#   6. Cross-arm core signature     — features persistent in BOTH arms
#                                     (paper-quality table)
#
# Persistence classes: persistent / resolving_early / established_late /
# P1_only / transient_4W / emergent_3mo / P1_3mo_only / persistent_directionswap.
# Same-direction sign required for "persistent".
#
# DEV LIMITATION — DO NOT EXPECT MEANINGFUL OUTPUT ON DEV:
#   dev_split.yaml subsets to 4W only (one M sample per group), by design,
#   to keep smoke tests fast. With only one age, EVERY classification will
#   be 'transient_4W' or 'none' — there's nothing to persist across. This
#   exercises the code paths but not the biology. To see persistence
#   working, the script needs ≥2 ages in the input data, ideally all 3.
#   On workstation (config/brain.yaml: P1 + 4W + 3mo), the full class
#   spectrum opens up.
#
# UPSTREAM REQUIREMENTS:
#   ✓ 08b_de_results.csv exists                    (gates views 1, 5, 6)
#   ✓ 08c_pathway_results.csv exists               (gates views 2, 6)
#   ✓ 08c_tf_activity.csv exists (08c run --tf)    (gates view 3)
#
# Run on dev (smoke test; expect 'transient_4W' for everything):
uv run python scripts/08g_cross_age.py --config config/dev_split.yaml
open results/dev/plots/08g_cross_age/01_gene_persistence/genes_persistence_class_barplot.png
cat results/dev/tables/08g_cross_age/08g_gene_persistence.csv
git add scripts/08g_cross_age.py && git commit -m "phase 8g: cross-age persistence"


# ============================================================================
# WORKSTATION RUN — same scripts, real config
# ============================================================================
# Brain (34 samples): swap every --config above to config/brain.yaml.
# Placenta (23 samples): config/placenta.yaml.
# Skip the dev_split_h5.py step entirely; --config <real>.yaml runs everything
# on real samples with real donor_ids. Phase 1 (CellBender) is required (GPU).
# Phase 7c (scANVI) requires the reference: block + staged Allen BCA h5ad.
# Pre-workstation TODO: refine pathways stress sets, train ABC CellTypist .pkl.
#
# IMPORTANT — DO NOT SKIP FLAGS:
#   8b: run with --subcluster $slug for EVERY focal cell type after main run
#   8c: ALWAYS pass --tf  (gates 8f view 5 TF concordance and 8g view 3;
#                          can't recover without re-running 8c)
#   8e: drop --n-perms 10 → use --n-perms 1000 (default) on workstation
#   8e: run --subcluster $slug for every focal cell type after main run
#   8f: drop --dev-test; pass real brain-config and placenta-config
#   8f: do NOT use the cp/sed duplicate trick on real data
#   8g: BRAIN ONLY (placenta exits cleanly with warning); needs all 3 ages
#       in the source data — which config/brain.yaml has by default
#
# After both tissues finish 8a–8e independently:
#   - 8f cross-tissue (requires both brain and placenta 8b/8c output)
#   - 8g cross-age (brain only)
