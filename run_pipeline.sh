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
# Phase 7b: subclustering (one cell type at a time)
# ============================================================================
# Re-runs HVG+scVI on the subset to resolve subtypes. Output 'subcluster' is
# INTEGER ids; 7d names them. Needs >=50 cells of the type.
uv run python scripts/07b_subcluster.py --config config/dev_split.yaml \
    --celltype "Excitatory neurons"
# repeat per lineage: "Inhibitory neurons", "Microglia", "Oligodendrocytes"...
open results/dev/plots/07b_subcluster/excitatory_neurons/umap_subclusters.png
cat results/dev/tables/07b_subcluster/07b_subcluster_excitatory_neurons_markers.csv
git add scripts/07b_subcluster.py && git commit -m "phase 7b: subclustering"


# ============================================================================
# Phase 7d: subcluster annotation (name the 7b integer ids)
# ============================================================================
# Track A: CellTypist majority_voting per cluster (if a model is configured).
# Track B: literature marker scoring from config/subcluster_markers.yaml,
# aggregated per cluster (mean score -> argmax). Already cluster-level by
# construction. Writes obs['subcluster_name'] back to the 7b h5ad.
uv run python scripts/07d_subcluster_annotate.py --config config/dev_split.yaml \
    --celltype "Excitatory neurons" --markers config/subcluster_markers.yaml
# placenta marker-only:  --celltype trophoblast --no-celltypist
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
uv run python scripts/08b_de.py --config config/dev_split.yaml
# subcluster DE:  --subcluster excitatory_neurons
# --no-expr-matrix to skip the offline-audit CSV; --expr-sig-only to trim it.
open results/dev/plots/08b_de/early_vs_relaxed_per_age/age-4W/Excitatory_neurons/volcano.png
cat results/dev/tables/08b_de/08b_de_results.csv
cat results/dev/tables/08b_de/08b_de_gene_expression_per_sample.csv
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
# TF activity (--tf or YAML run_tf_activity: true): ULM on same DE stats vs
# CollecTRI mouse network. Writes 08c_tf_activity.csv + barplot/volcano per
# celltype + TF*celltype heatmap per contrast. Needs network (omnipath).
Rscript scripts/fetch_genesets.R --out refs/msigdb_mouse.tsv      # one-time
uv run python scripts/08c_pathways.py --config config/dev_split.yaml --tf
# subcluster pathways:  --subcluster excitatory_neurons
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
# NEXT: 8e cell-cell communication (LIANA+ consensus),
#       8f cross-tissue (placenta -> brain cascades),
#       8g cross-age / persistence (derived from 8b/8c tables; lightest).
# ============================================================================


# ============================================================================
# WORKSTATION RUN — same scripts, real config
# ============================================================================
# Brain (34 samples): swap every --config above to config/brain.yaml.
# Placenta (23 samples): config/placenta.yaml.
# Skip the dev_split_h5.py step entirely; --config <real>.yaml runs everything
# on real samples with real donor_ids. Phase 1 (CellBender) is required (GPU).
# Phase 7c (scANVI) requires the reference: block + staged Allen BCA h5ad.
# Pre-workstation TODO: refine pathways stress sets, train ABC CellTypist .pkl.
