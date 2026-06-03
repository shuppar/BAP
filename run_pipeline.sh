# =============================================================================
# snRNA-seq prenatal-stress pipeline — run manual
# =============================================================================
# Run start to finish, copy-paste, top to bottom. Not executed as a script.
# Not a changelog (see git). Per phase: deps -> command -> key outputs -> commit.
#
# 3 groups: Early_Stress / Late_Stress / Relaxed(ref). Brain (P1/4W/3mo, n=34),
# placenta (E12.5/E18.5, n=23). Pools 1-4 = scVI batch_key. No dam ID -> pup is
# the statistical unit. See snRNAseq_project_summary.md for design + confounds.
#
# Dirs:  data tars ~/Downloads/BAP_data_1/processed_*_Pool{1..4}.tar
#        work       ~/Downloads/BAP_data_1/Analysis
# Configs: dev.yaml (3 samples x 500 cells, laptop test) | brain.yaml | placenta.yaml
# Swap --config to brain/placenta for the real run (workstation, GPU).
# =============================================================================


# --- Step 0: env (one-time) --------------------------------------------------
# uv (not conda — conda channels blocked at firewall). R for scDblFinder subprocess.
curl -LsSf https://astral.sh/uv/install.sh | sh
echo 'source $HOME/.local/bin/env' >> ~/.zshrc && source $HOME/.local/bin/env
brew install r
cd ~/Downloads/BAP_data_1/Analysis
uv sync


# --- Step 1: extract per-sample files from pool tars -------------------------
# Keep only filtered h5, raw h5 (for CellBender), metrics_summary.csv. Skips BAMs
# etc. -> ~80% disk saved.
mkdir -p data/Pool{1,2,3,4}
for P in 1 2 3 4; do
  tar -xvf ~/Downloads/BAP_data_1/processed_260411_Shiv_FLEX_260411_Shiv_FLEX_Pool$P.tar \
    -C data/Pool$P \
    --include='*sample_filtered_feature_bc_matrix.h5' \
    --include='*sample_raw_feature_bc_matrix.tar.gz' \
    --include='*metrics_summary.csv'
done
find data -name 'sample_filtered_feature_bc_matrix.h5' | wc -l   # expect 57


# --- Step 2: build configs from sample_metadata.csv --------------------------
# CSV is source of truth; YAMLs are generated (incl. annotation + reference blocks).
# Re-run whenever the CSV changes.
uv run python scripts/build_yaml.py
grep -c "^- id:" config/brain.yaml      # 34
grep -c "^- id:" config/placenta.yaml   # 23


# --- Phase 0: validate (MANDATORY first; ~1 min) -----------------------------
# Manifest, sex check (Y/Xist), fingerprints, confound tables. Catches sample
# swaps / missing metadata in minutes vs hours into compute. Review before Phase 2.
uv run python scripts/01_validate.py --config config/dev.yaml
ls results/dev/validation/                                # 6 files
git add scripts/01_validate.py config/ && git commit -m "phase 0: validation"


# --- Phase 2: per-sample QC --------------------------------------------------
# Per-sample MAD bounds + hard floors (min_counts 500, min_genes 200) + caps
# (snRNA: mt<=1%, hemo<=5%). Cohort-outlier flag catches failed preps (n>=5).
uv run python scripts/02_qc.py --config config/dev.yaml
open results/dev/plots/02_qc/E1-4WkM1_thresholds.png
cat results/dev/tables/summary_qc.csv                     # check cohort_outlier col
git add scripts/02_qc.py && git commit -m "phase 2: per-sample QC"


# --- Phase 3: doublets (scDblFinder per pool) --------------------------------
# Per pool, not per sample: doublets form within a capture. R via subprocess
# (not rpy2 — process isolation). One-time R deps:
R -e 'install.packages(c("optparse","Matrix","BiocManager"), repos="https://cloud.r-project.org"); BiocManager::install(c("scDblFinder","SingleCellExperiment","BiocParallel"), update=FALSE, ask=FALSE)'

uv run python scripts/03_doublets.py --config config/dev.yaml
cat results/dev/tables/summary_doublets.csv
git add scripts/03_doublets.py scripts/run_scdblfinder.R && git commit -m "phase 3: doublets"


# --- Phase 4: concat + lognorm + cell cycle + HVG ----------------------------
# Raw counts in .X (scVI needs them); lognorm in a layer, dropped after Phase 5
# to save disk. Cell cycle scored here (S/G2M/phase/cc_difference) for later use,
# NOT regressed. HVG = seurat_v3, batch_key=pool, w/ mito/ribo/hemo/sex exclusions.
uv run python scripts/04_integration_prep.py --config config/dev.yaml
open results/dev/plots/04_integration_prep/cell_cycle_scores.png   # 2 clouds at P1, ~1 at 3mo
open results/dev/plots/04_integration_prep/hvg_dispersion.png
git add scripts/04_integration_prep.py && git commit -m "phase 4: concat+lognorm+cell cycle+HVG"


# --- Phase 5: scVI integration -----------------------------------------------
# batch_key=pool ONLY. No categorical covariates — age/group/sex are biology and
# must stay in the latent (putting them in would integrate the stress signal away).
# Cell cycle conditioning OFF by default (proliferation may be real stress biology);
# enable via scvi.condition_cell_cycle in YAML if phase UMAP shows it driving clusters.
# GPU+BF16 auto on Ada, CPU fallback on Mac.
uv run python scripts/05_integration.py --config config/dev.yaml
open results/dev/plots/05_integration/umap_post_integration.png    # pools mix, biology splits
open results/dev/plots/05_integration/umap_post_phase.png          # cell cycle driving clusters?
git add scripts/05_integration.py && git commit -m "phase 5: scVI integration"


# --- Phase 6: clustering (multi-res Leiden + knee selection) -----------------
# Leiden at [0.2..2.0]; resolution auto-picked by geometric knee on the
# n_clusters-vs-resolution curve (one method, not silhouette too). Override if
# the knee plot disagrees with your eye.
uv add leidenalg
uv run python scripts/06_clustering.py --config config/dev.yaml
# uv run python scripts/06_clustering.py --config config/dev.yaml --resolution 0.6
open results/dev/plots/06_clustering/resolution_selection.png      # inspect FIRST
open results/dev/plots/06_clustering/cluster_composition_by_sample.png  # single-sample = batch artifact
git add scripts/06_clustering.py pyproject.toml uv.lock && git commit -m "phase 6: clustering"


# --- Phase 7: annotation (CellTypist + markers + composition) ----------------
# Two tracks: CellTypist reference transfer (per-age model from YAML) + marker
# scoring. Provisional labels from top marker score so composition has real names.
# Composition plots are PRELIMINARY — correct manual_annotation, then re-run.
uv add celltypist
# Models in brain.yaml (annotation.celltypist_models, per age). Only built-in
# mouse model is Developing_Mouse_Brain (P1); adult ages need custom .pkl or markers.
# dev.yaml inherits annotation+reference from brain.yaml via samples_from.
uv run python scripts/07_annotation.py --config config/dev.yaml
open results/dev/plots/07_annotation/marker_dotplot.png
open results/dev/plots/07_annotation/celltype_composition_by_group.png  # PRELIMINARY
cat results/dev/tables/annotation_summary.csv                           # fill manual_annotation
cat results/dev/tables/celltype_composition.csv                         # fractions + counts
git add scripts/07_annotation.py pyproject.toml uv.lock && git commit -m "phase 7: annotation"


# --- Phase 7b: subclustering (per cell type) ---------------------------------
# TODO — subcluster a coarse type (microglia, oligo lineage) by re-running
# HVG+scVI on the subset to resolve subtypes. Independent of 7c.


# --- Phase 7c: scANVI reference label transfer (cell type + region) ----------
# Tissue-agnostic via YAML `reference:` block. Brain -> region labels ONLY for
# region-restricted types (derived from reference: >=threshold of a type's ref
# cells in one region — no fuzzy regional claims). Placenta (region_key null) ->
# cell-type only. Reference is manual + workstation-only (multi-GB):
#   uv pip install "git+https://github.com/alleninstitute/abc_atlas_access.git"
#   # build labeled ref h5ad (cell-type + region cols on .obs), then:
#   uv run python scripts/prepare_reference.py --config config/brain.yaml \
#       --source /path/to/raw_reference.h5ad --max-cells-per-label 2000   # validates + stages
#   uv run python scripts/07c_label_transfer.py --config config/brain.yaml
# Skipped on dev (ref_h5ad null -> exits with guidance).
git add scripts/07c_label_transfer.py scripts/prepare_reference.py \
        scripts/_utils.py scripts/build_yaml.py config/ && git commit -m "phase 7c: scANVI label transfer + ref prep"


# --- Phase 8a: composition (propeller / speckle, via R subprocess) -----------
# Animal (donor) is the unit. propeller = arcsin/logit transform + limma
# moderated tests (good for small n). 2-group=t-test, 3-group=ANOVA F. Confounders
# (sex,pool) enter the limma design. NOT scCODA — its TF/arviz/numpy stack was a
# tar pit; propeller is a clean Bioconductor install (speckle+limma).
# R deps once:  R -e 'BiocManager::install(c("speckle","limma"))'
uv run python scripts/08a_composition.py --config config/dev.yaml
# --min-donors 2 to attempt n=2 groups (flagged unreliable_n<3)
open results/dev/plots/08a_composition/*/*/stacked_bar.png
open results/dev/plots/08a_composition/*/*/propeller_effects.png
cat results/dev/tables/composition_results.csv
git add scripts/08a_composition.py scripts/run_propeller.R && git commit -m "phase 8a: composition (propeller)"


# --- Phase 8b: pseudobulk DE (PyDESeq2) --------------------------------------
# Pseudobulk ONLY (never cell-level). Sum raw counts per donor per cell type.
# Design from contrast spec (~ sex + pool + group). Volcanoes label top genes
# (Ensembl->symbol via var['symbol']). flag/reliability/note carried per row.
uv run python scripts/08b_de.py --config config/dev.yaml          # --min-cells 5 on dev
# subcluster DE:  --subcluster microglia  (reads 08c_subclustered/microglia.h5ad)
open results/dev/plots/08b_de/*/*/*/volcano.png
cat results/dev/tables/de_results.csv
git add scripts/08b_de.py && git commit -m "phase 8b: pseudobulk DE + gene-labeled volcanoes"


# --- Phase 8c: pathway / GSEA (decoupler) ------------------------------------
# GSEA on DE Wald stats. Gene sets = mouse MSigDB via msigdbr (MH+M2+M5+M8),
# exported ONCE to refs/msigdb_mouse.tsv (decoupler's MSigDB-mouse fetch is broken).
# FDR corrected WITHIN each collection (kosher); FDR_pooled kept as reference.
# Plots: per-collection PANELS side by side (dotplot, volcano, celltype heatmap)
# + running-enrichment per top hit. M8 esp. useful for subclusters.
Rscript scripts/fetch_genesets.R --out refs/msigdb_mouse.tsv      # one-time
uv run python scripts/08c_pathways.py --config config/dev.yaml
# subcluster pathways:  --subcluster microglia
open results/dev/plots/08c_pathways/*/*/celltype_pathway_heatmap_panels.png
open results/dev/plots/08c_pathways/*/*/*/gsea_volcano_panels.png
cat results/dev/tables/pathway_results.csv
git add scripts/08c_pathways.py scripts/fetch_genesets.R && git commit -m "phase 8c: GSEA panels per collection"


# =============================================================================
# DEV-ONLY: composition/DE need >=2-3 donors per group, but dev has 1 sample
# per group. dev_pseudoreplicate.py SPLITS each sample's cells into N pseudo-
# donors so 8a/8b/8c code paths run. Numbers are MEANINGLESS (no real between-
# animal variance) — smoke test only. Order matters: pseudoreplicate -> 7b -> 8b.
#   uv run python dev_pseudoreplicate.py --n 3      # repo root, NOT committed/server
#   cp results/dev/h5ad/08_annotated/all_samples.h5ad.orig_backup <...>  # restore
# =============================================================================


# =============================================================================
# NEXT: Phase 8d trajectory, 8e communication (LIANA+), 8f cross-tissue,
#       8g cross-age/persistence (derived contrasts on 8b/8c tables).
#       Also TODO: subcluster annotation table (name the 7b integer subclusters).
# =============================================================================
