# Instructions for working on this project with Claude

Get broad context from the snRNAseq_project_summary.md file.

## Response style
- **Be brief.** No long preambles, no excessive caveats, no over-explaining.
- **No need to print your thoughts** unless it is helpful to either of us. This is very important, I don't want the chats to get so long because you keep on printing your thoughts
- **Don't restate what I just said.** Move to the substance.
- **Step by step.** Build one thing, verify it works, then move on. Don't write 5 files at once.
- **Be honest when something won't work** or when you're uncertain. Don't manufacture confidence.
- **Always give commands to run a specific script** (just mention where to run: Local Mac or remote WS), or rsync a specific file or folder.

## Code style
- **Parallel compute** wherever we can (training models, per-sample SoupX, etc.) — use multiprocessing or GPUs.
- **Simple > clever.** Plain Python scripts in `scripts/`, not a Python package. No Pydantic schemas, no abstract base classes, no dependency injection. Just functions and main().
- **One file per phase** (e.g. `01_validate.py`, `02_qc.py`). Each is runnable standalone.
- **Shared helpers in `scripts/_utils.py`** (leading underscore = not a phase entry point). Currently provides `load_config`, `add_lognorm`, `phase_paths`, `select_accelerator`. Add to it when something gets duplicated 2+ times.
- **Configs are plain YAML dicts.** No inheritance trees, no schema validation.
- **R is called as subprocess**, not via rpy2.
- **Idempotent steps where reasonable**, but don't over-engineer.
- **Raw counts in `.X`, lognorm computed on demand.** `04_integration_prep.py` computes the lognorm layer for Phase 5's pre-integration UMAP, then Phase 5 drops it before saving. Notebooks and downstream phases call `_utils.add_lognorm(adata)` after loading the integrated h5ad.
- **Time estimate.** Always give time estimate based on machine (WS or Mac) config.
- **Tmux command** should always be like this: tmux new -d -s soupx_placenta 'uv run python -u scripts/02_soupx.py --config config/placenta.yaml --n-jobs 6 2>&1 | tee logs/02_soupx_placenta_full.log' tail -f /home/poller/BAP-BrainPlacenta/logs/02_soupx_placenta_full.log

## Pipeline architecture decisions (don't re-litigate)
- **Language:** Python primary (Scanpy/scvi-tools). R subprocess for scDblFinder, propeller, SoupX.
- **Env:** uv + Python 3.12 on Mac (dev), workstation has GPU + R installed.
- **Subprocess R:** scripts/run-X.R called from Python, exchange via TSV/JSON.
- **Conda is blocked** at corporate firewall — don't suggest it.
- **scVI**: GPU phase, runs on workstation, not laptop.
- **Phase 1 ambient correction = SoupX (locked 2026-06-10).** CellBender abandoned 2026-06-05 (unresolvable `weakref.ReferenceType` pickle bug across all torch+pyro+numpy combos). SoupX via R subprocess is the active replacement; see §"Phase 1 = SoupX" below. CellBender Docker image was the only working option but we don't have Apptainer/Docker on the workstation.

## Dataset specifics
- **3 groups**, not 2: Early_Stress, Late_Stress, Relaxed (Relaxed = reference)
- **34 brain + 23 placenta** samples (after dropping duplicate CES2.3)
- **Ages:** brain P1/4W/3mo, placenta E12.5/E18.5
- **Pools (libraries):** Pool1-4, used as scVI batch_key
- **Known confounds** (see project doc §2): P1 Late Stress only in Pool3, placenta cross-age not comparable
- **No dam ID recorded** — treat each pup as independent observation, flag the caveat
- **Sex=TBD for all E12.5 placenta** — inferred from Y-chromosome via `01_validate.py`
- **`assigned_sex` is the source of truth for sex covariates**, not declared_sex (has unknowns) or inferred_sex (has ambiguous nulls; 10x Flex under-detects Xist). `sex_check.csv`: copies inferred_sex; replaces 'ambiguous' → 'F'; one declared/assigned mismatch (brain C5: declared F, Y-score 1.09 → assigned M, real swap to investigate). ALL downstream sex-stratified analyses MUST use `assigned_sex`. Source of truth: `sample_metadata.csv` `sex` column.

## Compute constraints
- **Laptop:** 12 GB RAM, Apple Silicon, no GPU. Dev only with subsetting (3 samples × 500 cells).
- **Workstation:** 258 GB RAM, 56 CPU cores, RTX 4500 Ada (24 GB VRAM). Production runs.
- **Network:** conda channels blocked; PyPI/CRAN/Bioconductor reachable. NVIDIA PyPI index (`https://pypi.nvidia.com`) reachable for cuML.

## Workflow conventions
- **`run_pipeline.sh` is a manual, not a history.** Tight start-to-finish walkthrough: what to run, in what order, what to inspect. NOT a changelog or bug/refactor diary (git history covers that). Self-contained — readable without opening the scripts — but terse: fragment comments over sentences, one block per phase (deps → command → key outputs → commit), one line of rationale only where a choice is non-obvious. Complete but not prolix.
- **Source of truth for samples:** `sample_metadata.csv`. YAML configs regenerate from it.
- **Dev runs:** `config/dev.yaml` (3 samples × 500 cells, ~1 GB RAM peak)
- **Full runs:** `config/brain.yaml` or `config/placenta.yaml`
- **Outputs:** `results/{tissue}/<phase>/` for prod, `results/dev/<phase>/` for laptop tests

## When asking me questions
- **One question at a time** is usually fine, **3 max**.
- **Single-select buttons** over multi-select where possible.
- **Don't ask for confirmation on small obvious things** — just do them.

## What I don't want
- Don't propose architectural pivots (containers, Nix, etc.) — we settled on uv + scripts.
- Don't add abstraction layers "for future flexibility." YAGNI.
- Don't write 200-line responses with 5 nested headers when 30 lines suffices.
- Don't apologize repeatedly when correcting something. Acknowledge briefly and move on.

## External identifiers: verify or flag
Several bugs came from writing external identifiers from memory and presenting them as verified. Root cause: plausible-looking names that don't exist. Rules:

- **Any external identifier must be verified against docs or flagged.** This covers: PyPI/conda package names, model names, gene symbols, API function names and keyword arguments, file/dataset names.
- If verified: fine, use it.
- If NOT verified: add an inline `# UNVERIFIED — check before prod` comment, OR wrap in a guard that fails loudly. Never write an unverified identifier as if confirmed.
- Cautionary examples:
  - `Mouse_Brain_Atlas` — does NOT exist; only `Developing_Mouse_Brain.pkl` ships with CellTypist.
  - `abc-atlas-access` (PyPI) — wrong; real package is `abc_atlas_access`, GitHub-only.
  - `score_genes(use_raw=False)` — runs on raw `.X`, not lognorm; need `layer="lognorm"`.
  - `Mlf1ip`/`Fam64a`/`Hn1` — outdated MGI symbols (now `Cenpu`/`Pimreg`/`Jpt1`).
  - `multi_class='ovr'` — removed in sklearn 1.7; CellTypist still hardcodes it (patched via sed; see §"CellTypist sklearn-1.7 patch").
  - `Cx3cr1` — NOT in the 10x Flex Mouse Transcriptome v2 panel; use P2ry12/Tmem119/Csf1r/Aif1 for microglia. (Also absent: Fam64a, Hmgb2, Hn1, Mlf1ip, Wdc.)
  - `datasplitter_kwargs` — correct scvi-tools 1.4.3 kwarg for num_workers/pin_memory (NOT data_loader_kwargs/dataloader_kwargs).

## No silent failures
Wrong-but-plausible output is worse than a crash, because it looks correct.

- **A correctness-critical step that can't run correctly must raise, not warn-and-continue.** Examples that now hard-fail:
  - Too few cell cycle genes match var_names → `raise` (likely Ensembl-vs-symbol mismatch).
  - No marker genes match var_names → `raise` (don't fall back to Leiden numbers).
  - `condition_cell_cycle: true` but `cc_difference` missing → `sys.exit`.
  - Required covariate column missing → `sys.exit`.
  - Phase 7 brain: any cell with `celltypist_class_predicted == "unset"` → `sys.exit` (a tier silently failed).
  - Phase 7 brain: lognorm layer missing when `apply_brain_marker_gate` runs → `sys.exit`.
  - Phase 7 brain: P1 scANVI result index ≠ P1 cells, or any P1 cell unlabeled → `sys.exit`.
- **Warn-and-skip is only acceptable when skipping an OPTIONAL output** and the skip is announced. Examples: CellTypist not installed → skip reference track; no markers → skip dotplot; an age has no CellTypist model → those cells get sentinel labels (announced).
- **Never leave NaN labels that surface as a "nan" category.** Fill with an explicit sentinel (`"no_subclass_model"`, `"no_region_model"`, `"unassigned_glia"`, etc.).
- When in doubt: fail loud and early (Phase 0 gate philosophy), not deep into a multi-hour run.

## Plots must carry biological meaning (not abstract designs)
Every figure should let a reader name the biology — the genes, cell types, or pathways affected — without cross-referencing a separate table.

- **Label the things that matter.** Volcano → name top significant genes on the plot. Composition → name cell types that shift. GSEA → name the pathways. Heatmaps → real gene/pathway/cell-type names on the axes.
- **Cap labels for readability** (e.g. top ~25 by significance); state how many more exist if truncated.
- **Gene identifiers must be human-readable.** If var_names are Ensembl IDs, map to symbols (var['symbol']) before labeling — never ship a plot of ENSMUSG IDs.
- **State the contrast and thresholds on the plot** (what-vs-what, padj/LFC cutoffs).
- Rule of thumb: if the figure can't tell you which genes/cell types/pathways are involved, it isn't done yet.

## Isolate fragile dependency stacks (don't pin them into the main env)
When a tool drags an incompatible dependency stack, isolate it — don't pin the main env backward to accommodate it.

- **scCODA abandoned** for composition (8a): TF/TFP/arviz/matplotlib/numpy/setuptools pins fought scanpy/scVI. Replaced with **propeller via R subprocess** — clean Bioconductor install, limma's empirical-Bayes moderation is better for small n anyway.
- **CellBender abandoned (2026-06-05).** torch 1.13.1/2.0.1/2.12, pyro 1.8.6, numpy <2/>=2, cellbender 0.3.0/0.3.2-master — all combos hit a `weakref.ReferenceType` pickle bug in checkpoint save (broadinstitute/CellBender #371/#386/#395). Only Docker image works; we don't have Apptainer/Docker. **Replaced with SoupX via R subprocess (2026-06-10);** see §"Phase 1 = SoupX".
- **Hard timebox: ~4 hours of debug per tool.** Past that, find an alternative or skip.

## Phase 1 = SoupX (locked 2026-06-10)
Ambient RNA correction is essential for this dataset — particularly for P1 brain, where lysed nucleated erythroblasts dump hemoglobin into the lysis buffer. (Note: the P1 erythrocyte-MISLABEL problem turned out to be a REFERENCE issue, not ambient — see Annotation conventions. SoupX correctly stripped Hb; the fix was switching P1 to the Rosenberg reference. SoupX is still essential for ambient cleanup generally.)

- **Why erythrocytes can appear at all in snRNA-seq:** Mature mammalian RBCs are anucleate, BUT P1 mouse brain still contains nucleated erythroblasts (basophilic / polychromatic / orthochromatic stages) in residual vasculature. A small number of real erythroid calls at P1 is biologically plausible.
- **Tool:** SoupX (CRAN), R subprocess. Avoids CellBender's pickle bug entirely.
- **Workflow:** cellranger filtered + raw counts → `SoupChannel` → `scran::quickCluster` → `setClusters` → `autoEstCont` (data-driven rho per cluster) → `adjustCounts` → corrected counts written as MTX + barcodes.tsv + features.tsv + soupx_summary.json. Python orchestrator assembles per-sample h5ad.
- **Manual rho fallback** (`--rho 0.10`) if `scran` install fails — bypasses clustering, uses fixed contamination fraction. ~10% is reasonable for snRNA-seq brain.
- **Scripts:** `scripts/run_soupx.R` (R worker) + `scripts/02_soupx.py` (Python orchestrator; parallel via `ProcessPoolExecutor`, default `--n-jobs 4`, RAM ~5-15 GB per concurrent sample).
- **Output:** `results/{tissue}/h5ad/02_soupx_corrected/{sample_id}.h5ad` + `tables/02_soupx/02_soupx_summary.csv` (per-sample rho_mean, pct_removed, n_cells).
- **Smoke test:** `--sample-ids E1 --n-jobs 1` on one sample before launching the full ~57-sample production run.
- **Downstream wiring:** `02_qc.py` has a prefer-soupx fallback — if `02_soupx_corrected/{id}.h5ad` exists, read that; else fall back to cellranger filtered h5. SoupX changes counts → invalidates Phase 2 onwards; full re-run from Phase 2 needed.

## renv Suggests workaround (2026-06-10)
renv installs `Suggests` dependencies by default — for SoupX this pulls Seurat → shiny → bslib → fs, and `fs` needs `libuv-dev` system library. Three layered fixes:

1. **Project-level setting (preferred, persistent):**
   ```r
   renv::settings$package.dependency.fields(
     c("Depends", "Imports", "LinkingTo"), persist = TRUE
   )
   ```
   Writes to `renv/settings.dcf`. All future installs skip Suggests in this project.
2. **Per-call argument** (only honored if project-level setting is not overriding): `renv::install("X", dependencies = c("Depends","Imports","LinkingTo"))`.
3. **Bundled libuv fallback:** `USE_BUNDLED_LIBUV=1 Rscript -e '...'` — even if Suggests sneak in, `fs` builds without the system library.
4. **Last resort (one-shot, sudo):** `sudo apt-get install -y libuv1-dev`.

## CellTypist sklearn-1.7 patch (2026-06-10)
CellTypist's `train.py` (lines 126 and 146) hardcodes `multi_class='ovr'` in `LogisticRegression(...)`, which sklearn 1.7+ removed. Two-pattern sed fixes it:

```bash
sed -i "s/multi_class = 'ovr', //;s/multi_class = 'ovr'//" \
  /home/poller/BAP-BrainPlacenta/.venv/lib/python3.12/site-packages/celltypist/train.py
```

- Removing this argument *improves* calibration: sklearn 1.7+ defaults to true multinomial softmax (was one-vs-rest in older versions).
- **Tech debt:** patch lives in `.venv/`; `uv sync` reverts it. Re-apply via `./setup-remote.sh --skip-references --skip-r --skip-cellbender` (Step 2.5 is idempotent).
- This patch is **required before CellTypist training** in any env using sklearn ≥1.7.

## CellTypist GPU training via cuML (locked 2026-06-10)
- **Install:** `cuml-cu12 cudf-cu12` from `https://pypi.nvidia.com`. Downgrades: numba 0.65→0.64, pyarrow 24→23, cuda-toolkit 13→12.9 (CUDA 13 driver is backward-compatible with 12.9 runtime).
- **Usage:** pass `use_GPU=True` to `celltypist.train()` (valid kwarg, default False).
- **Speedup observed:** ~9× for class (34 labels: CPU 123 min → GPU 14 min), >40× for subclass (334 labels: CPU >18 hours, never finished → GPU 27 min), ~3× for region (12 labels: CPU est. 60 min → GPU 19 min). Class+subclass+region full retrain ~60 min total on RTX 4500 Ada.
- **Output:** cuML LogReg fit, saved as `sklearn.linear_model._logistic.LogisticRegression` (CellTypist's design) — pkls are interchangeable with CPU-trained versions for inference.
- **L-BFGS line-search warning** on small problems is a problem-size issue, not a GPU issue. At full-data scale (92K cells × 12K features post-FS), L2 regularization (C=1.0) makes the optimization well-posed even when params > samples (subclass: 4M params, 92K samples).
- **Why CPU L-BFGS hung on subclass:** L-BFGS is not parallelizable across classes (each iteration computes one gradient over the full multinomial likelihood). Only inner BLAS matmul parallelizes. 334 classes × 92K cells single-threaded = days.

## Phase 8 conventions
- **Statistical unit is the animal (donor_id)**, never the cell. Composition = per-donor cell-type counts; DE = pseudobulk (sum raw counts per donor). No dam ID → each pup independent (anti-conservative; caveat carried in `flag`).
- **Corrected p-values everywhere** — propeller FDR (8a), DESeq2 padj (8b), per-collection BH FDR (8c). No raw p-value drives any significance call.
- **GSEA gene sets = mouse MSigDB via msigdbr** (MH+M2+M5+M8), exported once to refs/msigdb_mouse.tsv. NOT decoupler's get_resource(MSigDB,mouse) — it's broken. FDR corrected WITHIN each collection (sizes differ ~150x); pooled kept as reference.
- **Multi-database plots = side-by-side panels**, one subplot per collection, figure width scaled to #collections. Use constrained_layout (NOT with bbox_inches='tight' — they fight and clip).
- **--subcluster flag** on 8b/8c/8e runs on the 7b subcluster object (label col `subcluster_name` from 7d, else `subcluster` integer from 7b), writing `*_subcluster_{slug}` outputs to separate folders.
- **Subcluster runs are a LOOP, not one cell type.** Production runs every focal cell type through 8b, 8c, and 8e via the same `CELL_TYPES` array defined once in run_pipeline*.sh.
- **Subcluster naming (7d)**: 7b produces INTEGER subcluster ids; 7d names them via CellTypist majority_voting + literature-marker scoring from `config/subcluster_markers.yaml`, aggregated per cluster (mean score → argmax).
- **TF activity (8c) = REQUIRED in production**: ULM on DE Wald stats vs CollecTRI mouse network, per contrast×celltype, BH-FDR within celltype×contrast. Always pass `--tf`. Gates 8f view 5 and 8g view 3; cannot be recovered without re-running 8c. Needs network (omnipath).
- **8e cell-cell communication = LIANA+ in main env, no sidecar.** Three arms in one script: baseline `rank_aggregate` per group×age, differential via `li.multi.df_to_lr` reading 8b's Wald stats, per-donor for statistics. Covers all three group comparisons (ES-v-Rel, LS-v-Rel, ES-v-LS).
- **8f cross-tissue = six views, all reproducible from 8b/8c CSVs.** Two biologically aligned arms (E12.5 placenta Early → P1/4W/3mo brain Early; E18.5 placenta Late → same; P1 Late carries `confounded_with_pool` flag). Views: DEG overlap (hypergeom), RRHO (custom NumPy), pathway concordance from 8c GSEA, LR cross-tissue mechanistic hypotheses (placental ligand × brain receptor from liana mouseconsensus with `stress_axis` flag column), TF concordance from 8c TF activity, ORA of overlap genes vs MSigDB. The LR table is the publication-quality output.
- **NO cross-tissue cell-cell communication.** BBB makes literal placenta-cell-to-brain-cell signalling implausible — 8f view 4 (LR from DE) is the correctly framed endocrine/paracrine version.
- **8g cross-age persistence = brain only.** Placenta has incomplete cross-age factorial; 8g exits cleanly with `tissue: placenta`. Persistence classes: persistent, resolving_early, established_late, P1_only, transient_4W, emergent_3mo, P1_3mo_only, persistent_directionswap. Cross-arm core signature (view 6) = intersection of persistent calls in BOTH arms with consistent direction = paper-quality table.

## Phase 9 — two scientific arms (locked 2026-06-05)
Cross-species RRHO2 validation runs as TWO arms reported separately, NOT pooled.

- **ARM A — psychiatric/neurodevelopmental.** Nagy 2020 (GSE144136, MDD M), Maitra 2023 (GSE213982, MDD F+M, Mic1 = 38% of female-MDD DEGs), Velmeshev 2019 (UCSC autism), Herring 2022 (GSE168408, developmental PFC), Marsh 2022 (GSE198373, placenta).
- **ARM B — MS as stressed-cell signature reference.** Macnair 2025 (Zenodo 10.5281/zenodo.8338963, 632K nuclei), Absinta 2021 (GSE180759, MIMS-iron/MIMS-foamy), Jäkel 2019 (GSE118257, Oligo5/Oligo6).
- **CRITICAL framing for ARM B: NOT etiology.** Norwegian and Danish/Swedish registry data implicate maternal metabolic / nutritional / exposure variables for MS risk, not psychological stress. Valid claim: "mouse prenatal-stress microglia and OL share transcriptional features with disease-associated microglia/OL described in MS, consistent with a shared stressed-glia program."
- **Schirmer 2019 (PRJNA544731) DEFERRED.** Raw FASTQ only on SRA.
- **Mouse anchors DEFERRED** (Marques 2016, Falcão 2018, Velmeshev 2023, Braun 2023).
- **Subset RRHO is the paper-quality comparison.** Use `subset_labels` block in `config/cross_species_celltype_map.yaml` for headline figures.
- **Cross-species region matching uses `celltypist_region`** (brain). P1 region is now populated (parsed from Rosenberg labels), so P1 can participate in region-matched comparisons — not just adults.
- Loaders in `09_cross_species_validation.py` are STUBS — raise `NotImplementedError`. Smoke-test on Velmeshev first.

## Smoke-test policy
For every phase that takes >10 min on workstation: build a 1-sample (or 1-pool, 1-cluster) subset, run, verify outputs in expected paths, THEN launch full thing in tmux.
- Burned ~4 hours on CellBender's checkpoint bug at scale.
- Burned 30 min on missing `brain.yaml` CellTypist mappings.
- **Phases that NEED smoke tests:** 1 SoupX, 5 scVI, 7 P1-scANVI (Rosenberg transfer), 8c with --tf, 9 cross-species.
- **Phases too cheap to need smoke tests:** 0, 2, 3, 4, 6, 7 main (adult CellTypist), 7b, 7d, 8a, 8b, 8d, 8e (under ~30 min wall time per tissue).

## Dev workflow — split at h5 level, not at runtime
- **dev_split_h5.py** runs ONCE before Phase 0 on dev. Reads the 3 dev h5 files, writes 9 split h5 files (random barcode partition) into `data/dev_split/`, and emits `config/dev_split.yaml`.
- **No pipeline scripts are dev-aware.** All phase scripts run unchanged with `--config config/dev_split.yaml`.
- Pseudo-donors are random cell partitions → numbers MEANINGLESS, smoke test of code paths only.
- **Dev is 4W-only by design.** All 9 pseudo-donors are 4W M Pool1. 8g cannot be exercised meaningfully on dev (one age → every classification = `transient_4W` or `none`). Note: dev being 4W-only means the P1-scANVI branch is NOT exercised on dev — smoke-test it on real P1 data on the WS.

## Annotation conventions
- **Phase 7 uses per-cluster majority voting** (CellTypist convention), not per-cell argmax. Cells in one Leiden cluster share a label; low-purity (<60% majority OR runner-up >25%) clusters announced for manual review.
- **Phase 7d (subcluster naming) is already cluster-level** — scores aggregate per integer subcluster ID.
- **Brain Phase 7 = 4-tier, P1 via scANVI (UPDATED 2026-06-12; supersedes the 3-tier/Di-Bella scheme):**
  - **broad** (`celltypist_broad`, region-FREE, cross-age tier): derived from class per-age — ABC `class_to_broad_csv` for 4W/3mo, `config/rosenberg_class_to_broad.csv` for P1; trailing ` (region)` suffix stripped so all ages align. ~9 classes. Derived in Phase 7 by `derive_brain_broad`. THIS is the cross-age tier (8g persistence, 9 cross-species).
  - **class** (`celltypist_class`, canonical 8b/8c key): per-(Leiden cluster × age) MAJORITY vote. Region-TAGGED. ABC 34-vocab for 4W/3mo; Rosenberg ~18-vocab (e.g. CTX Glut, CB GABA) for P1. Two vocabularies coexist by design (per-age native).
  - **subclass** (`celltypist_subclass`, per-cell): ABC 334 for 4W/3mo; raw Rosenberg ~65 fine labels for P1.
  - **region** (`celltypist_region`, per-cell): ABC 12 for 4W/3mo; parsed from Rosenberg label prefix for P1 (CTX/CB/TH/HPF/OLF/STR/MB/non-regional).
  - **P1 = scANVI label transfer from Rosenberg 2018 P2-brain** (`run_scanvi_p1.py`, GPU subprocess called inside `07_annotation.py` for the P1 branch only — mirrors the R-subprocess pattern). Reference built by `prepare_rosenberg_reference.py` (→ `refs/rosenberg_p2brain_reference.h5ad` + 3 grouping CSVs in `config/`). Di Bella ABANDONED: cortex-only, mislabeled 42% of whole-brain P1 as erythrocyte (region-coverage failure, NOT ambient — SoupX had stripped Hb). 4W/3mo = CellTypist from ABC. `07c_label_transfer.py` DELETED.
  - YAML: `annotation.scanvi_p1.{ref_h5ad, labels_key, config_dir}` for P1; `annotation.celltypist_models.<age>.{class,subclass,region}` for adults (P1 entry null).
  - The gate/majority/audit path is SHARED across all ages — P1 differs only in where the per-cell label comes from (scANVI vs CellTypist); everything downstream is common code.
- **Placenta Phase 7 has no CellTypist model.** Curated literature markers + STAMP Spearman correlation against Liu 2024 reference (35 cell types, E9.5-E18.5). Tier 1+2 label-collapse (35 → ~21) + STRICT canonical-marker gates (Neutrophil / Lymphoid / Megakaryocyte) + Xist + Y-gene compartment scoring + EPC/TSC negative-control QC. Canonical key `celltype_majority`. NO broad tier yet — add a compartment grouping (trophoblast/decidua/immune/vascular/erythroid) when Phase 8f cross-tissue needs it.
- **⚠ Brain plot-vs-data provenance (2026-06-12):** Phase 7 brain ran ONCE; `celltypist_broad` was then patched region-free in place via `scripts/patch_broad_regionfree.py` (no Phase 7 re-run), and broad + per-age-class UMAPs replotted via `scripts/replot_brain_annotation.py`. The h5ad is authoritative; not every plot in `plots/07_annotation/` came from one clean run. A clean `07_annotation.py` re-run reconciles everything. The annotation/plotting SPLIT was discussed but DEFERRED (no time) — for now plots are still coupled to the annotation run.

## Brain marker gate (updated 2026-06-12)
STRICT canonical-marker gates for borderline brain calls. Runs over ALL (cluster×age) rows identically (P1 included). CellTypist/scANVI conf says "how sure among trained classes" — can't verify the cell expresses the biology the label requires. Gate catches false positives.

- **`BRAIN_GATE_CONFIG` in `scripts/07_annotation.py`** — five gated types:
  - microglia: ≥2 of {P2ry12, Tmem119, Csf1r, Aif1} → `unassigned_immune`. **Cx3cr1 REMOVED — absent from 10x Flex Mouse Transcriptome v2 panel.**
  - astrocyte: ≥2 of {Aqp4, Gja1, Slc1a3, Aldh1l1} → `unassigned_glia`
  - ol_lineage: ≥1 of {Mbp, Mog, Plp1, Mag, Pdgfra, Cspg4, Olig1, Olig2, Sox10} → `unassigned_glia`. **OPC markers ADDED — OPCs lack mature myelin genes, were wrongly demoted (P1 cluster 48, 11k cells).**
  - endothelial: ≥2 of {Cldn5, Pecam1, Cdh5} → `unassigned_vascular`
  - **erythroid (ADDED): ≥2 of {Hbb-bs, Hbb-bt, Hba-a1, Hba-a2, Alas2} → `unassigned_erythroid`.** Safety net, self-targeting: real placenta erythroid (Hb high, 91% express ≥2) pass; false brain calls (Hb=background) demote. With scANVI/Rosenberg, P1 has no erythroid labels → 0 demotions; gate is belt-and-suspenders.
- **`MARKER_PRESENCE_THRESHOLD = 0.20`** — marker "present" if ≥20% of cells in (cluster, age) have lognorm > 0.
- **Keyword matching is case-insensitive substring** on `winner_class`.
- **Audit CSV columns:** `markers_checked`, `markers_present`, `gate_outcome` (`no_gate`|`passed`|`demoted`), `gate_label`.
- **Production run 2026-06-12:** 4 demotions (Astro-Epen, weak astro markers), 0 erythroid demotions.

## Brain age-composition sanity (added 2026-06-10)
Diagnostic CSV `tables/07_annotation/07_annotation_age_composition_sanity.csv` flags developmentally-implausible (cluster × age) rows. **Informational only — does not modify labels.**

- `BRAIN_AGE_EXPECTATIONS` list in `scripts/07_annotation.py` — keyword + expected_ages + flag_name tuples.
- Current rules: radial glia / intermediate progenitor / IPC / neuroblast / glioblast / erythrocyte / erythroid progenitor → expected only at P1.
- Production run 2026-06-12: 0 flags (good).

## Plot format strategy (locked 2026-06-05)
- **Default: PNG @ 300 DPI** for all pipeline plots. Cell-level UMAPs with 600K dots aren't suitable for pure SVG.
- **Paper figures only** (~5-10 plots): refactor to PNG + PDF hybrid via `_utils.savefig(fig, path, dpi=300)`.
- Don't refactor all 13 plotting scripts. Targeted post-Phase 8 update only.
- **on-data UMAP labels require categorical dtype** — cast the color column to `category` before `sc.pl.umap(..., legend_loc="on data")`, else scanpy silently draws nothing.

## UMAP determinism (locked 2026-06-05)
- **Explicit `random_state`** in Phase 5 + Phase 6 (scanpy default has no seed).
- **Phase 5b UMAP seed sweep** (`05b_umap_sweep.py`, to be written): 5 seeds (42, 0, 7, 123, 2024), scVI training not re-run — just reprojections from saved latent.
- UMAP hyperparameters otherwise locked: `n_neighbors=15`, `min_dist=0.5`, `spread=1.0`, `init_pos='spectral'`, `metric='euclidean'`.

## Trajectory (8d) — no velocity, all ages equal
- **NO RNA velocity, NO CellRank.** 10x Flex is probe-based (exon-only), can't resolve spliced/unspliced. Without velocity CellRank only duplicates PAGA.
- **All ages treated identically in DPT** — pooled and per-age both run; age-split rows carry a `pool_age_confound` caveat in 'note' (not a validity gate).
- **PAGA edges are hypotheses, not transitions.** `trajectory_paga_edge_diagnostics.csv` flags ambient-driven / doublet-driven / shared-gene edges.

## Output organisation
- **Tables in per-phase subfolders**: `tables/<phase_dir>/<phase>_<name>.csv`. Helper: `_utils.phase_table_dir(cfg, label)`.
- **Plots also in per-phase subfolders**: `plots/<phase_dir>/...`.
- **Filenames carry the phase prefix** for identifiability when copied out.
- **Subcluster runs go to suffixed folders** (`plots/08e_communication_subcluster_excitatory_neurons/...`).

## Offline-audit CSVs (no workstation access post-run)
- **02_soupx** `02_soupx_summary.csv` — per-sample rho_mean, rho_min, rho_max, pct_removed, n_cells, n_clusters, elapsed_sec, mode (autoEst | manual). Lets you spot-check whether SoupX over-corrected any sample.
- **07** `07_annotation_class_per_cluster_age.csv` — augmented with `markers_checked`, `markers_present`, `gate_outcome`, `gate_label`. Reviewers can verify any gated label without re-running.
- **07** `07_annotation_age_composition_sanity.csv` — developmentally-implausible (cluster × age) rows.
- **8b** `08b_de_gene_expression_per_sample.csv` — per-sample mean lognorm of DE genes, keyed celltype × gene × sample_id.
- **8c** `08c_pathway_leading_edge.csv` — genes driving each significant pathway with log2FC + direction.
- **8c** `08c_tf_activity.csv` — TF activity scores + FDR per contrast×celltype×TF.
- **8d** `08d_trajectory_paga_edge_diagnostics.csv` — per cell-type-pair edge audit.
- **8e** `08e_lr_quantified.csv` — per-donor LR scores backing the Mann-Whitney comparisons.
- **8f** `08f_lr_cross_tissue.csv` — placental ligand × brain receptor mechanism hypotheses with `stress_axis` flag column. Publication-quality output.
- **8g** `08g_core_signature_genes.csv` — features persistent in BOTH stress arms with consistent direction. Paper-quality "core stress signature".

## Ask before strong scientific calls
- Don't drop an analysis, exclude samples/ages, or add complexity (extra tools, sidecar venvs) without checking it earns its place for THIS dataset. Surface the question; don't bake the decision in silently.
- Specifically don't propose: cross-tissue cell-cell communication (BBB makes it implausible — 8f view 4 is the correctly framed endocrine version); RNA velocity / CellRank (10x Flex probe-based, no spliced/unspliced); scCODA (dependency stack fights scanpy; propeller replaces it).

## Always state where files go
When presenting files in a Claude response, ALWAYS say which directory each one goes into (scripts/, config/, refs/, repo root, etc.). Don't leave the user to guess.

## Workstation infrastructure
- **SSH:** `ssh poller@172.17.213.147` (from Mac). User: `poller`.
- **Local repo on Mac:** `/Users/shuppar/Downloads/BAP_data_1/Analysis/`
- **Workstation project root:** `/home/poller/BAP-BrainPlacenta/` (NVMe).
- **Raw tars on workstation:** `/media/poller/PollerLab-1/BAP-data1/processed_260411_Shiv_FLEX_*.tar` (USB-HDD).
- **Extracted Cell Ranger output:** `/media/poller/PollerLab-1/BAP-data1/Analysis/data/Pool{1,2,3,4}/per_sample_outs/<sample_id>/` (USB-HDD). Reached from the project via the symlink `BAP-BrainPlacenta/data` → that path.
- **Disk layout rationale:** project on NVMe (fast intermediate h5ad writes); raw data on USB-HDD (read once per phase). Symlink bridges.
- **Pull-from-Mac convention:** code edits happen locally and rsync to the WS. Exclude list: `results/`, `data/`, `.venv/`, `.venv-cellbender/`, `__pycache__/`, `.git/`, `*.h5ad`, `logs/`, `.DS_Store`. Flags: `-av --progress --chmod=Fu+x`.
- **WS ↔ Mac MUST stay mirrored.** Any WS-side edit (sed, in-place script writes, one-liner config tweaks) must be rsync'd back to Mac the same turn. Reverse rsync:
  ```bash
  rsync -av --chmod=Fu+x \
    poller@172.17.213.147:/home/poller/BAP-BrainPlacenta/<path> \
    /Users/shuppar/Downloads/BAP_data_1/Analysis/<path>
  ```
- **tmux for any multi-minute job.** Detach `Ctrl-b d`, reattach `tmux attach -t <name>`, list `tmux ls`.
- **Pool → contents map:**
  - Pool1: 16 brain (3mo all groups + 4W males)
  - Pool2: 16 brain (4W females + P1 Early/Relaxed)
  - Pool3: 2 brain P1 Late + 14 placenta E12.5 (incl. duplicate CES2.3 to drop)
  - Pool4: 10 placenta (2 E12.5 Relaxed + all E18.5)
- **GPU operational hygiene:** pre-flight `nvidia-smi --query-gpu=memory.used` check before launching scVI/scANVI; refuse if non-display VRAM > 2 GB. `PYTORCH_CUDA_ALLOC_CONF=max_split_size_mb:512` for long scVI runs. Explicit cleanup between GPU phases: `del model; torch.cuda.empty_cache(); gc.collect()`.
