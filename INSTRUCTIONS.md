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
- **Parallel compute** wherever we can (training models, per-sample SoupX, per-slice propeller, etc.) — use multiprocessing or GPUs. See the "Parallelism is mandatory" section below — the standard is `_utils.parallel_map`.
- **Simple > clever.** Plain Python scripts in `scripts/`, not a Python package. No Pydantic schemas, no abstract base classes, no dependency injection. Just functions and main().
- **One file per phase** (e.g. `01_validate.py`, `02_qc.py`). Each is runnable standalone.
- **Shared helpers in `scripts/_utils.py`** (leading underscore = not a phase entry point). Currently provides `load_config`, `load_contrasts`, `phase_paths`, `phase_table_dir`, `add_lognorm`, `select_accelerator`, `iter_strata`, `parallel_map`, `unassigned_mask`. Add to it when something gets duplicated 2+ times.
- **Configs are plain YAML dicts.** No inheritance trees, no schema validation.
- **R is called as subprocess**, not via rpy2.
- **Idempotent steps where reasonable**, but don't over-engineer.
- **Raw counts in `.X`, lognorm computed on demand.** `04_integration_prep.py` computes the lognorm layer for Phase 5's pre-integration UMAP, then Phase 5 drops it before saving. Notebooks and downstream phases call `_utils.add_lognorm(adata)` after loading the integrated h5ad.
- **Time estimate.** Always give time estimate based on machine (WS or Mac) config.
- **Tmux command** should always be like this: tmux new -d -s soupx_placenta 'uv run python -u scripts/02_soupx.py --config config/placenta.yaml --n-jobs 6 2>&1 | tee logs/02_soupx_placenta_full.log' tail -f /home/poller/BAP-BrainPlacenta/logs/02_soupx_placenta_full.log
- **rsync command format**: Example -- rsync -av --chmod=Fu+x \
  /Users/shuppar/Downloads/BAP_data_1/Analysis/scripts/replot_brain_annotation.py \
  poller@172.17.213.147:/home/poller/BAP-BrainPlacenta/scripts/ (note: we are saving the WS results on MAC in appropriate subfolders in /home/poller/BAP-BrainPlacenta/results_WS (and not results)).

## Parallelism is mandatory for repeated work (not optional)
Any phase script that loops over samples / jobs / contrasts / cell-types and, per item,
launches a subprocess (R worker) or calls a heavy function MUST parallelize via
`_utils.parallel_map` and expose a `--n-jobs` CLI flag. A bare for-loop that starts one
subprocess per item is a performance bug — treat it like any other bug.

- `for item, result, err in parallel_map(fn, items, n_jobs=args.n_jobs, desc="..."):`
  yields (item, result, error); `err` is a string on failure (captured, not raised), so
  one bad item never aborts the batch — the caller decides what to do with it.
- Default `use_threads=True` (subprocess / I/O-bound, e.g. R workers). Use
  `use_threads=False` (processes) only for CPU-bound pure-Python work, and then `fn` and
  items must be picklable (top-level fn, no closures/lambdas).
- Default `--n-jobs` 8; on the workstation pass 16–24 for light R workers.
- Build the work list first (cheap pandas/IO), then hand it to parallel_map — don't
  interleave heavy serial setup inside the parallel loop.
- Reference implementation: `08a_composition.py` (collects per-slice propeller jobs, then
  runs them concurrently). Serial was ~3 h; parallel is single-digit minutes.
- **CPU-bound exception:** the shuffle null in `08b_disruption_shuffle_test.py` uses
  `use_threads=False` (process pool) — numpy permutation loops are CPU-bound pure Python
  and the GIL would serialize threads. The worker is a top-level function for picklability.

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
- **Significance must be obvious**, not a subtle mark. 8a heatmaps outline FDR<0.05 cells; volcanoes star/label significant genes. Rule of thumb: if the figure can't tell you which genes/cell types/pathways are involved, it isn't done yet.
- **Forest vs bar plots:** the 8b disruption mirror plot uses bars (showing MAGNITUDE of LOST/GAINED counts); the 8b shuffle-test plot uses bars for `|Δ|` with labels ALWAYS placed OUTSIDE the bars (never overlapping). The two visualizations answer different questions — bars in the disruption plot describe what was observed; the shuffle-test bars + within-stratum panel quantify whether the observation deviates from the k-preserving null.

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
- **Scripts:** `scripts/run_soupx.R` (R worker) + `scripts/02_soupx.py` (Python orchestrator; parallel via `_utils.parallel_map`, default `--n-jobs 4`, RAM ~5-15 GB per concurrent sample).
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
- **Speedup observed:** ~9× for class (34 labels: CPU 123 min → GPU 14 min), >40× for subclass (334 labels: CPU >18 hours, never finished → GPU 27 min), ~3× for region (12 labels). Class+subclass+region full retrain ~60 min total on RTX 4500 Ada.
- **Output:** cuML LogReg fit, saved as `sklearn.linear_model._logistic.LogisticRegression` (CellTypist's design) — pkls are interchangeable with CPU-trained versions for inference.
- **Why CPU L-BFGS hung on subclass:** L-BFGS is not parallelizable across classes (each iteration computes one gradient over the full multinomial likelihood). 334 classes × 92K cells single-threaded = days.

## Phase 8 conventions

**Cross-cutting (8a done; 8b–8g follow the same — wire via the shared helpers):**
- **Statistical unit is the animal (donor_id), never the cell.** Composition = per-donor cell-type counts; DE = pseudobulk (sum raw counts per donor). No dam ID → each pup independent (anti-conservative; caveat carried in `flag`).
- **Contrasts are declarative** (`config … contrasts:`, via `_utils.load_contrasts(cfg, kind="de")`) and shared across every 8x stage — never hard-code or synthesize a contrast in a phase script. The set: `early_vs_relaxed_per_age`, `late_vs_relaxed_per_age`, `omnibus_3group_per_age`, `early_vs_late_per_age` (brain only — it lives in the YAML), `within_group_across_age` (DE-only), `group_x_age_interaction` (DE-only). A phase consumes the contrasts it can handle and skips the rest with an announcement (8a skips across-age / interaction / sex-stratified contrasts).
- **Sex strata are declarative** (`config … strata: {sex: [combined, M, F]}`, via `_utils.iter_strata`). Applied to EVERY contrast in EVERY stage; `combined` = sex stays a covariate, `M`/`F` = subset the donors (sex then auto-drops from the design). Every output row/path carries a `sex` column. This SUPERSEDES the old `within_age_sex_stratified` contrast — remove that from `build_yaml.py` as part of wiring 8b. A stratum/group with only 1 donor is unavoidably skipped (e.g. P1 Relaxed females).
- **Drop, don't reassign, non-cell-types** — from numerator AND denominator, once up front: contaminants (`subcluster_name` starting `Contamination_` or `=="unresolved"`, from the 08c subcluster objects) and `unassigned*` gate labels (`_utils.unassigned_mask`). Record the dropped per-donor counts/fractions in a diagnostic CSV (`08a_dropped_cells_per_donor.csv` pattern) — never lose them silently. Plots show ONLY real, assigned cell types; the dropped mass lives in the diagnostic table.
- **Design `~ sex + pool + group`; drop any covariate that is constant OR perfectly aliased with `group` in a slice, and flag `confounded_with_pool`.** Canonical case: P1 `Late_Stress` is Pool3-only → pool ≡ group → rank-deficient → drop `pool`. Use the `aliased_with(df, factor, cov)` pattern from 8a. scVI batch correction does NOT touch count-level tests, so `pool` stays in the design wherever it's estimable.
- **`min_donors=2` to run a stratum; any group `<3` → `reliability=low_n`** (`config … composition: {min_donors: 2, reliable_donors: 3}`). Trust `ok` rows with a finite effect; `low_n` rows with NaN/inf effect are degenerate (a rare type absent in a tiny group) and inflate the raw FDR<0.05 count — read `reliability==ok` + finite effect for real hits.
- **Levels & granularities:** brain has region levels (`celltypist_region`) + `whole`; placenta `whole` only. Tiers: brain `celltypist_broad` (broad) + `celltypist_class` (class, canonical 8b/8c key) + subtype (focal coarse types exploded to `subcluster_name`); placenta `celltype_majority` + subtype. Fraction within a slice = cell type ÷ that donor's total cells in the slice (region → ÷ donor's cells in that region; whole → ÷ donor's total cells).

**Stage-specific:**
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

## Phase 8b follow-ups (locked 2026-06-15)
Three brain-only scripts operate on the master `08b_de_results.csv` (no DE re-run; the master CSV is already a deterministic function of the integrated h5ad + contrasts YAML). All three skip placenta cleanly — placenta has no `within_group_across_age` contrast (incomplete cross-age factorial; see project doc §2).

**`08b_developmental_disruption.py`** — classifies genes from the `within_group_across_age` contrast (per group, with all three pairwise age tests collapsed to the most-significant pair per gene) into 5 direction classes:
- `universal` — sig age-DE in all 3 groups (Relaxed, Early, Late). Developmental baseline.
- `relaxed_only` (`= LOST`) — sig only in Relaxed. Trajectory present in controls, lost under both stress regimens.
- `stress_shared` (`= GAINED`) — sig in BOTH Early AND Late, NOT in Relaxed. Trajectory absent in controls, induced by both stress regimens.
- `early_only` / `late_only` — sig in one stress group only.
Outputs `08b_developmental_disruption_summary.csv` (per sex × level × celltype: counts + mean `|LFC|` per group for the LOST class) + `08b_developmental_disruption_genes.csv` (long-form gene-level direction class assignments).

**`08b_followup_plots.py`** — two plot types per `(sex × level)`:
- Mirror disruption bar (LOST left red, GAINED right blue, cell types named on Y) + paired `|LFC|` boxplots showing effect-size collapse for LOST-trajectory genes under stress (Relaxed gray / Early red / Late blue).
- Stress-consistency stacked bars per age — Early-only (red) / Both-sig=convergent (gray) / Late-only (blue). Auto-skips when one of two stress contrasts is missing.

**`08b_disruption_shuffle_test.py`** — the k-preserving permutation null. Reports both permutation-based and analytic-binomial significance.
- **k-preserving null:** per gene, keep `k_i` = #groups in which it's sig (0/1/2/3) but randomize WHICH groups. Vectorized via per-row argsort of random scores on a `(n_genes, 3)` array. ~1 ms per shuffle on ~20K genes; 1000 shuffles per slice in single-digit seconds.
- **Per-category analytic binomial p-values** for all 6 disjoint sig-pattern categories (R-only / E-only / L-only / R∩E / R∩L / E∩L), modelling `obs ~ Binom(n_k_stratum, 1/3)`. Enrichment and depletion p-values reported separately; BH-corrected within each direction.
- **Within-stratum chi-square goodness-of-fit:** tests whether the three k=1 categories (R-only / E-only / L-only) are uniformly distributed and same for the three k=2 categories. Directly addresses "is Relaxed special?"
- **Output:** `08b_disruption_shuffle_test.csv` with columns: 6 obs counts, `n_k0..n_k3`, 6 categories × 4 p-values (enrichment + depletion + BH for each), `chi2_k1` + `p_chi2_k1` + `chi2_k2` + `p_chi2_k2`, permutation-null sanity check on the LOST/GAINED diff (`null_lost_mean`, `null_gained_mean`, `null_lost_p5/p95`, etc.). Plus the headline two-panel figure: Panel A = mirror `|Δ|` bar (LOST left, GAINED right; solid color for `Δ>0` enrichment, faded for `Δ<0` depletion; labels ALWAYS placed OUTSIDE the bars); Panel B = within-stratum 6-bar breakdown per cell type with dashed reference lines at `n_k1/3` and `n_k2/3`.

## Disruption analysis framing (locked 2026-06-15)
The headline biological finding for brain `within_group_across_age` is NOT just "LOST > GAINED" — that asymmetry is partly expected from marginal sig-rate maths (the AND requirement on GAINED). The cleanest disruption claim that survives the shuffle test is:

> "When age-DE signal is shared across two groups, Relaxed is almost always one of them. R∩E and R∩L are both massively enriched (↑***) over the k-preserving null in every brain broad cell type, while E∩L (=GAINED) is correspondingly depleted (↓***). Stress conditions independently lose age-trajectory genes from the Relaxed program rather than coordinately gaining new ones."

Supporting evidence at the broad-cell-type level (sex=combined, level=whole, brain main):
- 3 broad cell types significantly enriched LOST (`p_lost_BH < 0.05`): Olfactory ensheathing (Δ=+165 ↑**), Astrocytes/Ependymal (Δ=+130 ↑**), Vascular (Δ=+37 ↑*).
- ALL 7 broad cell types significantly depleted GAINED (`p_gained_dep_BH < 0.05`): E∩L counts fall to ~1/3 of null expectation.
- ALL 7 broad cell types massively enriched R∩E and R∩L over null (binomial enrichment p < 1e-3 for these k=2 categories).

Supporting evidence at the subcluster level: PAM_ATM_Microglia (Δ_LOST=+147 ↑**), BAM (Δ=+92 ↑**), Homeostatic_Microglia (Δ=+9 ↑**), Protoplasmic_Astrocyte (Δ=+334 ↑**), OPC (Δ=+72 ↑*), MFOL (Δ=+17 ↑*).

Methods caveat to include in paper: `within_group_across_age` is pool-confounded by design (each age uses different pools). The k-preserving shuffle test is partially robust to this — each gene's `k_i` is computed across the same pool structure, so the null preserves whatever pool-driven artifacts are present. The signal that exceeds this null is therefore biology beyond the pool structure.

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
- **Phases too cheap to need smoke tests:** 0, 2, 3, 4, 6, 7 main (adult CellTypist), 7b, 7d, 7e, 8a, 8b, 8b follow-ups, 8d, 8e (under ~30 min wall time per tissue, parallelized).

## Dev workflow — split at h5 level, not at runtime
- **dev_split_h5.py** runs ONCE before Phase 0 on dev. Reads the 3 dev h5 files, writes 9 split h5 files (random barcode partition) into `data/dev_split/`, and emits `config/dev_split.yaml`.
- **No pipeline scripts are dev-aware.** All phase scripts run unchanged with `--config config/dev_split.yaml`.
- Pseudo-donors are random cell partitions → numbers MEANINGLESS, smoke test of code paths only.
- **Dev is 4W-only by design.** All 9 pseudo-donors are 4W M Pool1. 8g cannot be exercised meaningfully on dev (one age → every classification = `transient_4W` or `none`). Note: dev being 4W-only means the P1-scANVI branch is NOT exercised on dev — smoke-test it on real P1 data on the WS. Same for 8b follow-ups: they need `within_group_across_age` which requires ≥2 ages.

## Annotation conventions
- **Phase 7 uses per-cluster majority voting** (CellTypist convention), not per-cell argmax. Cells in one Leiden cluster share a label; low-purity (<60% majority OR runner-up >25%) clusters announced for manual review.
- **Phase 7d (subcluster naming) is already cluster-level** — scores aggregate per integer subcluster ID.
- **Phase 7e (cell-type counts diagnostic)** — per-donor × cell-type counts CSV (brain: 3 granularities × `whole`+regions; placenta: `whole` only). Used to sanity-check 8a propeller inputs and for paper Table S?. Read off the annotated h5ad at `results/{tissue}/h5ad/08_annotated/all_samples.h5ad` (legacy "08_annotated" folder naming for what's logically Phase 7).
- **Brain Phase 7 = 4-tier, P1 via scANVI (UPDATED 2026-06-12; supersedes the 3-tier/Di-Bella scheme):**
  - **broad** (`celltypist_broad`, region-FREE, cross-age tier): derived from class per-age — ABC `class_to_broad_csv` for 4W/3mo, `config/rosenberg_class_to_broad.csv` for P1; trailing ` (region)` suffix stripped so all ages align. ~9 classes. THIS is the cross-age tier (8g persistence, 9 cross-species).
  - **class** (`celltypist_class`, canonical 8b/8c key): per-(Leiden cluster × age) MAJORITY vote. Region-TAGGED. ABC 34-vocab for 4W/3mo; Rosenberg ~18-vocab (e.g. CTX Glut, CB GABA) for P1. Two vocabularies coexist by design (per-age native).
  - **subclass** (`celltypist_subclass`, per-cell): ABC 334 for 4W/3mo; raw Rosenberg ~65 fine labels for P1.
  - **region** (`celltypist_region`, per-cell): ABC 12 for 4W/3mo; parsed from Rosenberg label prefix for P1 (CTX/CB/TH/HPF/OLF/STR/MB/non-regional).
  - **P1 = scANVI label transfer from Rosenberg 2018 P2-brain** (`run_scanvi_p1.py`, GPU subprocess called inside `07_annotation.py` for the P1 branch only). Di Bella ABANDONED: cortex-only, mislabeled 42% of whole-brain P1 as erythrocyte (region-coverage failure, NOT ambient — SoupX had stripped Hb). 4W/3mo = CellTypist from ABC. `07c_label_transfer.py` DELETED.
- **Placenta Phase 7 has no CellTypist model.** Curated literature markers + STAMP Spearman correlation against Liu 2024 reference (35 cell types, E9.5-E18.5). Tier 1+2 label-collapse (35 → ~21) + STRICT canonical-marker gates (Neutrophil / Lymphoid / Megakaryocyte) + Xist + Y-gene compartment scoring + EPC/TSC negative-control QC. Canonical key `celltype_majority`. NO broad tier yet — add a compartment grouping (trophoblast/decidua/immune/vascular/erythroid) when Phase 8f cross-tissue needs it.

## Brain marker gate (updated 2026-06-12)
STRICT canonical-marker gates for borderline brain calls. Runs over ALL (cluster×age) rows identically (P1 included). CellTypist/scANVI conf says "how sure among trained classes" — can't verify the cell expresses the biology the label requires. Gate catches false positives. Demoted labels become `unassigned_immune` / `unassigned_glia` / `unassigned_vascular` / `unassigned_erythroid` — which Phase 8 then DROPS (see Phase 8 conventions).

- **`BRAIN_GATE_CONFIG` in `scripts/07_annotation.py`** — five gated types:
  - microglia: ≥2 of {P2ry12, Tmem119, Csf1r, Aif1} → `unassigned_immune`. **Cx3cr1 REMOVED — absent from 10x Flex Mouse Transcriptome v2 panel.**
  - astrocyte: ≥2 of {Aqp4, Gja1, Slc1a3, Aldh1l1} → `unassigned_glia`
  - ol_lineage: ≥1 of {Mbp, Mog, Plp1, Mag, Pdgfra, Cspg4, Olig1, Olig2, Sox10} → `unassigned_glia`. **OPC markers ADDED.**
  - endothelial: ≥2 of {Cldn5, Pecam1, Cdh5} → `unassigned_vascular`
  - **erythroid: ≥2 of {Hbb-bs, Hbb-bt, Hba-a1, Hba-a2, Alas2} → `unassigned_erythroid`.**
- **`MARKER_PRESENCE_THRESHOLD = 0.20`** — marker "present" if ≥20% of cells in (cluster, age) have lognorm > 0.
- **Audit CSV columns:** `markers_checked`, `markers_present`, `gate_outcome` (`no_gate`|`passed`|`demoted`), `gate_label`.
- **Production run 2026-06-12:** 4 demotions (Astro-Epen, weak astro markers), 0 erythroid demotions.

## Brain age-composition sanity (added 2026-06-10)
Diagnostic CSV `tables/07_annotation/07_annotation_age_composition_sanity.csv` flags developmentally-implausible (cluster × age) rows. **Informational only — does not modify labels.**

- `BRAIN_AGE_EXPECTATIONS` list in `scripts/07_annotation.py` — keyword + expected_ages + flag_name tuples.
- Current rules: radial glia / intermediate progenitor / IPC / neuroblast / glioblast / erythrocyte / erythroid progenitor → expected only at P1.
- Production run 2026-06-12: 0 flags (good).

## 8b DE visualization blocklist (locked 2026-06-15)
The 8b master CSV reports DE on ALL genes. For visualization (heatmap top rows, bubble top rows, dotplot gene picks, volcano labels, RRHO labels), some genes are excluded by default — they overwhelm headline figures with technical / ambient signal, not stress biology.

- **`BLOCKLIST_FOR_VIZ` in `scripts/08b_de_summary.py`** — 17 genes:
  - Erythroid (residual ambient even after SoupX): Hbb-bs, Hbb-bt, Hba-a1, Hba-a2, Hbb-b1, Hbb-b2, Hbb-y, Hbb-bh1, Hbb-bh2, Alas2.
  - Sex-linked (escape sex-stratified design): Xist, Tsix, Ddx3y, Uty, Eif2s3y, Kdm5d, Eif2s3x.
- **`BLOCKLIST_PREFIXES = ("mt-",)`** — all mitochondrial transcripts (`mt-Co1`, `mt-Nd1`, etc.).
- **Master CSV is UNFILTERED.** The blocklist applies at TOP-N selection time only.
- **`--no-blocklist` CLI flag** to disable for QA / sensitivity analysis.
- Why these? In P1 brain, residual Hb in the master CSV showed up as ~650 sig age-DE genes in the `within_group_across_age` contrast — direction consistent across ALL groups including Relaxed. That's developmental erythroid signal (nucleated erythroblasts in residual vasculature, see Phase 1 = SoupX), NOT stress biology. Excluding them from visualization keeps the headline figures focused.

## Plot format strategy (locked 2026-06-05)
- **Default: PNG @ 300 DPI** for all pipeline plots. Cell-level UMAPs with 600K dots aren't suitable for pure SVG.
- **Paper figures only** (~5-10 plots): refactor to PNG + PDF hybrid via `_utils.savefig(fig, path, dpi=300)`.
- Don't refactor all plotting scripts. Targeted post-Phase 8 update only.
- **on-data UMAP labels require categorical dtype** — cast the color column to `category` before `sc.pl.umap(..., legend_loc="on data")`, else scanpy silently draws nothing.
- **constrained_layout fights `bbox_inches="tight"`** — `_utils.safe_fig` (and the local copy in `08b_followup_plots.py`/`08b_disruption_shuffle_test.py`) auto-detects and skips the bbox tightening when constrained_layout is active. Suptitle on a constrained_layout figure should NOT pass `y=` (let constrained_layout reserve title space automatically); footnotes go to `y=0.005` inside the figure not `y=-0.02`.

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
- **8b follow-ups share the 08b_de phase folder** — both tables and plots. Plots live under `plots/08b_de/summary/{disruption,consistency,shuffle_test}/{sex}/{level}.png`. Tables alongside the main 8b results CSV in `tables/08b_de/`.

## Offline-audit CSVs (no workstation access post-run)
- **02_soupx** `02_soupx_summary.csv` — per-sample rho_mean, rho_min, rho_max, pct_removed, n_cells, n_clusters, elapsed_sec, mode (autoEst | manual).
- **07** `07_annotation_class_per_cluster_age.csv` — augmented with `markers_checked`, `markers_present`, `gate_outcome`, `gate_label`.
- **07** `07_annotation_age_composition_sanity.csv` — developmentally-implausible (cluster × age) rows.
- **07e** `07e_celltype_counts.csv` — per-donor × cell-type counts (brain: 3 granularities × whole+regions; placenta: whole only). Sanity-check for 8a inputs; paper Table S?.
- **8a** `08a_dropped_cells_per_donor.csv` — per-donor contaminant + unassigned counts/fractions (the mass dropped from the tested composition).
- **8b** `08b_de_gene_expression_per_sample.csv` — per-sample mean lognorm of DE genes, keyed celltype × gene × sample_id.
- **8b follow-ups** `08b_developmental_disruption_summary.csv` — per (sex × level × celltype): 5 direction-class gene counts (universal / relaxed_only=LOST / stress_shared=GAINED / early_only / late_only) + mean `|LFC|` per group for the LOST class (the effect-size collapse columns).
- **8b follow-ups** `08b_developmental_disruption_genes.csv` — long-form gene-level direction class assignments per (sex × level × celltype × gene).
- **8b follow-ups** `08b_disruption_shuffle_test.csv` — k-preserving null per (sex × level × celltype): 6 disjoint sig-pattern category counts (R-only / E-only / L-only / R∩E / R∩L / E∩L), `n_k0..n_k3`, 6 categories × {enrichment p, depletion p} × {raw, BH}, within-stratum chi-square (`chi2_k1`, `chi2_k2` + p-values), permutation-null sanity check on the LOST/GAINED diff. Headline columns for the paper: `obs_lost`, `obs_gained`, `obs_r_only`, `obs_re_only`, `obs_el_only`, `p_lost_BH`, `p_gained_dep_BH`, `p_re_only_enr_BH`, `p_el_only_dep_BH`.
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
- **WS results saved on Mac under `results_WS/`** (not `results/`), in matching `{tissue}/{plots,tables}/<phase>` subfolders.
- **tmux for any multi-minute job.** Detach `Ctrl-b d`, reattach `tmux attach -t <name>`, list `tmux ls`.
- **Pool → contents map:**
  - Pool1: 16 brain (3mo all groups + 4W males)
  - Pool2: 16 brain (4W females + P1 Early/Relaxed)
  - Pool3: 2 brain P1 Late + 14 placenta E12.5 (incl. duplicate CES2.3 to drop)
  - Pool4: 10 placenta (2 E12.5 Relaxed + all E18.5)
- **GPU operational hygiene:** pre-flight `nvidia-smi --query-gpu=memory.used` check before launching scVI/scANVI; refuse if non-display VRAM > 2 GB. `PYTORCH_CUDA_ALLOC_CONF=max_split_size_mb:512` for long scVI runs. Explicit cleanup between GPU phases: `del model; torch.cuda.empty_cache(); gc.collect()`.
