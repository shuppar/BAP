# Instructions for working on this project with Claude

## Response style
- **Be brief.** No long preambles, no excessive caveats, no over-explaining.
- **Don't restate what I just said.** Move to the substance.
- **Step by step.** Build one thing, verify it works, then move on. Don't write 5 files at once.
- **No prolix bullet lists.** A short paragraph or 3-4 line list is usually enough.
- **Be honest when something won't work** or when you're uncertain. Don't manufacture confidence.

## Code style
- **Simple > clever.** Plain Python scripts in `scripts/`, not a Python package. No Pydantic schemas, no abstract base classes, no dependency injection. Just functions and main().
- **One file per phase** (e.g. `01_validate.py`, `02_qc.py`). Each is runnable standalone.
- **Shared helpers in `scripts/_utils.py`** (leading underscore = not a phase entry point). Currently provides `load_config`, `add_lognorm`, `phase_paths`, `select_accelerator`. Add to it when something gets duplicated 2+ times.
- **Configs are plain YAML dicts.** No inheritance trees, no schema validation.
- **R is called as subprocess**, not via rpy2.
- **Idempotent steps where reasonable**, but don't over-engineer.
- **Raw counts in `.X`, lognorm computed on demand.** `04_integration_prep.py` computes the lognorm layer for Phase 5's pre-integration UMAP, then Phase 5 drops it before saving. Notebooks and downstream phases call `_utils.add_lognorm(adata)` after loading the integrated h5ad.

## Pipeline architecture decisions (don't re-litigate)
- **Language:** Python primary (Scanpy/scvi-tools). R subprocess for scDblFinder + CellChat.
- **Env:** uv + Python 3.12 on Mac (dev), workstation has GPU + R installed.
- **Subprocess R:** scripts/run-X.R called from Python, exchange via TSV/JSON.
- **Conda is blocked** at corporate firewall — don't suggest it.
- **CellBender + scVI**: GPU phases, run on workstation, not laptop.

## Dataset specifics
- **3 groups**, not 2: Early_Stress, Late_Stress, Relaxed (Relaxed = reference)
- **34 brain + 23 placenta** samples (after dropping duplicate CES2.3)
- **Ages:** brain P1/4W/3mo, placenta E12.5/E18.5
- **Pools (libraries):** Pool1-4, used as scVI batch_key
- **Known confounds** (see project doc §2): P1 Late Stress only in Pool3, placenta cross-age not comparable
- **No dam ID recorded** — treat each pup as independent observation, flag the caveat
- **Sex=TBD** for all E12.5 placenta — infer from Y-chromosome

## Compute constraints
- **Laptop:** 12 GB RAM, Apple Silicon, no GPU. Dev only with subsetting (3 samples × 500 cells)
- **Workstation:** 258 GB RAM, 56 CPU cores, RTX 4500 Ada (24 GB VRAM). Production runs.
- **Network:** conda channels blocked; PyPI/CRAN/Bioconductor work.

## Workflow conventions
- **`run_pipeline.sh` is a manual, not a history.** Tight start-to-finish
  walkthrough: what to run, in what order, what to inspect. NOT a changelog or
  bug/refactor diary (git history covers that). Self-contained — readable without
  opening the scripts — but terse: fragment comments over sentences, one block
  per phase (deps → command → key outputs → commit), one line of rationale only
  where a choice is non-obvious. Complete but not prolix.
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

## External identifiers: verify or flag (added after a bug-hunt session)
Several bugs came from writing external identifiers from memory and presenting
them as verified. Root cause: plausible-looking names that don't exist. Rules:

- **Any external identifier must be verified against docs or flagged.** This
  covers: PyPI/conda package names, model names (e.g. CellTypist models),
  gene symbols, API function names and keyword arguments, file/dataset names.
- If verified: fine, use it.
- If NOT verified: add an inline `# UNVERIFIED — check before prod` comment, OR
  wrap the use in a guard that fails loudly with a helpful message. Never write
  an unverified identifier as if it were confirmed.
- Things that turned out wrong in practice (keep as cautionary examples):
  - `Mouse_Brain_Atlas` — does NOT exist; only `Developing_Mouse_Brain.pkl` ships with CellTypist.
  - `abc-atlas-access` (PyPI) — wrong; real package is `abc_atlas_access`, GitHub-only.
  - `score_genes(use_raw=False)` — runs on raw `.X`, not lognorm; need `layer="lognorm"`.
  - `Mlf1ip`/`Fam64a`/`Hn1` — outdated MGI symbols (now `Cenpu`/`Pimreg`/`Jpt1`).

## No silent failures (added same session)
Wrong-but-plausible output is worse than a crash, because it looks correct.

- **A correctness-critical step that can't run correctly must raise, not warn-and-continue.**
  Examples that now hard-fail (not warn):
  - Too few cell cycle genes match var_names → `raise` (likely Ensembl-vs-symbol mismatch).
  - No marker genes match var_names → `raise` (don't fall back to Leiden numbers).
  - `condition_cell_cycle: true` but `cc_difference` missing → `sys.exit`.
  - Required covariate column missing → `sys.exit` (don't silently drop it).
- **A warn-and-skip is only acceptable when skipping an OPTIONAL output** and the
  skip is clearly announced. Examples that legitimately warn-and-skip:
  - CellTypist not installed → skip reference track, use marker-based labels.
  - No curated markers present for a dotplot → skip that one plot.
  - An age has no CellTypist model → those cells get provisional labels (announced).
- **Never leave NaN labels that surface as a "nan" category** in plots/tables.
  Fill with an explicit sentinel (e.g. `"no_model"`) so the gap is visible.
- When in doubt: fail loud and early (the Phase 0 validation-gate philosophy),
  not deep into a multi-hour run or — worse — silently in the output.

## Plots must carry biological meaning (not abstract designs)
Every figure should let a reader name the biology — the genes, cell types, or
pathways affected — without cross-referencing a separate table. A shape with no
labels is a decoration, not a result.

- **Label the things that matter.** Volcano → name the top significant genes on
  the plot. Composition → name the cell types that shift. Pathway/GSEA → name the
  pathways, not just bars. Heatmaps → real gene/pathway/cell-type names on the axes.
- **Cap labels for readability** (e.g. top ~25 by significance), and say how many
  more exist if truncated. Readable beats exhaustive.
- **Gene identifiers must be human-readable.** If var_names are Ensembl IDs, map
  to symbols (var['symbol']) before labeling — never ship a plot of ENSMUSG IDs.
- **State the contrast and thresholds on the plot** (what-vs-what, padj/LFC cutoffs)
  so it's self-describing.
- Rule of thumb: if the figure can't tell you which genes/cell types/pathways are
  involved, it isn't done yet.

## Isolate fragile dependency stacks (don't pin them into the main env)
When a tool drags an incompatible dependency stack, isolate it — don't pin the
main env backward to accommodate it.
- **scCODA was abandoned** for composition (8a): its TF/TFP/arviz/matplotlib/
  numpy/setuptools pins cascaded endlessly and fought the scanpy/scVI stack.
  Replaced with **propeller (speckle+limma) via R subprocess** — clean Bioconductor
  install, and limma's empirical-Bayes moderation is better for small n anyway.
- **CellBender** stays in its own venv (.venv-cellbender) for the same reason.
- Lesson: if >2 dependency dominoes fall chasing one tool, stop and find the
  lighter alternative (often an R/Bioconductor tool via subprocess) rather than
  pinning the shared env into a corner.

## Phase 8 conventions
- **Statistical unit is the animal (donor_id)**, never the cell. Composition =
  per-donor cell-type counts; DE = pseudobulk (sum raw counts per donor). No dam
  ID → each pup independent (anti-conservative; caveat carried in `flag`).
- **Corrected p-values everywhere** — propeller FDR (8a), DESeq2 padj (8b),
  per-collection BH FDR (8c). No raw p-value drives any significance call.
- **GSEA gene sets = mouse MSigDB via msigdbr** (MH+M2+M5+M8), exported once to
  refs/msigdb_mouse.tsv. NOT decoupler's get_resource(MSigDB,mouse) — it's broken.
  FDR corrected WITHIN each collection (sizes differ ~150x); pooled kept as ref.
- **Multi-database plots = side-by-side panels**, one subplot per collection,
  figure width scaled to #collections so nothing shrinks/clips. Use
  constrained_layout (NOT with bbox_inches='tight' — they fight and clip).
- **--subcluster flag** on 8b/8c/8e runs on the 7b subcluster object (label
  col `subcluster_name` from 7d, else `subcluster` integer from 7b), writing
  `*_subcluster_{slug}` outputs to separate folders so they don't collide
  with main runs.
- **Subcluster runs are a LOOP, not one cell type.** Production runs every
  focal cell type through 8b, 8c, and 8e via the same `CELL_TYPES` array
  defined once in run_pipeline*.sh. Don't add a new cell type in one place
  and forget the others — every phase that takes `--subcluster` must loop.
- **Subcluster naming (7d)**: 7b produces INTEGER subcluster ids; 7d names them
  via CellTypist majority_voting (where a model exists) + literature-marker
  scoring from config/subcluster_markers.yaml, aggregated per cluster (mean
  score → argmax). Already cluster-level by construction.
- **TF activity (8c) = REQUIRED in production**: ULM on DE Wald stats vs
  CollecTRI mouse network, per contrast×celltype, BH-FDR within
  celltype×contrast. Always pass `--tf` (or YAML `pathways.run_tf_activity:
  true`). It gates 8f view 5 (TF concordance) and 8g view 3 (TF persistence),
  and cannot be recovered without re-running 8c. Needs network (omnipath).
- **8e cell-cell communication = LIANA+ in main env, no sidecar.** Three arms
  in one script: (1) baseline `rank_aggregate` per group×age, (2)
  differential via `li.multi.df_to_lr` reading 8b's Wald stats, (3) per-donor
  for statistics. Covers all three group comparisons (ES-v-Rel, LS-v-Rel,
  ES-v-LS). Clustered Δ heatmaps with pathway+celltype colour bars; rank-rank
  concordance scatters answer "do ES and LS hit the same programs?". The
  `--focus-celltypes` filter = union of YAML stress_focused_cell_types and
  top-10 cell types by max |Δ|.
- **8f cross-tissue = six views, all reproducible from 8b/8c CSVs.** Two
  biologically aligned arms (E12.5 placenta Early → P1/4W/3mo brain Early;
  E18.5 placenta Late → P1/4W/3mo brain Late; P1 Late carries
  confounded_with_pool flag). Views: DEG overlap (hypergeom), RRHO (custom
  NumPy), pathway concordance from 8c GSEA, LR cross-tissue mechanistic
  hypotheses (placental ligand × brain receptor from liana mouseconsensus
  with `stress_axis` flag column), TF concordance from 8c TF activity, ORA
  of overlap genes vs MSigDB. The LR table is the publication-quality output.
- **NO cross-tissue cell-cell communication.** Within-tissue interaction
  scoring (CellChat/CellPhoneDB framework) doesn't extend across the BBB —
  literal "placenta-cell-to-brain-cell signalling" is implausible and
  reviewers will flag it. The published placenta-brain axis (Goeden 2016,
  Vacher 2021, Wu 2017, Bonnin 2011) is endocrine/paracrine: placental
  ligand → circulation → brain receptor. 8f view 4 (LR from DE) is the
  correctly framed version. Don't add 8e-style cross-tissue CCC.
- **8g cross-age persistence = brain only.** Placenta has incomplete
  cross-age factorial (E12.5 = Early+Relaxed, E18.5 = Late+Relaxed; no
  factorial across ages), so 8g exits cleanly with `tissue: placenta`.
  Persistence classes: persistent, resolving_early, established_late,
  P1_only, transient_4W, emergent_3mo, P1_3mo_only, persistent_directionswap.
  Same-direction sign required for "persistent". Cross-arm core signature
  (view 6) = intersection of persistent calls in BOTH arms with consistent
  direction = paper-quality table.

## Dev workflow — split at h5 level, not at runtime
- **dev_split_h5.py** runs ONCE before Phase 0 on dev. Reads the 3 dev h5
  files, writes 9 split h5 files (random barcode partition) into
  `data/dev_split/`, and emits `config/dev_split.yaml` (9 pseudo-samples,
  donor_id suffixed _ps1/_ps2/_ps3, group/age/sex/pool inherited).
- **No pipeline scripts are dev-aware.** All phase scripts run unchanged with
  `--config config/dev_split.yaml`. Workstation uses `--config config/brain.yaml`
  (or placenta) with no splitter step.
- Pseudo-donors are random cell partitions of one animal → numbers MEANINGLESS,
  smoke test of code paths only. The split round-trips cleanly through
  `sc.read_10x_h5` (verified by writer self-check).
- The old `dev_pseudoreplicate.py` workflow (late mutation of donor_id before
  8a) is OBSOLETE — don't use it.
- **Dev is 4W-only by design.** All 9 pseudo-donors are 4W M Pool1 (one M
  sample per group × 3 splits). Consequence: 8g cannot be meaningfully
  exercised on dev — every classification ends up `transient_4W` or `none`
  because there are no other ages to compare against. That's correct
  behaviour, not a bug. 8g is workstation-only in practice; the dev run
  exercises code paths, not biology.

## Annotation conventions
- **Phase 7 uses per-cluster majority voting** (CellTypist convention), not
  per-cell argmax. Cells in one Leiden cluster share a label; low-purity
  (<60% majority) clusters announced in stdout for manual review.
- **Phase 7d (subcluster naming) is already cluster-level** by construction —
  scores aggregate per integer subcluster ID, then assign one name to all
  cells in the cluster.

## Trajectory (8d) — no velocity, all ages equal
- **NO RNA velocity, NO CellRank.** 10x Flex is probe-based (exon-only) and
  cannot resolve spliced/unspliced — velocity is not feasible (10x kb 25938615598477).
  Without velocity CellRank only duplicates PAGA, so it's dropped. PAGA +
  diffusion pseudotime (DPT) are the trajectory methods. Don't re-add velocity/CellRank.
- **All ages treated identically in DPT** — do NOT gate any analysis off for some
  ages and not others. DPT group comparison runs pooled across ages AND per age;
  age-split rows carry a pool_age_confound caveat in 'note' (not a validity gate).
- **PAGA edges are hypotheses, not transitions.** trajectory_paga_edge_diagnostics.csv
  flags ambient-driven / doublet-driven / shared-gene edges for offline audit.

## Output organisation
- **Tables in per-phase subfolders**: `tables/<phase_dir>/<phase>_<name>.csv`
  (e.g. `tables/08b_de/08b_de_results.csv`). Helper: `_utils.phase_table_dir(cfg, label)`.
- **Plots also in per-phase subfolders**: `plots/<phase_dir>/...`.
- **Filenames carry the phase prefix** so a file is identifiable by name alone
  when copied out of the directory.
- **Subcluster runs go to suffixed folders**, not the main phase folder, so
  re-running with `--subcluster <slug>` doesn't overwrite the main run's
  plots and tables. Convention used by 8b, 8c, 8e:
  `plots/08e_communication_subcluster_excitatory_neurons/...` and
  `tables/08e_communication_subcluster_excitatory_neurons/08e_*.csv`.

## Offline-audit CSVs (no workstation access post-run)
Beyond the standard result tables, persist enough to troubleshoot offline:
- **8b** `08b_de_gene_expression_per_sample.csv` — per-sample mean lognorm of DE
  genes, keyed celltype × gene × sample_id (+ n_cells backing each). Join target.
- **8c** `08c_pathway_leading_edge.csv` — genes driving each significant pathway
  with log2FC + direction, keyed contrast × level × celltype × pathway × gene.
  Join to the 8b matrix on (celltype, gene) for per-sample levels.
- **8c** `08c_tf_activity.csv` — TF activity scores + FDR per contrast×celltype×TF.
- **8d** `08d_trajectory_paga_edge_diagnostics.csv` — per cell-type-pair edge audit.
- **8e** `08e_lr_quantified.csv` — per-donor LR scores backing the Mann-Whitney
  group comparisons in the per-donor arm.
- **8f** `08f_lr_cross_tissue.csv` — placental ligand × brain receptor mechanism
  hypotheses with `stress_axis` flag column. The publication-quality output.
- **8g** `08g_core_signature_genes.csv` — features persistent in BOTH stress
  arms with consistent direction. The paper-quality "core stress signature".
- When adding analyses, ask: what's recoverable offline vs. lost? Persist the
  workstation-locked diagnostics as CSVs.

## Ask before strong scientific calls
- Don't drop an analysis, exclude samples/ages, or add complexity (extra tools,
  sidecar venvs) without checking it earns its place for THIS dataset (modality,
  n, confounds). Surface the question; don't bake the decision in silently.
- Specifically don't propose: cross-tissue cell-cell communication (BBB makes
  literal placenta-cell-to-brain-cell signalling implausible — 8f view 4 is
  the correctly framed endocrine version); RNA velocity / CellRank (10x Flex
  is probe-based — exon-only, no spliced/unspliced); scCODA (dependency stack
  fights scanpy; propeller via R is the replacement).

## Always state where files go
- When presenting files in a Claude response, ALWAYS say which directory each
  one goes into (scripts/, config/, repo root, etc.). Don't leave the user to
  guess. Especially important since the layout has multiple destinations
  (scripts/, config/, refs/, data/dev_split/, repo root).

## Workstation infrastructure (added at run start)
- **SSH:** `ssh poller@172.17.213.147` (from Mac). User: `poller`.
- **Local repo on Mac:** `/Users/shuppar/Downloads/BAP_data_1/Analysis/`
- **Workstation project root:** `/home/poller/BAP-BrainPlacenta/` (NVMe — fast).
  Renamed from `Analysis/` so the folder name reflects the study.
- **Raw tars on workstation:** `/media/poller/PollerLab-1/BAP-data1/processed_260411_Shiv_FLEX_260411_Shiv_FLEX_Pool{1,2,3,4}.tar` (USB-HDD).
- **Extracted Cell Ranger output:** `/media/poller/PollerLab-1/BAP-data1/Analysis/data/Pool{1,2,3,4}/per_sample_outs/<sample_id>/` (USB-HDD).
  Reached from the project via the symlink `BAP-BrainPlacenta/data` → that path.
  Matches `h5_path` / `raw_h5_path` in `sample_metadata.csv` when run from the
  project root.
- **Disk layout rationale:** project on NVMe (every intermediate h5ad write is
  fast); raw data stays on USB-HDD (read once per phase). Symlink bridges them.
- **Pull-from-Mac convention:** code edits happen locally and rsync to the WS.
  Standard exclude list: `results/`, `data/`, `.venv/`, `.venv-cellbender/`,
  `__pycache__/`, `.git/`, `*.h5ad`, `logs/`, `.DS_Store`.
  Standard rsync flags: `-av --progress --chmod=Fu+x` (the `--chmod=Fu+x`
  preserves the executable bit on `.sh` files — without it shell scripts land
  non-executable on the workstation).
- **tmux for any multi-minute job.** Detach: `Ctrl-b d` (Control key, not Cmd).
  Reattach: `tmux attach -t <name>`. List: `tmux ls`.
- **Pool → contents map** (confirmed against `sample_metadata.csv`):
  - Pool1: 16 brain (3mo all groups + 4W males)
  - Pool2: 16 brain (4W females + P1 Early/Relaxed)
  - Pool3: 2 brain P1 Late + 14 placenta E12.5 (incl. duplicate CES2.3 to drop)
  - Pool4: 10 placenta (2 E12.5 Relaxed + all E18.5)
