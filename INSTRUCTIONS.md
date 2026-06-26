# Instructions for working on this project with Claude

Get broad context from the `snRNAseq_project_summary.md` file.

## Response style
- **Be brief.** No long preambles, no excessive caveats, no over-explaining.
- **No need to print your thoughts** unless it is helpful to either of us. Keep chats from getting long by not printing reasoning unnecessarily.
- **Don't restate what I just said.** Move to the substance.
- **Step by step.** Build one thing, verify it works, then move on. Don't write 5 files at once.
- **Be honest when something won't work** or when you're uncertain. Don't manufacture confidence.
- **Always give commands to run a specific script** (mention where: Local Mac or remote WS), or rsync a specific file/folder.

## Code style
- **Parallel compute** wherever we can. See "Parallelism is mandatory" — the standard is `_utils.parallel_map`.
- **Simple > clever.** Plain Python scripts in `scripts/`, not a package. No Pydantic schemas, no ABCs, no dependency injection. Just functions and `main()`.
- **One file per phase** (e.g. `01_validate.py`). Each runnable standalone.
- **Shared helpers in `scripts/_utils.py`**: `load_config`, `load_contrasts`, `phase_paths`, `phase_table_dir`, `add_lognorm`, `select_accelerator`, `iter_strata`, `parallel_map`, `unassigned_mask`. Add to it when something is duplicated 2+ times.
- **Configs are plain YAML dicts.** No inheritance, no schema validation.
- **R is called as subprocess**, not via rpy2.
- **Idempotent steps where reasonable**, but don't over-engineer.
- **Raw counts in `.X`, lognorm computed on demand** via `_utils.add_lognorm(adata)`.
- **Time estimate.** Always give one based on machine (WS or Mac) config.
- **Tmux command** format: `tmux new -d -s soupx_placenta 'uv run python -u scripts/02_soupx.py --config config/placenta.yaml --n-jobs 6 2>&1 | tee logs/02_soupx_placenta_full.log'` then `tail -f /home/poller/BAP-BrainPlacenta/logs/02_soupx_placenta_full.log`
- **rsync format**: `rsync -av --chmod=Fu+x /Users/shuppar/Downloads/BAP_data_1/Analysis/scripts/<file> poller@172.17.213.147:/home/poller/BAP-BrainPlacenta/scripts/`. WS results saved on Mac under `/home/poller/BAP-BrainPlacenta/results_WS` (NOT `results`).

## Naming conventions (locked 2026-06-25)
- **Mouse pipeline:** numbered phase scripts `0X_*.py` / `08x_*.py`; R workers `run_*.R`.
- **Human cross-species (Phase 9):** scripts `h09X_...py` (`h09a`…`h09k`, plus `h09_summary_plots.py`, `h09k_diagnostics.py`, `h09k_plots.py`, `h09k_rrho_maps.py`); R workers `h_run_*.R` / `h_fetch_*.R` (e.g. `h_run_soupx_from_raw.R`, `h_run_singler.R`, `h_fetch_genesets.R`). Keeps the human arm visually separate from the numbered mouse pipeline.

## Parallelism is mandatory for repeated work (not optional)
Any phase script that loops over samples/jobs/contrasts/cell-types and, per item, launches a subprocess (R worker) or calls a heavy function MUST parallelize via `_utils.parallel_map` and expose `--n-jobs`. A bare for-loop starting one subprocess per item is a performance bug.
- `for item, result, err in parallel_map(fn, items, n_jobs=args.n_jobs, desc="..."):` yields (item, result, error).
- Default `use_threads=True` (subprocess/IO-bound). Use `use_threads=False` (processes) only for CPU-bound pure-Python work; then `fn` and items must be picklable.
- Default `--n-jobs` 8; on the workstation pass 16–24 for light R workers.
- Reference implementation: `08a_composition.py`.
- **CPU-bound exceptions (use_threads=False):** the shuffle null in `08b_disruption_shuffle_test.py`; the permutation nulls in `h09f_overlap_null.py` and `h09k_admati_2x2.py` (each chunk does a batch of RRHO shuffles; arrays pickled once per chunk, NOT per shuffle).
- **Plotting note:** 8e/8f/8g and all `h09*` summary plotting is intentionally SERIAL (matplotlib not thread-safe).
- **8f/8g parallelization caveat (2026-06-25):** single-threaded across slices; process-based parallel_map would pickle the 20M-row brain DE frame to every job. Pre-filter each slice serially, then parallel_map the small frames if ever needed.

## Pipeline architecture decisions (don't re-litigate)
- **Language:** Python primary (Scanpy/scvi-tools). R subprocess for scDblFinder, propeller, SoupX, fgsea, SingleR.
- **Env:** uv + Python 3.12. **Conda is blocked** at corporate firewall — don't suggest it.
- **scVI**: GPU phase, runs on workstation.
- **Phase 1 ambient correction = SoupX (locked 2026-06-10).** CellBender abandoned (pickle bug).

## Dataset specifics (mouse)
- **3 groups**, not 2: Early_Stress, Late_Stress, Relaxed (Relaxed = reference)
- **34 brain + 23 placenta** samples (after dropping duplicate CES2.3)
- **Ages:** brain P1/4W/3mo, placenta E12.5/E18.5
- **Pools (libraries):** Pool1-4, used as scVI batch_key
- **Known confounds:** P1 Late Stress only in Pool3; placenta cross-age not comparable
- **No dam ID recorded** — treat each pup as independent; flag the caveat
- **Sex=TBD for all E12.5 placenta** — inferred from Y-chromosome via `01_validate.py`
- **`assigned_sex` is the source of truth for sex covariates.** Source: `sample_metadata.csv` `sex` column.

## Compute constraints
- **Laptop:** 12 GB RAM, Apple Silicon, no GPU. Dev only with subsetting. **All Phase-9 plotting runs fine here offline** (small CSV/parquet inputs; the heavy compute is already done on WS).
- **Workstation:** 258 GB RAM, 56 CPU cores, RTX 4500 Ada (24 GB VRAM). Production runs.
- **Network:** conda channels blocked; PyPI/CRAN/Bioconductor reachable. NVIDIA PyPI reachable for cuML.

## Workflow conventions
- **`run_pipeline.sh` is a manual, not a history.**
- **Source of truth for samples:** `sample_metadata.csv`.
- **Dev runs:** `config/dev_split.yaml`. **Full runs:** `config/brain.yaml` / `config/placenta.yaml`.
- **Outputs:** `results/{tissue}/<phase>/` for prod, `results/dev/<phase>/` for laptop tests. **Phase-9 human outputs live under `data/human_validation/placenta/<dataset>/{h5ad,tables,plots}/`.**

## When asking me questions
- **One question at a time** usually; **3 max**.
- **Single-select buttons** over multi-select where possible.
- **Don't ask for confirmation on small obvious things** — just do them.

## What I don't want
- Don't propose architectural pivots (containers, Nix, etc.) — uv + scripts is settled.
- Don't add abstraction layers "for future flexibility." YAGNI.
- Don't write 200-line responses with 5 nested headers when 30 lines suffices.
- Don't apologize repeatedly when correcting something.

## External identifiers: verify or flag
Several bugs came from writing external identifiers from memory as if verified. **Any external identifier (PyPI/conda package, model name, gene symbol, API function/kwarg, file/dataset/COLUMN name) must be verified against docs/data or flagged.**
- Cautionary examples: `Mouse_Brain_Atlas` (doesn't exist); `abc-atlas-access` (real is `abc_atlas_access`); `score_genes(use_raw=False)` runs on raw `.X`; `Mlf1ip`/`Fam64a`/`Hn1` outdated (`Cenpu`/`Pimreg`/`Jpt1`); `multi_class='ovr'` removed in sklearn 1.7; `Cx3cr1` NOT in 10x Flex Mouse v2 panel (use P2ry12/Tmem119/Csf1r/Aif1); `datasplitter_kwargs` is the scvi-tools 1.4.3 kwarg.
- **Contrast names are NOT uniform across tissues:** brain `*_per_age` (age in `group_level`), placenta `*_E12.5`/`*_E18.5`. Match by family prefix (`contrast_family`), NEVER exact string.
- **`08b_developmental_disruption_genes.csv` stores its class in a column literally named `direction`** (NOT `direction_class`).
- **NEW (Phase 9):** verify external dataset structure on the WS before writing parsers — `unassigned_mask(obs, label_cols)` takes TWO args; Admati files are TRANSPOSED (metadata rows then gene rows, cells in columns); Admati gene IDs are human SYMBOLS (no Ensembl map); the `sn_..._allcells` file is trophoblast-ONLY despite the name; the figshare file ID for the sc all-cells matrix is **41003240** (5.9 GB). Always `cut -f1 | head` to confirm metadata-row vocabulary before parsing.

## No silent failures
Wrong-but-plausible output is worse than a crash. **A correctness-critical step that can't run correctly must raise, not warn-and-continue.** Warn-and-skip is only acceptable for an OPTIONAL output and the skip must be announced. Never leave NaN labels surfacing as a "nan" category.
- **Phase-9 guards:** `h09j` asserts each gene chunk's cell-column count matches the obs row count (catches ragged parsing); `h09k`/`h09e` skip a compartment loudly (`-- skipped`) when a group has too few donors; the human-DE helpers require both conditions ≥ MIN_DONORS.

## The internal-comma trap (CSV parsing)
Several Phase-8 CSVs have a `pair` column holding a Python list (internal comma inside quotes). **Never `cut`/`awk` these — always pandas.** For big mixed-dtype tables use `low_memory=False`; for the 6.7 GB leading-edge file read chunked (`chunksize=2_000_000`) with `usecols`.

## Plots must carry biological meaning (not abstract designs)
Every figure should let a reader name the biology without a side table.
- Label what matters (volcano → top sig genes; GSEA → pathways; heatmaps → real names). Cap labels (~top 25).
- Gene identifiers human-readable (`var['symbol']`). State contrast + thresholds on the plot.
- **Don't apply an effect-size floor to distribution plots** (volcano, rank-rank scatter, 8f LR scatter, h09 leading-edge scatters). Floors apply only to ranked/aggregating plots.
- **adjustText for label de-overlap** (`uv add adjustText`); `h09_summary_plots._label_points()` is the shared helper (leader lines + plain-annotate fallback). Reuse it in `h09k_plots.py`.
- **Phase-9 plotting is CSV-only and Mac-runnable.** Each plot script reads the saved tables (and `h09k_rankings.parquet` for RRHO maps) — never recompute DE. When a plot needs FDR-sized dots, the source table must carry FDR (h09k now saves `FDR_mouse`/`FDR_human`); if a CSV predates that, re-run the producing script once.

## Ask before strong scientific calls
- Don't drop an analysis, exclude samples/ages, or add complexity without checking it earns its place. Surface the question; don't bake it in silently.
- Don't propose: cross-tissue CCC; RNA velocity / CellRank; scCODA.
- When a result is sparse or surprising, interrogate before accepting OR fixing: cutoff artifact, silent-filter bug, or real biology? **Phase-9 example: the Admati 2×2 anomalies (peak=37; `discordant` trophoblast cell; loPE≫eoPE) were interrogated via `h09k_diagnostics.py` BEFORE interpretation — sign-flip = tail artifact (positive global Spearman), loPE>eoPE = mostly real (survives control-downsampling) not pure power. Always run the diagnostic before building the figure.**

## Phase 8 conventions
(unchanged — cross-cutting + 8f/8g-specific + stage-specific rules; statistical unit = animal; declarative contrasts/strata; drop contaminants+unassigned; `~ sex + pool + group`; `group_level` holds age for `*_per_age`; always filter sex+level before classifying; corrected p-values everywhere; GSEA = mouse MSigDB via msigdbr MH+M2+M5; TF activity REQUIRED in 8c production.)

## Phase 8f / 8g conventions (locked 2026-06-25)
(unchanged — six 8f views from 8b/8c CSVs; RRHO vectorized; `--logfc-cutoff 0.5`; plot-only `--plot-quantile 0.75`. 8g brain only; persistence classes per (celltype, level, feature, arm); B trajectory-shape / C persistence×disruption / View-7 8f-bridge modules. Interpretation locked: 0 persistent genes; IFN perinatal-transient; ECM durable; gliogenesis core.)

## Disruption analysis framing (8b, locked 2026-06-15)
(unchanged — "When age-DE signal is shared across two groups, Relaxed is almost always one of them.")

## Phase 9 — cross-species validation (PLACENTA ARM COMPLETE 2026-06-25; BRAIN ARM A NEXT)

**Two scientific arms reported separately, NOT pooled** (locked 2026-06-05): ARM A psychiatric/neurodevelopmental; ARM B MS as a stressed-glia signature reference (**NOT etiology** — valid claim is "shared stressed-glia program," not causation). These are the BRAIN arms and are **not yet started**. The PLACENTA arm below is complete.

### Method (locked, shared across all Phase-9 arms)
- **Compartment-level pseudobulk RRHO** is the bridge. Mouse and human trophoblast SUBTYPES have no 1:1 homology → compare at the COMPARTMENT level (trophoblast, decidua_stromal, vascular, immune, [erythroid]).
- **Mouse↔human 1:1 ortholog table** `refs/mouse_human_orthologs.tsv` (16,030 pairs, built by `h09e_build_ortholog_map.py` via pybiomart). Mouse stats renamed to human symbols, dedup, dropna before RRHO.
- **Ranking metric = signed DESeq2 Wald `stat`** (= log2FC/SE) — same as 8f/8c, so human and mouse arms are methodologically identical.
- **All RRHO / GSEA / leading-edge functions are lifted VERBATIM** from `08f_cross_tissue.py` (`rrho_matrix`, `classify_rrho_concordance`) and `08c_pathways.py` (`run_gsea_on_ranks` via `run_fgsea.R`, `run_tf_ulm`, `add_fdr`, `load_genesets_tsv`, `compute_leading_edge`). Because `h09*` scripts start with a LETTER they are importable (unlike `08f`/`08c` which start with a digit) — so later h09 scripts import earlier ones directly (`from h09e_cross_species_rrho import rrho_matrix, pseudobulk_de, ...`). Do not duplicate these functions; import them.
- **Human gene sets:** `refs/msigdb_human.tsv` (9,427 sets: H + C2:CP:REACTOME + C5:GO:BP; C8 omitted = human analog of the dropped mouse M8), built by `h_fetch_genesets.R` (msigdbr `db_species="HS"`). Human CollecTRI via `dc.op.collectri(organism="human")`.
- **GSEA concordance design:** two single-species fgsea runs → intersect FDR<0.05 same-sign (cleanest provenance; mirrors 8c). NOT a combined metric.
- **Permutation null:** gene-label shuffle, recompute the concordance-quadrant peak, empirical p (parallel, use_threads=False, chunked).

### ARM 1 — Gunter-Rahman (obesity, GSE271976) — COMPLETE
- Raw `_raw_feature_bc_matrix.h5` only → `h09a` does a knee/inflection cell-call (DropletUtils::barcodeRanks, **`--cutoff inflection`**) before SoupX. **Bug:** `read10xCounts` returns an HDF5 DelayedMatrix → coerce `as(m,"CsparseMatrix")` before SoupX.
- **HVG bug:** placental hormone genes (CGA/CGB/CSH/PSG/PRL) have extreme mean → break seurat_v3 loess → EXCLUDE before HVG (integration only).
- Annotation: marker-majority (`config/human_placenta_markers.yaml`, paper-exact subtypes) → `compartment`; SingleR (Vento-Tormo ref, `h_run_singler.R`, parallelized `--n-jobs 24` MulticoreParam — single-threaded SingleR hangs).
- **Findings:** decidua/vascular/trophoblast RRHO significant (permutation p≤1e-3), immune null (p=0.23, correct); HALLMARK_HYPOXIA concordant-up trophoblast+immune; leading-edge Jaccard ~0.4 (same genes); 0 concordant TFs (honest null). Shared hypoxia genes: SLC2A1, PGK1, PDK1, NDRG1, BNIP3L, ERRFI1.

### ARM 2 — Admati PE 2×2 (figshare 23264102) — COMPLETE
- **Download:** figshare AWS-WAF returns empty 202 with `x-amzn-waf-action: challenge` to plain curl/wget. **Bypass:** `curl -L -A "Mozilla/5.0 ... Chrome/120 ..." -o out.zip "https://ndownloader.figshare.com/files/41003240"` (browser User-Agent passes). File 41003240 = `sc_PE_allcells_with_metadata` (5.9 GB, all compartments, author-annotated).
- **Files:** `sn_..._allcells` is **trophoblast-only** (6,862 cells) despite the name; `sc_admati.zip` = 29 **filtered** Cell Ranger mtx triplets (NO raw → NO SoupX) + `PE_samples_metadata.xlsx`; the figshare sc all-cells txt is the substantive one.
- **h09j** streams the 5.9 GB transposed file (too wide to read whole: ~98k cells in columns) in gene-row chunks, accumulating pseudobulk via a sparse (cells×groups) indicator. Compartment from celltype PREFIX (`TB_`/`STROMAL_`/`VASCULAR_`/`IMMUNE_`). Condition from the 0/1 indicator rows (`early_control`/`late_control`/`early_PE`/`late_PE`). Powering: 10 eoPE / 3 early_control / 7 loPE / 6 late_control donors per compartment.
- **h09k** = the 2×2: mouse {E12.5 Early-vs-Relaxed, E18.5 Late-vs-Relaxed} × human {eoPE, loPE} × 4 compartments. Saves `h09k_rankings.parquet` (so replots/RRHO-maps never recompute), FDR-bearing concordant-pathways CSV, leading-edge CSV. Human DE design `~ condition` (sex not carried; GA confound handled by CONTRAST choice — eoPE vs early_control — not covariate).
- **Diagnostics (`h09k_diagnostics.py`) are mandatory before interpreting.** They established: (a) GA-matched diagonal NOT stronger (rejected); (b) the `discordant` E18.5×loPE×trophoblast cell is a TAIL ARTIFACT (global Spearman positive); (c) loPE≫eoPE is MOSTLY REAL (survives control-downsampling 6→3), partly power.
- **Findings — two conserved axes:** eoPE → HYPOXIA (broad; same genes as Gunter-Rahman: NDRG1/BNIP3L/ERRFI1/PLIN2/ANGPTL4/DDIT4); loPE → OXPHOS/electron-transport suppression (NDUFB4/CHCHD2/COX subunits) + insulin/peptide-hormone. Matches eoPE=hypoxic-placental, loPE=maternal-metabolic pathophysiology.
- **Caveats to state:** sc modality (scRNA vs mouse/Gunter-Rahman snRNA); no SoupX on Admati sc; RRHO peak magnitudes are n-sensitive (compare within column, not across); 3 early_controls make eoPE the thinner contrast.

### Phase 9 datasets
- ARM A brain (downloaded, not processed): Nagy 2020 GSE144136, Maitra 2023 GSE213982 (Mic1=38% female-MDD DEGs), Velmeshev 2019 (UCSC), Herring 2022 GSE168408 (age-anchor, raw). Smoke-test on Velmeshev first.
- ARM B MS (NOT etiology): Macnair 2025 (Zenodo 8338963), Absinta 2021 (GSE180759), Jäkel 2019 (GSE118257).
- ECHO-PATHWAYS (dbGaP phs003619/phs003620) — measured psychological stress, BULK, controlled — revision-stage upgrade.
- **Targeted angle:** does the brain stress signature (ECM/mesenchymal, IFN/immune, gliogenesis threads from 8f/8g) recover in human psychiatric/neurodevelopmental cortex? The brain cross-species result is the publication-tier pivot.

## Smoke-test policy
For every phase >10 min on WS: build a 1-sample/1-pool/1-cluster subset, run, verify, THEN launch full in tmux.
- NEED smoke tests: 1 SoupX, 5 scVI, 7 P1-scANVI, 8c `--tf`, **9 (h09a SoupX-from-raw, h09c scVI, h09k 2×2) — and confirm external file structure on WS before writing any parser.**
- **Phase-9 lesson:** smoke tests don't catch silent-empty joins or transposed-file mis-parses; the real checks are on production data (donor/condition/compartment census non-empty; gene-chunk column-count assertion; `head -1` the output CSV for expected columns).

## Annotation conventions
- **Phase 7 = per-cluster majority voting.** Brain 4-tier, P1 via scANVI (Rosenberg 2018). 4W/3mo via ABC CellTypist. Placenta = markers + STAMP vs Liu 2024 (`celltype_majority`).
- **Phase 9 human placenta:** Gunter-Rahman annotated by marker-majority (`human_placenta_markers.yaml`) + SingleR corroboration; Admati uses the authors' own published labels directly (mapped to compartments by prefix). Rejected training CellTypist on Marsh (Di Bella failure mode).

## Brain marker gate (updated 2026-06-12)
STRICT canonical-marker gates demote borderline calls to `unassigned_*` (Phase 8 DROPS them). microglia ≥2 of {P2ry12,Tmem119,Csf1r,Aif1} (Cx3cr1 REMOVED — off-panel); astrocyte ≥2 of {Aqp4,Gja1,Slc1a3,Aldh1l1}; ol_lineage ≥1 of {Mbp,Mog,Plp1,Mag,Pdgfra,Cspg4,Olig1,Olig2,Sox10}; endothelial ≥2 of {Cldn5,Pecam1,Cdh5}; erythroid ≥2 of {Hbb-bs,Hbb-bt,Hba-a1,Hba-a2,Alas2}. `MARKER_PRESENCE_THRESHOLD=0.20`.

## 8b DE visualization blocklist (locked 2026-06-15)
`BLOCKLIST_FOR_VIZ` in `08b_de_summary.py` (17 hemoglobin + sex-linked) + `BLOCKLIST_PREFIXES=("mt-",)`. `--no-blocklist` for QA.

## Plot format strategy (locked 2026-06-05)
- Default PNG @ 300 DPI (8e/8f/8g + h09 grids 140–150 DPI — many panels).
- Paper figures only: PNG + PDF hybrid via `_utils.savefig`.
- `constrained_layout` fights `bbox_inches="tight"` — `_utils.safe_fig` auto-detects.

## Output organisation
- Tables in per-phase subfolders: `tables/<phase_dir>/<phase>_<name>.csv`.
- Plots in `plots/<phase_dir>/...`.
- **Phase 9:** under `data/human_validation/placenta/<dataset>/`: `tables/h09{a-k}_*.csv` + `h09k_rankings.parquet`; `plots/h09_summary/`, `plots/h09k_admati_2x2/` (incl. `rrho_maps/`).

## Environment specifics
- **uv + renv (not conda).** Bootstrap `./setup-remote.sh`. `statsmodels>=0.14` pinned (8b shuffle + 8f/8g).
- **Phase-9 deps:** `pybiomart`, `decoupler` 2.1.6 + `omnipath`, `adjustText`, `pyarrow` (parquet). R: `SingleR` 2.14.0 (snapshotted). Mirror `renv.lock`/`pyproject.toml`/`uv.lock` back to Mac after env changes.
- **CellTypist sklearn-1.7 patch:** `sed -i "s/multi_class = 'ovr', //;s/multi_class = 'ovr'//" .venv/.../celltypist/train.py` — required before training; `uv sync` reverts it.
- **cuML via `https://pypi.nvidia.com`** for GPU CellTypist training.
- **renv Suggests workaround:** project-level `renv::settings$package.dependency.fields(c("Depends","Imports","LinkingTo"), persist=TRUE)`.

## Workstation infrastructure
- **SSH:** `ssh poller@172.17.213.147`. User `poller`.
- **Mac repo:** `/Users/shuppar/Downloads/BAP_data_1/Analysis/`
- **WS project root:** `/home/poller/BAP-BrainPlacenta/` (NVMe). **WS results mirrored to Mac under `results_WS/`.**
- **Raw tars / Cell Ranger:** `/media/poller/PollerLab-1/BAP-data1/...`; reached via symlink `BAP-BrainPlacenta/data`.
- **Human validation data:** `data/human_validation/{placenta,brain}/<dataset>/`.
- **Pull-from-Mac:** code edits local → rsync to WS (exclude `results/`, `data/`, `.venv*/`, `__pycache__/`, `.git/`, `*.h5ad`, `logs/`, `.DS_Store`; flags `-av --progress --chmod=Fu+x`).
- **WS ↔ Mac MUST stay mirrored.**
- **tmux for any multi-minute job.**
- **Pool → contents:** Pool1 = 16 brain (3mo + 4W males); Pool2 = 16 brain (4W females + P1 Early/Relaxed); Pool3 = 2 brain P1 Late + 14 placenta E12.5; Pool4 = 10 placenta.
- **GPU hygiene:** pre-flight `nvidia-smi`; `PYTORCH_CUDA_ALLOC_CONF=max_split_size_mb:512`; clear cache between GPU phases.

## Documentation
- Update `INSTRUCTIONS.md` and `snRNAseq_project_summary.md` in-place at session end (not via patch sheets).
- Draft handoff prompts at the close of each major phase.
- Standalone findings docs for cross-cutting results (e.g. `08f_08g_cross_tissue_persistence_findings_2026-06-25.md`; a Phase-9 placenta cross-species findings doc is worth writing before the brain arm).
