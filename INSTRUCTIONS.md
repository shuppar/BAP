# Instructions for working on this project with Claude

Get broad context from the `snRNAseq_project_summary.md` file.

## What "update the summary" and "update the instructions" mean (standing definitions)
These two phrases recur — here is exactly what I mean by each:

- **"Update the summary"** = edit `snRNAseq_project_summary.md` (the narrative/findings file)
  **in place** and return the **COMPLETE file** as a drop-in replacement. Surgical edits that
  change only what needs changing (status lines, the relevant phase section, publication
  strategy, next steps) and **preserve everything else verbatim**. Never a patch sheet, diff,
  or "add this snippet" fragment. Verify nothing was dropped (e.g. section-header count before
  vs after).
- **"Update the instructions"** = update THIS instructions block (it is pasted project text,
  NOT a `.md` file you can edit directly), so return the **complete updated instructions text**
  for me to paste, with new rules merged into the right existing sections and **no duplication**
  of rules already present.
- Do both **at session end** as a matter of course, plus draft a **handoff prompt** for the
  next chat at the close of each major phase. Write **complete files**, not paste-blocks, for
  documentation. Standalone findings docs (e.g. `08f_08g_..._findings.md`) for cross-cutting
  results are written BEFORE folding into the manuscript.

## Response style
- **Be brief.** No long preambles, no excessive caveats, no over-explaining.
- **No need to print your thoughts** unless it is helpful to either of us. Keep chats from getting long by not printing reasoning unnecessarily.
- **Don't restate what I just said.** Move to the substance.
- **Step by step.** Build one thing, verify it works, then move on. Don't write 5 files at once. Compile-check before handing a file over.
- **Be honest when something won't work** or when you're uncertain. Don't manufacture confidence.
- **Always give commands to run a specific script** (mention where: Local Mac or remote WS), or rsync a specific file/folder.
- **You (Claude) cannot read the WS filesystem** from your container. For anything on the WS, GIVE the command and ask me to paste the output — never pretend to have run it or guess the result.
- **Plots are instruments to surface trends.** When a figure is meant to reveal trends for further exploration, make sure it actually CAN: a single magnitude metric can structurally hide coherent small-effect programs (the peak-vs-pathway-lens lesson — peak-keyed views are dominated by bulk overlap and hide tight pathway programs like microglial IFN; provide the complementary pathway-keyed lens). Skip empty panels; never plot nulls as if they were signal.

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

## Commands & runs (don't make me ask)
- **Every multi-minute job goes in tmux** — not optional; the default for anything ≳2 min. Give the `tmux new -d -s <name> '...'` line + the `tail -f` follow line.
- **Always give a runtime estimate** with the command, tied to WS vs Mac config.
- **Every command states WHERE it runs** (Mac vs WS) + exact paths. rsync spells out source→dest explicitly.

## Parallelism is mandatory for repeated work (not optional)
Any phase script that loops over samples/jobs/contrasts/cell-types and, per item, launches a subprocess (R worker) or calls a heavy function MUST parallelize via `_utils.parallel_map` and expose `--n-jobs`. A bare for-loop starting one subprocess per item is a performance bug.
- `for item, result, err in parallel_map(fn, items, n_jobs=args.n_jobs, desc="..."):` yields (item, result, error).
- Default `use_threads=True` (subprocess/IO-bound). Use `use_threads=False` (processes) only for CPU-bound pure-Python work; then `fn` and items must be picklable.
- Default `--n-jobs` 8; on the workstation pass 16–24 for light R workers.
- Reference implementation: `08a_composition.py`.
- **CPU-bound exceptions (use_threads=False):** the shuffle null in `08b_disruption_shuffle_test.py`; the permutation nulls in `h09f_overlap_null.py`, `h09k_admati_2x2.py`, and `h10b_brain_rrho.py` (each chunk does a batch of RRHO shuffles; arrays pickled once per chunk, NOT per shuffle).
- **Plotting note:** 8e/8f/8g and all `h09*`/`h10*` summary plotting is intentionally SERIAL (matplotlib not thread-safe).
- **8f/8g parallelization caveat (2026-06-25):** single-threaded across slices; process-based parallel_map would pickle the 20M-row brain DE frame to every job. Pre-filter each slice serially, then parallel_map the small frames if ever needed.

## Pipeline architecture decisions (don't re-litigate)
Don't re-open locked decisions. If something genuinely needs revisiting, say WHY explicitly and ask — don't silently re-open it.
- **Language:** Python primary (Scanpy/scvi-tools). R subprocess for scDblFinder, propeller, SoupX, fgsea, SingleR.
- **Env:** uv + Python 3.12. **Conda is blocked** at corporate firewall — don't suggest it.
- **scVI**: GPU phase, runs on workstation.
- **Phase 1 ambient correction = SoupX (locked 2026-06-10).** CellBender abandoned (pickle bug).
- Don't propose architectural pivots (containers, Nix, etc.) — uv + scripts is settled. Don't add abstraction layers "for future flexibility." YAGNI.

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
- **Outputs:** `results/{tissue}/<phase>/` for prod, `results/dev/<phase>/` for laptop tests. **Phase-9 human outputs live under `data/human_validation/{placenta,brain}/<dataset>/{h5ad,tables,plots}/`.**

## When asking me questions
- **One question at a time** usually; **3 max**.
- **Single-select buttons** over multi-select where possible.
- **Don't ask for confirmation on small obvious things** — just do them.
- **Brainstorm → lock the decision (single-select) → then code.** Don't bake a strong scientific call in silently — surface it.

## Surprising results: interrogate, don't interpret-or-"fix"
- When a result is sparse, huge, or surprising, **run a diagnostic FIRST** (peak-vs-strength, quadrant decomposition, control downsampling, direction check) and establish artifact-vs-real BEFORE building interpretation or "fixing" the code. Established pattern: placenta Admati anomalies; the brain P1 / Maitra-213 / MS-microglia checks.
- A **blank/null cell often means "filtered for a correct reason"** (e.g. discordant direction), not "nothing there" — check the underlying ranking before concluding absence.

## What I don't want
- Don't propose architectural pivots (containers, Nix, etc.) — uv + scripts is settled.
- Don't add abstraction layers "for future flexibility." YAGNI.
- Don't write 200-line responses with 5 nested headers when 30 lines suffices.
- Don't apologize repeatedly when correcting something.

## Naming conventions (locked 2026-06-25)
- **Mouse pipeline:** numbered phase scripts `0X_*.py` / `08x_*.py`; R workers `run_*.R`.
- **Human cross-species — PLACENTA (Phase 9):** scripts `h09X_...py` (`h09a`…`h09k`, plus `h09_summary_plots.py`, `h09k_diagnostics.py`, `h09k_plots.py`, `h09k_rrho_maps.py`); R workers `h_run_*.R` / `h_fetch_*.R`.
- **Human cross-species — BRAIN (Phase 9):** scripts `h10X_...py` (`h10a` prep Velmeshev, `h10b` engine, `h10c/d/e` prep Maitra/Nagy/Macnair, `h10b_diagnostics.py`, `h10b_rrho_maps.py`, `h10_summary_plots.py`); shares the `h_run_*.R`/`h_fetch_*.R` workers. Keeps the human arms visually separate from the numbered mouse pipeline and from each other (h09=placenta, h10=brain).

## External identifiers: verify or flag
Several bugs came from writing external identifiers from memory as if verified. **Any external identifier (PyPI/conda package, model name, gene symbol, API function/kwarg, file/dataset/COLUMN name) must be verified against docs/data or flagged.** When unsure, give a WS command to check FIRST, parser SECOND.
- Cautionary examples: `Mouse_Brain_Atlas` (doesn't exist); `abc-atlas-access` (real is `abc_atlas_access`); `score_genes(use_raw=False)` runs on raw `.X`; `Mlf1ip`/`Fam64a`/`Hn1` outdated (`Cenpu`/`Pimreg`/`Jpt1`); `multi_class='ovr'` removed in sklearn 1.7; `Cx3cr1` NOT in 10x Flex Mouse v2 panel (use P2ry12/Tmem119/Csf1r/Aif1); `datasplitter_kwargs` is the scvi-tools 1.4.3 kwarg.
- **Contrast names are NOT uniform across tissues:** brain `*_per_age` (age in `group_level`), placenta `*_E12.5`/`*_E18.5`. Match by family prefix (`contrast_family`), NEVER exact string.
- **`08b_developmental_disruption_genes.csv` stores its class in a column literally named `direction`** (NOT `direction_class`).
- **Phase 9 (placenta):** `unassigned_mask(obs, label_cols)` takes TWO args; Admati files are TRANSPOSED (metadata rows then gene rows, cells in columns); Admati gene IDs are human SYMBOLS (no Ensembl map); the `sn_..._allcells` file is trophoblast-ONLY despite the name; the figshare file ID for the sc all-cells matrix is **41003240** (5.9 GB). Always `cut -f1 | head` to confirm metadata-row vocabulary before parsing.
- **Phase 9 (brain):** Velmeshev `exprMatrix.tsv.gz` is TRANSPOSED-wide (genes in rows, cells in header) + gene IDs are fused `ENSG|SYMBOL` (split on `|`); Maitra combined matrix has BOTH sexes — filter `F#`, and no condition is in any matrix file (needs `maitra_donor_meta.csv`); the M#→Nagy crosswalk is unrecoverable → use Nagy standalone; Nagy barcode subtype segment is OPTIONAL (Micro/Macro + Endo lack it — a too-strict regex silently drops microglia + endothelial); Macnair `np.where(str, str, np.nan)` raises DTypePromotionError on modern numpy → use `.map(dict)`; mouse 08b genes are mouse symbols, ortholog cols `mouse_symbol`/`human_symbol`; 08b celltype `OPC/Oligodendrocytes` is one class vs human Oli/OPC (merged-primary); 08b `Isocortex` level carries ONLY neurons by construction (T2=neurons by data, not choice).

## Confirm external file STRUCTURE before parsing
Files lie. Before writing any parser, confirm structure on the WS — files can be transposed, "allcells" can be a subset, gene IDs vary (Ensembl vs symbol vs fused), metadata can be inline in barcodes, mtx can be gene-major or cell-major. `head`/`zcat | head`/dim-check first, parser second. Smoke tests don't catch silent-empty joins or transposed-file mis-parses; the real checks are on production data (donor/condition/compartment census non-empty; gene-chunk column-count assertion; `head -1` the output CSV for expected columns).

## No silent failures
Wrong-but-plausible output is worse than a crash. **A correctness-critical step that can't run correctly must raise, not warn-and-continue.** Warn-and-skip is only acceptable for an OPTIONAL output and the skip must be announced. Never leave NaN labels surfacing as a "nan" category.
- **Phase-9 guards:** `h09j`/`h10a` assert each gene chunk's cell-column count matches the obs row count (catches ragged parsing); `h09k`/`h09e`/`h10b` skip a compartment/celltype loudly when a group has too few donors; the human-DE helpers require both conditions ≥ MIN_DONORS.

## The internal-comma trap (CSV parsing)
Several Phase-8 CSVs have a `pair` column holding a Python list (internal comma inside quotes). **Never `cut`/`awk` these — always pandas.** For big mixed-dtype tables use `low_memory=False`; for the 6.7 GB leading-edge file read chunked (`chunksize=2_000_000`) with `usecols`.

## Plots must carry biological meaning (not abstract designs)
Every figure should let a reader name the biology without a side table.
- Label what matters (volcano → top sig genes; GSEA → pathways; heatmaps → real names). Cap labels (~top 25).
- Gene identifiers human-readable (`var['symbol']`). State contrast + thresholds on the plot.
- **Don't apply an effect-size floor to distribution plots** (volcano, rank-rank scatter, 8f LR scatter, h09 leading-edge scatters). Floors apply only to ranked/aggregating plots.
- **adjustText for label de-overlap** (`uv add adjustText`); `h09_summary_plots._label_points()` is the shared helper (leader lines + plain-annotate fallback). Reuse it in `h09k_plots.py` and `h10_summary_plots.py`.
- **Per-panel scaling** when a shared color scale would crush secondary panels (the RRHO-map lesson — a global vmax driven by one peak-232 cell makes every other panel look empty; scale each panel to its own max, put the magnitude in the title).
- **Peak-keyed vs pathway-keyed lenses are both needed.** Peak/magnitude views answer "where is the bulk gene overlap" (→ dominated by neurons); pathway-keyed views (GSEA-FDR-gated, NOT peak-gated) answer "which coherent programs are conserved" (→ surface IFN/ECM/gliogenesis threads). Keep both; don't let a peak filter hide a tight pathway program.
- **Phase-9 plotting is CSV/parquet-only and Mac-runnable.** Each plot script reads the saved tables (and the rankings parquet for RRHO maps) — never recompute DE. When a plot needs FDR-sized dots, the source table must carry FDR; if a CSV predates that, re-run the producing script once.

## Ask before strong scientific calls
- Don't drop an analysis, exclude samples/ages, or add complexity without checking it earns its place. Surface the question; don't bake it in silently.
- Don't propose: cross-tissue CCC; RNA velocity / CellRank; scCODA.
- When a result is sparse or surprising, interrogate before accepting OR fixing: cutoff artifact, silent-filter bug, or real biology? **Phase-9 example: the Admati 2×2 anomalies were interrogated via `h09k_diagnostics.py` BEFORE interpretation; the brain P1/Maitra-213/MS-microglia anomalies via `h10b_diagnostics.py`. Always run the diagnostic before building the figure.**

## Phase 8 conventions
(unchanged — cross-cutting + 8f/8g-specific + stage-specific rules; statistical unit = animal; declarative contrasts/strata; drop contaminants+unassigned; `~ sex + pool + group`; `group_level` holds age for `*_per_age`; always filter sex+level before classifying; corrected p-values everywhere; GSEA = mouse MSigDB via msigdbr MH+M2+M5; TF activity REQUIRED in 8c production.)

## Phase 8f / 8g conventions (locked 2026-06-25)
(unchanged — six 8f views from 8b/8c CSVs; RRHO vectorized; `--logfc-cutoff 0.5`; plot-only `--plot-quantile 0.75`. 8g brain only; persistence classes per (celltype, level, feature, arm); B trajectory-shape / C persistence×disruption / View-7 8f-bridge modules. Interpretation locked: 0 persistent genes; IFN perinatal-transient; ECM durable; gliogenesis core.)

## Disruption analysis framing (8b, locked 2026-06-15)
(unchanged — "When age-DE signal is shared across two groups, Relaxed is almost always one of them.")

## Phase 9 — cross-species validation: PLACENTA arm (COMPLETE 2026-06-25)

**Two scientific arms reported separately, NOT pooled** (locked 2026-06-05): psychiatric/neurodevelopmental vs MS-as-stressed-glia-reference (**NOT etiology** — valid claim is "shared stressed-glia program," not causation).

### Method (locked, shared across all Phase-9 arms — placenta AND brain)
- **Compartment/broad-celltype-level pseudobulk RRHO** is the bridge. Placenta: mouse/human trophoblast SUBTYPES lack 1:1 homology → compare at COMPARTMENT level (trophoblast, decidua_stromal, vascular, immune, [erythroid]). Brain: types ARE homologous → broad 7-class.
- **Mouse↔human 1:1 ortholog table** `refs/mouse_human_orthologs.tsv` (16,030 pairs, `h09e_build_ortholog_map.py` via pybiomart). Mouse stats renamed to human symbols, dedup, dropna before RRHO.
- **Ranking metric = signed DESeq2 Wald `stat`** — same as 8f/8c, so human and mouse arms are methodologically identical.
- **All RRHO / GSEA / leading-edge / TF functions are lifted VERBATIM** from `08f`/`08c` (and re-exported through `h09e`/`h09g`/`h09h`/`h09k`). Later scripts IMPORT them, never duplicate (h09*/h10* start with a letter so they're importable).
- **Human gene sets:** `refs/msigdb_human.tsv` (H + C2:CP:REACTOME + C5:GO:BP), `h_fetch_genesets.R`. Human CollecTRI via `dc.op.collectri(organism="human")`.
- **GSEA concordance:** two single-species fgsea → intersect FDR<0.05 same-sign. NOT a combined metric.
- **Permutation null:** gene-label shuffle, recompute the concordance-quadrant peak, empirical p (parallel, use_threads=False, chunked).

### Placenta findings (brief)
ARM 1 Gunter-Rahman obesity: decidua/vascular/trophoblast RRHO significant (p≤1e-3), immune null; HALLMARK_HYPOXIA concordant-up trophoblast+immune; leading-edge Jaccard ~0.4; 0 concordant TFs. ARM 2 Admati PE 2×2: GA-matched diagonal REJECTED; two conserved axes — eoPE→HYPOXIA (same genes), loPE→OXPHOS/electron-transport suppression. Diagnostics mandatory (`h09k_diagnostics.py`).

## Phase 9 — cross-species validation: BRAIN arm (h10*, COMPUTE-COMPLETE 2026-06-29)
Four independent datasets, same engine, reported separately. Naming `h10*`. Method = the shared Phase-9 method above; bridge = broad 7-class (brain types ARE homologous).

### Datasets + units (each unit decision by its own structure, NOT a blanket rule)
- **Velmeshev 2019** (ASD, UCSC) — transposed-wide TSV, `ENSG|SYM`→symbol, align-by-position. Unit=sample (indiv×region PFC/ACC). Primary drops Neu-NRGN/Neu-mat; sensitivity variant adds them→ExN.
- **Maitra 2023** (MDD-F, GSE213982) — mtx triplet; filter `F#`; condition from `maitra_donor_meta.csv`. Unit=donor (all BA9). Mix dropped.
- **Nagy 2020** (MDD-M, GSE144136) — mtx triplet; barcode `{prefix}[_subtype].{donor}_{cond}_{batch}_{bc}` (subtype OPTIONAL); donor=`.NN` (17/17), condition inline (Suicide→MDD). Unit=donor. Mix dropped.
- **Macnair 2025** (MS, Zenodo 8338963) — discovery mtx triplet; `col_data` carries donor/`type_broad`/`diagnosis`/`matter`/`exclude_pseudobulk` (honor the flag); MS={SPMS,PPMS,RRMS} vs CTR; genes from `row_data` `symbol`; B/T cells dropped. Unit=donor (GM/WM has NO mouse analog → donor; matter kept descriptive only).

### Engine (h10b)
- Mouse rankings = **reuse 08b stat** (no recompute; consistent with paper Figs 2/3). Filter `contrast∈{early_vs_relaxed_per_age, late_vs_relaxed_per_age}`, `sex=combined`, `level∈{whole,Isocortex}`, 7-broad celltype map. Bridge mouse→human symbols.
- Human DE = PyDESeq2 from parquet read directly (mirrors `h09k.human_rankings`); `~ [sex+] diagnosis`; covar auto-drops when constant (Maitra-F/Nagy-M → `~ diagnosis`; Macnair keeps sex).
- Full 3×2 mouse grid (P1/4W/3mo × ES/LS) × {whole, Isocortex} × 7 broad; Isocortex→neurons only by data. Oli/OPC merged-primary (mouse `Oli_OPC` → human Oli AND OPC separately).
- `--tf` opt-in (CollecTRI human ULM); BRAIN TFs are NON-null (run with `--tf`). Saves `h10b_<ds>_rankings.parquet` so maps/plots never recompute.

### robust_class guard (locked 2026-06-29)
The argmax RRHO label flips on noise and can contradict the global Spearman sign. `robust_class` = directional label only if (margin ≥ 25% of peak AND Spearman sign agrees: concordant→r>0, discordant→r<0), else `ambiguous`. Columns `concordance_margin`, `spearman_agrees`. **Shuffle null runs on every cell regardless** — magnitude always tested, direction withheld when fragile. The RRHO MAP is ground truth; the label is summary-only.
> Backport this guard to the placenta arm at revision (currently margin-free argmax there; findings stand but the rule must be consistent across Fig 4).

### Diagnostics MANDATORY before interpretation
`h10b_diagnostics.py --dataset <ds>` — peak-vs-vector-strength (mechanical?), quadrant decomposition (tail artifact?), P1/4W/3mo per-celltype (age-specific vs across-the-board?). Run BEFORE building any interpretation or figure.

### Plotting (h10_summary_plots.py) — BOTH lenses
- **Peak-keyed** (overview, master heatmap, RRHO maps): filter `empirical_p<0.05`; robust_class as VISUAL encoding (directional solid/boxed, ambiguous light), not a gate. Per-panel map scaling.
- **Pathway-keyed** (thread-scanner: IFN/immune, ECM/mesench, gliogenesis, synaptic; Mic/IFN panels): gate on GSEA FDR<0.05 (the concordant-pathway definition), NOT on RRHO peak — this is the lens that surfaces small-but-coherent programs the peak views hide.
- Empty panels skipped. CSV/parquet-only, Mac-runnable, serial.

### Findings (locked, diagnostics-confirmed)
Neurons dominate by peak in all four. **MDD (both sexes) = headline: strong directional neuron-DOWN concordance at mouse 4W** (Maitra ExN 232 / InN 213, Spearman ~0.3; Nagy replicates). Female-MDD adds Mic, male-MDD adds Oli. **ASD weaker/bidirectional/ambiguous (P1).** **MS weak/glia-leaning.** **IFN/immune thread recovered as MICROGLIAL co-suppression in MDD/ASD (down_both, P1-prominent — the 8g signature) and INVERTS to up in MS** (MS microglia inflamed median +1.05 vs mouse stress −1.15 → legitimately discordant → the "MS≠etiology" directional control). ECM/mesenchymal (164 rows) + gliogenesis (34 rows) threads also present. **Interpretation:** mouse prenatal-stress neuronal programs converge most strongly on human MDD cortex; IFN thread microglial + inverts in MS; ECM/gliogenesis present; ASD weaker, MS distinct → disorder-specificity. The Fig 4 brain anchor.

### CLI (WS, from project root)
```bash
uv run python scripts/h10a_prep_velmeshev.py            # ~3-6 min (transposed stream, 2 variants)
uv run python scripts/h10c_prep_maitra.py               # ~3-6 min
uv run python scripts/h10d_prep_nagy.py                 # ~2-3 min
uv run python scripts/h10e_prep_macnair.py              # ~5-10 min (632K cells)
uv run python scripts/h10b_brain_rrho.py --dataset <ds> --tf --n-perm 5000 --n-jobs 16  # ~25-45 min, tmux
uv run python scripts/h10b_diagnostics.py --dataset <ds>   # seconds; MANDATORY
uv run python scripts/h10_summary_plots.py              # ~1-2 min, Mac-runnable
uv run python scripts/h10b_rrho_maps.py --dataset <ds>  # optional full grid
```

### Brain remaining work
(a) Fig 4 figure refinement — panels exist (both lenses); iterate any that don't tell the story (placenta Fig 4 was flagged "not yet representative" — apply that scrutiny). (b) Standalone brain findings doc before manuscript. (c) Velmeshev sensitivity variant (`--variant sensitivity`) into a quarantined subfolder + README. (d) Optional Herring 2022 age-anchor. (e) Backport robust_class to placenta (revision-stage). (f) ECHO-PATHWAYS dbGaP (measured stress, revision-stage upgrade).

## Smoke-test policy
For every phase >10 min on WS: build a 1-sample/1-pool/1-cluster subset, run, verify, THEN launch full in tmux.
- NEED smoke tests: 1 SoupX, 5 scVI, 7 P1-scANVI, 8c `--tf`, 9 placenta (h09a/h09c/h09k), 9 brain (h10a prep, h10b engine) — and confirm external file structure on WS before writing any parser.
- **Phase-9 lesson:** smoke tests don't catch silent-empty joins or transposed-file mis-parses; the real checks are on production data (donor/condition/compartment census non-empty; gene-chunk column-count assertion; `head -1` the output CSV for expected columns).

## Annotation conventions
- **Phase 7 = per-cluster majority voting.** Brain 4-tier, P1 via scANVI (Rosenberg 2018). 4W/3mo via ABC CellTypist. Placenta = markers + STAMP vs Liu 2024 (`celltype_majority`).
- **Phase 9 human placenta:** Gunter-Rahman by marker-majority + SingleR; Admati uses authors' own labels by prefix.
- **Phase 9 human brain:** all four datasets use authors' own published labels mapped to the broad 7-class (no re-clustering; pseudobulk on their annotations). Velmeshev fine clusters → broad via the locked map (primary drops Neu-NRGN/Neu-mat).

## Brain marker gate (updated 2026-06-12)
STRICT canonical-marker gates demote borderline calls to `unassigned_*` (Phase 8 DROPS them). microglia ≥2 of {P2ry12,Tmem119,Csf1r,Aif1} (Cx3cr1 REMOVED — off-panel); astrocyte ≥2 of {Aqp4,Gja1,Slc1a3,Aldh1l1}; ol_lineage ≥1 of {Mbp,Mog,Plp1,Mag,Pdgfra,Cspg4,Olig1,Olig2,Sox10}; endothelial ≥2 of {Cldn5,Pecam1,Cdh5}; erythroid ≥2 of {Hbb-bs,Hbb-bt,Hba-a1,Hba-a2,Alas2}. `MARKER_PRESENCE_THRESHOLD=0.20`.

## 8b DE visualization blocklist (locked 2026-06-15)
`BLOCKLIST_FOR_VIZ` in `08b_de_summary.py` (17 hemoglobin + sex-linked) + `BLOCKLIST_PREFIXES=("mt-",)`. `--no-blocklist` for QA.

## Plot format strategy (locked 2026-06-05)
- Default PNG @ 300 DPI (8e/8f/8g + h09/h10 grids 140–150 DPI — many panels).
- Paper figures only: PNG + PDF hybrid via `_utils.savefig`.
- `constrained_layout` fights `bbox_inches="tight"` — `_utils.safe_fig` auto-detects.

## Output organisation
- Tables in per-phase subfolders: `tables/<phase_dir>/<phase>_<name>.csv`.
- Plots in `plots/<phase_dir>/...`.
- **Phase 9 placenta:** under `data/human_validation/placenta/<dataset>/`: `tables/h09{a-k}_*.csv` + `h09k_rankings.parquet`; `plots/h09_summary/`, `plots/h09k_admati_2x2/`.
- **Phase 9 brain:** under `data/human_validation/brain/<dataset>/`: `tables/h10{a-e}_*` + `h10b_<ds>_*.{csv,parquet}`; `plots/{h10b_rrho_maps,h10_summary}/`; cross-dataset `data/human_validation/brain/_synthesis/plots/` (01_overview, 02_master, 03_thread_scanner, 04_microglia_ifn, 05_ifn_all_celltypes).

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
- Update the project summary (`snRNAseq_project_summary.md`) and these instructions in-place at session end (not via patch sheets) — see "What update means" at the top.
- Draft handoff prompts at the close of each major phase.
- Standalone findings docs for cross-cutting results (e.g. `08f_08g_cross_tissue_persistence_findings_2026-06-25.md`; a Phase-9 placenta and a Phase-9 brain cross-species findings doc each worth writing before the manuscript).
