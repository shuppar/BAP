# snRNA-seq Analysis Project: Prenatal Stress in Mouse Brain & Placenta

## 1. Study Overview

### Experimental Design

**Model:** Mouse model of prenatal stress with **three groups**:
- **Early Stress** — maternal stress applied during early gestation
- **Late Stress** — maternal stress applied during late gestation
- **Relaxed** — unstressed controls

**Tissues & Timepoints (actual, from sample manifest):**

| Tissue | Timepoints | Total samples |
|---|---|---|
| Brain | P1, 4 weeks, 3 months | 34 |
| Placenta | E12.5, E18.5 | 23 (after removing duplicate CES2.3) |

**Sample allocation (brain, n=34):**

| Age | Early Stress | Late Stress | Relaxed |
|---|---|---|---|
| P1 | 2M + 2F | 2M (no F) | 2M + 2F |
| 4-week | 2M + 2F | 2M + 2F | 2M + 2F |
| 3-month | 2M + 2F | 2M + 2F | 2M + 2F |

**Sample allocation (placenta, n=23):**

| Age | Early Stress | Late Stress | Relaxed |
|---|---|---|---|
| E12.5 | 9 (sex TBD) | — | 6 (sex TBD; was 7, dropped 1 duplicate) |
| E18.5 | — | 2M + 2F | 2M + 2F |

**Critical design notes:**
- **Three groups, not two** — this changes every contrast (was binary stress vs control, now Early-vs-Relaxed and Late-vs-Relaxed as primary contrasts)
- **P1 Late Stress has no females** (only 2 males) — these samples are also confounded with Pool 3 (see §2)
- **Placenta has no complete factorial** — E12.5 has Early+Relaxed only, E18.5 has Late+Relaxed only. Cross-age placenta comparisons are not analyzable.
- **Placenta E12.5 sex is undetermined** at sampling — will be inferred from Y-chromosome expression
- **No dam ID recorded** — pups from the same litter cannot be identified as such, so the litter random effect cannot be modeled (see §2)

**Sequencing pools (= libraries = scVI batch_key):**

| Pool | Composition |
|---|---|
| Pool 1 | 16 brain samples (3-month + part of 4-week + male 4-week) |
| Pool 2 | 16 brain samples (P1 Early+Relaxed + part of 4-week + female 4-week) |
| Pool 3 | 2 brain samples (P1 Late Stress only) + 14 placenta E12.5 samples |
| Pool 4 | 10 placenta samples (2 E12.5 Relaxed + all E18.5) |

**Data type:** 10x Genomics **Flex** chemistry, processed through Cell Ranger multi
- Input files: `.h5` files (filtered + raw matrices per sample, in `per_sample_outs/<sample_id>/`)
- **Pre-prep nuclei counts: 0.6M–25M per sample** (cell-counter measurement after nuclei isolation)
- **Cell Ranger-called cells: ~10K–30K per sample** (post-capture, post-Cell-Ranger cell-calling) — TO BE CONFIRMED from `metrics_summary.csv`
- Total cells post-QC: estimated ~500–700K brain, ~300–400K placenta

**Brain region:** Whole brain (not microdissected)

---

## 1b. Implementation Status

**All phases through 8g implemented and smoke-tested on dev (laptop, 3 samples × 3 pseudo-donors × 500 cells/donor).** Workstation execution pending; cross-tissue (8f) and cross-age (8g) views designed to operate on completed 8b/8c CSVs so they need no re-runs.

> **Workstation target** (where production runs go; see §3 for detail):
> Linux box, **258 GB RAM, 56 CPU cores, 1× NVIDIA RTX 4500 Ada (24 GB VRAM)**.
> GPU and CPU compute on the same host. Conda blocked at firewall — use `uv` + `renv` only. R + Rscript on PATH. **CellBender abandoned 2026-06-05 (pickle bug); replaced by SoupX via R subprocess (2026-06-10).** Everything runs from the main uv-managed `.venv/` (no sidecar venvs). cuML installed via NVIDIA PyPI for GPU LogReg (CellTypist training, locked 2026-06-10). Runbook: `run_pipeline_WS.sh`.

| Phase | Status | Script(s) |
|---|---|---|
| 0 Validation | ✓ done | `01_validate.py` |
| 1 Ambient RNA (SoupX) | ✓ scripts written, smoke-test pending (2026-06-10) | `02_soupx.py` + `run_soupx.R` |
| 2 Per-sample QC | ✓ done (TBD: wire prefer-soupx fallback) | `02_qc.py` |
| 3 Doublet detection | ✓ done (per-pool via R subprocess) | `03_doublets.py` + `run_scdblfinder.R` |
| 4 Concat + HVG + cell cycle | ✓ done | `04_integration_prep.py` |
| 5 scVI integration | ✓ done (CPU dev; GPU workstation, BF16) | `05_integration.py` |
| 6 Clustering (Leiden, igraph) | ✓ done | `06_clustering.py` |
| 7 Annotation (per-cluster majority) | ✓ done | `07_annotation.py` |
| 7b Subclustering | ✓ done | `07b_subcluster.py` |
| 7d Subcluster annotation | ✓ done | `07d_subcluster_annotate.py` + `config/subcluster_markers.yaml` |
| 7c scANVI label transfer | ✓ done (workstation-only; exits cleanly on dev) | `07c_label_transfer.py` + `prepare_reference.py` |
| 8a Composition (propeller) | ✓ done (R subprocess) | `08a_composition.py` + `run_propeller.R` |
| 8b Pseudobulk DE (PyDESeq2) | ✓ done + per-sample expr matrix offline audit | `08b_de.py` |
| 8c GSEA + leading-edge + TF activity | ✓ done (decoupler ULM + CollecTRI) | `08c_pathways.py` + `fetch_genesets.R` |
| 8d Trajectory (PAGA + DPT) | ✓ done + edge diagnostics offline audit | `08d_trajectory.py` |
| 8e Cell-cell communication | ✓ done (LIANA+ baseline + differential + per-donor) | `08e_communication.py` + `_08e_plots_baseline.py` + `_08e_plots_differential.py` + `_08e_plots_perdonor.py` |
| 8f Cross-tissue | ✓ done (six views; placenta→brain cascades) | `08f_cross_tissue.py` |
| 8g Cross-age / persistence | ✓ done (six views; derived from 8b/8c) | `08g_cross_age.py` |

**Dev pre-step (one-time, NOT part of workstation pipeline):**

| Tool | Purpose | Where it lives |
|---|---|---|
| `dev_split_h5.py` | Reads 3 dev h5 files, writes 9 split h5 files + `config/dev_split.yaml` (random barcode partition, donor_id suffixed _ps1/_ps2/_ps3). All downstream phases run unchanged with `--config config/dev_split.yaml`. Pseudo-donors are random partitions of one animal → numbers MEANINGLESS, smoke test only. | repo root |

**Key implementation notes (may differ from the original plan in §5/§6):**

1. **Flat scripts layout, no `src/snrna/` package.** Per `INSTRUCTIONS.md`. Shared helpers live in `scripts/_utils.py` (`load_config`, `add_lognorm`, `phase_paths`, `phase_table_dir`, `select_accelerator`). Phase scripts are numbered standalone files in `scripts/`.

2. **CellBender skipped on laptop dev.** Phase 1 is GPU-only; the dev pipeline goes directly from raw Cell Ranger h5 → Phase 2. CellBender will run on the workstation in production.

3. **Raw counts in `.X`, lognorm not persisted.** Phase 4 builds `.layers["lognorm"]` for Phase 5's pre-integration UMAP, then Phase 5 drops it before writing the final h5ad. Notebooks recompute via `_utils.add_lognorm(adata)`. Saves ~50% disk on the 600K-cell production object.

4. **QC: per-sample MAD + hard floors + hard caps + cohort-outlier flag.** Each cell must pass per-sample MAD bounds AND the absolute floors `min_counts=500`, `min_genes=200` AND hard caps `pct_mt≤1.0`, `pct_hemo≤5.0`. A cohort-outlier flag (sample median UMI/genes >3 cohort-MADs below median) is added when n≥5, catching failed-prep samples whose own MAD made their bounds permissive.

5. **scDblFinder per pool, not per sample.** Per-pool matrix combined, passed to scDblFinder with `samples=` arg so simulated doublets respect within-sample boundaries. Cells classified `doublet` are removed (not just flagged).

6. **HVG exclusion lists implemented.** `var["use_for_scvi"]` = `highly_variable AND NOT hvg_excluded`. Generic exclusions: mito, ribo, hemoglobin, sex-linked (Xist, Ddx3y, Uty, Eif2s3y, Kdm5d, Tsix). Placenta adds Prl*, Psg*, Cgb*, Cga.

7. **scVI uses BF16 mixed precision on GPU, FP32 on CPU.** Auto-detected via `_utils.select_accelerator()`. Dev runs on tiny data (<5K cells) auto-cap `max_epochs` at 50.

8. **Leiden uses igraph backend** with `flavor="igraph", n_iterations=2, directed=False` (scanpy's future default, faster than leidenalg).

9. **Phase 7 annotation uses per-cluster majority voting**, not per-cell argmax. Cells in one Leiden cluster share a label; low-purity (<60% majority) clusters announced in stdout for manual review. Phase 7d (subcluster naming) is already cluster-level by construction.

10. **scANVI implemented but workstation-only** (`07c_label_transfer.py`). Requires `reference:` block in YAML + a labeled reference h5ad. Exits cleanly on dev with `"no 'reference:' block"` — by design. The same Allen Brain Cell Atlas build also trains a CellTypist `.pkl` for adult 4W/3mo brain (no built-in adult mouse CellTypist model exists).

11. **scCODA abandoned for composition (8a).** Its TF/TFP/arviz dependency stack fought the scanpy/scvi-tools stack. Replaced with **propeller (speckle+limma) via R subprocess** — clean Bioconductor install, limma's empirical-Bayes moderation is better for small n anyway.

12. **TF activity (8c) via decoupler ULM + CollecTRI mouse network.** Runs on the same DE Wald-stat vector GSEA uses. BH-FDR within celltype×contrast. Plots: per-celltype barplot + volcano, per-contrast TF×celltype heatmap (significant TFs only, capped at top N by |score|). Enable with `--tf` or YAML `pathways.run_tf_activity: true`. Needs network (omnipath).

13. **No RNA velocity, no CellRank.** 10x Flex is probe-based (exon-only), can't resolve spliced/unspliced. Without velocity CellRank only duplicates PAGA. Trajectory = PAGA + diffusion pseudotime. All ages treated identically in DPT (no gating off "steady-state" ages).

14. **Dev workflow: `dev_split_h5.py`** runs once before Phase 0, writes 9 split h5 files + `config/dev_split.yaml`. NO pipeline scripts are dev-aware — same scripts run on workstation with `--config config/brain.yaml`. The split round-trips cleanly through `sc.read_10x_h5` (verified by writer self-check). The old `dev_pseudoreplicate.py` (late mutation of donor_id before 8a) is OBSOLETE.

15. **Output organisation: per-phase subfolders.** Tables in `tables/<phase_dir>/<phase>_<name>.csv`, plots in `plots/<phase_dir>/...`. Helper: `_utils.phase_table_dir(cfg, label)`. Cross-phase reads (e.g. 8c reads 8b's `de_results`) updated to use the per-phase paths.

16. **Offline-audit CSVs.** Built into 8b/8c/8d to make results troubleshootable without workstation access: per-sample DE-gene expression matrix (8b), pathway leading-edge genes with log2FC (8c), TF activity (8c), PAGA edge diagnostics (8d). Join keys: (celltype, gene) bridges 8c leading-edge ↔ 8b expression matrix.

17. **8e cell-cell communication: LIANA+ in main env, no sidecar.** liana-py 1.7.3 via `uv add liana`. Three arms in one script: (1) baseline `rank_aggregate` per group×age, (2) differential via `li.multi.df_to_lr` reading 8b's Wald stats, (3) per-donor `rank_aggregate` → Mann-Whitney across donors. Covers all three group comparisons (ES-v-Rel, LS-v-Rel, ES-v-LS) explicitly. Clustered Δ heatmaps (seaborn clustermap, pathway+celltype colour bars), rank-rank concordance scatters answering "do ES and LS hit the same programs?", optional `--zscore-rows` for pattern view alongside absolute Δ.

18. **8e subcluster handling: separate output trees.** `--subcluster <slug>` loads `08c_subclustered/{slug}.h5ad` and writes to `plots/08e_communication_subcluster_{slug}/` and `tables/08e_communication_subcluster_{slug}/` — no collision with main run. Cell-type column auto-detected via priority list `[subcluster_name, subcluster, manual_annotation, scanvi_celltype, celltypist_majority, provisional_celltype]`.

19. **8f cross-tissue: six views, all reproducible from 8b/8c CSVs.** Two biologically aligned arms (E12.5 placenta Early → P1/4W/3mo brain Early; E18.5 placenta Late → P1/4W/3mo brain Late; P1 Late carries `confounded_with_pool` flag). Views: (1) hypergeometric DEG overlap, (2) RRHO (custom NumPy ~30 lines), (3) pathway concordance from 8c GSEA, (4) LR cross-tissue mechanistic hypotheses (placental ligand × brain receptor from liana mouseconsensus; `stress_axis` column flags GR/MR/CRH/cytokine genes from a curated 60-gene list — Goeden/Vacher/Wu/Bonnin axes), (5) TF concordance from 8c TF activity, (6) ORA of overlap genes vs MSigDB. **NO cross-tissue CCC** — within-tissue interaction scoring doesn't extend across BBB; view 4 (LR from DE) is the correctly framed placenta-brain endocrine/paracrine version. Effect-size scatters (Spearman ρ) sit alongside heatmaps for genes, pathways, and TFs.

20. **8g cross-age persistence: brain-only by design.** Placenta has incomplete cross-age factorial (E12.5 = Early+Relaxed, E18.5 = Late+Relaxed; no factorial across ages), so 8g exits cleanly with `tissue: placenta`. Six views: (1-3) persistence classification at gene/pathway/TF level using classes `persistent / resolving_early / established_late / P1_only / transient_4W / emergent_3mo / P1_3mo_only / persistent_directionswap` — same-direction sign required for "persistent"; (4) effect-size trajectories P1→4W→3mo for top features; (5) Early-vs-Late overlap per age (hypergeometric + full-list Spearman ρ); (6) cross-arm core signature (features persistent in BOTH arms, same direction — the paper-quality table).

21. **`--tf` flag is required in production.** It gates: 8f view 5 (TF concordance) and 8g view 3 (TF persistence). Without it both silently skip and can't be recovered without re-running 8c. Always pass `--tf` to 8c, or set `pathways.run_tf_activity: true` in YAML.

22. **Dev limitation acknowledged: 8g cannot be meaningfully exercised on dev.** `dev_split.yaml` subsets to 4W only (one M sample per group), by design, to keep smoke tests fast. With one age, every classification is `transient_4W` or `none` — code paths work, but the biology is meaningless. To exercise 8g properly needs ≥2 ages in input data; this lives on the workstation.

23. **Phase 1 = SoupX via R subprocess (locked 2026-06-10).** CellBender abandoned 2026-06-05 (pickle bug across all torch+pyro+numpy combos). SoupX is the active replacement: CRAN package, R-Bioconductor stack only, no torch dependencies. Per-sample workflow: cellranger filtered + raw counts → `SoupChannel` → `scran::quickCluster` → `setClusters` → `autoEstCont` (data-driven rho per cluster) → `adjustCounts`. Manual rho fallback (`--rho 0.10`) if scran installation fails. Scripts: `02_soupx.py` (Python orchestrator, parallel via `ProcessPoolExecutor`) + `run_soupx.R` (per-sample worker). Output: `results/{tissue}/h5ad/02_soupx_corrected/{sample_id}.h5ad`. Need to wire `02_qc.py` with prefer-soupx fallback; full re-run from Phase 2 needed once SoupX completes (corrected counts change HVG selection, scVI training, clustering, annotation).

24. **Brain CellTypist models retrained on GPU via cuML (locked 2026-06-10).** Three per-age models built from the ABC adult atlas (WMB-10Xv3): `refs/celltypist_brain_adult_class.pkl` (34 labels), `_subclass.pkl` (334 labels), `_region.pkl` (12 anatomical_division_label categories). Use `use_GPU=True` in `celltypist.train()` for ~9× speedup on class, >40× on subclass (CPU L-BFGS hangs on 334-class multinomial; not parallelizable across classes). Total retrain ~60 min on RTX 4500 Ada. cuML installed via `cuml-cu12 cudf-cu12` from `https://pypi.nvidia.com`; numba 0.65→0.64, pyarrow 24→23, cuda-toolkit 13→12.9 (CUDA 13 driver is backward-compatible with 12.9 runtime). Output classifiers save as `sklearn.linear_model._logistic.LogisticRegression` — pkls are interchangeable with CPU-trained versions for inference.

25. **CellTypist sklearn-1.7+ patch (locked 2026-06-10).** CellTypist's `train.py` hardcodes `multi_class='ovr'` in `LogisticRegression(...)` (lines 126, 146), which sklearn 1.7+ removed. Two-pattern sed fixes it; the patch lives in `.venv/` and `uv sync` reverts it. Re-apply via `./setup-remote.sh` Step 2.5 (idempotent). Removing the argument also improves calibration: sklearn 1.7+ defaults to true multinomial softmax (was one-vs-rest), giving non-bimodal `class_conf` distributions (median ~0.999 post-fix vs 0.000 with old SGD-trained pkls).

26. **Brain Phase 7 = three-tier CellTypist (locked 2026-06-10).** `class` (canonical, used by 8b/8c) = per-(Leiden cluster × age) MAJORITY vote of per-cell predictions, P1 uses `Developing_Mouse_Brain.pkl` (Di Bella 2021, built-in), 4W/3mo uses `refs/celltypist_brain_adult_class.pkl`. `subclass` = per-cell raw labels for 4W/3mo only (consumed by 7b/7d at subcluster level), P1 → sentinel `"no_subclass_model"`. `region` = per-cell raw labels for 4W/3mo only (consumed by Phase 9 for region-matched cross-species comparison), P1 → sentinel `"no_region_model"`. YAML schema: `annotation.celltypist_models.<age>.{class, subclass, region}` (nested per age).

27. **Brain marker gate + age-composition sanity (added 2026-06-10).** `BRAIN_GATE_CONFIG` in `07_annotation.py` — STRICT canonical-marker gates for borderline CellTypist calls: microglia ≥2 of {Cx3cr1, P2ry12, Tmem119, Csf1r, Aif1} → else demote to `unassigned_immune`; astrocyte ≥2 of {Aqp4, Gja1, Slc1a3, Aldh1l1} → `unassigned_glia`; OL lineage ≥1 of {Mbp, Mog, Plp1, Mag} → `unassigned_glia`; endothelial ≥2 of {Cldn5, Pecam1, Cdh5} → `unassigned_vascular`. Marker presence threshold = ≥20% of cells in the (cluster, age) group with lognorm > 0. Audit CSV gains `markers_checked`, `markers_present`, `gate_outcome` (`no_gate` | `passed` | `demoted`), `gate_label` columns. First production run: 4 demotions (496 cells, all "30 Astro-Epen" with only Aldh1l1 expressed). Separate `07_annotation_age_composition_sanity.csv` flags developmentally-implausible (cluster × age) rows (RG/IPC/neuroblast/glioblast/erythrocyte outside P1) — informational only, does not modify labels.

28. **P1 erythrocyte caveat (documented 2026-06-10).** First Phase 7 brain run produced ~81K P1 cells (~45% of P1) labeled "Blood: Erythrocyte" / "Blood: Erythroid progenitor" across clusters 0, 19, 20, 45. Two contributors: (a) Di Bella 2021's labels cover nucleated erythroblasts (basophilic / polychromatic / orthochromatic stages), which still exist in P1 mouse brain vasculature — a small number of these calls is biologically real; (b) the bulk inflation is ambient hemoglobin contamination from lysed erythroblasts during nuclei prep, with CellTypist's LogReg pattern-matching contaminated nuclei to erythroid lineage. CellBender would have prevented this; QC `pct_hemo ≤ 5%` hard cap catches the worst but not all. SoupX (Phase 1) is the engineered fix; an erythrocyte gate (`Hbb-bs / Hbb-bt / Hba-a1 / Hba-a2 / Alas2` ≥ 3 at ≥30% cells) is deferred pending SoupX results.

29. **renv Suggests workaround (locked 2026-06-10).** renv installs `Suggests` dependencies by default — for SoupX this drags Seurat → shiny → bslib → fs → libuv-dev (a system library we don't have). Project-level fix: `renv::settings$package.dependency.fields(c("Depends", "Imports", "LinkingTo"), persist = TRUE)` — writes to `renv/settings.dcf`, all future installs skip Suggests in this project. Belt-and-suspenders: set `USE_BUNDLED_LIBUV=1` in the env before any `Rscript` invocation so even if a stray Suggests sneaks in, `fs` builds without the system library.

30. **SoupX rho ranges — full production run completed 2026-06-10.** All 34 brain + 23 placenta samples returned `status=ok`, no failures, no `rho > 0.30` outliers. Per-tissue ranges: brain 3mo 0.026-0.054 (median ~0.034), brain 4W 0.014-0.039 (median ~0.023), brain P1 0.017-0.102 (median ~0.073) — **L2_P1 outlier at rho=0.017** (every other P1 sample is 0.053-0.102; smallest n_cells=11,023; could be cleaner dissection or autoEstCont underestimating on a small cluster structure; not blocking, re-run with `--rho 0.075 --sample-ids L2_P1` if needed). Placenta E12.5 0.025-0.111 (median ~0.046), E18.5 0.034-0.085 with **LCP\* (E18.5 Relaxed) clustering at 0.075-0.085** — see §2 "Pool 4 ambient signature" caveat. Total ~785K brain + ~549K placenta pre-QC nuclei. Summary CSVs: `results/{tissue}/tables/02_soupx/02_soupx_summary.csv`.

---

## 2. Critical Considerations & Caveats

### Statistical Power Limitations

**n=2 per sex per condition per age is a real limitation:**
- Cannot reliably test sex × condition interactions
- Pseudobulk DE with n=2 vs n=2 detects only large effect sizes (logFC > 1.5–2)
- Compositional analysis (scCODA) has weak variance estimation
- Single-cell-level DE treating cells as replicates is **incorrect** — reviewers will flag this

**Decision: Sex-agnostic analysis with sex as a covariate**
- Pool sexes per group (effective n=4 vs n=4 per age per condition pair)
- DESeq2 design: `~ sex + group` (sex as nuisance variable; group is the 3-level factor Early/Late/Relaxed)
- Per-sex stratified analyses as secondary/exploratory only
- Sex × group interactions reported as exploratory, underpowered

**No dam ID — litter random effect cannot be modeled:**
- Prenatal stress was applied to **dams**, not pups; pups from the same litter share gestational environment
- Two pups from one stressed dam are **not fully independent observations** for testing the stress effect
- Without dam IDs we cannot fit `~ sex + group + (1|dam)` and must treat each pup as an independent biological replicate
- This is anti-conservative: p-values are slightly optimistic, effect sizes are accurate
- **Mitigation:** include `pool` (= harvest+library batch) as a covariate where it captures some of the dam structure (pups harvested same day often from same/few dams), use FDR < 0.05 (not 0.1), require reasonable effect sizes for top hits
- **Explicit caveat in methods section:** "Dam identity was not recorded; pseudobulk DE treats each pup as an independent biological replicate, which may inflate statistical significance for traits that aggregate at the litter level."

### 10x Flex-Specific Considerations

1. **Probe-based capture (not poly-A)** — biased toward exonic sequences
2. **RNA velocity is complicated** — may not have enough unspliced reads
   - Pipeline auto-checks spliced/unspliced ratio; falls back to PAGA + pseudotime if Flex data insufficient
3. **Ambient RNA from probe leakage** — CellBender still needed
4. **Doublet rates** comparable to standard 3' chemistry

### snRNA-seq-Specific Considerations

1. **High intronic read fraction** is normal (nuclear RNA includes nascent transcripts)
2. **Mitochondrial % should be near zero** in good nuclei prep — high %mt indicates cytoplasmic contamination
3. **Ambient RNA is severe** — nuclei prep releases cytoplasmic RNA into the soup
   - Placenta especially: hemoglobin contamination (Hbb, Hba)
   - Brain: Malat1, mitochondrial genes
4. **Cell type annotation references** built from scRNA-seq may not transfer perfectly

### RNA Velocity & Cell-Cell Communication on snRNA-seq

Both are **defensible and published** but require careful framing:

**RNA Velocity:**
- Validated against matched scRNA-seq (correlations 0.94–0.99 in Alzheimer microglia data)
- Use **veloVI** or scVelo dynamical model (not steady-state)
- Caveat: velocity magnitudes noisier on nuclei due to nuclear export variation
- **For this study:** Use at P1 only (proliferating populations). Skip for 4W/3-month (steady-state).
- Cross-check with PAGA + pseudotime

**Cell-Cell Communication:**
- Published precedent in liver, brain, kidney snRNA-seq
- Use **LIANA+** as wrapper (runs CellChat, CellPhoneDB, NATMI, etc. for consensus)
- Particularly well-motivated for **placenta** (signaling-heavy tissue)
- Acknowledge in methods: inferring signaling potential from nuclear transcriptional state

### Pool/Library Confounding — REAL AND CONSEQUENTIAL

The actual pool structure has several confounds that constrain what's analyzable.

**Pool × age × group matrix (brain):**

| Age | Group | Pool 1 | Pool 2 | Pool 3 |
|---|---|---|---|---|
| 3-month | Early Stress | 4 | 0 | 0 |
| 3-month | Late Stress | 4 | 0 | 0 |
| 3-month | Relaxed | 4 | 0 | 0 |
| 4-week | Early Stress | 2 | 2 | 0 |
| 4-week | Late Stress | 1 | 3 | 0 |
| 4-week | Relaxed | 1 | 3 | 0 |
| P1 | Early Stress | 0 | 4 | 0 |
| P1 | Late Stress | 0 | 0 | 2 |
| P1 | Relaxed | 0 | 4 | 0 |

**Pool × age × group matrix (placenta):**

| Age | Group | Pool 3 | Pool 4 |
|---|---|---|---|
| E12.5 | Early Stress | 9 | 0 |
| E12.5 | Relaxed | 5 | 2 |
| E18.5 | Late Stress | 0 | 4 |
| E18.5 | Relaxed | 0 | 4 |

**What's analyzable cleanly:**
- 4-week brain: Early vs Late vs Relaxed (all 3 groups in both Pool 1 and Pool 2 — best-balanced age)
- 3-month brain: Early vs Late vs Relaxed (single pool, no pool confound to worry about within this age)
- P1 brain: Early vs Relaxed (both in Pool 2)
- E18.5 placenta: Late Stress vs Relaxed (all in Pool 4, no within-age confound)

**What's confounded and must be flagged in any analysis:**
- **P1 brain Late Stress vs anything:** the only 2 P1 Late Stress samples are in Pool 3, separately from the P1 Early/Relaxed in Pool 2. Cannot disentangle Pool 3 batch effect from Late Stress effect at P1.
- **Brain age trajectories (P1 → 4W → 3-month):** each age is dominated by a different pool (3-month entirely Pool 1; 4-week split Pool 1/2; P1 split Pool 2/3). scVI integration can partially correct, but interpretation requires caution.
- **Pool 3 mixes tissues:** has 2 brain samples and 14 placenta samples in the same library. Ambient RNA contamination across tissues is possible. The 2 P1 Late Stress brain samples may carry placental gene signatures.
- **Placenta cross-age (E12.5 vs E18.5):** each age = one pool + different conditions. **Cannot do meaningful E12.5 vs E18.5 comparison.**
- **E12.5 placenta Early vs Relaxed:** mostly within Pool 3, with 2 Relaxed in Pool 4. Acceptable but the 2 Pool-4 Relaxed samples will look like a small batch effect.

**Mitigations:**
- Always include `pool` as a covariate in DESeq2 designs (`~ sex + pool + group`)
- For scVI integration: `batch_key=pool` is correct
- Tag underpowered/confounded contrasts as `confounded_with_pool` in output tables; downweight in interpretation
- For the §8g persistence analysis, note that any age-specific signal at P1 vs 4W vs 3-month is partially confounded with pool

**Pool 4 ambient signature differs from other pools (SoupX 2026-06-10):**
- All four E18.5 Relaxed (LCP*) samples cluster at rho 0.075-0.085, distinctly higher than E12.5 Relaxed (CES*, range 0.037-0.103, median 0.047) and E18.5 Late (LSP*, range 0.034-0.067 within the same Pool 4).
- Both LCP* and LSP* live in Pool 4, but LSP* is more spread → suggests a real prep/age effect on top of any pool effect, not a pure batch artifact.
- Already controlled by `~ sex + pool + group` in the DE design; no extra mitigation needed.
- Recorded here so reviewers (and future-us) see we noticed it.

---

## 3. Compute Environment

**Available machine:** 258 GB RAM, 56 CPU cores, **1× NVIDIA RTX 4500 Ada Generation (24 GB VRAM)**, GPU on the same box as the CPU compute. Likely a shared lab workstation (Xorg + gnome-shell using ~290 MB VRAM; other Python processes occasionally present).

**Implications:**
- CellBender on GPU: ~1–2 hrs per sample (vs. 12–24 hrs CPU)
- scVI on GPU with BF16 mixed precision: ~2–3 hrs per integrated tissue object (vs. 24–48 hrs CPU)
- CellBender GPU-parallel: 2 samples simultaneously (start conservative; bump to 3 if VRAM peak stays under 20 GB)
- Total wall time per tissue: ~1–1.5 days unattended; both tissues + cross-tissue analysis in ~3–4 days of compute spread across a week of calendar time
- No need for external GPU rental; on-box GPU handles the full pipeline

**System RAM budget (per-tissue, ~720K input cells × ~20K genes, ~550–600K post-QC):**
- Single Flex h5 sample loaded sparse: 0.5–1.5 GB
- CellBender peak: 8–15 GB per sample (2 parallel on GPU → ~30 GB system RAM)
- Concatenated post-QC sparse object: 30–55 GB resident
- scVI training peak (object + model + optimizer + pinned-memory workers): 60–90 GB
- Downstream (DE, composition, pathway, LIANA+): <30 GB
- **Worst-case peak: ~90 GB during scVI** → ~2.5× headroom on 258 GB

**GPU VRAM budget (24 GB total, ~23.5 GB usable after display):**
- CellBender per sample: 4–8 GB (2 samples parallel fits comfortably; 3 fits if individual peaks stay <8 GB)
- scVI training at `batch_size=1024` with BF16: 8–12 GB
- scANVI: similar to scVI
- **Do NOT use rapids-singlecell** at this cell count — neighbors/UMAP on 600K cells needs 15–22 GB VRAM with no margin and locks the GPU. CPU runs the same operations in 10–30 min on scVI's 30-dim latent, which is fine.

**Ada-specific optimizations:**
- **BF16 mixed precision** for scVI (`precision: "bf16-mixed"` in trainer kwargs) — Ada has strong BF16 throughput; ~1.3–1.5× speedup and ~halved activation memory, no accuracy loss. Requires scvi-tools ≥1.0.
- **Sustained boost clocks** — the 4500 Ada's 210 W TDP runs cooler than a 4090 (450 W) and won't thermally throttle on multi-day jobs. Good fit for CellBender's long sequential runs.

**Key engineering principles:**
- Sparse matrices throughout — one accidental `.toarray()` or `sc.pp.scale()` on the full matrix can densify to 48+ GB instantly. Check `issparse(adata.X)` before any operation that might densify.
- Don't carry redundant layers (`counts` + `X` + `.raw` + `scaled` = 4× memory); drop `.raw` after extraction; scale inline on HVG subset only (scVI doesn't need scaled input).
- Process samples sequentially through Phase 1–3, write `.h5ad` checkpoint, free RAM, then concatenate.
- Process tissues sequentially (brain end-to-end, then placenta); only re-load both together for cross-tissue analysis (§8f), where you're working with pseudobulk summaries.
- Save checkpoints after every phase — scVI's PyTorch state doesn't always release cleanly even after `del model; torch.cuda.empty_cache()`. Kernel restart + reload between major phases is safest.
- **Serialize GPU phases** — don't run CellBender and scVI on the GPU at the same time.
- Wrap long jobs in `tmux`; set `CUDA_VISIBLE_DEVICES=0` explicitly; check `nvidia-smi` at the top of `run.py` to refuse launch if non-display VRAM is already in use.

---

## 4. Analysis Strategy

### Ecosystem: Python (Scanpy + scvi-tools)

**Rationale:**
- scVI/scANVI superior for complex batch structure (multiple libraries × ages × conditions × tissues)
- Better Python tooling for snRNA-seq quirks (CellBender, scDblFinder via rpy2)
- Pseudobulk DE via PyDESeq2 / decoupler
- Scales better for ~1.5M nuclei
- Can bridge to R for specific steps (DESeq2, CellChat) via rpy2 if needed

### Joint Integration, Per-Age Analysis

**The key analytical structure:**

| Step | Scope | Rationale |
|---|---|---|
| Integration (scVI) | All ages together | One consistent latent space |
| Clustering + annotation | Joint object | Canonical cell type labels |
| UMAP | Computed once on full data | Stable coordinates across all figures |
| Composition analysis | Per age | Cell type proportions are age-specific |
| Differential expression | Per age | `~ sex + pool + group` within each age (3-group factor) |
| Trajectory analysis | Cross-age | Maturation lineages (oligo, microglia, astrocyte) — caveat: pool-age confounded |
| Velocity | P1 only (brain) | Proliferating populations only |
| Cell-cell communication | Per age, per group | Differential signaling |

**Per-age UMAPs use the same coordinates** (just subset the integrated object). No recomputing.

**Note on the 3-group factor:** group ∈ {Early Stress, Late Stress, Relaxed}. Primary contrasts are Early-vs-Relaxed and Late-vs-Relaxed (Relaxed as reference). Early-vs-Late is a secondary contrast (compares the two stress timings to each other).

### Reference Atlases (per age)

**Brain:**
- **P1:** Di Bella et al. 2021 (Nature) — developing mouse cortex; Rosenberg et al. 2018 — P2/P11 backup
- **4W & 3-month:** Allen Brain Cell Atlas (Yao et al. 2023) — definitive adult mouse reference
- Cross-validate with CellTypist running multiple references

**Placenta:**
- Marsh & Blelloch 2020 (primary)
- Han et al. for E12-stage references
- Additional published mouse placenta atlases to be identified during execution

### Pre-Registration

**Before looking at results, document:**
- Focal cell types (a priori based on stress literature: microglia, oligo lineage, exc/inh neurons)
- Primary vs. secondary contrasts
- Significance thresholds
- Validation criteria

Protects against fishing accusations given the n=2 constraint.

---

## 5. Detailed Analysis Pipeline

### Phase 0: Validation (Mandatory First Run)

**Modules:** `src/metadata.py` + `src/sex_check.py` + `src/gene_sets.py`

This phase produces **no analysis output, only validation**. It must pass before any compute-heavy phase runs. Catches 80% of problems in 5 minutes that would otherwise surface 12 hours into CellBender.

```bash
python run.py --config config/brain.yaml --step validate
```

1. **Sample manifest validation** (`src/metadata.py`):
   - All h5 paths exist and are readable
   - All required metadata fields present per sample (no missing condition, sex, age, library, donor_id)
   - No duplicate sample IDs
   - Lock sample order alphabetically by ID (reproducibility — concatenation order affects downstream `obs_names` suffixes)
   - **Print categorical balance matrix:** n per group × age × sex — warn if any cell is <2
   - **Print library × condition × sex contingency table** — this is the confound check from §2
2. **Per-sample fingerprint** (`src/io.py`):
   - For each sample: n_cells, median UMI, median genes, top 20 genes, %mt, %hemo, %ribo
   - Single CSV: `tables/sample_fingerprints.csv`
   - Flag outliers (>3 MADs from cohort median on any metric)
3. **Y-chromosome inferred sex check** (`src/sex_check.py`):
   - Score Y-linked markers (Ddx3y, Uty, Eif2s3y, Kdm5d) and Xist per cell
   - Aggregate to sample-level inferred sex
   - Compare to metadata-declared sex; write mismatch report
   - Add `sex_check_passed` boolean to per-sample manifest (gates downstream filtering)
4. **Gene-set sanity check** (`src/gene_sets.py`):
   - Verify mouse Mt-, Rps/Rpl, Hbb/Hba/Alas2 patterns resolve to >0 genes in at least one h5
   - Pre-cache canonical gene lists (mt, ribo, hemo, sex-linked, stress-relevant) for reuse across phases
5. **Empty-droplet check:**
   - Confirm raw matrices (not just filtered) are loadable and contain ≥50K below-knee droplets per sample for CellBender to learn ambient profile

**Output:** Single HTML report `results/{tissue}/00_validation_report.html` summarizing all the above. Mandatory review before launching Phase 1.

**Auto-saved plots:**
- `manifest_balance_matrix.png` (group × age × sex heatmap)
- `library_confound_check.png` (library × condition × sex contingency)
- `sex_check_scatter.png` (Y-score vs. Xist-score per sample, colored by declared sex)
- `sample_fingerprints_heatmap.png` (QC metric outliers across cohort)

### Phase 1: Load + Metadata Attachment + Ambient RNA Correction

**Module:** `src/io.py` + `src/ambient.py` + `src/metadata.py`

1. Load 10x Flex `.h5` per sample (raw + filtered matrices)
2. **Attach metadata to `.obs`** (hard-validated): `sample_id`, `donor_id`, `condition`, `sex` (declared), `sex_inferred` (from Phase 0), `age` or `stage`, `library`, `batch_run_date` (if available)
3. **Gene-name harmonization:** Ensembl IDs as primary key (`adata.var_names`), gene symbols in `adata.var['symbol']`. Resolves the `Hbb-b1` vs `Hbb-bs` problem and avoids mid-pipeline ambiguity
4. **Preserve empty droplets:** keep ≥50K below-knee droplets in raw matrix for CellBender
5. Run **CellBender** (`remove-background`) on raw matrix
   - GPU mode, `--epochs 150`, `--cells-posterior-reg 50`
   - GPU-parallel 2 samples (bump to 3 if VRAM peak <20 GB)
   - Re-does cell vs. empty droplet call (better than Cell Ranger defaults)
6. Save: `results/{tissue}/h5ad/02_ambient_corrected/{sample_id}.h5ad`

**Auto-saved plots:**
- `{sample_id}_barcode_rank.png` (original Cell Ranger knee)
- `{sample_id}_cellbender_posterior.png` (CellBender cell probability)
- `{sample_id}_ambient_genes.png` (top contaminants: Hbb in placenta, Malat1 in brain)
- `summary_ambient_fraction.png` (across samples — flags outliers)

### Phase 2: Per-Sample QC

**Module:** `src/qc.py`

1. Compute metrics: `n_genes`, `total_counts`, `pct_counts_mt`, `pct_counts_ribo`, `pct_counts_hemo`
2. **Automated MAD-based thresholding** (5 MADs default, per sample)
3. Apply per-sample manual overrides from `sample_overrides.yaml` if present
4. Flag sex chromosome genes (Xist, Ddx3y, Uty, Eif2s3y, Kdm5d) — useful as sex-mixup QC

**snRNA-specific defaults:**
```yaml
qc:
  pct_mt_max: 1.0           # nuclei should be near-zero
  pct_hemo_max: 5.0         # critical for placenta
  min_genes: auto           # MAD-based
  n_mads: 5
  pct_counts_in_top_20_genes_max: 50
```

**Auto-saved plots:**
- `{sample_id}_violin_pre/post.png`
- `{sample_id}_scatter_counts_vs_genes.png` (with threshold lines)
- `{sample_id}_histogram_thresholds.png` (auto-cutoffs marked)
- `summary_qc_table.csv` + `summary_cells_per_sample.png`

### Phase 3: Doublet Detection

**Module:** `src/doublets.py`

1. Run **scDblFinder** per library (before merging — doublets form within a capture)
2. Optionally cross-check with **scrublet**
3. Remove doublets

**Auto-saved plots:**
- `{sample_id}_doublet_score_dist.png`
- `{sample_id}_doublets_on_umap.png` (preliminary per-sample UMAP)
- `summary_doublet_rate_per_sample.png`

### Phase 4: Concatenation + HVG Selection

**Module:** `src/integration.py` (part 1)

1. Concatenate all samples for the tissue
2. Log-normalize
3. Identify HVGs (2000–4000 for brain, 2000 for placenta)
4. **Flag for HVG exclusion:** hemoglobin, sex-linked genes, pregnancy genes (placenta) — prevent dominance of integration

**Auto-saved plots:**
- `hvg_dispersion.png`
- `pre_integration_umap_by_library.png` (shows batch effect)
- `pre_integration_umap_by_condition.png` (diagnostic for library/condition confounding)
- `pre_integration_umap_by_age.png`

### Phase 5: Integration with scVI/scANVI

**Module:** `src/integration.py` (part 2)

1. Setup scVI:
   - `batch_key=library`
   - `categorical_covariate_keys=[age, condition, sex]` (preserve as biology)
   - `continuous_covariate_keys=[pct_counts_mt]`
2. Train scVI (CPU mode, `max_epochs=200`, `batch_size=512`, early stopping)
3. If reference atlas provided → use scANVI for label-aware integration
4. Save trained model + integrated AnnData

**Auto-saved plots:**
- `scvi_loss_curve.png`
- `post_integration_umap_by_library.png` (should mix well now)
- `post_integration_umap_by_condition.png`
- `post_integration_umap_by_age.png`
- `integration_metrics_scib.png` (batch correction vs. bio conservation)

### Phase 6: Clustering

**Module:** `src/clustering.py`

1. Build neighbor graph on scVI latent space
2. **Multi-resolution Leiden clustering** (0.2, 0.4, 0.6, 0.8, 1.0, 1.5, 2.0)
3. Compute **clustree** diagram (clusters splitting across resolutions)
4. Default resolution via silhouette/modularity; YAML-overridable

**Auto-saved plots:**
- `clustree.png` (essential for choosing resolution)
- `umap_leiden_res{X}.png` (one per resolution)
- `cluster_qc_metrics.png` (per-cluster mean QC — catches junk clusters)
- `cluster_composition_by_sample.png` (catches single-sample clusters = batch artifacts)

### Phase 7: Annotation

**Module:** `src/annotation.py`

**Two-track approach:**

1. **Reference-based:** CellTypist or scANVI label transfer
   - Brain P1: Di Bella 2021 + Rosenberg 2018
   - Brain 4W/3-month: Allen Brain Cell Atlas
   - Placenta: Marsh & Blelloch 2020 + others
2. **Marker-based:** `rank_genes_groups` + curated marker dotplots from YAML

Final annotations = manual reconciliation (notebook step).

**Auto-saved plots:**
- `umap_celltypist_predictions.png` + confidence scores
- `marker_dotplot_curated.png`
- `marker_heatmap_top10_per_cluster.png`
- `umap_final_annotations.png`
- `annotation_confusion_matrix.png` (reference vs. marker agreement)

### Phase 8: Downstream Biology

**Architectural principle:** All contrasts (composition, DE, pathway, communication) are defined **declaratively in YAML**, not hard-coded. The downstream engine iterates over the contrast specification, runs the same pseudobulk/scCODA/LIANA+ machinery for each, and writes uniformly-structured output tables. Adding a new contrast = edit YAML, not Python.

**Contrast families (see §6 config block for full spec):**

| Family | Example | Status |
|---|---|---|
| Within-age, Early vs Relaxed | 4W brain Early vs Relaxed | ✓ Primary |
| Within-age, Late vs Relaxed | 4W brain Late vs Relaxed | ✓ Primary (⚠ P1: pool-confounded) |
| Within-age, Early vs Late | 4W brain Early vs Late | ✓ Secondary |
| 3-group omnibus per age | 4W brain Early vs Late vs Relaxed (F-test) | ✓ Primary |
| Within-group, across ages | Early Stress P1 vs 4W vs 3-month | ⚠ Pool-confounded |
| Group × age interaction | (P1 E − P1 R) vs. (4W E − 4W R) | ⚠ Underpowered + partly confounded |
| Persistent vs. resolving | DE in P1 ∩ 4W ∩ 3-month vs. P1-only | Set operation; confound caveat applies |
| Within-age, sex-stratified | 4W male Early vs Relaxed | ⚠ Exploratory (n=2 per cell) |
| Sex × group × age interaction | Full three-way | ⚠ Underpowered, reported only |

All exploratory and confounded contrasts produce output tagged `flag: underpowered_exploratory` or `flag: confounded_with_pool` in tables and figures.

**Reference level for `group`:** Relaxed (so positive logFC = upregulated in stress).

#### 8a. Composition (`src/composition.py`)
- **scCODA** (Bayesian, accounts for sum-to-1 constraint)
- **propeller** cross-check (frequentist)
- Iterates over all contrast families above (animal as statistical unit)
- Plots: stacked bars, boxplots, scCODA forest plots, heatmaps per contrast
- Master table: `composition_results.csv` with `[contrast_name, contrast_family, celltype, log2FC, FDR, flag]`

#### 8b. Differential Expression (`src/de.py`)
**Critical:** Pseudobulk only, never single-cell-level. Animal (pup, via `donor_id`) is the statistical unit. Dam-level random effects cannot be modeled (no dam ID — see §2).

1. Sum counts per `donor_id` per cell type → pseudobulk matrix
2. Filter: ≥10 cells in ≥3 samples per group
3. **PyDESeq2** (primary) / edgeR via rpy2 (cross-check)
4. **Iterate over all contrast families** from the declarative YAML spec
5. Master combined table: `de_results.csv` with `[contrast_name, contrast_family, celltype, gene, logFC, padj, direction, flag]`

**Plots per contrast:** volcano per cell type, MA plots, DEG count heatmap, top DEG dotplot, pseudobulk PCA (sanity check — should separate by biological factor, not library).

#### 8c. Pathway Analysis (`src/pathways.py`)
- GSEA on ranked DE statistics (decoupler + MSigDB Hallmark, Reactome, GO BP)
- **Stress-relevant gene sets:** GR target genes, HPA axis, neuroinflammation, synaptic, mitochondrial, ER-stress, oxidative phosphorylation
- TF activity inference (CollecTRI)
- Per-cell pathway scoring (decoupler ulm/mlm)
- Runs against every contrast in §8b
- Master table: `pathway_results.csv` with `[contrast_name, contrast_family, celltype, pathway, NES, FDR, flag]`

#### 8d. Trajectory (`src/trajectory.py`) — Brain primarily
- **PAGA** (robust, reviewer-friendly) — always run
- **Diffusion pseudotime** anchored at progenitor clusters
- **veloVI / scVelo** — P1 only, with auto-diagnostic for Flex feasibility
- **CellRank 2** — fate probabilities (works without velocity)
- **Cross-group comparison:** pseudotime distributions, PAGA structure, developmental "speed" — compare each stress group (Early, Late) against Relaxed
- **Cross-age trajectories:** P1 → 4W → 3-month maturation, compared between groups (pool-confounded — see §2)
- **Focal lineages:** oligodendrocyte (OPC → mature OL), microglia developmental states, astrocyte maturation
- Sex-stratified trajectories as secondary (exploratory)

#### 8e. Cell-Cell Communication (`src/communication.py`)
- **LIANA+** (consensus across CellChat, CellPhoneDB, NATMI, etc.)
- Iterates over all contrast families (within-age, across-age, sex-stratified)
- **Differential communication** computed per contrast
- **Placenta focus:** trophoblast ↔ decidua ↔ fetal vasculature
- **Brain focus:** neuron ↔ glia, microglia ↔ neuron (neuroinflammation)
- Master table: `communication_results.csv` with `[contrast_name, contrast_family, source_celltype, target_celltype, ligand, receptor, score_delta, FDR, flag]`

#### 8f. Cross-Tissue Link (`src/cross_tissue.py`) — THE UNIQUE ANGLE
- DEG overlap: placenta vs. brain at biologically-aligned timepoints
- Pathway concordance heatmaps
- **RRHO** (rank-rank hypergeometric overlap) — standard for transcriptional signature comparison
- Placental-derived signals with known brain receptors
- **Biologically-aligned cross-tissue mappings (given actual sample availability):**
  - **Early Stress arm:** E12.5 placenta (Early-vs-Relaxed) → P1 brain (Early-vs-Relaxed) → 4W brain (Early-vs-Relaxed) → 3-month brain (Early-vs-Relaxed)
  - **Late Stress arm:** E18.5 placenta (Late-vs-Relaxed) → P1 brain (Late-vs-Relaxed) → 4W brain (Late-vs-Relaxed) → 3-month brain (Late-vs-Relaxed)
  - The placenta timing aligns with the stress exposure window (Early stress → mid-gestation placenta; Late stress → late-gestation placenta), giving two parallel temporal cascades to compare
- **Caveat:** the P1 Late Stress brain samples are pool-confounded (see §2); the Late Stress cascade should be interpreted with this in mind
- Iterates across all contrast families (so the cross-tissue link is also computed for sex-stratified contrasts as exploratory)

#### 8g. Cross-Age & Persistence Analysis (`src/cross_age.py`)
**The "developmental trajectory of the stress signature" view.** Operates on §8b–c output tables (no re-running of DE). Run separately for Early-vs-Relaxed and Late-vs-Relaxed contrasts.

1. **Persistence sets** (gene-level), separately for Early-vs-Relaxed and Late-vs-Relaxed:
   - Persistent: DE at P1 AND 4W AND 3-month (same direction)
   - Early-only: DE at P1 only
   - Emergent: DE at 3-month only
   - Transient: DE at 4W only
   - **Caveat:** P1 Late Stress is pool-confounded (see §2), so Late-vs-Relaxed persistence calls involving P1 carry a `confounded_with_pool` flag
2. **Pathway persistence:** same logic, but at the pathway level (more robust given n=2)
3. **Trajectory of effect size:** for each celltype × gene, plot logFC across ages — visualize whether stress effect attenuates, amplifies, or reverses
4. **Cross-age compositional drift:** does the stress group's cell type composition diverge from Relaxed over development?
5. **Within-group developmental contrasts:** P1 vs 4W and 4W vs 3-month *within* the Early Stress group, *within* the Late Stress group, *within* the Relaxed group. Compare: "does stress group development differ from Relaxed development?"
6. **Sex-stratified persistence (exploratory):** same as #1 but per sex
7. **Early vs Late comparison:** at each age, do the Early-Stress and Late-Stress signatures overlap or are they distinct? Hypergeometric test on the DE gene sets.

**Plots:**
- `persistence_venn_early.png` (P1 ∩ 4W ∩ 3-month Early-vs-Relaxed DEGs per celltype)
- `persistence_venn_late.png` (same for Late-vs-Relaxed)
- `effect_size_trajectory.png` (logFC across ages, one line per top gene, faceted by Early/Late)
- `pathway_persistence_heatmap.png` (pathways × ages, colored by NES, faceted by celltype × contrast)
- `stress_vs_relaxed_developmental_divergence.png` (UMAP centroid drift across ages, separate for Early and Late)
- `early_vs_late_overlap_venn.png` (overlap of Early-vs-Relaxed and Late-vs-Relaxed DEGs per age)

**Master table:** `cross_age_results.csv` with `[celltype, gene_or_pathway, contrast (Early_vs_Relaxed | Late_vs_Relaxed), P1_logFC, 4W_logFC, 3mo_logFC, persistence_class, sex_specific, flag]`

### Phase 9: Reporting

**Module:** `src/figures.py` + `notebooks/05_figures.ipynb`

- Publication-ready figures
- **HTML summary report per tissue** with all QC + biological findings
- Export canonical final AnnData

### Phase 10: Reproducibility

- Run logs (timestamp, params, versions, RAM peak)
- Config snapshots saved per run
- Fixed random seeds (scVI, UMAP, Leiden)
- **`uv.lock` (Python) and `renv.lock` (R)** — both committed to git; pin every transitive dependency to exact versions and hashes (see §12)
- Git for code; data outside git
- `manifest.json` with file checksums

---

## 6. Pipeline Architecture

### Directory Structure

> **Implementation note:** the original plan in this section described an `src/snrna/` Python package. We deliberately switched to a flat `scripts/` layout (see `INSTRUCTIONS.md` and §1b). The structure below reflects what's actually in the repo. The function descriptions from the original `src/` plan are preserved as inline comments to document scope of each phase script.

```
Analysis/
├── config/
│   ├── brain.yaml              # 34 samples, generated from sample_metadata.csv
│   ├── placenta.yaml           # 23 samples, generated
│   └── dev.yaml                # hand-maintained: samples_from brain.yaml + subset.sample_ids
├── scripts/
│   ├── _utils.py               # shared helpers: load_config, add_lognorm, phase_paths, select_accelerator
│   ├── build_yaml.py           # regenerate brain.yaml + placenta.yaml from sample_metadata.csv
│   ├── 01_validate.py          # Phase 0: manifest validation, sex check, fingerprints
│   ├── 02_qc.py                # Phase 2: per-sample QC (MAD + floors + caps + cohort flag)
│   ├── 03_doublets.py          # Phase 3: scDblFinder per pool (Python driver)
│   ├── run_scdblfinder.R       # Phase 3: R subprocess called by 03_doublets.py
│   ├── 04_integration_prep.py  # Phase 4: concat + log-norm + HVG + exclusion
│   ├── 05_integration.py       # Phase 5: scVI training + UMAPs
│   ├── 06_clustering.py        # Phase 6: multi-res Leiden + clustree (TBD)
│   ├── 07_annotation.py        # Phase 7: CellTypist + markers (TBD)
│   ├── 02_ambient.py           # Phase 1: CellBender wrapper, workstation-only (TBD)
│   └── ...                     # additional phase scripts as we go
├── notebooks/
│   ├── 01_qc.ipynb             # load 03_qc_filtered/*.h5ad and inspect inline (TBD)
│   ├── 02_integration.ipynb    # load 06_integrated/all_samples.h5ad and replot (TBD)
│   └── ...                     # one per phase, thin: just load + plot
├── data/                       # raw Cell Ranger h5 files (gitignored)
│   ├── Pool1/ Pool2/ Pool3/ Pool4/
│   └── per_sample_outs/{sample_id}/sample_filtered_feature_bc_matrix.h5
├── results/                    # all outputs (gitignored)
│   └── {tissue}/               # brain/ placenta/ dev/
│       ├── h5ad/
│       │   ├── 03_qc_filtered/             # one .h5ad per sample
│       │   ├── 04_doublets_removed/        # one .h5ad per sample
│       │   ├── 05_integration_ready/all_samples.h5ad   # concatenated, post-HVG
│       │   └── 06_integrated/
│       │       ├── all_samples.h5ad        # post-scVI, with X_scVI + UMAPs
│       │       └── scvi_model/             # trained scVI model dir
│       ├── plots/
│       │   ├── 02_qc/                      # per-sample violin/scatter/threshold plots
│       │   ├── 03_doublets/                # per-pool score hist + per-sample rate bar
│       │   ├── 04_integration_prep/        # HVG dispersion, cells per sample, exclusion summary
│       │   └── 05_integration/             # pre/post UMAPs, loss curve
│       ├── tables/
│       │   ├── summary_qc.csv              # per-sample post-QC counts + cohort outlier flags
│       │   ├── summary_doublets.csv
│       │   ├── summary_integration_prep.csv
│       │   └── scvi_training_history.csv
│       └── validation/                     # Phase 0 outputs (flat dir, not under plots/)
├── sample_metadata.csv         # canonical source of truth for samples
├── run_pipeline.sh             # running log of commands (not executable)
├── pyproject.toml              # Python deps (uv reads this)
├── uv.lock                     # Python lockfile — COMMIT TO GIT
├── .python-version             # pins Python to 3.12
├── .gitignore
└── README.md
```

### Original phase-by-phase scope (still valid as a logical map)

The original `src/` module names from the planning phase map to the current flat scripts as follows. Function-level scope from the original plan is preserved:

| Original `src/` module | Current implementation | Status |
|---|---|---|
| `io.py`, `metadata.py`, `sex_check.py`, `gene_sets.py` | inlined in `01_validate.py` and `02_qc.py` | done |
| `qc.py` | `02_qc.py` | done |
| `ambient.py` | `02_ambient.py` (TBD, workstation) | not started |
| `doublets.py` | `03_doublets.py` + `run_scdblfinder.R` | done |
| `integration.py` | `04_integration_prep.py` + `05_integration.py` (split into prep + scVI training) | done |
| `clustering.py` | `06_clustering.py` (TBD) | not started |
| `annotation.py` | `07_annotation.py` (TBD) | not started |
| `contrasts.py`, `composition.py`, `de.py`, `pathways.py`, `trajectory.py`, `communication.py`, `cross_age.py`, `cross_tissue.py` | Phase 8 sub-scripts (TBD) | not started |
| `provenance.py`, `plotting.py`, `utils.py` | partially in `_utils.py`; rest deferred until needed | partial |


### YAML Configuration

**Example `config/brain.yaml`** — generated from `sample_metadata.csv` by `scripts/build_yaml.py`. Re-run that script any time the CSV changes. Full file has all 34 brain samples; structure shown below.

```yaml
tissue: brain
group_reference: Relaxed              # +logFC = upregulated in stress
results_dir: results/brain

samples:
  - id: E1
    donor_id: m_E1
    h5: data/Pool1/.../E1/sample_filtered_feature_bc_matrix.h5
    raw_h5: data/Pool1/.../E1/sample_raw_feature_bc_matrix.tar.gz  # for CellBender
    age: 3mo
    group: Early_Stress
    sex: F
    pool: Pool1
    library: Pool1
  - id: L1-P1M1
    donor_id: m_L1
    h5: data/Pool3/.../L1-P1M1/sample_filtered_feature_bc_matrix.h5
    raw_h5: data/Pool3/.../L1-P1M1/sample_raw_feature_bc_matrix.tar.gz
    age: P1
    group: Late_Stress
    sex: M
    pool: Pool3                       # ⚠ confounded — see §2 pool table
    library: Pool3
  # ... 32 more samples

qc:
  pct_mt_max: 1.0                     # snRNA: near zero
  pct_hemo_max: 5.0                   # critical for placenta
  n_mads: 5                           # MAD bounds on n_genes & log-counts
  min_counts: 500                     # hard UMI floor
  min_genes: 200                      # hard gene floor

sex_markers:
  y_linked: [Ddx3y, Uty, Eif2s3y, Kdm5d]
  x_linked: [Xist]

random_seed: 42
```

**Example `config/dev.yaml`** — hand-maintained, uses indirection:

```yaml
tissue: brain
results_dir: results/dev
samples_from: config/brain.yaml       # pull sample records from here

subset:
  enabled: true
  sample_ids: [E1-4WkM1, L1-4WkM1, S1-4WkM1]   # one per group, all 4W M, Pool1
  max_cells_per_sample: 500                     # used by 02_qc.py

qc: {pct_mt_max: 1.0, pct_hemo_max: 5.0, n_mads: 5, min_counts: 500, min_genes: 200}
sex_markers: {y_linked: [Ddx3y, Uty, Eif2s3y, Kdm5d], x_linked: [Xist]}
random_seed: 42
```

**Not yet in YAML** (planned for future phases, will be added when the corresponding scripts land):

```yaml
age_groups:
  P1: {reference_atlas: /refs/di_bella_2021.h5ad, backup_reference: /refs/rosenberg_2018.h5ad, run_velocity: true}
  4W: {reference_atlas: /refs/abc_atlas.h5ad, run_velocity: false}
  3mo: {reference_atlas: /refs/abc_atlas.h5ad, run_velocity: false}

integration:
  batch_key: pool
  categorical_covariates: [age, group, sex]
  continuous_covariates: [pct_counts_mt]
  n_hvg: 3000                         # 2000 for placenta

scvi:
  n_layers: 2
  n_latent: 30
  max_epochs: 400
  batch_size: 1024
  early_stopping_patience: 30
  # accelerator/precision auto-selected by _utils.select_accelerator()
```

### Future YAML expansion — declarative contrast spec (TBD in code)

The downstream phases (8a composition, 8b DE, 8c pathway, 8e communication, 8g cross-age) will be driven by a single declarative contrast block. This is the planned shape — not yet wired up since Phase 6 hasn't started — but the spec is preserved here verbatim so it can be dropped into the YAML when needed:

```yaml
# Declarative contrast specification — drives §8a (composition), §8b (DE),
# §8c (pathway), §8e (communication), and §8g (cross-age) downstream engines.
contrasts:

  # ------ PRIMARY ------

  early_vs_relaxed_per_age:
    description: "Early Stress vs Relaxed, within each age — primary"
    design: "~ sex + pool + group"
    group_by: age
    test: group
    levels: [Early_Stress, Relaxed]                # Early vs Relaxed
    flag: primary

  late_vs_relaxed_per_age:
    description: "Late Stress vs Relaxed, within each age — primary (⚠ P1 confounded with Pool 3)"
    design: "~ sex + pool + group"
    group_by: age
    test: group
    levels: [Late_Stress, Relaxed]
    flag: primary
    confound_warnings:
      P1: "Late Stress at P1 is in Pool 3 only; Relaxed at P1 is in Pool 2 — pool-vs-group is fully confounded at this age"

  omnibus_3group_per_age:
    description: "F-test: do the three groups differ at each age?"
    design: "~ sex + pool + group"
    group_by: age
    test: group_omnibus                            # tests if ANY group differs
    flag: primary

  # ------ SECONDARY ------

  early_vs_late_per_age:
    description: "Early Stress vs Late Stress — secondary, do the two stress timings produce the same signature?"
    design: "~ sex + pool + group"
    group_by: age
    test: group
    levels: [Early_Stress, Late_Stress]
    flag: secondary

  within_group_across_age:
    description: "Developmental trajectory within each group — pool-confounded with age, interpret cautiously"
    design: "~ sex + age"
    group_by: group
    test: age
    pairwise: [[P1, 4W], [4W, 3mo], [P1, 3mo]]
    flag: confounded_with_pool

  # ------ EXPLORATORY ------

  within_age_sex_stratified:
    description: "Stress contrasts within each sex, within each age — exploratory"
    design: "~ pool + group"
    group_by: [age, sex]
    test: group
    levels: [Early_Stress, Relaxed]
    flag: underpowered_exploratory

  group_x_age_interaction:
    description: "Does the stress effect change with age?"
    design: "~ sex + pool + group * age"
    test: "group:age"
    flag: underpowered_exploratory

  sex_x_group_x_age_interaction:
    description: "Full three-way — reported only, not interpreted"
    design: "~ sex * group * age + pool"
    test: "sex:group:age"
    flag: underpowered_exploratory

  # ------ POST-HOC SET OPERATIONS (run by src/cross_age.py) ------

  persistence_early:
    description: "Persistent / early-only / emergent / transient classification of Early-vs-Relaxed DEGs"
    source_contrast: early_vs_relaxed_per_age
    ages_required: [P1, 4W, 3mo]
    flag: derived

  persistence_late:
    description: "Same for Late-vs-Relaxed DEGs (P1 carries confound flag)"
    source_contrast: late_vs_relaxed_per_age
    ages_required: [P1, 4W, 3mo]
    flag: derived

stress_focused_cell_types:
  - microglia
  - oligodendrocyte_lineage
  - excitatory_neurons
  - inhibitory_neurons
  - astrocytes

compute:
  n_threads: 56
  use_gpu: true
  cuda_device: 0
  gpu_model: "RTX 4500 Ada"

cellbender:
  epochs: 150
  cells_posterior_reg: 50
  parallel_samples: 2
  use_cuda: true

scvi:
  accelerator: gpu
  devices: 1
  precision: "bf16-mixed"
  max_epochs: 400
  batch_size: 1024
  early_stopping: true
  early_stopping_patience: 30
  n_layers: 2
  n_latent: 30
  dataloader_kwargs:
    num_workers: 4
    pin_memory: true

scanpy:
  n_jobs: 56

qc:
  pct_mt_max: 1.0
  pct_hemo_max: 5.0
  min_genes: auto
  n_mads: 5
  pct_counts_in_top_20_genes_max: 50

random_seed: 42
```

**Example `config/placenta.yaml`:** same structure but with `tissue: placenta`, `age_groups: {E12.5, E18.5}`, only the contrasts that are analyzable (E12.5 Early-vs-Relaxed, E18.5 Late-vs-Relaxed — see §2 for what's NOT analyzable due to incomplete factorial).

### CLI Usage

```bash
# MANDATORY FIRST RUN — Phase 0 validation, no compute, ~5 min
python run.py --config config/brain.yaml --step validate

# Run one phase
python run.py --config config/brain.yaml --step qc

# Run multiple phases
python run.py --config config/brain.yaml --steps qc,integration,annotation

# Run everything (will refuse to start if --step validate hasn't passed)
python run.py --config config/brain.yaml --step all

# Resume from last valid checkpoint
python run.py --config config/brain.yaml --step downstream --resume

# Dry run (config syntax check only — no validation of data)
python run.py --config config/brain.yaml --step all --dry-run
```

**Each step:**
1. Checks if input `.h5ad` exists (resumes from last valid checkpoint)
2. Loads config
3. Calls relevant `src/` module
4. Writes output `.h5ad`, plots, tables, log
5. Updates manifest

---

## 7. Open Questions Before Code Writing

1. **Library structure (CRITICAL):**
   - How are samples distributed across 10x Flex multiplexes?
   - Are stress + control mixed within libraries, or in separate runs?
   - Are M + F mixed within libraries?
   - Total library count per tissue?

2. **Additional data to incorporate:**
   - Behavioral data on offspring (anxiety, HPA function)?
   - Corticosterone / HPA hormone measurements?
   - Archived tissue available for RNAscope/IHC validation?
   - Parallel bulk RNA-seq or proteomics?

3. **Detailed pipeline requirements** (user to provide):
   - Any tool preferences/swaps
   - Specific plot styles
   - Additional analyses not covered

---

## 8. Publication Strategy (Deferred Discussion)

**Realistic IF 12–15 targets:**
- Nature Communications, Molecular Psychiatry, Biological Psychiatry, Genome Biology
- Pure bioinformatics → likely IF 6–10 (Communications Biology, iScience)
- To push to IF 12–15: typically need RNAscope/IHC validation + behavioral data

**Strongest publication framing:**
- **Cross-tissue developmental cascade:** E12.5/E18.5 placenta → P1/4W/3-month brain transcriptional programming under prenatal stress; two parallel arms (Early Stress vs Late Stress) showing distinct trajectories
- Multi-age × multi-tissue × two-stress-windows design is genuinely unique
- Placenta-brain axis is a hot topic with few snRNA-seq papers

**To strengthen publication:**
- RNAscope/IHC validation of top 2–3 findings (~$5–10K, 2–3 months)
- Behavioral validation of stress model
- External dataset comparison (human prenatal stress, postmortem psychiatric)

---

## 9. Remote Workflow & Repo Layout

### Environment

- **Remote machine** accessed via VPN + SSH from local Mac
- **VSCode Remote-SSH** is the primary IDE (host alias: `remote-snRNA`)
- **uv + renv** for environment management (see §12 — conda channels are blocked by corporate firewall, so the standard conda/mamba workflow isn't available)
- **tmux** for long-running jobs (CellBender, scVI) so they survive SSH drops
- **HTML reports** auto-generated per phase for quick remote sanity checks (no need to download plots)

### GPU operational hygiene

The GPU is on the same workstation that may be shared (Xorg + occasional other Python jobs visible in `nvidia-smi`). For multi-day CellBender + scVI runs:

- **Pre-flight check in `run.py`:** query `nvidia-smi --query-gpu=memory.used --format=csv,noheader,nounits`; warn or refuse to launch if non-display VRAM is >2 GB
- **Monitor during runs:** `nvidia-smi -l 5` in a side tmux pane to watch VRAM and utilization
- **Fragmentation control:** `export PYTORCH_CUDA_ALLOC_CONF=max_split_size_mb:512` before launching long scVI runs
- **Explicit cleanup between GPU phases:** `del model; torch.cuda.empty_cache(); gc.collect()` — PyTorch doesn't always release cleanly
- **Long-job launch pattern:** `CUDA_VISIBLE_DEVICES=0 tmux new -d -s scvi 'python run.py --config config/brain.yaml --step integration'`
- **Coordinate with co-users** of the workstation before launching — a co-launched ClearMap or other GPU job mid-run will OOM both. Consider a simple `/tmp/gpu_in_use_by_<user>` lock file convention.
- **Clean up orphan processes** before launching: stale Python processes occasionally hold tiny VRAM allocations (a few MB each) but indicate the GPU is shared, so confirm no one else has an active job.

### Code organisation: scripts + notebooks (no duplication)

Single codebase, two invocation modes:

| File | Role | Plots |
|---|---|---|
| `src/*.py` | Real implementation — all logic lives here as functions | Saved to disk |
| `run.py` | Batch CLI runner — imports from `src/`, runs in tmux for unattended jobs | Saved to disk |
| `notebooks/*.ipynb` | Thin interactive viewers — **load saved `.h5ad` checkpoints** and inspect inline in VSCode | Inline |

**Notebook pattern (default):** load the checkpoint produced by `run.py`, don't re-run.
```python
adata = sc.read_h5ad("results/brain/h5ad/03_qc_filtered/all_samples.h5ad")
sc.pl.violin(adata, ['n_genes', 'pct_counts_mt'])
```

Notebooks only re-call `src/` functions when interactively tweaking parameters.

### Typical iteration loop

1. `tmux` → `python run.py --config config/brain.yaml --step qc` (runs unattended)
2. Reconnect via VSCode → open `notebooks/01_qc.ipynb` → inspect inline
3. If thresholds need tweaking: edit `config/brain.yaml`, rerun the step
4. Or open `results/brain/report.html` for a quick all-in-one summary
5. Move to next phase

---

## 10. Next Steps

All scripts through 8g are implemented and smoke-tested on dev. The remaining work is:

**Pre-workstation prep (manual, in any order):**
1. Refine the pathways stress gene sets in `08c_pathways.py` (currently a SCAFFOLD with UNVERIFIED placeholders) with real literature lists. The `STRESS_AXIS_GENES` curated list in `08f_cross_tissue.py` (60 genes covering GR/MR, CRH, HPA, cytokines, GABA, serotonin, BDNF, steroidogenesis) is a good starting point and overlaps with what 8c should use.
2. Build the Allen Brain Cell Atlas reference: `uv pip install abc_atlas_access` from GitHub, stage a labeled ref h5ad via `prepare_reference.py`, train a CellTypist `.pkl` on the ABC subset for adult 4W/3mo brain (no built-in adult mouse CellTypist model exists), and point `annotation.celltypist_models.{4W,3mo}` at the trained model.
3. Identify mouse placenta reference atlas (Marsh & Blelloch 2020 is the primary candidate); stage labeled h5ad similarly.

**Workstation execution (when ready):**
4. Run `02_ambient.py` (CellBender, GPU) — script TBD in `.venv-cellbender`.
5. Run brain pipeline end-to-end with `--config config/brain.yaml`; per-phase commands and required flags documented in `run_pipeline_WS.sh` (workstation-optimised runbook).
6. Run placenta pipeline end-to-end with `--config config/placenta.yaml`; skip 8g (placenta has incomplete cross-age factorial; script exits cleanly).
7. Run 8f cross-tissue once both brain and placenta finish (needs both tissues' 8b/8c output).

**Phase 9 reporting (deferred until after workstation results):**
8. Assemble publication figures from real workstation outputs; the per-phase outputs are organised under `plots/<phase>/...` and `tables/<phase>/...`, so a thin figure-assembly notebook can pull headline panels into a paper-ready layout.

**Phase 10 provenance (small, can be done before or alongside workstation runs):**
9. `manifest.json` + `provenance.py` snapshot of config + git hash + package versions + output checksums per run. Lock files (uv.lock, renv.lock) are already in place; this just adds an audit trail per execution.

**Open questions** (most can be resolved on first workstation pass):
- Behavioral data on offspring (anxiety, HPA function) — to incorporate alongside transcriptomics?
- Corticosterone / HPA hormone measurements available?
- Archived tissue for RNAscope/IHC validation of top hits?

---

## 11. Summary of Key Decisions

| Decision Point | Choice | Rationale |
|---|---|---|
| **Ecosystem** | Python (Scanpy + scvi-tools) | Better batch integration, scales to ~1.5M nuclei |
| **Compute** | RTX 4500 Ada (24 GB) + 258 GB RAM + 56 cores | On-box GPU + ample system RAM; ~1.5 days/tissue |
| **Validation gate** | Phase 0 mandatory before any compute | Catches sample swaps, missing metadata, confounds in 5 min |
| **Sex assignment** | Y-chromosome inferred + declared metadata, mismatch flagged | Detects sample swaps; also needed for E12.5 placenta (sex unknown at sampling) |
| **Donor tracking** | `donor_id` distinct from `sample_id` in `.obs`; **no dam ID available** | Pup is the statistical unit; dam random effect cannot be modeled (see §2 caveat) |
| **Group factor** | 3-level: Early_Stress / Late_Stress / Relaxed (reference) | Reflects actual experimental design (not binary stress/control) |
| **Batch_key for integration** | `pool` (Pool 1–4) | Reflects actual library/multiplexing structure |
| **Gene naming** | Ensembl IDs as `var_names`, symbols in `var['symbol']` | Avoids `Hbb-b1` vs `Hbb-bs` ambiguity |
| **Ambient RNA** | CellBender per sample (GPU, 150 epochs) | Re-does cell calling, handles snRNA-seq contamination |
| **Doublets** | scDblFinder per pool | Doublets form within capture |
| **QC thresholds** | Per-sample MAD (5 MADs) + hard floors (min_counts=500, min_genes=200) + hard caps (mt≤1%, hemo≤5%) | Per-sample is standard practice; floors catch debris that MAD misses on skewed distributions |
| **Cohort outlier flag** | Sample median UMI/genes >3 cohort-MADs below median → flag (n≥5 only) | Catches failed-prep samples whose own MAD made bounds permissive |
| **Integration** | scVI, batch_key=pool, BF16 mixed precision on GPU | State-of-the-art for complex batch structure; Ada-optimized. scANVI deferred until Phase 7 atlas labels available |
| **Lognorm layer** | Computed in Phase 4, dropped after Phase 5; recomputed on demand via `_utils.add_lognorm` | Saves ~50% disk at 600K-cell scale; project doc §3 advice |
| **HVG selection** | seurat_v3, batch_key=pool, n=3000 brain / 2000 placenta | seurat_v3 works on raw counts (compatible with scVI workflow); batch-aware avoids pool-driven HVGs |
| **HVG exclusions** | mito/ribo/hemo/sex-linked always; +Prl/Psg/Cgb/Cga for placenta | Technical artifacts shouldn't drive integration; pregnancy hormones dominate placenta variance |
| **Code layout** | Flat `scripts/` + `_utils.py`; no `src/` package | Per `INSTRUCTIONS.md`: simple > clever. Helpers added when duplication appears |
| **Object scope** | All ages in one AnnData per tissue | Consistent labels, cross-age trajectories |
| **UMAP** | Computed once on full data | Stable coordinates for all figures |
| **Contrast specification** | Declarative YAML, not hard-coded | Adding contrasts = config edit, not code change |
| **Composition + DE + pathway + LIANA+** | All driven by declarative contrasts | Uniform output schema across analyses |
| **DE method** | Pseudobulk + PyDESeq2, pup as statistical unit | Animal is the unit; cell-level DE is incorrect |
| **DE design (primary)** | `~ sex + pool + group` per age | Sex + pool as covariates; group is the 3-level factor of interest |
| **Pool-confounded contrasts** | Flagged `confounded_with_pool` in output | Honest about P1 Late-Stress + cross-age limitations (see §2) |
| **Cross-age analysis** | §8g persistence per Early-vs-Relaxed and Late-vs-Relaxed | Identifies developmental trajectory; flagged for pool-confounding |
| **Cross-tissue** | E12.5 placenta → Early-Stress brain cascade; E18.5 → Late-Stress brain cascade | Aligns placenta sampling window with stress exposure timing |
| **Sex-stratified contrasts** | Run but flagged `underpowered_exploratory` | Honest about n=2 limitations |
| **Trajectory** | PAGA + diffusion pseudotime only; NO velocity, NO CellRank | 10x Flex is probe-based (exon-only) — velocity not feasible; without velocity CellRank duplicates PAGA |
| **All ages equal in DPT** | Pool + per-age comparisons; no gating off "steady-state" ages | User principle: don't drop analyses for some ages, just tag caveats |
| **Cell-cell communication** | LIANA+ consensus, driven by contrasts | Multi-method robustness across all contrast families |
| **Provenance** | `manifest.json` + git hash + env lock from Phase 0 | Reproducibility wired in from line 1, not bolted on |
| **Environment** | uv + renv (not conda) | Conda channels blocked at firewall; PyPI/CRAN/Bioconductor reachable |
| **Python↔R bridge** | R called as subprocess (not rpy2) | Robust to crashes; easier to debug; no build-against-system-R fragility |
| **CellBender isolation** | Separate `.venv-cellbender/` venv | CellBender pins PyTorch hard; conflicts with scvi-tools if shared env |
| **Leiden backend** | igraph (`flavor="igraph", n_iterations=2, directed=False`) | Faster, scanpy's future default |
| **Phase 7 annotation** | Per-cluster majority voting, NOT per-cell argmax | CellTypist convention; cells in one cluster share an identity by definition |
| **Subcluster naming (7d)** | CellTypist majority + literature markers (`config/subcluster_markers.yaml`) aggregated per integer cluster | Cluster-level by construction; readable names for 7b integer ids |
| **TF activity (8c)** | decoupler ULM on DE Wald stats vs CollecTRI mouse network | Per contrast×celltype; BH-FDR within celltype×contrast; GR/HPA-axis readout for stress |
| **Composition tool** | propeller (R subprocess), NOT scCODA | scCODA's TF/arviz/numpy pins fought scanpy stack; propeller is cleaner + limma moderation is better for small n |
| **Dev workflow** | `dev_split_h5.py` writes 9 split h5 + `config/dev_split.yaml`; ALL phase scripts run unchanged | No dev-aware code in phase scripts; workstation uses `config/brain.yaml` with no splitter step |
| **Output layout** | Tables in `tables/<phase>/<phase>_<name>.csv`; plots in `plots/<phase>/...` | Phase-grouped, filename-identifiable when copied out |
| **Offline-audit CSVs** | 8b per-sample expr matrix, 8c leading-edge + TF activity, 8d PAGA edge diagnostics, 02_soupx summary (rho + pct_removed) | Troubleshoot results after workstation access ends |
| **Ask before strong scientific calls** | Don't drop analyses, exclude samples, or add tools without surfacing the trade-off | User principle: don't bake decisions in silently |
| **Phase 1 ambient correction** | **SoupX via R subprocess** (replaces abandoned CellBender, locked 2026-06-10) | CellBender's pickle bug unfixable without Docker; SoupX is the established R-Bioconductor alternative. `scran::quickCluster` + `autoEstCont` for data-driven rho per cluster; manual rho fallback available |
| **Brain CellTypist training** | cuML GPU LogReg (locked 2026-06-10) | ~9× speedup on class (123 min → 14 min); >40× on subclass (>18h CPU hang → 27 min GPU). CPU L-BFGS not parallelizable across 334 classes; cuML handles it on RTX 4500 Ada |
| **Brain CellTypist schema** | Three-tier per age: class (per-cluster majority, used by 8b/8c) + subclass (per-cell, 4W/3mo only, used by 7b/7d) + region (per-cell, 4W/3mo only, used by Phase 9) | Class is the canonical 8b/8c key. Subclass + region only at adult ages because P1 has cell types absent from ABC's adult taxonomy. P1 carries sentinel labels `no_subclass_model` / `no_region_model` |
| **Brain marker gate (Phase 7)** | STRICT canonical-marker gates for microglia / astrocyte / OL lineage / endothelial; demote failing labels to `unassigned_immune` / `unassigned_glia` / `unassigned_vascular` | CellTypist conf and marker gates measure different things; gate catches LogReg-confident calls that lack underlying biology. First run: 4 demotions (496 cells), 11 passes — sensible behaviour |
| **CellTypist sklearn-1.7 patch** | Sed-patch CellTypist's `train.py` to remove `multi_class='ovr'` (sklearn 1.7+ removed it; CellTypist still hardcodes) | Required before any retrain. Re-apply via `setup-remote.sh` if `uv sync` reverts the venv. Removing the kwarg also improves calibration (sklearn 1.7+ defaults to true multinomial softmax) |
| **renv Suggests workaround** | `renv::settings$package.dependency.fields(c("Depends","Imports","LinkingTo"), persist=TRUE)` + `USE_BUNDLED_LIBUV=1` | Prevents renv pulling Seurat→shiny→bslib→fs→libuv chain when installing SoupX. Persists across sessions |

---

## 12. Environment & Deployment

### Why uv + renv (not conda)

The workstation has **conda channels blocked at the corporate firewall** (verified: connections to `repo.anaconda.com` and `conda.anaconda.org` time out). However, **PyPI, CRAN, and Bioconductor are all reachable**. This rules out the entire conda family (conda, mamba, micromamba, pixi all hit the same blocked channels) and pushes us toward language-specific package managers.

The chosen stack:

| Concern | Tool | Why |
|---|---|---|
| Python deps | **uv** (Astral) | Single static binary, no admin needed, reads pyproject.toml, lockfile with hashes baked in. Roughly 5–10× faster than pip for env solves. |
| R deps | **renv** (Posit) | Project-local R library, reproducible via renv.lock. Standard tool in R community. |
| Python↔R interop | **subprocess + JSON/TSV** | Robust isolation; R crash doesn't kill Python; debuggable independently. |
| CellBender isolation | **separate uv venv** | CellBender's PyTorch pin conflicts with scvi-tools; running it in a sidecar venv avoids resolver hell. |

### What runs where

| Layer | Tool | Lock file (commit to git) |
|---|---|---|
| Main Python env | `uv` → `.venv/` | `uv.lock` |
| CellBender Python env | `uv` → `.venv-cellbender/` | (pinned inline in `setup-remote.sh`) |
| R packages | `renv` → `.renv-cache/` | `renv.lock` |
| Python version | `.python-version` → uv downloads CPython 3.11 if needed | (file itself is in git) |

### Bootstrap on the workstation

One command, idempotent:

```bash
./setup-remote.sh
```

This script:
1. Installs `uv` to `~/.local/bin` if not present (single curl-piped-to-sh, no admin)
2. Runs `uv sync` — creates `.venv/`, installs all Python deps from `uv.lock`
3. Installs R via `apt` if not present (requires sudo)
4. Runs `Rscript scripts/install-r-packages.R` — installs CRAN + Bioconductor packages, writes `renv.lock`
5. Creates `.venv-cellbender/` with CellBender + its pinned PyTorch

First run: 20–40 min (R compilation dominates). Subsequent runs from `uv.lock` + `renv.lock`: ~5 min.

### Day-to-day on the Mac (development)

```bash
# Edit code locally, push to remote
make sync                  # rsync ./ -> remote:~/snrna-project/ (excludes data, venvs)

# Run a phase on the remote
ssh remote-snRNA
cd ~/snrna-project
make validate              # Phase 0
make qc                    # Phase 2
make integration           # Phase 5 (GPU-heavy, run in tmux)
```

`make help` lists all available targets.

### Python ↔ R subprocess pattern

For each R-based step (scDblFinder, edgeR cross-check, CellChat), there's a paired pattern:

**Python side (`src/doublets.py`):**
```python
import subprocess, pandas as pd
result = subprocess.run([
    "Rscript", "scripts/run-scdblfinder.R",
    "--input", str(mtx_path),
    "--barcodes", str(barcodes_path),
    "--genes", str(genes_path),
    "--output", str(output_tsv),
    "--library", library_id,
], check=True, capture_output=True, text=True)
scores = pd.read_csv(output_tsv, sep="\t")
adata.obs = adata.obs.merge(scores, on="barcode", how="left")
```

**R side (`scripts/run-scdblfinder.R`):**
- Loads counts via Matrix::readMM
- Runs scDblFinder
- Writes TSV with `[barcode, doublet_score, doublet_class, library]`
- Python reads the TSV and joins to `.obs`

Rationale (rather than rpy2):
- **Process isolation:** an R crash, segfault, or memory issue doesn't kill the Python pipeline
- **Debuggable:** each R script runs standalone — `Rscript scripts/run-scdblfinder.R --input ... --output ...` works without any Python involvement
- **No build-against-system-R fragility:** rpy2 must be compiled against the exact R version installed, and breaks subtly when R is updated
- **Cross-language contract is explicit:** the TSV/JSON schema is the API; no implicit shared memory

The cost is some boilerplate (writing matrices to disk and reading back results), but for a pipeline that runs maybe a dozen times total, this is a rounding error.

### Updating dependencies

```bash
# Python: add a package
uv add scvelo              # updates pyproject.toml + uv.lock
git add pyproject.toml uv.lock && git commit -m "deps: add scvelo"

# R: add a package interactively
R
> renv::install("bioc::DropletUtils")   # bioc:: prefix for Bioconductor
> renv::snapshot()                       # updates renv.lock
git add renv.lock && git commit -m "deps: add DropletUtils"
```

### Long-term reproducibility

Anyone who clones the repo and runs `./setup-remote.sh` on a Linux x86_64 machine with network access to PyPI/CRAN/Bioconductor gets a near-identical environment. The lockfiles pin every transitive dependency. For paper revisions in 18 months, this gives you:

- **`uv.lock`** — exact Python package versions with hashes
- **`renv.lock`** — exact R package versions
- **`pyproject.toml`** + `scripts/install-r-packages.R` — top-level intent, human-readable
- **Git history** — track every env change

The only thing not pinned is system libraries (glibc, CUDA runtime). For absolute reproducibility, you would containerize — but given Apptainer isn't currently available on the workstation and conda channels are blocked, the uv + renv approach is the pragmatic best.

### Fallback if PyPI/CRAN/Bioc ever get blocked

If the firewall tightens further:
- **PyPI blocked → uv broken.** Fallback: build wheels off-site, transfer with `uv pip install <wheel.whl>` or set up a local devpi mirror.
- **CRAN blocked → R install broken.** Fallback: use Posit Package Manager (`packagemanager.posit.co`) which is sometimes whitelisted when CRAN isn't, or transfer R package tarballs.
- **All blocked + Apptainer added later → containerize.** Build off-site via GitHub Actions (Linux x86_64 native), transfer the `.sif` file.

For now, the current network policy is workable.

### Files in this layer

| File | Purpose |
|---|---|
| `pyproject.toml` | Python deps + project metadata |
| `uv.lock` | Generated by `uv sync`; commit to git |
| `.python-version` | Pins to 3.11 |
| `setup-remote.sh` | One-shot bootstrap on the workstation |
| `scripts/install-r-packages.R` | R-side installer (CRAN + Bioconductor + GitHub for CellChat) |
| `renv.lock` | Generated by `renv::snapshot()`; commit to git |
| `Makefile` | Convenience commands: `make setup`, `make sync`, `make validate`, etc. |
| `.gitignore` | Excludes venvs, caches, data, results |
| `scripts/run_soupx.R` | R worker for Phase 1 SoupX (added 2026-06-10) |
| `scripts/02_soupx.py` | Python orchestrator for Phase 1 SoupX (added 2026-06-10) |

### GPU stack additions (2026-06-10)

- **cuML installation:** `uv pip install cuml-cu12 cudf-cu12 --index-url https://pypi.nvidia.com --extra-index-url https://pypi.org/simple/`. Required for GPU-accelerated CellTypist training (`use_GPU=True`). Forces numba 0.65→0.64, pyarrow 24→23, cuda-toolkit 13→12.9. CUDA 13 driver is backward-compatible with the 12.9 runtime, so no driver change needed.
- **CellTypist sklearn-1.7 patch** (per §1b note 25) is applied via sed in `setup-remote.sh` Step 2.5. Idempotent. Re-run `./setup-remote.sh --skip-references --skip-r --skip-cellbender` whenever `uv sync` touches the venv.

---

## 13. Ambient RNA Correction Strategy (locked 2026-06-10)

### Why we need it

snRNA-seq nuclei preps release cytoplasmic and erythroid RNA into the lysis buffer, contaminating every droplet with hemoglobin (Hbb-bs / Hbb-bt / Hba-a1 / Hba-a2 / Alas2) and brain-specific abundant transcripts (Malat1, mitochondrial). For P1 brain specifically, this is severe: nucleated erythroblasts in residual vasculature rupture readily, and `pct_hemo ≤ 5%` QC catches only the worst cases. First Phase 7 brain run produced ~81K P1 cells (~45% of P1) labeled "Blood: Erythrocyte" — most of those are contaminated neurons / glia, not real erythroblasts.

### Tool choice: SoupX (R), not CellBender (Python)

CellBender was abandoned 2026-06-05 after an unresolvable `weakref.ReferenceType` pickle bug across torch 1.13.1 / 2.0.1 / 2.12 + pyro 1.8.6 + numpy <2/>=2 + cellbender 0.3.0 / 0.3.2-master. Only the official Docker image works; we don't have Apptainer/Docker on the workstation. SoupX is the established R-Bioconductor alternative — same purpose, completely different toolstack (no torch/pyro/numpy pickling).

### Workflow

| Step | Tool | Purpose |
|---|---|---|
| Load filtered counts | `DropletUtils::read10xCounts` | per-sample filtered .h5 (cells) |
| Load raw counts | `DropletUtils::read10xCounts` | per-sample raw (all droplets, soup profile source) |
| Build SoupChannel | `SoupX::SoupChannel(tod, toc)` | filtered + raw matrices |
| Quick clusters | `scran::quickCluster` | per-sample rough clustering for autoEstCont |
| Estimate contamination | `SoupX::autoEstCont` | data-driven per-cluster rho |
| Adjust counts | `SoupX::adjustCounts(sc, roundToInt=TRUE)` | corrected count matrix |
| Write outputs | MTX + barcodes.tsv + features.tsv + JSON | Python orchestrator assembles per-sample h5ad |

### Manual rho fallback

If `scran` install fails, `02_soupx.py --rho 0.10` bypasses clustering and uses a fixed 10% contamination fraction. Less data-driven but functional; reasonable starting value for snRNA-seq brain (typically 5–15% post-autoEst).

### Output organisation

- `results/{tissue}/h5ad/02_soupx_corrected/{sample_id}.h5ad` — corrected counts + sample metadata in obs + SoupX summary in uns
- `results/{tissue}/tables/02_soupx/02_soupx_summary.csv` — per-sample rho_mean / rho_min / rho_max / n_cells / pct_removed / n_clusters / elapsed_sec / mode (autoEst | manual)
- `results/{tissue}/logs/02_soupx/{sample_id}.log` — R subprocess stdout/stderr

### Sanity bounds (for spot-checking SoupX outputs)

- **`rho_mean`:** 0.02–0.10 typical for adult brain; 0.05–0.15 for P1 (more ambient given blood). > 0.30 = SoupX over-correcting, investigate.
- **`pct_removed`:** 2–8% for adult; 5–15% for P1. >30% = over-correction.
- **Post-SoupX `pct_hemo` mean/median:** should drop substantially in non-erythroid clusters. P1 cluster 0 / 19 / 20 / 45 hemoglobin signal should mostly disappear.

### Smoke-test gate

Per smoke-test policy: ONE sample first (`--sample-ids E1 --n-jobs 1`), verify sanity bounds, then production run with `--n-jobs 6`. Smoke test = ~5–15 min; full 57-sample production with `--n-jobs 6` on 56 cores = ~100 min wall-clock.

### Downstream wiring (TBD)

`02_qc.py` needs a prefer-soupx fallback: if `02_soupx_corrected/{id}.h5ad` exists, read that; else fall back to cellranger filtered h5. SoupX changes counts, so Phase 2 onwards must be re-run end-to-end (HVG selection, scVI training, clustering, annotation all change).
