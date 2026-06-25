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
- **Three groups, not two** — this changes every contrast (Early-vs-Relaxed and Late-vs-Relaxed as primary contrasts, Early-vs-Late secondary)
- **P1 Late Stress has no females** (only 2 males) — these samples are also confounded with Pool 3 (see §2)
- **Placenta has no complete factorial** — E12.5 has Early+Relaxed only, E18.5 has Late+Relaxed only. Cross-age placenta comparisons are not analyzable.
- **Placenta E12.5 sex is undetermined** at sampling — inferred from Y-chromosome expression
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
- Total cells post-QC: ~785K brain, ~549K placenta pre-QC nuclei (SoupX-corrected); ~661K brain + ~397K placenta after Phase 7 annotation and Phase 8 contaminant/unassigned drop

**Brain region:** Whole brain (not microdissected)

---

## 1b. Implementation Status

**All phases through 8a implemented and run on the workstation; 8b main DE done; 8b follow-ups (disruption + shuffle null) done for brain main + 3 focal subclusters; 8c–8d smoke-tested; 8e PRODUCTION-COMPLETE (compute + all plotting, both tissues × broad+subtype node schemes, 2026-06-24); 8f–8g implemented and smoke-tested, production wiring in progress.** Cross-tissue (8f) and cross-age (8g) views operate on completed 8b/8c CSVs so they need no re-runs.

> **Workstation target** (where production runs go; see §3 for detail):
> Linux box, **258 GB RAM, 56 CPU cores, 1× NVIDIA RTX 4500 Ada (24 GB VRAM)**.
> GPU and CPU compute on the same host. Conda blocked at firewall — use `uv` + `renv` only. R + Rscript on PATH. **CellBender abandoned 2026-06-05 (pickle bug); replaced by SoupX via R subprocess (2026-06-10).** Everything runs from the main uv-managed `.venv/` (no sidecar venvs). cuML installed via NVIDIA PyPI for GPU LogReg (CellTypist training, locked 2026-06-10). Runbook: `run_pipeline_WS.sh`.

| Phase | Status | Script(s) |
|---|---|---|
| 0 Validation | ✓ done | `01_validate.py` |
| 1 Ambient RNA (SoupX) | ✓ done (full production run 2026-06-10) | `02_soupx.py` + `run_soupx.R` |
| 2 Per-sample QC | ✓ done | `02_qc.py` |
| 3 Doublet detection | ✓ done (per-pool via R subprocess) | `03_doublets.py` + `run_scdblfinder.R` |
| 4 Concat + HVG + cell cycle | ✓ done | `04_integration_prep.py` |
| 5 scVI integration | ✓ done (GPU workstation, BF16) | `05_integration.py` |
| 6 Clustering (Leiden, igraph) | ✓ done | `06_clustering.py` |
| 7 Annotation (per-cluster majority; P1 via scANVI) | ✓ done | `07_annotation.py` + `run_scanvi_p1.py` |
| 7b Subclustering | ✓ done | `07b_subcluster.py` |
| 7d Subcluster annotation | ✓ done | `07d_subcluster_annotate.py` + `config/subcluster_markers.yaml` |
| 7e Cell-type counts diagnostic | ✓ done | `07e_celltype_counts.py` |
| 8a Composition (propeller) | ✓ done, both tissues (strata + parallel + unassigned-drop, 2026-06-15) | `08a_composition.py` + `run_propeller.R` |
| 8b Pseudobulk DE (PyDESeq2) | ✓ done — brain main + 7 focal subclusters | `08b_de.py` |
| 8b Summary plots (heatmap/upset/bar/bubble/rrho/dotplot/grid/venn) | ✓ done | `08b_de_summary.py` |
| 8b Developmental disruption analysis | ✓ done (brain main + 3 focal subclusters) | `08b_developmental_disruption.py` |
| 8b Follow-up plots (disruption + consistency) | ✓ done (brain main + 3 focal subclusters) | `08b_followup_plots.py` |
| 8b Disruption shuffle test (k-preserving null + within-stratum binomial) | ✓ done (brain main + 3 focal subclusters) | `08b_disruption_shuffle_test.py` |
| 8c GSEA + leading-edge + TF activity | ✓ smoke-tested (decoupler ULM + CollecTRI) | `08c_pathways.py` + `fetch_genesets.R` |
| 8d Trajectory (PAGA + DPT) | ✓ smoke-tested | `08d_trajectory.py` |
| 8e Cell-cell communication (compute) | ✓ DONE — both tissues, broad + subtype node schemes, 3 arms, n_perms=1000 (2026-06-24) | `08e_comms.py` + `_08e_*` workers |
| 8e Communication (plotting) | ✓ DONE — both tissues × both schemes; 7 plot families incl. focal-fan grids + per-pathway LR detail (2026-06-24) | `08e_comms_summary.py` + `_08e_plots_{baseline,differential,stats,pathway}.py` |
| 8f Cross-tissue | ✓ smoke-tested (six views) | `08f_cross_tissue.py` |
| 8g Cross-age / persistence | ✓ smoke-tested (brain only) | `08g_cross_age.py` |
| 9 Cross-species RRHO2 | loaders are stubs | `09_cross_species_validation.py` |

**Key implementation notes (selected; full history in git):**

1. **Flat scripts layout, no `src/snrna/` package.** Shared helpers in `scripts/_utils.py`: `load_config`, `load_contrasts`, `phase_paths`, `phase_table_dir`, `add_lognorm`, `select_accelerator`, **`iter_strata`, `parallel_map`, `unassigned_mask`** (last three added during 8a). Phase scripts are numbered standalone files.

2. **Raw counts in `.X`, lognorm not persisted.** Phase 4 builds `.layers["lognorm"]`, Phase 5 drops it. Recompute via `_utils.add_lognorm(adata)`.

3. **QC: per-sample MAD + hard floors + hard caps + cohort-outlier flag.** `min_counts=500`, `min_genes=200`, `pct_mt≤1.0`, `pct_hemo≤5.0`.

4. **scDblFinder per pool**, `samples=` arg so simulated doublets respect within-sample boundaries.

5. **HVG exclusion lists.** mito/ribo/hemoglobin/sex-linked always; placenta adds Prl*/Psg*/Cgb*/Cga.

6. **scVI uses BF16 mixed precision on GPU.** Auto-detected via `_utils.select_accelerator()`.

7. **Leiden uses igraph backend** (`flavor="igraph", n_iterations=2, directed=False`).

8. **Phase 7 = per-cluster majority voting.** P1 brain annotated by scANVI from Rosenberg 2018 (Di Bella abandoned — cortex-only, mislabeled 42% of whole-brain P1 as erythrocyte). 4W/3mo via ABC CellTypist. Brain marker gate demotes unsupported calls to `unassigned_*`; Phase 8 then DROPS those.

9. **scCODA abandoned for composition (8a).** Replaced with **propeller (speckle+limma) via R subprocess** — clean Bioconductor install, limma's empirical-Bayes moderation better for small n.

10. **Phase 1 = SoupX via R subprocess (locked 2026-06-10).** CellBender abandoned 2026-06-05 (pickle bug). SoupX rho ranges all healthy (no `rho>0.30` outliers); full 57-sample run done.

11. **Brain CellTypist models retrained on GPU via cuML (locked 2026-06-10).** Three per-age models from ABC adult atlas. sklearn-1.7 patch (sed) required before training.

12. **Phase 8a finalized (2026-06-15):** propeller composition now runs the full sex-strata schema (combined + M + F, via `iter_strata`) on every per-age group/omnibus contrast, parallelized across slices via `_utils.parallel_map` (`--n-jobs`; ~3 h serial → minutes). Contaminants + `unassigned*` cells dropped from numerator and denominator (recorded in `08a_dropped_cells_per_donor.csv`); `pool` dropped + flagged `confounded_with_pool` only where perfectly aliased with `group` (P1 Late, and one placenta E12.5 sex stratum). `min_donors=2` to run, `<3` → `low_n`. Change heatmaps **bold-outline FDR<0.05 cells**; makeup bars stay descriptive. This is the template for 8b–8g.

13. **Phase 8b follow-ups locked (2026-06-15) — disruption + shuffle null.** Three brain-only scripts consume the master `08b_de_results.csv` (no DE re-run) and produce the headline disruption analysis. Findings: every brain broad cell type shows **R∩E and R∩L enriched (↑\*\*\*) and E∩L depleted (↓\*\*\*)** under a k-preserving null — i.e., when age-DE signal is shared across two groups, Relaxed is almost always one of them. See INSTRUCTIONS.md §"Disruption analysis framing" and §5 below for the full claim.

14. **Phase 8e production-complete (compute + plotting, locked 2026-06-24).** Split compute/plot like 8b/8c. Compute (`08e_comms.py`) ran both tissues × broad+subtype node schemes, three arms (baseline/differential/perdonor), n_perms=1000. **Key result: placenta differential = 447 sig LR pairs (`interaction_padj`); brain differential = NULL at every level/contrast** (receptor-side DE too weak at coarse grouping — confirmed real, not a bug). Plotting (`08e_comms_summary.py`) produces 7 families; the readable graph is the **focal-fan small-multiples grid** (one cell type pinned per panel) — hairball all-pairs graphs kept only as supplementary. Per-pathway **LR dotplot + ranked lollipop** answer "which LR pairs changed, between which cells." A **slice-specific adaptive effect floor** (each plot computes its own cutoff from its own pairs; `--q-stat 0.25` differential, `--q-delta 0.75` baseline) controls density, applied to all Δ/aggregating/per-pathway plots but NOT distribution plots (volcano, rank-rank). Sex NOT stratified (n≈2/group degenerate). See §5 + INSTRUCTIONS.md §"Phase 8e plotting".

---

## 2. Critical Considerations & Caveats

### Statistical Power Limitations

**n=2 per sex per condition per age is a real limitation:**
- Cannot reliably test sex × condition interactions
- Pseudobulk DE with n=2 vs n=2 detects only large effect sizes (logFC > 1.5–2)
- Single-cell-level DE treating cells as replicates is **incorrect** — reviewers will flag this

**Decision: combined analysis primary, sex strata run systematically but flagged**
- **Primary (`combined` stratum):** pool sexes per group (effective n=4 vs n=4 per age per condition pair), with sex as a covariate. DESeq2 design: **`~ sex + pool + group`** (sex + pool as nuisance variables; group is the 3-level factor Early/Late/Relaxed).
- **Sex strata (`M`, `F`):** run for every contrast via the declarative strata mechanism (see §8), but flagged `low_n` / `underpowered_exploratory` — within a single sex, n drops to ~2/group, so these are exploratory.
- Sex × group interactions reported as exploratory, underpowered.

**No dam ID — litter random effect cannot be modeled:**
- Prenatal stress was applied to **dams**, not pups; pups from one stressed dam are **not fully independent** for testing the stress effect
- Without dam IDs we cannot fit `~ sex + pool + group + (1|dam)` and must treat each pup as an independent biological replicate
- This is anti-conservative: p-values are slightly optimistic, effect sizes are accurate
- **Mitigation:** include `pool` (= harvest+library batch) as a covariate where estimable (pups harvested same day often from same/few dams), use FDR < 0.05, require reasonable effect sizes
- **Explicit caveat in methods:** "Dam identity was not recorded; pseudobulk DE treats each pup as an independent biological replicate, which may inflate statistical significance for traits that aggregate at the litter level."

### 10x Flex-Specific Considerations

1. **Probe-based capture (not poly-A)** — biased toward exonic sequences
2. **RNA velocity is not feasible** — probe-based (exon-only); no spliced/unspliced. Trajectory = PAGA + diffusion pseudotime only (no velocity, no CellRank).
3. **Ambient RNA from probe leakage** — handled by SoupX (Phase 1)
4. **Doublet rates** comparable to standard 3' chemistry

### snRNA-seq-Specific Considerations

1. **High intronic read fraction** is normal
2. **Mitochondrial % should be near zero** in good nuclei prep
3. **Ambient RNA is severe** — placenta especially (hemoglobin); brain (Malat1, mito). Corrected by SoupX.
4. **Cell type annotation references** built from scRNA-seq may not transfer perfectly

### Cell-Cell Communication on snRNA-seq

Defensible and published but framed carefully. Use **LIANA+** (consensus). Particularly well-motivated for **placenta** (signaling-heavy). Acknowledge: inferring signaling potential from nuclear transcriptional state. **No cross-tissue CCC** (BBB) — the placenta→brain link is the LR-from-DE endocrine framing in 8f view 4.

### Pool/Library Confounding — REAL AND CONSEQUENTIAL

**What's analyzable cleanly:**
- 4-week brain: Early vs Late vs Relaxed (all 3 groups in both Pool 1 and Pool 2 — best-balanced age)
- 3-month brain: Early vs Late vs Relaxed (single pool, no within-age pool confound)
- P1 brain: Early vs Relaxed (both in Pool 2)
- E18.5 placenta: Late Stress vs Relaxed (all in Pool 4)

**What's confounded and must be flagged (`confounded_with_pool`):**
- **P1 brain Late Stress vs anything:** the only 2 P1 Late Stress samples are in Pool 3, separate from P1 Early/Relaxed in Pool 2. Pool ≡ group at this age → `pool` is dropped from the design and the contrast flagged (see §8 pool-aliasing rule).
- **Brain age trajectories (P1 → 4W → 3-month):** each age dominated by a different pool. scVI integration partially corrects; interpretation requires caution. **8b follow-ups partially mitigate:** the k-preserving shuffle null in `08b_disruption_shuffle_test.py` preserves each gene's pool-driven `k_i` value, so signal exceeding that null is biology beyond pool structure.
- **Pool 3 mixes tissues:** 2 brain + 14 placenta samples in one library. The 2 P1 Late Stress brain samples may carry placental signatures.
- **Placenta cross-age (E12.5 vs E18.5):** each age = one pool + different conditions. **Cannot do meaningful E12.5 vs E18.5 comparison.**
- **E12.5 placenta Early vs Relaxed:** mostly Pool 3, with 2 Relaxed in Pool 4. Acceptable; in one sex stratum the Pool-4 Relaxed cause pool to alias with group → pool dropped + flagged there.

**Mitigations:**
- Include `pool` as a covariate in DESeq2 designs where estimable (`~ sex + pool + group`); drop + flag where aliased with group.
- For scVI integration: `batch_key=pool`.
- Tag underpowered/confounded contrasts; downweight in interpretation.

**Pool 4 ambient signature differs (SoupX 2026-06-10):** E18.5 Relaxed (LCP*) cluster at rho 0.075-0.085, higher than E18.5 Late (LSP*, 0.034-0.067, same Pool 4) → real prep/age effect, not pure batch artifact. Controlled by `~ sex + pool + group`. Also note placenta E18.5 has a large contaminant fraction (~18% Late vs ~14% Relaxed = DSC↔trophoblast bleed), surfaced by the 8a dropped-cells diagnostic and dropped from the test.

---

## 3. Compute Environment

**Available machine:** 258 GB RAM, 56 CPU cores, **1× NVIDIA RTX 4500 Ada (24 GB VRAM)**, GPU on the same box. Likely a shared lab workstation.

**Implications:**
- SoupX per-sample (R): parallel via `parallel_map`, ~5-15 GB RAM per concurrent sample
- scVI on GPU with BF16: ~2-3 hrs per integrated tissue object
- Propeller (8a) and other R-subprocess loops: parallelized across cores via `parallel_map` (`--n-jobs 24`)
- 8b shuffle test (8b follow-up): CPU-bound numpy permutations, parallelized as processes via `parallel_map(use_threads=False)`. Brain main + 3 subclusters complete in ~5-10 min at `--n-jobs 16`.
- Total wall time per tissue: ~1-1.5 days unattended

**System RAM budget (per-tissue, ~720K input cells):** worst-case peak ~90 GB during scVI → ~2.5× headroom on 258 GB.

**GPU VRAM budget (24 GB):** scVI at `batch_size=1024` BF16: 8-12 GB. **Do NOT use rapids-singlecell** at this cell count — neighbors/UMAP on 600K cells needs 15-22 GB VRAM with no margin. CPU runs the same in 10-30 min on scVI's 30-dim latent.

**Key engineering principles:** sparse matrices throughout (check `issparse` before any densifying op); don't carry redundant layers; process samples sequentially through Phase 1-3 then concatenate; process tissues sequentially; save checkpoints after every phase; serialize GPU phases; wrap long jobs in `tmux`; set `CUDA_VISIBLE_DEVICES=0`.

---

## 4. Analysis Strategy

### Ecosystem: Python (Scanpy + scvi-tools), R via subprocess

**Rationale:** scVI/scANVI superior for complex batch structure; better Python tooling for snRNA-seq quirks; pseudobulk DE via PyDESeq2; scales to ~1.5M nuclei. R called as subprocess (not rpy2) for scDblFinder, propeller, SoupX.

### Joint Integration, Per-Age Analysis

| Step | Scope | Rationale |
|---|---|---|
| Integration (scVI) | All ages together | One consistent latent space |
| Clustering + annotation | Joint object | Canonical cell type labels |
| UMAP | Computed once on full data | Stable coordinates across all figures |
| Composition analysis | Per age (+ per region, brain) | Cell type proportions are age/region-specific |
| Differential expression | Per age | `~ sex + pool + group` within each age (3-group factor) |
| Disruption + shuffle (8b follow-ups) | Per (sex × level × celltype), brain only | Needs ≥2 ages = `within_group_across_age` contrast |
| Trajectory analysis | Cross-age | Maturation lineages — caveat: pool-age confounded |
| Cell-cell communication | Per age, per group | Differential signaling |

**Note on the 3-group factor:** group ∈ {Early Stress, Late Stress, Relaxed}. Primary contrasts Early-vs-Relaxed and Late-vs-Relaxed (Relaxed reference). Early-vs-Late secondary. Sex handled by the declarative strata (combined + M + F) — see §8.

### Reference Atlases (per age)

**Brain:**
- **P1:** Rosenberg et al. 2018 (Science, GSE110823, P2 brain) via scANVI label transfer — PRIMARY (Di Bella 2021 cortex-only abandoned: mislabeled 42% of whole-brain P1 as erythrocyte)
- **4W & 3-month:** Allen Brain Cell Atlas (Yao et al. 2023) — definitive adult mouse reference

**Placenta:**
- Curated literature markers + STAMP Spearman correlation against Liu 2024 reference (35 cell types, E9.5-E18.5). No CellTypist model.

### Pre-Registration

Document before looking at results: focal cell types (microglia, oligo lineage, exc/inh neurons, astrocytes), primary vs. secondary contrasts, significance thresholds, validation criteria. Protects against fishing accusations given the n=2 constraint.

---

## 5. Detailed Analysis Pipeline

(Phases 0–7 as implemented; see §1b for status and §"Annotation conventions"/§"Phase 1 = SoupX" in INSTRUCTIONS.md for locked details. Phase 0 validation is a mandatory gate: manifest validation, Y-chromosome sex check, per-sample fingerprints, gene-set sanity, empty-droplet check.)

### Phase 7e: Cell-type Counts Diagnostic (`07e_celltype_counts.py`) — DONE

Per-donor × cell-type count CSV. Brain: 3 granularities (`celltypist_broad` / `_class` / `_subclass`) × (`whole` + per `celltypist_region`). Placenta: `celltype_majority` × `whole` only. Long-form columns: `tissue, donor_id, sample_id, age, group, sex, pool, granularity, level, region, celltype, n_cells, category`. Used to sanity-check 8a propeller inputs (the rows are exactly the propeller numerator) and for paper Table S?. Outputs `tables/07_annotation/07e_celltype_counts.csv`. Production run 2026-06-15: 32,428 rows brain (668,969 cells, 0.33% unassigned), 478 rows placenta (449,467 cells, 3.87% unassigned).

### Phase 8: Downstream Biology

**Architectural principle:** All contrasts (composition, DE, pathway, communication) are defined **declaratively in YAML** (`contrasts:` + `strata:`), not hard-coded. The downstream engine iterates over the contrast × sex-stratum specification, runs the same machinery for each, and writes uniformly-structured output tables tagged with `sex`. Adding a new contrast = edit YAML.

**Reference level for `group`:** Relaxed (positive logFC = upregulated in stress).

**Cross-cutting Phase-8 rules (see INSTRUCTIONS.md "Phase 8 conventions"):** statistical unit = animal (donor_id); sex strata combined/M/F on every contrast; contaminants + `unassigned*` dropped (num + denom) and recorded in a diagnostic CSV; `~ sex + pool + group` with pool dropped+flagged where aliased; `min_donors=2`, `<3` → `low_n`; per-item parallelism via `_utils.parallel_map` + `--n-jobs`.

#### 8a. Composition (`08a_composition.py` + `run_propeller.R`) — DONE
- **propeller (speckle+limma)** via R subprocess — NOT scCODA (abandoned; dependency stack fought scanpy).
- Per-donor cell-type counts (animal as statistical unit). Grid per test: sex stratum {combined, M, F} × contrast × age × level × granularity, **parallelized across slices** with `parallel_map`.
- **Levels:** `whole` + (brain) each `celltypist_region`; placenta `whole` only. **Denominator = cells in the slice** (region → ÷ donor's cells in that region; whole → ÷ donor's total cells).
- **Granularities:** brain `celltypist_broad` + `celltypist_class` + subtype (focal coarse types exploded to `subcluster_name`); placenta `celltype_majority` + subtype.
- **Cleaning:** contaminants (`Contamination_*`/`unresolved`) + `unassigned*` dropped from numerator and denominator; dropped per-donor counts/fractions → `08a_dropped_cells_per_donor.csv`.
- **Pool:** `~ sex + pool + group`; pool dropped + flagged `confounded_with_pool` where aliased with group. `min_donors=2`, `<3` → `low_n`.
- **Outputs:** master `08a_composition_results.csv` [tissue, sex, contrast, flag, age, level, granularity, category, test_type, prop_ratio, log2_prop_ratio, statistic, pvalue, fdr, reliability, note]; change heatmaps (`heatmaps/{granularity}/{sex}/{contrast}_{age}.png`; color = log2(prop ratio); cols = whole + regions; **FDR<0.05 cells bold-outlined**); descriptive stacked makeup bars; dropped-cells diagnostic. **Read `reliability==ok` + finite `log2_prop_ratio` for real hits.**
- **Verified hits:** brain — P1 cortical OPC/Oligo, Astrocytes/Ependymal, Radial_Glia, Choroid_Plexus DOWN (Early); 3mo CTXsp OEC / Isocortex Vascular DOWN. placenta — E18.5 Late: B cell / Monocyte / Lymphoid / Dendritic DOWN, Yolk-sac epithelial UP; E12.5 B cell DOWN.

#### 8b. Differential Expression (`08b_de.py`) — DONE
Pseudobulk only, never single-cell-level. Animal (pup, via `donor_id`) is the statistical unit. Sum counts per donor per cell type; filter ≥10 cells in ≥3 samples per group; **PyDESeq2**; iterate over all contrasts × sex strata; design `~ sex + pool + group` (pool dropped+flagged where aliased); contaminant + `unassigned_mask` drop before pseudobulk; master `de_results.csv` [contrast, sex, contrast_family, celltype, gene, logFC, padj, direction, flag]. **Plots:** volcano per cell type with **gene symbols** (`var['symbol']`), top genes labeled, contrast + thresholds on plot. Offline-audit per-sample expr matrix CSV. **Sig thresholds LOCKED:** `padj<0.05 & |log2FC|>1`.

**8b summary plots (`08b_de_summary.py`):** 8 plot types — heatmap, upset, bar, bubble, rrho, dotplot, grid (volcano grid with gene labels), venn. Master CSV is unfiltered; a `BLOCKLIST_FOR_VIZ` constant (17 hemoglobin + sex-linked genes + `mt-*` prefix) keeps headline figures focused on stress biology. `--no-blocklist` for QA.

**Brain main DE production result (2026-06-15):** 20.3M rows total, 763K sig hits. Subclusters: immune 7,963 sig; opc_oligodendrocytes 30,029; astrocytes_ependymal 23,932.

#### 8b follow-ups (locked 2026-06-15) — DONE
Three brain-only scripts that consume the master `08b_de_results.csv` (no DE re-run needed). Skip placenta cleanly — placenta has no `within_group_across_age` contrast (incomplete cross-age factorial; see §2).

**8b-i Developmental disruption (`08b_developmental_disruption.py`):** Classifies genes from the `within_group_across_age` contrast (per group; pairwise age tests collapsed to most-significant per gene) into 5 direction classes:
- `universal` — sig in all 3 groups (developmental baseline)
- `relaxed_only` (`= LOST`) — sig only in Relaxed
- `stress_shared` (`= GAINED`) — sig in BOTH Early AND Late, NOT in Relaxed
- `early_only` / `late_only` — sig in one stress group

Outputs: `08b_developmental_disruption_summary.csv` (counts + mean `|LFC|` per group for LOST class) + `08b_developmental_disruption_genes.csv` (long-form gene-level direction class assignments).

**8b-ii Follow-up plots (`08b_followup_plots.py`):** Two plot types per `(sex × level)`:
- Mirror disruption bar (LOST left red, GAINED right blue) + paired `|LFC|` boxplots showing effect-size collapse for LOST-trajectory genes under stress (Relaxed gray / Early red / Late blue).
- Stress-consistency stacked bars per age — Early-only (red) / Both-sig=convergent (gray) / Late-only (blue). Auto-skips when one of two stress contrasts is missing.

**8b-iii Shuffle null test (`08b_disruption_shuffle_test.py`):** k-preserving permutation null + analytic binomial tests + within-stratum chi-square. **This is the inferential validation of the disruption claim.**
- For each gene, keep `k_i` = #groups in which it's sig (0/1/2/3) but randomize WHICH groups. Vectorized via per-row argsort.
- Per-category analytic binomial p-values for all 6 disjoint sig-pattern categories (R-only / E-only / L-only / R∩E / R∩L / E∩L) modelling `obs ~ Binom(n_k_stratum, 1/3)`. Enrichment + depletion p-values, BH-corrected within each direction.
- Within-stratum chi-square goodness-of-fit (k=1 trio and k=2 trio each tested vs uniform).
- Headline figure: 2-panel. Panel A = mirror `|Δ|` bar (LOST left, GAINED right; solid color for Δ>0, faded for Δ<0). Panel B = within-stratum 6-bar breakdown per cell type with dashed reference lines at `n_k1/3` and `n_k2/3`.
- **Headline biological finding (brain combined×whole):** ALL 7 broad cell types show massively enriched R∩E and R∩L (k=2 trio, ↑\*\*\*) and depleted E∩L (↓\*\*\*). 3 cell types additionally show enriched LOST (Olfactory ensheathing ↑\*\*, Astrocytes/Ependymal ↑\*\*, Vascular ↑\*). Subcluster level: PAM_ATM_Microglia ↑\*\*, BAM ↑\*\*, Protoplasmic_Astrocyte ↑\*\*, OPC ↑\*, MFOL ↑\*, Homeostatic_Microglia ↑\*\*. See INSTRUCTIONS.md §"Disruption analysis framing" for the paper-quality claim.

#### 8c. Pathway Analysis (`08c_pathways.py`)
GSEA on ranked DE statistics (decoupler + mouse MSigDB via msigdbr, FDR within collection); stress-relevant gene sets (GR targets, HPA, neuroinflammation, synaptic, mito, ER-stress, oxphos); **TF activity (CollecTRI, REQUIRED — `--tf`)**; per-cell pathway scoring. Runs against every contrast × sex stratum.

#### 8d. Trajectory (`08d_trajectory.py`) — Brain primarily
PAGA (always) + diffusion pseudotime anchored at progenitors. **No velocity, no CellRank** (Flex probe-based). Focal lineages: oligo (OPC→OL), microglia, astrocyte. PAGA edges = hypotheses (edge diagnostics CSV).

#### 8e. Cell-Cell Communication (`08e_comms.py` compute + `08e_comms_summary.py` plotting) — PRODUCTION-COMPLETE (2026-06-24)
LIANA-py 1.7.3, split compute/plot (mirrors 8b/8c). Covers ES-v-Rel, LS-v-Rel, ES-v-LS. Placenta focus: trophoblast ↔ decidua ↔ fetal vasculature. Brain focus: neuron ↔ glia, microglia ↔ neuron.

**Three arms** — `baseline` (`rank_aggregate` per group×age, pooled cells) is DESCRIPTIVE only (pseudoreplication, not a stress test) but carries `specificity_fdr` (BH within group×age×level, from `cellphone_pvals` perms); `differential` (`df_to_lr` on 8b Wald stats) is the PRIMARY inferential arm, sig col `interaction_padj`, animal-respecting via 8b; `perdonor` (`rank_aggregate` per donor → MW-U) is corroboration (null at n≈4, expected). **Placenta differential = 447 sig pairs; brain differential = NULL at every level/contrast** (receptor-side DE too weak at coarse grouping — confirmed real, not a bug).

**Node schemes** (`--node-scheme {broad,subtype}`): broad = canonical key (dir `08e_communication/`); subtype = `comms_subtype` (focal-subcluster substates ≥300 cells as nodes, smaller + non-focal collapsed to parent broad; dir `08e_communication_subtype/`; focal map in `config/stress_pathways_8e.yaml`). **Levels:** brain whole + regional (`celltypist_region`); placenta whole only; differential arm is whole-only. **Sex:** combined-only — NOT stratified (n≈2/group per sex makes all three arms degenerate); stated as a methods limitation.

**Plotting (7 families, CSV-only no recompute):** `01_overview`, `02_baseline` (descriptive landscape; hairball all-pairs graphs kept as supplementary only), `03_differential` (PRIMARY, FDR-backed: volcano/dotplot/Δ-network/Δ-chord/Δ-celltype-heatmap/sender-receiver up-down bars/stress-signature heatmap), `04_sender_receiver`, `05_per_donor` (null corroboration), `06_by_pathway` (per stress pathway: chord+network graphs + **LR dotplot** + **ranked LR lollipop** + **cross-scheme companion** broad-vs-subtype), `07_focal_grids`. **The readable graph format = focal-fan small-multiples grid** (one cell type pinned per panel, edges fan to all others) — per-group descriptive + 4 Δ variants (count/magnitude × baseline/differential), sqrt-scaled widths, `unassigned_*` dropped, across whole + per-pathway + regional + subtype.

**Pathway → gene set:** `graph_pathways` in `config/stress_pathways_8e.yaml` (9 brain + 11 placenta MH Hallmarks) × that pathway's 8c leading-edge genes (FDR<0.05, level=='whole'); LR pair in-pathway if ligand OR receptor in set. **Slice-specific adaptive effect floor (THE density control):** each Δ/aggregating/per-pathway plot keeps pairs with |effect| ≥ the q-quantile WITHIN that exact slice (`--q-stat 0.25` differential on |interaction_stat|; `--q-delta 0.75` baseline on |Δ score|); NOT applied to distribution plots (volcano, rank-rank) or per-group descriptive graphs; `top_n` (30 dotplot/25 ranked) is the hard backstop. **specificity_fdr edge-filter** (≤0.05 in EITHER group) on baseline Δ graphs.

#### 8f. Cross-Tissue Link (`08f_cross_tissue.py`) — THE UNIQUE ANGLE
Six views, reproducible from 8b/8c CSVs. Two biologically aligned arms: E12.5 placenta (Early) → P1/4W/3mo brain (Early); E18.5 placenta (Late) → same (P1 Late flagged `confounded_with_pool`). Views: DEG overlap (hypergeom), RRHO (custom NumPy), pathway concordance, **LR cross-tissue hypotheses (placental ligand × brain receptor, `stress_axis` flag — the publication-quality output)**, TF concordance, ORA. **No cross-tissue CCC** (BBB).

#### 8g. Cross-Age & Persistence (`08g_cross_age.py`) — brain only
Operates on 8b–c output tables. Persistence classes (persistent / resolving_early / established_late / P1_only / transient_4W / emergent_3mo / P1_3mo_only / persistent_directionswap), separately for Early-vs-Relaxed and Late-vs-Relaxed (P1 Late flagged). Cross-arm core signature (view 6) = features persistent in BOTH arms, same direction = paper-quality table. Placenta exits cleanly (incomplete cross-age factorial).

### Phase 9: Cross-Species RRHO2 Validation
Two arms reported separately (see INSTRUCTIONS.md "Phase 9"): ARM A psychiatric/neurodevelopmental; ARM B MS as a stressed-glia signature reference (NOT etiology). Loaders are stubs; smoke-test on Velmeshev first.

### Phase 10: Reproducibility
Run logs, config snapshots, fixed seeds, `uv.lock` + `renv.lock` committed, git for code, `manifest.json` with checksums.

---

## 6. Pipeline Architecture

### Directory Structure (flat `scripts/`, no `src/` package)

```
Analysis/
├── config/        brain.yaml, placenta.yaml (from sample_metadata.csv), dev_split.yaml
├── scripts/       _utils.py + numbered phase scripts + run_*.R workers
├── notebooks/     thin per-phase viewers (load checkpoint, plot inline)
├── data/          raw Cell Ranger h5 (gitignored; symlink to USB-HDD on WS)
├── results/       all outputs (gitignored); {tissue}/{h5ad,plots,tables,validation}/
├── refs/          reference atlases, msigdb_mouse.tsv, celltypist pkls
├── sample_metadata.csv      canonical source of truth for samples
├── run_pipeline_WS.sh       workstation runbook (manual, not changelog)
├── pyproject.toml + uv.lock + .python-version (3.12)
├── setup-remote.sh + scripts/install-r-packages.R + renv.lock
└── Makefile
```

### YAML Configuration

`config/brain.yaml` and `config/placenta.yaml` are generated from `sample_metadata.csv` by `scripts/build_yaml.py` (re-run when the CSV changes). `config/dev_split.yaml` is emitted by `dev_split_h5.py`. Configs include: `tissue`, `group_reference: Relaxed`, `results_dir`, `samples` (with `donor_id`, `h5`, `raw_h5`, `age`, `group`, `sex`, `pool`/`library`), `qc`, `sex_markers`, `random_seed`, plus the integration/scvi/cellbender→soupx blocks, and the declarative Phase-8 blocks below.

**Declarative Phase-8 blocks (shared across 8a–8g):**

```yaml
strata:
  sex: [combined, M, F]          # iter_strata; applied to EVERY contrast

composition:
  min_donors: 2                  # run a stratum if both groups have >=2 donors
  reliable_donors: 3             # <3 in any group -> reliability=low_n

contrasts:
  early_vs_relaxed_per_age:      # PRIMARY
    design: "~ sex + pool + group"
    group_by: age
    test: group
    levels: [Early_Stress, Relaxed]
    flag: primary
  late_vs_relaxed_per_age:       # PRIMARY (P1 pool-confounded -> pool dropped+flagged)
    design: "~ sex + pool + group"
    group_by: age
    test: group
    levels: [Late_Stress, Relaxed]
    flag: primary
  omnibus_3group_per_age:        # PRIMARY (F-test; placenta never has 3 groups/age -> skipped there)
    design: "~ sex + pool + group"
    group_by: age
    test: group_omnibus
    flag: primary
  early_vs_late_per_age:         # SECONDARY (brain only)
    design: "~ sex + pool + group"
    group_by: age
    test: group
    levels: [Early_Stress, Late_Stress]
    flag: secondary
  within_group_across_age:       # DE-only (8b); pool-confounded with age.
                                 # Also the input contrast for 8b follow-ups
                                 # (disruption + shuffle null) — brain only.
    design: "~ sex + age"
    group_by: group
    test: age
    pairwise: [[P1, 4W], [4W, 3mo], [P1, 3mo]]
    flag: confounded_with_pool
  group_x_age_interaction:       # DE-only (8b, DESeq2 LRT); underpowered
    design: "~ sex + pool + group * age"
    test: "group:age"
    flag: underpowered_exploratory
  # NOTE: the old `within_age_sex_stratified` contrast is REMOVED — superseded by
  # the strata mechanism (combined/M/F applied to every contrast).

stress_focused_cell_types: [microglia, oligodendrocyte_lineage, excitatory_neurons, inhibitory_neurons, astrocytes]
```

`config/placenta.yaml`: same structure with `tissue: placenta`, ages E12.5/E18.5, `whole` level only, and only the analyzable contrasts (E12.5 Early-vs-Relaxed, E18.5 Late-vs-Relaxed; no omnibus, no early_vs_late, no cross-age — see §2). 8b follow-ups skip placenta cleanly.

### CLI Usage

```bash
python run.py --config config/brain.yaml --step validate     # MANDATORY first
python run.py --config config/brain.yaml --step qc
python run.py --config config/brain.yaml --step all          # refuses if validate hasn't passed
# Phase-8 scripts run standalone with --n-jobs, e.g.:
uv run python scripts/08a_composition.py --config config/brain.yaml --n-jobs 24
# 8b follow-ups (brain only):
uv run python scripts/08b_developmental_disruption.py --config config/brain.yaml
uv run python scripts/08b_followup_plots.py --config config/brain.yaml
uv run python scripts/08b_disruption_shuffle_test.py --config config/brain.yaml --n-perm 1000 --n-jobs 16
```

---

## 7. Open Questions Before Next Stage

1. Behavioral data on offspring (anxiety, HPA function) — to incorporate alongside transcriptomics?
2. Corticosterone / HPA hormone measurements available?
3. Archived tissue for RNAscope/IHC validation of top hits?

---

## 8. Publication Strategy (Deferred Discussion)

**Realistic IF 12–15 targets:** Nature Communications, Molecular Psychiatry, Biological Psychiatry, Genome Biology. Pure bioinformatics → likely IF 6–10.

**Strongest framing:** cross-tissue developmental cascade — E12.5/E18.5 placenta → P1/4W/3-month brain transcriptional programming under prenatal stress; two parallel arms (Early vs Late) with distinct trajectories. Multi-age × multi-tissue × two-stress-windows design is genuinely unique.

**Headline figures to date (from 8b follow-ups, locked 2026-06-15):**
- Mirror disruption bar + effect-size collapse boxes (`08b_followup_plots.py:plot_disruption`) — descriptive panel
- k-preserving shuffle null + within-stratum 6-bar breakdown (`08b_disruption_shuffle_test.py:plot_shuffle`) — inferential panel

Together these support the framing: *"prenatal stress predominantly disrupts existing developmental programs rather than inducing new ones; when age-DE signal is shared across two groups, Relaxed is almost always one of them."*

**To strengthen:** RNAscope/IHC validation of top 2–3 findings; behavioral validation; external dataset comparison (Phase 9 cross-species).

---

## 9. Remote Workflow & Repo Layout

- **Remote machine** via VPN + SSH from local Mac; VSCode Remote-SSH (host alias `remote-snRNA`).
- **uv + renv** (conda blocked). **tmux** for long jobs. **HTML reports** per phase.
- **WS ↔ Mac mirrored:** code edits local → rsync to WS (`--chmod=Fu+x`); any WS-side edit rsync'd back. WS results saved on Mac under `results_WS/`.
- **Single codebase, two modes:** `scripts/*.py` (real implementation, runnable standalone) + thin `notebooks/*.ipynb` (load saved `.h5ad` checkpoints, inspect inline).

(Full SSH/path/pool details in INSTRUCTIONS.md "Workstation infrastructure".)

---

## 10. Next Steps

**8c production run (current priority):**
1. Run 8c (`uv run python scripts/08c_pathways.py --config config/brain.yaml --tf`) for brain main + 7 focal subclusters. Same for placenta where applicable.
2. Verify GSEA leading-edge tables match the LOST/GAINED direction-class genes from 8b follow-ups (expected: stress-relevant gene sets — GR targets, HPA, neuroinflammation, ER-stress — should show up in the leading-edges of LOST-class genes).

**8e production-complete (2026-06-24):** cell-cell communication compute + all plotting done, both tissues × broad+subtype node schemes. Brain differential null (receptor-side DE too weak); placenta 447 sig LR pairs is the real CCC signal. Centerpiece plots = focal-fan grids + per-pathway LR dotplot/ranked, all with slice-specific adaptive effect floors.

**8d production run:** trajectory, both tissues (still smoke-tested only).

**8f/8g (need both tissues complete through 8c):** cross-tissue cascade + cross-age persistence. 8g brain-only; 8f draws on both tissues' 8b/8c CSVs.

**Phase 9 (after 8):** implement the stub loaders (smoke-test on Velmeshev), run ARM A + ARM B RRHO2 separately.

**Open questions** (resolve on first full pass): behavioral data, corticosterone, archived tissue for validation.

---

## 11. Summary of Key Decisions

| Decision Point | Choice | Rationale |
|---|---|---|
| **Ecosystem** | Python (Scanpy + scvi-tools); R via subprocess | Better batch integration, scales to ~1.5M nuclei |
| **Compute** | RTX 4500 Ada (24 GB) + 258 GB RAM + 56 cores | On-box GPU + ample RAM; ~1.5 days/tissue |
| **Validation gate** | Phase 0 mandatory before any compute | Catches sample swaps, missing metadata, confounds in 5 min |
| **Sex assignment** | `assigned_sex` (Y-chromosome inferred + declared, mismatch flagged) | Detects swaps; needed for E12.5 placenta (sex unknown at sampling) |
| **Donor tracking** | `donor_id` distinct from `sample_id`; **no dam ID** | Pup is the statistical unit; dam random effect cannot be modeled |
| **Group factor** | 3-level: Early_Stress / Late_Stress / Relaxed (reference) | Reflects actual design |
| **Sex schema (Phase 8)** | combined + M + F strata on EVERY contrast (`iter_strata`); `combined` primary, M/F flagged low_n | Declarative; replaces the old `within_age_sex_stratified` contrast |
| **Batch_key for integration** | `pool` (Pool 1–4) | Reflects actual multiplexing |
| **Ambient RNA** | SoupX per sample (R subprocess) | CellBender pickle bug unfixable without Docker |
| **Doublets** | scDblFinder per pool | Doublets form within capture |
| **Integration** | scVI, batch_key=pool, BF16 on GPU | SOTA for complex batch; Ada-optimized |
| **Annotation (brain)** | 4-tier; P1 via scANVI (Rosenberg 2018); 4W/3mo via ABC CellTypist | Di Bella cortex-only mislabeled 42% of P1 as erythrocyte |
| **Annotation (placenta)** | Markers + STAMP vs Liu 2024; `celltype_majority` | No mouse placenta CellTypist model |
| **DE method** | Pseudobulk + PyDESeq2, pup as statistical unit | Cell-level DE is incorrect |
| **DE design** | `~ sex + pool + group` per age | Sex + pool covariates; group is the factor of interest |
| **DE sig thresholds** | `padj<0.05 & |log2FC|>1` | LOCKED across all 8b outputs |
| **Composition tool** | **propeller (R subprocess), NOT scCODA** | scCODA's deps fought scanpy; limma moderation better for small n |
| **Composition denominator** | cells in the slice (region → ÷ donor's region cells; whole → ÷ donor's total) | Region columns answer region-internal proportions |
| **Non-cell-types (Phase 8)** | **DROP contaminants + `unassigned*`** (num + denom), record in diagnostic CSV; not reassigned | Not real cell types — testing them pollutes FDR; plots stay clean |
| **Pool when aliased with group** | **drop `pool` + flag `confounded_with_pool`** (P1 Late; one placenta E12.5 sex stratum) | Rank-deficient otherwise; scVI batch correction doesn't cover count-level tests |
| **Composition reliability** | `min_donors=2` to run; `<3` → `low_n`; trust `ok` + finite effect | low_n NaN/inf rows are degenerate (rare type absent in tiny group) |
| **Contrast specification** | Declarative YAML (`contrasts:` + `strata:`), not hard-coded | Adding contrasts = config edit |
| **Parallelism** | `_utils.parallel_map` + `--n-jobs`, mandatory for repeated subprocess/heavy loops | 8a serial ~3 h → minutes; one-subprocess-per-item loop is a bug |
| **Significance on plots** | outline FDR<0.05 (heatmaps), label sig genes (volcano), thresholds on figure | Reader can name the biology without a side table |
| **Trajectory** | PAGA + diffusion pseudotime; NO velocity, NO CellRank | Flex probe-based — no spliced/unspliced |
| **Cell-cell communication** | LIANA+ consensus, contrast-driven; NO cross-tissue CCC | BBB makes literal cross-tissue signalling implausible (8f view 4 is endocrine) |
| **Cross-tissue** | E12.5→Early brain cascade; E18.5→Late brain cascade | Aligns placenta sampling window with stress timing |
| **Cross-age (8g)** | brain only; persistence per Early/Late arm | Placenta has incomplete cross-age factorial |
| **8e arms** | baseline=descriptive (not a stress test); differential=primary (FDR via `interaction_padj`); per-donor=corroboration (null at n≈4) | Only the differential arm carries animal-respecting stress inference |
| **8e brain differential** | NULL everywhere — reported, not fixed | Receptor-side DE too weak at coarse grouping; confirmed real |
| **8e sex** | combined-only, NOT stratified | n≈2/group per sex → all three arms degenerate |
| **8e readable graph** | focal-fan small-multiples grid (one cell type per panel) | Hairball all-pairs graphs unreadable at 20+ nodes; kept only as supplementary |
| **8e Δ grids** | 4 variants (count/magnitude × baseline/differential), sqrt-scaled, unassigned dropped | sqrt is the one width-scaling that reads across all metric ranges (ratios 13–197) |
| **8e effect floor** | slice-specific adaptive: \|effect\| ≥ q-quantile within each plot's own pairs; q-stat=0.25, q-delta=0.75 | Each plot thins by its own distribution; not a global constant |
| **8e floor scope** | Δ/aggregating/per-pathway plots only; NOT volcano/rank-rank or per-group descriptive graphs | Flooring a distribution plot hides the points that give it meaning |
| **8e LR-pair detail** | per-pathway LR dotplot + ranked lollipop (LR pair is the unit) | Graphs answer topology; the dotplot/ranked answer "which LR pairs changed, between which cells" |
| **8b follow-ups: disruption framing** | NOT "LOST > GAINED" alone (partly maths artifact); the clean claim is "R∩E + R∩L enriched, E∩L depleted across all cell types under k-preserving null" | LOST > GAINED is partly explained by `n_k1 > n_k2` and the AND requirement on GAINED; the within-stratum R∩E/R∩L pattern is robust to those |
| **8b shuffle-test null** | k-preserving (each gene keeps total #sig groups but randomizes WHICH groups) | Symmetric independence null overestimated marginals 3× and falsely showed obs as depleted; k-preserving null correctly captures gene-level overlap structure |
| **8b shuffle-test inferences** | Per-category analytic binomial (`Binom(n_k_stratum, 1/3)`) for the 6 disjoint sig patterns + within-stratum chi-square | Analytic is exact; permutation kept as sanity check on LOST/GAINED diff |
| **8b shuffle-test plot layout** | 2-panel matching disruption plot: Δ mirror bar (Panel A, narrow) + within-stratum 6-bar breakdown (Panel B, wide) | Visual continuity with the descriptive disruption plot; labels OUTSIDE bars to prevent overlap |
| **Environment** | uv + renv (not conda) | Conda channels blocked at firewall |
| **Python↔R bridge** | R as subprocess (not rpy2) | Process isolation; debuggable; no build-against-system-R fragility |
| **Provenance** | `manifest.json` + git hash + uv.lock/renv.lock | Reproducibility wired in from Phase 0 |

---

## 12. Environment & Deployment

**Why uv + renv (not conda):** conda channels blocked at the corporate firewall; PyPI/CRAN/Bioconductor reachable. Python deps via **uv** (`.venv/`, `uv.lock`); R deps via **renv** (`.renv-cache/`, `renv.lock`); Python↔R via subprocess + TSV/JSON. Bootstrap: `./setup-remote.sh` (idempotent; installs uv, `uv sync`, R packages via `scripts/install-r-packages.R`, applies the CellTypist sklearn-1.7 sed patch). cuML from `https://pypi.nvidia.com` for GPU CellTypist training.

**Lock files committed:** `uv.lock`, `renv.lock`, `.python-version` (3.12). The only thing not pinned is system libraries (glibc, CUDA runtime). For absolute reproducibility you would containerize, but Apptainer isn't available and conda is blocked, so uv + renv is the pragmatic best.

**New explicit dep (2026-06-15):** `statsmodels>=0.14` added to `pyproject.toml` for BH correction in `08b_disruption_shuffle_test.py` (the script falls back to `np.nan` if missing but explicit pin guarantees the import path).

**renv Suggests workaround, CellTypist sklearn-1.7 patch, cuML install** — see the matching sections in INSTRUCTIONS.md.

---

## 13. Ambient RNA Correction (SoupX, locked 2026-06-10)

snRNA-seq preps release cytoplasmic + erythroid RNA into the soup. SoupX (R subprocess) replaces CellBender (abandoned 2026-06-05, pickle bug). Workflow: cellranger filtered + raw → `SoupChannel` → `scran::quickCluster` → `autoEstCont` (data-driven rho) → `adjustCounts` → MTX/barcodes/features/JSON → Python assembles per-sample h5ad. Manual rho fallback (`--rho 0.10`). Sanity bounds: `rho_mean` 0.02-0.10 adult / 0.05-0.15 P1; `pct_removed` 2-8% adult / 5-15% P1; >0.30 = over-correction. Full 57-sample run done; summary in `results/{tissue}/tables/02_soupx/02_soupx_summary.csv`. (Note: the P1 erythrocyte-mislabel was a REFERENCE issue, fixed by switching P1 to Rosenberg — not ambient; SoupX correctly stripped Hb.)

`02_qc.py` has a prefer-soupx fallback; SoupX changes counts → Phase 2 onward re-run from corrected counts.
