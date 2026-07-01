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

**All phases through 8g implemented and run on the workstation; 8b main DE done; 8b follow-ups (disruption + shuffle null) done for brain main + 3 focal subclusters; 8c–8d smoke-tested; 8e PRODUCTION-COMPLETE; 8f PRODUCTION-COMPLETE; 8g PRODUCTION-COMPLETE.** **Phase 9 (cross-species) PLACENTA ARM PRODUCTION-COMPLETE (2026-06-25): two independent human-placenta validations (Gunter-Rahman obesity + Admati PE), full RRHO → permutation null → concordant GSEA → leading-edge → plots.** **Phase 9 BRAIN ARM COMPUTE-COMPLETE (2026-06-29): four independent human validations (Velmeshev ASD, Maitra MDD-female, Nagy MDD-male, Macnair MS) through the same engine (`h10*`) — prep + RRHO + permutation null + concordant GSEA + leading-edge + TF + diagnostics + comprehensive plots. Headline: mouse prenatal-stress neuronal programs converge strongly + directionally on human MDD cortex (both sexes); the 8f/8g IFN/immune thread is recovered as microglial co-suppression in MDD/ASD and inverts (up) in MS.** **Phase 9 BRAIN ARM 5 (Hwang/Girgenti 2025 PTSD/MDD dlPFC, `h10f`) IN PROGRESS (2026-07-01): a fifth, on-target *trauma* arm being built by reclustering the deposited counts (Zenodo 15186498, open CC-BY) through the SAME pipeline as the other four brain arms and the mouse data — own scVI+Leiden+marker annotation, own PyDESeq2 rankings — for cross-Fig-4 consistency; the authors' Supplementary Tables are retained purely as a validation/cross-check layer, never as pipeline input. Motivating test: PTSD microglia are suppressed (SPP1 down) while MDD microglia are up, so PTSD should be concordant with the mouse microglial IFN/immune down-thread and MDD divergent — a PTSD-vs-MDD internal directional control within one dataset.** Cross-tissue (8f) and cross-age (8g) views operate on completed 8b/8c CSVs so they need no re-runs.

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
| 8a Composition (propeller) | ✓ done, both tissues | `08a_composition.py` + `run_propeller.R` |
| 8b Pseudobulk DE (PyDESeq2) | ✓ done — brain main + 7 focal subclusters | `08b_de.py` |
| 8b Summary plots | ✓ done | `08b_de_summary.py` |
| 8b Developmental disruption | ✓ done (brain main + 3 focal subclusters) | `08b_developmental_disruption.py` |
| 8b Follow-up plots | ✓ done | `08b_followup_plots.py` |
| 8b Disruption shuffle test | ✓ done | `08b_disruption_shuffle_test.py` |
| 8c GSEA + leading-edge + TF activity | ✓ smoke-tested (decoupler ULM + CollecTRI) | `08c_pathways.py` + `fetch_genesets.R` |
| 8d Trajectory (PAGA + DPT) | ✓ smoke-tested | `08d_trajectory.py` |
| 8e Cell-cell communication (compute) | ✓ DONE — both tissues, 3 arms, n_perms=1000 | `08e_comms.py` + `_08e_*` workers |
| 8e Communication (plotting) | ✓ DONE — 7 plot families | `08e_comms_summary.py` + `_08e_plots_*.py` |
| 8f Cross-tissue | ✓ PRODUCTION-COMPLETE — six views | `08f_cross_tissue.py` |
| 8g Cross-age / persistence | ✓ PRODUCTION-COMPLETE — comprehensive + B/C/View-7 | `08g_cross_age.py` |
| **9 Cross-species — PLACENTA, Gunter-Rahman (obesity)** | **✓ PRODUCTION-COMPLETE (2026-06-25)** | `h09a`–`h09h`, `h09_summary_plots.py` |
| **9 Cross-species — PLACENTA, Admati (PE 2×2)** | **✓ PRODUCTION-COMPLETE (2026-06-25)** | `h09j`, `h09k`, `h09k_diagnostics.py`, `h09k_plots.py`, `h09k_rrho_maps.py` |
| **9 Cross-species — BRAIN (4 datasets)** | **✓ COMPUTE-COMPLETE (2026-06-29)** | `h10a`–`h10e`, `h10b_brain_rrho.py` (engine), `h10b_diagnostics.py`, `h10b_rrho_maps.py`, `h10_summary_plots.py` |
| **9 Cross-species — BRAIN ARM 5 (Hwang/Girgenti PTSD/MDD, `h10f`)** | **⏳ IN PROGRESS (2026-07-01) — recluster + own DE, then RRHO via `h10b`** | `h10f_prep_hwang.py` (+ reuse `h10b` engine, `h10f_validate_vs_tables.py`) |

**Naming convention (locked 2026-06-25):** human-side scripts are `h09X_...` (e.g. `h09a_prep_human_placenta.py`); R workers are `h_run_*.R` / `h_fetch_*.R` (e.g. `h_run_soupx_from_raw.R`, `h_fetch_genesets.R`). This keeps the human cross-species arm visually separate from the numbered mouse pipeline.

**Key implementation notes (selected; Phase 9 additions at the end; full history in git):**

1. **Flat scripts layout, no `src/snrna/` package.** Shared helpers in `scripts/_utils.py`. Phase scripts are numbered standalone files.
2. **Raw counts in `.X`, lognorm not persisted.** Recompute via `_utils.add_lognorm(adata)`.
3. **QC: per-sample MAD + hard floors + hard caps + cohort-outlier flag.**
4. **scDblFinder per pool**, `samples=` arg.
5. **HVG exclusion lists.** mito/ribo/hemoglobin/sex-linked always; placenta adds Prl*/Psg*/Cgb*/Cga.
6. **scVI uses BF16 mixed precision on GPU.**
7. **Leiden uses igraph backend.**
8. **Phase 7 = per-cluster majority voting.** P1 brain via scANVI from Rosenberg 2018. 4W/3mo via ABC CellTypist.
9. **scCODA abandoned for composition (8a). propeller (speckle+limma) via R subprocess.**
10. **Phase 1 = SoupX via R subprocess (locked 2026-06-10).**
11. **Brain CellTypist models retrained on GPU via cuML (locked 2026-06-10).**
12. **Phase 8a finalized (2026-06-15):** full sex-strata schema, parallelized; contaminants + `unassigned*` dropped; this is the template for 8b–8g.
13. **Phase 8b follow-ups locked (2026-06-15) — disruption + shuffle null.**
14. **Phase 8e production-complete (compute + plotting, locked 2026-06-24).** placenta differential = 447 sig LR pairs; brain differential = NULL.
15. **Phase 8f production-complete (2026-06-25).** Six views. Headline: cross-tissue IFN/complement/cytokine co-suppression + ECM-flavoured LR table.
16. **Phase 8g production-complete (2026-06-25).** 0 persistent genes; pathways persist; IFN perinatal-transient; ECM/mesenchymal durable; gliogenesis cross-arm core.
17. **Phase 9 placenta cross-species arm production-complete (2026-06-25).** Two independent human-placenta validations (Gunter-Rahman obesity, Admati PE). Full method, key findings, and per-script detail in §5 "Phase 9" below. Core engineering: a mouse↔human 1:1 ortholog bridge (`refs/mouse_human_orthologs.tsv`, 16,030 pairs via pybiomart); compartment-level pseudobulk RRHO (mouse/human trophoblast subtypes have no 1:1 homology, so the **compartment** is the bridge level); all RRHO/GSEA/leading-edge functions lifted **verbatim** from the mouse arm (`08f`/`08c`) so the human and mouse arms are methodologically identical. Human MSigDB (`refs/msigdb_human.tsv`, 9,427 sets H/C2:REACTOME/C5:GOBP) via `h_fetch_genesets.R`.

---

## 2. Critical Considerations & Caveats

### Statistical Power Limitations

**n=2 per sex per condition per age is a real limitation:**
- Cannot reliably test sex × condition interactions
- Pseudobulk DE with n=2 vs n=2 detects only large effect sizes (logFC > 1.5–2)
- Single-cell-level DE treating cells as replicates is **incorrect**

**Decision: combined analysis primary, sex strata run systematically but flagged**
- **Primary (`combined` stratum):** pool sexes per group, with sex as a covariate. Design: **`~ sex + pool + group`**.
- **Sex strata (`M`, `F`):** flagged `low_n` / `underpowered_exploratory`.

**No dam ID — litter random effect cannot be modeled** — each pup treated as independent (anti-conservative; explicit methods caveat).

### 10x Flex-Specific Considerations
1. Probe-based capture (not poly-A) — biased toward exonic sequences
2. RNA velocity not feasible — Trajectory = PAGA + diffusion pseudotime only
3. Ambient RNA from probe leakage — handled by SoupX
4. Doublet rates comparable to standard 3'

### snRNA-seq-Specific Considerations
1. High intronic read fraction is normal
2. Mitochondrial % near zero in good nuclei prep
3. Ambient RNA severe (placenta hemoglobin; brain Malat1/mito) — SoupX-corrected
4. scRNA-built references may not transfer perfectly

### Cell-Cell Communication on snRNA-seq
Use **LIANA+** (consensus). No cross-tissue CCC (BBB) — placenta→brain link is the LR-from-DE endocrine framing in 8f view 4.

### Pool/Library Confounding — REAL AND CONSEQUENTIAL
- **Cleanly analyzable:** 4-week brain (all 3 groups in Pool1+Pool2), 3-month brain (single pool), P1 brain Early-vs-Relaxed (Pool2), E18.5 placenta Late-vs-Relaxed (Pool4).
- **Confounded/flagged (`confounded_with_pool`):** P1 brain Late Stress (Pool3-only); brain age trajectories; Pool3 mixes tissues; placenta cross-age (not comparable).
- **Mitigations:** `pool` covariate where estimable, drop+flag where aliased; `batch_key=pool` for scVI; k-preserving shuffle null preserves pool-driven structure.

---

## 3. Compute Environment

**Available machine:** 258 GB RAM, 56 CPU cores, 1× NVIDIA RTX 4500 Ada (24 GB VRAM), GPU on the same box.

**Implications:** SoupX parallel via `parallel_map`; scVI on GPU BF16 ~2-3 hrs/tissue; propeller/R loops parallelized; shuffle tests CPU-bound numpy. Total ~1-1.5 days/tissue unattended. **Phase 9 human prep** (e.g. streaming the 5.9 GB Admati sc matrix, 16-cell 2×2 RRHO with 5k-shuffle nulls) all fits comfortably; the heavy outputs are small CSV/parquet so **all Phase-9 plotting is trivial and can run on the Mac offline**.

---

## 4. Analysis Strategy

### Ecosystem: Python (Scanpy + scvi-tools), R via subprocess
scVI/scANVI for batch structure; PyDESeq2 pseudobulk DE; scales to ~1.5M nuclei. R as subprocess for scDblFinder, propeller, SoupX, fgsea, SingleR.

### Joint Integration, Per-Age Analysis
Integration/clustering/annotation/UMAP on the joint object; composition/DE/communication per age; trajectory and cross-age cross-cutting.

### Reference Atlases
- **Brain P1:** Rosenberg 2018 (scANVI). **4W/3mo:** ABC CellTypist.
- **Placenta:** curated markers + STAMP vs Liu 2024 (`celltype_majority`).
- **Phase 9 human placenta references:** Vento-Tormo 2018 (first-trimester, SingleR corroboration for Gunter-Rahman annotation); authors' own labels used directly for Admati.

---

## 5. Detailed Analysis Pipeline

(Phases 0–7 as implemented; see §1b for status. Phase 0 validation is a mandatory gate.)

### Phase 7e: Cell-type Counts Diagnostic — DONE
Per-donor × cell-type count CSV used to sanity-check 8a propeller inputs.

### Phase 8: Downstream Biology
**Architectural principle:** All contrasts declarative in YAML (`contrasts:` + `strata:`). Reference level for `group`: Relaxed.

#### 8a. Composition — DONE
propeller (speckle+limma) via R subprocess; per-donor counts; sex strata; contaminants+unassigned dropped; pool dropped+flagged where aliased.

#### 8b. Differential Expression — DONE
Pseudobulk only; animal = unit; `~ sex + pool + group`; `padj<0.05 & |log2FC|>1`. Brain main: 20.3M rows, 763K sig.

#### 8b follow-ups — DONE
Developmental disruption (5 direction classes; class column literally named **`direction`**) + shuffle null. Headline: when age-DE is shared across two groups, Relaxed is almost always one of them.

#### 8c. Pathway Analysis
fgsea-multilevel via `run_fgsea.R`; decoupler ULM + CollecTRI TF activity; leading-edge table. **These functions are the basis of the Phase-9 human GSEA/TF — lifted verbatim.**

#### 8d. Trajectory — Brain primarily
PAGA + diffusion pseudotime.

#### 8e. Cell-Cell Communication — PRODUCTION-COMPLETE
LIANA-py; three arms; placenta differential 447 sig, brain differential NULL.

#### 8f. Cross-Tissue Link — PRODUCTION-COMPLETE
Six views; LR table is the publication-quality output; IFN/complement/cytokine co-suppression concordant placenta→brain.

#### 8g. Cross-Age & Persistence — brain only — PRODUCTION-COMPLETE
0 persistent genes; pathways persist; IFN perinatal-transient; ECM/mesenchymal durable; gliogenesis cross-arm core.

**Unified 8f+8g spine:** (1) ECM/mesenchymal — cross-tissue + durable; (2) IFN/immune — cross-tissue but perinatal-transient; (3) gliogenesis — durable, regimen-independent.

---

### Phase 9: Cross-Species Validation — PLACENTA ARM PRODUCTION-COMPLETE (2026-06-25)

Two **independent** human-placenta validations of the mouse prenatal-stress placental signature, run as separate arms and reported separately. Both bridge to mouse via a 1:1 ortholog table and compare at the **compartment** level (trophoblast, decidua_stromal, vascular, immune, [erythroid]) because mouse and human trophoblast subtypes have no 1:1 homology. The mouse anchor for both is **E18.5 Late-Stress vs Relaxed** (and, for Admati, also E12.5 Early-Stress vs Relaxed). RRHO2-style pseudobulk concordance is the primary method; ranking metric = signed DESeq2 Wald `stat` (matches 8f/8c).

**Framing (locked):** "convergent cell-type stress programs," NOT etiologic causal claims. Honest nulls reported as-is.

#### ARM 1 — Gunter-Rahman 2025 (maternal obesity), GSE271976 — COMPLETE
Human **term** placenta snRNA-seq, maternal obese vs lean (20 samples: 10 maternal-facing, 10 fetal-facing), raw `_raw_feature_bc_matrix.h5` only.

Pipeline (`h09a`→`h09h`, then `h09_summary_plots.py`):
- **h09a** `h09a_prep_human_placenta.py` + `h_run_soupx_from_raw.R` — knee/inflection cell-call (DropletUtils::barcodeRanks, `--cutoff inflection`) → SoupX → per-sample h5ad. Bug fixed: `read10xCounts` returns HDF5 DelayedMatrix → coerce `as(m,"CsparseMatrix")` before SoupX. 78,103 cells.
- **h09b** QC + scDblFinder (reuses mouse `run_scdblfinder.R`) → 66,289 singlets (~5% of paper's 62,864 — validates).
- **h09c** HVG (seurat_v3) → scVI (batch=sample, continuous covariate pct_mt, BF16) → Leiden, 24 clusters. Bug: placental hormone genes (CGA/CGB/CSH/PSG) break seurat_v3 loess → EXCLUDE before HVG.
- **h09d** `h09d_annotate.py` + `h_run_singler.R` + `config/human_placenta_markers.yaml` — marker-majority (paper-exact subtypes) → `subtype`+`compartment`; SingleR (Vento-Tormo ref) corroboration, parallelized `--n-jobs 24`. 21/24 clusters agree.
- **h09e** `h09e_cross_species_rrho.py` — compartment-level RRHO, mouse E18.5 Late-vs-Relaxed ↔ human obese-vs-lean. Mouse DE recomputed fresh (`~sex+group`), human (`~side+condition`). RRHO funcs copied verbatim from `08f`. Ortholog bridge built by `h09e_build_ortholog_map.py` (pybiomart, 16,030 1:1). Compartment map in `config/cross_species_celltype_map.yaml` under `placenta_compartments:`.
- **h09f** `h09f_overlap_null.py` — concordant-up leading-edge genes + 10,000-shuffle permutation null (parallel, use_threads=False).
- **h09g** `h09g_pathways_tf.py` + `refs/msigdb_human.tsv` (built by `h_fetch_genesets.R`) — two single-species fgsea → intersect FDR<0.05 same-sign; TF via decoupler ULM + CollecTRI human. Functions lifted verbatim from `08c`.
- **h09h** `h09h_leading_edge.py` — per-species leading-edge + overlap per concordant pathway (`compute_leading_edge` verbatim from `08c`).
- **h09_summary_plots.py** — CSV-only: RRHO grid, NES concordance scatter, pathway dotplots, hypoxia leading-edge. adjustText for label de-overlap.

**Findings (Gunter-Rahman):**
- **RRHO concordance, permutation-backed:** decidua_stromal (peak 8.0, empirical p=1e-4), vascular (5.4, p=3e-4), trophoblast (5.3, p=1e-3) all significant; **immune weak and correctly null** (peak 2.3, p=0.23). RRHO maps: decidua = cleanest up-up hotspot; vascular = bidirectional; trophoblast = moderate bidirectional; immune = diffuse.
- **Concordance lives in the tails** (global Spearman ≤0.03) — RRHO is designed for exactly this; don't claim genome-wide correlation.
- **Concordant pathways headline = HALLMARK_HYPOXIA concordantly UP in trophoblast and immune** — independently recovers Gunter-Rahman's EVT-hypoxia thesis from the mouse side. decidua = coordinated cell-cycle/proliferation DOWN (E2F/G2M/MYC/DNA-replication); vascular = ribosome-biogenesis/rRNA DOWN + peptide-hormone UP.
- **Gene-level conservation:** leading-edge Jaccard ~0.36–0.43 — the **same genes** drive the pathways. Shared trophoblast hypoxia genes: SLC2A1, PGK1, PDK1, NDRG1, BNIP3L, ERRFI1, FBP1, DUSP1, FOS, ATF3, etc. immune hypoxia (64 shared) incl. VEGFA, LDHA, CXCR4.
- **0 concordant TFs** all compartments — honest null (ULM underpowered on a single n=1 ranking vector).

#### ARM 2 — Admati 2023 (preeclampsia), figshare 23264102 — COMPLETE (2×2)
Human placenta, eoPE + loPE + GA-matched controls, both scRNA and snRNA. **Downloaded via browser User-Agent curl to bypass the figshare AWS-WAF challenge** (`x-amzn-waf-action: challenge`); plain curl/wget returns an empty 202.

Files: `sc_admati.zip` (29 filtered Cell Ranger mtx triplets + `PE_samples_metadata.xlsx`, NO raw/annotation); `sn_PE_TB_allcells_with_metadata` (snRNA trophoblast-only, 6,862 cells); `sc_PE_allcells_with_metadata` (figshare file 41003240, 5.9 GB, **all compartments, author-annotated** — the substantive file).

- **sn trophoblast arm (`h09i`) — DEFERRED.** 3 PE vs 2 control donors → underpowered. Designed as a targeted confirmation (sign-test + GSEA of the conserved hypoxia genes in eoPE) but crashed on duplicate gene symbols (`var_names_make_unique` needed before DESeq2). Set aside in favour of the powered sc arm.
- **sc all-compartments 2×2 arm — the real PE validation.**
  - **h09j** `h09j_prep_admati_sc.py` — STREAM-parses the 5.9 GB transposed file (23 metadata rows then genes; ~98k cells in columns, too wide for a whole-frame read) in gene-row chunks, accumulating pseudobulk via a sparse (cells × groups) indicator. Compartment from celltype prefix (`TB_`→trophoblast, `STROMAL_`→decidua_stromal, `VASCULAR_`→vascular, `IMMUNE_`→immune). Output: per-(donor×compartment) pseudobulk parquet + group meta. **Powering: 10 eoPE / 3 early_control / 7 loPE / 6 late_control donors per compartment.** No erythroid (already dropped from the mouse RRHO). No SoupX (filtered matrices only, no raw) — flagged.
  - **h09k** `h09k_admati_2x2.py` — the 2×2: mouse {E12.5 Early-vs-Relaxed, E18.5 Late-vs-Relaxed} × human {eoPE vs early_control, loPE vs late_control} × 4 compartments. Per cell: RRHO + permutation null + concordant GSEA + leading-edge (all machinery imported from h09e/h09g/h09h). Mouse E12.5 ranking computed fresh (`~sex+group`). Saves `h09k_rankings.parquet` (so replots/RRHO-maps never recompute DE), `h09k_rrho_2x2_summary.csv`, `h09k_concordant_pathways_2x2.csv` (with FDR), `h09k_leading_edge_2x2.csv`.
  - **h09k_diagnostics.py** — interrogates two anomalies before interpreting (see findings).
  - **h09k_plots.py** — CSV-only: 2×2 peak grid; **diverging dotplot** (the two conserved axes; dot size = −log10 FDR); per-subtype named pathway dotplots; named gene panels.
  - **h09k_rrho_maps.py** — the per-compartment 2×2 RRHO rank-rank heatmaps (reads the saved rankings).

**Findings (Admati PE 2×2):**
- **GA-matched diagonal hypothesis REJECTED.** Diagonal mean concordance peak (9.41) is *not* higher than off-diagonal (10.62). Confirmed by diagnostics: mouse E12.5 vs E18.5 trophoblast rankings are positively correlated (r=+0.28), and eoPE/loPE stat distributions are comparable. Structure is organized by **PE subtype × biology**, not by gestational stage. (The stressor-timing confound on the mouse axis was accepted going in; the result moots it.)
- **Two distinct conserved axes — the headline:**
  - **eoPE → HYPOXIA** (concordant UP, broad across all 4 compartments). Shared trophoblast hypoxia leading-edge genes overlap the Gunter-Rahman set: NDRG1, BNIP3L, ERRFI1, PLIN2, ANGPTL4, DDIT4, PFKFB3, ADM, CITED2, PLAUR, plus PAM, TIPARP, TNFAIP3, ACKR3, INHA, KLF6. eoPE's top leading-edge pathways: HALLMARK_HYPOXIA, RESPONSE_TO_OXYGEN, TNFA-NFKB, GLYCOLYSIS. → eoPE is the third independent human stressor recovering the conserved hypoxia program (mouse + obesity + eoPE).
  - **loPE → OXIDATIVE PHOSPHORYLATION / electron-transport SUPPRESSION** (concordant DOWN, broad) + insulin/peptide-hormone response. Shared trophoblast OXPHOS genes (all down in both): NDUFB4, CHCHD2, NDUFS5/B5/B6/B1/C1/C2/V2/AB1, COX20/14/6A1/6B1/6C/7C, UQCRQ, LDHB, MPC2. loPE's top leading-edge pathways: PEPTIDE_HORMONE, INSULIN, RESPIRATORY_ELECTRON_TRANSPORT, AEROBIC_RESPIRATION, TRANSLATION. → matches loPE's maternal/metabolic pathophysiology.
- **loPE RRHO peaks far exceed eoPE's** (e.g. trophoblast loPE peak 24–37 vs eoPE 5–7). Diagnostics show this is **partly power** (control-n: downsampling loPE controls 6→3 drops the peak ~30%, to 17–26) but **mostly real** (even matched at 3 controls loPE stays 17–26 ≫ eoPE 5–7). Interpretation: the OXPHOS/translation program loPE shares is a *large, coherent* gene set → stronger RRHO tail overlap; report peak magnitudes with the n-sensitivity caveat (compare within column, not across).
- **`E18.5 × loPE × trophoblast` is classified `discordant`** but this is a **tail artifact** — the global Spearman is positive (+0.09) and both mouse ages correlate positively with loPE; one strong off-diagonal quadrant flips the label. Do NOT headline the discordant cell.

#### Datasets catalogued for Phase 9 (placenta + future brain)
- **Gunter-Rahman GSE271976** (term, obesity, snRNA) — used.
- **Admati 2023 figshare 23264102** (eoPE/loPE, sc + sn) — sc all-cells file 41003240 used; WAF bypass required.
- **Vento-Tormo 2018** (first-trimester reference, CELLxGENE) — SingleR ref.
- **Marsh 2022 GSE198373** (mid-gestation) — downloaded, DROPPED from Phase 9 (no stress contrast; E12.5 anchor only).
- **ECHO-PATHWAYS** (dbGaP phs003619 CANDLE + phs003620 GAPPS) — only measured-psychological-stress placenta data; BULK; controlled-access (2–6 mo). PI action item; checklist at `refs/dbgap_application_checklist.md`. Revision-stage upgrade.
- **Brain (PROCESSED — see brain section below):** Velmeshev 2019, Maitra 2023 GSE213982, Nagy 2020 GSE144136, Macnair 2025 (Zenodo 8338963). **Hwang/Girgenti 2025 PTSD/MDD dlPFC (Zenodo 15186498, open CC-BY) — ADDED 2026-07-01 as brain ARM 5 (`h10f`, in progress; reclustering).** Herring 2022 GSE168408 deferred (age-anchor, raw, 24GB tar).

### Phase 9: Cross-Species Validation — BRAIN ARM (h10*) — COMPUTE-COMPLETE (2026-06-29)

Tests whether the mouse prenatal-stress brain signature (8f/8g threads: ECM/mesenchymal, IFN/immune, gliogenesis) shows convergent gene-expression/pathway changes in human psychiatric/neurodevelopmental/demyelinating cortex. **Similarity, NOT etiology** — no human prenatal-stress paradigm exists.

**Four independent datasets, same engine (`h10*`), reported separately:**

| Dataset | Disorder | Region | Donors | Unit | Modality | Condition source |
|---|---|---|---|---|---|---|
| Velmeshev 2019 (UCSC) | ASD (peds 4–22y) | PFC+ACC | 31 | sample (indiv×region) | scRNA | `meta.tsv` `diagnosis` |
| Maitra 2023 (GSE213982, F) | MDD female | dlPFC BA9 | 18/20 | donor | snRNA | `maitra_donor_meta.csv` |
| Nagy 2020 (GSE144136) | MDD male | dlPFC BA9 | 17/17 | donor | scRNA | barcode-inline |
| Macnair 2025 (Zenodo 8338963) | MS (stressed-glia ref) | GM+WM | 26/54 | donor | snRNA | `col_data` `diagnosis` |

(Nagy is the male MDD arm standalone; the GSE213982 combined matrix folds Nagy males in as M1–M34 but the M#→Nagy crosswalk is unrecoverable — Nagy standalone keeps condition inline.) Herring 2022 deferred (age-anchor only).

**Method (engine identical to placenta h09, functions lifted VERBATIM from h09e/h09g/h09h/h09k):**
- **Bridge = broad 7-class** {ExN, InN, Ast, Oli, OPC, Mic, Endo} (brain types ARE homologous, unlike trophoblast). Finer human subtypes = exploratory secondary (Velmeshev sensitivity variant only).
- **Two tiers:** T1 whole-mouse-brain × all 7; T2 mouse `Isocortex` (ABC region token) neurons × human cortical neurons. Isocortex carries ONLY ExN/InN by construction → T2 = neurons; glia T1-only (the microglial IFN thread can't get cortical restriction — stated).
- **Full 3×2 mouse grid, no anchor pre-selection:** 3 ages (P1/4W/3mo) × {early_vs_relaxed, late_vs_relaxed} = 6 rankings/celltype. NO diagonal hypothesis; age-matched cells highlighted not pre-collapsed.
- **Mouse rankings = reuse 08b Wald `stat`** (no recompute, identical to paper Figs 2/3), bridged to human symbols via `mouse_human_orthologs.tsv`. **Human rankings = PyDESeq2 from the h10a/c/d/e pseudobulk parquet** (`~ [sex+] diagnosis`; sex auto-drops when constant).
- **RRHO → 5000-shuffle permutation null → concordant GSEA (two single-species fgsea, intersect FDR<0.05 same-sign) → leading-edge → TF (CollecTRI human ULM).** `refs/msigdb_human.tsv`.
- **`robust_class` guard (NEW):** the argmax label flips on noise; trust the directional label only if (margin ≥25% of peak AND Spearman sign agrees: concordant→r>0, discordant→r<0) else `ambiguous`. Shuffle null runs on EVERY cell — magnitude always tested, direction withheld when fragile. The RRHO MAP is ground truth.

**Diagnostics (`h10b_diagnostics.py`) MANDATORY before interpretation** — peak-vs-vector-strength, quadrant decomposition, age-specificity. Established the findings are real, not mechanical.

**Findings (diagnostics-confirmed):**
- **Neurons dominate every dataset by RRHO peak** (ExN/InN top-2 human celltypes in all four).
- **MDD (both sexes) = headline: strong, directional, neuron-specific DOWN-concordance at mouse 4W.** Maitra `4W ExN→ExN` peak 232 (concordant_down, margin 17, Spearman +0.33, p=2e-4); `4W InN→InN` peak 213 (margin 56, r +0.29) — the strongest, most coherent cross-species cells in the whole project. Nagy replicates (4W ExN 105, InN 76, Oli 61). 4W wins ONLY neuron rows (glia peak 3mo) → celltype+age-specific, NOT mechanical.
- **Sex nuance:** female-MDD (Maitra) adds Mic (peak 26); male-MDD (Nagy) adds Oli (61).
- **ASD (Velmeshev): weaker, bidirectional, perinatal.** P1 ExN/InN peaks ~49 but ALL `ambiguous` (near-tied quadrants, flat Spearman) — significant magnitude, untrustworthy direction. Genuinely different/weaker than MDD.
- **MS (Macnair): weak, glia-leaning** by peak. Distinct from psychiatric arms — stressed-glia framing holds.
- **TFs did NOT null** (unlike placenta): 37/156/137/47 concordant TFs (Velm/Mait/Nagy/Macn); many recur across grid cells.

**The 8f/8g THREADS recovered in human (pathway-keyed lens — see plotting note):**
- **IFN/immune co-suppression is CONSERVED and MICROGLIAL.** Maitra/Nagy/Velmeshev microglia: 40–48 concordant immune pathways each, almost entirely `down_both` — the exact 8g signature (INTERFERON_ALPHA/GAMMA_RESPONSE, INTERFERON_SIGNALING, TNFA_SIGNALING_VIA_NFKB), NES strongly negative both species, P1-prominent (matches 8g "IFN perinatal-transient"). **This is the 8f/8g IFN/immune thread recovered in human MDD + ASD cortex.**
- **MS INVERTS the IFN direction — a clean directional control.** Macnair microglia show only 4 concordant immune rows, all `up_both`. Direct check of the MS microglial ranking: canonical IFN/inflammatory genes are strongly UP in MS (median Wald +1.05; IRF7 +2.7, NFKB1 +3.9, BST2 +4.2, C1QA +2.9 — classic neuroinflammation), while mouse prenatal-stress microglia go DOWN (median −1.15). Opposite directions → legitimately non-concordant. **MS microglia ARE inflamed — just opposite to prenatal stress.** Confirms stress-suppresses / MS-activates; a mechanistic validation of the "MS ≠ etiology" framing.
- **ECM/mesenchymal (164 concordant rows) and gliogenesis (34 rows) threads also present** — EPITHELIAL_MESENCHYMAL_TRANSITION, MESENCHYMAL_CELL_DIFFERENTIATION, COLLAGEN_FORMATION; GLIOGENESIS, GLIAL_CELL/ASTROCYTE_DIFFERENTIATION. The full 8f/8g spine is recoverable.

**Key plotting lesson (locked):** peak-keyed views (overview/master/maps) are dominated by neurons (bulk gene overlap) and STRUCTURALLY HIDE small-but-coherent programs like microglial IFN (small peak, tight pathway). Pathway-keyed views (thread-scanner, Mic/IFN) gate on GSEA FDR — NOT RRHO peak — and surface the threads. **Both lenses are kept; they answer different questions** ("where is bulk overlap" vs "which coherent programs are conserved").

**Interpretation (locked):** *Mouse prenatal-stress neuronal programs converge most strongly on human MDD cortex (directional, both sexes, mouse 4W signature); the IFN/immune thread is recovered as microglial co-suppression in MDD/ASD and inverts in MS; ECM and gliogenesis threads also present — with a weaker/bidirectional ASD signature and a distinct/inverted MS signature establishing disorder-specificity.* The Fig 4 brain anchor.

**Scripts:** `h10a_prep_velmeshev.py` (transposed TSV, ENSG|SYM split, primary+sensitivity), `h10c/d/e_prep_*` (mtx triplets → pseudobulk parquet), `h10b_brain_rrho.py` (engine: `--dataset --variant --n-perm --n-jobs --tf`), `h10b_diagnostics.py`, `h10b_rrho_maps.py`, `h10_summary_plots.py` (comprehensive: peak-keyed overview/master/maps + pathway-keyed thread-scanner/Mic-IFN; filter `empirical_p<0.05` for peak views, GSEA-FDR for pathway views; CSV/parquet-only, Mac-runnable).

**Outputs:** `data/human_validation/brain/<ds>/{tables,plots}/` + `_synthesis/plots/` (01_overview, 02_master, 03_thread_scanner, 04_microglia_ifn, 05_ifn_all_celltypes).

**Caveats:** Velmeshev/Nagy scRNA vs mouse/Maitra/Macnair snRNA (modality); RRHO peaks n-sensitive (compare within column); Macnair no SoupX + no GM/WM bridge (mouse whole-brain, unit=donor); "age-matched" is a soft guess (strongest MDD signal lands at mouse 4W regardless of human age); no dam ID (8b caveat carries).

> **OPEN:** the placenta arm (h09e/h09k) used margin-free argmax labels; findings stand (only the known E18.5×loPE×trophoblast cell was fragile, already handled) but the `robust_class` margin guard should be backported at revision for Fig-4-wide consistency.

### Phase 9: Cross-Species Validation — BRAIN ARM 5 (Hwang/Girgenti PTSD/MDD, `h10f`) — IN PROGRESS (2026-07-01)

A fifth brain arm and the single most on-target *stress* dataset available: **Hwang/Girgenti et al. 2025 *Nature*, "Single-cell transcriptomic and chromatin dynamics of the human brain in PTSD"** — snRNA-seq dlPFC, **935,371 discovery nuclei × 27,982 genes, 105 donors (39 CON / 36 MDD / 36 PTSD)**, single region. PTSD is the closest adult analog to prenatal stress, and the paper reports **microglial pro-inflammatory SPP1 signalling *decreased* in PTSD but *increased* in MDD** — i.e. neuroimmune/microglial *suppression* specific to the PTSD brain. Since the mouse prenatal-stress microglia go DOWN on the IFN/immune thread, the prediction is clean: **PTSD concordant, MDD divergent — a PTSD-vs-MDD internal directional control within one dataset/pipeline** (tighter than the cross-dataset MS inversion).

**Data provenance & the deposit gotcha (locked).** Zenodo 15186498 (open, CC-BY-4.0). Files: `RNA_count_mat.npz` (**scipy sparse COO, uint32, cells×genes, nnz≈2.85B — NOT dense**, despite the repo's `get_RNA_count_matrix.py` `toarray`+`savez` script which did not generate the deposit), `RNA_cellnames.txt` (`{MS####XX channel}-{16bp barcode}`, 105 channels), `RNA_genenames.txt` (9-col CSV, col1 = ENSG). **The per-cell obs (celltype `anno`, `Condition`, `Channel`) was NOT deposited** — it lives in the undeposited `RNA_FINAL.zarr`. So the snRNA download is counts + barcodes + genes only; diagnosis and celltype must be recovered externally.

**Decision (locked 2026-07-01): recluster from scratch, do NOT shortcut off the tables.** Build `h10f` through the **same pipeline as the other four brain arms and the mouse data** — own scVI (batch = channel, BF16) + Leiden + marker annotation (their `lake_genes` panel → broad-7), own PyDESeq2 `~ sex + Condition` rankings, then the `h10b` RRHO engine. Consistency across Fig 4 outranks the (available) shortcut of ranking directly off the authors' supplementary DE tables. Diagnosis is recovered from **Supplementary Table 1** (channel `MS####XX` → Condition/Sex; 105 channels all join cleanly). ENSG→symbol via `refs/mouse_human_orthologs.tsv` `human_ensembl` (no extra pybiomart call). **No SoupX** (filtered counts only, no raw — flagged, like Admati). Two human contrasts (PTSD-vs-CON, MDD-vs-CON) run against the same mouse 08b rankings; `robust_class` guard native.

**Supplementary Tables = validation/cross-check layer ONLY (never pipeline input).** The authors' published per-gene results are retained at `data/human_validation/brain/hwang_ptsd/supp/` purely to sanity-check our own reclustered analysis:
- **Table 6** — DESeq2 pseudobulk, PTSD+MDD × celltype, **sig-only (FDR ≤ 0.05)**, `Genename/Condition/Celltype/log2FC/FDR`. Direction/effect-size cross-check against our DESeq2, including the **MIC PTSD-down / MDD-up SPP1 check** (the headline internal control). Matched estimator (DESeq2) but truncated.
- **Table 3 (MAST)** and **Table 4 (Wilcox)** — **full transcriptome** (FDR spans 0→~1), `Genename`+`Geneid`(ENSG) + `PTSD/MDD {MAST,Wilcox} log2FC`+`FDR`, all 16 author celltypes. Full-ranking concordance references (different estimator than our DESeq2 — for rank cross-check only).
- **Table 2** — cells per celltype/subtype: annotation sanity-check against our cluster proportions.
- **Table 1** — sample metadata: the diagnosis source + demographic covariates (Sex/Ancestry/AgeDeath/PMI/RIN/Smoking/Antidepressant).
- Author celltype scheme (from `utils/lists.py`): coarse `class` {EXC, INH, OLG, OPC, END, AST, MG} + fine `subclass` {CUX2/RORB/FEZF2/OPRK1 (EXC), LAMP5/KCNG1/VIP/SST/PVALB (INH)} → maps onto broad-7. `h10f_validate_vs_tables.py` reports per-celltype sign-concordance + the SPP1 panel.

**Scripts (planned):** `h10f_prep_hwang.py` (sparse-load counts → CSR; join Condition/Sex from Table 1 by channel prefix; ENSG→symbol; QC-light → HVG → scVI → Leiden → marker-annotate broad-7 → pseudobulk per donor×celltype → emit `h10f_hwang_rankings.parquet` mirroring `h10b_<ds>_rankings.parquet`), then `h10b_brain_rrho.py --dataset hwang --dx-contrast {PTSD,MDD} --tf`, `h10b_diagnostics.py`, `h10f_validate_vs_tables.py`, `h10_summary_plots.py`. Smoke-test on 2 channels before the full scVI run (GPU, ~2–4 hr).

### Phase 10: Reproducibility
Run logs, config snapshots, fixed seeds, `uv.lock` + `renv.lock` committed, git, `manifest.json` with checksums.

---

## 6. Pipeline Architecture

### Directory Structure (flat `scripts/`, no `src/` package)
```
Analysis/
├── config/        brain.yaml, placenta.yaml, dev_split.yaml,
│                  human_placenta_markers.yaml, cross_species_celltype_map.yaml
├── scripts/       _utils.py + numbered phase scripts + h09* human scripts + run_*.R / h_*.R workers
├── notebooks/     thin per-phase viewers
├── data/          raw Cell Ranger h5 (gitignored); data/human_validation/{placenta,brain}/...
├── results/       all outputs (gitignored)
├── refs/          reference atlases, msigdb_mouse.tsv, msigdb_human.tsv,
│                  mouse_human_orthologs.tsv, celltypist pkls, vento_tormo_2018/
├── sample_metadata.csv
├── run_pipeline_WS.sh
├── pyproject.toml + uv.lock + .python-version (3.12)
├── setup-remote.sh + scripts/install-r-packages.R + renv.lock
└── Makefile
```

### YAML Configuration
`config/brain.yaml` / `config/placenta.yaml` from `sample_metadata.csv`. Declarative Phase-8 blocks (`strata`, `composition`, `contrasts`, `stress_focused_cell_types`). Phase-9 cross-species compartment map in `config/cross_species_celltype_map.yaml` (`placenta_compartments:` block; brain block untouched, for ARM A).

### CLI Usage (Phase 9 placenta)
```bash
# Gunter-Rahman arm (WS):
uv run python scripts/h09a_prep_human_placenta.py ...      # → h09b → h09c → h09d
uv run python scripts/h09e_cross_species_rrho.py
uv run python scripts/h09f_overlap_null.py --n-perm 10000 --n-jobs 24
uv run python scripts/h09g_pathways_tf.py
uv run python scripts/h09h_leading_edge.py
uv run python scripts/h09_summary_plots.py                 # plots (Mac-runnable)

# Admati PE 2×2 arm (WS):
uv run python scripts/h09j_prep_admati_sc.py
uv run python scripts/h09k_admati_2x2.py --n-perm 5000 --n-jobs 24
uv run python scripts/h09k_diagnostics.py
uv run python scripts/h09k_plots.py                        # plots (Mac-runnable)
uv run python scripts/h09k_rrho_maps.py                    # plots (Mac-runnable)
```

---

## 7. Open Questions Before Next Stage
1. Behavioral data on offspring — incorporate alongside transcriptomics?
2. Corticosterone / HPA hormone measurements available?
3. Archived tissue for RNAscope/IHC validation? (priority targets: ECM thread Fn1/F13a1/Col8a1; Glul; perinatal IFN/ISGs; **and now the conserved placental hypoxia genes NDRG1/BNIP3L/ANGPTL4/DDIT4** validated across mouse + obesity + eoPE.)

---

## 8. Publication Strategy

**Realistic IF 12–15 targets:** Nature Communications, Molecular Psychiatry, Biological Psychiatry, Genome Biology.

**Strongest framing:** cross-tissue developmental cascade (E12.5/E18.5 placenta → P1/4W/3mo brain), two parallel stress-window arms, **plus cross-species validation of the placental signature in two independent human cohorts.**

**Four-figure plan:**
- **Fig 1** atlas + composition.
- **Fig 2 (disruption):** mirror disruption bar + effect-size collapse + k-preserving shuffle null. Claim: prenatal stress disrupts existing developmental programs rather than inducing new ones.
- **Fig 3 (cross-tissue cascade + persistence, 8f+8g):** three threads — ECM/mesenchymal, IFN/immune (perinatal-transient), gliogenesis.
- **Fig 4 (cross-species placenta validation):** the mouse prenatal-stress placental signature recapitulated in human placenta. **Two independent stressors:** (a) Gunter-Rahman obesity — decidua/vascular/trophoblast concordant (permutation-backed), HALLMARK_HYPOXIA conserved with ~40% gene-level overlap; (b) Admati PE — **two conserved axes:** eoPE recovers the same hypoxia program (third independent stressor; same genes), loPE recovers OXPHOS/electron-transport suppression. Honest nuance: GA-matched diagonal not stronger; immune weak; TF null. The convergence of mouse stress + human obesity + human eoPE on one hypoxia gene set (NDRG1/BNIP3L/ANGPTL4/DDIT4/SLC2A1/PGK1) is the headline.

**Decisive tier driver — NOW RESOLVED:** the **brain** cross-species arm (Phase 9 brain, COMPUTE-COMPLETE 2026-06-29). The mouse brain stress signature DOES recover in human cortex: strong directional neuronal-down concordance with human MDD (both sexes, mouse 4W), the 8f/8g IFN/immune thread recovered as microglial co-suppression in MDD/ASD (inverting in MS), plus ECM and gliogenesis threads. Both placenta and brain cross-species arms are now done → the cross-species package supports the upper-tier framing. **A fifth brain arm (Hwang/Girgenti PTSD/MDD, `h10f`, in progress) adds the closest adult *trauma* analog and a PTSD-vs-MDD internal directional control (predicted: PTSD microglia suppressed → concordant with mouse; MDD microglia up → divergent), strengthening the Fig 4 brain panel.**

**To strengthen:** RNAscope/IHC of top findings; behavioral validation; brain cross-species; ECHO-PATHWAYS measured-stress upgrade at revision.

---

## 9. Remote Workflow & Repo Layout
Remote via VPN+SSH from Mac; uv+renv; tmux; HTML reports. Code edits local → rsync to WS; WS results mirrored to Mac under `results_WS/`. **Phase-9 plotting can run entirely on the Mac** (small CSV/parquet inputs).

---

## 10. Next Steps

**Mouse arm COMPLETE through 8g. Phase 9 PLACENTA cross-species COMPLETE (Gunter-Rahman + Admati).**

**Immediate:**
1. **Phase 9 BRAIN ARM 5 — Hwang/Girgenti PTSD/MDD (`h10f`) — IN PROGRESS (2026-07-01).** Recluster the deposited counts (Zenodo 15186498) through the standard pipeline; diagnosis from Supp Table 1, own scVI/Leiden/marker annotation → own PyDESeq2 → `h10b` RRHO (PTSD-vs-CON + MDD-vs-CON). Validate against Supp Tables 1/2/3/6 (`h10f_validate_vs_tables.py`), especially the MIC PTSD-down/MDD-up SPP1 internal control. Folds into Fig 4 as the trauma arm + PTSD/MDD directional control.
2. **Phase 9 BRAIN ARM (4 datasets) — COMPUTE-COMPLETE (2026-06-29).** Four datasets (Velmeshev/Maitra/Nagy/Macnair) through `h10*`. **Remaining:** (a) Fig 4 brain figure refinement — the comprehensive plots exist (peak-keyed + pathway-keyed); iterate any panel that doesn't tell the story (the placenta Fig 4 was flagged "not yet representative" — apply that scrutiny); (b) write a standalone brain cross-species findings doc (like the 8f/8g one) before folding into the manuscript; (c) run the Velmeshev sensitivity variant (`--variant sensitivity`, Neu-NRGN/Neu-mat→ExN) into a quarantined subfolder + README; (d) optional: Herring 2022 age-anchor; (e) backport the `robust_class` margin guard to the placenta arm at revision.
2. **Figure refinement for Fig 4** — current placenta plots are functional but not yet "representative"; redesign deferred (options discussed: multi-stressor conserved-gene heatmap; compartment×pathway bubble grid; story-driven schematic).
3. **Documentation** — this summary + INSTRUCTIONS.md kept current (done 2026-06-25).

**Optional / revision-stage:**
- Fix + run `h09i` sn trophoblast targeted confirmation (dup-gene bug).
- ECHO-PATHWAYS dbGaP application (measured psychological stress).
- EMT/Myogenesis leading-edge confirmation; subcluster-level immune persistence; behavioral/corticosterone/archived-tissue validation.

---

## 11. Summary of Key Decisions

| Decision Point | Choice | Rationale |
|---|---|---|
| Ecosystem | Python (Scanpy + scvi-tools); R subprocess | Batch integration, scale |
| DE method | Pseudobulk + PyDESeq2, animal = unit | Cell-level DE incorrect |
| DE design | `~ sex + pool + group` per age | Sex+pool covariates |
| Composition | propeller (R), not scCODA | limma moderation for small n |
| Ambient RNA | SoupX per sample (R subprocess) | CellBender pickle bug |
| Trajectory | PAGA + DPT; no velocity | Flex probe-based |
| CCC | LIANA+; no cross-tissue CCC | BBB |
| Cross-tissue (8f) | E12.5→Early / E18.5→Late; placenta-whole × brain-{whole+regions} | sampling-window alignment |
| Cross-age (8g) | brain only; persistence per arm; region-resolved | placenta incomplete factorial |
| **Cross-species (Phase 9)** | **compartment-level pseudobulk RRHO; mouse↔human 1:1 ortholog bridge; functions lifted verbatim from 8f/8c** | **trophoblast subtypes lack 1:1 homology; method must match mouse arm** |
| **Phase-9 human GSEA** | **two single-species fgsea → intersect FDR<0.05 same-sign** | **cleanest provenance; mirrors 8c** |
| **Phase-9 mouse anchor** | **E18.5 Late-vs-Relaxed (both arms); + E12.5 Early-vs-Relaxed for Admati 2×2** | **all human stressors commensurable to one mouse signature** |
| **Admati PE arm** | **sc all-compartments (author-annotated) primary; sn trophoblast deferred (3v2)** | **powered, all compartments; modality caveat (sc) stated** |
| **Hwang PTSD arm (h10f)** | **recluster from deposited counts (own scVI/annotation/DE); authors' Supp Tables = cross-check only, NOT input** | **consistency with the other 4 brain arms + mouse outranks the tables-shortcut; obs (celltype) not deposited anyway** |
| **GA-matched 2×2 framing** | **tested, REJECTED — structure is by PE subtype not stage** | **diagonal not stronger; eoPE→hypoxia, loPE→OXPHOS** |
| Environment | uv + renv (not conda) | Conda blocked at firewall |
| Python↔R bridge | R as subprocess | Process isolation |

---

## 12. Environment & Deployment

**uv + renv (not conda).** Bootstrap `./setup-remote.sh`. Lock files committed (`uv.lock`, `renv.lock`, `.python-version` 3.12).

**Phase-8 deps:** `statsmodels>=0.14` (BH correction in 8b shuffle + 8f/8g).

**Phase-9 deps (added this arm):** `pybiomart` (ortholog map), `decoupler` 2.1.6 + `omnipath` (TF activity, CollecTRI human), `adjustText` (label de-overlap in plots), `pyarrow` (parquet for h09j/h09k rankings). R side: `SingleR` 2.14.0 (added + snapshotted) alongside SoupX/scran/DropletUtils/Matrix/rhdf5/msigdbr/optparse/BiocParallel. Mirror `renv.lock`/`pyproject.toml`/`uv.lock` back to Mac after env changes.

**CellTypist sklearn-1.7 patch; cuML via NVIDIA PyPI; renv Suggests workaround** — see INSTRUCTIONS.md.

---

## 13. Ambient RNA Correction (SoupX, locked 2026-06-10)
SoupX (R subprocess) replaces CellBender. cellranger filtered+raw → SoupChannel → scran::quickCluster → autoEstCont → adjustCounts. **Phase-9 Gunter-Rahman reuses this** (`h_run_soupx_from_raw.R`, with a knee/inflection cell-call up front since only raw matrices ship). **Admati sc has NO raw matrices → no SoupX (flagged); Admati uses author-annotated counts directly.**
