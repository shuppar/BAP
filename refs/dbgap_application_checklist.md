# dbGaP / Synapse Application Checklist — Stage 2 controlled-access datasets

This is a checklist of what YOU (the PI or PI-delegated requester) need to do
to apply for the controlled-access human validation datasets recommended in
the cross-species validation plan. Claude cannot do these for you — they
require your eRA Commons account and your institutional signing official.

Typical approval time: **2–6 months**. Apply as early as possible in the
analysis cycle so the data is in hand by the time the mouse pipeline produces
DE tables to validate against.

---

## What you need before starting any application

- **eRA Commons account** — apply at https://commons.era.nih.gov/ if you don't
  have one. Requires institution sign-off (1–2 weeks).
- **Institutional signing official (SO)** — find them via your university's
  research administration office. The SO must approve every dbGaP request.
- **NIH-style data use plan** — short paragraph stating the research question,
  what you'll do with the data, who has access, and security/retention.
- **Local data security review** — many institutions require an IRB review for
  using human genomic data. Start the IRB process in parallel with dbGaP.

---

## Dataset 1 — ECHO-PATHWAYS placental bulk RNA-seq (HIGHEST PRIORITY)

**Why:** the only large dataset directly linking maternal prenatal stress to
placental transcriptomics. RRHO2 against your mouse E18 trophoblast/decidua
pseudobulk is the headline cross-species comparison.

- **Accessions:** `phs003619.v1.p1` (CANDLE, n=794) + `phs003620` (GAPPS, n=289)
- **Total N:** 1,065 mother-child pairs (1,029 with childhood-trauma data,
  874 with prenatal-SLE data)
- **Modality:** placental bulk RNA-seq + measured maternal prenatal stressful
  life events (SLEs) + childhood traumatic events (CTEs)
- **Key paper:** Cao-Lei et al. 2024 *Mol Psychiatry*
- **Where to apply:** https://dbgap.ncbi.nlm.nih.gov/
- **Citations needed in application:** Cao-Lei 2024 + your study aims

## Dataset 2 — Pique-Regi/Garcia-Flores term placenta scRNA-seq

**Why:** the cleanest single-cell-level human term placenta reference (matches
mouse E18 most closely). The Marsh 2022 mid-gestation open dataset is a
partial substitute but not term.

- **Accession:** `phs001886.v5.p1`
- **N:** 42 term placentas (24 term-in-labor + 18 term-no-labor)
- **Modality:** scRNA-seq of chorioamniotic membrane, basal plate, villi
- **Key papers:** Pique-Regi 2019 *eLife* (v1) + Garcia-Flores 2024 *Sci Transl Med* (v5)
- **Where to apply:** https://dbgap.ncbi.nlm.nih.gov/

## Dataset 3 — PsychENCODE brainSCOPE adult cortex atlas

**Why:** 388 donors, 2.8M PFC nuclei, harmonized across SZ/ASD/BD/AD/control.
Broadest cross-disorder adult reference for your 12W mouse timepoint.

- **Accession:** Synapse, requires PsychENCODE Data Use Certificate
- **N:** 388 donors, 2.8M nuclei
- **Modality:** snRNA-seq + snATAC-seq + snMultiome (uniformly processed)
- **Key paper:** Emani et al. 2024 *Science* (brainSCOPE)
- **Where to apply:** https://www.synapse.org/ + SAGE Bionetworks
- **Citation:** Emani 2024 + brainSCOPE consortium description

## Dataset 4 — Hwang/Girgenti PTSD/MDD dlPFC atlas

**Why:** most recent and uniquely valuable — PTSD/MDD are the most plausible
adult outcomes of prenatal stress. Cell-type-resolved DEGs in micro, astro,
OL lineage, both Ex and In neurons.

- **Accession:** Zenodo 15186498 (processed, open) + Synapse (raw, controlled)
- **N:** 111 donors, >2M dlPFC nuclei
- **Modality:** snRNA-seq + snATAC-seq
- **Key paper:** Hwang/Girgenti et al. 2025 *Nature*
- **Where to apply:** Processed data on Zenodo is open (no application).
  Raw data via PsychENCODE/Synapse — same channel as Dataset 3.
- **Note:** start with the Zenodo processed data immediately (open access);
  apply for raw only if you need it for re-analysis.

---

## Practical sequencing — recommended application order

1. **Week 0–1:** confirm eRA Commons, identify SO, draft data use plan.
2. **Week 1–2:** submit `phs003619` + `phs003620` (ECHO-PATHWAYS). Highest priority.
3. **Week 1–2:** submit `phs001886.v5.p1` (term placenta). Same dbGaP channel,
   same paperwork; bundle both in the same SO meeting to save round-trips.
4. **Week 2–4:** Synapse data-use certificate for PsychENCODE brainSCOPE.
5. **Week 4+:** monitor approvals. Typical timeline 2–6 months for dbGaP;
   Synapse is usually faster (weeks).

While these are pending, run RRHO2 on the four open Stage-1 datasets
(Nagy/Maitra, Velmeshev, Herring, Marsh). That gives a defensible validation
even before controlled-access data lands.

---

## Caveats

- **No published human dataset combines placenta snRNA-seq with measured
  maternal psychological stress at large scale.** ECHO-PATHWAYS is the
  closest — but it's BULK RNA-seq, not single-cell. The RRHO comparison from
  your mouse data will be: mouse pseudobulk per cell type ↔ ECHO bulk
  (across thousands of placentas). That's defensible and standard practice.
- **Sex stratification:** Maitra 2023 found striking sex-specific cell-type
  contributions to MDD DEGs (female: Mic1 microglia 38%; male: ExN10_L46
  53%). If your mouse study is mixed-sex, plan to stratify the human
  comparison by sex when applying.
- **Several accessions reported in literature have minor inconsistencies.**
  The validation doc flags Suryawanshi 2018 GEO ID, Wang Q 2022 CNGB ID, and
  the exact raw-data accession for Hwang/Girgenti 2025 as unverified.
  Confirm against each paper's Data Availability statement before submitting.
