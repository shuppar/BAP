# Phase 8f Cross-Tissue (Placenta → Brain): Production Run, Findings & Queries

**Date:** 2026-06-25
**Scope:** Production wiring of `08f_cross_tissue.py`, interrogation of all six views, and the resulting Fig 3 interpretation. Brain + placenta, reproducible from completed 8b/8c CSVs (no DE re-runs).

---

## 1. What 8f does

Six analytical views, all operating on completed 8b (`08b_de_results.csv`) and 8c (`08c_pathway_results.csv`, `08c_tf_activity.csv`, `08c_pathway_leading_edge.csv`) tables. The cross-tissue join is always **placenta-WHOLE × brain-{whole + 13 regions}** (placenta is whole-only), across two biologically aligned arms:

- **Early arm:** E12.5 placenta (Early) → P1 / 4W / 3mo brain (Early)
- **Late arm:** E18.5 placenta (Late) → P1 / 4W / 3mo brain (Late; P1 Late flagged `confounded_with_pool`)

| View | Output | Status |
|---|---|---|
| 1. DEG overlap (hypergeometric) | `08f_deg_overlap.csv` | sparse (2 FDR<0.05 at \|log2FC\|>1) |
| 2. RRHO (rank-rank) | `08f_rrho_summary.csv` | rich (4,430 rows) |
| 3. Pathway concordance | `08f_pathway_concordance.csv` | **richest / headline (244k rows)** |
| 4. LR cross-tissue (ligand×receptor) | `08f_lr_cross_tissue.csv` | moderate (see §3) |
| 5. TF concordance | `08f_tf_concordance.csv` | partial corroboration |
| 6. Overlap-gene ORA | `08f_overlap_enrichment.csv` | supplementary |

---

## 2. Methodology & fixes applied this session

Four silent-correctness bugs were found and fixed before the production run (the same class that bit 8c summary plots):

1. **`_extract_age` parsed a non-existent format.** It expected `age-4W`; the 8b/8c `group_level` column holds the age **directly** (`P1`, `4W`, `3mo`, `E12.5`, `E18.5`). The old parser returned `None` for every row → empty join that still exits cleanly. Fixed to return `group_level` as-is.
2. **No `sex` filter.** All three strata (combined/M/F) were piling into one celltype → duplicate gene rows corrupting RRHO `set_index("gene")` and inflating overlap counts. Added `--sex` (default `combined`).
3. **No `level` filter.** whole + 13 regions mixed. Added asymmetric level handling: placenta pinned to `whole`; brain iterates `whole` + regions, with `brain_level` threaded through every output row, FDR grouping, and plot path.
4. **Mixed-dtype reads.** All `pd.read_csv` now use `low_memory=False`.

Other engineering:
- **RRHO vectorized.** The per-cell-type-pair RRHO matrix was an O(bins²) double loop rebuilding a Python `set` and calling `scipy.hypergeom.sf` per cell (~1,600×/matrix). Replaced with a 2D membership-matrix product + one array-valued `hypergeom.sf`. **Bit-identical output (max abs diff 0.0), ~36× faster** on a 6k-gene case (more on larger sets). This was necessary because the 14× region multiplier made the serial version take hours.
- **DEG cutoff.** Initial production run at `|log2FC|>1` (8b-locked). Re-run at `|log2FC|>0.5` to populate the overlap/LR views (RRHO uses ranked stats and pathway/TF use FDR, so those are cutoff-independent). The 0.5 run is the basis for the LR numbers below.
- **Plot-only quantile density floor** (`--plot-quantile`, default 0.75): floors LR bar + pathway/TF dotplots within each slice; **never** floors distribution plots (scatters) or heatmaps. CSVs always complete.
- **Parallelization deliberately NOT added.** The script is single-threaded across slices, but the run completes on its own and the per-slice parallelization would require pre-filtering to dodge pickling the 20M-row brain DE frame — judged too risky for a just-validated join. Revisit only if reruns become frequent.

`contrast_family` lift into `_utils.py` was **not** needed for 8f: the `ARMS` dict carries explicit per-arm contrast strings (`early_vs_relaxed_E12.5`, `late_vs_relaxed_E18.5`, `*_per_age`) that match the CSVs exactly.

---

## 3. Findings

### 3.1 LR cross-tissue table (View 4) — moderate, NOT a stress-hormone relay

At `|log2FC|>0.5`: **93 candidate LR pairs** (59 discordant, 26 concordant_down, 8 concordant_up = 34 concordant), concentrated in the **Early→P1** window (23 of 34 concordant).

- **Only 1 pair touches the canonical stress axis** (GR/CRH/cytokine/serotonin/BDNF curated set), and it is **discordant**: `Ntf3 → Ntrk3` (placental Mesenchyme → brain Excitatory neurons, Late/P1). Cannot be framed as a clean relay.
- **What carries the signal is ECM / angiogenic / growth-factor signalling.** Top concordant ligands: **Fn1 (8 pairs)**, Fgf2/Fgf11, Efnb2 (ephrin-B2), Vegfa/Pgf, Col9a3, Adam15, F13a1, Hp.
- The dominant relationship is actually **anti-correlated** (59/93 discordant), which complicates a simple cascade narrative.

**Read:** the LR view is a reasonable mechanistic *supplementary* table — a candidate developmental-signalling/ECM cascade anchored on Fn1 — but it is thin and stress-axis-weak. It should not carry Fig 3 alone.

### 3.2 Pathway concordance (View 3) — STRONG, persistent, the real Fig 3

After dropping M8 cell-identity sets (per the 8c convention), multiple **Hallmark/Reactome/GO immune pathways are concordant placenta→brain across ALL THREE brain ages (P1, 4W, 3mo), in BOTH arms.**

Dominant pattern = **coordinated co-SUPPRESSION** (concordant_down, all 3 ages):

| Pathway | Arm | Collection | n rows |
|---|---|---|---|
| REACTOME_CYTOKINE_SIGNALING_IN_IMMUNE_SYSTEM | Early | M2 | 156 |
| HALLMARK_INTERFERON_GAMMA_RESPONSE | Late | MH | 120 |
| HALLMARK_INTERFERON_ALPHA_RESPONSE | Late | MH | 117 |
| HALLMARK_COMPLEMENT | Early | MH | 100 |
| HALLMARK_ALLOGRAFT_REJECTION | Late | MH | 63 |
| HALLMARK_TNFA_SIGNALING_VIA_NFKB | Early | MH | 60 |
| HALLMARK_IL2_STAT5_SIGNALING | Late | MH | 60 |
| HALLMARK_INFLAMMATORY_RESPONSE | Late | MH | 55 |
| REACTOME_INTERFERON_SIGNALING | Early | M2 | 40 |
| HALLMARK_COMPLEMENT | Late | MH | 39 |

Exception — **`TNFA_SIGNALING_VIA_NFKB` is bidirectional** in the Early arm (concordant_up 176 rows / 5 ct_pairs; concordant_down 60 / 5). Direction is cell-type-pair-specific → report as context-dependent, not a clean up or down.

### 3.3 TF concordance (View 5) — partial corroboration

Same persistence cut on `08f_tf_concordance.csv` (concordance classes `concordant_activated` / `concordant_repressed`):

- **Persistent (all 3 ages):** **Jun & Jund (AP-1) concordant_repressed** in the Late arm (48 / 12 rows) — AP-1 is a canonical stress-responsive IEG family. Clean and interpretable.
- **Nr3c1 (glucocorticoid receptor) concordant_activated** (Early 3mo/4W; Late P1) — the stress axis the LR view missed surfaces here at the regulator level. Modest but on-target.
- **NFkB family (Nfkb1/Rela/Rel) bidirectional, low-n** — matches the bidirectional TNFα pathway result.
- **Interferon TFs (Irf1/Irf9/Stat1) appear concordant_repressed but only at 1–2 ages, low-n.** The IFN regulatory layer is NOT persistent here — so the IFN-suppression claim rests on pathway + leading-edge evidence, not TF activity.

### 3.4 Leading-edge confirmation (8c `08c_pathway_leading_edge.csv`)

Pulled top recurrent leading-edge genes per persistent pathway, **whole/combined, per tissue**. This disambiguates real regulated program vs composition:

- **IFN-γ — real, shared.** Brain ∩ placenta share `Ifih1, Ddx60, Herc6, Parp14, Tap1, B2m, Xaf1, Nmi, Rnf213`; placenta adds textbook ISGs `Stat1, Isg15, Usp18, Irf9, Ifit3, Gbp3`.
- **IFN-α — real, shared, cleanest.** Both tissues driven by `Bst2, Ifitm3, Ifih1, Adar, Irf7, Irf9, Ifi27, B2m, Rtp4, Herc6`; placenta adds `Oas1a, Isg15, Usp18, Cmpk2`. `Irf7`/`Stat1` in the leading edge IS the regulatory fingerprint, so the IFN story is self-supporting at the effector level despite the noisy TF-activity layer.
- **Complement — messier; rename.** Leading edge is NOT the classic cascade (`C1q*`, `C3`, `C4`). Brain: `Clu, Ctsl, Fyn, Tnfaip3, Pla2g4a`; placenta: `Mmp15, Plaur, Fn1, Pdgfb, Mmp14, Timp2`. This is the ECM/protease arm of the Hallmark set (and `Fn1` ties back to the LR result). Report as matrix-remodeling/protease, not classical complement.
- **Cytokine signalling — proteasome-driven.** `Psm*` subunits, `Elob`, `Skp1`, ubiquitin machinery in both tissues. The protein-degradation arm of the Reactome set, not cytokine ligands/receptors.
- **Microglial-ageing (M8) — real state shift, NOT composition.** Leading edge is metabolic/redox/proteostasis genes (`Cst3, Prdx1, Cd81, Gsn, Tmed3, Ndufb9` brain; `Gpx1, Prdx5, Sparc, Lum, Cd63, Eef1a1` placenta), NOT pan-myeloid identity markers (`Csf1r/Aif1/P2ry12` absent). So it's a genuine shared low-grade metabolic/redox state shift — but it reads as a broad metabolic/stress-proteostasis signature, not specifically "microglial ageing." Keep as a hypothesis with the M8 label treated as a loose descriptor.

---

## 4. Defensible Fig 3 framing

> **Prenatal stress imposes a persistent, regimen-shared suppression of interferon signalling concordant between placenta and brain from birth to 3 months, driven by shared interferon-stimulated gene effectors (Irf7, Stat1, Isg15, Bst2, Ifih1, B2m…), with cell-type-context-dependent TNFα–NFkB modulation and a secondary ECM/protease thread.**

- **Primary, robust:** persistent cross-tissue IFN-α/γ co-suppression (pathway + leading-edge; both arms, all 3 ages, both tissues). Independent of the thin LR table.
- **Secondary:** TNFα–NFkB context-dependence; AP-1 (Jun/Jund) repression; Nr3c1 activation at the TF level.
- **Mechanistic supplementary:** the LR/ECM cascade (Fn1-anchored).
- **Caveats to carry:** rename "complement"→ECM/protease and "cytokine signalling"→proteasome per their leading edges; TNFα direction is cell-type-pair-specific; 59/93 LR pairs are discordant; `within_group`/region contrasts inherit the pool-age confound.

---

## 5. Queries run (reproducible)

All run on WS from `/home/poller/BAP-BrainPlacenta`. **Do not `cut`/`awk` these CSVs** — the `pair` column contains an internal comma (`"['Early_Stress', 'Late_Stress']"`) that shifts positional fields. Use pandas.

**Production run (looser cutoff):**
```bash
rm -rf results/brain/plots/08f_cross_tissue results/brain/tables/08f_cross_tissue
uv run python -u scripts/08f_cross_tissue.py \
  --brain-config config/brain.yaml --placenta-config config/placenta.yaml \
  --logfc-cutoff 0.5 2>&1 | tee logs/08f_cross_tissue_logfc05.log
```

**LR table inspection:**
```python
import pandas as pd
lr = pd.read_csv('08f_lr_cross_tissue.csv')
lr['direction'].value_counts()
lr[lr['stress_axis'].notna() & (lr['stress_axis']!='')]          # real stress-axis hits (note .notna())
lr[lr['direction'].str.startswith('concordant')]['ligand'].value_counts().head(10)
```

**Pathway concordance — persistent inflammatory cut (drop M8):**
```python
pw = pd.read_csv('08f_pathway_concordance.csv', low_memory=False)
pw = pw[pw['collection'].isin(['MH','M2','M5'])]
conc = pw[pw['concordance_class'].isin(['concordant_up','concordant_down'])]
kw = ['INFLAMMAT','TNF','IL6','IL2_STAT','JAK_STAT','INTERFERON','IFN','COMPLEMENT',
      'NFKB','CYTOKINE','INNATE','TGF_BETA','IL1','CHEMOKINE','TOLL','ALLOGRAFT']
inf = conc[conc['pathway'].str.upper().str.contains('|'.join(kw), na=False)]
rec = (inf.groupby(['arm','pathway','concordance_class'])
          .agg(n_rows=('brain_age','size'),
               ages=('brain_age', lambda s: ','.join(sorted(s.unique()))),
               n_ages=('brain_age','nunique'),
               n_ctpairs=('brain_celltype','nunique'),
               collection=('collection','first')).reset_index())
rec[rec['n_ages']==3].sort_values(['concordance_class','n_rows'], ascending=[True,False])
```

**TF concordance — same persistence cut:**
```python
tf = pd.read_csv('08f_tf_concordance.csv', low_memory=False)
conc = tf[tf['concordance_class'].str.startswith('concordant', na=False)]
fam = ['IRF','STAT','NFKB','REL','JUN','FOS','NR3C','CEBP','SPI1','RUNX','IKZF']
sel = conc[conc['TF'].str.upper().str.contains('|'.join(fam), na=False)]
(sel.groupby(['arm','TF','concordance_class'])
    .agg(n_rows=('brain_age','size'),
         ages=('brain_age', lambda s: ','.join(sorted(s.unique()))),
         n_ages=('brain_age','nunique'),
         n_ctpairs=('brain_celltype','nunique')).reset_index()
    .sort_values(['n_ages','n_rows'], ascending=False))
```

**Leading-edge check (6.7 GB file — chunked, pandas only, never awk):**
```python
PATHS = {'brain':'results/brain/tables/08c_pathways/08c_pathway_leading_edge.csv',
         'placenta':'results/placenta/tables/08c_pathways/08c_pathway_leading_edge.csv'}
TARGETS = ['INTERFERON_GAMMA','INTERFERON_ALPHA','HALLMARK_COMPLEMENT',
           'CYTOKINE_SIGNALING_IN_IMMUNE_SYSTEM','MICROGLIAL_CELL_AGEING']
cols = ['sex','level','collection','pathway','gene','log2FC','direction']
for tis, f in PATHS.items():
    keep = []
    for ch in pd.read_csv(f, usecols=cols, chunksize=2_000_000, low_memory=False):
        ch = ch[(ch['level']=='whole') & (ch['sex']=='combined')]
        m = ch['pathway'].str.contains('|'.join(TARGETS), na=False)
        if m.any(): keep.append(ch[m])
    df = pd.concat(keep, ignore_index=True)
    for t in TARGETS:
        sub = df[df['pathway'].str.contains(t, na=False)]
        if not sub.empty:
            print(tis, t, sub['gene'].value_counts().head(15).index.tolist())
```

**Leading-edge file schema:** `tissue, sex, contrast, flag, group_level, pair, level, celltype, n_donors_total, reliability, note, collection, pathway, NES, pathway_FDR, leading_edge_rank, gene, log2FC, rank_stat, direction` (one row per leading-edge gene).

---

## 6. Open items / next

- **8g audit (next).** Check `08g_cross_age.py` for the same four issues as 8f (group_level-as-age parse, sex/level filtering, exact-string vs `contrast_family` contrast matching, `low_memory=False`) before any production run. Pointed question for 8g: do the **same ISGs** (Irf7, Stat1, Isg15, Bst2, Ifih1…) **persist across ages WITHIN the brain** — pairing within-brain persistence (8g) with the cross-tissue concordance (8f) on the same interferon program.
- **DEG-overlap View 1 is sparse even at 0.5** (2 FDR<0.05) — pathway concordance is the better vehicle; do not lean on View 1.
- **Possible follow-ups:** confirm the TNFα-up vs -down split is driven by distinct, coherent cell types (not noise); subcluster-level pathway concordance for the immune compartment.
