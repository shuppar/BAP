# Handoff v2 — workstation setup complete, ready for Phase 0

This handoff supersedes `HANDOFF_to_new_chat.md`. Use it when opening a fresh
chat to continue from where today left off.

## One-paragraph summary

> Workstation is fully bootstrapped. All Python + R packages installed,
> CellBender sidecar venv built, ABC atlas downloaded and processed
> (`refs/abc_brain_ref.h5ad` + `refs/celltypist_brain_adult.pkl`), MSigDB gene
> sets fetched, `brain.yaml` regenerated with all 3 CellTypist models (P1
> built-in + 4W/3mo from the bespoke ABC build pointing to the same .pkl).
> Human Stage-1 validation datasets currently downloading in a background
> tmux. Phase 9 cross-species scaffold + cell-type map YAML written but
> per-dataset loaders are stubs. The next step is **Phase 0 validation on
> brain**, then **Phase 1 CellBender brain** (the longest single step,
> ~1 day on GPU).

## Workstation infra (already done)

- **Project root:** `/home/poller/BAP-BrainPlacenta/` (NVMe, fast).
- **Raw data:** `/media/poller/PollerLab-1/BAP-data1/Analysis/data/` (USB-HDD),
  reached via the `BAP-BrainPlacenta/data` symlink. 57 unique samples
  (34 brain + 23 placenta; CES2.3 dropped as duplicate).
- **SSH:** `ssh poller@172.17.213.147`.
- **Local repo on Mac:** `/Users/shuppar/Downloads/BAP_data_1/Analysis/`.
- **Rsync convention (from Mac):**
  ```bash
  rsync -av --progress --chmod=Fu+x \
    --exclude='results/' --exclude='data/' --exclude='.venv/' \
    --exclude='.venv-cellbender/' --exclude='__pycache__/' --exclude='.git/' \
    --exclude='.DS_Store' --exclude='*.h5ad' --exclude='logs/' \
    /Users/shuppar/Downloads/BAP_data_1/Analysis/ \
    poller@172.17.213.147:/home/poller/BAP-BrainPlacenta/
  ```
  `--chmod=Fu+x` is mandatory — without it `.sh` files lose exec bit.

## What was changed in this session

All these files were updated on the Mac and rsynced to the workstation:

| File | What changed |
|---|---|
| `scripts/install-r-packages.R` | Bioc bumped 3.19 → 3.21 (R 4.5 compat); removed CellChat (unused, transitive deps fail); added msigdbr; wrapped `renv::snapshot` in tryCatch |
| `scripts/prepare_brain_reference.py` | NEW. Downloads ABC atlas WMB-10Xv3, subsamples to ~300 cells/subclass (~92K cells), trains CellTypist .pkl |
| `scripts/build_yaml.py` | `REFERENCE_CONFIG.brain` now points at `refs/abc_brain_ref.h5ad` + `refs/celltypist_brain_adult.pkl` (same .pkl for both 4W and 3mo — ABC atlas is adult P56) |
| `setup-remote.sh` | Steps 3 (abc_atlas_access install), 7 (idempotent MSigDB fetch), 8 (brain reference build); `--skip-references` flag |
| `run_pipeline_WS.sh` | Pre-flight now checks for both reference files |
| `scripts/download_human_validation.sh` | NEW. Downloads 4 open Stage-1 human datasets (idempotent) |
| `scripts/queue_human_downloads.sh` | NEW. Polls for abc_ref tmux to end, then triggers downloads |
| `refs/dbgap_application_checklist.md` | NEW. Stage-2 dbGaP/Synapse application guide |
| `scripts/09_cross_species_validation.py` | NEW. Phase 9 scaffold with RRHO2 utility, pseudobulk DE, mouse→human gene mapping. **Per-dataset loaders are stubs.** |
| `config/cross_species_celltype_map.yaml` | NEW. Mouse → human cell-type mapping per dataset. Best-effort labels; verify against actual data after download. |
| `INSTRUCTIONS.md` | Added workstation infrastructure section + rsync convention |

## What's done on the workstation right now

- ✅ `uv sync` (.venv) and abc_atlas_access installed
- ✅ R 4.5.2 + Bioc 3.21 + scDblFinder/edgeR/SingleCellExperiment/DropletUtils/speckle/limma/msigdbr
- ✅ CellBender sidecar `.venv-cellbender/` (Python 3.10)
- ✅ `refs/msigdb_mouse.tsv` (44 MB)
- ✅ `refs/abc_brain_ref.h5ad` (4.6 GB, 92,463 cells × 32,285 genes)
- ✅ `refs/celltypist_brain_adult.pkl` (40 MB, ~330 subclass labels)
- ✅ `refs/abc_atlas/` cache (~50 GB of regional h5ads; can delete after Phase 9 if disk pressure)
- ✅ `config/brain.yaml` has all 3 CellTypist models + ref_h5ad + labels_key=cell_type + region_key=region
- 🔄 `human_dl` tmux: downloading Maitra 2023 (~700 KB/s); will continue to Nagy, Velmeshev, Herring, Marsh

## Verified facts not to re-litigate

- **Project name on workstation is `BAP-BrainPlacenta`**, renamed from `Analysis`. Same files, different folder name.
- **CellChat was intentionally removed** from `install-r-packages.R`. Phase 8e uses LIANA+ in Python; CellChat is dead code.
- **One CellTypist .pkl for both 4W and 3mo.** ABC atlas is adult P56, biologically close to both. Saves compute, no biology lost.
- **WMB-10Xv3** was chosen as the ABC chemistry (closest to user's 10x Flex). WMB-10Xv2 deliberately excluded.
- **Naive uppercase mouse→human ortholog mapping** in Phase 9 (works for ~85% of orthologs). HGNC HCOP TSV is a TODO for higher accuracy.
- **placenta.yaml deliberately has `celltypist_models: {}` and `ref_h5ad: null`** — no placenta reference built yet. Phase 7 placenta will fall back to marker-only. Building a Marsh & Blelloch placenta reference is a separate decision (probably faster than ABC, ~hours).
- **Stage-2 controlled-access datasets** (ECHO-PATHWAYS, PsychENCODE, Hwang/Girgenti, Pique-Regi) need the user to apply via dbGaP/Synapse personally. Claude can't help with submission. Reference: `refs/dbgap_application_checklist.md`.

## Known small issues / TODOs

- `queue_human_downloads.sh` polls for tmux session, not the script — interactive `tmux new -s X <cmd>` sessions persist after command ends. Workaround: `tmux kill-session -t abc_ref` manually when builds finish, or use `tmux new -s X -d "<cmd>"` (detached mode).
- Phase 9 per-human-dataset loaders are STUBS that raise NotImplementedError with docstrings explaining what's needed. Fill in after each dataset finishes downloading.
- `cross_species_celltype_map.yaml` uses best-effort labels/regex patterns from the published papers — verify against actual downloaded `meta.tsv` files before running Phase 9.
- Some compile-time Ubuntu sys-libs were missing for `ragg`/`fs` (transitive deps of plotting packages). We worked around by skipping snapshot validation. If you ever want to install ragg properly, the missing lib is likely `libfreetype-dev` (the new name; setup-remote.sh installs `libfreetype6-dev` which is the old name).

## Today's failure modes (so the next Claude doesn't repeat)

- Bioc 3.19 → 3.21 (R 4.5 needs Bioc 3.21).
- `fs`/`ragg` install failure from leftover CellChat transitive deps → dropped CellChat entirely.
- Orphan packages in `renv/library` from earlier failed run → `rm -rf renv renv.lock` before retry.
- `.Rprofile` referenced deleted `renv/activate.R` → remove `.Rprofile` and let renv::init regenerate.
- `renv::snapshot` validation fails when broken transitive deps remain → wrap in tryCatch, verify pipeline packages explicitly.
- ABC atlas API: `directory='WMB-10Xv3'` (chemistry only), `file_name='WMB-10Xv3-{region}/raw'` (region in file_name, not directory).
- AnnData column collision after metadata join: drop conflicting cols from obs first.
- rsync without `--chmod=Fu+x` strips exec bit on `.sh` files.

## Suggested next-chat opening message

Paste this as the first message of the new chat:

> I'm continuing workstation execution of the snRNA-seq prenatal stress
> pipeline. Workstation: 258 GB RAM, 56 cores, RTX 4500 Ada (24 GB VRAM),
> R 4.5.2 + Bioc 3.21, uv-managed Python venv. Project root
> `/home/poller/BAP-BrainPlacenta/` (NVMe), raw data on USB-HDD via symlink.
> Per `HANDOFF_v2.md` (attached), all environment + R + CellBender +
> ABC brain reference + CellTypist .pkl + MSigDB are built and wired into
> `config/brain.yaml`. Human Stage-1 validation downloads are running in a
> background tmux. Phase 9 scaffold exists with per-dataset loader stubs.
>
> Today's specific goal: **[FILL IN — e.g. "Phase 0 brain validation +
> Phase 1 CellBender brain (kick off in tmux)", or "QC + doublets brain
> after Phase 1 done"]**.
>
> Working from `run_pipeline_WS.sh`. When something fails I'll paste the
> full traceback and the command that produced it.

## Files to upload to the new chat

**Project files (paste into the new project's Files):**
- `snRNAseq_project_summary.md` (unchanged)
- `INSTRUCTIONS.md` (the updated one from today — has workstation infra section)
- `HANDOFF_v2.md` (this file)
- `Human_Datasets_for_Cross-Species_Validation...md` (unchanged)

**First-message attachments (upload at chat start):**
- `run_pipeline_WS.sh` (updated)
- `_utils.py`
- All phase scripts: `01_validate.py`, `02_qc.py`, `03_doublets.py`,
  `04_integration_prep.py`, `05_integration.py`, `06_clustering.py`,
  `07_annotation.py`, `07b_subcluster.py`, `07c_label_transfer.py`,
  `07d_subcluster_annotate.py`, `08a_composition.py`, `08b_de.py`,
  `08c_pathways.py`, `08d_trajectory.py`, `08e_communication.py`,
  `_08e_plots_*.py`, `08f_cross_tissue.py`, `08g_cross_age.py`,
  `09_cross_species_validation.py`
- R scripts: `run_scdblfinder.R`, `run_propeller.R`, `fetch_genesets.R`
- Configs: `config/brain.yaml`, `config/placenta.yaml`,
  `config/subcluster_markers.yaml`, `config/cross_species_celltype_map.yaml`
- `sample_metadata.csv`

**Probably skip (save context):**
- `prepare_brain_reference.py` (done, won't be re-run)
- `prepare_reference.py` (validator, not used directly)
- `install-r-packages.R`, `setup-remote.sh` (done, won't re-run)
- `download_human_validation.sh`, `queue_human_downloads.sh` (running)
- `build_yaml.py` (done unless schema changes)

## Cadence reminders for next chat

- Use tmux for any multi-minute job. `Ctrl-b d` to detach.
- Wrap GPU phases (1 and 5) in tmux. Phase 1 CellBender is ~1 day per tissue.
- Don't try to run the full pipeline in one session; checkpoint between phases.
- Recommended order: brain through 8e in tmux while placenta also through 8e
  in a parallel tmux (different terminal). Then 8f cross-tissue. Then 8g
  cross-age (brain only). Then Phase 9 cross-species (after human downloads
  done AND mouse 8b done).
