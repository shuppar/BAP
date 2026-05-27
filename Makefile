# Makefile for the snRNA-seq prenatal-stress pipeline.
#
# Usage examples:
#   make setup           # one-time: bootstrap env on the workstation
#   make sync            # rsync local code -> remote workstation
#   make validate        # run Phase 0 validation (brain)
#   make qc              # run QC phase (brain)
#   make qc-placenta     # run QC phase (placenta)
#   make gpu-check       # verify GPU is visible to PyTorch
#   make lock            # regenerate Python + R lockfiles
#   make clean           # remove venvs (keeps lockfiles)
#
# Configuration: edit the variables below or override on command line:
#   make sync REMOTE=user@host PROJECT_DIR=/scratch/myproj

# =============================================================================
# Configuration
# =============================================================================

# Remote workstation SSH alias (configure in ~/.ssh/config).
# The doc uses 'remote-snRNA' but adjust to whatever you set up.
REMOTE       ?= remote-snRNA
PROJECT_DIR  ?= ~/snrna-project

# Default tissue config for the analysis-running targets
TISSUE       ?= brain
CONFIG       := config/$(TISSUE).yaml

# Python invocation: use `uv run` so it works whether or not the venv is activated
PY           := uv run python

# =============================================================================
# Targets — local (run on Mac)
# =============================================================================

.PHONY: help
help:
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | \
		awk 'BEGIN {FS = ":.*?## "}; {printf "  \033[36m%-20s\033[0m %s\n", $$1, $$2}'

.PHONY: sync
sync: ## rsync code from local Mac -> remote workstation (excludes data, venvs, caches)
	rsync -avzP \
		--exclude '.venv*' \
		--exclude '.renv-cache' \
		--exclude '__pycache__' \
		--exclude '.pytest_cache' \
		--exclude '.ruff_cache' \
		--exclude '.mypy_cache' \
		--exclude '.ipynb_checkpoints' \
		--exclude 'results' \
		--exclude 'data' \
		--exclude '*.h5ad' \
		--exclude '*.h5' \
		--exclude '.DS_Store' \
		--exclude '.git/objects' \
		./ $(REMOTE):$(PROJECT_DIR)/

.PHONY: sync-back
sync-back: ## rsync results/plots/tables from remote -> local (for figures, reports)
	rsync -avzP \
		$(REMOTE):$(PROJECT_DIR)/results/$(TISSUE)/plots/ \
		./results/$(TISSUE)/plots/
	rsync -avzP \
		$(REMOTE):$(PROJECT_DIR)/results/$(TISSUE)/tables/ \
		./results/$(TISSUE)/tables/
	rsync -avzP \
		$(REMOTE):$(PROJECT_DIR)/results/$(TISSUE)/report.html \
		./results/$(TISSUE)/report.html

.PHONY: ssh
ssh: ## SSH into the remote workstation
	ssh $(REMOTE)

# =============================================================================
# Targets — remote (run on workstation, via Make or manually)
# =============================================================================

.PHONY: setup
setup: ## ONE-TIME bootstrap: install uv, Python deps, R, R packages, CellBender venv
	./setup-remote.sh

.PHONY: setup-skip-r
setup-skip-r: ## bootstrap without R (use if R already configured)
	./setup-remote.sh --skip-r

.PHONY: gpu-check
gpu-check: ## verify PyTorch sees the GPU and report VRAM
	$(PY) -c "import torch; \
		print('CUDA available:', torch.cuda.is_available()); \
		print('Device:', torch.cuda.get_device_name(0) if torch.cuda.is_available() else 'NONE'); \
		print('VRAM GB:', round(torch.cuda.get_device_properties(0).total_memory / 1e9, 1) if torch.cuda.is_available() else 0); \
		print('BF16 supported:', torch.cuda.is_bf16_supported() if torch.cuda.is_available() else False)"

.PHONY: r-check
r-check: ## verify R packages load correctly
	Rscript -e 'for (p in c("scDblFinder","edgeR","CellChat","jsonlite")) \
		cat(sprintf("%-15s %s\n", p, if (requireNamespace(p, quietly=TRUE)) "OK" else "MISSING"))'

# =============================================================================
# Pipeline phases
# =============================================================================

.PHONY: validate
validate: ## Phase 0: manifest validation, sex check, fingerprints (NO compute)
	$(PY) run.py --config $(CONFIG) --step validate

.PHONY: validate-placenta
validate-placenta:
	$(PY) run.py --config config/placenta.yaml --step validate

.PHONY: qc
qc: ## Phase 2: per-sample QC
	$(PY) run.py --config $(CONFIG) --step qc

.PHONY: ambient
ambient: ## Phase 1: CellBender ambient RNA correction (uses .venv-cellbender)
	$(PY) run.py --config $(CONFIG) --step ambient

.PHONY: doublets
doublets: ## Phase 3: scDblFinder (Python invokes R via subprocess)
	$(PY) run.py --config $(CONFIG) --step doublets

.PHONY: integration
integration: ## Phase 5: scVI/scANVI (GPU-heavy, several hours)
	$(PY) run.py --config $(CONFIG) --step integration

.PHONY: annotation
annotation: ## Phase 7: CellTypist + marker-based annotation
	$(PY) run.py --config $(CONFIG) --step annotation

.PHONY: downstream
downstream: ## Phase 8: DE, composition, pathway, communication (all contrasts)
	$(PY) run.py --config $(CONFIG) --step downstream

.PHONY: cross-age
cross-age: ## Phase 8g: persistence / cross-age analysis
	$(PY) run.py --config $(CONFIG) --step cross_age

.PHONY: cross-tissue
cross-tissue: ## Phase 8f: placenta -> brain link (runs after both tissues complete)
	$(PY) run.py --step cross_tissue

# =============================================================================
# Lock file management
# =============================================================================

.PHONY: lock
lock: ## regenerate uv.lock (Python) and renv.lock (R) — commit both to git
	uv lock
	Rscript -e 'renv::snapshot(prompt = FALSE)'

.PHONY: lock-python
lock-python: ## regenerate Python lock file only
	uv lock

.PHONY: lock-r
lock-r: ## regenerate R lock file only
	Rscript -e 'renv::snapshot(prompt = FALSE)'

# =============================================================================
# Cleaning
# =============================================================================

.PHONY: clean
clean: ## remove venvs but keep lockfiles (next `make setup` rebuilds from lock)
	rm -rf .venv .venv-cellbender .renv-cache

.PHONY: clean-results
clean-results: ## remove ALL results (DESTRUCTIVE — confirms first)
	@read -p "Delete results/? Type 'yes' to confirm: " ans && [ "$$ans" = "yes" ] && rm -rf results/

# =============================================================================
# Dev tools
# =============================================================================

.PHONY: lint
lint: ## run ruff linter on src/
	uv run ruff check src/ tests/

.PHONY: format
format: ## run ruff formatter on src/
	uv run ruff format src/ tests/

.PHONY: test
test: ## run pytest
	uv run pytest tests/
