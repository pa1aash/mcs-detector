# =====================================================================
# Reproducibility entry points for the mcs-detector code + data release.
#
# TWO PATHS, deliberately separated:
#   * FAST / CACHED (no Geant4, no new simulation): regenerate every figure,
#     run the analytic unit tests, and re-audit the headline numbers against the
#     committed data. Runs in minutes on any laptop.
#   * FULL SIMULATION (needs Geant4 11.3.2 + the built mcs_sim binary + many
#     core-hours): the Monte-Carlo campaign that produces the raw .root output.
#     Raw .root is gitignored/regenerable; the committed reduced outputs are
#     sufficient for the fast path above.
#
# Environment (both paths):  conda env create -f environment.yml && conda activate g4highland
# =====================================================================
PYTHON   ?= python
PDFLATEX ?= pdflatex

.PHONY: all figures test audit campaign clean help

## --- FAST / CACHED PATH (no Geant4, no new simulation) ---------------

all: figures test audit    ## fast reproduction: figures from committed data + unit tests + ledger audit

figures:                   ## regenerate all publication figures from committed data (deterministic geometry; NO Geant4)
	@for t in rectilinear schwarzp gyroid diamond voronoi; do \
	  raw=data/geom_stats/voxel/$${t}_f40_c2.5_camp_vox.raw; \
	  if test ! -f $$raw; then echo "== voxels $$t"; $(PYTHON) geom/make_campaign_voxels.py $$t 0.40 || exit 1; fi; \
	done
	@cd figs/src && for f in fig_geometry.py fig_pt.py fig_validation.py \
	  fig_width_invariance.py fig_neff_collapse_3energy.py fig_boundary.py \
	  fig_rocking.py fig_impact.py fig_affordability.py fig_systematics.py; do \
	  echo "== $$f"; $(PYTHON) $$f || exit 1; done
	@if command -v $(PDFLATEX) >/dev/null 2>&1; then \
	  echo "== fig_mechanism.tex"; cd figs/src && \
	  $(PDFLATEX) -interaction=nonstopmode -halt-on-error -output-directory=.. fig_mechanism.tex >/dev/null && \
	  rm -f ../fig_mechanism.aux ../fig_mechanism.log; \
	else echo "== fig_mechanism.tex skipped (no pdflatex; the committed figs/fig_mechanism.pdf is used)"; fi

test:                      ## analytic theory-identity unit tests (no data needed)
	pytest -q analysis/lib/test_theory.py

audit:                     ## independently re-verify every headline number against committed data (claims ledger)
	$(PYTHON) tools/check_ledger.py

## --- FULL SIMULATION PATH (Geant4 11.3.2 + compute box; ~310 core-hours) ---

campaign:                  ## FULL Monte-Carlo campaign (needs the built sim/build/mcs_sim + g4highland; produces raw .root)
	@echo "FULL SIMULATION PATH --- needs Geant4 11.3.2 and the built mcs_sim binary."
	@echo "  1. build the engine:   make -C sim            (CMake + Geant4)"
	@echo "  2. run the campaign:   bash sim/run_campaign.sh   (Highland/floor validation)"
	@echo "                         python sim/campaign.py     (resumable topology x fill x energy grid)"
	@echo "  3. high-stat reruns:   bash sim/escalate_500.sh ; bash sim/escalate_voronoi_3e7.sh"
	@echo "  Raw .root lands in data/runs/ (gitignored). Then re-run the analysis/ scripts to"
	@echo "  regenerate data/analysis/*.json, and 'make figures' + 'make audit' as above."
	@echo "  Total ~310 core-hours; use a compute-optimized instance, not a laptop."

clean:                     ## remove regenerated figure build artifacts
	rm -f figs/fig_mechanism.aux figs/fig_mechanism.log

help:                      ## list targets
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | sort | awk 'BEGIN{FS=":.*?## "}{printf "  \033[36m%-10s\033[0m %s\n", $$1, $$2}'
