# =============================================================================
# Methylation pipeline
# Usage:   make <target> SOURCE=<data source>    (data sources: see config.R)
#   make qc            SOURCE=ppmi_p140   sample sheet, QC/normalization, QC plots
#   make preprocessing SOURCE=ppmi_p140   cell counts, sources of variation, ComBat
#   make all           SOURCE=ppmi_p140   qc then preprocessing
# Run from the repository root with the methylation conda environment active.
# Each step logs to logs/<step>_<SOURCE>_<timestamp>.log, and make stops at the
# first step that fails. For long runs:
#   nohup make all SOURCE=ppmi_p140 > logs/make_ppmi_p140.log 2>&1 &
# =============================================================================

SOURCES := ppmi_p140 psomagen
SOURCE  ?=
RSCRIPT ?= Rscript
LOG_DIR := logs
STAMP   := $(shell date +%m%d_%H%M)

# Sample sheet script for each data source
SHEET_ppmi_p140 := qc/01_build_sample_sheet_ppmi.R
SHEET_psomagen  := qc/01_build_sample_sheet_psomagen.R

# Run an R script with the data source (and any extra arguments) and log it.
# Scripts are run through source(), which reads the whole file before running
# it, so editing a script during a run can't break the run
#   $(call run_r,<script>,<log name>,<extra arguments>)
define run_r
	@echo "[$$(date +%H:%M)] $(1) $(SOURCE) $(3) -> $(LOG_DIR)/$(2)_$(SOURCE)_$(STAMP).log"
	@$(RSCRIPT) -e 'source("$(1)")' $(SOURCE) $(3) > $(LOG_DIR)/$(2)_$(SOURCE)_$(STAMP).log 2>&1
endef

# Steps run one after another, never in parallel (they share memory-heavy data)
.NOTPARALLEL:
.PHONY: help all qc preprocessing check-source

# Default target: show usage instead of starting a run
help:
	@sed -n '2,10p' Makefile | sed 's/^# \{0,1\}//'

all: check-source
	@$(MAKE) --no-print-directory qc SOURCE=$(SOURCE) STAMP=$(STAMP)
	@$(MAKE) --no-print-directory preprocessing SOURCE=$(SOURCE) STAMP=$(STAMP)

qc: check-source
	$(call run_r,$(SHEET_$(SOURCE)),01_sample_sheet)
	$(call run_r,qc/02_qc.R,02_qc)
	$(call run_r,qc/03_qc_plots.R,03_qc_plots)
	@echo "[$$(date +%H:%M)] qc complete for $(SOURCE)"

preprocessing: check-source
	$(call run_r,preprocessing/1_cell_counts.R,1_cell_counts)
	$(call run_r,preprocessing/2_variation_sources.R,2_variation_sources_raw,raw)
	$(call run_r,preprocessing/3_combat.R,3_combat)
	$(call run_r,preprocessing/2_variation_sources.R,2_variation_sources_combat,combat)
	@echo "[$$(date +%H:%M)] preprocessing complete for $(SOURCE)"

check-source:
	@if [ -z "$(filter $(SOURCE),$(SOURCES))" ] || [ -z "$(SOURCE)" ]; then \
	  echo "Set SOURCE to one of: $(SOURCES)  (e.g. make qc SOURCE=ppmi_p140)"; exit 1; \
	fi
	@mkdir -p $(LOG_DIR)
