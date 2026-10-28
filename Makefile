SOURCES  := ppmi_p140 psomagen
SOURCE   ?=
ANCESTRY ?=
RSCRIPT  ?= Rscript
LOG_DIR  := logs
STAMP    := $(shell date +%m%d_%H%M)
comma    := ,
MODULES  := qc preprocessing clocks ewas mqtl
MQTL_LOG  = $(subst $(comma),_,$(ANCESTRY))

SHEET_ppmi_p140 := qc/1_build_sample_sheet_ppmi.R
SHEET_psomagen  := qc/1_build_sample_sheet_psomagen.R

log_file = $(LOG_DIR)/$(firstword $(subst /, ,$(1)))/$(2)_$(SOURCE)_$(STAMP).log

define run_r
	@echo "[$$(date +%H:%M)] $(1) $(SOURCE) $(3) -> $(call log_file,$(1),$(2))"
	@$(RSCRIPT) -e 'source("$(1)")' $(SOURCE) $(3) > $(call log_file,$(1),$(2)) 2>&1
endef

# Step 0: check the inputs of the given modules before running them
define check_inputs
	@$(RSCRIPT) R/check_inputs.R $(SOURCE) $(1) $(ANCESTRY)
endef

.NOTPARALLEL:
.PHONY: help all check qc preprocessing clocks ewas mqtl check-source check-ancestry

help:
	@echo "Usage: make <check|qc|preprocessing|clocks|ewas|mqtl|all> SOURCE=<$(subst $() ,|,$(SOURCES))> [ANCESTRY=<label>]"
	@echo "See 'Running the pipeline' in README.md"

all: check-source
	$(call check_inputs,qc$(comma)preprocessing$(comma)clocks$(comma)ewas$(if $(ANCESTRY),$(comma)mqtl))
	@$(MAKE) --no-print-directory qc SOURCE=$(SOURCE) STAMP=$(STAMP)
	@$(MAKE) --no-print-directory preprocessing SOURCE=$(SOURCE) STAMP=$(STAMP)
	@$(MAKE) --no-print-directory clocks SOURCE=$(SOURCE) STAMP=$(STAMP)
	@$(MAKE) --no-print-directory ewas SOURCE=$(SOURCE) STAMP=$(STAMP)
	@if [ -n "$(ANCESTRY)" ]; then \
	  $(MAKE) --no-print-directory mqtl SOURCE=$(SOURCE) ANCESTRY=$(ANCESTRY) STAMP=$(STAMP); \
	fi

qc: check-source
	$(call check_inputs,qc)
	$(call run_r,$(SHEET_$(SOURCE)),1_sample_sheet)
	$(call run_r,qc/2_qc.R,2_qc)
	$(call run_r,qc/3_qc_plots.R,3_qc_plots)
	@echo "[$$(date +%H:%M)] qc complete for $(SOURCE)"

preprocessing: check-source
	$(call check_inputs,preprocessing)
	$(call run_r,preprocessing/1_cell_counts.R,1_cell_counts)
	$(call run_r,preprocessing/2_variation_sources.R,2_variation_sources_raw,raw)
	$(call run_r,preprocessing/3_combat.R,3_combat)
	$(call run_r,preprocessing/2_variation_sources.R,2_variation_sources_combat,combat)
	$(call run_r,preprocessing/4_cpg_annotation.R,4_cpg_annotation)
	@echo "[$$(date +%H:%M)] preprocessing complete for $(SOURCE)"

clocks: check-source
	$(call check_inputs,clocks)
	$(call run_r,clocks/1_estimate.R,1_estimate)
	$(call run_r,clocks/2_compare.R,2_compare)
	$(call run_r,clocks/3_associations.R,3_associations)
	$(call run_r,clocks/4_plots.R,4_plots)
	@echo "[$$(date +%H:%M)] clocks complete for $(SOURCE)"

ewas: check-source
	$(call check_inputs,ewas)
	$(call run_r,ewas/1_case_control.R,1_case_control)
	$(call run_r,ewas/2_age.R,2_age)
	$(call run_r,ewas/3_exposure.R,3_exposure)
	$(call run_r,ewas/4_plots.R,4_plots)
	@echo "[$$(date +%H:%M)] ewas complete for $(SOURCE)"

mqtl: check-source check-ancestry
	$(call check_inputs,mqtl)
	$(call run_r,mqtl/1_mqtl_sample_map.R,$(MQTL_LOG)_1_sample_map,$(ANCESTRY))
	@echo "[$$(date +%H:%M)] mqtl/2_mqtl_genotypes.sh $(SOURCE) $(ANCESTRY) -> $(call log_file,mqtl/,$(MQTL_LOG)_2_genotypes)"
	@bash mqtl/2_mqtl_genotypes.sh $(SOURCE) $(ANCESTRY) > $(call log_file,mqtl/,$(MQTL_LOG)_2_genotypes) 2>&1
	$(call run_r,mqtl/3_mqtl.R,$(MQTL_LOG)_3_mqtl,$(ANCESTRY))
	$(call run_r,mqtl/4_mqtl_permutations.R,$(MQTL_LOG)_4_permutations,$(ANCESTRY))
	$(call run_r,mqtl/5_mqtl_plots.R,$(MQTL_LOG)_5_plots,$(ANCESTRY))
	@echo "[$$(date +%H:%M)] mqtl complete for $(SOURCE) $(ANCESTRY)"

check: check-source
	$(call check_inputs,$(if $(ANCESTRY),all,qc$(comma)preprocessing$(comma)clocks$(comma)ewas))

check-ancestry:
	@if [ -z "$(ANCESTRY)" ]; then \
	  echo "Set ANCESTRY for the mQTL analysis (e.g. make mqtl SOURCE=psomagen ANCESTRY=AFR)"; exit 1; \
	fi

check-source:
	@if [ -z "$(filter $(SOURCE),$(SOURCES))" ] || [ -z "$(SOURCE)" ]; then \
	  echo "Set SOURCE to one of: $(SOURCES)  (e.g. make qc SOURCE=ppmi_p140)"; exit 1; \
	fi
	@mkdir -p $(addprefix $(LOG_DIR)/,$(MODULES))
