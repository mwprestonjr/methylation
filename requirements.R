# =============================================================================
# requirements.R - R/Bioconductor packages not managed by conda
# Author: GP2 Subtypes and Mechanisms - M.E.
# Date: April 2026
# Updated: Oct 15, 2026
# Description: Installs the packages environment.yml doesn't provide, at the
#              versions in use. Run once, after creating the conda
#              environment:
#                conda env create -f environment.yml
#                conda activate methylation
#                Rscript requirements.R
# =============================================================================

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# Bioconductor packages not on conda. version = "3.20" pins them to the
# Bioconductor release the conda Bioconductor packages come from;
# update = FALSE leaves the conda-installed packages alone
BiocManager::install(c(
    "minfiData",              # needed by maxprobes
    "sva",                    # ComBat batch correction
    "FlowSorted.Blood.EPIC",  # cell type deconvolution (IDOL reference)
    "sesame",                 # per-sample QC stats (qc/2_qc.R)
    "sesameData",             # sesame manifests and reference data
    "methylclock",            # epigenetic clocks (clocks/)
    "methylclockData"         # clock coefficients (clocks/)
), version = "3.20", ask = FALSE, update = FALSE)

# Download sesame reference data (manifests, idat signatures) into the
# ExperimentHub cache; needed once per machine before sesame can read idats
sesameData::sesameDataCache()

# Reinstall preprocessCore without threading (required for Linux/GCP
# environments; see qc/README.md)
BiocManager::install("preprocessCore",
                     configure.args = "--disable-threading",
                     force = TRUE, update = FALSE)

# CRAN and GitHub packages, pinned to the versions in use
install.packages("remotes")

# MatrixEQTL (mqtl/): conda-forge only has 2.3
remotes::install_version("MatrixEQTL", version = "2.4", upgrade = "never")

# Cross-reactive probe lists (R/array_profiles.R); depends on minfiData
remotes::install_github("markgene/maxprobes@c2120dba972e12115280ef274ff80550cee5b264",
                        upgrade = "never")

# DunedinPACE pace of aging (clocks/)
remotes::install_github("danbelsky/DunedinPACE@4b569983543e51d1022aecec9a25e694bb3a336a",
                        upgrade = "never")
