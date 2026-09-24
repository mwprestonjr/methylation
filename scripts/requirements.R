# =============================================================================
# requirements.R - R/Bioconductor packages not managed by conda
# Author: GP2 Subtypes and Mechanisms - M.E.
# Date: April 2026
# =============================================================================

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c(
    "minfi",
    "minfiData",
    "IlluminaHumanMethylationEPICanno.ilm10b4.hg19",
    "IlluminaHumanMethylationEPICmanifest",
    "IlluminaHumanMethylationEPICv2anno.20a1.hg38",   # EPICv2 annotation (Scripts 02, 03)
    "IlluminaHumanMethylationEPICv2manifest",          # EPICv2 manifest, loaded by minfi
    "DMRcate",
    "limma",
    "bumphunter",
    "sva",                    # ComBat batch correction
    "FlowSorted.Blood.EPIC"   # cell type deconvolution
), ask = FALSE)

# Install preprocessCore without threading (required for Linux/GCP environments)
BiocManager::install("preprocessCore",
                     configure.args = "--disable-threading",
                     force = TRUE)

# CRAN packages
install.packages("remotes")

# GitHub package for cross-reactive probe removal - depends on minfiData
remotes::install_github("markgene/maxprobes")
# Methylation clocks (Horvath, Hannum, PhenoAge, etc.) and pace of aging
# Run from an activated env (compilers must be on PATH). Deps needing cmake/libuv
# come from conda first:
#   mamba install -n methylation -c conda-forge r-nloptr r-httpuv r-deriv r-doby r-lme4
BiocManager::install(c("methylclock", "methylclockData"), ask = FALSE)
remotes::install_github("danbelsky/DunedinPACE")
