# =============================================================================
# Array profiles
# Author: GP2 Subtypes and Mechanisms - M.E., M.P.
# Date: Sept 29, 2026
# Description: Everything in the pipeline that differs by array type. Script 01
#              detects the array per chip (detect_array); Script 02 looks up
#              the profile for the array in the sample sheet. Everything else
#              (bead counts, sex probes, SNP probes, getSex, Funnorm) reads the
#              annotation attached to the rgSet and works for any array here.
#              To support another array, add a detect_array range and a profile.
# =============================================================================

# Each profile has:
#   anno_pkg   minfi annotation package to load
#   annotation value for annotation(rgSet); set explicitly because minfi
#              can't detect EPICv2 by itself
#   xreactive  function(ann) returning cross-reactive probe names, where ann
#              is getAnnotation(rgSet)
ARRAY_PROFILES <- list(
  EPICv1 = list(
    anno_pkg   = "IlluminaHumanMethylationEPICanno.ilm10b4.hg19",
    annotation = c(array = "IlluminaHumanMethylationEPIC",
                   annotation = "ilm10b4.hg19"),
    xreactive  = function(ann) {
      unlist(maxprobes::xreactive_probes(array_type = "EPIC"))
    }
  ),
  EPICv2 = list(
    anno_pkg   = "IlluminaHumanMethylationEPICv2anno.20a1.hg38",
    annotation = c(array = "IlluminaHumanMethylationEPICv2",
                   annotation = "20a1.hg38"),
    # maxprobes only has an EPICv1 list; map it to EPICv2 names via EPICv1_Loci.
    # NOTE: probes new on EPICv2 (no v1 equivalent) are not screened
    xreactive  = function(ann) {
      cross_reactive_v1 <- unlist(maxprobes::xreactive_probes(array_type = "EPIC"))
      ann$Name[ann$EPICv1_Loci %in% cross_reactive_v1]
    }
  )
)

# Array type of a sample from the number of bead types in its green idat
# (EPICv2: 1,105,209; EPICv1: 1,051,815-1,052,641; 450k: 622,399). Anything
# else (e.g. a genotyping chip) is "Unknown"
detect_array <- function(basename) {
  n_beads <- nrow(illuminaio::readIDAT(paste0(basename, "_Grn.idat"))$Quants)
  dplyr::case_when(dplyr::between(n_beads, 1100000, 1110000) ~ "EPICv2",
                   dplyr::between(n_beads, 1045000, 1060000) ~ "EPICv1",
                   dplyr::between(n_beads,  615000,  625000) ~ "450k",
                   TRUE                                      ~ "Unknown")
}
