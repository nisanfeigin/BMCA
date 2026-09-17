# =============================================================================
# Record the exact software environment used to produce the figures
# =============================================================================
#
# Run this once, in the same R session / on the same machine you use for the
# final figure run, and commit the resulting files. Nature's code checklist
# asks for dependency versions and the versions the code was tested on.
#
#   Rscript R/capture_session_info.R
# =============================================================================

pkgs <- c(
  "biomaRt", "cowplot", "data.table", "dendextend", "doParallel", "dplyr",
  "enrichR", "forcats", "foreach", "ggdendro", "ggh4x", "ggnewscale",
  "ggplot2", "ggpubr", "ggraph", "ggrastr", "ggrepel", "Hmisc", "igraph",
  "janitor", "lubridate", "Matrix", "NMF", "openxlsx", "patchwork",
  "pheatmap", "purrr", "RColorBrewer", "readr", "readxl", "reshape2",
  "rlang", "rstatix", "scalop", "scales", "stringr", "survival", "survminer",
  "tibble", "tidyr", "tidytext", "tidyverse", "viridis"
)

invisible(lapply(pkgs, function(p) {
  suppressPackageStartupMessages(
    try(library(p, character.only = TRUE), silent = TRUE)
  )
}))

writeLines(capture.output(sessionInfo()), "sessionInfo.txt")

versions <- vapply(pkgs, function(p) {
  tryCatch(as.character(packageVersion(p)), error = function(e) "NOT INSTALLED")
}, character(1))

write.csv(
  data.frame(package = pkgs, version = unname(versions)),
  "package_versions.csv", row.names = FALSE
)

cat("Wrote sessionInfo.txt and package_versions.csv\n")
