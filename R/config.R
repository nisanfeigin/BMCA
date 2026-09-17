# =============================================================================
# Shared configuration for the BMCA analysis scripts
# =============================================================================
#
# Every script in R/ is written with paths relative to the project root and
# begins with source("R/config.R"). This file replaces the machine-specific
# setwd() calls that were used during development.
#
# Set the project root in one of three ways (checked in this order):
#   1. open the BMCA.Rproj file in RStudio, or
#   2. export BMCA_ROOT=/path/to/project before starting R, or
#   3. start R with the project root as the working directory.
# =============================================================================

.bmca_root <- Sys.getenv("BMCA_ROOT", unset = getwd())
if (.bmca_root != getwd()) setwd(.bmca_root)

# ---- Directory layout -------------------------------------------------------
DATA_DIR  <- "BMCA/data"   # .RDS objects produced by the preprocessing pipeline
CSV_DIR   <- "BMCA/csv"    # gene signatures, published reference sets, clinical tables
PLOTS_DIR <- "BMCA/plots"  # figure panels written by the scripts

# ---- Sanity check -----------------------------------------------------------
.missing <- c(DATA_DIR, CSV_DIR)[!dir.exists(c(DATA_DIR, CSV_DIR))]
if (length(.missing) > 0) {
  stop(
    "Cannot find: ", paste(.missing, collapse = ", "), "\n",
    "Working directory is: ", getwd(), "\n",
    "Set the project root as described in R/config.R (see README.md).",
    call. = FALSE
  )
}

# ---- Output directories -----------------------------------------------------
for (.d in file.path(PLOTS_DIR, c("", "fig1", "fig3", "fig4", "fig5"))) {
  if (!dir.exists(.d)) dir.create(.d, recursive = TRUE)
}
rm(.d, .missing, .bmca_root)

# ---- Figure export defaults -------------------------------------------------
# Raster panels are exported at 300 dpi; vector panels are written as PDF.
FIG_DPI <- 300
