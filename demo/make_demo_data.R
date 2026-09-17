# =============================================================================
# Demo data generator
# =============================================================================
#
# Creates a small SIMULATED dataset with the same structure as the real BMCA
# objects, so that run_demo.R can be executed by anyone without access to the
# controlled-access patient data.
#
# The simulated data is not derived from any patient sample. It reproduces the
# schema (column names, factor levels, object shapes) and a plausible signal,
# so the demo exercises the same code paths as the manuscript analyses.
#
# Output: demo/data/  (~3 MB, written in a few seconds)
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(tidyr)
})

set.seed(42)

dir.create("demo/data", recursive = TRUE, showWarnings = FALSE)

# ---- Design -----------------------------------------------------------------
n_patients   <- 24
cell_types   <- c("Malignant", "Myeloid", "T_NK", "B_Plasma", "Endothelial",
                  "Fibroblast", "Pericyte", "Astrocyte", "Oligodendrocyte",
                  "Neuron")
primary_sites <- c("melanoma", "NSCLC", "breast", "colorectal")

patients <- tibble(
  patient      = sprintf("P%02d", seq_len(n_patients)),
  primary_site = rep(primary_sites, length.out = n_patients),
  location     = rep(c("Intracranial", "Extracranial"), each = n_patients / 2),
  seq_tech     = rep(c("sn", "sc"), length.out = n_patients),
  study        = rep(c("demo_study_1", "demo_study_2"), length.out = n_patients)
)

# ---- Per-cell metadata ------------------------------------------------------
# Brain-resident cell types are only recovered in single-nucleus intracranial
# samples, mirroring the real atlas; fibroblasts are enriched extracranially.
compose <- function(location, seq_tech) {
  w <- c(Malignant = 40, Myeloid = 16, T_NK = 12, B_Plasma = 4,
         Endothelial = 7, Fibroblast = 8, Pericyte = 4,
         Astrocyte = 0, Oligodendrocyte = 0, Neuron = 0)
  if (location == "Intracranial") {
    w["Fibroblast"] <- 3
    w["Myeloid"]    <- 22
    if (seq_tech == "sn") {
      w["Astrocyte"] <- 6; w["Oligodendrocyte"] <- 5; w["Neuron"] <- 3
    }
  }
  w / sum(w)
}

metadata <- patients %>%
  rowwise() %>%
  do({
    p <- .
    n <- sample(600:1400, 1)
    w <- compose(p$location, p$seq_tech)
    w <- w * runif(length(w), 0.7, 1.3)          # per-patient variability
    ct <- sample(cell_types, n, replace = TRUE, prob = w / sum(w))
    tibble(
      CellID       = sprintf("%s_cell%05d", p$patient, seq_len(n)),
      patient      = p$patient,
      study        = p$study,
      primary_site = p$primary_site,
      location     = p$location,
      seq_tech     = p$seq_tech,
      cell_type    = ct
    )
  }) %>%
  ungroup()

# a coarse cell-state annotation, as in metadata_states_*.RDS
state_pool <- c("Hypoxia", "Stress", "Cell_Cycle", "Interferon", "Inflammatory")
metadata <- metadata %>%
  mutate(state = sample(state_pool, n(), replace = TRUE))

saveRDS(metadata, "demo/data/demo_metadata.rds")

# ---- Patient x cell-type abundance matrix -----------------------------------
celltype_patient_matrix <- metadata %>%
  count(patient, cell_type) %>%
  pivot_wider(names_from = cell_type, values_from = n, values_fill = 0)

saveRDS(celltype_patient_matrix, "demo/data/demo_celltype_patient_matrix.rds")

# ---- Pseudobulk expression matrix (fibroblasts) -----------------------------
# Structure matches BMCA/data/UMI_pb_<cell_type>.RDS: a named list of
# gene x sample count matrices, one per cohort.
n_genes <- 2000
genes   <- sprintf("GENE%04d", seq_len(n_genes))

samples <- patients$patient
base_mu <- 2^rnorm(n_genes, mean = 5, sd = 2)

# 60 genes are genuinely higher intracranially, 60 extracranially
up_ic <- sample(n_genes, 60)
up_ec <- sample(setdiff(seq_len(n_genes), up_ic), 60)

counts <- sapply(seq_along(samples), function(j) {
  mu <- base_mu
  if (patients$location[j] == "Intracranial") {
    mu[up_ic] <- mu[up_ic] * 4
  } else {
    mu[up_ec] <- mu[up_ec] * 4
  }
  mu <- mu * runif(n_genes, 0.6, 1.6)            # gene-level noise
  rpois(n_genes, lambda = mu * runif(1, 0.8, 1.2))
})
dimnames(counts) <- list(genes, samples)

saveRDS(list(lin_sn = counts), "demo/data/demo_UMI_pb_Fibroblast.rds")

# ---- Ground truth, for checking the demo output -----------------------------
saveRDS(
  list(up_intracranial = genes[up_ic], up_extracranial = genes[up_ec]),
  "demo/data/demo_truth.rds"
)

cat("Demo data written to demo/data/\n")
cat("  cells:    ", nrow(metadata), "\n")
cat("  patients: ", n_patients, "\n")
cat("  genes:    ", n_genes, "\n")
