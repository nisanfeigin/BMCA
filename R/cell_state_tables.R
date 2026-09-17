# =============================================================================
# Cell-state and cell-type abundance tables (shared upstream step)
# =============================================================================
#
# Builds the per-patient cell-type and cell-state abundance matrices used by the
# figure scripts, and compares lymphocyte meta-programs against published
# lymphocyte subset signatures by Jaccard similarity.
#
# Part of: A brain metastasis cell atlas reveals multicellular adaptation to
#          the neural microenvironment (Feigin et al.)
#
# Run from the project root (see README.md). All paths below are relative to it.
#
# Inputs:
#   BMCA/csv/all_MPs.xlsx
#   BMCA/csv/lymphocyte_subset_signatures_merged.rds
#   BMCA/data/celltype_patient_matrix_proccesed.RDS
#   BMCA/data/ecotypes_tib_proccesed.RDS
#   BMCA/data/metadata_all_studies.rds
#
# Outputs:
#   BMCA/data/  (3 files)
#
# Original filename: states_tibs.R
# =============================================================================

# ---- Packages ---------------------------------------------------------------
pkgs <- c(
  "dendextend", "dplyr", "ggdendro", "ggplot2", "ggrepel", "Hmisc",
  "janitor", "Matrix", "purrr", "RColorBrewer", "readxl", "reshape2",
  "scales", "tidyr", "tidyverse", "viridis"
)
invisible(lapply(pkgs, library, character.only = TRUE))

# ---- Project paths ----------------------------------------------------------
source("R/config.R")

# =============================================================================
sheets <- c("GD", "CD8", "CD4")

mps <- map(sheets, ~ read_excel("BMCA/csv/all_MPs.xlsx", sheet = .x))
names(mps) <- sheets
ext_sigs <- readRDS("BMCA/csv/lymphocyte_subset_signatures_merged.rds")

jaccard <- function(a, b) {
  a <- unique(a[!is.na(a)]); b <- unique(b[!is.na(b)])
  length(intersect(a, b)) / length(union(a, b))
}

ext_long <- tibble(sig = names(ext_sigs),
                   genes = map(ext_sigs, ~ unique(na.omit(.x))))

make_heatmap <- function(tib, title) {
  mp_long <- tibble(mp = names(tib),
                    genes = map(tib, ~ unique(na.omit(.x))))

  m <- outer(
    seq_len(nrow(mp_long)), seq_len(nrow(ext_long)),
    Vectorize(function(i, j) jaccard(mp_long$genes[[i]], ext_long$genes[[j]]))
  )
  rownames(m) <- mp_long$mp
  colnames(m) <- ext_long$sig

  # drop columns whose max Jaccard is below 0.1
  m <- m[, apply(m, 2, max) >= 0.1, drop = FALSE]

  # order both axes by hclust
  row_ord <- rownames(m)[hclust(dist(m))$order]
  col_ord <- colnames(m)[hclust(dist(t(m)))$order]

  as.data.frame(m) |>
    rownames_to_column("mp") |>
    pivot_longer(-mp, names_to = "sig", values_to = "jaccard") |>
    mutate(mp  = factor(mp,  levels = row_ord),
           sig = factor(sig, levels = col_ord)) |>
    ggplot(aes(sig, mp, fill = jaccard)) +
    geom_tile(color = "grey90", linewidth = 0.1) +
    scale_fill_viridis_c(option = "magma", name = "Jaccard") +
    labs(title = title, x = "External signature", y = "Meta-program") +
    theme_minimal(base_size = 8) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 4),
          panel.grid = element_blank())
}
p_gd  <- make_heatmap(mps$GD,  "GD")
p_cd8 <- make_heatmap(mps$CD8, "CD8")
p_cd4 <- make_heatmap(mps$CD4, "CD4")

p_gd; p_cd8; p_cd4

cts <- c("Endothelial", "Fibroblast", "Pericyte", "B_Plasma", "T_NK", "Myeloid",
         "Malignant",       "Oligodendrocyte", "Neuron",          "Astrocyte")
metadata_all <- as_tibble(readRDS("BMCA/data/metadata_all_studies.rds")) %>%
  filter(!(study == "anoop")) %>%
  filter(cell_type %in% cts)

tab_tib <- metadata_all %>%
  tabyl(patient, cell_type) %>%
  as_tibble()
tab_tib_t <- metadata_all %>%
  tabyl(patient, t_sub_anno, ) %>%
  as_tibble() %>%
  dplyr::select(-Other, -UA,  -NA_, -DP)
tab_tib <- tab_tib %>% left_join(tab_tib_t, by = "patient")
saveRDS(tab_tib, "BMCA/data/celltype_patient_matrix_new.RDS")

# -----------------------------------------------------------------------------
celltype_patient_matrix <- tab_tib
metadata_all_sc <- metadata_all %>%
  filter(seq_tech == "sc") %>%
  dplyr::select(patient) %>%
  distinct()
celltype_patient_matrix <- unique(celltype_patient_matrix)
celltype_patient_matrix <- celltype_patient_matrix %>%
  mutate(Neuron = if_else(patient %in% metadata_all_sc$patient, 0L, Neuron)) %>%
  mutate(Astrocyte = if_else(patient %in% metadata_all_sc$patient, 0L, Astrocyte))

celltype_patient_matrix_mal <- celltype_patient_matrix %>%
  rowwise() %>%
  mutate(across(-patient, ~ .x / sum(c_across(-patient)) * 100)) %>%
  ungroup() %>%
  dplyr::select(patient, Malignant)

celltype_patient_matrix_tme<- celltype_patient_matrix %>%
  dplyr::select(-Malignant) %>%
  rowwise() %>%
  ungroup() %>%
  dplyr::select(-CD8, -CD4, -NK, -GD) %>%
  rowwise() %>%
  mutate(across(-patient, ~ .x / sum(c_across(-patient)) * 100)) %>%
  ungroup()

celltype_patient_matrix_tme

rowSums(celltype_patient_matrix_tme[,-1])

celltype_patient_matrix_T_NK<- celltype_patient_matrix %>%
  dplyr::select(patient, CD4, CD8, GD, NK) %>%
  rowwise() %>%
  mutate(across(-patient, ~ .x / sum(c_across(-patient)) * 100)) %>%
  ungroup()

rowSums(celltype_patient_matrix_T_NK[,-1])

celltype_patient_matrix <- celltype_patient_matrix_mal %>%
  left_join(celltype_patient_matrix_tme, by = "patient") %>%
  left_join(celltype_patient_matrix_T_NK, by = "patient")

saveRDS(celltype_patient_matrix, "BMCA/data/celltype_patient_matrix_proccesed.RDS")

celltype_patient_matrix <- readRDS("BMCA/data/celltype_patient_matrix_proccesed.RDS")

summary_tib <- metadata_all %>%
  count(patient, state) %>%                 # count number of cells per patient/state
  pivot_wider(names_from = state,           # make one column per state
              values_from = n,              # values = counts
              values_fill = 0)

summary_tib <- summary_tib %>%
  dplyr::select(
    -matches("LQ|MT|RP|NA|Tech|Mesenchymal|Doublet|not_assigned|Not_assigned|Not_Assigned", ignore.case = FALSE),
    -matches("^(?!Astrocyte|Oligodendrocyte|Neuron).*UA_\\d+$", perl = TRUE)
  )
# first, identify the patients that are Extracranial
extracranial_patients <- metadata_all %>%
  filter(location == "Extracranial") %>%
  distinct(patient) %>%
  pull(patient)

# now set relevant columns to NA for those patients
summary_tib <- summary_tib %>%
  mutate(across(
    matches("^(Neuron|Oligodendrocyte|Astrocyte)"),
    ~ if_else(patient %in% extracranial_patients, NA, .x)
  ))

summary_tib <- summary_tib %>%
  filter(patient %in% celltype_patient_matrix$patient)

cell_types <- colnames(summary_tib)[colnames(summary_tib) != "patient"]
cell_types <- unique(sub("_.*", "", cell_types))
cell_types[cell_types == "B"] <- "B_Plasma"
cell_types

# Loop through each cell type
for (ct in cell_types) {
  print(ct)

  # Find patients where the value in celltype_patient_matrix is < 3
  patients_to_mask <- celltype_patient_matrix %>%
    filter(.data[[ct]] < 3) %>%
    pull(patient)

  # Mask all columns in summary_tib that belong to this cell type for these patients
  summary_tib <- summary_tib %>%
    mutate(across(
      starts_with(ct),
      ~ if_else(patient %in% patients_to_mask, NA_integer_, .x)
    ))
}

# Create a list to store tibbles per cell type
summary_list <- lapply(cell_types, function(ct) {
  summary_tib %>%
    dplyr::select(patient, starts_with(ct))
})

# Name the list elements with the cell types
names(summary_list) <- cell_types

summary_list <- lapply(summary_list, function(sub_tib) {
  sub_tib %>%
    rowwise() %>%
    mutate(across(-patient, ~ if(sum(c_across(-patient), na.rm = TRUE) < 50) NA_real_ else .x)) %>%
    ungroup()
  })


 summary_list <- lapply(summary_list, function(sub_tib) {
  sub_tib %>%
    rowwise() %>%
    mutate(across(-patient, ~ .x / sum(c_across(-patient), na.rm = TRUE) * 100)) %>%
    ungroup()
})
names(summary_list)


for(ct_name in names(summary_list)) {

  sub_tib <- summary_list[[ct_name]]

  # Compute row sums for each patient, ignoring NAs
  plot_data <- sub_tib %>%
    rowwise() %>%
    mutate(row_sum = sum(c_across(-patient), na.rm = TRUE)) %>%
    ungroup()

  # Create bar plot
  p <- ggplot(plot_data, aes(x = patient, y = row_sum)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    labs(title = paste("Total", ct_name, "per patient"),
         x = "Patient",
         y = "Total") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))

  print(p)
}

patients_md <- metadata_all %>%
  dplyr::select(patient, primary_site, location, study, seq_tech, subtype, Donor) %>%
  distinct()


# Left join all tibbles in summary_list by "patient"
combined_summary <- purrr::reduce(summary_list, left_join, by = "patient")
final_matrix <- combined_summary %>%
  left_join(patients_md, by = "patient") %>%
  filter(patient %in% celltype_patient_matrix$patient)

final_matrix <- unique(final_matrix)

patients <- intersect(final_matrix$patient, celltype_patient_matrix$patient)

final_matrix <- final_matrix %>%
  filter(patient %in% patients)

celltype_patient_matrix <- celltype_patient_matrix %>%
  filter(patient %in% patients)

ecotypes_tib <- celltype_patient_matrix %>%
  left_join(final_matrix, by = "patient") %>%
  filter(patient %in% patients) %>%
  distinct()

major_cells <- c("B_Plasma", "Endothelial", "Fibroblast", "Pericyte", "Myeloid",
                 "Neuron", "Astrocyte", "Oligodendrocyte")

for (cell in major_cells) {
  # Find all columns that start with the cell type name and underscore
  sub_cols <- grep(paste0("^", cell, "_"), names(ecotypes_tib), value = TRUE)

  # Set subtypes to NA if major cell type < 3
  ecotypes_tib <- ecotypes_tib %>%
    mutate(across(all_of(sub_cols), ~ if_else(.data[[cell]] < 3, NA_real_, .x)))
}


# Pivot T/NK columns to long format
tnk_long <- ecotypes_tib %>%
  dplyr::select(patient, CD4, CD8, GD, NK) %>%
  pivot_longer(cols = c(CD4, CD8, GD, NK),
               names_to = "CellType",
               values_to = "Value")

# Calculate counts above thresholds for annotation
annotations <- tnk_long %>%
  group_by(CellType) %>%
  dplyr::summarise(
    above5 = sum(Value > 5, na.rm = TRUE),
    above10 = sum(Value > 10, na.rm = TRUE)
  )

# Plot histograms with annotations
ggplot(tnk_long, aes(x = Value, fill = CellType)) +
  geom_histogram(binwidth = 5, color = "white", alpha = 0.7) +
  facet_wrap(~CellType, scales = "free_x") +
  geom_text(
    data = annotations,
    aes(
      x = Inf, y = Inf,
      label = paste0(">", 5, ": ", above5, "\n>", 10, ": ", above10)
    ),
    hjust = 1.1, vjust = 1.1,
    inherit.aes = FALSE
  ) +
  labs(title = "Distribution of T/NK Subtypes with Threshold Counts",
       x = "Percentage",
       y = "Number of Patients") +
  theme_minimal()

T_NK_cells <- c("CD4", "CD8", "GD", "NK")

for (cell in T_NK_cells) {
  # Find all columns that start with the cell type name and underscore
  sub_cols <- grep(paste0("^", cell, "_"), names(ecotypes_tib), value = TRUE)

  # Set subtypes to NA if major cell type < 3
  ecotypes_tib <- ecotypes_tib %>%
    mutate(across(all_of(sub_cols), ~ if_else(.data[[cell]] < 3, NA_real_, .x)))
}

ecotypes_tib <- ecotypes_tib %>%
  mutate(
    across(
      .cols = matches("Microglia|Neuron_|Oligodendrocyte|Astrocyte"),
      .fns = ~ if_else(location == "Extracranial", NA_real_, .)
    ),
    across(
      .cols = matches("Epithelial"),
      .fns = ~ if_else(location == "Intracranial", NA_real_, .)
    ),
    across(
      .cols = matches("Skin"),
      .fns = ~ if_else(primary_site != "melanoma", NA_real_, .)
    )
  )

ecotypes_tib <- ecotypes_tib %>%

  # 1. Remove cols with "oubl", "Tech", or "tech"
  dplyr::select(-matches("oubl|[Tt]ech")) %>%

  # 2. Sum Malignant_Stress_2 into Malignant_Stress, remove _2
  mutate(Malignant_Stress = Malignant_Stress + Malignant_Stress_2) %>%
  dplyr::select(-Malignant_Stress_2) %>%
  # 2. Sum B_Plasma_Stress_2 into B_Plasma_Stress, remove _2
  mutate(B_Plasma_Stress = B_Plasma_Stress_1 + B_Plasma_Stress_2) %>%
  dplyr::select(-B_Plasma_Stress_2, -B_Plasma_Stress_1) %>%

  # 3. Sum B_Plasma_CC_1 + B_Plasma_CC_2 into B_Plasma_CC, remove both
  mutate(B_Plasma_CC = B_Plasma_CC + B_Plasma_CC_1 + B_Plasma_CC_2) %>%
  dplyr::select(-B_Plasma_CC_1, -B_Plasma_CC_2) %>%

  # 4. Sum B_Plasma_Plasma_1 + B_Plasma_Plasma_2 into B_Plasma_Plasma, remove both
  mutate(B_Plasma_Plasma = B_Plasma_Plasma + B_Plasma_Plasma_1 + B_Plasma_Plasma_2) %>%
  dplyr::select(-B_Plasma_Plasma_1, -B_Plasma_Plasma_2) %>%

  # 5. Sum all Fibroblast_Pericyte_like* into Fibroblast_Pericyte_like
  mutate(Fibroblast_Pericyte_like = rowSums(dplyr::select(., matches("Fibroblast_Pericyte_like")))) %>%
  dplyr::select(-matches("Fibroblast_Pericyte_like_")) %>%

  # 6. Sum Malignant_Hypoxia_2 into Malignant_Hypoxia, remove _2
  mutate(Malignant_Hypoxia = Malignant_Hypoxia + Malignant_Hypoxia_2) %>%
  dplyr::select(-Malignant_Hypoxia_2) %>%

  # 7. Sum Malignant_Respiration_1 + _2 into Malignant_Respiration, remove both
  mutate(Malignant_Respiration = rowSums(dplyr::select(., matches("Malignant_Respiration")))) %>%
  dplyr::select(-Malignant_Respiration_1, -Malignant_Respiration_2) %>%

  # 8. Sum Myeloid_Macrophage_1 + _2 into Myeloid_Macrophage, remove both
  mutate(Myeloid_Macrophage = rowSums(dplyr::select(., matches("Myeloid_Macrophage")))) %>%
  dplyr::select(-Myeloid_Macrophage_1, -Myeloid_Macrophage_2) %>%

  # 9. Add CD8_CT_total as sum of all CD8_CT* cols
  mutate(CD8_CT = rowSums(dplyr::select(., matches("CD8_CT")))) %>%
  dplyr::select(-CD8_CT1, -CD8_CT2) %>%

  # BONUS: Fix Pericyte_stress vs Pericyte_Stress (capitalization inconsistency)
  mutate(Pericyte_Stress = Pericyte_Stress + Pericyte_stress) %>%
  dplyr::select(-Pericyte_stress) %>%

  # BONUS: Fix Fibroblast_stress vs Fibroblast_Stress
  mutate(Fibroblast_Stress = Fibroblast_Stress + Fibroblast_stress) %>%
  dplyr::select(-Fibroblast_stress) %>%

  # BONUS: Sum Myeloid_INFR_1 + _2 into Myeloid_INFR
  mutate(Myeloid_INFR = rowSums(dplyr::select(., matches("Myeloid_INFR")))) %>%
  dplyr::select(-Myeloid_INFR_1, -Myeloid_INFR_2)

saveRDS(ecotypes_tib, "BMCA/data/ecotypes_tib_proccesed.RDS")

# -----------------------------------------------------------------------------
ecotypes_tib <- readRDS("BMCA/data/ecotypes_tib_proccesed.RDS")

can_types <- "all"
primary_sites <-  c("colorectal", "breast", "NSCLC", "melanoma")

ecotypes_tib <- ecotypes_tib %>%
  filter(primary_site %in% primary_sites)

mean_threshold <- 1

ecotypes_states <- ecotypes_tib %>%
  dplyr::select(patient,
    where(~ is.numeric(.x) && sum(!is.na(.x)) >= 5)
  ) %>%
  dplyr::select(patient,
    where(~ is.numeric(.x) && mean(.x, na.rm = TRUE) >= mean_threshold)
  )

metadata_p <-  metadata_all %>%
  dplyr::select(patient, primary_site, location, study, seq_tech, subtype, Donor) %>%
  distinct()
ecotypes_states <- ecotypes_states %>% left_join(metadata_p, by = "patient")
metadata_states <- ecotypes_states %>%
  dplyr::select(patient, primary_site, location, study, seq_tech, subtype, Donor)

metadata_states <- metadata_states %>%
  filter(patient %in% ecotypes_states$patient)
ecotypes_states <- ecotypes_states %>%
  dplyr::select(-study, -location, -seq_tech, -subtype, -primary_site, -Donor)
saveRDS(ecotypes_states, "BMCA/data/ecotypes_tib_proccesed_filtered.RDS")
saveRDS(metadata_states, "BMCA/data/metadata_states_proccesed_filtered.RDS")

# -----------------------------------------------------------------------------

