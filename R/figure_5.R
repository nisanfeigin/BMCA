# =============================================================================
# Figure 5 - Multicellular ecotypes, co-occurrence networks and clinical associations
# =============================================================================
#
# Co-occurrence network analysis of ecotype/state features split by melanoma vs.
# carcinoma, brain-adaptation enrichment across cell types, associations with sex
# and age, and Kaplan-Meier survival analysis (OS and PFS).
#
# Part of: A brain metastasis cell atlas reveals multicellular adaptation to
#          the neural microenvironment (Feigin et al.)
#
# Run from the project root (see README.md). All paths below are relative to it.
#
# Inputs:
#   BMCA/csv/Clinical information for paired samples.csv
#   BMCA/csv/hodge2019_markers_neurons.csv
#   BMCA/csv/oligo_marques.csv
#   BMCA/csv/Sample table with date-BrM-0730.csv
#   BMCA/csv/sigs_all_Astrocyte_2.csv
#   BMCA/data/community_fraction_tib.RDS
#   BMCA/data/ecotypes_tib_proccesed_filtered.RDS
#   BMCA/data/endo_cs.RDS
#   BMCA/data/fib_cs.RDS
#   BMCA/data/metadata_all_studies.rds
#   BMCA/data/metadata_states_proccesed_filtered.RDS
#   BMCA/data/patient_means_MNR_NS.RDS
#   BMCA/data/results_list_not_combined_ecotypes_lenient.RDS
#   BMCA/data/scores_for_MHC.RDS
#   BMCA/data/UMI_pb_<cell_type>.RDS   # path built at run time
#
# Outputs:
#   BMCA/data/  (3 files)
#   BMCA/plots/fig5/  (12 files)
#
# Original filename: Figure_5_4.R
# =============================================================================

# ---- Packages ---------------------------------------------------------------
pkgs <- c(
  "dendextend", "dplyr", "enrichR", "ggdendro", "ggnewscale", "ggplot2",
  "ggpubr", "ggraph", "ggrepel", "Hmisc", "igraph", "lubridate", "Matrix",
  "openxlsx", "patchwork", "purrr", "RColorBrewer", "readxl", "reshape2",
  "scales", "stringr", "survival", "survminer", "tibble", "tidyr",
  "tidyverse", "viridis"
)
invisible(lapply(pkgs, library, character.only = TRUE))

# ---- Project paths ----------------------------------------------------------
source("R/config.R")

# =============================================================================
metadata <- as_tibble(readRDS("BMCA/data/metadata_all_studies.rds"))
# -----------------------------------------------------------------------------

# =============================================================================
# Co-occurrence Network Analysis of Ecotype/State Features
# =============================================================================
# Two cancer flows: melanoma vs. carcinoma (all other primary sites)
# Single threshold tier:
#   Lenient : global mean + 0*SD  &  per-study mean + 1*SD
#
# One plot per seq_tech × location group, melanoma | carcinoma side by side.
# =============================================================================


# =============================================================================
# PARAMETERS  —  adjust these without touching the rest of the code
# =============================================================================

# Input / output files
PATH_METADATA <- "BMCA/data/metadata_states_proccesed_filtered.RDS"
PATH_ECOTYPES <- "BMCA/data/ecotypes_tib_proccesed_filtered.RDS"
PATH_MERGED   <- "BMCA/data/ecotypes_tib_merged_states.RDS"

# Lenient threshold  (GLOBAL_SD_MULT = 0, STUDY_SD_MULT = 1)
LENIENT_GLOBAL_SD_MULT <- 0
LENIENT_STUDY_SD_MULT  <- 1

# Permutation test
N_SHUFFLES  <- 200
PERM_SEED   <- 1
MIN_SAMPLES       <- 10   # minimum total patients in the comparison
MIN_FEATURE_COUNT <- 3 # minimum patients "high" for each individual feature
# Edge significance cutoff
P_ADJ_THRESH <- 0.05
P_ADJ_METHOD <- "BH"
# P-value method: "permutation" (p_adj_shuff) or "fisher" (p_adj)
P_VALUE_METHOD <- "permutation"   # <── switch here

# Edge visual parameters
EDGE_WIDTH  <- 1.2
EDGE_ALPHA  <- 0.7

# Metadata columns (everything else treated as a numeric feature)
META_COLS <- c("patient", "study", "location", "seq_tech",
               "subtype", "primary_site", "Donor")

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================
# ── Select p-value column based on global method switch ──────────────────────
get_padj <- function(pairs_df) {
  if (P_VALUE_METHOD == "permutation") {
    pairs_df$p_adj_shuff
  } else {
    pairs_df$p_adj
  }
}

# ── Cell-type label from feature name prefix ──────────────────────────────────
get_cell_type <- function(feature_name) {
  case_when(
    str_detect(feature_name, "^Malignant")       ~ "Malignant",
    str_detect(feature_name, "^CD8")             ~ "CD8",
    str_detect(feature_name, "^CD4")             ~ "CD4",
    str_detect(feature_name, "^GD")              ~ "GD",
    str_detect(feature_name, "^NK")              ~ "NK",
    str_detect(feature_name, "^Myeloid")         ~ "Myeloid",
    str_detect(feature_name, "^B_Plasma")        ~ "B/Plasma",
    str_detect(feature_name, "^Fibroblast")      ~ "Fibroblast",
    str_detect(feature_name, "^Pericyte")        ~ "Pericyte",
    str_detect(feature_name, "^Endothelial")     ~ "Endothelial",
    str_detect(feature_name, "^Astrocyte")       ~ "Astrocyte",
    str_detect(feature_name, "^Oligodendrocyte") ~ "Oligodendrocyte",
    str_detect(feature_name, "^Neuron")          ~ "Neuron",
    TRUE                                         ~ "Other"
  )
}

# ── Compute per-study thresholds ──────────────────────────────────────────────
compute_study_thresholds <- function(dat_long, sd_mult) {
  dat_long %>%
    group_by(study, feature) %>%
    summarise(
      study_mean = mean(value, na.rm = TRUE),
      study_sd   = sd(value,   na.rm = TRUE),
      .groups    = "drop"
    ) %>%
    mutate(study_cut = study_mean + sd_mult * study_sd)
}

# ── Compute global thresholds ─────────────────────────────────────────────────
compute_global_thresholds <- function(dat_long, sd_mult) {
  dat_long %>%
    group_by(feature) %>%
    summarise(
      global_mean = mean(value, na.rm = TRUE),
      global_sd   = sd(value,   na.rm = TRUE),
      .groups     = "drop"
    ) %>%
    mutate(global_cut = global_mean + sd_mult * global_sd)
}

# ── Flag "high" patients ──────────────────────────────────────────────────────
flag_high <- function(merged_subset, feature_cols,
                      study_sd_mult, global_sd_mult) {

  dat_long <- merged_subset %>%
    dplyr::select(patient, study, seq_tech, location, all_of(feature_cols)) %>%
    pivot_longer(all_of(feature_cols), names_to = "feature", values_to = "value")

  study_thr  <- compute_study_thresholds(dat_long, study_sd_mult)
  global_thr <- compute_global_thresholds(dat_long, global_sd_mult)

  dat_long %>%
    left_join(study_thr  %>% dplyr::select(study, feature, study_cut),  by = c("study", "feature")) %>%
    left_join(global_thr %>% dplyr::select(feature, global_cut),        by = "feature") %>%
    mutate(
      high = case_when(
        is.na(value)                            ~ NA,
        value > study_cut & value > global_cut ~ TRUE,
        TRUE                                    ~ FALSE
      )
    )
}
run_cooccurrence <- function(dat_high_sub) {

  presence_wide <- dat_high_sub %>%
    dplyr::select(patient, feature, high) %>%
    group_by(patient, feature) %>%
    summarise(high = any(high == TRUE), .groups = "drop") %>%
    pivot_wider(names_from = feature, values_from = high, values_fill = FALSE)

  features <- setdiff(colnames(presence_wide), "patient")
  nF       <- length(features)

  if (nF < 2) return(tibble())

  # ── Remove features that are "high" in fewer than MIN_FEATURE_COUNT patients ──
  feature_counts <- colSums(presence_wide[, features], na.rm = TRUE)
  features <- features[feature_counts >= MIN_FEATURE_COUNT]
  nF <- length(features)

  if (nF < 2) return(tibble())

  # Fisher exact tests
  pair_results <- vector("list", nF * (nF - 1) / 2)
  k <- 1

  for (i in seq_len(nF - 1)) {
    for (j in seq(i + 1, nF)) {

      f1  <- features[i]
      f2  <- features[j]
      sub <- presence_wide %>%
        dplyr::select(all_of(c("patient", f1, f2))) %>%
        filter(!is.na(.data[[f1]]) & !is.na(.data[[f2]]))

      if (nrow(sub) < MIN_SAMPLES) next   # ← minimum total patients

      both  <- sum( sub[[f1]] &  sub[[f2]])
      only1 <- sum( sub[[f1]] & !sub[[f2]])
      only2 <- sum(!sub[[f1]] &  sub[[f2]])
      none  <- sum(!sub[[f1]] & !sub[[f2]])

      ft <- fisher.test(
        matrix(c(both, only1, only2, none), nrow = 2),
        alternative = "greater"
      )

      pair_results[[k]] <- tibble(
        feature1   = f1,
        feature2   = f2,
        N          = nrow(sub),
        both_count = both,
        odds_ratio = unname(ft$estimate),
        pvalue     = ft$p.value
      )
      k <- k + 1
    }
  }

  pairdf <- bind_rows(pair_results) %>%
    mutate(p_adj = p.adjust(pvalue, method = P_ADJ_METHOD))

  if (nrow(pairdf) == 0) return(tibble())

  # Permutation test
  mat0 <- presence_wide %>%
    dplyr::select(all_of(features)) %>%
    mutate(across(everything(), ~ as.integer(!is.na(.) & (. == TRUE)))) %>%
    as.matrix()

  feat1_idx <- match(pairdf$feature1, colnames(mat0))
  feat2_idx <- match(pairdf$feature2, colnames(mat0))
  both_obs  <- mapply(function(i, j) sum(mat0[, i] & mat0[, j]),
                      feat1_idx, feat2_idx)

  shuffle_once <- function() {
    m_sh <- apply(mat0, 2, sample)
    mapply(function(i, j) sum(m_sh[, i] & m_sh[, j]), feat1_idx, feat2_idx)
  }

  set.seed(PERM_SEED)
  shuff_counts <- replicate(N_SHUFFLES, shuffle_once())

  pairdf$pvalue_shuff <- vapply(seq_along(both_obs), function(k) {
    if (is.na(both_obs[k])) return(NA_real_)
    mean(shuff_counts[k, ] >= both_obs[k])
  }, numeric(1))

  pairdf$p_adj_shuff <- p.adjust(pairdf$pvalue_shuff, method = P_ADJ_METHOD)

  pairdf
}

# ── Build one network panel ───────────────────────────────────────────────────
build_network_plot <- function(pairs, panel_title) {

  sig <- if (nrow(pairs) > 0) {
    pairs %>% filter(get_padj(pairs) < P_ADJ_THRESH)
  } else tibble()

  if (nrow(sig) == 0) {
    return(ggplot() +
             annotate("text", x = 0.5, y = 0.5,
                      label = paste0("No significant edges\n", panel_title),
                      hjust = 0.5, vjust = 0.5, size = 4) +
             theme_void() + labs(title = panel_title))
  }

  edge_df <- sig %>%
    mutate(n1 = pmin(feature1, feature2),
           n2 = pmax(feature1, feature2)) %>%
    distinct(n1, n2, .keep_all = TRUE)

  g <- graph_from_data_frame(edge_df %>% dplyr::select(n1, n2), directed = FALSE)
  V(g)$cell_type <- get_cell_type(V(g)$name)
  V(g)$degree    <- igraph::degree(g)

  set.seed(PERM_SEED)
  layout_fr <- create_layout(g, layout = "fr")

  cell_types <- sort(unique(V(g)$cell_type))
  pal <- colorRampPalette(brewer.pal(min(8, length(cell_types)), "Set2"))(length(cell_types))
  names(pal) <- cell_types

  edge_list <- igraph::as_data_frame(g, what = "edges") %>%
    mutate(
      from_idx = match(from, V(g)$name),
      to_idx   = match(to,   V(g)$name),
      x        = layout_fr$x[from_idx],
      y        = layout_fr$y[from_idx],
      xend     = layout_fr$x[to_idx],
      yend     = layout_fr$y[to_idx]
    )

  ggraph(layout_fr) +
    geom_segment(
      data = edge_list,
      aes(x = x, y = y, xend = xend, yend = yend),
      colour    = "steelblue",
      linewidth = EDGE_WIDTH,
      alpha     = EDGE_ALPHA
    ) +
    geom_node_point(aes(color = cell_type, size = degree)) +
    geom_node_text(aes(label = name), repel = TRUE, size = 2.5) +
    scale_color_manual(values = pal, name = "Cell Type") +
    scale_size_continuous(range = c(3, 10), name = "Degree") +
    theme_void() +
    labs(
      title    = panel_title,
      subtitle = paste0("Significant edges: ", nrow(edge_df),
                        "  |  Nodes: ", vcount(g))
    ) +
    theme(
      plot.title      = element_text(size = 12, face = "bold"),
      legend.position = "bottom",
      legend.box      = "vertical"
    )
}

# =============================================================================
# 1. LOAD & MERGE DATA
# =============================================================================

metadata        <- readRDS(PATH_METADATA)
ecotypes_states <- readRDS(PATH_ECOTYPES)

ecotypes_states$CD8_CT_total <- NULL
clusters_fibro <- readRDS("BMCA/data/fib_cs.RDS") %>%
  pivot_wider(names_from = cluster, values_from = score)
clusters_endo <- readRDS("BMCA/data/endo_cs.RDS") %>%
  pivot_wider(names_from = cluster, values_from = score)
clusters_all <- full_join(clusters_fibro, clusters_endo, by = "patient")
clusters_all <- clusters_all %>%
  rename_with(~ .x %>%
                str_replace("^fib_",  "Fibroblast_") %>%
                str_replace("^endo_", "Endothelial_"),
              -patient)

ecotypes_states <- ecotypes_states %>%
  left_join(clusters_all, by = "patient")

merged <- ecotypes_states %>%
  left_join(metadata, by = "patient")

saveRDS(merged, PATH_MERGED)

feature_cols <- setdiff(names(merged), META_COLS)
message("Loaded: ", nrow(merged), " patients | ", length(feature_cols), " features.")

# =============================================================================
# 2. SPLIT INTO CANCER FLOWS
# =============================================================================

merged_melanoma  <- merged %>% filter(primary_site == "melanoma")
merged_carcinoma <- merged %>% filter(primary_site != "melanoma")

message("Melanoma  patients : ", nrow(merged_melanoma))
message("Carcinoma patients : ", nrow(merged_carcinoma))

# =============================================================================
# 3. FLAG "HIGH" PATIENTS — lenient tier only, both flows
# =============================================================================

message("\nComputing high-patient flags (lenient tier only)...")

dat_high <- list(
  melanoma  = flag_high(merged_melanoma,  feature_cols,
                        LENIENT_STUDY_SD_MULT, LENIENT_GLOBAL_SD_MULT),
  carcinoma = flag_high(merged_carcinoma, feature_cols,
                        LENIENT_STUDY_SD_MULT, LENIENT_GLOBAL_SD_MULT)
)

# =============================================================================
# 4. LOOP OVER seq_tech × location — one combined side-by-side plot per group
# =============================================================================

groups       <- merged %>% distinct(seq_tech, location) %>% drop_na()
results_list <- list()

for (g in seq_len(nrow(groups))) {

  seqg       <- groups$seq_tech[g]
  locg       <- groups$location[g]
  group_name <- paste0(seqg, "_", locg)
  message("\n══ Group: ", group_name, " ══")

  panel_plots <- list()

  for (flow in c("melanoma", "carcinoma")) {

    message("  ── Flow: ", flow)

    sub <- dat_high[[flow]] %>%
      filter(seq_tech == seqg, location == locg)

    n_patients <- n_distinct(sub$patient)

    if (n_patients == 0) {
      message("    No patients — empty panel.")
      panel_plots[[flow]] <- ggplot() +
        annotate("text", x = 0.5, y = 0.5,
                 label = paste0("No data\n(", flow, ")"),
                 hjust = 0.5, vjust = 0.5, size = 4) +
        theme_void() +
        labs(title = paste0(str_to_title(flow), " — no data"))
      next
    }

    pairs <- run_cooccurrence(sub)

    panel_plots[[flow]] <- build_network_plot(
      pairs       = pairs,
      panel_title = paste0(str_to_title(flow), "  (n = ", n_patients, ")")
    )

    results_list[[group_name]][[flow]] <- pairs

    message("    Pairs tested: ", nrow(pairs),
            "  | significant: ",
            if (nrow(pairs) > 0) sum(get_padj(pairs) < P_ADJ_THRESH, na.rm = TRUE) else 0)
  }

  combined <- (panel_plots[["melanoma"]] | panel_plots[["carcinoma"]]) +
    plot_annotation(
      title    = paste0("Co-occurrence Network  |  ", group_name),
      subtitle = "Lenient threshold: global mean  &  per-study mean + 1 SD",
      theme = theme(
        plot.title    = element_text(size = 15, face = "bold"),
        plot.subtitle = element_text(size = 10, color = "grey40")
      )
    )

  print(combined)
  message("  \u2713 Plot printed for: ", group_name)
}
message("\n\u2713 All groups processed. Results stored in results_list.")
saveRDS(results_list, "BMCA/data/results_list_not_combined_ecotypes_lenient.RDS")

# =============================================================================
# 5. COMBINE sn + sc MODALITIES — one tibble per location × cancer_type
# =============================================================================

results_list <- readRDS("BMCA/data/results_list_not_combined_ecotypes_lenient.RDS")

combine_modalities <- function(results_list, location, cancer_type) {

  get_pairs <- function(modality) {
    key <- paste0(modality, "_", location)
    if (!key %in% names(results_list))          return(NULL)
    if (!cancer_type %in% names(results_list[[key]])) return(NULL)
    results_list[[key]][[cancer_type]] %>%
      mutate(n1 = pmin(feature1, feature2),
             n2 = pmax(feature1, feature2),
             modality_src = modality)
  }

  all_pairs <- bind_rows(get_pairs("sn"), get_pairs("sc"))

  if (is.null(all_pairs) || nrow(all_pairs) == 0) return(tibble())

  sig_pairs <- all_pairs %>%
    group_by(n1, n2) %>%
    filter(any(get_padj(pick(everything())) < P_ADJ_THRESH)) %>%
    ungroup()

  if (nrow(sig_pairs) == 0) return(tibble())

  sig_pairs %>%
    group_by(n1, n2) %>%
    summarise(
      modality = {
        sn_sig <- any(modality_src == "sn" & get_padj(pick(everything())) < P_ADJ_THRESH)
        sc_sig <- any(modality_src == "sc" & get_padj(pick(everything())) < P_ADJ_THRESH)
        if      (sn_sig & sc_sig) "both"
        else if (sn_sig)          "sn"
        else                      "sc"
      },
      mean_odds_ratio = mean(odds_ratio,   na.rm = TRUE),
      mean_pval_shuff = mean(pvalue_shuff, na.rm = TRUE),
      min_padj_shuff  = min(p_adj_shuff,   na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(location = location, cancer_type = cancer_type)
}

locations    <- c("Intracranial", "Extracranial")
cancer_types <- c("melanoma", "carcinoma")

combined_results_list <- setNames(
  lapply(locations, function(loc) {
    setNames(
      lapply(cancer_types, function(ct) combine_modalities(results_list, loc, ct)),
      cancer_types
    )
  }),
  locations
)

for (loc in locations) {
  for (ct in cancer_types) {
    tib <- combined_results_list[[loc]][[ct]]
    message(sprintf("%-15s | %-10s → %d edges  |  modality: [%s]",
                    loc, ct, nrow(tib),
                    paste(sort(unique(tib$modality)), collapse = ", ")))
  }
}

# =============================================================================
# 6. PLOT — combined_results_list
# Edge colour = modality (sn / sc / both); node fill = cell type
# =============================================================================


# ── Cell-type levels (order = legend order) ──────────────────────────────────
major_cts <- c("Neuron", "Astrocyte", "Oligodendrocyte",
               "Fibroblast", "Pericyte", "Endothelial",
               "CD4", "CD8", "NK", "GD", "T_NK",
               "B_Plasma", "Myeloid", "Malignant", "Other")

# ── Manual palette (NPG-derived where a hue fits) ────────────────────────────
CT_COLORS <- c(
  Neuron          = "#3E8E41",   # green      ┐
  Astrocyte       = "#EFC000",   # yellow     ├ brain lineages
  Oligodendrocyte = "#00A087",   # teal       ┘
  Fibroblast      = "#7E6148",   # brown
  Pericyte        = "#B09C85",   # tan
  Endothelial     = "#8B1A1A",   # dark red
  NK             = "#9ECAE1",   # T/NK lightest ┐
  CD8             = "#4292C6",   #               │
  CD4              = "#2171B5",   #               ├ blue ramp
  GD              = "#08519C",   #               │
  T_NK            = "#08306B",   # T/NK darkest  ┘
  B_Plasma        = "#7D4E9E",   # purple
  Myeloid         = "#E8871A",   # orange
  Malignant       = "#E64B35",   # red
  Other           = "#BEBEBE"    # unmapped — should stay empty
)

MODALITY_COLORS <- c(sn = "#a8d1f5", sc = "#778873", both = "#6a0dad")

stopifnot(
  "palette must cover every level" = all(major_cts %in% names(CT_COLORS)),
  "palette has keys outside major_cts" = all(names(CT_COLORS) %in% major_cts)
)

if (!exists("PERM_SEED")) PERM_SEED <- 42L

# ── Cell-type resolution ─────────────────────────────────────────────────────
# Longest-prefix match against major_cts first (fixes T_NK and any T_NK_*),
# then fall back to get_cell_type(), then normalise "B/Plasma" -> "B_Plasma".
resolve_ct <- function(x, levels = major_cts) {

  x  <- as.character(x)
  lv <- setdiff(levels, "Other")
  lv <- lv[order(nchar(lv), decreasing = TRUE)]

  hit <- vapply(x, function(nm) {
    m <- lv[nm == lv | startsWith(nm, paste0(lv, "_"))]
    if (length(m)) m[1] else NA_character_
  }, character(1), USE.NAMES = FALSE)

  if (anyNA(hit)) {
    fb <- tryCatch(as.character(get_cell_type(x)),
                   error = function(e) rep(NA_character_, length(x)))
    fb <- str_replace_all(fb, "/", "_")
    fb[!fb %in% names(CT_COLORS)] <- NA_character_
    hit <- ifelse(is.na(hit), fb, hit)
  }

  hit[is.na(hit)] <- "Other"
  factor(hit, levels = levels)
}

# ── Display labels: strip cell-type prefix, no "*" ───────────────────────────
make_label <- function(node_name, cell_type) {

  nm <- str_replace_all(as.character(node_name), "/", "_")
  ct <- str_replace_all(as.character(cell_type), "/", "_")
  ct[is.na(ct)] <- ""

  lab        <- nm
  has_prefix <- nzchar(ct) & startsWith(nm, paste0(ct, "_"))
  is_bare    <- nzchar(ct) & nm == ct

  lab[has_prefix] <- substring(nm[has_prefix], nchar(ct[has_prefix]) + 2L)

  lab <- str_replace_all(lab, fixed("MHC.I"),    "MHCI")
  lab <- str_replace_all(lab, fixed("targets."), "targets/")
  lab <- str_replace_all(lab, "_", " ")

  lab[is_bare]    <- str_replace_all(ct[is_bare], "_", "-")   # T_NK -> T-NK
  lab[!nzchar(lab)] <- nm[!nzchar(lab)]                       # never blank
  lab
}

empty_panel <- function(msg) {
  ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = msg,
             hjust = 0.5, vjust = 0.5, size = 3) +
    theme_void()
}

# ── Panel builder ────────────────────────────────────────────────────────────
build_combined_network_plot <- function(edge_tib, panel_title = "",
                                        ct_levels = major_cts) {

  if (is.null(edge_tib) || !is.data.frame(edge_tib) || nrow(edge_tib) == 0)
    return(empty_panel(paste0("No edges\n", panel_title)))

  if (!all(c("n1", "n2", "modality") %in% names(edge_tib)))
    return(empty_panel(paste0("Missing columns\n", panel_title)))

  edge_tib_vis <- edge_tib %>%
    filter(!is.na(n1), !is.na(n2), n1 != n2) %>%
    mutate(
      modality        = as.character(modality),
      edge_colour_val = unname(MODALITY_COLORS[modality]),
      edge_colour_val = ifelse(is.na(edge_colour_val), "grey50", edge_colour_val)
    )

  if (nrow(edge_tib_vis) == 0)
    return(empty_panel(paste0("No edges\n", panel_title)))

  g <- graph_from_data_frame(
    d        = edge_tib_vis %>% dplyr::select(n1, n2, edge_colour_val, modality),
    directed = FALSE
  )

  V(g)$cell_type <- resolve_ct(V(g)$name)          # raw names in, always
  V(g)$degree    <- igraph::degree(g)
  V(g)$label_txt <- make_label(V(g)$name, as.character(V(g)$cell_type))

  unmapped <- unique(V(g)$name[V(g)$cell_type == "Other"])
  if (length(unmapped))
    warning(panel_title, " — unmapped nodes: ",
            paste(unmapped, collapse = ", "), call. = FALSE)

  set.seed(PERM_SEED)
  layout_fr <- create_layout(g, layout = "fr")

  edge_list <- igraph::as_data_frame(g, what = "edges") %>%
    mutate(
      from_idx = match(from, V(g)$name),
      to_idx   = match(to,   V(g)$name),
      x        = layout_fr$x[from_idx],
      y        = layout_fr$y[from_idx],
      xend     = layout_fr$x[to_idx],
      yend     = layout_fr$y[to_idx]
    )

  modality_legend_df <- tibble(
    x = rep(NA_real_, 3), y = rep(NA_real_, 3),
    modality = factor(names(MODALITY_COLORS), levels = names(MODALITY_COLORS))
  )

  ggraph(layout_fr) +
    geom_segment(
      data = edge_list,
      aes(x = x, y = y, xend = xend, yend = yend, colour = edge_colour_val),
      linewidth = 1.2, alpha = 0.75
    ) +
    scale_colour_identity(guide = "none") +
    geom_node_point(aes(fill = cell_type),
                    shape = 21, size = 5, colour = "white",
                    show.legend = TRUE) +
    geom_node_text(aes(label = label_txt), repel = TRUE, size = 2,
                   max.overlaps = Inf, seed = PERM_SEED) +
    scale_fill_manual(
      values   = CT_COLORS,
      limits   = ct_levels,
      drop     = FALSE,
      na.value = "grey80",
      labels   = function(x) str_replace_all(x, "_", "-"),
      name     = "Cell type",
      guide    = guide_legend(
        override.aes   = list(shape = 21, size = 3.5, colour = "white"),
        title.position = "top", nrow = 4, order = 1
      )
    ) +
    ggnewscale::new_scale_colour() +
    geom_point(data = modality_legend_df,
               aes(x = x, y = y, colour = modality), na.rm = TRUE) +
    scale_colour_manual(
      values = MODALITY_COLORS,
      breaks = names(MODALITY_COLORS),
      labels = c(sn = "Single-nucleus only", sc = "Single-cell only",
                 both = "Both"),
      name   = "Modality",
      guide  = guide_legend(override.aes = list(size = 3),
                            title.position = "top", order = 2)
    ) +
    theme_void() +
    labs(title = NULL) +
    theme(
      legend.position = "bottom",
      legend.box      = "horizontal",
      legend.title    = element_text(size = 9, face = "bold"),
      legend.text     = element_text(size = 8)
    )
}

# ── Driver ───────────────────────────────────────────────────────────────────
network_figs <- list()

for (loc in c("Extracranial", "Intracranial")) {

  panel_plots <- list()

  for (ct in c("melanoma", "carcinoma")) {

    tib <- tryCatch(combined_results_list[[loc]][[ct]],
                    error = function(e) NULL)

    panel_plots[[ct]] <- build_combined_network_plot(
      edge_tib    = tib,
      panel_title = paste0(str_to_title(ct), " | ", loc)
    )
  }

  network_figs[[loc]] <- (panel_plots[["melanoma"]] | panel_plots[["carcinoma"]]) +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom")

  print(network_figs[[loc]])
  message("\u2713 Plotted: ", loc)
}

panel_plots[["carcinoma"]]

ggsave("BMCA/plots/fig5/carcinoma_graph.pdf",
       plot = panel_plots[["carcinoma"]], width = 200, height = 150,
       units = "mm", dpi = 300)

# -----------------------------------------------------------------------------
naa

metadata <- as_tibble(readRDS("BMCA/data/metadata_all_studies.rds"))

metadata_p <- metadata %>%
  dplyr::select(patient, primary_site, study,
                immunotherapy, Location, Sex, Age,
                Os_status,
                Os_month, PFS) %>%
  distinct() %>%
  filter(study == "lin_sn") %>%
  filter(primary_site != "melanoma") %>%
  distinct()

ecotypes_states_sub <- ecotypes_states %>%
  filter(patient %in% metadata_p$patient)
community_fraction_tib <- ecotypes_states_sub %>%
  transmute(
    patient,
    Brain_adap = rowMeans(
      across(c(Malignant_Neuronal_related_2,
               Malignant_Neuronal_related_3,
               Myeloid_Scavenger_Immunosuppressive,
               Fibroblast_C1,
               Myeloid_Stress,
               Malignant_Chromatin_Remodeling),
             ~ as.numeric(scale(.x))),   # (x - mean) / sd, per column
      na.rm = TRUE
    )
  ) %>%
  mutate(Brain_adap = Brain_adap - mean(Brain_adap, na.rm = TRUE))
saveRDS(community_fraction_tib,"BMCA/data/community_fraction_tib.RDS")

community_fraction_tib <- readRDS("BMCA/data/community_fraction_tib.RDS")
bp_df <- community_fraction_tib %>%
  inner_join(
    metadata_p %>% dplyr::select(patient, primary_site) %>% distinct(),
    by = "patient"
  ) %>%
  filter(!is.na(Brain_adap), !is.na(primary_site)) %>%
  mutate(primary_site = fct_reorder(primary_site, Brain_adap, .fun = median))

# every pairwise combination of the primary_site levels
pairs <- combn(levels(bp_df$primary_site), 2, simplify = FALSE)

high_co <- 0.5
low_co <- 0.15
bp_df <- bp_df %>%
  mutate(primary_site = fct_relevel(primary_site,
                                    "breast", "NSCLC", "colorectal"))
p <- ggplot(bp_df, aes(x = primary_site, y = Brain_adap)) +
  geom_boxplot(aes(fill = primary_site),
               outlier.shape = NA, width = 0.6, alpha = 0.33,
               colour = "grey25", linewidth = 0.4) +
  geom_jitter(width = 0.2, size = 1, alpha = 0.45, colour = "black") +
  geom_hline(yintercept = high_co, linetype = "dashed",
             colour = "black", linewidth = 0.6) +
  geom_hline(yintercept = low_co, linetype = "dashed",
             colour = "black", linewidth = 0.6) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
  labs(x = NULL, y = "BrM adaptation score", fill = "Primary site") +
  theme_classic(base_size = 7) +
  theme(
    axis.text.x     = element_text(angle = 30, hjust = 1, colour = "black"),
    axis.text.y     = element_text(colour = "black"),
    legend.position = "none"
  )

p

ggsave("BMCA/plots/fig5/brm_enrich_boxes.pdf", plot = p, width = 30, height = 50,
       units = "mm", dpi = 300)


grp <- ifelse(bp_df$primary_site == "colorectal", "colorectal", "other")

wilcox.test(Brain_adap ~ grp, data = bp_df)   # Wilcoxon
t.test(Brain_adap ~ grp, data = bp_df)        # t-test

# -----------------------------------------------------------------------------

Sample_clic_data <- read_csv("BMCA/csv/Sample table with date-BrM-0730.csv")
Sample_clic_data <- Sample_clic_data %>%
  filter(Cohort == "PUMC") %>%
  dplyr::select(Sample_ID, Anatomical_location, Sex, Age_at_sampling,
                'Treatment_received (time)', 'Sampling date', 'The date of progression',
                'The date of death', 'Last follow-up', 'OS status')
# See the actual non-NA values and their format

# Also useful:

# Parse only valid dd/mm/yyyy strings; everything else -> NA (no warnings)
parse_dmy <- function(x) {
  x <- trimws(as.character(x))
  x[!grepl("^\\d{1,2}/\\d{1,2}/\\d{4}$", x)] <- NA
  dmy(x)
}

Sample_clic_data %>%
  summarise(across(c(`Sampling date`, `The date of progression`,
                     `The date of death`, `Last follow-up`, `OS status`),
                   list(class   = ~ class(.x)[1],
                        n_nonNA = ~ sum(!is.na(.x))),
                   .names = "{.col}__{.fn}")) %>%
  tidyr::pivot_longer(everything()) %>%
  print(n = Inf)

Sample_clic_data %>%
  dplyr::select(`The date of progression`, `The date of death`, `Last follow-up`) %>%
  tidyr::pivot_longer(everything()) %>%
  filter(!is.na(value)) %>%
  distinct(name, value) %>%
  print(n = 40)

Sample_clic_data %>%
  summarise(across(c(`Sampling date`, `The date of progression`,
                     `The date of death`, `Last follow-up`, `OS status`),
                   list(class   = ~ class(.x)[1],
                        n_nonNA = ~ sum(!is.na(.x))),
                   .names = "{.col}__{.fn}")) %>%
  tidyr::pivot_longer(everything(),
                      values_transform = list(value = as.character)) %>%
  print(n = Inf)

# =============================================================================
# CLINICAL DATE PARSING + SURVIVAL ENDPOINT DERIVATION
# =============================================================================

DATE_SENTINELS <- c("", "/", "//", "-", "--", ".", "NA", "na", "n/a", "N/A",
                    "no relapse", "No relapse", "unknown", "Unknown",
                    "alive", "Alive", "not applicable", "NR", "nr")

# --- inspect first, then fill these in ---------------------------------------
print(count(Sample_clic_data, `OS status`))

DEAD_LABELS  <- c("Dead", "dead", "DOD", "Deceased", "1")
ALIVE_LABELS <- c("Alive", "alive", "NED", "AWD", "0")

# ── Date parser: dmy only, sentinels to NA, loud on anything else ────────────
parse_date_any <- function(x, colname = "") {

  if (inherits(x, "Date"))    return(x)
  if (inherits(x, "POSIXct")) return(as.Date(x))
  if (is.numeric(x))          return(as.Date(x, origin = "1899-12-30"))

  x <- str_squish(as.character(x))
  x[x %in% DATE_SENTINELS] <- NA_character_

  out <- suppressWarnings(
    as.Date(parse_date_time(x, orders = "dmy", quiet = TRUE))
  )

  bad <- !is.na(x) & is.na(out)
  if (any(bad))
    warning(colname, " — unparsed: ",
            paste(unique(x[bad]), collapse = " | "), call. = FALSE)

  out
}

norm_status <- function(x) {
  s <- str_squish(as.character(x))
  s[s %in% DATE_SENTINELS] <- NA_character_
  case_when(
    s %in% DEAD_LABELS  ~ "dead",
    s %in% ALIVE_LABELS ~ "alive",
    is.na(s)            ~ NA_character_,
    TRUE                ~ "unmapped"
  )
}

# ── Derive endpoints ─────────────────────────────────────────────────────────
Sample_clic_data <- Sample_clic_data %>%
  mutate(
    date_sampling    = parse_date_any(`Sampling date`,             "Sampling date"),
    date_progression = parse_date_any(`The date of progression`,   "Progression"),
    date_death       = parse_date_any(`The date of death`,         "Death"),
    date_lastfu      = parse_date_any(`Last follow-up`,            "Last follow-up"),

    os_status_norm = norm_status(`OS status`),

    # a death is known to have occurred if either the date or the status says so
    death_known = !is.na(date_death) | os_status_norm %in% "dead",
    # ...but a death with no date has no usable time
    death_dateless = death_known & is.na(date_death),

    # ---- OS ----
    OS_event   = as.integer(!is.na(date_death)),
    OS_enddate = coalesce(date_death, date_lastfu),
    OS_month   = as.numeric(OS_enddate - date_sampling) / 30.44,

    # ---- PFS ----
    PFS_event   = as.integer(!is.na(date_progression) | !is.na(date_death)),
    PFS_enddate = pmin(date_progression, date_death, na.rm = TRUE),
    PFS_enddate = coalesce(PFS_enddate, date_lastfu),
    PFS_month   = as.numeric(PFS_enddate - date_sampling) / 30.44,

    # ---- validity flags ----
    flag_no_sampling = is.na(date_sampling),
    flag_no_anchor   = is.na(OS_enddate),
    flag_negative    = (!is.na(OS_month)  & OS_month  < 0) |
      (!is.na(PFS_month) & PFS_month < 0),
    flag_pfs_gt_os   = !is.na(PFS_month) & !is.na(OS_month) & PFS_month > OS_month + 1e-6,

    OS_evaluable  = !flag_no_sampling & !is.na(OS_month)  & OS_month  >= 0 & !death_dateless,
    PFS_evaluable = !flag_no_sampling & !is.na(PFS_month) & PFS_month >= 0 & !death_dateless
  )

# ── QC report ────────────────────────────────────────────────────────────────
qc <- Sample_clic_data %>%
  summarise(
    n_rows            = n(),
    n_sampling_parsed = sum(!is.na(date_sampling)),
    n_no_anchor       = sum(flag_no_anchor),
    n_death_dateless  = sum(death_dateless, na.rm = TRUE),
    n_status_unmapped = sum(os_status_norm %in% "unmapped"),
    n_negative        = sum(flag_negative, na.rm = TRUE),
    n_pfs_gt_os       = sum(flag_pfs_gt_os, na.rm = TRUE),
    n_OS_evaluable    = sum(OS_evaluable),
    n_OS_events       = sum(OS_event[OS_evaluable]),
    n_PFS_evaluable   = sum(PFS_evaluable),
    n_PFS_events      = sum(PFS_event[PFS_evaluable]),
    median_OS_month   = median(OS_month[OS_evaluable],  na.rm = TRUE),
    median_PFS_month  = median(PFS_month[PFS_evaluable], na.rm = TRUE)
  )
print(as.data.frame(qc))

# status vs date agreement — the check that matters most
Sample_clic_data %>%
  count(os_status_norm, has_death_date = !is.na(date_death)) %>%
  print()

# every row that will be dropped, and why
Sample_clic_data %>%
  filter(!OS_evaluable) %>%
  transmute(Sample_ID, `OS status`, date_sampling, date_death, date_lastfu,
            reason = case_when(
              flag_no_sampling ~ "no sampling date",
              death_dateless   ~ "death known, date missing",
              flag_no_anchor   ~ "no death or follow-up date",
              flag_negative    ~ "negative interval",
              TRUE             ~ "other"
            )) %>%
  print(n = Inf)

surv_df <- filter(Sample_clic_data, OS_evaluable)
Sample_clic_data$PFS <- Sample_clic_data$PFS_month
metadata$samp_date <- NULL
metadata$brain_adaptation_score <- NULL
metadata$Os_month <- NULL
metadata$Os_status <- NULL
metadata$PFS <- NULL
metadata$`Dissamination Model` <- NULL
metadata <- metadata[,1:31]


Sample_clic_data <- Sample_clic_data %>%
  dplyr::select(Sample_ID, 'OS_month', 'PFS', `OS status`, Anatomical_location,
                Sex, Age_at_sampling)

colnames(Sample_clic_data) <- c("patient", "Os_month", "PFS", "Os_status", "Location",
                                "Sex", "Age")

Sample_clic_data$patient <- paste0("lin_sn_", Sample_clic_data$patient)

intersect(Sample_clic_data$patient, metadata_p$patient)

Sample_clic_data

metadata <- metadata %>% left_join(Sample_clic_data, by = "patient")

metadata$Location[metadata$Location == "Intracranial"] <- "Not specified"


metadata_p <- metadata %>%
  dplyr::select(patient, primary_site, study,
                Location, Sex, Age,
                Os_status, Os_month, PFS) %>%
  filter(patient %in% community_fraction_tib$patient) %>%
  filter(study == "lin_sn") %>%
  filter(primary_site %in% c("breast", "NSCLC")) %>%
  distinct() %>%
  filter(!is.na(Sex)) %>%          # fixed: is.na(), not isNA()
  distinct(patient, .keep_all = TRUE)   # force exactly 1 row per patient
community_fraction_tib_sub <- community_fraction_tib %>%
  left_join(metadata_p, by = "patient") %>%
  filter(primary_site %in% c("breast", "NSCLC"))

community_fraction_tib_sub

community_fraction_tib_sub <- community_fraction_tib_sub %>%
  mutate(`BrM Remodeling Score` = case_when(
    Brain_adap > high_co ~ "High",
    Brain_adap < low_co  ~ "Low",
  ))


# ── Setup ─────────────────────────────────────────────────────────────────────
score_levels <- c("Low", "Intermediate", "High")
score_cols   <- c(Low = "#3B6EA5", Intermediate = "grey75", High = "#C0392B")

dat <- community_fraction_tib_sub %>%
  mutate(score = factor(`BrM Remodeling Score`, levels = score_levels)) %>%
  filter(!is.na(score))

# ── Enrichment helper ─────────────────────────────────────────────────────────
# Fisher's exact test for the global group × score association, plus a per-cell
# hypergeometric test asking whether each score class is enriched within a group.
enrich_test <- function(data, group_var) {
  d <- data %>%
    filter(!is.na(.data[[group_var]])) %>%
    mutate(grp = droplevels(as.factor(.data[[group_var]])))

  tab <- table(d$grp, d$score)

  fisher_p <- tryCatch(
    fisher.test(tab, simulate.p.value = TRUE, B = 1e5)$p.value,
    error = function(e) NA_real_
  )

  N <- sum(tab)
  hyper <- expand.grid(grp   = rownames(tab),
                       score = colnames(tab),
                       stringsAsFactors = FALSE) %>%
    mutate(
      observed  = mapply(function(g, s) tab[g, s],   grp, score),  # k in cell
      group_n   = mapply(function(g, s) sum(tab[g, ]), grp, score), # n drawn
      score_tot = mapply(function(g, s) sum(tab[, s]), grp, score), # K successes
      expected  = group_n * score_tot / N,
      # P(>= observed) : enrichment of this score class in this group
      enrich_p  = phyper(observed - 1, score_tot, N - score_tot, group_n,
                         lower.tail = FALSE),
      # P(<= observed) : depletion
      deplete_p = phyper(observed, score_tot, N - score_tot, group_n,
                         lower.tail = TRUE)
    ) %>%
    group_by(score) %>%
    mutate(enrich_padj = p.adjust(enrich_p, method = "BH")) %>%
    ungroup() %>%
    as_tibble()

  list(table = tab, fisher_p = fisher_p, hyper = hyper)
}

# ── Categorical variables: stacked proportions + enrichment ───────────────────
cat_vars <- c("Location", "Sex")

for (var in cat_vars) {

  res <- enrich_test(dat, var)

  cat("\n================ ", var, " ================\n")
  cat("Contingency table (group x score):\n"); print(res$table)
  cat(sprintf("\nFisher's exact test (global association): p = %.4g\n",
              res$fisher_p))
  cat("\nPer-cell hypergeometric enrichment:\n")
  print(res$hyper %>%
          arrange(enrich_p) %>%
          dplyr::select(grp, score, observed, expected, enrich_p, enrich_padj))

  df_plot <- dat %>%
    filter(!is.na(.data[[var]])) %>%
    mutate(group = droplevels(as.factor(.data[[var]])))

  p <- ggplot(df_plot, aes(x = group, fill = score)) +
    geom_bar(position = "fill", alpha = 0.9, colour = "white", linewidth = 0.3) +
    scale_fill_manual(values = score_cols, name = "BrM adaptation score") +
    scale_y_continuous(labels = scales::percent_format(), expand = c(0, 0)) +
    labs(
      title    = sprintf("Fisher's exact p = %.3g", res$fisher_p),
      x        = var,
      y        = "Proportion of patients"
    ) +
    theme_classic(base_size = 11, base_family = "Helvetica") +
    theme(
      axis.title  = element_text(face = "bold", size = 10),
      axis.text   = element_text(size = 9, color = "grey10"),
      axis.text.x = element_text(angle = 30, hjust = 1),
      plot.title  = element_text(hjust = 0.5, size = 10)
    )

  print(p)
}

# ── Same, stratified by primary_site ──────────────────────────────────────────
for (var in cat_vars) {
  for (site in unique(dat$primary_site)) {

    d_site <- dat %>% filter(primary_site == site)
    if (nrow(d_site) < 3 || n_distinct(d_site[[var]], na.rm = TRUE) < 2) next

    res <- enrich_test(d_site, var)
    cat(sprintf("\n---- %s | %s ----\n", site, var))
    print(res$table)
    cat(sprintf("Fisher's exact p = %.4g\n", res$fisher_p))
    print(res$hyper %>% arrange(enrich_p) %>%
            dplyr::select(grp, score, observed, expected, enrich_p))
  }

  df_plot <- dat %>%
    filter(!is.na(.data[[var]])) %>%
    mutate(group = droplevels(as.factor(.data[[var]])))
  p <- ggplot(df_plot, aes(x = group, fill = score)) +
    geom_bar(position = "fill", alpha = 0.9, colour = "white", linewidth = 0.3) +
    facet_wrap(~ primary_site, scales = "free_x") +
    scale_fill_manual(values = score_cols, name = "BrM adaptation score") +
    scale_y_continuous(labels = scales::percent_format(), expand = c(0, 0)) +
    labs(x = var, y = "Proportion of patients") +
    theme_classic(base_size = 7, base_family = "Helvetica") +
    theme(
      axis.text.x      = element_text(angle = 30, hjust = 1),
      strip.background = element_rect(fill = "grey95", colour = NA),
      legend.position  = "right"
    )
  print(p)
}

ggsave("BMCA/plots/fig5/Sex_remo.pdf", plot = p, width = 100, height = 75,
       units = "mm", dpi = 300)

# ── Location, excluding "Not specified" ───────────────────────────────────────
dat_loc <- dat %>% filter(Location != "Not specified")
res_loc  <- enrich_test(dat_loc, "Location")
cat("\n================  Location (excl. 'Not specified')  ================\n")
print(res_loc$table)
cat(sprintf("Fisher's exact p = %.4g\n", res_loc$fisher_p))
print(res_loc$hyper %>% arrange(enrich_p) %>%
        dplyr::select(grp, score, observed, expected, enrich_p, enrich_padj))

# ── Age: continuous predictor, categorical outcome ────────────────────────────
# A hypergeometric test needs a categorical predictor, so it doesn't apply to
# Age directly. The categorical analog of the old Age~score correlation is to
# ask whether Age distribution differs across score classes (Kruskal-Wallis).
df_age <- dat %>% filter(!is.na(Age))

p_age <- ggplot(df_age, aes(x = score, y = Age)) +
  geom_boxplot(outlier.size = 0.8, linewidth = 0.4, alpha = 0.85) +
  geom_jitter(width = 0.15, size = 1, alpha = 0.5, colour = "grey20") +
  stat_compare_means(method = "kruskal.test", size = 3.5) +
  facet_wrap(~ primary_site) +
  labs(x = "BrM adaptation score", y = "Age") +
  theme_classic(base_size = 7, base_family = "Helvetica") +
  theme(
    legend.position = "none",
    strip.background = element_rect(fill = "grey95", colour = NA)
  )
print(p_age)
ggsave("BMCA/plots/fig5/age_remo.pdf", plot = p_age, width = 70, height = 75,
       units = "mm", dpi = 300)


# ── Categorical variables: boxplot ────────────────────────────────────────────
cat_vars <- c("Location", "Sex")

for (var in cat_vars) {

  df_plot <- community_fraction_tib_sub %>%
    filter(!is.na(.data[[var]]), !is.na(Brain_adap), !is.na(primary_site)) %>%
    mutate(group = as.factor(.data[[var]]))

  p <- ggplot(df_plot, aes(x = group, y = Brain_adap, fill = group)) +
    geom_boxplot(outlier.size = 0.8, linewidth = 0.4, alpha = 0.8) +
    geom_jitter(width = 0.15, size = 1, alpha = 0.5, color = "grey20") +
    stat_compare_means(
      method       = "wilcox.test",
      label        = "p.format",
      size         = 3.5,
      comparisons  = combn(levels(as.factor(df_plot$group)), 2, simplify = FALSE)
    ) +
    facet_wrap(~ primary_site, scales = "free_x") +
    labs(
      title = "",
      x     = var,
      y     = "BrM Remodeling Score",
      fill  = var
    ) +
    theme_classic(base_size = 11, base_family = "Helvetica") +
    theme(
      axis.title       = element_text(face = "bold", size = 10),
      axis.text        = element_text(size = 9, color = "grey10"),
      axis.text.x      = element_text(angle = 30, hjust = 1),
      legend.position  = "none",
      plot.title       = element_text(hjust = 0.5, face = "bold", size = 11),
      strip.background = element_rect(fill = "grey95", colour = NA),
      strip.text       = element_text(face = "bold", size = 10)
    )

  print(p)
}

for (var in cat_vars) {
  df_plot <- community_fraction_tib_sub %>%
    filter(!is.na(.data[[var]]), !is.na(Brain_adap)) %>%
    mutate(group = as.factor(.data[[var]]))

  p <- ggplot(df_plot, aes(x = group, y = Brain_adap, fill = group)) +
    geom_boxplot(outlier.size = 0.8, linewidth = 0.4, alpha = 0.8) +
    geom_jitter(width = 0.15, size = 1, alpha = 0.5, color = "grey20") +
    stat_compare_means(
      method       = "wilcox.test",
      label        = "p.format",
      size         = 3.5,
      comparisons  = combn(levels(as.factor(df_plot$group)), 2, simplify = FALSE)
    ) +
    labs(
      title = "",
      x     = var,
      y     = "BrM Remodeling Score",
      fill  = var
    ) +
    theme_classic(base_size = 11, base_family = "Helvetica") +
    theme(
      axis.title       = element_text(face = "bold", size = 10),
      axis.text        = element_text(size = 9, color = "grey10"),
      axis.text.x      = element_text(angle = 30, hjust = 1),
      legend.position  = "none",
      plot.title       = element_text(hjust = 0.5, face = "bold", size = 11)
    )

  print(p)
}

df_plot <- community_fraction_tib_sub %>%
  filter(!is.na(.data[["Location"]]), !is.na(Brain_adap)) %>%
  mutate(group = as.factor(.data[["Location"]])) %>%
  filter(Location != "Not specified")

p <- ggplot(df_plot, aes(x = group, y = Brain_adap)) +
  geom_boxplot(outlier.size = 0.8, linewidth = 0.4, alpha = 0.8) +
  geom_jitter(width = 0.15, size = 1, alpha = 0.5, color = "grey20") +
  stat_compare_means(
    method    = "wilcox.test",
    ref.group = ".all.",          # each group vs the pooled rest
    label     = "p.format",
    label.y   = max(df_plot$Brain_adap, na.rm = TRUE) * 1.08,
    size      = 2
  ) +
  labs(
    title = "",
    x     = "Location",
    y     = "BrM adaptation score"
  ) +
  theme_classic(base_size = 7, base_family = "Helvetica") +
  theme(
    axis.text.x      = element_text(angle = 30, hjust = 1),
    legend.position  = "none",
    plot.title       = element_text(hjust = 0.5)
  )
print(p)

ggsave("BMCA/plots/fig5/loc_remo.pdf", plot = p, width = 60, height = 70,
       units = "mm", dpi = 300)

# ── Age: scatter plot with correlation ────────────────────────────────────────
df_age <- community_fraction_tib_sub %>%
  filter(!is.na(Age), !is.na(Brain_adap))

ggplot(df_age, aes(x = Age, y = Brain_adap)) +
  geom_point(size = 2, alpha = 0.7, color = "grey30") +
  geom_smooth(method = "lm", se = TRUE, color = "#C0392B", linewidth = 0.8) +
  stat_cor(method = "pearson", label.x.npc = "left", size = 3.5) +
  facet_wrap(~ primary_site) +
  labs(
    title = "",
    x     = "Age",
    y     = "BrM Remodeling Score"
  ) +
  theme_classic(base_size = 11, base_family = "Helvetica") +
  theme(
    axis.title  = element_text(face = "bold", size = 10),
    axis.text   = element_text(size = 9, color = "grey10"),
    plot.title  = element_text(hjust = 0.5, face = "bold", size = 11)
  )


km_data <- community_fraction_tib_sub %>%
  mutate(
    PFS_clean = suppressWarnings(as.numeric(PFS)),
    brain_adap_group = case_when(
      Brain_adap < 0.15   ~ "Low BrM Remodeling",
      Brain_adap > 0.5 ~ "High BrM Remodeling",
      TRUE              ~ NA_character_
    )
  ) %>%
  filter(!is.na(PFS_clean), !is.na(brain_adap_group)) %>%
  mutate(event = 1)  # ← all patients with PFS have experienced the event

fit <- survfit(Surv(PFS_clean, event) ~ brain_adap_group, data = km_data)

p <- ggsurvplot(
  fit,
  data          = km_data,
  pval          = TRUE,
  conf.int      = FALSE,
  risk.table    = FALSE,
  palette       = c("#C0392B", "#2980B9"),
  xlab          = "Time (months)",
  ylab          = "Progression Free Survival",
  legend        = c(0.8, 0.85),        # x,y position inside plot; or "top"/"bottom"/"right"
  legend.title  = "Brain adaptation",
  legend.labs   = c("High", "Low"),  # order = alphabetical factor levels
  font.legend   = 7,
  ggtheme       = theme_classic(base_size = 7)
)
p
ggsave("BMCA/plots/fig5/PFS_kaplan.pdf", plot = p$plot, width = 100, height = 60,
       units = "mm", dpi = 300)


# Check this first:

km_data_os <- community_fraction_tib_sub %>%
  distinct(patient, .keep_all = TRUE) %>%
  mutate(
    Os_month_clean = suppressWarnings(as.numeric(Os_month)),
    event_os = case_when(
      Os_status == "deceased" ~ 1,
      Os_status == "live"     ~ 0,
      TRUE                    ~ NA_real_   # catches NA and "/"
    ),
    brain_adap_group = case_when(
      Brain_adap < 0.15   ~ "Low BrM Remodeling",
      Brain_adap > 0.5    ~ "High BrM Remodeling",
      TRUE                ~ NA_character_
    )
  ) %>%
  filter(!is.na(Os_month_clean), !is.na(event_os), !is.na(brain_adap_group))

fit_os <- survfit(Surv(Os_month_clean, event_os) ~ brain_adap_group, data = km_data_os)

p_os <- ggsurvplot(
  fit_os,
  data          = km_data_os,
  pval          = TRUE,
  conf.int      = FALSE,
  risk.table    = FALSE,
  palette       = c("#C0392B", "#2980B9"),
  xlab          = "Time (months)",
  ylab          = "Overall Survival",
  legend        = c(0.8, 0.85),
  legend.title  = "BrM Remodeling",
  legend.labs   = c("High BrM Remodeling", "Low BrM Remodeling"),
  font.legend   = 7,
  ggtheme       = theme_classic(base_size = 7)
)
p_os
ggsave("BMCA/plots/fig5/OS_kaplan.pdf", plot = p_os$plot, width = 100, height = 60,
       units = "mm", dpi = 300)

km_data_os_allevent <- community_fraction_tib_sub %>%
  distinct(patient, .keep_all = TRUE) %>%
  mutate(
    Os_month_clean = suppressWarnings(as.numeric(Os_month)),
    event_os = 1,   # treat every loss to follow-up as death
    brain_adap_group = case_when(
      Brain_adap < 0.15   ~ "Low BrM Remodeling",
      Brain_adap > 0.5    ~ "High BrM Remodeling",
      TRUE                ~ NA_character_
    )
  ) %>%
  filter(!is.na(Os_month_clean), !is.na(brain_adap_group))

fit_os_allevent <- survfit(Surv(Os_month_clean, event_os) ~ brain_adap_group, data = km_data_os_allevent)

p_os_allevent <- ggsurvplot(
  fit_os_allevent,
  data          = km_data_os_allevent,
  pval          = TRUE,
  conf.int      = FALSE,
  risk.table    = FALSE,
  palette       = c("#C0392B", "#2980B9"),
  xlab          = "Time (months)",
  ylab          = "Overall Survival (all censored = death)",
  legend        = c(0.8, 0.85),
  legend.title  = "BrM Remodeling",
  legend.labs   = c("High BrM Remodeling", "Low BrM Remodeling"),
  font.legend   = 11,
  censor        = TRUE,
  ggtheme       = theme_classic(base_size = 11)
)
p_os_allevent

community_fraction_tib_sub_high <- community_fraction_tib_sub %>%
  filter(Brain_adap>0.5)


# -----------------------------------------------------------------------------
metadata_p <- metadata %>%
  dplyr::select(patient, primary_site, study) %>%
  distinct() %>%
  filter(patient %in% community_fraction_tib$patient) %>%
  filter(primary_site %in% c("breast", "NSCLC"))


ecotypes_states_mal <- ecotypes_states %>%
  dplyr::select(patient, Myeloid_Microglia, Myeloid_Scavenger_Immunosuppressive,
                Malignant, Myeloid_Stress, Malignant_Chromatin_Remodeling,
                Malignant_Neuronal_related_2, Malignant_Neuronal_related_3,
                Malignant_secretory_CRC_enriched, Malignant_PDAC_related,
                Fibroblast_C1, Fibroblast_C2, Endothelial_C1, Endothelial_C2,
                starts_with("Oligo"),
                starts_with("Neuron"),
                starts_with("Astro")) %>%
  filter(patient %in% community_fraction_tib$patient)

ecotypes_states_mal$brain_tme <- ecotypes_states_mal$Myeloid_Microglia +
  ecotypes_states_mal$Oligodendrocyte +
  ecotypes_states_mal$Astrocyte +
  ecotypes_states_mal$Neuron

ecotypes_states_mal$brain_tme_to_purity <- ecotypes_states_mal$brain_tme/
  (100-ecotypes_states_mal$Malignant)

patient_means_MNR_NS <- readRDS("BMCA/data/patient_means_MNR_NS.RDS") %>%
  filter(patient %in% community_fraction_tib$patient) %>%
  dplyr::select(patient, NTF4, BDNF)

colnames(patient_means_MNR_NS) <- c("patient", "Malignant_MNR_NTF4",
                                    "Malignant_MNR_BDNF")

scores_for_MHC <- readRDS("BMCA/data/scores_for_MHC.RDS")

APC_cell_types <- c("Myeloid", "B_Plasma", "Endothelial")

mhcii_tib <- APC_cell_types %>%
  lapply(function(ct) {
    scores_for_MHC[[ct]] %>%
      dplyr::select(patient, MHCII) %>%
      dplyr::rename(!!paste0("MHCII_", ct) := MHCII)   # dplyr::rename so := works
  }) %>%
  purrr::reduce(dplyr::full_join, by = "patient") %>%   # purrr::reduce, not IRanges::reduce
  dplyr::mutate(
    MHCII_mean_APCs = rowMeans(
      dplyr::across(dplyr::all_of(paste0("MHCII_", APC_cell_types))),
      na.rm = TRUE
    )
  ) %>%
  dplyr::arrange(patient)

ecotypes_meta <- metadata_p %>%
  left_join(community_fraction_tib, by = "patient") %>%
  left_join(ecotypes_states_mal, by = "patient") %>%
  #left_join(patient_means_MNR_NS, by = "patient") %>%
  left_join(mhcii_tib %>% dplyr::select(patient, MHCII_mean_APCs), by = "patient")
patients_w_ba <- ecotypes_meta %>% filter(Brain_adap > 0) %>% pull(patient)
saveRDS(patients_w_ba, "BMCA/data/patients_w_ba.RDS")

# -----------------------------------------------------------------------------
ecotypes_tib_proccesed_filtered <- readRDS("BMCA/data/ecotypes_tib_proccesed_filtered.RDS")
metadata_states_proccesed_filtered <- readRDS("BMCA/data/metadata_states_proccesed_filtered.RDS")
patients_wo_ba <- ecotypes_meta %>% filter(Brain_adap < 0) %>% pull(patient)
ecotypes_tib_proccesed_filtered <- ecotypes_tib_proccesed_filtered %>% filter(patient %in% patients_wo_ba)
metadata_states_proccesed_filtered <- metadata_states_proccesed_filtered %>% filter(patient %in% patients_wo_ba)

# --- 1. Define states per theme ---
theme_states <- list(
  `Antigen Presentation` = list(
    Malignant   = c("Malignant_INFR_MHC", "Malignant_MHC.I", "Malignant_EpiSen_MHC"),
    Myeloid     = c("Myeloid_DC"),
    B_Plasma    = c("B_Plasma_Antigen_presentation"),
    GD          = c("GD_MHCII_Dysfunctional"),
    CD8         = c("CD8_CT_MHCII_INFG"),
    Endothelial = c("Endothelial_HEV")
  ),
  `Hypoxia` = list(
    Malignant   = c("Malignant_Hypoxia"),
    Myeloid     = c("Myeloid_Hypoxia"),
    Fibroblast  = c("Fibroblast_Hypoxia")
  )
)

# --- 2. Join metadata ---
df_full <- ecotypes_tib_proccesed_filtered %>%
  left_join(metadata_states_proccesed_filtered %>% dplyr::select(patient, location),
            by = "patient") %>%
  filter(!is.na(location))

# --- 3. Compute delta + Wilcoxon per theme x cell type ---
results <- imap_dfr(theme_states, function(cell_list, theme_name) {

  imap_dfr(cell_list, function(cols, cell_type) {

    tmp <- df_full %>%
      dplyr::select(patient, location, all_of(cols)) %>%
      rowwise() %>%
      mutate(mean_val = mean(c_across(all_of(cols)), na.rm = TRUE)) %>%
      ungroup() %>%
      filter(!is.na(mean_val))

    if (length(unique(tmp$location)) < 2) return(NULL)

    means <- tmp %>%
      group_by(location) %>%
      summarise(mean_val = mean(mean_val, na.rm = TRUE), .groups = "drop")

    intra <- means %>% filter(location == "Intracranial") %>% pull(mean_val)
    extra <- means %>% filter(location == "Extracranial") %>% pull(mean_val)

    if (length(intra) == 0 || length(extra) == 0) return(NULL)

    test <- wilcox.test(mean_val ~ location, data = tmp, exact = FALSE)

    tibble(
      theme     = theme_name,
      cell_type = cell_type,
      delta     = intra - extra,
      p_value   = test$p.value
    )
  })
}) %>%
  mutate(
    p_adj = p.adjust(p_value, method = "BH"),
    stars = case_when(
      p_adj < 0.001 ~ "***",
      p_adj < 0.01  ~ "**",
      p_adj < 0.05  ~ "*",
      TRUE          ~ ""
    )
  ) %>%
  # Order cell_type within each theme by delta
  group_by(theme) %>%
  mutate(cell_type = fct_reorder(cell_type, delta)) %>%
  ungroup()
results <- results %>% filter(theme == 'Antigen Presentation')
# --- 4. Plot ---
ggplot(results, aes(x = cell_type, y = delta)) +
  geom_bar(stat = "identity", fill = "steelblue", color = "black", linewidth = 0.4) +
  geom_hline(yintercept = 0, linewidth = 0.5) +
  geom_text(
    aes(
      label = stars,
      y     = delta + ifelse(delta >= 0, 0.3, -1)
    ),
    size = 5, vjust = 0.5
  ) +
  facet_wrap(~ theme, scales = "free_x") +
  labs(
    title = "",
    x     = "Antigen-Presentation Transcriptional Programs per Cell Type",
    y     = "Mean Intracranial − Mean Extracranial"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x        = element_text(angle = 45, hjust = 1),
    panel.grid.major.x = element_blank(),
    plot.title         = element_text(hjust = 0.5, face = "bold"),
    strip.text         = element_text(face = "bold", size = 12)
  )

# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
ecotypes_meta <- metadata_p %>%
  left_join(community_fraction_tib, by = "patient") %>%
  left_join(ecotypes_states_mal, by = "patient") %>%
 # left_join(patient_means_fibro, by = "patient") %>%
  #left_join(patient_means_MNR_NS, by = "patient") %>%
  left_join(mhcii_tib %>% dplyr::select(patient, MHCII_mean_APCs), by = "patient")


ecotypes_meta <- ecotypes_meta %>%
  filter(study == "lin_sn" & primary_site %in% c("NSCLC", "breast")) %>%
  dplyr::select(-primary_site, -study)


colnames(ecotypes_meta) <- c("patient",
                             "Brain adap. score (mean)",
                             "Microglia", "Myelo. immunosupp.",
                             "Malignant", "Myelo. stress",
                             "Mal. chromatin", "MNR2",
                             "MNR3", "Mal. epi. diff. 1",
                             "Mal. epi. diff. 2", "Fib_C1", "Fib_C2", "Endo_C1",
                             "Endo_C2",
                             "Oligodendrocyte", "Mature", "Oligo1", "Oligo2",
                             "Stress response", "Neuron", "Excitatory", "Inhibitory",
                             "Astrocyte", "Reactive", "Non-Reactive", "brain_tme",
                             "brain_tme_to_purity",
                             "MHCII (APCs)"
                             )

ecotypes_meta <- ecotypes_meta %>%
  dplyr::select(patient,
                "Myelo. stress",
                'Mal. chromatin',
                MNR2,
                "Myelo. immunosupp.",
                MNR3,
                "Mal. epi. diff. 1",
                "Mal. epi. diff. 2",
                'Fib_C1',
                'MHCII (APCs)')

# ----- Select numeric columns -----
mat <- ecotypes_meta %>%
  dplyr::select(2:ncol(.)) %>%
  as.data.frame()

n          <- ncol(mat)
feat_names <- colnames(mat)

# ----- Spearman correlation + p-value -----
cor_mat <- matrix(NA, n, n, dimnames = list(feat_names, feat_names))
p_mat   <- matrix(NA, n, n, dimnames = list(feat_names, feat_names))

for (i in seq_len(n)) {
  for (j in seq_len(n)) {
    x    <- mat[[i]]
    y    <- mat[[j]]
    keep <- complete.cases(x, y)
    if (sum(keep) >= 2) {
      test          <- cor.test(x[keep], y[keep], method = "spearman")
      cor_mat[i, j] <- test$estimate
      p_mat[i, j]   <- test$p.value
    }
  }
}

# ----- BH correction on upper triangle only -----
upper_idx           <- which(upper.tri(p_mat))
p_adjusted          <- p.adjust(p_mat[upper_idx], method = "BH")

p_mat_adj           <- matrix(NA, n, n, dimnames = list(feat_names, feat_names))
p_mat_adj[upper_idx] <- p_adjusted
p_mat_adj[lower.tri(p_mat_adj)] <- t(p_mat_adj)[lower.tri(p_mat_adj)]

# ----- Significance labels (BH-adjusted) -----
sig_mat <- matrix("", n, n, dimnames = list(feat_names, feat_names))
sig_mat[!is.na(p_mat_adj) & p_mat_adj < 0.05]  <- "*"
sig_mat[!is.na(p_mat_adj) & p_mat_adj < 0.01]  <- "**"
sig_mat[!is.na(p_mat_adj) & p_mat_adj < 0.001] <- "***"
diag(sig_mat) <- ""

# ----- Hierarchical clustering -----
hc      <- hclust(dist(cor_mat))
ord     <- hc$order
cor_ord <- cor_mat[ord, ord]
sig_ord <- sig_mat[ord, ord]

# ----- Build plot_df -----
plot_df <- expand.grid(
  Var1 = feat_names[ord],
  Var2 = feat_names[ord],
  stringsAsFactors = FALSE
) %>%
  mutate(
    Correlation = as.vector(cor_ord),
    Sig         = as.vector(sig_ord)
  )

plot_df$Var1 <- factor(plot_df$Var1, levels = feat_names[ord])
plot_df$Var2 <- factor(plot_df$Var2, levels = feat_names[ord])


# ----- Plot -----
p <- ggplot(plot_df, aes(x = Var2, y = Var1, fill = Correlation)) +
  geom_tile(color = "grey85", linewidth = 0.3) +
  geom_text(
    data     = subset(plot_df, Sig != ""),
    aes(label = Sig),
    size     = 2,
    fontface = "bold"
  ) +
  scale_fill_gradient2(
    low      = "darkblue",
    mid      = "white",
    high     = "darkred",
    midpoint = 0,
    limits   = c(-1, 1)
  ) +
  coord_fixed() +
  theme_minimal(base_size = 7) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.title  = element_blank(),
    panel.grid  = element_blank()
  )

p

ggsave(file.path( "BMCA/plots/fig5/heatmap_brm_priamry.pdf"), plot = p,
       units = "mm", width = 70, height = 70, dpi = 300)


ecotypes_meta <- metadata_p %>%
  left_join(community_fraction_tib, by = "patient") %>%
  left_join(ecotypes_states_mal, by = "patient") %>%
  left_join(patient_means_fibro, by = "patient") %>%
  #left_join(patient_means_MNR_NS, by = "patient") %>%
  left_join(mhcii_tib %>% dplyr::select(patient, MHCII_mean_APCs), by = "patient")


ecotypes_meta <- ecotypes_meta %>%
  filter(study == "lin_sn" & primary_site %in% c("NSCLC", "breast")) %>%
  dplyr::select(-primary_site, -study, -Fibroblast_NTRK2, -Malignant,
                -Fibroblast_C2, -Endothelial_C1, -Endothelial_C2,
                -Malignant_secretory_CRC_enriched, -Malignant_PDAC_related,
                -MHCII_mean_APCs, -brain_tme, -brain_tme_to_purity, -Fibroblast_NTRK2,
                -Fibroblast_NTRK3)


colnames(ecotypes_meta) <- c("patient",
                             "Brain adap. score",
                             "Microglia", "Myelo. immunosupp.",
                             "Myelo. stress",
                             "Mal. chromatin", "MNR2",
                             "MNR3", "Neuronal rel. CAFs",
                             "Oligodendrocyte", "Mature", "Oligo1", "Oligo2",
                             "Stress response", "Neuron", "Excitatory", "Inhibitory",
                             "Astrocyte", "Non-Reactive", "Reactive"
)

# ----- Select numeric columns -----
mat <- ecotypes_meta %>%
  dplyr::select(2:ncol(.)) %>%
  as.data.frame()

n          <- ncol(mat)
feat_names <- colnames(mat)

# ----- Spearman correlation + p-value -----
cor_mat <- matrix(NA, n, n, dimnames = list(feat_names, feat_names))
p_mat   <- matrix(NA, n, n, dimnames = list(feat_names, feat_names))

for (i in seq_len(n)) {
  for (j in seq_len(n)) {
    x    <- mat[[i]]
    y    <- mat[[j]]
    keep <- complete.cases(x, y)
    if (sum(keep) >= 2) {
      test          <- cor.test(x[keep], y[keep], method = "spearman")
      cor_mat[i, j] <- test$estimate
      p_mat[i, j]   <- test$p.value
    }
  }
}

# ----- BH correction on upper triangle only -----
upper_idx           <- which(upper.tri(p_mat))
p_adjusted          <- p.adjust(p_mat[upper_idx], method = "BH")

p_mat_adj           <- matrix(NA, n, n, dimnames = list(feat_names, feat_names))
p_mat_adj[upper_idx] <- p_adjusted
p_mat_adj[lower.tri(p_mat_adj)] <- t(p_mat_adj)[lower.tri(p_mat_adj)]

# ----- Significance labels (BH-adjusted) -----
sig_mat <- matrix("", n, n, dimnames = list(feat_names, feat_names))
sig_mat[!is.na(p_mat_adj) & p_mat_adj < 0.05]  <- "*"
sig_mat[!is.na(p_mat_adj) & p_mat_adj < 0.01]  <- "**"
sig_mat[!is.na(p_mat_adj) & p_mat_adj < 0.001] <- "***"
diag(sig_mat) <- ""

# ----- Hierarchical clustering -----
hc      <- hclust(dist(cor_mat))
ord     <- hc$order
cor_ord <- cor_mat[ord, ord]
sig_ord <- sig_mat[ord, ord]

# ----- Build plot_df -----
plot_df <- expand.grid(
  Var1 = feat_names[ord],
  Var2 = feat_names[ord],
  stringsAsFactors = FALSE
) %>%
  mutate(
    Correlation = as.vector(cor_ord),
    Sig         = as.vector(sig_ord)
  )

plot_df$Var1 <- factor(plot_df$Var1, levels = feat_names[ord])
plot_df$Var2 <- factor(plot_df$Var2, levels = feat_names[ord])

# ---- group definitions (order matters) ----
col_groups <- tribble(
  ~feature,                            ~col_group,
  "Oligodendrocyte",                   "Brain cell types",
  "Neuron",                            "Brain cell types",
  "Microglia",                 "Brain cell types",
  "Astrocyte",                         "Brain cell types",
  "Inhibitory",                 "Neuron subtypes",
  "Excitatory",                 "Neuron subtypes",
  "Reactive",                "Astrocyte subtypes",
  "Non-Reactive",             "Astrocyte subtypes",
  "Stress response",                 "Oligodendrocyte subtypes",
  "Mature",                 "Oligodendrocyte subtypes",
  "Oligo2",                 "Oligodendrocyte subtypes",
  "Oligo1",                 "Oligodendrocyte subtypes"
)

row_groups <- tribble(
  ~feature,                               ~row_group,
  "MNR2",         "Features",
  "Neuronal rel. CAFs",                     "Features",
  "Mal. chromatin",       "Features",
  "MNR3",         "Features",
  "Myelo. immunosupp.",  "Features",
  "Myelo. stress",                       "Features",
  "Brain adap. score",                           "Brain adaptation"
)

col_group_levels <- c("Brain cell types", "Neuron subtypes",
                      "Astrocyte subtypes", "Oligodendrocyte subtypes")
row_group_levels <- c("Features", "Brain adaptation")   # top band, bottom band

row_feats <- row_groups$feature
col_feats <- col_groups$feature

# ---- features must exist in cor_mat ----
missing <- setdiff(c(row_feats, col_feats), rownames(cor_mat))
if (length(missing) > 0)
  stop("These aren't in cor_mat: ", paste(missing, collapse = ", "))

# ---- subset precomputed matrices to (rows x cols) ----
cor_sub <- cor_mat[row_feats, col_feats, drop = FALSE]
sig_sub <- sig_mat[row_feats, col_feats, drop = FALSE]

plot_df <- as.data.frame(cor_sub) %>%
  rownames_to_column("row_feat") %>%
  pivot_longer(-row_feat, names_to = "col_feat", values_to = "Correlation") %>%
  left_join(
    as.data.frame(sig_sub) %>%
      rownames_to_column("row_feat") %>%
      pivot_longer(-row_feat, names_to = "col_feat", values_to = "Sig"),
    by = c("row_feat", "col_feat")
  ) %>%
  left_join(row_groups, by = c("row_feat" = "feature")) %>%
  left_join(col_groups, by = c("col_feat" = "feature")) %>%
  mutate(
    row_feat  = factor(row_feat,  levels = rev(row_feats)),  # first listed -> top
    col_feat  = factor(col_feat,  levels = col_feats),
    row_group = factor(row_group, levels = row_group_levels),
    col_group = factor(col_group, levels = col_group_levels)
  )

# ---- plot ----
ggplot(plot_df, aes(x = col_feat, y = row_feat, fill = Correlation)) +
  geom_tile(colour = "grey85", linewidth = 0.3) +
  geom_text(aes(label = Sig), size = 2, fontface = "bold") +
  scale_fill_gradient2(low = "darkblue", mid = "white", high = "darkred",
                       midpoint = 0, limits = c(-1, 1)) +
  facet_grid(row_group ~ col_group,
             scales = "free", space = "free", switch = "y",
             labeller = labeller(col_group = label_wrap_gen(width = 12))) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 7) +
  theme(
    axis.text.x       = element_text(angle = 45, hjust = 1, vjust = 1),
    panel.grid        = element_blank(),
    panel.spacing     = unit(4, "pt"),
    strip.text.y.left = element_blank(),
    strip.placement   = "outside",
    strip.background  = element_rect(fill = "grey95", colour = NA)
  )

p <- ggplot(plot_df, aes(x = col_feat, y = row_feat, fill = Correlation)) +
  geom_tile(colour = "grey85", linewidth = 0.3) +
  geom_text(aes(label = Sig), size = 2, fontface = "bold") +
  scale_fill_gradient2(low = "darkblue", mid = "white", high = "darkred",
                       midpoint = 0, limits = c(-1, 1)) +
  facet_grid(row_group ~ col_group,
             scales = "free", space = "free", switch = "y",
             labeller = labeller(col_group = label_wrap_gen(width = 12))) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 7) +
  theme(
    axis.text.x       = element_text(angle = 45, hjust = 1, vjust = 1),
    panel.grid        = element_blank(),
    panel.spacing     = unit(4, "pt"),
    strip.text.y.left = element_blank(),
    strip.placement   = "outside",
    strip.background  = element_rect(fill = "grey95", colour = NA)
  )   # <- note: your snippet was missing this closing paren on theme() + the plot

p
ggsave(file.path( "BMCA/plots/fig5/heatmap_brain_enrch_brain_cells.pdf"), plot = p,
       units = "mm", width = 120, height = 60, dpi = 300)
# -----------------------------------------------------------------------------
library(Matrix); library(dplyr); library(tibble); library(purrr)
library(ggplot2); library(ggrepel); library(patchwork)

## ---- config -----------------------------------------------------------------
data_dir   <- "BMCA/data"
fig_dir    <- "BMCA/plots"
csv_dir    <- "BMCA/csv"
study      <- "lin_sn"
labels     <- c("Neuron", "Astrocyte", "Oligodendrocyte")
cell_types <- c("Neuron", "Astrocyte", "Oligodendrocyte")
hi_co      <- 0.5     # high Brain_adap
low_co     <- 0.15    # low Brain_adap
LOG2FC_CO  <- 2
P_CO       <- 0.05
PADJ_CO    <- 0.05
LABEL_BY   <- "padj"  # "p" | "padj"  -> volcano y-axis and significance call
N_TOP      <- 10      # top genes per direction, ranked by log2FC only

list_of_pbs <- list()
for (ct in labels) {
  print(ct)
  pb <- readRDS(paste0("BMCA/data/UMI_pb_", ct, ".RDS"))[["lin_sn"]]
  list_of_pbs[[ct]] <- pb
}


# 1. CPM-normalise each pseudobulk sample (column), then log2
norm_pb <- lapply(list_of_pbs, function(m) {
  m   <- as.matrix(m)
  lib <- colSums(m)
  cpm <- t(t(m) / lib) * 1e6      # each column sums to 1e6
  log2(cpm + 1)                    # use log2(cpm/10 + 1) to down-weight low counts
})

# 2. Per-gene mean log2-CPM within each cell type  -> genes x 13
ct_mean <- sapply(norm_pb, rowMeans)

# 3. Build per-gene df for each cell type: avg expr + log2FC vs other cell types
cell_types <- colnames(ct_mean)
plot_df <- bind_rows(lapply(cell_types, function(ct) {
  others <- setdiff(cell_types, ct)
  data.frame(
    gene      = rownames(ct_mean),
    cell_type = ct,
    avg_expr  = ct_mean[, ct],                                   # x  (this cell type)
    # avg_expr = rowMeans(ct_mean),                              # <- grand-mean MA-plot version
    log2fc    = ct_mean[, ct] - rowMeans(ct_mean[, others, drop = FALSE])  # y
  )
}))

# 4. One faceted overview figure
ggplot(plot_df, aes(avg_expr, log2fc)) +
  geom_point(size = 0.00001, alpha = 0.1) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey50") +
  facet_wrap(~ cell_type, scales = "free_x") +
  labs(x = "Mean expression (log2 CPM)",
       y = "log2FC vs other cell types") +
  theme_bw(base_size = 20)

top_genes_n <- 500

for (ct in cell_types) {
  print(ct)
  d <- dplyr::filter(plot_df, cell_type == ct)

  # x threshold: value leaving 1000 genes to its right
  x_cut <- sort(d$avg_expr, decreasing = TRUE)[top_genes_n]

  # y threshold: value leaving 1000 genes above it
  y_cut <- sort(d$log2fc, decreasing = TRUE)[top_genes_n]
  dim(d %>% filter(log2fc > y_cut))
  # top 20 genes by log2FC
  top20 <- d %>% dplyr::slice_max(log2fc, n = 50)

  p <- ggplot(d, aes(avg_expr, log2fc)) +
    geom_point(size = 1, alpha = 0.35, colour = "grey40") +
    geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey60") +
    geom_vline(xintercept = x_cut, linewidth = 0.4,
               linetype = "dashed", colour = "firebrick") +
    geom_vline(xintercept = 2, linewidth = 0.4,
               linetype = "dashed", colour = "firebrick") +
    geom_hline(yintercept = y_cut, linewidth = 0.4,
               linetype = "dashed", colour = "steelblue") +
    geom_point(data = top20, size = 0.9, colour = "firebrick") +
    ggrepel::geom_text_repel(
      data = top20, aes(label = gene),
      size = 2.5, max.overlaps = Inf,
      segment.size = 0.2, segment.colour = "grey50",
      box.padding = 0.3, min.segment.length = 0
    ) +
    labs(title = ct,
         x = "Mean expression (log2 CPM)",
         y = "log2FC vs other cell types") +
    theme_bw(base_size = 10)

  print(p)
}

top_genes <- lapply(cell_types, function(ct) {
  d <- dplyr::filter(plot_df, cell_type == ct)
  d <- d[order(d$log2fc, decreasing = TRUE), ]
  d <- head(d, top_genes_n)                    # top 1000 by log2FC first
  d <- dplyr::filter(d, avg_expr > 3)   # then drop low-expression
  d$gene
})
names(top_genes) <- labels

hodge_short <- readr::read_csv("BMCA/csv/hodge2019_markers_neurons.csv")
oligo_marques <- readr::read_csv("BMCA/csv/oligo_marques.csv")
astro_nmf <- readr::read_csv("BMCA/csv/sigs_all_Astrocyte_2.csv")

## ---- neurotransmitter synthesis genes, one colour each ----------------------
nt_genes <- c(
)
names(nt_genes)
genes_neuro <- unlist(hodge_short, use.names = FALSE)
genes_oligo <- c()
genes_astro <- unlist(astro_nmf, use.names = FALSE)

genes_neuro <- unique(c(genes_neuro, top_genes[["Neuron"]], names(nt_genes)))
genes_oligo <- unique(c(genes_oligo, top_genes[["Oligodendrocyte"]]))
genes_astro <- unique(c(genes_astro, top_genes[["Astrocyte"]]))

marker_genes_by_pair <- list(
  Neuron      = genes_neuro,
  Astrocyte  = genes_astro,
  Oligodendrocyte = genes_oligo
)

## ---- Brain_adap + (optional) strata per patient ----------------------------
#metadata <- as_tibble(readRDS(file.path(data_dir, "metadata_all_studies.rds")))
cohort <- readRDS(file.path(data_dir, "community_fraction_tib.RDS")) %>%
  left_join(metadata_p, by = "patient") %>%
  filter(primary_site != "colorectal", !is.na(Brain_adap)) %>%
  mutate(Brain_adap = as.numeric(Brain_adap))

uo_co <- 0.5 #high brain adap
low_co <- 0.15 #low brain adap

## short prefixes for the output tibble column names
ct_short <- c(Neuron = "neuro", Astrocyte = "astro", Oligodendrocyte = "oligo")

## ---- clean marker lists: uppercase, drop blanks -----------------------------
clean_genes <- function(x) {
  x <- toupper(as.character(unlist(x, use.names = FALSE)))
  x <- x[!is.na(x) & nzchar(x) & x != "NA"]
  unique(x)
}
marker_genes_by_pair <- lapply(marker_genes_by_pair, clean_genes)

## ---- core DE for one cell type ---------------------------------------------
run_de_ct <- function(ct) {
  pb <- list_of_pbs[[ct]]

  keep_pt <- intersect(colnames(pb), cohort$patient)
  pb <- pb[, keep_pt, drop = FALSE]
  bs <- cohort$Brain_adap[match(keep_pt, cohort$patient)]

  hi <- which(bs >= hi_co)
  lo <- which(bs <= low_co)
  message(sprintf("%-16s matched=%2d  high=%2d  low=%2d",
                  ct, length(keep_pt), length(hi), length(lo)))
  if (length(hi) < 2 || length(lo) < 2) { message("  -> skipping"); return(NULL) }

  ## 1. CPM
  cpm <- Matrix::t(Matrix::t(pb) / Matrix::colSums(pb)) * 1e6

  ## 2. filter to this cell type's own marker genes
  g_keep <- intersect(rownames(cpm), marker_genes_by_pair[[ct]])
  message(sprintf("  -> %d/%d marker genes found",
                  length(g_keep), length(marker_genes_by_pair[[ct]])))
  cpm <- as.matrix(cpm[g_keep, , drop = FALSE])

  ## 3. log2(x + 1)
  lg <- log2(cpm + 1)
  lg <- lg[matrixStats::rowMaxs(lg) > 0, , drop = FALSE]   # drop all-zero genes only

  ## 4. log2FC + Wilcoxon
  mh <- rowMeans(lg[, hi, drop = FALSE])
  ml <- rowMeans(lg[, lo, drop = FALSE])
  pv <- apply(lg, 1, function(v)
    suppressWarnings(wilcox.test(v[hi], v[lo], exact = FALSE)$p.value))

  res <- tibble(
    gene      = rownames(lg),
    cell_type = ct,
    n_hi      = length(hi),
    n_lo      = length(lo),
    mean_hi   = mh,
    mean_lo   = ml,
    log2FC    = mh - ml,
    p         = pv,
    padj      = p.adjust(pv, method = "BH")
  )

  yv <- if (LABEL_BY == "padj") res$padj else res$p
  co <- if (LABEL_BY == "padj") PADJ_CO else P_CO
  res$sig <- !is.na(yv) & yv < co & abs(res$log2FC) >= LOG2FC_CO
  res$direction <- ifelse(!res$sig, "ns", ifelse(res$log2FC > 0, "high_BA", "low_BA"))

  ## ---- top N per side by log2FC ONLY (ignores p / padj) ---------------------
  top_hi <- res$gene[order(res$log2FC, decreasing = TRUE)][seq_len(min(N_TOP, sum(res$log2FC > 0)))]
  top_lo <- res$gene[order(res$log2FC, decreasing = FALSE)][seq_len(min(N_TOP, sum(res$log2FC < 0)))]
  res$top_side <- ifelse(res$gene %in% top_hi, "high_BA",
                         ifelse(res$gene %in% top_lo, "low_BA", NA_character_))

  res <- arrange(res, desc(log2FC))
  message(sprintf("  -> %d sig by %s+FC (%d up high-BA, %d up low-BA) | labelling %d + %d by FC",
                  sum(res$sig), LABEL_BY,
                  sum(res$direction == "high_BA"), sum(res$direction == "low_BA"),
                  length(top_hi), length(top_lo)))
  res
}

de_res <- lapply(labels, run_de_ct)
names(de_res) <- labels
de_res <- de_res[!vapply(de_res, is.null, logical(1))]
de_all <- bind_rows(de_res)

## ---- 50 x 6 tibble: top N genes per cell type per direction -----------------
FC_CO <- 2   # |log2FC| threshold for the gene table

pad_to <- function(x, n) c(x, rep(NA_character_, max(0, n - length(x))))[seq_len(n)]

gene_sets <- map(names(de_res), function(ct) {
  r   <- de_res[[ct]]
  pre <- ct_short[[ct]]
  hi  <- r %>% filter(log2FC >=  FC_CO) %>% arrange(desc(log2FC)) %>% pull(gene)
  lo  <- r %>% filter(log2FC <= -FC_CO) %>% arrange(log2FC)       %>% pull(gene)
  set_names(list(hi, lo), paste0(pre, c("_high", "_low")))
}) %>% flatten()

n_max <- max(lengths(gene_sets), 1)
top50_tbl <- as_tibble(map(gene_sets, pad_to, n = n_max))

print(map_int(gene_sets, length))          # how many genes per column
print(top50_tbl, n = 15)

readr::write_csv(top50_tbl,
                 file.path(csv_dir, paste0("genes_absFC", FC_CO, "_brainadap_", study, ".csv")))

## ---- neurotransmitter synthesis genes, one colour each ----------------------
nt_genes <- c(
)

nt_cols <- c(                    # histamine - teal
)

## legend labels: "GAD1 (GABA)"
nt_labels <- setNames(sprintf("%s  (%s)", names(nt_genes), unname(nt_genes)),
                      names(nt_genes))

FORCE_LABEL <- list(
  Neuron = names(nt_genes)
)

## ---- volcano ----------------------------------------------------------------
volcano_ct <- function(res) {
  ct <- res$cell_type[1]
  res$y <- -log10(if (LABEL_BY == "padj") res$padj else res$p)
  res$y[!is.finite(res$y)] <- NA
  hline <- -log10(if (LABEL_BY == "padj") PADJ_CO else P_CO)

  want   <- FORCE_LABEL[[ct]]
  if (is.null(want)) want <- character(0)
  forced <- intersect(want, res$gene)
  if (length(want) && length(setdiff(want, forced))) {
    message(sprintf("  %s: forced genes absent from results -> %s",
                    ct, paste(setdiff(want, forced), collapse = ", ")))
  }

  lab_top <- filter(res, !is.na(top_side), !gene %in% forced)
  lab_frc <- filter(res, gene %in% forced)
  lab_frc$gene <- factor(lab_frc$gene, levels = names(nt_cols))

  p <- ggplot(res, aes(log2FC, y)) +
    geom_vline(xintercept = c(-LOG2FC_CO, LOG2FC_CO), linetype = "dashed",
               colour = "grey60", linewidth = 0.3) +
    geom_hline(yintercept = hline, linetype = "dashed",
               colour = "grey60", linewidth = 0.3) +
    geom_point(data = filter(res, !sig), colour = "grey85", size = 0.9, alpha = 0.6) +
    geom_point(data = filter(res,  sig), colour = "black",  size = 1.4, alpha = 0.9) +
    geom_point(data = lab_top, shape = 21, size = 1.6, stroke = 0.3,
               colour = "black", fill = NA) +
    geom_text_repel(data = lab_top, aes(label = gene), size = 1.9,
                    max.overlaps = Inf, min.segment.length = 0,
                    segment.size = 0.15, segment.colour = "grey50",
                    box.padding = 0.15, point.padding = 0.1, force = 2)

  if (nrow(lab_frc)) {
    p <- p +
      geom_point(data = lab_frc, aes(colour = gene),
                 shape = 19, size = 2.6, alpha = 0.95) +
      geom_point(data = lab_frc, shape = 21, size = 2.6, stroke = 0.4,
                 colour = "black", fill = NA) +
      geom_text_repel(data = lab_frc, aes(label = gene, colour = gene),
                      size = 3.1, fontface = "bold",
                      max.overlaps = Inf, min.segment.length = 0,
                      segment.size = 0.35, box.padding = 0.5,
                      point.padding = 0.25, force = 10, seed = 42,
                      show.legend = FALSE) +
      scale_colour_manual(values = nt_cols, labels = nt_labels,
                          breaks = levels(droplevels(lab_frc$gene)),
                          name = "NT synthesis")
  }

  p +
    labs(title    = ct,
         x = "log2FC  (high - low Brain_adap)",
         y = paste0("-log10(", LABEL_BY, ")")) +
    theme_classic(base_size = 10) +
    theme(legend.position = "right",
          legend.key.size = unit(0.4, "cm"),
          legend.text  = element_text(size = 7),
          legend.title = element_text(size = 8))
}

plots <- lapply(de_res, volcano_ct)
plots

plots_tog <- plots[[1]] + plots[[2]] + plots[[3]]

ggsave("BMCA/plots/fig5/plots_tog.pdf", plots_tog,
       width = 250, height = 80, units = "mm",
       dpi = 300)

saveRDS(de_res, file.path(data_dir, paste0("de_brainadap_by_ct_", study, ".RDS")))
enriched_genes_BA <- top50_tbl
write.csv(enriched_genes_BA, file.path(csv_dir, "top_per_cellstate_by_resid_ct.csv"), row.names = FALSE)
# -----------------------------------------------------------------------------

# ---- 1. Databases --------------------------------------------------
databases <- c("GO_Biological_Process_2026", "MSigDB_Hallmark_2020",
               "GO_Cellular_Component_2026",
               "Reactome_2022")

# ---- Fixed color map: one color per database, defined once ---------
db_colors <- setNames(
  brewer.pal(max(3, length(databases)), "Set1")[seq_along(databases)],
  databases
)

out_dir <- "BMCA/plots"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- 2. Wide tibble -> named list of gene vectors ------------------
gene_lists <- enriched_genes_BA %>%
  as.list() %>%
  map(~ .x[!is.na(.x) & .x != ""])
map_int(gene_lists, length)   # sanity check

# ---- 3. Run enrichr ONCE per cell type -----------------------------
enrich_results <- map(gene_lists, ~ enrichr(.x, databases))

# ---- 4. Manual plot titles -----------------------------------------
# Keys = names(enrich_results); values = display titles. Edit freely.
plot_titles <- c(
  "neuro_high" = "neuro_high",
  "neuro_low"  = "neuro_low",
  "astro_high" = "astro_high",
  "astro_low"  = "astro_low",
  "oligo_high" = "oligo_high",
  "oligo_low"  = "oligo_low"
)

# ---- 5. Combined data frame across all cell types ------------------
n_terms <- 5

plot_df <- imap_dfr(enrich_results, function(enrich_list, cell_type) {
  bind_rows(enrich_list, .id = "Database") %>%
    mutate(
      CellType      = cell_type,
      Database      = factor(Database, levels = databases),
      Count         = as.numeric(gsub("/.*", "", Overlap)),
      neg_log10_fdr = -log10(Adjusted.P.value),
      gene_labels   = str_replace_all(Genes, ";", ", ")
    ) %>%
    arrange(Adjusted.P.value) %>%
    slice_head(n = n_terms)
}) %>%
  mutate(CellType = factor(CellType, levels = names(enrich_results)))

GENE_SIZE      <- 12
GENES_2ROW_CO  <- 4     # more than this -> 2 rows
GENES_3ROW_CO  <- 10    # more than this -> 3 rows
GENES_4ROW_CO  <- 14    # more than this -> 4 rows

wrap_genes <- function(x, co2 = GENES_2ROW_CO, co3 = GENES_3ROW_CO,
                       co4 = GENES_4ROW_CO) {
  vapply(x, function(s) {
    if (is.na(s) || !nzchar(s)) return("")
    g <- trimws(strsplit(s, ",")[[1]])
    g <- g[nzchar(g)]
    n <- length(g)
    nrow_out <- if (n > co4) 4L else if (n > co3) 3L else if (n > co2) 2L else 1L
    if (nrow_out == 1L) return(paste(g, collapse = ", "))
    idx <- cut(seq_len(n), breaks = nrow_out, labels = FALSE)   # even split
    paste(vapply(split(g, idx), paste, character(1), collapse = ", "),
          collapse = "\n")
  }, character(1), USE.NAMES = FALSE)
}

make_plot <- function(df_ct, ct) {
  ttl <- if (ct %in% names(plot_titles)) plot_titles[[ct]] else ct
  df_ct$gene_labels_wrapped <- wrap_genes(df_ct$gene_labels)

  ggplot(df_ct,
         aes(x = reorder(Term, neg_log10_fdr),
             y = neg_log10_fdr, fill = Database)) +
    geom_bar(stat = "identity", alpha = 0.8) +
    geom_text(aes(label = gene_labels_wrapped),
              y = 0, hjust = 0, vjust = 0.5,
              size = GENE_SIZE / .pt, lineheight = 0.9,
              color = "black", fontface = "italic") +
    coord_flip(clip = "off") +
    scale_x_discrete(labels = function(x) str_wrap(x, width = WRAP_W)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    scale_fill_manual(values = db_colors, drop = FALSE, limits = databases) +
    theme_minimal(base_size = 14) +
    theme(
      text         = element_text(size = 14),
      axis.text    = element_text(size = 14),
      axis.title   = element_text(size = 14),
      axis.text.y  = element_text(size = 14, lineheight = 0.85),
      legend.text  = element_text(size = 14),
      legend.title = element_text(size = 14),
      plot.title   = element_text(size = 14, face = "bold")
    ) +
    labs(title = ttl,
         x = "Enrichment Term",
         y = "-log10(FDR Adjusted P-value)",
         fill = "Database")
}

plot_list <- plot_df %>%
  group_split(CellType) %>%
  set_names(map_chr(., ~ as.character(.x$CellType[1]))) %>%
  imap(~ make_plot(.x, .y))

plot_list

hodge_short <- readr::read_csv("BMCA/csv/hodge2019_markers_neurons.csv")
oligo_marques <- readr::read_csv("BMCA/csv/oligo_marques.csv")
astro_nmf <- readr::read_csv("BMCA/csv/sigs_all_Astrocyte_2.csv")


sigs <- as_tibble(cbind(hodge_short, oligo_marques, astro_nmf))
enriched_genes_BA


# --- turn each column into a clean gene set (drop NA, uppercase, dedupe) ---
to_sets <- function(tb) {
  lapply(tb, function(col) {
    g <- toupper(col[!is.na(col)])
    unique(g)
  })
}

sets_sigs <- to_sets(sigs)
sets_ba   <- to_sets(enriched_genes_BA)

jaccard <- function(a, b) {
  u <- length(union(a, b))
  if (u == 0) return(NA_real_)
  length(intersect(a, b)) / u
}

# --- pairwise Jaccard: sigs (rows) x enriched_genes_BA (cols) ---
jac <- expand_grid(sig = names(sets_sigs), ba = names(sets_ba)) |>
  mutate(jaccard = map2_dbl(sig, ba, ~ jaccard(sets_sigs[[.x]], sets_ba[[.y]])))

# --- order both axes by hierarchical clustering on the Jaccard matrix ---
m <- jac |>
  pivot_wider(names_from = ba, values_from = jaccard) |>
  tibble::column_to_rownames("sig") |>
  as.matrix()

ord_sig <- rownames(m)[hclust(dist(m))$order]
ord_ba  <- colnames(m)[hclust(dist(t(m)))$order]

jac <- jac |>
  mutate(sig = factor(sig, levels = ord_sig),
         ba  = factor(ba,  levels = ord_ba))

n_shared <- function(a, b) length(intersect(a, b))

jac <- expand_grid(sig = names(sets_sigs), ba = names(sets_ba)) |>
  mutate(jaccard = map2_dbl(sig, ba, ~ jaccard(sets_sigs[[.x]], sets_ba[[.y]])),
         shared  = map2_int(sig, ba, ~ n_shared(sets_sigs[[.x]], sets_ba[[.y]])))

p <- ggplot(jac, aes(ba, sig, fill = jaccard)) +
  geom_tile(colour = "grey92") +
  geom_text(aes(label = shared), size = 2.6) +
  scale_fill_gradient(low = "white", high = "#1b5e20",
                      na.value = "grey95", name = "Jaccard") +
  coord_fixed() +
  labs(x = NULL, y = NULL) +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p

library(tidyverse); library(cowplot); library(patchwork)

BASE <- 6

LINEAGE <- tribble(
  ~lineage,          ~sig_pre,           ~ba_pre,
  "Neuron",          "Neuron",           "neuro",
  "Oligodendrocyte", "Olig",  "oligo",
  "Astrocyte",       "Astro",        "astro"
)

jac_l <- jac %>%
  mutate(sig_pre = str_extract(sig, "^[^_]+"),
         ba_pre  = str_extract(ba,  "^[^_]+")) %>%
  inner_join(LINEAGE, by = c("sig_pre", "ba_pre")) %>%
  mutate(lineage = factor(lineage, levels = LINEAGE$lineage),
         ba      = factor(ba, levels = c(paste0(unique(ba_pre), "_high"),
                                         paste0(unique(ba_pre), "_low"))))
OLIGO_KEEP <- c("OPC", "MOL")     # matched against the part after the prefix

jac_l <- jac_l %>%
  filter(lineage != "Oligodendrocyte" |
           str_replace(sig, "^[^_]+_", "") %in% OLIGO_KEEP) %>%
  droplevels()

JMAX <- max(jac_l$jaccard, na.rm = TRUE)

make_hm <- function(d, ttl) {
  ggplot(d, aes(factor(sig), fct_rev(ba), fill = jaccard)) +
    geom_tile(colour = "grey90", linewidth = 0.1) +
    geom_text(aes(label = paste0(sprintf("%.2f", jaccard), "\n(", shared, ")"),
                  colour = jaccard > 0.5 * JMAX),
              size = BASE / .pt, lineheight = 0.85, show.legend = FALSE) +
    scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "grey20")) +
    scale_fill_gradient(low = "white", high = "darkgreen",
                        name = "Jaccard", limits = c(0, JMAX)) +
    scale_x_discrete(labels = function(x) str_replace(x, "^[^_]+_", ""),
                     expand = c(0, 0)) +
    scale_y_discrete(labels = function(x) str_replace(x, "^[^_]+_", ""),
                     expand = c(0, 0)) +
    labs(x = NULL, y = NULL, title = ttl) +
    theme_cowplot(font_size = BASE, rel_small = 1) +
    theme(axis.text.x     = element_text(angle = 45, hjust = 1, colour = "black"),
          axis.text.y     = element_text(colour = "black"),
          axis.ticks      = element_blank(),
          axis.line       = element_blank(),
          legend.key.size = unit(0.25, "cm"))
}

hm_list <- jac_l %>%
  group_split(lineage) %>%
  set_names(map_chr(., ~ as.character(.x$lineage[1]))) %>%
  imap(~ make_hm(.x, .y))

fig_jac <- wrap_plots(hm_list, nrow = 1) +
  plot_layout(guides = "collect", widths = map_dbl(hm_list, ~ n_distinct(.x$data$sig)))

ggsave("BMCA/plots/fig5/dge_jaccard_by_lineage.pdf", fig_jac,
       width = 180, height = 40, units = "mm", dpi =  300)
# -----------------------------------------------------------------------------

df <- read_csv("BMCA/csv/Clinical information for paired samples.csv") %>%
  dplyr::select(Met, Model)
df$patient <- paste0("lin_sn_", df$patient)
intersect(df$patient, metadata$patient)
metadata <- metadata %>% left_join(df, by = "patient")
metadata_p <- metadata %>% dplyr::select(patient, primary_site) %>% distinct() %>%
  filter(patient %in% community_fraction_tib$patient)
community_fraction_tib <- community_fraction_tib %>% left_join(df, by = "patient")
community_fraction_tib <- community_fraction_tib %>% left_join(metadata_p, by = "patient")


plot_df <- community_fraction_tib %>%
  filter(!is.na(`Dissamination Model`), `Dissamination Model` != "NA") %>%
  filter(primary_site != "colorectal")

ggplot(plot_df, aes(x = `Dissamination Model`, y = Brain_adap)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, colour = "grey40") +
  geom_jitter(aes(colour = primary_site),
              width = 0.15, height = 0, size = 2.2, alpha = 0.85) +
  stat_compare_means(method = "t.test", aes(label = paste0("p = ", after_stat(p.format)))) +
  labs(x = "Dissemination model",
       y = "Brain adaptation score",
       colour = "Primary site") +
  theme_classic(base_size = 12)

p <- ggplot(plot_df, aes(x = `Dissamination Model`, y = Brain_adap)) +
  stat_summary(fun = median, geom = "crossbar", width = 0.5,
               colour = "grey40", fatten = 1.5) +
  geom_jitter(aes(colour = primary_site),
              width = 0.12, height = 0, size = 3, alpha = 0.9) +
  stat_compare_means(method = "t.test",
                     aes(label = paste0("p = ", after_stat(p.format))),
                     label.x = 1.4) +          # <- move right (1 = Early, 2 = Late)
  labs(x = "Dissemination model", y = "Brain adaptaion score",
       colour = "Primary site") +
  theme_classic(base_size = 12)
p
ggsave("BMCA/plots/fig5/diss_model_brain_remodeling.pdf", p, width = 5, height = 5, dpi = 600)


plot_df <- plot_df %>%
  mutate(`Brain Met Remodeling Level` = case_when(
    Brain_adap > 0.5  ~ "High",
    Brain_adap < 0.15 ~ "Low",
    TRUE              ~ "Medium"   # the 0.15–0.5 range; use NA_character_ to leave blank
  ))


# ---- Fisher's exact test on the 2x3 table -------------------------
tab <- table(plot_df$`Dissamination Model`, plot_df$`Brain Met Remodeling Level`)
tab   # inspect the contingency table first

ft <- fisher.test(tab)
ft$p.value

# ---- Plot with the p-value in the subtitle ------------------------
ggplot(plot_df, aes(x = `Dissamination Model`, fill = `Brain Met Remodeling Level`)) +
  geom_bar(position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = c("Low"    = "#4575b4",
                               "Medium" = "#fee090",
                               "High"   = "#d73027")) +
  theme_minimal(base_size = 12) +
  labs(x = "Dissemination Model",
       y = "Proportion of patients",
       fill = "Brain Met Remodeling",
       subtitle = paste0("Fisher's exact test, p = ", signif(ft$p.value, 3)))

