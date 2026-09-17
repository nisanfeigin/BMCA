# =============================================================================
# Figure 1 - Stromal compartment: fibroblast, pericyte and endothelial analyses
# =============================================================================
#
# Pseudobulk differential expression between intracranial and extracranial samples
# for the stromal compartment, fibroblast origin/contamination assessment against
# published fibroblast references, endothelial heatmaps, pathway enrichment, and
# the myeloid-fibroblast TGFB2 correlation analyses.
#
# Part of: A brain metastasis cell atlas reveals multicellular adaptation to
#          the neural microenvironment (Feigin et al.)
#
# Run from the project root (see README.md). All paths below are relative to it.
#
# Inputs:
#   BMCA/csv/all_MPS.xlsx
#   BMCA/csv/barnett_endo.csv
#   BMCA/csv/canonical_markers2.csv
#   BMCA/csv/cords_fib.csv
#   BMCA/csv/endo_joyce.csv
#   BMCA/csv/gao_fib.csv
#   BMCA/csv/luo_fib.csv
#   BMCA/data/centered_cpm_lin_sn.RDS
#   BMCA/data/expmat_cpm_lin_sn.RDS
#   BMCA/data/metadata_all_studies.rds
#   BMCA/data/score_tib_Mesenchimal.RDS
#   BMCA/data/UMI_pb_<cell_type>.RDS   # path built at run time
#   BMCA/data/UMI_pb_Fibroblast.RDS
#   BMCA/data/UMI_pb_Myeloid.RDS
#
# Outputs:
#   BMCA/data/  (7 files)
#   BMCA/plots/  (20 files)
#   BMCA/plots/fig1/  (13 files)
#
# Original filename: 22_perycites_scores.R
# =============================================================================

# ---- Packages ---------------------------------------------------------------
pkgs <- c(
  "biomaRt", "cowplot", "dendextend", "doParallel", "dplyr", "enrichR",
  "foreach", "ggdendro", "ggplot2", "ggpubr", "ggrepel", "Matrix", "NMF",
  "openxlsx", "parallel", "patchwork", "pheatmap", "purrr",
  "RColorBrewer", "readxl", "reshape2", "rlang", "scales", "scalop",
  "stringr", "tibble", "tidyr", "tidytext", "tidyverse", "viridis"
)
invisible(lapply(pkgs, library, character.only = TRUE))

# ---- Project paths ----------------------------------------------------------
source("R/config.R")

# =============================================================================
# Path to
metadata <- as_tibble(readRDS("BMCA/data/metadata_all_studies.rds"))

IC_cells <- metadata %>%
  filter(location == "Intracranial") %>%
  filter(cell_type == "Fibroblast" | cell_type == "Pericyte") %>%
  pull(CellID)

EC_cells <- metadata %>%
  filter(location == "Extracranial") %>%
  filter(cell_type == "Fibroblast" | cell_type == "Pericyte") %>%
  pull(CellID)

Mesenchimal_cells <- metadata %>%
  filter(cell_type == "Fibroblast" | cell_type == "Pericyte") %>%
  pull(CellID)

score_tib <- readRDS(paste0("BMCA/data/score_tib_Mesenchimal.RDS"))
fib <-  score_tib %>%
  filter(annotation %in% c("sn_CAF1_Fibroblasts", "sn_CAF4_Fibroblasts")) %>%
  pull(CellID)

peri <-  score_tib %>%
  filter(annotation %in% c("sn_Pericyte-like_Fibroblasts")) %>%
  pull(CellID)

score_tib_mes <-  score_tib %>%
  filter(annotation %in% c("sn_UA_2", "sn_Basement_Membrane", "sc_CC",
                           "sc_Stress", "sc_Hypoxia", "sc_Smooth_muscle", "sn_UA_1")) %>%
  filter(CellID %in% Mesenchimal_cells)

score_tib_mes <- score_tib_mes %>%
  mutate(
    Pericyte_avg   = rowMeans(cbind(sc_Pericyte, sn_Pericyte), na.rm = TRUE),
    Fibroblast_avg = rowMeans(cbind(sc_Fibroblasts, sn_Fibroblasts), na.rm = TRUE)
  ) %>%
  left_join(metadata %>% dplyr::select(CellID, location, seq_tech), by = "CellID")


score_tib_mes_sn_ic <-score_tib_mes %>%
  filter(seq_tech == "sn" & location == "Intracranial")

score_tib_mes_sn_ec <-score_tib_mes %>%
  filter(seq_tech == "sn" & location == "Extracranial")

score_tib_mes_sc_ic <-score_tib_mes %>%
  filter(seq_tech == "sc" & location == "Intracranial")

score_tib_mes_sc_ec <-score_tib_mes %>%
  filter(seq_tech == "sc" & location == "Extracranial")

ggplot(score_tib_mes_sn_ic,
       aes(x = Pericyte_avg, y = Fibroblast_avg)) +
  stat_density_2d(
    na.rm = TRUE,
    contour_var = "ndensity",
    bins = 10,
    color = "Black"
  ) +
  scale_fill_viridis_d() +
  theme_classic() +
  labs(
    x = "Average Pericyte score",
    y = "Average Fibroblast score",
    title = "Pericyte vs Fibroblast scores sn_Intracranial"
  ) +
  geom_abline(slope = 2, intercept = -1, color = "darkred", linetype = "dashed", linewidth = 1)


ggplot(score_tib_mes_sn_ec,
       aes(x = Pericyte_avg, y = Fibroblast_avg)) +
  stat_density_2d(
    na.rm = TRUE,
    contour_var = "ndensity",
    bins = 10,
    color = "Black"
  ) +
  scale_fill_viridis_d() +
  theme_classic() +
  labs(
    x = "Average Pericyte score",
    y = "Average Fibroblast score",
    title = "Pericyte vs Fibroblast scores sn_Extracranial"
  )

ggplot(score_tib_mes_sc_ic,
       aes(x = Pericyte_avg, y = Fibroblast_avg)) +
  stat_density_2d(
    na.rm = TRUE,
    contour_var = "ndensity",
    bins = 10,
    color = "Black"
  ) +
  scale_fill_viridis_d() +
  theme_classic() +
  labs(
    x = "Average Pericyte score",
    y = "Average Fibroblast score",
    title = "Pericyte vs Fibroblast scores sc_Intracranial"
  ) +
  geom_abline(slope = 2, intercept = 1.6, color = "darkred", linetype = "dashed", linewidth = 1)


ggplot(score_tib_mes_sc_ec,
       aes(x = Pericyte_avg, y = Fibroblast_avg)) +
  stat_density_2d(
    na.rm = TRUE,
    contour_var = "ndensity",
    bins = 10,
    color = "Black"
  ) +
  scale_fill_viridis_d() +
  theme_classic() +
  labs(
    x = "Average Pericyte score",
    y = "Average Fibroblast score",
    title = "Pericyte vs Fibroblast scores sc_Extracranial"
  ) +
  geom_abline(slope = 1, intercept = -1.5, color = "darkred", linetype = "dashed", linewidth = 1)


score_tib <-  score_tib %>%
  filter(annotation %in% c("sc_Fibroblasts", "sn_Fibroblasts", "sc_Pericyte", "sn_Pericyte")) %>%
  filter(CellID %in% IC_cells)


score_tib2 <- score_tib %>%
  mutate(
    Pericyte_avg   = rowMeans(cbind(sc_Pericyte, sn_Pericyte), na.rm = TRUE),
    Fibroblast_avg = rowMeans(cbind(sc_Fibroblasts, sn_Fibroblasts), na.rm = TRUE)
  )

ggplot(score_tib2,
       aes(x = Pericyte_avg, y = Fibroblast_avg)) +
  geom_point(alpha = 0.3, size = 0.6) +
  theme_minimal(base_size = 12) +
  labs(
    x = "Average Pericyte score (sc + sn)",
    y = "Average Fibroblast score (sc + sn)",
    title = "Pericyte vs Fibroblast scores"
  )

ggplot(score_tib,
       aes(x = sn_Pericyte, y = sn_Fibroblasts)) +
  geom_point(alpha = 0.3, size = 0.6) +
  theme_minimal(base_size = 12) +
  labs(
    x = "Average Pericyte score (sn)",
    y = "Average Fibroblast score (sn)",
    title = "Pericyte vs Fibroblast scores -sn"
  )

ggplot(score_tib,
       aes(x = sc_Pericyte, y = sc_Fibroblasts)) +
  geom_point(alpha = 0.3, size = 0.6) +
  theme_minimal(base_size = 12) +
  labs(
    x = "Average Pericyte score (sc)",
    y = "Average Fibroblast score (sc)",
    title = "Pericyte vs Fibroblast scores - sc"
  )

#7A5633"   # dark brown

ggplot(score_tib2,
       aes(x = Pericyte_avg, y = Fibroblast_avg)) +
  stat_density_2d(
    na.rm = TRUE,
    contour_var = "ndensity",
    bins = 10,
    color = "Black"
  ) +
  scale_fill_viridis_d() +
  theme_classic() +
  labs(
    x = "Average Pericyte score (sc + sn)",
    y = "Average Fibroblast score (sc + sn)",
    title = "Pericyte vs Fibroblast scores"
  )

ggplot(score_tib,
       aes(x = sn_Pericyte, y = sn_Fibroblasts)) +
  stat_density_2d(
    na.rm = TRUE,
    contour_var = "ndensity",
    bins = 10,
    color = "Black"
  ) +
  scale_fill_viridis_d() +
  theme_classic()   +
  labs(
    x = "Average Pericyte score (sn)",
    y = "Average Fibroblast score (sn)",
    title = "Pericyte vs Fibroblast scores -sn"
  )

ggplot(score_tib,
       aes(x = sc_Pericyte, y = sc_Fibroblasts)) +
  stat_density_2d(
    na.rm = TRUE,
    contour_var = "ndensity",
    bins = 10,
    color = "Black"
  ) +
  scale_fill_viridis_d() +
  theme_classic()   +
  labs(
    x = "Average Pericyte score (sc)",
    y = "Average Fibroblast score (sc)",
    title = "Pericyte vs Fibroblast scores -sc"
  )


score_tib_mes_annotated <- score_tib_mes %>%
  mutate(
    annotation_new = case_when(

      ## sn + Intracranial
      seq_tech == "sn" & location == "Intracranial" ~
        if_else(
          Fibroblast_avg > (2 * Pericyte_avg - 1),
          "Fibroblast",
          "Pericyte"
        ),

      ## sc + Intracranial
      seq_tech == "sc" & location == "Intracranial" ~
        if_else(
          Fibroblast_avg > (2 * Pericyte_avg + 1.6),
          "Fibroblast",
          "Pericyte"
        ),

      ## sc + Extracranial
      seq_tech == "sc" & location == "Extracranial" ~
        if_else(
          Fibroblast_avg > (1 * Pericyte_avg - 1.5),
          "Fibroblast",
          "Pericyte"
        ),

      ## sn + Extracranial → all Fibroblast
      seq_tech == "sn" & location == "Extracranial" ~
        "Fibroblast",

      ## fallback (if any unexpected rows exist)
      TRUE ~ NA_character_
    )
  )


score_tib_mes_annotated <- score_tib_mes %>%
  mutate(
    annotation_new = case_when(

      ## sn + Intracranial
      seq_tech == "sn" & location == "Intracranial" ~
        if_else(
          Fibroblast_avg > (2 * Pericyte_avg - 1),
          "Fibroblast",
          "Pericyte"
        ),

      ## sc + Intracranial
      seq_tech == "sc" & location == "Intracranial" ~
        if_else(
          Fibroblast_avg > (2 * Pericyte_avg + 1.6),
          "Fibroblast",
          "Pericyte"
        ),

      ## sc + Extracranial
      seq_tech == "sc" & location == "Extracranial" ~
        if_else(
          Fibroblast_avg > (1 * Pericyte_avg - 1.5),
          "Fibroblast",
          "Pericyte"
        ),

      ## sn + Extracranial → all Fibroblast
      seq_tech == "sn" & location == "Extracranial" ~
        "Fibroblast",

      ## fallback (if any unexpected rows exist)
      TRUE ~ NA_character_
    )
  )

score_tib_mes_annotated <- score_tib_mes_annotated %>% dplyr::select(CellID, annotation_new)
metadata <- metadata %>% left_join(score_tib_mes_annotated, by = "CellID")
saveRDS(metadata, "BMCA/data/metadata_all_studies_bu.rds")
# -----------------------------------------------------------------------------
metadata <- as_tibble(readRDS("BMCA/data/metadata_all_studies.rds"))

library(tidyverse); library(cowplot)

BASE <- 12          # presentation, not print

comp_df <- metadata %>%
  filter(!doublet_flag,
         !is.na(cell_type), !is.na(location),
         location %in% c("Intracranial", "Extracranial")) %>%
  count(location, patient, cell_type, name = "n_cells") %>%
  group_by(location, patient) %>%
  mutate(frac = n_cells / sum(n_cells)) %>%          # within patient -> sums to 1
  group_by(location, cell_type) %>%
  summarise(frac = mean(frac), n_pat = n(), .groups = "drop") %>%
  group_by(location) %>%
  mutate(frac = frac / sum(frac)) %>%                # renormalise: absent types aren't NA
  ungroup() %>%
  mutate(cell_type = factor(cell_type, levels = major_cts),
         location  = factor(location, levels = c("Extracranial", "Intracranial")))

n_lab <- metadata %>%
  filter(!doublet_flag, location %in% c("Intracranial", "Extracranial")) %>%
  distinct(location, patient) %>%
  count(location, name = "n") %>%
  mutate(lab = paste0(location, "\n(n = ", n, ")"))

p_comp <- ggplot(comp_df, aes(location, frac, fill = cell_type)) +
  geom_col(width = 0.6, colour = "white", linewidth = 0.2) +
  scale_fill_manual(values = CT_COLORS, drop = FALSE, name = "Cell type",
                    labels = function(x) str_replace_all(x, "_", "/")) +
  scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
  scale_x_discrete(labels = setNames(n_lab$lab, n_lab$location)) +
  labs(x = NULL, y = "Mean fraction of cells per patient") +
  theme_cowplot(font_size = BASE) +
  theme(axis.text = element_text(colour = "black"),
        legend.key.size = unit(0.5, "cm"))

ggsave("BMCA/plots/celltype_composition_IC_vs_EC.pdf", p_comp,
       width = 120, height = 150, units = "mm")

cell_types <- c("Endothelial",     "Fibroblast",      "Pericyte",
                "B_Plasma", "Myeloid", "CD8", "CD4", "NK",
                "Malignant",       "Oligodendrocyte",
                "Neuron",          "Astrocyte")

metadata_p <- metadata %>%
  filter(study == "lin_sn") %>%
  dplyr::select(patient, primary_site, location, Donor) %>%
  distinct()

list_of_pbs <- list()
for (ct in cell_types) {
  print(ct)
  pb <- readRDS(paste0("BMCA/data/UMI_pb_", ct, ".RDS"))[["lin_sn"]]
  #pb <- pb[, intersect(colnames(pb), intracranial), drop = FALSE]  # <- filter here
  list_of_pbs[[ct]] <- pb
}

fib_pb <- list_of_pbs[["Fibroblast"]]


# intracranial patients only
intracranial <- metadata_p %>%
  dplyr::filter(location == "Intracranial") %>%
  dplyr::pull(patient) %>%
  unique()


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

top_genes_n <- 1000

for (ct in cell_types) {
  print(ct)
  d <- dplyr::filter(plot_df, cell_type == ct)

  # x threshold: value leaving 1000 genes to its right
  x_cut <- sort(d$avg_expr, decreasing = TRUE)[1000]

  # y threshold: value leaving 1000 genes above it
  y_cut <- sort(d$log2fc, decreasing = TRUE)[top_genes_n]

  # top 20 genes by log2FC
  top20 <- d %>% dplyr::slice_max(log2fc, n = 20)

  p <- ggplot(d, aes(avg_expr, log2fc)) +
    geom_point(size = 0.4, alpha = 0.35, colour = "grey40") +
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

genes_of_interest <- c("NTRK2", "NTRK3")
ct <- "Fibroblast"
for (ct in cell_types) {
  d <- dplyr::filter(plot_df, cell_type == ct)

  # x threshold: value leaving 1000 genes to its right
  x_cut <- sort(d$avg_expr, decreasing = TRUE)[1000]

  # y threshold: value leaving 1000 genes above it
  y_cut <- sort(d$log2fc, decreasing = TRUE)[top_genes_n]

  # top 20 genes by log2FC
  top20 <- d %>% dplyr::slice_max(log2fc, n = 20)

  # genes of interest present in this cell type
  goi <- d %>% dplyr::filter(gene %in% genes_of_interest)

  p <- ggplot(d, aes(avg_expr, log2fc)) +
    geom_point(size = 0.4, alpha = 0.35, colour = "grey40") +
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
    # genes of interest on top
    geom_point(data = goi, size = 1.6, colour = "black") +
    ggrepel::geom_text_repel(
      data = goi, aes(label = gene),
      size = 3, fontface = "bold", colour = "black",
      max.overlaps = Inf, segment.size = 0.3,
      box.padding = 0.5, min.segment.length = 0
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

names(top_genes) <- cell_types
saveRDS(top_genes, "BMCA/data/top_genes_dge_by_pb.RDS")

## ============================================================
## Design matrix for the 12 paired comparisons (24 columns)
## ============================================================
attrs    <- c("donor", "primary", "location")
GAP_PAIR <- 0

build_block <- function(focal) {
  others <- setdiff(attrs, focal)
  combos <- expand.grid(o1 = c(0, 1), o2 = c(0, 1))
  map_dfr(seq_len(4), function(k) {
    v <- setNames(rep(NA_real_, 3), attrs)
    v[others[1]] <- combos$o1[k]
    v[others[2]] <- combos$o2[k]
    tibble(block = focal, pair = k, role = c("test", "ref"),
           donor    = c(if_else(focal == "donor",    1, v["donor"]),
                        if_else(focal == "donor",    0, v["donor"])),
           primary  = c(if_else(focal == "primary",  1, v["primary"]),
                        if_else(focal == "primary",  0, v["primary"])),
           location = c(if_else(focal == "location", 1, v["location"]),
                        if_else(focal == "location", 0, v["location"])))
  })
}

design <- map_dfr(attrs, build_block) %>%
  mutate(block    = factor(block, levels = c("donor", "primary", "location")),
         role     = factor(role,  levels = c("test", "ref")),
         col_id   = row_number(),
         x_within = (pair - 1) * (2 + GAP_PAIR) + if_else(role == "test", 1, 2))

## ---- Helper: one block -> one heatmap object ---------------
make_block_hm <- function(focal, title) {
  hm_df <- design %>%
    filter(block == focal) %>%
    dplyr::select(pair, role, col_id, x_within, donor, primary, location) %>%
    pivot_longer(c(donor, primary, location),
                 names_to = "attr", values_to = "on") %>%
    mutate(attr = factor(attr, levels = c("location", "primary", "donor")))

  ggplot(hm_df, aes(x_within, attr)) +
    geom_tile(aes(fill = factor(on)), colour = "grey45",
              width = 0.95, height = 0.9) +
    scale_fill_manual(values = c(`0` = "grey85", `1` = "black"), guide = "none") +
    scale_y_discrete(labels = c(donor = "Donor", primary = "Primary",
                                location = "Location")) +
    labs(x = NULL, y = NULL, title = title) +
    theme_minimal_grid() +
    theme(axis.text.x  = element_blank(),
          axis.ticks.x = element_blank(),
          panel.grid   = element_blank(),
          axis.text.y  = element_blank(),
          plot.title   = element_blank())
}

## ---- Three separate objects --------------------------------
p_donor    <- make_block_hm("donor",    "Donor")
p_primary  <- make_block_hm("primary",  "Primary")
p_location <- make_block_hm("location", "Location")

p_donor
p_primary
p_location

## ---- Save individually -------------------------------------
ggsave("fib_design_donor.pdf",    p_donor,    width = 3.2, height = 1.6)
ggsave("fib_design_primary.pdf",  p_primary,  width = 3.2, height = 1.6)
ggsave("fib_design_location.pdf", p_location, width = 3.2, height = 1.6)
# -----------------------------------------------------------------------------

COR_METHOD <- "pearson"
MIN_CPM    <- 1
MIN_FRAC   <- 0.10
USE_FILTER <- TRUE          # FALSE = use the gene list as-is, no expression filter

samples <- colnames(fib_pb)
meta <- metadata_p %>% filter(patient %in% samples) %>% distinct(patient, .keep_all = TRUE)
meta <- meta[match(samples, meta$patient), ]
stopifnot(identical(meta$patient, samples))
meta <- meta %>% mutate(donor = if_else(is.na(Donor), patient, Donor))

cpm    <- as.matrix(fib_pb)
logcpm <- log2(cpm + 1)

genes <- top_genes[["Fibroblast"]]
if (USE_FILTER) {
  expressed <- rowSums(cpm >= MIN_CPM) >= (MIN_FRAC * ncol(cpm))
  genes <- intersect(genes, rownames(cpm)[expressed])
} else {
  genes <- intersect(genes, rownames(cpm))
}
message("genes used: ", length(genes), " / ", length(top_genes[["Fibroblast"]]))

cor_mat <- cor(logcpm[genes, ], method = COR_METHOD)

ann <- meta %>% dplyr::select(patient, primary_site, location, donor)
pair_tbl <- as.data.frame(cor_mat) %>%
  rownames_to_column("s1") %>%
  pivot_longer(-s1, names_to = "s2", values_to = "sim") %>%
  filter(s1 < s2) %>%
  left_join(ann, by = c("s1" = "patient")) %>%
  rename(donor1 = donor, primary1 = primary_site, loc1 = location) %>%
  left_join(ann, by = c("s2" = "patient")) %>%
  rename(donor2 = donor, primary2 = primary_site, loc2 = location) %>%
  mutate(donor    = as.integer(donor1   == donor2),
         primary  = as.integer(primary1 == primary2),
         location = as.integer(loc1      == loc2),
         triple   = paste(donor, primary, location, sep = "-"))

SIM_RANGE <- range(pair_tbl$sim)
message("COR_METHOD = ", COR_METHOD,
        " | sim range: ", paste(round(SIM_RANGE, 3), collapse = " to "))

## ============================================================
## Design matrix for the 12 paired comparisons (24 columns)
## ============================================================
attrs    <- c("donor", "primary", "location")
GAP_PAIR <- 0

build_block <- function(focal) {
  others <- setdiff(attrs, focal)
  combos <- expand.grid(o1 = c(0, 1), o2 = c(0, 1))
  map_dfr(seq_len(4), function(k) {
    v <- setNames(rep(NA_real_, 3), attrs)
    v[others[1]] <- combos$o1[k]
    v[others[2]] <- combos$o2[k]
    tibble(block = focal, pair = k, role = c("test", "ref"),
           donor    = c(if_else(focal == "donor",    1, v["donor"]),
                        if_else(focal == "donor",    0, v["donor"])),
           primary  = c(if_else(focal == "primary",  1, v["primary"]),
                        if_else(focal == "primary",  0, v["primary"])),
           location = c(if_else(focal == "location", 1, v["location"]),
                        if_else(focal == "location", 0, v["location"])))
  })
}

design <- map_dfr(attrs, build_block) %>%
  mutate(block    = factor(block, levels = c("donor", "primary", "location")),
         role     = factor(role,  levels = c("test", "ref")),
         col_id   = row_number(),
         x_within = (pair - 1) * (2 + GAP_PAIR) + if_else(role == "test", 1, 2))

## ============================================================
## COMBINED figure: all 24 columns on one axis
## ============================================================
GAP_BLOCK <- 2.0
block_w   <- 8

design <- design %>%
  mutate(block_i = as.integer(block),
         x_pos   = x_within + (block_i - 1) * (block_w + GAP_BLOCK))

block_mid <- design %>% group_by(block) %>%
  summarise(x = mean(range(x_pos)), .groups = "drop")
sep_x_c <- design %>% group_by(block) %>%
  summarise(hi = max(x_pos), .groups = "drop") %>%
  slice(1:2) %>% mutate(x = hi + GAP_BLOCK / 2) %>% pull(x)

col_map <- design %>%
  mutate(triple = paste(donor, primary, location, sep = "-")) %>%
  dplyr::select(block, pair, role, x_pos, x_within, triple)

dat <- pair_tbl %>% inner_join(col_map, by = "triple")

n_tbl <- col_map %>%
  left_join(count(dat, x_pos, name = "n"), by = "x_pos") %>%
  mutate(n = replace_na(n, 0L))

stat <- dat %>%
  group_by(block, pair) %>%
  summarise(ok = n_distinct(role) == 2,
            p  = if (first(ok)) t.test(sim ~ role)$p.value else NA_real_,
            x  = mean(x_pos), .groups = "drop") %>%
  mutate(label = if_else(is.na(p), "n/a", paste0("p=", signif(p, 2))))

X_LIM_C <- c(0.5, max(design$x_pos) + 0.5)

p_box_c <- ggplot(dat, aes(x_pos, sim, group = x_pos, fill = role)) +
  geom_vline(xintercept = sep_x_c, colour = "grey80", linewidth = 0.4) +
  geom_boxplot(width = 0.9, linewidth = 0.3, alpha = 0.5) +
  geom_jitter(width = 0.18, height = 0, size = 0.5, alpha = 0.2, colour = "darkgreen") +
  geom_text(data = stat, aes(x, SIM_RANGE[2] + 0.05, label = label),
            inherit.aes = FALSE, size = 2.6) +
  geom_text(data = block_mid, aes(x, SIM_RANGE[2] + 0.12, label = str_to_title(block)),
            inherit.aes = FALSE, fontface = "bold", size = 3.4) +
  scale_x_continuous(limits = X_LIM_C, expand = c(0, 0)) +
  scale_fill_manual(values = c(test = "black", ref = "grey80"), guide = "none") +
  coord_cartesian(ylim = c(SIM_RANGE[1] - 0.02, SIM_RANGE[2] + 0.16), clip = "off") +
  labs(x = NULL, y = paste0(str_to_title(COR_METHOD), " correlation")) +
  theme_cowplot() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

hm_df <- design %>%
  dplyr::select(x_pos, donor, primary, location) %>%
  pivot_longer(-x_pos, names_to = "attr", values_to = "on") %>%
  mutate(attr = factor(attr, levels = c("location", "primary", "donor")))

p_hm_c <- ggplot(hm_df, aes(x_pos, attr)) +
  geom_tile(aes(fill = factor(on)), colour = "grey45", width = 1, height = 0.9) +
  geom_text(data = n_tbl, aes(x_pos, 0.4, label = n),
            inherit.aes = FALSE, size = 2, colour = "grey35") +
  scale_fill_manual(name = NULL, values = c(`0` = "grey85", `1` = "black"),
                    labels = c(`0` = "Different", `1` = "Same")) +
  scale_x_continuous(limits = X_LIM_C, expand = c(0, 0)) +
  scale_y_discrete(labels = c(donor = "Donor", primary = "Primary", location = "Location"),
                   expand = expansion(add = c(0.8, 0.5))) +
  labs(x = NULL, y = NULL) +
  theme_minimal_grid() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        panel.grid = element_blank(), axis.text.y = element_text(face = "bold"),
        legend.position = "right")

fig <- p_box_c / p_hm_c + plot_layout(heights = c(3, 1))
fig

## ============================================================
## REDUCED figure: only the biologically valid contrasts
## ============================================================
GAP <- 1

contrasts_def <- tribble(
  ~contrast,       ~role,   ~triple,
  "Co-migration",  "focal", "1-1-0",
  "Co-migration",  "ref",   "0-1-0",
  "Location",      "focal", "0-1-1",
  "Location",      "ref",   "0-1-0",
  "Primary",       "focal", "0-1-1",
  "Primary",       "ref",   "0-0-1"
) %>%
  separate(triple, c("donor", "primary", "location"),
           sep = "-", convert = TRUE, remove = FALSE) %>%
  mutate(contrast = factor(contrast, levels = c("Co-migration", "Location", "Primary")),
         ci       = as.integer(contrast),
         within   = if_else(role == "focal", 1, 2),
         x_pos    = within + (ci - 1) * (2 + GAP))

dat_r <- pair_tbl %>% inner_join(contrasts_def, by = "triple")

sep_x   <- c(3, 6)
X_LIM   <- c(0.5, 8.5)
mid_lab <- contrasts_def %>% group_by(contrast) %>%
  summarise(x = mean(x_pos), .groups = "drop")

stat_r <- dat_r %>%
  group_by(contrast) %>%
  summarise(p     = t.test(sim ~ role)$p.value,
            focal = mean(sim[role == "focal"]),
            ref   = mean(sim[role == "ref"]),
            .groups = "drop") %>%
  mutate(ratio = focal / ref,
         p_txt = if_else(p < 0.001, "p < 0.001", paste0("p = ", signif(p, 2))),
         label = paste0("ratio = ", signif(ratio, 3), "\n", p_txt))

ratio_lab <- stat_r %>% left_join(mid_lab, by = "contrast")

n_r <- contrasts_def %>%
  left_join(count(dat_r, x_pos, name = "n"), by = "x_pos") %>%
  mutate(n = replace_na(n, 0L))

BASE <- 6

p_box <- ggplot(dat_r, aes(x_pos, sim, group = x_pos, fill = role)) +
  geom_vline(xintercept = sep_x, colour = "grey85", linewidth = 0.4) +
  geom_boxplot(width = 0.7, outlier.shape = NA, linewidth = 0.3, alpha = 0.55) +
  geom_jitter(width = 0.16, height = 0, size = 0.6, alpha = 0.3, colour = "black") +
  geom_text(data = ratio_lab, aes(x, SIM_RANGE[2] + 0.09, label = label),
            inherit.aes = FALSE, size = BASE / .pt, colour = "grey20", lineheight = 0.9) +
  geom_text(data = mid_lab, aes(x, SIM_RANGE[2] + 0.17, label = contrast),
            inherit.aes = FALSE, size = BASE / .pt) +
  scale_x_continuous(limits = X_LIM, expand = c(0, 0)) +
  scale_fill_manual(values = c(focal = "#3B6FB6", ref = "grey80"), guide = "none") +
  coord_cartesian(ylim = c(SIM_RANGE[1] - 0.02, SIM_RANGE[2] + 0.21), clip = "off") +
  labs(x = NULL, y = paste0(str_to_title(COR_METHOD), " correlation")) +
  theme_cowplot(font_size = BASE) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

p_hm <- ggplot(hm_r, aes(x_pos, attr)) +
  geom_tile(aes(fill = factor(on)), colour = "grey45", width = 0.9, height = 0.9) +
  scale_fill_manual(name = NULL, values = c(`0` = "grey85", `1` = "black"),
                    labels = c(`0` = "Different", `1` = "Same")) +
  scale_x_continuous(limits = c(0.5, 8.5), expand = c(0, 0)) +
  scale_y_discrete(labels = c(donor = "Donor", primary = "Primary", location = "Location")) +
  labs(x = NULL, y = NULL) +
  theme_minimal_grid(font_size = BASE) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        panel.grid = element_blank(),
        legend.position = "right",
        legend.key.size = unit(0.25, "cm"))

fig_r <- p_box / p_hm + plot_layout(heights = c(3, 1))
fig_r
ggsave("BMCA/plots/fig1/comparisons.pdf", fig_r, width = 100, height = 80,
       units = "mm",
       dpi = 300)


## ============================================================
## GLOBAL CUTOFFS (single source of truth)
## ============================================================
LOG2FC_CO   <- 2
P_CO        <- 0.05
TOP_GENES_N <- 1000
EXPR_FLOOR  <- 3      # avg_expr filter for top_genes

cell_types <- c("Endothelial", "Fibroblast", "Pericyte",
                "B_Plasma", "Myeloid", "CD8", "CD4", "NK",
                "Malignant", "Oligodendrocyte", "Neuron", "Astrocyte")

brain_lineages <- c("Neuron", "Oligodendrocyte", "Astrocyte")

## ---- helper: raw (un-centered) log2-CPM matrix for a cell type ----------
load_logcpm <- function(ct) {
  print(ct)
  pb  <- as.matrix(readRDS(paste0("BMCA/data/UMI_pb_", ct, ".RDS"))[["lin_sn"]])
  lib <- colSums(pb)
  log2(t(t(pb) / lib) * 1e6 + 1)          # genes x samples
}

## ---- sample -> location lookup -----------------------------------------
loc_vec <- metadata_p$location
names(loc_vec) <- metadata_p$patient
intra_set <- metadata_p$patient[metadata_p$location == "Intracranial"]
extra_set <- metadata_p$patient[metadata_p$location == "Extracranial"]

## ============================================================
## 1. Define top_genes per cell type (IC samples only)
##    top by log2FC vs other cell types, then expr floor
## ============================================================
# intracranial-only pseudobulks for the cross-cell-type comparison
list_of_pbs <- list()
for (ct in cell_types) {
  print(ct)
  m  <- as.matrix(readRDS(paste0("BMCA/data/UMI_pb_", ct, ".RDS"))[["lin_sn"]])
  m  <- m[, intersect(colnames(m), intra_set), drop = FALSE]
  list_of_pbs[[ct]] <- m
}

# CPM-normalise + log2 each pseudobulk, per-gene mean -> genes x cell_type
norm_pb <- lapply(list_of_pbs, function(m) {
  lib <- colSums(m)
  log2(t(t(m) / lib) * 1e6 + 1)
})
ct_mean <- sapply(norm_pb, rowMeans)

# per-gene df: avg expr in ct + log2FC vs mean of other cell types
plot_df <- bind_rows(lapply(colnames(ct_mean), function(ct) {
  others <- setdiff(colnames(ct_mean), ct)
  tibble(gene      = rownames(ct_mean),
         cell_type = ct,
         avg_expr  = ct_mean[, ct],
         log2fc    = ct_mean[, ct] - rowMeans(ct_mean[, others, drop = FALSE]))
}))

# top TOP_GENES_N by log2FC, then drop low-expression
top_genes <- lapply(colnames(ct_mean), function(ct) {
  d <- plot_df %>% filter(cell_type == ct) %>% arrange(desc(log2fc))
  d <- head(d, TOP_GENES_N) %>% filter(avg_expr > EXPR_FLOOR)
  d$gene
})
names(top_genes) <- colnames(ct_mean)

## ============================================================
## 2. IC vs EC differential expression -> de_tibbles
##    (exclude brain-resident lineages: no EC samples)
## ============================================================
ct_test <- setdiff(names(top_genes), brain_lineages)

de_tibbles <- list()
for (ct in ct_test) {
  message(ct)

  logmat <- load_logcpm(ct)
  g      <- intersect(top_genes[[ct]], rownames(logmat))
  logmat <- logmat[g, , drop = FALSE]
  centered <- logmat - rowMeans(logmat)          # center for the reported log2fc

  loc        <- loc_vec[colnames(centered)]
  intra_cols <- which(loc == "Intracranial")
  extra_cols <- which(loc == "Extracranial")
  print(length(intra_cols))
  print(length(extra_cols))

  if (length(intra_cols) < 5 || length(extra_cols) < 5) {
    message("  skip ", ct, " (intra=", length(intra_cols),
            ", extra=", length(extra_cols), ")")
    next
  }

  mean_intra <- rowMeans(centered[, intra_cols, drop = FALSE])
  mean_extra <- rowMeans(centered[, extra_cols, drop = FALSE])
  log2fc     <- mean_intra - mean_extra

  pval <- vapply(seq_len(nrow(centered)), function(i)
    suppressWarnings(wilcox.test(centered[i, intra_cols],
                                 centered[i, extra_cols])$p.value),
    numeric(1))
  padj <- p.adjust(pval, method = "BH")

  de_tibbles[[ct]] <- tibble(
    gene       = rownames(centered),
    cell_type  = ct,
    mean_intra = mean_intra,
    mean_extra = mean_extra,
    log2fc     = log2fc,
    p_value    = pval,
    padj       = padj)
}

## =======================0=====================================
## 3. Contamination POTENTIAL (replaces contamination_score)
##    potential(x in ct) = max_L mean_IC log2CPM(x in L)
##                         - mean_EC log2CPM(x in ct)
##    L in {Neuron, Oligodendrocyte, Astrocyte}; UN-CENTERED log2-CPM
## ============================================================
# brain-lineage reference: per-gene mean over IC samples, then max across lineages
brain_mean <- list()
for (L in brain_lineages) {
  m  <- load_logcpm(L)
  ic <- intersect(colnames(m), intra_set)
  brain_mean[[L]] <- rowMeans(m[, ic, drop = FALSE])
}
all_genes <- Reduce(union, lapply(brain_mean, names))
brain_mat <- sapply(brain_lineages, function(L) brain_mean[[L]][all_genes])
rownames(brain_mat) <- all_genes
brain_ref <- apply(brain_mat, 1, max, na.rm = TRUE)
brain_ref[!is.finite(brain_ref)] <- NA           # genes absent from all 3 lineages

# per cell type: EC baseline (raw log2-CPM) + potential for every gene
for (ct in names(de_tibbles)) {
  df <- de_tibbles[[ct]]
  m  <- load_logcpm(ct)
  ec <- intersect(colnames(m), extra_set)
  ic <- intersect(colnames(m), intra_set)
  ec_baseline <- rowMeans(m[, ec, drop = FALSE])
  ic_baseline <- rowMeans(m[, ic, drop = FALSE])

  df$brain_ref_log2cpm       <- unname(brain_ref[df$gene])
  df$ec_baseline_log2cpm     <- unname(ec_baseline[df$gene])
  df$ic_baseline_log2cpm     <- unname(ic_baseline[df$gene])
  df$contamination_potential <- df$brain_ref_log2cpm - df$ic_baseline_log2cpm
  de_tibbles[[ct]] <- df
}


LFC_CO  <- 2
PADJ_CO <- 0.05

de_all <- bind_rows(de_tibbles)

n_sig <- de_all %>%
  filter(abs(log2fc) > LFC_CO, padj < PADJ_CO) %>%
  mutate(dir = if_else(log2fc > 0, "IC-up", "EC-up")) %>%
  count(cell_type, dir, name = "n") %>%
  complete(cell_type = unique(de_all$cell_type),
           dir       = c("IC-up", "EC-up"),
           fill      = list(n = 0L)) %>%
  mutate(dir       = factor(dir, levels = c("IC-up", "EC-up")),
         cell_type = fct_reorder(cell_type, n, .fun = sum))

BASE <- 6

p_nsig <- ggplot(n_sig, aes(n, cell_type, fill = dir)) +
  geom_col(width = 0.7, position = position_dodge(width = 0.75)) +
  geom_text(aes(label = n), position = position_dodge(width = 0.75),
            hjust = -0.3, size = BASE / .pt, colour = "grey25") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
  scale_fill_manual(name = NULL,
                    values = c(`IC-up` = "#C0504D", `EC-up` = "#3B6FB6")) +
  labs(x = paste0("Genes with |log2FC| > ", LFC_CO, " and FDR < ", PADJ_CO),
       y = NULL) +
  theme_cowplot(font_size = BASE) +
  theme(legend.position = "top",
        legend.key.size = unit(0.25, "cm"))

p_nsig
ggsave("BMCA/plots/fig1/DGE_bars.pdf", p_nsig, width = 100, height = 100,
       dpi = 300, units = "mm")
## ============================================================
## 4. Volcano plots, significant genes coloured by potential
## ============================================================
# shared potential range so panels are comparable
pot_all   <- unlist(lapply(de_tibbles, function(df)
  df$contamination_potential[!is.na(df$padj) &
                               df$padj < P_CO &
                               abs(df$log2fc) > LOG2FC_CO]))
pot_range <- range(pot_all, na.rm = TRUE)

plots <- list()
ct <- "Fibroblast"
for (ct in names(de_tibbles)) {
  df  <- de_tibbles[[ct]]
  sig <- !is.na(df$padj) & df$padj < P_CO & abs(df$log2fc) > LOG2FC_CO
  d_bg  <- df[!sig, , drop = FALSE]
  d_sig <- df[sig,  , drop = FALSE]

  p <- ggplot(df, aes(log2fc, -log10(padj))) +
    geom_point(data = d_bg, colour = "black", size = 1, alpha = 0.5) +
    geom_point(data = d_sig, aes(colour = contamination_potential),
               size = 2.2, alpha = 0.95) +
    scale_colour_gradient2(
      low = "blue", mid = "darkgrey", high = "red",
      midpoint = 0, limits = pot_range,
      name = "Contamination\npotential\n(log2-CPM)") +
    geom_vline(xintercept = c(-LOG2FC_CO, LOG2FC_CO),
               linetype = "dashed", linewidth = 0.4, colour = "grey50") +
    geom_hline(yintercept = -log10(P_CO),
               linetype = "dashed", linewidth = 0.4, colour = "grey50") +
    geom_text_repel(
      data = d_sig, aes(label = gene),
      size = 2.4, max.overlaps = Inf, colour = "grey20",
      segment.size = 0.2, segment.colour = "grey60",
      box.padding = 0.3, min.segment.length = 0) +
    labs(title = ct,
         x = "log2FC (intracranial - extracranial)",
         y = expression(-log[10]~italic(padj))) +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(face = "bold"))

  plots[[ct]] <- p
  print(p)
}

ggsave("DGE_fibros.pdf", p, width = 200, height = 200,
       dpi = 300, units = "mm")
## ============================================================
## ============================================================
## 5. Per cell type: log2FC (IC-EC) vs contamination potential
##    significant genes coloured IC vs EC, all sig genes labelled
## ============================================================
scatter_plots <- list()
for (ct in names(de_tibbles)) {
  df <- de_tibbles[[ct]]
  d  <- df[!is.na(df$contamination_potential), , drop = FALSE]
  if (!nrow(d)) next

  sig <- !is.na(d$padj) & d$padj < P_CO & abs(d$log2fc) > LOG2FC_CO
  d$class <- dplyr::case_when(
    sig & d$log2fc > 0 ~ "IC-enriched",
    sig & d$log2fc < 0 ~ "EC-enriched",
    TRUE               ~ "ns")
  d$class <- factor(d$class, levels = c("IC-enriched", "EC-enriched", "ns"))

  d_lab <- d[d$class != "ns", , drop = FALSE]      # label every significant gene

  pal <- c(`IC-enriched` = "#C53B33",   # red  = intracranial
           `EC-enriched` = "#2C6FB0",   # blue = extracranial
           `ns`          = "grey75")

  r <- cor(d$log2fc, d$contamination_potential,
           method = "spearman", use = "complete.obs")

  p <- ggplot(d, aes(contamination_potential, log2fc)) +
    geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey70") +
    geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey70") +
    #geom_smooth(method = "loess", se = TRUE, colour = "black",
     #           linewidth = 0.6, fill = "grey80", alpha = 0.4) +
    geom_point(aes(colour = class, size = class != "ns", alpha = class != "ns")) +
    scale_colour_manual(values = pal, name = NULL, drop = FALSE) +
    scale_size_manual(values = c(`FALSE` = 0.8, `TRUE` = 1.8), guide = "none") +
    scale_alpha_manual(values = c(`FALSE` = 0.4, `TRUE` = 0.95), guide = "none") +
    ggrepel::geom_text_repel(
      data = d_lab, aes(label = gene, colour = class),
      size = 2.4, max.overlaps = Inf, show.legend = FALSE,
      segment.size = 0.2, segment.colour = "grey60",
      box.padding = 0.3, min.segment.length = 0) +
    labs(title = paste0(ct),
         x = "Contamination potential\n(max brain lineage - EC baseline, log2-CPM)",
         y = "log2FC (intracranial - extracranial)") +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(face = "bold"))

  scatter_plots[[ct]] <- p
  print(p)
}

de_top <- de_tibbles   # <-- run this while the TOP-GENES version is in memory


## ---- switches ----------------------------------------------
STUDY     <- "lin_sn"
LFC_CO    <- 2
PADJ_CO   <- 0.05
IC_REGEX  <- "brain|intracranial"   # <-- set from count() first
MIN_CPM   <- 1
MIN_FRAC  <- 0.10
CT        <- "Fibroblast"
ANN_TYPES <- c("Neuron", "Oligodendrocyte", "Astrocyte")

ic_samples <- metadata %>%
  filter(study == STUDY, str_detect(tolower(location), IC_REGEX)) %>%
  pull(patient) %>% unique()

## ---- helper: IC-only log2CPM matrix for one cell type -------
get_ic_logcpm <- function(ct) {
  cpm <- as.matrix(readRDS(paste0("BMCA/data/UMI_pb_", ct, ".RDS"))[[STUDY]])
  cpm <- cpm[, intersect(colnames(cpm), ic_samples), drop = FALSE]
  log2(cpm + 1)
}

## ---- gene set + correlation --------------------------------
genes <- de_tibbles[[CT]] %>%
  filter(log2fc > LFC_CO, padj < PADJ_CO) %>%
  pull(gene)

cpm_ct <- as.matrix(readRDS(paste0("BMCA/data/UMI_pb_", CT, ".RDS"))[[STUDY]])
cpm_ct <- cpm_ct[, intersect(colnames(cpm_ct), ic_samples), drop = FALSE]

keep   <- rowSums(cpm_ct >= MIN_CPM) >= (MIN_FRAC * ncol(cpm_ct))
genes  <- intersect(genes, rownames(cpm_ct)[keep])
logcpm <- log2(cpm_ct[genes, , drop = FALSE] + 1)
genes  <- genes[apply(logcpm, 1, var) > 0]
logcpm <- logcpm[genes, , drop = FALSE]

message(CT, ": ", length(genes), " genes x ", ncol(cpm_ct), " IC samples")

cor_mat <- cor(t(logcpm), method = "pearson")
hc      <- hclust(as.dist(1 - cor_mat), method = "average")
ord     <- hc$labels[hc$order]

## ---- annotation: mean log2CPM in brain cell types ----------
ann_df <- map_dfr(ANN_TYPES, function(act) {
  m <- get_ic_logcpm(act)
  tibble(gene      = ord,
         cell_type = act,
         expr      = if_else(ord %in% rownames(m),
                             rowMeans(m[intersect(ord, rownames(m)), , drop = FALSE])[ord],
                             NA_real_),
         n_samp    = ncol(m))
}) %>%
  mutate(gene      = factor(gene, levels = ord),
         cell_type = factor(cell_type, levels = rev(ANN_TYPES)))

walk(ANN_TYPES, ~ message(.x, ": ",
                          sum(is.na(filter(ann_df, cell_type == .x)$expr)), " genes missing"))

## ---- heatmap -----------------------------------------------
cor_df <- as.data.frame(cor_mat) %>%
  rownames_to_column("g1") %>%
  pivot_longer(-g1, names_to = "g2", values_to = "r") %>%
  mutate(g1 = factor(g1, levels = ord),
         g2 = factor(g2, levels = rev(ord)))

lab_size <- if (length(ord) <= 60) 5 else 0

p_gcor <- ggplot(cor_df, aes(g1, g2, fill = r)) +
  geom_tile() +
  scale_fill_gradient2(low = "#3B4CC0", mid = "white", high = "#B40426",
                       midpoint = 0, limits = c(-1, 1), name = "Pearson r") +
  scale_x_discrete(drop = FALSE, expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(x = NULL, y = NULL,
       title = paste0(CT, " — IC-up genes, ", STUDY,
                      ", n = ", ncol(cpm_ct), " IC samples")) +
  theme_cowplot(font_size = 5) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = lab_size),
        axis.ticks.y = element_blank(), axis.line = element_blank())

p_ann <- ggplot(ann_df, aes(gene, cell_type, fill = expr)) +
  geom_tile() +
  scale_fill_viridis_c(option = "magma", na.value = "grey90",
                       name = "Mean\nlog2(CPM+1)") +
  scale_x_discrete(drop = FALSE, expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(x = NULL, y = NULL) +
  theme_cowplot(font_size = 5) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = lab_size),
        axis.text.y = element_text(size = 6),
        axis.ticks  = element_blank(), axis.line = element_blank())

## ---- annotation: contamination potential -------------------
pot_df <- de_tibbles[[CT]] %>%
  dplyr::select(gene, contamination_potential) %>%
  filter(gene %in% ord) %>%
  distinct(gene, .keep_all = TRUE) %>%
  right_join(tibble(gene = ord), by = "gene") %>%
  mutate(gene = factor(gene, levels = ord),
         row  = "Contam. potential")

pot_range <- range(pot_df$contamination_potential, na.rm = TRUE)
message("pot_range: ", paste(round(pot_range, 2), collapse = " to "),
        " | missing: ", sum(is.na(pot_df$contamination_potential)))

p_pot <- ggplot(pot_df, aes(gene, row, fill = contamination_potential)) +
  geom_tile() +
  scale_fill_gradient2(low = "#2166AC", mid = "#E8E8E8", high = "#B2182B",
                       midpoint = 0, limits = pot_range, na.value = "grey90",
                       name = "Contamination\npotential\n(log2-CPM)") +
  scale_x_discrete(drop = FALSE, expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(x = NULL, y = NULL) +
  theme_cowplot(font_size = 5) +
  theme(axis.text.x = element_blank(), axis.ticks = element_blank(),
        axis.line = element_blank(),
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.25, "cm"))

fig_gcor <- p_gcor / p_pot / p_ann +
  plot_layout(heights = c(10, 0.6, 1.2), guides = "collect")

fig_gcor

ggsave("BMCA/plots/fig1/fibroblast_heatmap.pdf", fig_gcor, width = 150, height = 130,
       dpi = 300, units = "mm")

K <- 2

clu <- cutree(hc, k = K)
gene_clusters <- split(names(clu), clu)          # list of 2 gene vectors
names(gene_clusters) <- paste0("C", names(gene_clusters))

sapply(gene_clusters, length)

n_max <- max(lengths(gene_clusters))
top <- gene_clusters %>%
  map(~ c(.x, rep(NA_character_, n_max - length(.x)))) %>%
  as_tibble()

write.csv(top, "BMCA/csv/fib_top_genes_by_group.csv", row.names = FALSE)

gao_fib <- readr::read_csv("BMCA/csv/gao_fib.csv")
library(dplyr); library(tidyr)

wide <- gao_fib %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC), .by_group = TRUE) %>%
  mutate(rank = row_number()) %>%
  ungroup() %>%
  dplyr::select(cluster, rank, gene) %>%
  pivot_wider(names_from = cluster, values_from = gene)  %>%
  dplyr::select(-rank)
gao_fib <- wide[1:30,]
luo_fib <- readr::read_csv("BMCA/csv/luo_fib.csv")
library(dplyr); library(tibble)

# gene-name columns are every 3rd: positions 1, 4, 7, ...
gene_cols <- seq(1, ncol(luo_fib), by = 3)

# cluster IDs (c1, c2, ...) live in the first data row at those columns
clust_names <- as.character(unlist(luo_fib[1, gene_cols]))

luo_wide <- luo_fib[-1, gene_cols]          # drop the sub-header row
colnames(luo_wide) <- clust_names
luo_fib <- as_tibble(luo_wide)[,1:8]

cords_fib <- readr::read_csv("BMCA/csv/cords_fib.csv")
library(dplyr); library(tidyr)

cords_fib <- cords_fib %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC), .by_group = TRUE) %>%
  mutate(rank = row_number()) %>%
  ungroup() %>%
  dplyr::select(cluster, rank, gene) %>%
  pivot_wider(names_from = cluster, values_from = gene) %>%
  dplyr::select(-rank)
cords_fib <- cords_fib[1:30,]

fib_tables <- list(Gao = gao_fib, Luo = luo_fib, Cords = cords_fib)

jaccard  <- function(a, b) length(intersect(a, b)) / length(union(a, b))
top_sets <- lapply(top, function(x) unique(na.omit(x)))

for (nm in names(fib_tables)) {
  fib <- fib_tables[[nm]]
  fib_sets <- lapply(fib, function(x) unique(na.omit(x)))

  # rows = top groups, cols = fib clusters
  jac <- expand_grid(row = names(top_sets), col = names(fib_sets)) %>%
    mutate(jac_val = purrr::map2_dbl(row, col,
                                     ~ jaccard(top_sets[[.x]], fib_sets[[.y]])))

  # order by hierarchical clustering
  m <- jac %>%
    pivot_wider(names_from = col, values_from = jac_val) %>%
    tibble::column_to_rownames("row") %>% as.matrix()
  row_ord <- rownames(m)[hclust(dist(m))$order]
  col_ord <- colnames(m)[hclust(dist(t(m)))$order]

  jac <- jac %>%
    mutate(row = factor(row, levels = row_ord),
           col = factor(col, levels = col_ord))

  p <- ggplot(jac, aes(col, row, fill = jac_val)) +
    geom_tile(colour = "grey90", linewidth = 0.1) +
    scale_fill_gradient(low = "white", high = "darkgreen",
                        name = "Jaccard", limits = c(0, NA)) +
    coord_fixed() +
    labs(x = NULL, y = NULL,
         title = paste0(nm, "et al.")) +
    theme_minimal(base_size = 9) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid = element_blank())
  print(p)

  ggsave(sprintf("BMCA/plots/jaccard_top_vs_%s.pdf", nm), p,
         width = 6, height = 5)
}


# -----------------------------------------------------------------------------

# ---- 1. Databases --------------------------------------------------
databases <- c("GO_Biological_Process_2026", "MSigDB_Hallmark_2020",
               "GO_Cellular_Component_2026", "GO_Molecular_Function_2026",
               "Reactome_2022")

# ---- Fixed color map: one color per database, defined once ---------
db_colors <- setNames(
  brewer.pal(max(3, length(databases)), "Set1")[seq_along(databases)],
  databases
)

out_dir <- "BMCA/plots"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- 2. Wide tibble -> named list of gene vectors ------------------
gene_lists <- top %>%
  as.list() %>%
  map(~ .x[!is.na(.x) & .x != ""])
map_int(gene_lists, length)   # sanity check

# ---- 3. Run enrichr ONCE per cell type -----------------------------
enrich_results <- map(gene_lists, ~ enrichr(.x, databases))

# ---- 4. Manual plot titles -----------------------------------------
# Keys = names(enrich_results); values = display titles. Edit freely.
plot_titles <- c(
  "C1" = "C1",
  "C2"  = "C2"
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

GENE_SIZE      <- 5
GENES_2ROW_CO  <- 4     # more than this -> 2 rows
GENES_3ROW_CO  <- 10    # more than this -> 3 rows
GENES_4ROW_CO  <- 14    # more than this -> 4 rows
WRAP_W         <- 40    # y-axis term label wrap width

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
    theme_minimal(base_size = 6) +
    labs(title = ttl,
         x = "Enrichment Term",
         y = "-log10(FDR Adjusted P-value)",
         fill = "Database")
}

plot_list <- plot_df %>%
  group_split(CellType) %>%
  set_names(map_chr(., ~ as.character(.x$CellType[1]))) %>%
  imap(~ make_plot(.x, .y))

plot_list[[1]]


fig <- plot_list[[1]] + plot_list[[2]]

ggsave("BMCA/plots/fig1/enrichment_C1_C2.pdf", fig,
       width = 200, height = 70, units = "mm")

## =============================================================
## Extracranial (EC) counterpart of the IC enrichment analysis
## Paste directly after the enrichr block (after enrichment_C1_C2.pdf)
## Reuses from earlier in the script:
##   STUDY, LFC_CO, PADJ_CO, IC_REGEX, MIN_CPM, MIN_FRAC, CT,
##   metadata, de_tibbles, databases, db_colors,
##   n_terms, wrap_genes(), make_plot(), plot_titles
## =============================================================

## ---- switches ----------------------------------------------
EC_CLUSTER <- FALSE   # TRUE  -> split EC genes into K_EC co-expression clusters
# FALSE -> treat all EC-up genes as a single set
K_EC       <- 2      # only used when EC_CLUSTER = TRUE

## ---- 1. EC samples (complement of ic_samples) ---------------
ec_samples <- metadata %>%
  filter(study == STUDY,
         !is.na(location),
         !str_detect(tolower(location), IC_REGEX)) %>%
  pull(patient) %>% unique()

message("EC samples: ", length(ec_samples),
        " | IC samples: ", length(ic_samples))

## ---- 2. EC-up gene set, expressed in EC samples -------------
cpm_ec <- as.matrix(readRDS(paste0("BMCA/data/UMI_pb_", CT, ".RDS"))[[STUDY]])
cpm_ec <- cpm_ec[, intersect(colnames(cpm_ec), ec_samples), drop = FALSE]

stopifnot(ncol(cpm_ec) > 0)

genes_ec <- de_tibbles[[CT]] %>%
  filter(log2fc < -LFC_CO, padj < PADJ_CO) %>%   # note the sign flip
  pull(gene)

keep_ec   <- rowSums(cpm_ec >= MIN_CPM) >= (MIN_FRAC * ncol(cpm_ec))
genes_ec  <- intersect(genes_ec, rownames(cpm_ec)[keep_ec])
logcpm_ec <- log2(cpm_ec[genes_ec, , drop = FALSE] + 1)
genes_ec  <- genes_ec[apply(logcpm_ec, 1, var) > 0]
logcpm_ec <- logcpm_ec[genes_ec, , drop = FALSE]

message(CT, " EC: ", length(genes_ec), " genes x ", ncol(cpm_ec), " EC samples")

## ---- 3. Cluster EC genes (same scheme as the IC side) -------
if (EC_CLUSTER && length(genes_ec) >= K_EC * 5) {

  cor_mat_ec <- cor(t(logcpm_ec), method = "pearson")
  hc_ec      <- hclust(as.dist(1 - cor_mat_ec), method = "average")
  ord_ec     <- hc_ec$labels[hc_ec$order]

  clu_ec           <- cutree(hc_ec, k = K_EC)
  gene_clusters_ec <- split(names(clu_ec), clu_ec)
  names(gene_clusters_ec) <- paste0("E", names(gene_clusters_ec))

} else {

  if (EC_CLUSTER) message("Too few EC genes to cluster into ", K_EC,
                          " groups - using a single set instead")
  gene_clusters_ec <- list(E1 = genes_ec)

}

print(sapply(gene_clusters_ec, length))

## ---- 4. Save EC gene groups --------------------------------
n_max_ec <- max(lengths(gene_clusters_ec))
top_ec <- gene_clusters_ec %>%
  map(~ c(.x, rep(NA_character_, n_max_ec - length(.x)))) %>%
  as_tibble()

write.csv(top_ec, "BMCA/csv/fib_EC_top_genes_by_group.csv", row.names = FALSE)

## ---- 5. enrichR on the EC groups ---------------------------
gene_lists_ec <- gene_clusters_ec %>%
  map(~ .x[!is.na(.x) & .x != ""])
print(map_int(gene_lists_ec, length))

enrich_results_ec <- map(gene_lists_ec, ~ enrichr(.x, databases))

## ---- 6. Titles (make_plot() reads plot_titles from globalenv) ----
plot_titles <- c(
  plot_titles,
  setNames(paste0(names(gene_clusters_ec), " (EC-up)"),
           names(gene_clusters_ec))
)

## ---- 7. Same plot_df construction as the IC side -----------
plot_df_ec <- imap_dfr(enrich_results_ec, function(enrich_list, grp) {
  bind_rows(enrich_list, .id = "Database") %>%
    mutate(
      CellType      = grp,
      Database      = factor(Database, levels = databases),
      Count         = as.numeric(gsub("/.*", "", Overlap)),
      neg_log10_fdr = -log10(Adjusted.P.value),
      gene_labels   = str_replace_all(Genes, ";", ", ")
    ) %>%
    arrange(Adjusted.P.value) %>%
    slice_head(n = n_terms)
}) %>%
  mutate(CellType = factor(CellType, levels = names(enrich_results_ec)))

stopifnot(nrow(plot_df_ec) > 0)

## ---- 8. Plots ----------------------------------------------
plot_list_ec <- plot_df_ec %>%
  group_split(CellType) %>%
  set_names(map_chr(., ~ as.character(.x$CellType[1]))) %>%
  imap(~ make_plot(.x, .y))

fig_ec <- wrap_plots(plot_list_ec, nrow = 1)
fig_ec

ggsave(sprintf("BMCA/plots/fig1/enrichment_EC_%s.pdf",
               paste(names(plot_list_ec), collapse = "_")),
       fig_ec,
       width = 200, height = 70, units = "mm", device = cairo_pdf)

## ---- 9. Optional: IC and EC on one sheet -------------------
fig_ic_ec <- wrap_plots(c(plot_list, plot_list_ec), nrow = 2)

ggsave("BMCA/plots/fig1/enrichment_IC_vs_EC.pdf", fig_ic_ec,
       width = 200, height = 140, units = "mm", device = cairo_pdf)

# -----------------------------------------------------------------------------
## =============================================================
## Fibroblast volcano - points coloured by IC cluster (C1/C2),
## EC side black; labels = genes in the top 3 enriched terms of
## each of the three enrichment analyses (C1, C2, EC)
## Run after both enrichr blocks (needs enrich_results + enrich_results_ec)
## =============================================================

library(tidyverse); library(ggrepel); library(cowplot)

## ---- switches ----------------------------------------------
BASE      <- 7        # base font size (pt)
LAB_SIZE  <- 7        # gene label size (pt)
LFC_CO    <- 2
PADJ_CO   <- 0.05
TOP_TERMS <- 3        # top N terms per enrichment analysis
MAX_LAB   <- Inf      # cap labels per set; set e.g. 20 if too crowded

CLUST_COL <- c(C1 = "#1F7D53", C2 = "#6A1E55", EC = "black")

## ---- 1. genes belonging to the top N terms -----------------
## enrichR returns gene symbols upper-cased and ";"-separated,
## so everything is matched on toupper().
genes_from_top_terms <- function(x, n_terms = TOP_TERMS) {
  df <- if (is.data.frame(x[[1]])) bind_rows(x, .id = "Database")
  else bind_rows(map(x, ~ bind_rows(.x, .id = "Database")), .id = "Group")
  df %>%
    arrange(Adjusted.P.value) %>%
    slice_head(n = n_terms) %>%
    pull(Genes) %>%
    str_split(";") %>% unlist() %>% trimws() %>%
    toupper() %>% unique() %>% setdiff("")
}

lab_sets <- list(
  C1 = genes_from_top_terms(enrich_results[["C1"]]),
  C2 = genes_from_top_terms(enrich_results[["C2"]]),
  EC = genes_from_top_terms(enrich_results_ec)   # pools E1/E2 into one "EC"
)

message("genes in top ", TOP_TERMS, " terms - ",
        paste(names(lab_sets), lengths(lab_sets), sep = ": ", collapse = " | "))

lab_map <- imap_dfr(lab_sets, ~ tibble(gene_up = .x, set = .y)) %>%
  distinct(gene_up, .keep_all = TRUE)

## ---- 2. volcano data + point classes -----------------------
clust_map <- tibble(
  gene    = unlist(gene_clusters, use.names = FALSE),
  cluster = rep(names(gene_clusters), lengths(gene_clusters))
)

volc_df <- de_top[["Fibroblast"]] %>%
  left_join(clust_map, by = "gene") %>%
  mutate(
    gene_up   = toupper(gene),
    neglog10p = -log10(pmax(padj, .Machine$double.xmin)),
    sig       = padj < PADJ_CO & abs(log2fc) > LFC_CO,
    pt_class  = case_when(
      sig & log2fc > 0 & cluster == "C1" ~ "C1",
      sig & log2fc > 0 & cluster == "C2" ~ "C2",
      sig & log2fc < 0                   ~ "EC",
      TRUE                               ~ "Other"
    ),
    pt_class = factor(pt_class, levels = c("C1", "C2", "EC", "Other"))
  ) %>%
  arrange(pt_class == "Other" | pt_class == "EC", decreasing = TRUE)

## ---- 3. labels, kept inside the correct quadrant -----------
lab_df <- volc_df %>%
  filter(sig) %>%
  inner_join(lab_map, by = "gene_up") %>%
  filter((set %in% c("C1", "C2") & log2fc > 0) |
           (set == "EC" & log2fc < 0)) %>%
  mutate(lab_col = unname(CLUST_COL[set])) %>%
  group_by(set) %>%
  slice_min(padj, n = MAX_LAB, with_ties = FALSE) %>%
  ungroup()

message("labelled: ",
        paste(names(table(lab_df$set)), table(lab_df$set),
              sep = "=", collapse = " | "))

lab_ic <- filter(lab_df, log2fc > 0)
lab_ec <- filter(lab_df, log2fc < 0)

X_LIM <- max(LFC_CO * 1.1, max(abs(volc_df$log2fc), na.rm = TRUE) * 1.05)

## ---- 4. plot -----------------------------------------------
p_volc <- ggplot(volc_df, aes(log2fc, neglog10p)) +
  geom_vline(xintercept = c(-LFC_CO, LFC_CO), linetype = "dashed",
             colour = "grey60", linewidth = 0.25) +
  geom_hline(yintercept = -log10(PADJ_CO), linetype = "dashed",
             colour = "grey60", linewidth = 0.25) +
  geom_point(aes(fill = pt_class, size = pt_class, alpha = pt_class),
             shape = 21, stroke = 0) +
  geom_text_repel(data = lab_ic,
                  aes(label = gene, colour = lab_col),
                  size = LAB_SIZE / .pt, fontface = "italic",
                  xlim = c(0, NA),
                  min.segment.length = 0, segment.size = 0.2,
                  segment.colour = "grey60",
                  box.padding = 0.15, max.overlaps = Inf,
                  show.legend = FALSE) +
  geom_text_repel(data = lab_ec,
                  aes(label = gene, colour = lab_col),
                  size = LAB_SIZE / .pt, fontface = "italic",
                  xlim = c(NA, 0),
                  min.segment.length = 0, segment.size = 0.2,
                  segment.colour = "grey60",
                  box.padding = 0.15, max.overlaps = Inf,
                  show.legend = FALSE) +
  scale_colour_identity() +
  scale_fill_manual(name = NULL,
                    values = c(CLUST_COL, Other = "grey70"),
                    breaks = c("C1", "C2", "EC")) +
  scale_size_manual(values = c(C1 = 1.5, C2 = 1.5, EC = 1.5, Other = 1),
                    guide = "none") +
  scale_alpha_manual(values = c(C1 = 1, C2 = 1, EC = 1, Other = 0.4),
                     guide = "none") +
  scale_x_continuous(limits = c(-X_LIM, X_LIM)) +
  labs(x = expression(log[2]~"FC (IC vs EC)"),
       y = expression(-log[10]~"FDR")) +
  theme_cowplot(font_size = BASE, rel_small = 1) +
  theme(legend.position = c(0.02, 0.97),
        legend.justification = c(0, 1),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.25, "cm"),
        legend.background = element_blank())

p_volc

ggsave("BMCA/plots/fig1/volcano_fibroblast_enrichlabels.pdf", p_volc,
       width = 100, height = 100, units = "mm", dpi = 300)
gene_cols <- list(EC_enriched = genes_ec, IC_C1 = gene_clusters$C1, IC_C2 = gene_clusters$C2)
write.csv(as_tibble(map(gene_cols, ~ c(.x, rep(NA_character_, max(lengths(gene_cols)) - length(.x))))),
          "BMCA/csv/fib_EC_IC_C1_C2_genes.csv", row.names = FALSE, na = "")
# -----------------------------------------------------------------------------
de_top[["Fibroblast"]]

library(tidyverse); library(ggrepel); library(cowplot)

BASE     <- 10
LFC_CO   <- 2        # where the dashed lines sit
LFC_LAB  <- LFC_CO   # threshold used to pick labels (see note below)
PADJ_CO  <- 0.05

CLUST_COL <- c(C1 = "#1F7D53", C2 = "#6A1E55")

clust_map <- tibble(
  gene    = unlist(gene_clusters, use.names = FALSE),
  cluster = rep(names(gene_clusters), lengths(gene_clusters))
)

volc_df <- de_top[["Fibroblast"]] %>%
  left_join(clust_map, by = "gene") %>%
  mutate(cluster   = factor(replace_na(cluster, "Other"),
                            levels = c("C1", "C2", "Other")),
         neglog10p = -log10(pmax(padj, .Machine$double.xmin))) %>%
  arrange(desc(cluster == "Other"))

# --- everything in the two top quadrants ---
lab_df <- volc_df %>%
  filter(padj < PADJ_CO, abs(log2fc) > LFC_LAB) %>%
  mutate(lab_col = if_else(log2fc > 0,
                           unname(CLUST_COL[as.character(cluster)]) %|% "black",
                           "black"))

X_LIM <- max(LFC_CO * 1.1, max(abs(volc_df$log2fc), na.rm = TRUE) * 1.05)

p_volc <- ggplot(volc_df, aes(log2fc, neglog10p)) +
  geom_vline(xintercept = c(-LFC_CO, LFC_CO), linetype = "dashed",
             colour = "grey60", linewidth = 0.25) +
  geom_hline(yintercept = -log10(PADJ_CO), linetype = "dashed",
             colour = "grey60", linewidth = 0.25) +
  geom_point(aes(fill = cluster, size = cluster, alpha = cluster),
             shape = 21, stroke = 0) +
  geom_text_repel(data = lab_df,
                  aes(label = gene, colour = lab_col),
                  size = BASE / .pt, fontface = "italic",
                  min.segment.length = 0, segment.size = 0.2,
                  segment.colour = "grey60",
                  box.padding = 0.15, max.overlaps = Inf,
                  show.legend = FALSE) +
  scale_colour_identity() +
  scale_fill_manual(name = NULL,
                    values = c(CLUST_COL, Other = "black"),
                    breaks = c("C1", "C2")) +
  scale_size_manual(values = c(C1 = 1.5, C2 = 1.5, Other = 1),
                    guide = "none") +
  scale_alpha_manual(values = c(C1 = 1, C2 = 1, Other = 0.4),
                     guide = "none") +
  scale_x_continuous(limits = c(-X_LIM, X_LIM)) +
  labs(x = expression(log[2]~"FC (IC vs EC)"),
       y = expression(-log[10]~"FDR")) +
  theme_cowplot(font_size = BASE, rel_small = 1) +
  theme(legend.position = c(0.02, 0.97),
        legend.justification = c(0, 1),
        legend.key.size = unit(0.25, "cm"),
        legend.background = element_blank())


ggsave("volcano_fibroblast.pdf", p_volc,
       width = 250, height = 200, units = "mm")
ggsave("BMCA/plots/fig1/volcano_fibroblast.pdf", p_volc,
       width = 188, height = 150, units = "mm")

library(tidyverse); library(cowplot)

BASE   <- 6
STUDY  <- "lin_sn"
TARGET <- "TGFB2"
COR_METHOD <- "pearson"
N_LAB  <- 50
HL   <- "HALLMARK_TGF_BETA_SIGNALING"

## --- hallmark gene set ---------------------------------------
hm_genes <- msigdbr(species = "Homo sapiens", category = "H") %>%
  filter(gs_name == HL) %>%
  pull(gene_symbol) %>%
  unique()

length(hm_genes)     # ~54

mye_all <- log2(as.matrix(readRDS("BMCA/data/UMI_pb_Myeloid.RDS")[[STUDY]]) + 1)
fib_all <- log2(as.matrix(readRDS("BMCA/data/UMI_pb_Fibroblast.RDS")[[STUDY]]) + 1)

loc_sets <- list(Intracranial = intra_set, Extracranial = extra_set)
s <- c("A", "b")
c(s, hm_genes)
rank_one <- function(samps_loc, loc) {
  s <- Reduce(intersect, list(colnames(mye_all), colnames(fib_all), samps_loc))
  message(loc, ": ", length(s), " shared samples")
  if (length(s) < 5) return(tibble())

  mye <- mye_all[, s, drop = FALSE]
  fib <- fib_all[, s, drop = FALSE]
  if (!TARGET %in% rownames(fib) || stats::var(fib[TARGET, ]) == 0) {
    warning(loc, ": ", TARGET, " missing or invariant"); return(tibble())
  }
  g <- union(intersect(top_genes[["Myeloid"]], rownames(mye)),
             intersect(hm_genes, rownames(mye)))
  g <- g[apply(mye[g, , drop = FALSE], 1, stats::var) > 0]

  tibble(gene = g,
         r    = as.numeric(cor(t(mye[g, , drop = FALSE]),
                               fib[TARGET, ], method = COR_METHOD))) %>%
    filter(!is.na(r)) %>%
    arrange(r) %>%
    mutate(rank = row_number(), location = loc, n_samp = length(s))
}

make_rank_plot <- function(d, ttl) {
  lab_df <- bind_rows(slice_head(d, n = N_LAB), slice_tail(d, n = N_LAB))
  ggplot(d, aes(rank, r)) +
    geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.25) +
    geom_point(size = 0.5, alpha = 0.6, colour = "grey30") +
    ggrepel::geom_text_repel(data = lab_df, aes(label = gene),
                             size = BASE / .pt, fontface = "italic",
                             min.segment.length = 0, segment.size = 0.2,
                             segment.colour = "grey60",
                             box.padding = 0.15, max.overlaps = Inf) +
    scale_y_continuous(limits = c(-1, 1)) +
    scale_x_continuous(expand = expansion(mult = 0.03)) +
    labs(x = "Rank (ascending correlation)",
         y = paste0(str_to_title(COR_METHOD), " r vs fibroblast ", TARGET),
         title = paste0(ttl, " (", d$n_samp[1], " samples, ",
                        nrow(d), " genes)")) +
    theme_cowplot(font_size = BASE, rel_small = 1)
}

## --- IC ---
cor_ic <- rank_one(intra_set, "Intracranial")
p_ic   <- make_rank_plot(cor_ic, "Intracranial")
ggsave("BMCA/plots/myeloid_vs_fib_TGFB2_rank_IC.pdf", p_ic,
       width = 89, height = 70, units = "mm")

## --- EC ---
cor_ec <- rank_one(extra_set, "Extracranial")
p_ec   <- make_rank_plot(cor_ec, "Extracranial")
ggsave("BMCA/plots/myeloid_vs_fib_TGFB2_rank_EC.pdf", p_ec,
       width = 89, height = 70, units = "mm")
library(enrichR); library(tidyverse); library(patchwork); library(RColorBrewer)

R_CO    <- 0.5
n_terms <- 5

databases <- c("GO_Biological_Process_2026", "MSigDB_Hallmark_2020",
               "GO_Cellular_Component_2026", "GO_Molecular_Function_2026",
               "Reactome_2022")
db_colors <- setNames(brewer.pal(max(3, length(databases)), "Set1")[seq_along(databases)],
                      databases)

## --- gene sets from the correlation ranking -------------------
make_sets <- function(cor_df) {
  list(Positive = cor_df$gene[cor_df$r >  R_CO],
       Negative = cor_df$gene[cor_df$r < -R_CO]) %>%
    keep(~ length(.x) >= 5)          # enrichr is meaningless below this
}

sets_ic <- make_sets(cor_ic)
sets_ec <- make_sets(cor_ec)

map_int(sets_ic, length)
map_int(sets_ec, length)

## --- run + plot for one location ------------------------------
run_enrich <- function(gene_sets, loc) {
  if (!length(gene_sets)) { message(loc, ": no set passes |r| > ", R_CO); return(NULL) }

  enrich_results <- map(gene_sets, ~ enrichr(.x, databases))

  plot_df <- imap_dfr(enrich_results, function(el, set_name) {
    bind_rows(el, .id = "Database") %>%
      mutate(Database      = factor(Database, levels = databases),
             SetName       = set_name,
             neg_log10_fdr = -log10(Adjusted.P.value),
             gene_labels   = str_replace_all(Genes, ";", ", ")) %>%
      arrange(Adjusted.P.value) %>%
      slice_head(n = n_terms)
  }) %>%
    mutate(SetName = factor(SetName, levels = names(gene_sets)))

  pl <- plot_df %>%
    group_split(SetName) %>%
    set_names(map_chr(., ~ as.character(.x$SetName[1]))) %>%
    imap(~ make_plot(.x, paste0(loc, " — ", .y)))

  fig <- wrap_plots(pl, nrow = 1)
  ggsave(sprintf("BMCA/plots/enrich_myeloid_TGFB2_%s.pdf", loc), fig,
         width = 200, height = 70, units = "mm")
  fig
}

fig_ic <- run_enrich(sets_ic, "Intracranial")
fig_ec <- run_enrich(sets_ec, "Extracranial")

library(tidyverse); library(cowplot)

BASE  <- 6
R_CO  <- 0.5
N_LAB <- 12

cmp_df <- bind_rows(cor_ic, cor_ec) %>%
  dplyr::select(gene, location, r) %>%
  pivot_wider(names_from = location, values_from = r) %>%
  filter(!is.na(Intracranial), !is.na(Extracranial)) %>%
  mutate(delta = Intracranial - Extracranial,
         hl    = case_when(
           Intracranial >  R_CO & Extracranial <  R_CO ~ "IC-specific",
           Extracranial >  R_CO & Intracranial <  R_CO ~ "EC-specific",
           Intracranial >  R_CO & Extracranial >  R_CO ~ "Shared",
           TRUE                                        ~ "Other") %>%
           factor(levels = c("IC-specific", "EC-specific", "Shared", "Other"))) %>%
  arrange(hl == "Other" | hl == "Shared")     # highlighted points on top

count(cmp_df, hl)

lab_df <- cmp_df %>% slice_max(abs(delta), n = N_LAB)

p_cmp <- ggplot(cmp_df, aes(Extracranial, Intracranial)) +
  geom_hline(yintercept = 0, colour = "grey85", linewidth = 0.25) +
  geom_vline(xintercept = 0, colour = "grey85", linewidth = 0.25) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              colour = "grey60", linewidth = 0.25) +
  geom_hline(yintercept = c(-R_CO, R_CO), linetype = "dotted",
             colour = "grey70", linewidth = 0.25) +
  geom_vline(xintercept = c(-R_CO, R_CO), linetype = "dotted",
             colour = "grey70", linewidth = 0.25) +
  geom_point(aes(fill = hl, size = hl, alpha = hl), shape = 21, stroke = 0) +
  ggrepel::geom_text_repel(data = lab_df, aes(label = gene),
                           size = BASE / .pt, fontface = "italic",
                           min.segment.length = 0, segment.size = 0.2,
                           segment.colour = "grey60",
                           box.padding = 0.15, max.overlaps = Inf) +
  scale_fill_manual(name = NULL,
                    values = c(`IC-specific` = "#1F7D53",
                               `EC-specific` = "#6A1E55",
                               Shared        = "#E8871A",
                               Other         = "black"),
                    breaks = c("IC-specific", "EC-specific", "Shared")) +
  scale_size_manual(values = c(`IC-specific` = 2, `EC-specific` = 2,
                               Shared = 2, Other = 2), guide = "none") +
  scale_alpha_manual(values = c(`IC-specific` = 1, `EC-specific` = 1,
                                Shared = 1, Other = 0.35), guide = "none") +
  scale_x_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.5)) +
  scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.5)) +
  coord_fixed() +
  labs(x = paste0("r vs fibroblast TGFB2 — extracranial"),
       y = paste0("r vs fibroblast TGFB2 — intracranial")) +
  theme_cowplot(font_size = BASE, rel_small = 1) +
  theme(legend.position = c(0.02, 0.98),
        legend.justification = c(0, 1),
        legend.key.size = unit(0.25, "cm"),
        legend.background = element_blank())
p_cmp
ggsave("BMCA/plots/myeloid_TGFB2_cor_IC_vs_EC.pdf", p_cmp,
       width = 89, height = 89, units = "mm")

library(tidyverse); library(cowplot); library(msigdbr)

BASE <- 6
HL   <- "HALLMARK_TGF_BETA_SIGNALING"

## --- hallmark gene set ---------------------------------------
hm_genes <- msigdbr(species = "Homo sapiens", category = "H") %>%
  filter(gs_name == HL) %>%
  pull(gene_symbol) %>%
  unique()

length(hm_genes)     # ~54

make_rank_plot <- function(d, ttl) {
  d <- d %>% mutate(in_hm = gene %in% hm_genes) %>% arrange(in_hm)
  lab_df <- filter(d, in_hm)

  ggplot(d, aes(rank, r)) +
    geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.25) +
    geom_point(aes(colour = in_hm, size = in_hm, alpha = in_hm), stroke = 0) +
    ggrepel::geom_text_repel(data = lab_df, aes(label = gene),
                             size = BASE / .pt, fontface = "italic",
                             colour = "#B2182B",
                             min.segment.length = 0, segment.size = 0.2,
                             segment.colour = "grey60",
                             box.padding = 0.15, max.overlaps = Inf) +
    scale_colour_manual(name = NULL,
                        values = c(`FALSE` = "grey75", `TRUE` = "#B2182B"),
                        breaks = TRUE, labels = "TGF-\u03b2 signalling") +
    scale_size_manual(values  = c(`FALSE` = 0.5, `TRUE` = 1.3), guide = "none") +
    scale_alpha_manual(values = c(`FALSE` = 0.4, `TRUE` = 1),   guide = "none") +
    scale_y_continuous(limits = c(-1, 1)) +
    scale_x_continuous(expand = expansion(mult = 0.03)) +
    labs(x = "Rank (ascending correlation)",
         y = "Pearson r vs fibroblast TGFB2",
         title = paste0(ttl, " (", d$n_samp[1], " samples, ",
                        nrow(d), " genes, ", sum(d$in_hm), " hallmark)")) +
    theme_cowplot(font_size = BASE, rel_small = 1) +
    theme(legend.position = c(0.02, 0.98),
          legend.justification = c(0, 1),
          legend.key.size = unit(0.25, "cm"),
          legend.background = element_blank())
}

p_ic <- make_rank_plot(cor_ic, "Intracranial")
ggsave("BMCA/plots/myeloid_vs_fib_TGFB2_rank_IC_hallmark.pdf", p_ic,
       width = 89, height = 70, units = "mm")

p_ec <- make_rank_plot(cor_ec, "Extracranial")
ggsave("BMCA/plots/myeloid_vs_fib_TGFB2_rank_EC_hallmark.pdf", p_ec,
       width = 89, height = 70, units = "mm")


mye_mps <- read_excel("BMCA/csv/all_MPS.xlsx", sheet = "Myeloid")
library(tidyverse); library(cowplot)

BASE <- 6

## --- MP gene sets ---------------------------------------------
mye_sets <- mye_mps %>%
  as.list() %>%
  map(~ .x[!is.na(.x) & .x != ""])

map_int(mye_sets, length)

## --- make sure MP genes are actually in the ranking -----------
## rank_one() only tests top_genes[["Myeloid"]]; widen it once:
ALL_MP <- unique(unlist(mye_sets, use.names = FALSE))

rank_one <- function(samps_loc, loc) {
  s <- Reduce(intersect, list(colnames(mye_all), colnames(fib_all), samps_loc))
  message(loc, ": ", length(s), " shared samples")
  if (length(s) < 5) return(tibble())

  mye <- mye_all[, s, drop = FALSE]
  fib <- fib_all[, s, drop = FALSE]
  if (!TARGET %in% rownames(fib) || stats::var(fib[TARGET, ]) == 0) {
    warning(loc, ": ", TARGET, " missing or invariant"); return(tibble())
  }

  g <- union(intersect(top_genes[["Myeloid"]], rownames(mye)),
             intersect(ALL_MP,                 rownames(mye)))
  g <- g[apply(mye[g, , drop = FALSE], 1, stats::var) > 0]

  tibble(gene = g,
         r    = as.numeric(cor(t(mye[g, , drop = FALSE]),
                               fib[TARGET, ], method = COR_METHOD))) %>%
    filter(!is.na(r)) %>%
    arrange(r) %>%
    mutate(rank = row_number(), location = loc, n_samp = length(s))
}

cor_ic <- rank_one(intra_set, "Intracranial")
cor_ec <- rank_one(extra_set, "Extracranial")

## --- one plot per MP ------------------------------------------
mp_plot <- function(d, mp_genes, mp_name, loc) {
  d <- d %>% mutate(hit = gene %in% mp_genes) %>% arrange(hit)
  ggplot(d, aes(rank, r)) +
    geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.25) +
    geom_point(aes(colour = hit, size = hit, alpha = hit), stroke = 0) +
    ggrepel::geom_text_repel(data = filter(d, hit), aes(label = gene),
                             size = BASE / .pt, fontface = "italic",
                             colour = "#B2182B",
                             min.segment.length = 0, segment.size = 0.2,
                             segment.colour = "grey60",
                             box.padding = 0.15, max.overlaps = Inf) +
    scale_colour_manual(values = c(`FALSE` = "grey78", `TRUE` = "#B2182B"),
                        guide = "none") +
    scale_size_manual(values  = c(`FALSE` = 0.4, `TRUE` = 1.2), guide = "none") +
    scale_alpha_manual(values = c(`FALSE` = 0.35, `TRUE` = 1),  guide = "none") +
    scale_y_continuous(limits = c(-1, 1)) +
    scale_x_continuous(expand = expansion(mult = 0.03)) +
    labs(x = "Rank (ascending correlation)",
         y = paste0("r vs fibroblast ", TARGET),
         title = paste0(str_replace_all(mp_name, "_", " "), " — ", loc,
                        " (", sum(d$hit), "/", length(mp_genes), " genes)")) +
    theme_cowplot(font_size = BASE, rel_small = 1)
}

## --- write one multi-page PDF per location --------------------
save_all <- function(d, loc) {
  if (!nrow(d)) return(invisible(NULL))
  pl <- imap(mye_sets, ~ mp_plot(d, .x, .y, loc))
  pdf(sprintf("BMCA/plots/myeloid_MPs_vs_%s_%s.pdf", TARGET, loc),
      width = 89 / 25.4, height = 70 / 25.4, onefile = TRUE)
  walk(pl, print)
  dev.off()
  pl
}

pl_ic <- save_all(cor_ic, "Intracranial")
pl_ec <- save_all(cor_ec, "Extracranial")

library(tidyverse); library(cowplot)

BASE <- 6

mye_sets <- mye_mps %>%
  as.list() %>%
  map(~ .x[!is.na(.x) & .x != ""])

ALL_MP <- unique(unlist(mye_sets, use.names = FALSE))

## --- coverage check: which MP genes can be scored at all? -----
missing_tbl <- tibble(gene = ALL_MP,
                      in_matrix = gene %in% rownames(mye_all))
count(missing_tbl, in_matrix)
setdiff(ALL_MP, rownames(mye_all))        # genes absent from the pseudobulk

## --- per-location correlation over ALL MP genes ---------------
cor_mp <- function(samps_loc, loc) {
  s <- Reduce(intersect, list(colnames(mye_all), colnames(fib_all), samps_loc))
  message(loc, ": ", length(s), " samples")

  mye <- mye_all[, s, drop = FALSE]
  y   <- fib_all[TARGET, s]

  g  <- intersect(ALL_MP, rownames(mye))
  v  <- apply(mye[g, , drop = FALSE], 1, stats::var)
  message("  ", sum(v == 0), " MP genes invariant here (dropped)")
  g  <- g[v > 0]

  tibble(gene = g,
         r    = as.numeric(cor(t(mye[g, , drop = FALSE]), y, method = COR_METHOD)),
         location = loc, n_samp = length(s))
}

cor_long <- bind_rows(cor_mp(intra_set,  "Intracranial"),
                      cor_mp(extra_set,  "Extracranial"))

## --- mean r per MP --------------------------------------------
mp_long <- imap_dfr(mye_sets, ~ tibble(MP = .y, gene = .x)) %>%
  inner_join(cor_long, by = "gene", relationship = "many-to-many") %>%
  group_by(location, MP) %>%
  summarise(mean_r = mean(r, na.rm = TRUE),
            n_gene = n(),
            n_tot  = length(mye_sets[[first(MP)]]),
            .groups = "drop") %>%
  mutate(location = factor(location, levels = c("Intracranial", "Extracranial")),
         MP       = factor(MP, levels = names(mye_sets)))

## coverage per MP — anything well below n_tot is under-sampled
mp_long %>% filter(location == "Intracranial") %>%
  dplyr::select(MP, n_gene, n_tot) %>% print(n = Inf)

## --- heatmap ---------------------------------------------------
CO <- max(abs(mp_long$mean_r), na.rm = TRUE)

p_mp <- ggplot(mp_long, aes(MP, fct_rev(location), fill = mean_r)) +
  geom_tile(colour = "grey90", linewidth = 0.15) +
  geom_text(aes(label = sprintf("%.2f", mean_r),
                colour = abs(mean_r) > 0.6 * CO),
            size = BASE / .pt, show.legend = FALSE) +
  scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "grey20")) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0, limits = c(-CO, CO),
                       name = paste0("Mean r vs\nfibroblast ", TARGET)) +
  scale_x_discrete(expand = c(0, 0),
                   labels = function(x) str_replace_all(x, "_", " ")) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(x = NULL, y = NULL) +
  theme_cowplot(font_size = BASE, rel_small = 1) +
  theme(axis.text.x     = element_text(angle = 45, hjust = 1, colour = "black"),
        axis.text.y     = element_text(colour = "black"),
        axis.ticks      = element_blank(),
        axis.line       = element_blank(),
        legend.key.size = unit(0.25, "cm"))

ggsave("BMCA/plots/myeloid_MP_meanR_vs_TGFB2.pdf", p_mp,
       width = 180, height = 45, units = "mm")
# -----------------------------------------------------------------------------


## ---- switches ----------------------------------------------
STUDY     <- "lin_sn"
LFC_CO    <- 2
PADJ_CO   <- 0.05
IC_REGEX  <- "brain|intracranial"   # <-- set from count() first
MIN_CPM   <- 1
MIN_FRAC  <- 0.10
CT        <- "Endothelial"
ANN_TYPES <- c("Neuron", "Oligodendrocyte", "Astrocyte")

ic_samples <- metadata %>%
  filter(study == STUDY, str_detect(tolower(location), IC_REGEX)) %>%
  pull(patient) %>% unique()

## ---- helper: IC-only log2CPM matrix for one cell type -------
get_ic_logcpm <- function(ct) {
  cpm <- as.matrix(readRDS(paste0("BMCA/data/UMI_pb_", ct, ".RDS"))[[STUDY]])
  cpm <- cpm[, intersect(colnames(cpm), ic_samples), drop = FALSE]
  log2(cpm + 1)
}

## ---- gene set + correlation --------------------------------
genes <- de_tibbles[[CT]] %>%
  filter(log2fc > LFC_CO, padj < PADJ_CO) %>%
  pull(gene)

cpm_ct <- as.matrix(readRDS(paste0("BMCA/data/UMI_pb_", CT, ".RDS"))[[STUDY]])
cpm_ct <- cpm_ct[, intersect(colnames(cpm_ct), ic_samples), drop = FALSE]

keep   <- rowSums(cpm_ct >= MIN_CPM) >= (MIN_FRAC * ncol(cpm_ct))
genes  <- intersect(genes, rownames(cpm_ct)[keep])
logcpm <- log2(cpm_ct[genes, , drop = FALSE] + 1)
genes  <- genes[apply(logcpm, 1, var) > 0]
logcpm <- logcpm[genes, , drop = FALSE]

message(CT, ": ", length(genes), " genes x ", ncol(cpm_ct), " IC samples")

cor_mat <- cor(t(logcpm), method = "pearson")
hc      <- hclust(as.dist(1 - cor_mat), method = "average")
ord     <- hc$labels[hc$order]

## ---- annotation: mean log2CPM in brain cell types ----------
ann_df <- map_dfr(ANN_TYPES, function(act) {
  m <- get_ic_logcpm(act)
  tibble(gene      = ord,
         cell_type = act,
         expr      = if_else(ord %in% rownames(m),
                             rowMeans(m[intersect(ord, rownames(m)), , drop = FALSE])[ord],
                             NA_real_),
         n_samp    = ncol(m))
}) %>%
  mutate(gene      = factor(gene, levels = ord),
         cell_type = factor(cell_type, levels = rev(ANN_TYPES)))

walk(ANN_TYPES, ~ message(.x, ": ",
                          sum(is.na(filter(ann_df, cell_type == .x)$expr)), " genes missing"))

## ---- heatmap -----------------------------------------------
cor_df <- as.data.frame(cor_mat) %>%
  rownames_to_column("g1") %>%
  pivot_longer(-g1, names_to = "g2", values_to = "r") %>%
  mutate(g1 = factor(g1, levels = ord),
         g2 = factor(g2, levels = rev(ord)))

lab_size <- if (length(ord) <= 60) 5 else 0

p_gcor <- ggplot(cor_df, aes(g1, g2, fill = r)) +
  geom_tile() +
  scale_fill_gradient2(low = "#3B4CC0", mid = "white", high = "#B40426",
                       midpoint = 0, limits = c(-1, 1), name = "Pearson r") +
  scale_x_discrete(drop = FALSE, expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(x = NULL, y = NULL,
       title = paste0(CT, " — IC-up genes, ", STUDY,
                      ", n = ", ncol(cpm_ct), " IC samples")) +
  theme_cowplot(font_size = 5) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = lab_size),
        axis.ticks.y = element_blank(), axis.line = element_blank())

p_ann <- ggplot(ann_df, aes(gene, cell_type, fill = expr)) +
  geom_tile() +
  scale_fill_viridis_c(option = "magma", na.value = "grey90",
                       name = "Mean\nlog2(CPM+1)") +
  scale_x_discrete(drop = FALSE, expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(x = NULL, y = NULL) +
  theme_cowplot(font_size = 5) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = lab_size),
        axis.text.y = element_text(size = 6),
        axis.ticks  = element_blank(), axis.line = element_blank())

## ---- annotation: contamination potential -------------------
pot_df <- de_tibbles[[CT]] %>%
  dplyr::select(gene, contamination_potential) %>%
  filter(gene %in% ord) %>%
  distinct(gene, .keep_all = TRUE) %>%
  right_join(tibble(gene = ord), by = "gene") %>%
  mutate(gene = factor(gene, levels = ord),
         row  = "Contam. potential")

pot_range <- range(pot_df$contamination_potential, na.rm = TRUE)
message("pot_range: ", paste(round(pot_range, 2), collapse = " to "),
        " | missing: ", sum(is.na(pot_df$contamination_potential)))

p_pot <- ggplot(pot_df, aes(gene, row, fill = contamination_potential)) +
  geom_tile() +
  scale_fill_gradient2(low = "#2166AC", mid = "#E8E8E8", high = "#B2182B",
                       midpoint = 0, limits = pot_range, na.value = "grey90",
                       name = "Contamination\npotential\n(log2-CPM)") +
  scale_x_discrete(drop = FALSE, expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(x = NULL, y = NULL) +
  theme_cowplot(font_size = 5) +
  theme(axis.text.x = element_blank(), axis.ticks = element_blank(),
        axis.line = element_blank(),
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.25, "cm"))

fig_gcor <- p_gcor / p_pot / p_ann +
  plot_layout(heights = c(10, 0.6, 1.2), guides = "collect")

fig_gcor

ggsave("BMCA/plots/fig1/endothelial_heatmap.pdf", fig_gcor, width = 150, height = 130,
       dpi = 300, units = "mm")

K <- 4

clu <- cutree(hc, k = K)
gene_clusters <- split(names(clu), clu)          # list of 2 gene vectors
names(gene_clusters) <- paste0("C", names(gene_clusters))

sapply(gene_clusters, length)

n_max <- max(lengths(gene_clusters))
top <- gene_clusters %>%
  map(~ c(.x, rep(NA_character_, n_max - length(.x)))) %>%
  as_tibble()

write.csv(top, "BMCA/csv/endo_top_genes_by_group.csv", row.names = FALSE)


TOP_N  <- 30
FDR_CO <- 0.05

## ---- raw inputs + their column mappings --------------------
endo_raw <- list(
  Joyce   = list(df  = readr::read_csv("BMCA/csv/endo_joyce.csv"),
                 clu = "cluster",   gene = "gene",  lfc = "avg_log2FC",    padj = "p.adj"),
  Barnett = list(df  = readr::read_csv("BMCA/csv/barnett_endo.csv"),
                 clu = "cell_type", gene = "names", lfc = "logfoldchanges", padj = "pvals_adj")
)

## ---- long -> named list of top-N gene vectors --------------
make_sets <- function(df, clu, gene, lfc, padj,
                      top_n = TOP_N, fdr = FDR_CO) {
  df %>%
    dplyr::select(cluster = all_of(clu), gene = all_of(gene),
                  lfc = all_of(lfc), padj = all_of(padj)) %>%
    filter(!is.na(padj), padj < fdr, lfc > 0) %>%
    group_by(cluster) %>%
    slice_max(lfc, n = top_n, with_ties = FALSE) %>%
    ungroup() %>%
    { split(.$gene, .$cluster) } %>%
    lapply(function(x) unique(na.omit(x)))
}

## ---- plotting -----------------------------------------------
top_sets <- gene_clusters

make_jac_plot <- function(ref_sets, nm) {
  jac <- expand_grid(row = names(top_sets), col = names(ref_sets)) %>%
    mutate(jac_val = purrr::map2_dbl(row, col,
                                     ~ jaccard(top_sets[[.x]], ref_sets[[.y]])),
           n_ovl   = purrr::map2_int(row, col,
                                     ~ length(intersect(top_sets[[.x]], ref_sets[[.y]]))))

  m <- jac %>%
    dplyr::select(row, col, jac_val) %>%
    pivot_wider(names_from = col, values_from = jac_val) %>%
    tibble::column_to_rownames("row") %>% as.matrix()

  row_ord <- if (nrow(m) > 2) rownames(m)[hclust(dist(m))$order]      else rownames(m)
  col_ord <- if (ncol(m) > 2) colnames(m)[hclust(dist(t(m)))$order]   else colnames(m)

  jac <- jac %>%
    mutate(row = factor(row, levels = row_ord),
           col = factor(col, levels = col_ord))

  ggplot(jac, aes(col, row, fill = jac_val)) +
    geom_tile(colour = "grey90", linewidth = 0.1) +
    geom_text(aes(label = ifelse(n_ovl > 0, n_ovl, ""),
                  colour = jac_val > 0.5 * max(jac_val)),
              size = 2.5, show.legend = FALSE) +
    scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "grey20")) +
    scale_fill_gradient(low = "white", high = "darkgreen",
                        name = "Jaccard", limits = c(0, NA)) +
    coord_fixed() +
    labs(x = NULL, y = NULL, title = paste0(nm, " et al.")) +
    theme_minimal(base_size = 20) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid  = element_blank())
}

## ---- loop ----------------------------------------------------
endo_sets_all <- list()
jac_plots     <- list()

for (nm in names(endo_raw)) {
  s <- endo_raw[[nm]]
  ref_sets <- make_sets(s$df, s$clu, s$gene, s$lfc, s$padj)

  message(nm, ": ", length(ref_sets), " clusters | genes per cluster: ",
          paste(lengths(ref_sets), collapse = ", "))

  p <- make_jac_plot(ref_sets, nm)
  print(p)

  ggsave(sprintf("BMCA/plots/jaccard_top_vs_%s_endo.pdf", nm), p,
         width = 6, height = 5)

  endo_sets_all[[nm]] <- ref_sets
  jac_plots[[nm]]     <- p
}


# -----------------------------------------------------------------------------

# ---- 1. Databases --------------------------------------------------
databases <- c("GO_Biological_Process_2026", "MSigDB_Hallmark_2020",
               "GO_Cellular_Component_2026", "GO_Molecular_Function_2026",
               "Reactome_2022")

# ---- Fixed color map: one color per database, defined once ---------
db_colors <- setNames(
  brewer.pal(max(3, length(databases)), "Set1")[seq_along(databases)],
  databases
)

out_dir <- "BMCA/plots"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- 2. Wide tibble -> named list of gene vectors ------------------
gene_lists <- top %>%
  dplyr::select(-C3, -C4) %>%
  as.list() %>%
  map(~ .x[!is.na(.x) & .x != ""])
map_int(gene_lists, length)   # sanity check

# ---- 3. Run enrichr ONCE per cell type -----------------------------
enrich_results <- map(gene_lists, ~ enrichr(.x, databases))

# ---- 4. Manual plot titles -----------------------------------------
# Keys = names(enrich_results); values = display titles. Edit freely.
plot_titles <- c(
  "C1" = "C1",
  "C2"  = "C2"
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

GENE_SIZE      <- 5
GENES_2ROW_CO  <- 4     # more than this -> 2 rows
GENES_3ROW_CO  <- 10    # more than this -> 3 rows
GENES_4ROW_CO  <- 14    # more than this -> 4 rows
WRAP_W         <- 40    # y-axis term label wrap width

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
    theme_minimal(base_size = 6) +
    labs(title = ttl,
         x = "Enrichment Term",
         y = "-log10(FDR Adjusted P-value)",
         fill = "Database")
}

plot_list <- plot_df %>%
  group_split(CellType) %>%
  set_names(map_chr(., ~ as.character(.x$CellType[1]))) %>%
  imap(~ make_plot(.x, .y))

plot_list[[1]]


fig <- plot_list[[1]] + plot_list[[2]]

ggsave("BMCA/plots/fig1/enrichment_endo.pdf", fig,
       width = 200, height = 70, units = "mm")

gene_clusters$C3 <- NULL
gene_clusters$C4 <- NULL

## ---- pooled IC matrix, all studies -------------------------
ic_all  <- metadata %>%
  filter(str_detect(tolower(location), IC_REGEX)) %>%
  pull(patient) %>% unique()

pb_list  <- readRDS(paste0("BMCA/data/UMI_pb_", CT, ".RDS"))
common_g <- Reduce(intersect, lapply(pb_list, rownames))

pool <- map(names(pb_list), function(st) {
  cpm <- as.matrix(pb_list[[st]])[common_g, , drop = FALSE]
  cpm[, intersect(colnames(cpm), ic_all), drop = FALSE]
})
names(pool) <- names(pb_list)

## duplicate patient IDs across studies?
dups <- unlist(lapply(pool, colnames))
if (any(duplicated(dups))) message("DUPLICATE IDs: ",
                                   paste(unique(dups[duplicated(dups)]), collapse = ", "))

cpm_all <- do.call(cbind, pool)
mat_log <- log2(cpm_all + 1)


message("pooled: ", ncol(mat_log), " IC samples | ",
        length(common_g), " common genes")
map_int(gene_clusters, ~ length(intersect(.x, rownames(mat_log))))

## ---- score -------------------------------------------------
score_set <- function(mat_log, genes, center = TRUE) {
  g <- intersect(genes, rownames(mat_log))
  m <- mat_log[g, , drop = FALSE]
  if (center) m <- m - rowMeans(m)          # relative expression, Tirosh-style
  colMeans(m)
}

cluster_scores <- imap_dfc(gene_clusters, ~ tibble(!!.y := score_set(mat_log, .x))) %>%
  mutate(patient = colnames(mat_log), .before = 1)

scores_long <- cluster_scores %>%
  pivot_longer(-patient, names_to = "cluster", values_to = "score") %>%
  mutate(cluster = paste0("Endothelial_", cluster)) %>%
  filter(cluster %in% c("Endothelial_C1", "Endothelial_C2"))

saveRDS(scores_long, "BMCA/data/endo_cs.RDS")
library(tidyverse); library(ggrepel); library(cowplot)

BASE     <- 5
LFC_CO   <- 2        # where the dashed lines sit
LFC_LAB  <- LFC_CO   # threshold used to pick labels
PADJ_CO  <- 0.05

C1_COL <- "#9E2A3A"

clust_map <- tibble(
  gene    = unlist(gene_clusters, use.names = FALSE),
  cluster = rep(names(gene_clusters), lengths(gene_clusters))
)

volc_df <- de_top[["Endothelial"]] %>%
  left_join(clust_map, by = "gene") %>%
  mutate(hl        = factor(if_else(cluster %in% "C1", "C1", "Other"),
                            levels = c("C1", "Other")),
         neglog10p = -log10(pmax(padj, .Machine$double.xmin))) %>%
  arrange(desc(hl == "Other"))          # C1 points drawn on top

# --- everything in the two top quadrants ---
lab_df <- volc_df %>%
  filter(padj < PADJ_CO, abs(log2fc) > LFC_LAB) %>%
  mutate(lab_col = if_else(log2fc > 0 & hl == "C1", C1_COL, "black"))

X_LIM <- max(LFC_CO * 1.1, max(abs(volc_df$log2fc), na.rm = TRUE) * 1.05)

p_volc <- ggplot(volc_df, aes(log2fc, neglog10p)) +
  geom_vline(xintercept = c(-LFC_CO, LFC_CO), linetype = "dashed",
             colour = "grey60", linewidth = 0.25) +
  geom_hline(yintercept = -log10(PADJ_CO), linetype = "dashed",
             colour = "grey60", linewidth = 0.25) +
  geom_point(aes(fill = hl, size = hl, alpha = hl),
             shape = 21, stroke = 0) +
  geom_text_repel(data = lab_df,
                  aes(label = gene, colour = lab_col),
                  size = BASE / .pt, fontface = "italic",
                  min.segment.length = 0, segment.size = 0.2,
                  segment.colour = "grey60",
                  box.padding = 0.15, max.overlaps = Inf,
                  show.legend = FALSE) +
  scale_colour_identity() +
  scale_fill_manual(name = NULL,
                    values = c(C1 = C1_COL, Other = "grey75"),
                    breaks = "C1") +
  scale_size_manual(values = c(C1 = 1.3, Other = 0.7), guide = "none") +
  scale_alpha_manual(values = c(C1 = 1, Other = 0.4), guide = "none") +
  scale_x_continuous(limits = c(-X_LIM, X_LIM)) +
  labs(x = expression(log[2]~"FC (IC vs EC)"),
       y = expression(-log[10]~"FDR")) +
  theme_cowplot(font_size = BASE, rel_small = 1) +
  theme(legend.position = c(0.02, 0.97),
        legend.justification = c(0, 1),
        legend.key.size = unit(0.25, "cm"),
        legend.background = element_blank())


ggsave("BMCA/plots/fig1/volcano_endothelial.pdf", p_volc,
       width = 188, height = 150, units = "mm")


# -----------------------------------------------------------------------------


## ============================================================
## GLOBAL CUTOFFS
## ============================================================
LOG2FC_CO  <- 2
P_CO       <- 0.05
EXPR_FLOOR <- 3      # mean log2-CPM floor (replaces the top_genes restriction)

cell_types <- c("Endothelial", "Fibroblast", "Pericyte",
                "B_Plasma", "Myeloid", "CD8", "CD4", "NK",
                "Malignant", "Oligodendrocyte", "Neuron", "Astrocyte")
brain_lineages <- c("Neuron", "Oligodendrocyte", "Astrocyte")

load_logcpm <- function(ct) {
  pb  <- as.matrix(readRDS(paste0("BMCA/data/UMI_pb_", ct, ".RDS"))[["lin_sn"]])
  lib <- colSums(pb)
  log2(t(t(pb) / lib) * 1e6 + 1)
}

loc_vec <- metadata_p$location; names(loc_vec) <- metadata_p$patient
intra_set <- metadata_p$patient[metadata_p$location == "Intracranial"]
extra_set <- metadata_p$patient[metadata_p$location == "Extracranial"]

## ============================================================
## 1. IC vs EC DE over ALL genes (expression floor only)
## ============================================================
ct_test <- setdiff(cell_types, brain_lineages)

## ---- Reference IC/EC significant sets from the TOP-GENES version --------
## Run this on the top-genes de_tibbles BEFORE building the all-genes one.
ref_class <- lapply(de_tibbles, function(df) {
  sig <- !is.na(df$padj) & df$padj < P_CO & abs(df$log2fc) > LOG2FC_CO
  tibble(gene = df$gene[sig],
         ref  = ifelse(df$log2fc[sig] > 0, "IC-enriched", "EC-enriched"))
})
# ref_class[[ct]]$gene = genes that were significant in the top-genes run

de_tibbles <- list()
for (ct in ct_test) {
  message(ct)

  logmat <- load_logcpm(ct)

  # expression floor: keep genes with mean log2-CPM > EXPR_FLOOR (all samples)
  keep_g <- rowMeans(logmat) > EXPR_FLOOR
  logmat <- logmat[keep_g, , drop = FALSE]

  centered <- logmat - rowMeans(logmat)

  loc        <- loc_vec[colnames(centered)]
  intra_cols <- which(loc == "Intracranial")
  extra_cols <- which(loc == "Extracranial")
  if (length(intra_cols) < 2 || length(extra_cols) < 2) {
    message("  skip ", ct); next
  }

  mean_intra <- rowMeans(centered[, intra_cols, drop = FALSE])
  mean_extra <- rowMeans(centered[, extra_cols, drop = FALSE])
  log2fc     <- mean_intra - mean_extra

  pval <- vapply(seq_len(nrow(centered)), function(i)
    suppressWarnings(wilcox.test(centered[i, intra_cols],
                                 centered[i, extra_cols])$p.value),
    numeric(1))
  padj <- p.adjust(pval, method = "BH")

  de_tibbles[[ct]] <- tibble(
    gene = rownames(centered), cell_type = ct,
    mean_intra = mean_intra, mean_extra = mean_extra,
    log2fc = log2fc, p_value = pval, padj = padj)
}

## ============================================================
## 2. Contamination potential (unchanged; now over all kept genes)
## ============================================================
brain_mean <- list()
for (L in brain_lineages) {
  m  <- load_logcpm(L)
  ic <- intersect(colnames(m), intra_set)
  brain_mean[[L]] <- rowMeans(m[, ic, drop = FALSE])
}
all_genes <- Reduce(union, lapply(brain_mean, names))
brain_mat <- sapply(brain_lineages, function(L) brain_mean[[L]][all_genes])
rownames(brain_mat) <- all_genes
brain_ref <- apply(brain_mat, 1, max, na.rm = TRUE)
brain_ref[!is.finite(brain_ref)] <- NA

for (ct in names(de_tibbles)) {
  df <- de_tibbles[[ct]]
  m  <- load_logcpm(ct)
  ec <- intersect(colnames(m), extra_set)
  ic <- intersect(colnames(m), intra_set)
  ec_baseline <- rowMeans(m[, ec, drop = FALSE])
  ic_baseline <- rowMeans(m[, ic, drop = FALSE])

  df$brain_ref_log2cpm       <- unname(brain_ref[df$gene])
  df$ec_baseline_log2cpm     <- unname(ec_baseline[df$gene])
  df$ic_baseline_log2cpm     <- unname(ic_baseline[df$gene])
  df$contamination_potential <- df$brain_ref_log2cpm - df$ic_baseline_log2cpm
  de_tibbles[[ct]] <- df
}

# ============================================================
## 3. Per cell type: contamination potential (x) vs log2FC (y)
## ============================================================
pal <- c(`IC-enriched` = "#C53B33", `EC-enriched` = "#2C6FB0", `ns` = "grey75")

scatter_plots <- list()
for (ct in names(de_tibbles)) {
  df <- de_tibbles[[ct]]
  d  <- df[!is.na(df$contamination_potential), , drop = FALSE]
  if (!nrow(d)) next

  sig <- !is.na(d$padj) & d$padj < P_CO & abs(d$log2fc) > LOG2FC_CO

  # was this gene in the top-genes significant set? (and which direction)
  rc      <- ref_class[[ct]]
  ref_dir <- setNames(rc$ref, rc$gene)[d$gene]      # NA if not in top-genes set

  d$class <- factor(dplyr::case_when(
    sig & !is.na(ref_dir) ~ ref_dir,                # in top-genes set -> IC/EC colour
    sig &  is.na(ref_dir) ~ "New",                  # significant only in all-genes -> black
    TRUE                  ~ "ns"),
    levels = c("IC-enriched", "EC-enriched", "New", "ns"))

  d_lab <- d[d$class != "ns", , drop = FALSE]        # label all significant genes

  pal <- c(`IC-enriched` = "#C53B33",   # red
           `EC-enriched` = "#2C6FB0",   # blue
           `New`         = "black",     # newly significant in all-genes run
           `ns`          = "grey75")

  r <- cor(d$log2fc, d$contamination_potential,
           method = "spearman", use = "complete.obs")

  p <- ggplot(d, aes(contamination_potential, log2fc)) +
    geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey70") +
    geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey70") +
  #  geom_smooth(method = "loess", se = TRUE, colour = "grey40",
   #             linewidth = 0.6, fill = "grey80", alpha = 0.4) +
    geom_point(aes(colour = class, size = class != "ns", alpha = class != "ns")) +
    scale_colour_manual(values = pal, name = NULL, drop = FALSE) +
    scale_size_manual(values = c(`FALSE` = 0.6, `TRUE` = 1.8), guide = "none") +
    scale_alpha_manual(values = c(`FALSE` = 0.3, `TRUE` = 0.95), guide = "none") +
    ggrepel::geom_text_repel(
      data = d_lab, aes(label = gene, colour = class),
      size = 2.4, max.overlaps = Inf, show.legend = FALSE,
      segment.size = 0.2, segment.colour = "grey60",
      box.padding = 0.3, min.segment.length = 0) +
    labs(title = paste0(ct, "  (Spearman r = ", signif(r, 2), ")"),
         x = "Contamination potential\n(max brain lineage - EC baseline, log2-CPM)",
         y = "log2FC (intracranial - extracranial)") +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(face = "bold"))

  scatter_plots[[ct]] <- p
  print(p)
}
## ============================================================
## Ref-significant genes: sharing count as a bar plot
##   one bar per (gene, cell type it was significant in)
##   height = # cell types (full pool) with log2FC > FC_COUNT_CO
##   colour = the cell type the gene was significant in
##   x-axis  = gene (all labelled), sorted by height
## ============================================================
FC_COUNT_CO <- 2
# de_top = the TOP-GENES de_tibbles (with contamination_potential added)
# capture it while that version is in memory:  de_top <- de_tibbles

## ---- 1. reference IC-up calls (top-genes version) ----------------------
ref_sig <- bind_rows(lapply(names(de_top), function(ct) {
  df <- de_top[[ct]]
  s  <- !is.na(df$padj) & df$padj < P_CO & df$log2fc > LOG2FC_CO
  df[s, c("gene", "cell_type", "log2fc", "contamination_potential")]
}))

## ---- 2. n_fc over the FULL gene pool of ALL cell types -----------------
## count, per gene, in how many cell types log2FC(IC-EC) > FC_COUNT_CO
## computed straight from the pseudobulks so it spans every gene
loc_vec   <- metadata_p$location; names(loc_vec) <- metadata_p$patient
brain_lineages <- c("Neuron", "Oligodendrocyte", "Astrocyte")
ct_test   <- setdiff(cell_types, brain_lineages)

load_logcpm <- function(ct) {
  pb  <- as.matrix(readRDS(paste0("BMCA/data/UMI_pb_", ct, ".RDS"))[["lin_sn"]])
  lib <- colSums(pb)
  log2(t(t(pb) / lib) * 1e6 + 1)
}

fc_long <- bind_rows(lapply(ct_test, function(ct) {
  m   <- load_logcpm(ct)
  loc <- loc_vec[colnames(m)]
  ic  <- which(loc == "Intracranial"); ec <- which(loc == "Extracranial")
  if (length(ic) < 1 || length(ec) < 1) return(NULL)
  cen <- m - rowMeans(m)
  tibble(gene      = rownames(m),
         cell_type = ct,                              # <-- add this
         log2fc    = rowMeans(cen[, ic, drop = FALSE]) -
           rowMeans(cen[, ec, drop = FALSE]))
}))

fc_count <- fc_long %>%
  group_by(gene) %>%
  summarise(n_fc = sum(log2fc > FC_COUNT_CO, na.rm = TRUE), .groups = "drop")

## ---- 3. one row per (gene, significant cell type) ----------------------
plot_tbl <- ref_sig %>%
  left_join(fc_count, by = "gene") %>%
  filter(!is.na(contamination_potential), !is.na(n_fc)) %>%
  mutate(cell_type = factor(cell_type),
         bar_id    = paste(gene, cell_type, sep = "__")) %>%
  arrange(desc(n_fc), gene) %>%
  mutate(bar_id = factor(bar_id, levels = bar_id))     # lock x order by height

## ---- 4. bar plot: one bar per (gene, cell type) ------------------------
p_share <- ggplot(plot_tbl, aes(bar_id, n_fc, fill = cell_type)) +
  geom_col(width = 0.8) +
  scale_x_discrete(labels = setNames(plot_tbl$gene, plot_tbl$bar_id)) +
  scale_y_continuous(breaks = seq(0, max(plot_tbl$n_fc), by = 1),
                     expand = expansion(mult = c(0, 0.05))) +
  scale_fill_brewer(palette = "Paired", name = "Significant in") +
  labs(title = paste0("Ref-significant genes: IC-enrichment sharing across cell types"),
       x = NULL,
       y = paste0("# cell types with log2FC > ", FC_COUNT_CO)) +
  theme_bw(base_size = 11) +
  theme(plot.title  = element_text(face = "bold"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10),
        panel.grid.major.x = element_blank())

print(p_share)

## optional: save wide so all gene labels are legible
ggsave("ref_sig_sharing_bars.pdf", p_share,
       width = max(8, nrow(plot_tbl) * 0.12), height = 4.5, device = cairo_pdf)

goi <- c("NTRK2", "NTRK3")

ntrk <- fc_long %>%
  filter(gene %in% goi) %>%
  arrange(gene, desc(log2fc))

# which cell types clear the log2FC > 1 bar
ntrk_above <- ntrk %>% filter(log2fc > FC_COUNT_CO)

cat("Cell types with log2FC >", FC_COUNT_CO, ":\n")
ntrk_above %>%
  group_by(gene) %>%
  summarise(n_fc       = n(),
            cell_types = paste(cell_type, collapse = ", "),
            .groups = "drop") %>%
  print()

# full per-cell-type values (so you see the ones that just missed too)
cat("\nAll cell types, log2FC (IC-EC):\n")
ntrk %>%
  mutate(above = log2fc > FC_COUNT_CO,
         log2fc = round(log2fc, 2)) %>%
  print(n = Inf)


# ---- 1. Databases --------------------------------------------------
databases <- c("GO_Biological_Process_2026", "MSigDB_Hallmark_2020",
               "GO_Cellular_Component_2026",
               "Reactome_2022")

# ---- Fixed color map: one color per database, defined once ---------
db_colors <- setNames(
  brewer.pal(max(3, length(databases)), "Set1")[seq_along(databases)],
  databases
)

# ---- 2. Wide tibble -> named list of gene vectors ------------------
gene_lists <- dge_tbl %>%
  as.list() %>%
  map(~ .x[!is.na(.x) & .x != ""])
map_int(gene_lists, length)   # sanity check

# sizes of the examined gene lists, per cell type / direction
list_sizes <- map_int(gene_lists, length)

# ---- 3. Run enrichr ONCE per cell type / direction -----------------
enrich_results <- map(gene_lists, ~ enrichr(.x, databases))

# ---- 4. Facet names, derived from the *_IC / *_EC column names -----
facet_labels <- setNames(
  names(gene_lists) %>%
    str_replace("_IC$", " (Intracranial)") %>%
    str_replace("_EC$", " (Extracranial)"),
  names(gene_lists)
)

# ---- 5. Build ONE combined data frame, tag direction --------------
n_terms <- 10
plot_df <- imap_dfr(enrich_results, function(enrich_list, cell_type) {
  bind_rows(enrich_list, .id = "Database") %>%
    mutate(
      CellType      = cell_type,
      Database      = factor(Database, levels = databases),
      Count         = as.numeric(gsub("/.*", "", Overlap)),
      pct_examined  = 100 * Count / list_sizes[[cell_type]],
      neg_log10_fdr = -log10(Adjusted.P.value),
      gene_labels = {
        g <- str_split(Genes, ";")[[1]]
        half <- ceiling(length(g) / 2)
        paste(
          paste(g[seq_len(half)], collapse = ", "),
          paste(g[-seq_len(half)], collapse = ", "),
          sep = "\n"
        )
      }    ) %>%
    arrange(Adjusted.P.value) %>%
    slice_head(n = n_terms)
}) %>%
  mutate(
    CellType  = factor(CellType, levels = names(enrich_results)),
    Direction = if_else(str_detect(CellType, "_IC$"), "IC", "EC")
  )

# --- Base replacements for tidytext's reorder_within / scale_x_reordered ---
reorder_within <- function(x, by, within, fun = mean, sep = "___", ...) {
  new_x <- paste(x, within, sep = sep)
  stats::reorder(new_x, by, FUN = fun)
}
scale_x_reordered <- function(..., sep = "___", width = 25) {
  reg <- paste0(sep, ".+$")
  ggplot2::scale_x_discrete(
    labels = function(x) stringr::str_wrap(gsub(reg, "", x), width = width),
    ...
  )
}


# ---- 6. Plot builder, called once per direction -------------------
make_enrich_plot <- function(df, plot_title) {
  ggplot(df,
         aes(x = reorder_within(Term, neg_log10_fdr, CellType),
             y = neg_log10_fdr, fill = Database)) +
    geom_bar(stat = "identity", alpha = 0.8) +
    geom_text(aes(label = gene_labels),
              y = 0.05, hjust = 0, vjust = 0.5, size = 9 / .pt,
              color = "black", lineheight = 0.85) +
    coord_flip(clip = "off") +
    scale_x_reordered(width = 25) +
    scale_fill_manual(values = db_colors, drop = FALSE) +
    facet_wrap(~ CellType, scales = "free_y",
               labeller = as_labeller(facet_labels)) +
    theme_minimal(base_size = 10) +
    theme(
      text          = element_text(size = 10),
      axis.text     = element_text(size = 10),
      axis.title    = element_text(size = 10),
      axis.text.y   = element_text(size = 10, lineheight = 0.85),
      strip.text    = element_text(size = 10),
      legend.text   = element_text(size = 10),
      legend.title  = element_text(size = 10),
      plot.title    = element_text(size = 10, face = "bold")
    ) +
    labs(title = plot_title,
         x = "Enrichment Term",
         y = "-log10(FDR Adjusted P-value)",
         fill = "Database")
}

plot_IC <- make_enrich_plot(filter(plot_df, Direction == "IC"),
                            "Intracranial-enriched")
plot_EC <- make_enrich_plot(filter(plot_df, Direction == "EC"),
                            "Extracranial-enriched")

plot_IC
plot_EC

top <- dge_tbl %>%
  dplyr::select(Fibroblast_IC, Fibroblast_EC)

gao_fib <- readr::read_csv("BMCA/csv/gao_fib.csv")
library(dplyr); library(tidyr)

wide <- gao_fib %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC), .by_group = TRUE) %>%
  mutate(rank = row_number()) %>%
  ungroup() %>%
  dplyr::select(cluster, rank, gene) %>%
  pivot_wider(names_from = cluster, values_from = gene)  %>%
  dplyr::select(-rank)
gao_fib <- wide[1:30,]
luo_fib <- readr::read_csv("BMCA/csv/luo_fib.csv")
library(dplyr); library(tibble)

# gene-name columns are every 3rd: positions 1, 4, 7, ...
gene_cols <- seq(1, ncol(luo_fib), by = 3)

# cluster IDs (c1, c2, ...) live in the first data row at those columns
clust_names <- as.character(unlist(luo_fib[1, gene_cols]))

luo_wide <- luo_fib[-1, gene_cols]          # drop the sub-header row
colnames(luo_wide) <- clust_names
luo_fib <- as_tibble(luo_wide)[,1:8]

cords_fib <- readr::read_csv("BMCA/csv/cords_fib.csv")
library(dplyr); library(tidyr)

cords_fib <- cords_fib %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC), .by_group = TRUE) %>%
  mutate(rank = row_number()) %>%
  ungroup() %>%
  dplyr::select(cluster, rank, gene) %>%
  pivot_wider(names_from = cluster, values_from = gene) %>%
  dplyr::select(-rank)
cords_fib <- cords_fib[1:30,]

fib_tables <- list(Gao = gao_fib, Luo = luo_fib, Cords = cords_fib)
fib_sets
jaccard  <- function(a, b) length(intersect(a, b)) / length(union(a, b))
top_sets <- lapply(top, function(x) unique(na.omit(x)))

for (nm in names(fib_tables)) {
  fib <- fib_tables[[nm]]
  fib_sets <- lapply(fib, function(x) unique(na.omit(x)))

  # rows = top groups, cols = fib clusters
  # rows = top groups, cols = fib clusters
  jac <- expand_grid(row = names(top_sets), col = names(fib_sets)) %>%
    mutate(
      jac_val   = purrr::map2_dbl(row, col,
                                  ~ jaccard(top_sets[[.x]], fib_sets[[.y]])),
      n_overlap = purrr::map2_int(row, col,
                                  ~ length(intersect(top_sets[[.x]], fib_sets[[.y]])))
    )
  # order by hierarchical clustering
  m <- jac %>%
    dplyr::select(row, col, jac_val) %>%
    pivot_wider(names_from = col, values_from = jac_val) %>%
    tibble::column_to_rownames("row") %>% as.matrix()
  row_ord <- rownames(m)[hclust(dist(m))$order]
  col_ord <- colnames(m)[hclust(dist(t(m)))$order]

  jac <- jac %>%
    mutate(row = factor(row, levels = row_ord),
           col = factor(col, levels = col_ord))

  p <- ggplot(jac, aes(col, row, fill = jac_val)) +
    geom_tile(colour = "grey90", linewidth = 0.1) +
    geom_text(aes(label = n_overlap,
                  colour = jac_val > max(jac_val, na.rm = TRUE) / 2),
              size = 2.6, show.legend = FALSE) +
    scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "grey20")) +
    scale_fill_gradient(low = "white", high = "darkred",
                        name = "Jaccard", limits = c(0, NA)) +
    coord_fixed() +
    labs(x = NULL, y = NULL,
         title = paste0(nm, " et al.")) +
    theme_minimal(base_size = 9) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid = element_blank())

  print(p)

  #ggsave(sprintf("BMCA/plots/jaccard_top_vs_%s.pdf", nm), p,
   #      width = 6, height = 5)
}
# -----------------------------------------------------------------------------


# genes for fibroblast (may be <1000)
g <- top_genes$Fibroblast

fib <- norm_pb$Fibroblast[g, ]
fib_c <- fib - rowMeans(fib)     # subtract each gene's mean across samples
cc <- cor(fib_c, method = "pearson")

# --- hierarchical clustering on 1 - correlation ---
hc <- hclust(as.dist(1 - cc), method = "average")
ord <- hc$labels[hc$order]

# --- long format, ordered by dendrogram ---
cc_long <- cc %>%
  as.data.frame() %>%
  rownames_to_column("sample_x") %>%
  pivot_longer(-sample_x, names_to = "sample_y", values_to = "pearson") %>%
  mutate(
    sample_x = factor(sample_x, levels = ord),
    sample_y = factor(sample_y, levels = ord)
  )

p <- ggplot(cc_long, aes(sample_x, sample_y, fill = pearson)) +
  geom_tile() +
  scale_fill_gradient2(
    low = "darkblue", mid = "white", high = "darkred",
    midpoint = 0,
    limits = c(-0.8, 0.8),
    oob = scales::squish,        # cap: values beyond limits clamp to end colour
    name = "Pearson"
  ) +
  coord_fixed() +
  labs(x = NULL, y = NULL,
       title = "Fibroblast — sample-sample Pearson (top genes)") +
  theme_minimal(base_size = 8) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 0),
    axis.text.y = element_text(size = 0),
    panel.grid = element_blank()
  )
p

metadata <- metadata %>%
  distinct()
metadata_p <- metadata %>%
  filter(patient %in% rownames(cc)) %>%
  dplyr::select(patient, primary_site, location, Donor, Location, "Dissamination Model") %>%
  distinct()
library(dplyr); library(tibble)

patient_to_donor <- metadata_p %>%
  filter(location == "Intracranial", !is.na(Donor)) %>%
  distinct(patient, Donor) %>%
  deframe()                       # names = patient, values = Donor

# invert: names = original Donor value, values = patient name
donor_to_patient <- setNames(names(patient_to_donor), patient_to_donor)

# recode the Donor column (apply to whichever table holds real Donor values)
metadata_p <- metadata_p %>%
  mutate(Donor = unname(donor_to_patient[Donor]))
metadata_p$Donor[metadata_p$Donor == "lin_sn_B2234330"] <- NA

metadata_p$location[metadata_p$location == "Intracranial"] <- "BrM"
metadata_p$location[metadata_p$location == "Extracranial"] <- "Primary"

colnames(metadata_p) <- c("Patient",             "Primary site",
                          "Location (BrM/Primary)",
                          "Paired Patient",               "Location in the brain",
                          "Dissamination model")

metadata_p

# --- dendrogram order (already computed) ---
hc  <- hclust(as.dist(1 - cc), method = "average")
ord <- hc$labels[hc$order]

# --- align metadata to the heatmap sample order ---
meta_ord <- metadata_p %>%
  filter(Patient %in% ord) %>%
  mutate(Patient = factor(Patient, levels = ord)) %>%
  arrange(Patient)

# ================= main heatmap =================
cc_long <- cc %>%
  as.data.frame() %>%
  rownames_to_column("sample_x") %>%
  pivot_longer(-sample_x, names_to = "sample_y", values_to = "pearson") %>%
  mutate(sample_x = factor(sample_x, levels = ord),
         sample_y = factor(sample_y, levels = ord))

p_heat <- ggplot(cc_long, aes(sample_x, sample_y, fill = pearson)) +
  geom_tile() +
  scale_fill_gradient2(low = "darkblue", mid = "white", high = "darkred",
                       midpoint = 0, limits = c(-0.8, 0.8),
                       oob = scales::squish, name = "Pearson") +
  coord_fixed() +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 8) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size = 0),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank(),
        plot.margin = margin(1, 1, 1, 1))

p_heat
# ================= annotation strip builder =================
anno_strip <- function(df, col, palette_name = NULL, cols = NULL) {
  g <- ggplot(df, aes(x = Patient, y = 1, fill = .data[[col]])) +
    geom_tile(height = 1) +
    coord_fixed(ratio = 4) +   # each tile drawn 4x taller than one sample is wide
    labs(x = NULL, y = col, fill = col) +
    theme_minimal(base_size = 12) +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          panel.grid = element_blank(),
          axis.title.y = element_text(angle = 0, hjust = 1, vjust = 0.5, size = 12),
          plot.margin = margin(0, 1, 0, 1))
  if (!is.null(cols))          g <- g + scale_fill_manual(values = cols, na.value = "grey90")
  else if (!is.null(palette_name)) g <- g + scale_fill_brewer(palette = palette_name, na.value = "grey90")
  g
}
a_site   <- anno_strip(meta_ord, "Primary site",        palette_name = "Set3")
a_loc2   <- anno_strip(meta_ord, "Location (BrM/Primary)",            palette_name = "Paired")
n_Paired_Patient <- dplyr::n_distinct(meta_ord$'Paired Patient'[!is.na(meta_ord$'Paired Patient')])
Paired_Patient_cols <- setNames(
  grDevices::hcl.colors(n_Paired_Patient, "Spectral"),   # or "Set3", "Dark2", etc.
  sort(unique(na.omit(meta_ord$'Paired Patient')))
)
a_Paired_Patient <- anno_strip(meta_ord, "Paired Patient", cols = Paired_Patient_cols)

a_site
a_loc2
a_Paired_Patient
p_heat

lk <- theme(legend.key.size = unit(5, "mm"),
            legend.text = element_text(size = 8),
            legend.title = element_text(size = 8))

ggsave("BMCA/plots/a_site_body.pdf",
       a_site  + theme(legend.position = "none"),
       width = 6, height = 0.5)
ggsave("BMCA/plots/a_site_legend.pdf",
       get_legend(a_site  + lk), width = 2,
       height = 2)
ggsave("BMCA/plots/a_loc2_body.pdf",
       a_loc2  + theme(legend.position = "none"),
       width = 6, height = 1)
ggsave("BMCA/plots/a_loc2_legend.pdf",
       get_legend(a_loc2  + lk), width = 2,
       height = 2)
ggsave("BMCA/plots/a_Paired Patient_body.pdf",
       a_Paired_Patient + theme(legend.position = "none"),
       width = 6, height = 1)
ggsave("BMCA/plots/a_Paired Patient_legend.pdf",
       get_legend(a_Paired_Patient + lk), width = 2,
       height = 2)
ggsave("BMCA/plots/p_heat_body.pdf",
       p_heat  + theme(legend.position = "none"),
       width = 6, height = 6)
ggsave("BMCA/plots/p_heat_legend.pdf",
       get_legend(p_heat  + lk), width = 2,
       height = 2)

library(dplyr); library(tidyr); library(ggplot2)

samples <- rownames(cc)
meta_m <- metadata_p %>% distinct(Patient, .keep_all = TRUE)
meta_m <- meta_m[match(samples, meta_m$Patient), ]
stopifnot(all(meta_m$Patient == samples))

attrs <- c(
  "Primary site"           = "Primary site",
  "Location (BrM/Primary)" = "Location (BrM/Primary)",
  "Paired Patient"                  = "Paired Patient"
)
not_shared <- c("Not specified")

cc_self <- cc; diag(cc_self) <- NA

# top-3 to same-attr AND top-3 to different-attr, per sample
top_long <- lapply(names(attrs), function(lbl) {
  v <- meta_m[[attrs[[lbl]]]]
  valid <- !is.na(v) & !(v %in% not_shared)

  do.call(rbind, lapply(seq_along(samples), function(k) {
    if (!valid[k]) return(NULL)
    same <- setdiff(which(valid & v == v[k]), k)
    diff <- setdiff(which(valid & v != v[k]), k)
    mk <- function(idx, grp) {
      if (length(idx) == 0) return(NULL)
      vals <- head(sort(cc_self[k, idx], decreasing = TRUE), 3)
      tibble(sample = samples[k], attribute = lbl, group = grp, cor = vals)
    }
    rbind(mk(same, "same"), mk(diff, "different"))
  }))
}) %>% bind_rows() %>%
  mutate(attribute = factor(attribute, levels = names(attrs)),
         group = factor(group, levels = c("same", "different")))

# ---- Wilcoxon same vs different, per attribute, BH across attributes ----
pvals <- top_long %>%
  group_by(attribute) %>%
  summarise(
    p = tryCatch(
      wilcox.test(cor[group == "same"], cor[group == "different"])$p.value,
      error = function(e) NA_real_),
    y = max(cor, na.rm = TRUE) + 0.03,
    .groups = "drop") %>%
  mutate(p_bh = p.adjust(p, method = "BH"),
         label = ifelse(is.na(p_bh), "ns",
                        paste0("p=", signif(p_bh, 2))))

p <- ggplot(top_long, aes(attribute, cor, fill = group)) +
  geom_boxplot(outlier.shape = NA, width = 0.6,
               position = position_dodge(0.7)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.12,
                                             dodge.width = 0.7),
             size = 0.5, alpha = 0.5) +
  geom_text(data = pvals, aes(attribute, y, label = label),
            inherit.aes = FALSE, size = 2.8) +
  scale_fill_manual(
    values = c(same = "firebrick", different = "grey70"),
    labels = c(same = "Within group", different = "Outside of group")
  ) +  labs(x = NULL, y = "Top Pearson correlation", fill = NULL,
       title = "") +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))
p
ggsave("BMCA/plots/fibro_top3_sameVSdiff_boxplot.pdf", p, width = 7, height = 4)

pvals   # inspect raw + BH p-values

library(dplyr); library(tibble)

# 3 groups from the dendrogram
grp <- cutree(hc, k = 3)              # named vector: sample -> 1/2/3
table(grp)                            # group sizes

# full normalized fibroblast matrix: genes x samples
fib_all <- as.matrix(fib_c)
fib_all <- fib_all[, names(grp)]      # align columns to grp order

# per-group log2FC (group mean - rest mean) on log2-CPM
top <- lapply(sort(unique(grp)), function(g) {
  in_g  <- names(grp)[grp == g]
  out_g <- names(grp)[grp != g]
  l2fc  <- rowMeans(fib_all[, in_g, drop = FALSE]) -
    rowMeans(fib_all[, out_g, drop = FALSE])
  names(sort(l2fc, decreasing = TRUE))[1:100]
})
names(top) <- paste0("group", sort(unique(grp)))


# -----------------------------------------------------------------------------


centered_cpm <- readRDS(paste0("BMCA/data/centered_cpm_lin_sn.RDS"))

metadata_fibro <- metadata %>%
  filter(cell_type == "Fibroblast") %>%
  filter(CellID %in% colnames(centered_cpm))
centered_cpm <- centered_cpm[,metadata_fibro$CellID]
pca <- prcomp(t(centered_cpm), center = F, scale = T, rank = 50)

metadata_fibro$PC1 <- pca$x[,1]
metadata_fibro$PC2 <- pca$x[,2]

ggplot(data = metadata_fibro, aes(x = PC1, y = PC2)) +
  geom_point(size = 0.5) +
  xlab("PC1") +
  ylab("PC2") +
  theme_classic() +
  theme(
    axis.text = element_text(size=12),
    axis.title = element_text(size=12),
    plot.title = element_text(size=14),
    legend.text = element_text(size=10),
    legend.title = element_text(size = 12),
    legend.key = element_blank(),
    panel.background = element_rect(fill = "white"),
    axis.line = element_line(colour = "black"),
    aspect.ratio = 1
  )

ggplot(data = metadata_fibro, aes(x = PC1, y = PC2, color = primary_site)) +
  geom_point(size = 0.5) +
  xlab("PC1") +
  ylab("PC2") +
  theme_classic() +
  theme(
    axis.text = element_text(size=12),
    axis.title = element_text(size=12),
    plot.title = element_text(size=14),
    legend.text = element_text(size=10),
    legend.title = element_text(size = 12),
    legend.key = element_blank(),
    panel.background = element_rect(fill = "white"),
    axis.line = element_line(colour = "black"),
    aspect.ratio = 1
  )

ggplot(data = metadata_fibro, aes(x = PC1, y = PC2, color = location)) +
  geom_point(size = 0.5) +
  xlab("PC1") +
  ylab("PC2") +
  theme_classic() +
  theme(
    axis.text = element_text(size=12),
    axis.title = element_text(size=12),
    plot.title = element_text(size=14),
    legend.text = element_text(size=10),
    legend.title = element_text(size = 12),
    legend.key = element_blank(),
    panel.background = element_rect(fill = "white"),
    axis.line = element_line(colour = "black"),
    aspect.ratio = 1
  )


# -----------------------------------------------------------------------------
## Chunk 2 ##
# -----------------------------------------------------------------------------

umap <- uwot::tumap(X = pca$x[, 1:50])

## save umap coordinates
metadata_fibro$UMAP1 <- umap[,1]
metadata_fibro$UMAP2 <- umap[,2]

## plot basic umap
ggplot(data = metadata_fibro, aes(x = UMAP1, y = UMAP2)) +
  geom_point(size = 0.1) +
  xlab("UMAP1") +
  ylab("UMAP2") +
  theme_classic() +
  theme(axis.text = element_text(size=20), axis.title = element_text(size=20),
        plot.title = element_text(size=24), legend.text = element_text(size=16),
        legend.key=element_blank(), panel.background = element_rect(fill = "white"),
        axis.line = element_line(colour = "black"), aspect.ratio = 1)

ggplot(data = metadata_fibro, aes(x = UMAP1, y = UMAP2, color = primary_site)) +
  geom_point(size = 0.01) +
  xlab("UMAP1") +
  ylab("UMAP2") +
  theme_classic() +
  guides(colour = guide_legend(override.aes = list(size = 1))) +
  theme(axis.text = element_text(size=5), axis.title = element_text(size=5),
        legend.key = element_blank(), panel.background = element_rect(fill = "white"),
        axis.line = element_line(colour = "black"), aspect.ratio = 1)

ggplot(data = metadata_fibro, aes(x = UMAP1, y = UMAP2, color = location)) +
  geom_point(size = 0.00001) +
  xlab("UMAP1") +
  ylab("UMAP2") +
  theme_classic() +
  guides(colour = guide_legend(override.aes = list(size = 1))) +
  theme(axis.text = element_text(size=5), axis.title = element_text(size=5),
        legend.key = element_blank(), panel.background = element_rect(fill = "white"),
        axis.line = element_line(colour = "black"), aspect.ratio = 1)

# -----------------------------------------------------------------------------
## Chunk 3 ##
# -----------------------------------------------------------------------------

metadata_fibro$km_cluster <- as.factor(
  paste0(
    "KM",
    kmeans(
      x = pca$x[, 1:50],
      centers = 5,
      iter.max = 1000
    )$cluster
  )
)

ggplot(data = metadata_fibro, aes(x = UMAP1, y = UMAP2, color = km_cluster)) +
  geom_point(size = 0.01) +
  xlab("UMAP1") +
  ylab("UMAP2") +
  theme_classic() +
  guides(colour = guide_legend(override.aes = list(size = 1))) +
  theme(axis.text = element_text(size=5), axis.title = element_text(size=5),
        legend.key = element_blank(), panel.background = element_rect(fill = "white"),
        axis.line = element_line(colour = "black"), aspect.ratio = 1)

# -----------------------------------------------------------------------------
## Chunk 4 ##
# -----------------------------------------------------------------------------

cell_type_markers <- read_csv("BMCA/csv/canonical_markers2.csv")

centered_cpm_markers <- centered_cpm[rownames(centered_cpm) %in% cell_type_markers$gene,]

cluster_mat <- sapply(
  levels(metadata_fibro$km_cluster),
  function(cl) {
    rowMeans(centered_cpm_markers[, metadata_fibro$km_cluster == cl, drop = FALSE])
  }
)


cluster_mat[1:5, 1:5]
pheatmap(
  cluster_mat,
  clustering_distance_rows = "correlation",
  clustering_distance_cols = "correlation",
  clustering_method = "average"
)

# Genes to exclude
exclude_genes <- c("COL1A1", "COL1A2", "COL3A1", "DCN", "LUM")

# Subset the matrix
cluster_mat_sub <- cluster_mat[!rownames(cluster_mat) %in% exclude_genes, ]

# Convert to long format
df_long <- melt(cluster_mat_sub)

# Overlayed histogram per column
ggplot(df_long, aes(x = value, fill = cluster)) +
  geom_histogram(alpha = 0.3, position = "identity", bins = 100) +
  theme_minimal() +
  labs(title = "Overlayed histograms per column (excluding specific genes)",
       x = "Value", y = "Count") +
  theme(legend.position = "none")  # remove legend if too cluttered

# Genes to exclude
exclude_genes <- c("COL1A1", "COL1A2", "COL3A1", "DCN", "LUM")

# Subset matrix to exclude those genes
cluster_mat_sub <- cluster_mat[!rownames(cluster_mat) %in% exclude_genes, ]

# Logical matrix: TRUE if value < 2
lt2_mat <- cluster_mat_sub < 2

# Columns where all remaining rows are < 2
cols_all_below2 <- colnames(cluster_mat_sub)[colSums(lt2_mat) == nrow(cluster_mat_sub)]

# Result
cols_all_below2
length(cols_all_below2)  # how many columns meet this criteria

metadata_fibro <- metadata_fibro %>% filter(km_cluster %in%cols_all_below2)
# Count how many cells per patient
table_patients <- table(metadata_fibro$patient)

# Bar plot
barplot(table_patients,
        main = "Number of cells per patient",
        xlab = "Patient",
        ylab = "Count",
        las = 2,         # rotate x labels
        col = "steelblue")

expmat_cpm <- tryCatch(
  {
    readRDS(paste0("BMCA/data/expmat_cpm_lin_sn.RDS"))
  },
  error = function(e) {
    message("⚠️ Failed to load file for study: ", s)
    return(NULL)
  }
)

cols <- intersect(colnames(expmat_cpm), metadata_fibro$CellID)

# Skip if no overlapping cells
if (length(cols) == 0) {
  message("⚠️ No matching cells found for study: ", s)
  next
}

expmat_cpm_fil <- expmat_cpm[, cols]

metadata_fibro <- metadata_fibro %>%
  filter(CellID %in% cols)

# Split column indices by patient
patient_list <- split(metadata_fibro$CellID, metadata_fibro$patient)

# Keep only patients with >= 50 cells
patient_list <- patient_list[lengths(patient_list) >= 30]
# Skip study if no patients left
if (length(patient_list) == 0) {
  message("⚠️ No patients with ≥50 cells for study: ", s, " (cell type: ", ct, ")")
  next
}

# Initialize list to collect average expression per patient
patient_expr_list <- list()

for (p in names(patient_list)) {
  cell_ids <- patient_list[[p]]
  mat_subset <- expmat_cpm_fil[, cell_ids, drop = FALSE]

  # Average across columns (cells), keeping it sparse
  avg_expr <- Matrix::rowMeans(mat_subset)

  # Turn into a dgCMatrix with one column
  avg_mat <- Matrix::Matrix(avg_expr, sparse = TRUE)
  colnames(avg_mat) <- p  # name the column as patient ID
  rownames(avg_mat) <- rownames(expmat_cpm_fil)

  patient_expr_list[[p]] <- avg_mat
}

# Combine all patient columns into one matrix
patient_avg <- do.call(cbind, patient_expr_list)
# 1️⃣ Compute average expression per gene (rowMeans works on sparse matrices)
gene_avg <- rowMeans(patient_avg)

# 2️⃣ Order genes by average expression (descending)
top_genes <- names(sort(gene_avg, decreasing = TRUE))[1:1000]

expmat_patient_avg_Fib <- patient_avg[top_genes,]
# -----------------------------------------------------------------------------
metadata <- as_tibble(readRDS("BMCA/data/metadata_all_studies.rds"))

anno_patient <- metadata %>%
  distinct(patient, primary_site, Donor, location, study, seq_tech) %>%
  filter(study == "lin_sn") %>%
  distinct()
patients <- intersect(anno_patient$patient, colnames(expmat_patient_avg_Fib))

expmat_patient_avg <- log2(expmat_patient_avg_Fib + 1)
expmat_patient_avg <- expmat_patient_avg - rowMeans(expmat_patient_avg)
expmat_patient_avg_mat <- as.matrix(expmat_patient_avg)
cor_mat <- cor(
  expmat_patient_avg_mat,
  method = "spearman",
  use = "pairwise.complete.obs"
)

anno_row <- anno_patient |>
  as.data.frame()

rownames(anno_row) <- anno_row$patient
anno_row$patient <- NULL
anno_row <- anno_row[rownames(cor_mat), , drop = FALSE]
stopifnot(identical(rownames(anno_row), rownames(cor_mat)))
anno_row[] <- lapply(anno_row, as.factor)


pheatmap(
  cor_mat,
  color = colorRampPalette(c("darkblue", "white", "darkred"))(100),
  clustering_distance_rows = as.dist(1 - cor_mat),
  clustering_distance_cols = as.dist(1 - cor_mat),
  clustering_method = "average",
  show_rownames = FALSE,
  show_colnames = FALSE,
  annotation_row = anno_row,
  border_color = NA
)

# Compute HC from correlation matrix
hc <- hclust(as.dist(1 - cor_mat), method = "average")

# Cut into 2 clusters
clusters <- cutree(hc, k = 2)

# Add cluster info to annotation
anno_row$cluster <- factor(clusters[match(rownames(anno_row), names(clusters))])
# make cluster factor
cluster_factor <- anno_row$cluster
mat_sub <- expmat_patient_avg_mat[, patients]
# run Wilcoxon test for each gene
de_results <- apply(mat_sub, 1, function(x) {
  # x = expression of one gene across patients
  wilcox.test(x ~ cluster_factor)$p.value
})

# compute fold change between clusters
fc <- rowMeans(mat_sub[, cluster_factor == 2]) - rowMeans(mat_sub[, cluster_factor == 1])

# assemble results table
res <- data.frame(
  gene = rownames(mat_sub),
  logFC = fc,
  pval  = de_results
)


# adjust p-values
res$padj <- p.adjust(res$pval, method = "BH")

ggplot(res, aes(x = logFC, y = -log10(padj))) +
  geom_point(alpha = 0.6) +
  geom_vline(xintercept = c(-1,1), linetype = "dashed", color = "red") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  theme_minimal() +
  labs(
    x = "Log2 fold change (cluster2 vs cluster1)",
    y = "-log10(adj. p-value)",
    title = "Volcano plot of DE genes"
  )


anno_row_1 <- anno_row %>%
  filter(cluster == "1") %>%
  select(Donor, location, primary_site)

expmat_patient_avg_mat_1 <- expmat_patient_avg_mat[,rownames(anno_row_1)]- rowMeans(expmat_patient_avg_mat[,rownames(anno_row_1)])
cor_mat_1 <- cor(
  expmat_patient_avg_mat_1,
  method = "spearman",
  use = "pairwise.complete.obs"
)
pheatmap(
  cor_mat_1,
  color = colorRampPalette(c("darkblue", "white", "darkred"))(100),
  clustering_distance_rows = as.dist(1 - cor_mat_1),
  clustering_distance_cols = as.dist(1 - cor_mat_1),
  clustering_method = "average",
  show_rownames = FALSE,
  show_colnames = FALSE,
  annotation_row = anno_row_1,
  border_color = NA
)

# symmetric breaks around 0
max_abs <- max(abs(cor_mat_1), na.rm = TRUE)
breaks <- seq(-max_abs, max_abs, length.out = 101)

# Define colors only for location
ann_colors <- list(
  location = c(
    Extracranial = "orange",
    Intracranial = "purple"
  )
)

pheatmap(
  cor_mat_1,
  color = colorRampPalette(c("darkblue", "white", "darkred"))(100),
  breaks = breaks,
  clustering_distance_rows = as.dist(1 - cor_mat_1),
  clustering_distance_cols = as.dist(1 - cor_mat_1),
  clustering_method = "average",
  show_rownames = FALSE,
  show_colnames = FALSE,
  annotation_row = anno_row_1,
  annotation_colors = ann_colors,  # only location controlled
  border_color = NA
)

pheatmap(
  cor_mat_1,
  color = colorRampPalette(c("darkblue", "white", "darkred"))(100),
  breaks = breaks,
  clustering_distance_rows = as.dist(1 - cor_mat_1),
  clustering_distance_cols = as.dist(1 - cor_mat_1),
  clustering_method = "average",
  show_rownames = FALSE,
  show_colnames = FALSE,
  annotation_row = anno_row_1,
  annotation_colors = ann_colors,
  border_color = NA,
  # Add these three lines:
  filename = "BMCA/plots/fig1/correlation_heatmap.pdf",
  width = 11,
  height = 10
)
r <- 3
pheatmap(
  cor_mat_1,
  color = colorRampPalette(c("darkblue", "white", "darkred"))(100),
  breaks = breaks,
  clustering_distance_rows = as.dist(1 - cor_mat_1),
  clustering_distance_cols = as.dist(1 - cor_mat_1),
  clustering_method = "average",
  show_rownames = FALSE,
  show_colnames = FALSE,
  annotation_row = anno_row_1,
  annotation_colors = ann_colors,  # only location controlled
  border_color = NA,
  cutree_rows = r,   # <-- number of row clusters
  cutree_cols = r    # <-- number of column clusters
)

# -----------------------------------------------------------------------------

# Set grouping variable
group_var <- "primary_site"  # could be "location" or "primary_site"

# Process and plot
as.data.frame(as.table(cor_mat_1)) %>%
  setNames(c("Sample1", "Sample2", "Correlation")) %>%
  filter(Sample1 != Sample2) %>%  # remove self-correlations
  left_join(anno_row_1 %>% mutate(Sample = rownames(cor_mat_1)) %>% select(Sample, all_of(group_var)),
            by = c("Sample1" = "Sample")) %>%
  rename(Group1 = all_of(group_var)) %>%
  left_join(anno_row_1 %>% mutate(Sample = rownames(cor_mat_1)) %>% select(Sample, all_of(group_var)),
            by = c("Sample2" = "Sample")) %>%
  rename(Group2 = all_of(group_var)) %>%
  filter(Group1 == Group2 & !is.na(Group1)) %>%  # keep only within-group correlations
  group_by(Sample1) %>%
  slice_max(order_by = Correlation, n = 5) %>%   # keep up to 5 highest per sample
  ungroup() %>%
  ggplot(aes(x = Group1, y = Correlation)) +
  geom_boxplot(fill = "lightblue") +
  theme_minimal() +
  labs(x = group_var, y = "Correlation", title = paste("Top 5 within-group correlations by", group_var)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


n = 3
# Function to get top N correlations within a group
get_top_cor <- function(cor_mat, anno_df, group_var, top_n = n) {
  as.data.frame(as.table(cor_mat)) %>%
    setNames(c("Sample1", "Sample2", "Correlation")) %>%
    filter(Sample1 != Sample2) %>%  # remove self-correlations
    left_join(anno_df %>% mutate(Sample = rownames(cor_mat)) %>% select(Sample, all_of(group_var)),
              by = c("Sample1" = "Sample")) %>%
    rename(Group1 = all_of(group_var)) %>%
    left_join(anno_df %>% mutate(Sample = rownames(cor_mat)) %>% select(Sample, all_of(group_var)),
              by = c("Sample2" = "Sample")) %>%
    rename(Group2 = all_of(group_var)) %>%
    filter(Group1 == Group2 & !is.na(Group1)) %>%  # within-group only
    group_by(Sample1) %>%
    slice_max(order_by = Correlation, n = top_n) %>%
    ungroup() %>%
    mutate(GroupType = group_var)  # label for plotting
}

# Get top n correlations for each grouping
top_location <- get_top_cor(cor_mat_1, anno_row_1, "location", top_n = n)
top_primary <- get_top_cor(cor_mat_1, anno_row_1, "primary_site", top_n = n)
top_donor   <- get_top_cor(cor_mat_1, anno_row_1, "Donor", top_n = n)

# Combine all
top_all <- bind_rows(top_location, top_primary, top_donor)

# Plot combined boxplot
ggplot(top_all, aes(x = GroupType, y = Correlation, fill = GroupType)) +
  geom_boxplot() +
  theme_minimal() +
  labs(x = "Grouping variable", y = "Top n correlations per sample",
       title = "Comparison of top correlations across groups") +
  scale_fill_manual(values = c("location" = "skyblue",
                               "primary_site" = "lightgreen",
                               "Donor" = "salmon")) +
  theme(legend.position = "none")


# Use the same top_all from before
# top_all <- bind_rows(top_location, top_primary, top_donor)

# Plot with p-values
ggplot(top_all, aes(x = GroupType, y = Correlation, fill = GroupType)) +
  geom_boxplot() +
  theme_minimal() +
  labs(x = "Grouping variable", y = "Top 3 correlations per sample",
       title = "Comparison of top correlations across groups") +
  scale_fill_manual(values = c("location" = "skyblue",
                               "primary_site" = "lightgreen",
                               "Donor" = "salmon")) +
  theme(legend.position = "none") +
  stat_compare_means(method = "t.test",
                     label = "p.signif",
                     comparisons = list(c("location", "primary_site"),
                                        c("location", "Donor"),
                                        c("primary_site", "Donor")))
# 1. Assign your plot to an object
p_corr <- ggplot(top_all, aes(x = GroupType, y = Correlation, fill = GroupType)) +
  geom_boxplot() +
  theme_minimal() +
  labs(x = "", y = "Average Correlation Within A Group",
       title = "") +
  scale_fill_manual(values = c("location" = "skyblue",
                               "primary_site" = "lightgreen",
                               "Donor" = "salmon")) +
  theme(legend.position = "none") +
  stat_compare_means(method = "t.test",
                     label = "p.signif",
                     comparisons = list(c("location", "primary_site"),
                                        c("location", "Donor"),
                                        c("primary_site", "Donor")))

# 2. Save to PDF
ggsave(
  filename = "BMCA/plots/fig1/Top_Correlations_Comparison.pdf",
  plot = p_corr,
  width = 7,       # Width in inches
  height = 6,      # Height in inches
  device = "pdf"
)
## -------------------------------
## Hierarchical clustering
## -------------------------------
dist_mat <- as.dist(1 - cor_mat_1)
hc <- hclust(dist_mat, method = "average")
clusters <- cutree(hc, k = r)


## -------------------------------
## Align expression matrix
## -------------------------------
expmat <- expmat_patient_avg_mat_1[, names(clusters)]

## -------------------------------
## DGE: cluster vs all others
## -------------------------------
dge_list <- lapply(1:r, function(k) {

  in_clust  <- clusters == k
  out_clust <- clusters != k

  logFC <- rowMeans(expmat[, in_clust, drop = FALSE]) -
    rowMeans(expmat[, out_clust, drop = FALSE])

  pval <- apply(expmat, 1, function(g)
    wilcox.test(g[in_clust], g[out_clust], exact = FALSE)$p.value
  )

  data.frame(
    gene  = rownames(expmat),
    logFC = logFC,
    pval  = pval,
    padj  = p.adjust(pval, method = "BH"),
    cluster = k,
    stringsAsFactors = FALSE
  )
})


names(dge_list) <- paste0("Cluster_", 1:r)
res_2 <- dge_list[[1]]
top20_per_cluster <- lapply(dge_list, function(df) {
  # keep only significant genes
  sig_genes <- df[df$padj < 0.05, ]

  # order by logFC descending and take top 20
  head(sig_genes[order(sig_genes$logFC, decreasing = TRUE), ], 20)
})

