# =============================================================================
# Figure 4 - TME composition and cell-state shifts between intracranial and extracranial tumours
# =============================================================================
#
# Cell-type composition comparisons (unpaired and donor-paired), myeloid/lymphoid
# ratio, cell-state delta-z heatmaps, IFNG signalling in T/NK cells, and Jaccard
# comparisons against published myeloid meta-program sets.
#
# Part of: A brain metastasis cell atlas reveals multicellular adaptation to
#          the neural microenvironment (Feigin et al.)
#
# Run from the project root (see README.md). All paths below are relative to it.
#
# Inputs:
#   BMCA/csv/joyce_mps.csv
#   BMCA/csv/tyler_myeloid_MP.csv
#   BMCA/data/celltype_patient_matrix_proccesed.RDS
#   BMCA/data/ecotypes_tib_proccesed_filtered.RDS
#   BMCA/data/md_tib_Myeloid.RDS
#   BMCA/data/metadata_all_studies.rds
#   BMCA/data/metadata_states_proccesed_filtered.RDS
#   BMCA/data/score_tib_Myeloid.RDS
#   BMCA/data/UMI_pb_T_NK.RDS
#
# Outputs:
#   BMCA/plots/fig4/  (12 files)
#
# Original filename: figure_4_1.R
# =============================================================================

# ---- Packages ---------------------------------------------------------------
pkgs <- c(
  "dplyr", "ggh4x", "ggplot2", "ggpubr", "ggrepel", "patchwork", "purrr",
  "readr", "readxl", "reshape2", "rstatix", "scales", "stringr", "tidyr",
  "tidyverse", "viridis"
)
invisible(lapply(pkgs, library, character.only = TRUE))

# ---- Project paths ----------------------------------------------------------
source("R/config.R")

# =============================================================================
metadata <- as_tibble(readRDS("BMCA/data/metadata_all_studies.rds"))

metadata_p <- metadata %>%
  dplyr::select(patient, primary_site, location, study, seq_tech) %>%
  distinct()

metadata_mes <- metadata %>%
  filter(cell_type %in% c("Fibroblast", "Pericyte") &
           location == "Intracranial")

metadata_p_fil <- metadata_p %>%
  filter(primary_site %in% c("melanoma", "NSCLC", "colorectal", "breast"))


# Count primary_site within each location × seq_tech panel
plot_df <- metadata_p_fil %>%
  count(location, seq_tech, primary_site, name = "n")
plot_df$seq_tech[plot_df$seq_tech == "sc"] = "Single-cell"
plot_df$seq_tech[plot_df$seq_tech == "sn"] = "Single-nucleus"

p <- ggplot(plot_df,
            aes(x = reorder(primary_site, -n), y = n, fill = primary_site)) +
  geom_col(width = 0.72) +
  geom_text(aes(label = n), vjust = -0.35, size = 2) +
  facet_grid(location ~ seq_tech,
             scales = "free_x", space = "free_x") +
  scale_fill_brewer(palette = "Set2", guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(x = "Primary site", y = "Number of samples") +
  theme_classic(base_size = 7, base_family = "Helvetica") +
  theme(
    axis.text.x     = element_text(angle = 45, hjust = 1, colour = "black"),
    axis.text.y     = element_text(colour = "black"),
    axis.line       = element_line(linewidth = 0.3, colour = "black"),
    axis.ticks      = element_line(linewidth = 0.3, colour = "black"),
    strip.background = element_rect(fill = "grey92", colour = NA),
    strip.text      = element_text(size = 7),
    panel.spacing   = unit(0.6, "lines"),
    plot.margin     = margin(4, 6, 4, 4)
  )
p

ggsave("BMCA/plots/fig4/primary_site_faceted.pdf", p,
       width = 60, height = 65, units = "mm", dpi = 300)
# -----------------------------------------------------------------------------
celltype_patient_matrix <- readRDS("BMCA/data/celltype_patient_matrix_proccesed.RDS")

celltype_patient_matrix <- celltype_patient_matrix %>%
  left_join(metadata_p, by = "patient")
celltype_patient_matrix <- celltype_patient_matrix%>%
  filter(!(study == "anoop")) %>%
  filter(primary_site %in% c("melanoma", "NSCLC", "breast", "colorectal"))


celltype_patient_matrix %>%
  count(primary_site, location, seq_tech) %>%
  ggplot(aes(x = location, y = n, fill = primary_site)) +
  geom_col(position = "stack", width = 0.7) +
  facet_wrap(~ seq_tech, scales = "free_x") +
  scale_fill_brewer(palette = "Set2") +
  labs(
    x = "Location",
    y = "Number of patients",
    fill = "Primary site",
    title = "Distribution of primary site by location",
    subtitle = "Faceted by sequencing technology"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_rect(fill = "grey90"),
    strip.text = element_text(face = "bold"),
    legend.position = "bottom"
  )

# -----------------------------------------------------------------------------
# ── Groupings ────────────────────────────────────────────────────────────────
immune_cells  <- c("Myeloid", "B_Plasma", "T_NK", "CD4", "CD8", "GD", "NK")
stromal_cells <- c("Fibroblast", "Pericyte", "Endothelial")
cell_features <- c("Malignant", immune_cells, stromal_cells)

carcinoma_sites <- c("NSCLC", "breast", "colorectal")
# ── Long format + recode primary_site ────────────────────────────────────────
long_data <- celltype_patient_matrix %>%
  select(patient, primary_site, location, seq_tech, all_of(cell_features)) %>%
  mutate(
    primary_site = if_else(primary_site %in% carcinoma_sites, "Carcinoma", primary_site),
    primary_site = recode(primary_site, "melanoma" = "Melanoma")
  ) %>%
  pivot_longer(cols = all_of(cell_features),
               names_to = "cell_type",
               values_to = "value") %>%
  mutate(
    cell_group = case_when(
      cell_type == "Malignant"          ~ "Malignant",
      cell_type %in% immune_cells       ~ "Immune",
      cell_type %in% stromal_cells      ~ "Stromal"
    ),
    cell_group   = factor(cell_group,   levels = c("Malignant", "Immune", "Stromal")),
    cell_type    = factor(cell_type,    levels = cell_features),
    seq_tech     = factor(seq_tech,     levels = c("sc", "sn")),
    primary_site = factor(primary_site, levels = c("Carcinoma", "Melanoma"))
  )

# ── Wilcoxon + BH correction ─────────────────────────────────────────────────
stats_raw <- long_data %>%
  group_by(cell_type, cell_group, primary_site, seq_tech) %>%
  filter(n_distinct(location) == 2) %>%
  wilcox_test(value ~ location) %>%
  ungroup()

stats_corrected <- stats_raw %>%
  mutate(p.adj = p.adjust(p, method = "BH")) %>%
  mutate(significance = case_when(
    p.adj < 0.001 ~ "***",
    p.adj < 0.01  ~ "**",
    p.adj < 0.05  ~ "*",
    TRUE          ~ "ns"
  ))

# ── Summary stats (SE) ────────────────────────────────────────────────────────
plot_data <- long_data %>%
  group_by(primary_site, seq_tech, location, cell_type, cell_group) %>%
  summarise(
    mean_val = mean(value, na.rm = TRUE),
    se_val   = sd(value, na.rm = TRUE) / sqrt(sum(!is.na(value))),
    n        = n(),
    .groups  = "drop"
  )

# ── Bracket positions ─────────────────────────────────────────────────────────
bracket_data <- plot_data %>%
  group_by(cell_type, cell_group, primary_site, seq_tech) %>%
  summarise(
    y_pos = max(mean_val + se_val, na.rm = TRUE) * 1.08,
    .groups = "drop"
  ) %>%
  left_join(
    stats_corrected %>% select(cell_type, primary_site, seq_tech, p.adj, significance),
    by = c("cell_type", "primary_site", "seq_tech")
  ) %>%
  filter(!is.na(significance), significance != "ns")

p <- ggplot(plot_data, aes(x = seq_tech, y = mean_val, fill = location)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.65) +
  geom_errorbar(
    aes(ymin = pmax(mean_val - se_val, 0), ymax = mean_val + se_val),
    position = position_dodge(width = 0.8),
    width = 0.2,
    linewidth = 0.45,
    color = "grey20"
  ) +
  geom_text(
    data = bracket_data,
    aes(x = seq_tech, y = y_pos, label = significance),
    inherit.aes = FALSE,
    size = 2,
    fontface = "bold",
    color = "grey20"
  ) +
  facet_nested(
    primary_site ~ cell_group + cell_type,
    scales    = "free_y",
    space     = "fixed",
    nest_line = element_line(color = "grey40", linewidth = 0.4)
  ) +
  scale_fill_manual(values = c("Intracranial" = "#C0392B", "Extracranial" = "#2980B9")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +  # ← clean base, room for stars
  labs(
    x    = "Sequencing technology",
    y    = "Mean abundance (%)",
    fill = NULL   # ← cleaner: no legend title, just colour swatches
  ) +
  theme_classic(base_size = 7, base_family = "Helvetica") +  # ← Nature uses Helvetica
  theme(
    # ── Facet strips ──────────────────────────────────────────────────────────
    strip.text.x         = element_text(color = "grey10"),
    strip.text.y         = element_text(color = "grey10", angle = 0),
    strip.background     = element_rect(fill = "grey95", color = NA),
    ggh4x.facet.nestline = element_line(colour = "grey50", linewidth = 0.3),

    # ── Axes ──────────────────────────────────────────────────────────────────
    axis.text.x          = element_text(color = "grey10"),
    axis.text.y          = element_text(color = "grey10"),
    axis.title.x         = element_text(face = "bold", color = "grey10",
                                        margin = margin(t = 6)),
    axis.title.y         = element_text(face = "bold", color = "grey10",
                                        margin = margin(r = 6)),
    axis.line            = element_line(linewidth = 0.4, color = "grey30"),
    axis.ticks           = element_line(linewidth = 0.3, color = "grey30"),

    # ── Legend ────────────────────────────────────────────────────────────────
    legend.position      = "bottom",
    legend.text          = element_text(color = "grey10"),
    legend.key.size      = unit(0.35, "cm"),
    legend.spacing.x     = unit(0.2, "cm"),

    # ── Panel ─────────────────────────────────────────────────────────────────
    panel.grid.major.y   = element_line(linewidth = 0.25, color = "grey90"),  # ← subtle y grid
    panel.grid.major.x   = element_blank(),
    panel.grid.minor     = element_blank(),
    panel.spacing        = unit(0.3, "lines"),

    # ── Overall ───────────────────────────────────────────────────────────────
    plot.margin          = margin(8, 8, 8, 8)
  )

p

ggsave("BMCA/plots/fig4/cell_types_bars_all.pdf", p,
       width = 180, height = 80, units = "mm", dpi = 300)
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
metadata_donor <- metadata %>%
  filter(!is.na(Donor)) %>%
  dplyr::select(patient, Donor) %>%
  distinct()


immune_cells  <- c("Myeloid", "B_Plasma", "T_NK", "CD4", "CD8", "GD", "NK")
stromal_cells <- c("Fibroblast", "Pericyte", "Endothelial")
cell_features <- c("Malignant", immune_cells, stromal_cells)

celltype_patient_matrix_donor <- celltype_patient_matrix %>%
  left_join(metadata_donor, by = "patient") %>%
  filter(!is.na(Donor))
# Long format, carry Donor through; one value per Donor × location × cell_type
paired_data <- celltype_patient_matrix_donor %>%
  select(Donor, location, all_of(cell_features)) %>%
  pivot_longer(all_of(cell_features), names_to = "cell_type", values_to = "value") %>%
  group_by(Donor, location, cell_type) %>%
  summarise(value = mean(value, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    cell_type = factor(cell_type, levels = cell_features),
    location  = factor(location,  levels = c("Extracranial", "Intracranial"))
  )

# Keep only donors with BOTH locations for each cell_type (complete pairs)
complete_pairs <- paired_data %>%
  group_by(cell_type, Donor) %>%
  filter(n_distinct(location) == 2) %>%
  ungroup()

# Paired Wilcoxon signed-rank per cell_type (arrange ensures correct pairing)
stats_paired <- complete_pairs %>%
  arrange(cell_type, Donor, location) %>%
  group_by(cell_type) %>%
  wilcox_test(value ~ location, paired = TRUE) %>%
  ungroup() %>%
  mutate(
    p.adj = p.adjust(p, method = "BH"),
    significance = case_when(
      p.adj < 0.001 ~ "***", p.adj < 0.01 ~ "**",
      p.adj < 0.05  ~ "*",   TRUE         ~ "ns"
    )
  )

# Star positions (top of each free-y panel)
star_pos <- complete_pairs %>%
  group_by(cell_type) %>%
  summarise(y_pos = max(value, na.rm = TRUE) * 1.05, .groups = "drop") %>%
  left_join(stats_paired %>% select(cell_type, significance, p.adj),
            by = "cell_type") %>%
  filter(significance != "ns")

p <- ggplot(complete_pairs, aes(x = location, y = value)) +
  geom_line(aes(group = Donor), colour = "grey75", linewidth = 0.3, alpha = 0.6) +
  geom_point(aes(colour = location), size = 1, alpha = 0.85) +
  geom_text(data = star_pos, aes(x = 1.5, y = y_pos, label = significance),
            inherit.aes = FALSE, size = 2, fontface = "bold") +
  facet_wrap(~ cell_type, scales = "free_y", nrow = 1) +
  scale_colour_manual(values = c("Extracranial" = "#2980B9",
                                 "Intracranial" = "#C0392B"), guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.1))) +
  labs(x = NULL, y = "Abundance (%)") +
  theme_classic(base_size = 7, base_family = "Helvetica") +
  theme(
    axis.text   = element_text(colour = "black"),
    axis.text.x = element_text(angle = 30, hjust = 1),
    axis.line   = element_line(linewidth = 0.3),
    axis.ticks  = element_line(linewidth = 0.3),
    strip.background = element_rect(fill = "grey95", colour = NA),
    panel.spacing = unit(0.5, "lines")
  )

p

ggsave("BMCA/plots/fig4/cell_types_paired.pdf", p,
       width = 200, height = 50, units = "mm", dpi = 300)
# -----------------------------------------------------------------------------

# ── Column groups ─────────────────────────────────────────────────────────────
immune_cells   <- c("Myeloid", "B_Plasma", "T_cell", "NK")   # T_cell = T_NK*(1-NK)
tsub_cells     <- c("CD4", "CD8", "GD")
stromal_cells  <- c("Fibroblast", "Pericyte", "Endothelial")
heat_features  <- c("Malignant", immune_cells, tsub_cells, stromal_cells)
carcinoma_sites <- c("NSCLC", "breast", "colorectal")

star <- function(p) dplyr::case_when(
  is.na(p)  ~ "", p < 0.001 ~ "***", p < 0.01 ~ "**", p < 0.05 ~ "*", TRUE ~ "")

# ── Unpaired helper (rows 1 & 2) ──────────────────────────────────────────────
delta_z_unpaired <- function(df, group_label) {
  long <- df %>%
    mutate(T_cell = T_NK * (1 - NK)) %>%
    select(patient, location, seq_tech, all_of(heat_features)) %>%
    pivot_longer(all_of(heat_features), names_to = "cell_type", values_to = "value") %>%
    group_by(cell_type, seq_tech) %>%
    mutate(z = (value - mean(value, na.rm = TRUE)) / sd(value, na.rm = TRUE)) %>%
    ungroup()
  deltas <- long %>%
    group_by(cell_type, location) %>%
    summarise(mean_z = mean(z, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = location, values_from = mean_z) %>%
    mutate(delta = Intracranial - Extracranial) %>% select(cell_type, delta)
  stats <- long %>% filter(!is.na(z)) %>%
    group_by(cell_type) %>% filter(n_distinct(location) == 2) %>%
    wilcox_test(z ~ location) %>% ungroup() %>%
    mutate(p.adj = p.adjust(p, method = "BH"))
  deltas %>% left_join(stats %>% select(cell_type, p.adj), by = "cell_type") %>%
    mutate(group = group_label)
}

# ── Paired helper (row 3): donor×location collapse, paired Wilcoxon ───────────
delta_z_paired <- function(df, group_label) {
  long <- df %>%
    mutate(T_cell = T_NK * (1 - NK)) %>%
    select(Donor, location, seq_tech, all_of(heat_features)) %>%
    pivot_longer(all_of(heat_features), names_to = "cell_type", values_to = "value") %>%
    group_by(cell_type, seq_tech) %>%
    mutate(z = (value - mean(value, na.rm = TRUE)) / sd(value, na.rm = TRUE)) %>%
    ungroup() %>%
    group_by(cell_type, Donor, location) %>%
    summarise(z = mean(z, na.rm = TRUE), .groups = "drop")
  deltas <- long %>%
    group_by(cell_type, location) %>%
    summarise(mean_z = mean(z, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = location, values_from = mean_z) %>%
    mutate(delta = Intracranial - Extracranial) %>% select(cell_type, delta)
  stats <- long %>% filter(!is.na(z)) %>%
    group_by(cell_type, Donor) %>% filter(n_distinct(location) == 2) %>% ungroup() %>%
    arrange(cell_type, Donor, location) %>%
    group_by(cell_type) %>% wilcox_test(z ~ location, paired = TRUE) %>% ungroup() %>%
    mutate(p.adj = p.adjust(p, method = "BH"))
  deltas %>% left_join(stats %>% select(cell_type, p.adj), by = "cell_type") %>%
    mutate(group = group_label)
}

# ── Build rows ────────────────────────────────────────────────────────────────
mel <- celltype_patient_matrix %>%
  filter(primary_site == "melanoma") %>% delta_z_unpaired("Melanoma")
carc_all <- celltype_patient_matrix %>%
  filter(primary_site %in% carcinoma_sites) %>% delta_z_unpaired("Carcinoma (all samples)")
paired_all <- celltype_patient_matrix_donor %>%
  group_by(Donor) %>% filter(n_distinct(location) == 2) %>% ungroup() %>%
  delta_z_paired("Carcinoma (paired samples)")

# ── Assemble ──────────────────────────────────────────────────────────────────
heat <- bind_rows(mel, carc_all, paired_all) %>%
  mutate(
    sig        = star(p.adj),
    cell_type  = factor(cell_type, levels = heat_features),
    group      = factor(group, levels = c("Melanoma",
                                          "Carcinoma (all samples)",
                                          "Carcinoma (paired samples)")),
    row_group  = factor(if_else(group == "Melanoma", "Melanoma", "Carcinoma"),
                        levels = c("Melanoma", "Carcinoma")),
    cell_group = factor(case_when(
      cell_type == "Malignant"     ~ "Malignant",
      cell_type %in% immune_cells  ~ "Immune",
      cell_type %in% tsub_cells    ~ "T subtypes",
      cell_type %in% stromal_cells ~ "Stromal"
    ), levels = c("Malignant", "Immune", "T subtypes", "Stromal"))
  )

disp <- c(Malignant = "Malignant", Myeloid = "Myeloid", B_Plasma = "B/Plasma",
          T_cell = "T", NK = "NK", CD4 = "CD4", CD8 = "CD8", GD = "\u03b3\u03b4",
          Fibroblast = "Fibroblast", Pericyte = "Pericyte", Endothelial = "Endothelial")

lim <- max(abs(heat$delta), na.rm = TRUE)


y_labels <- c(
  "Melanoma" = "Melanoma",
  "Carcinoma (all samples)" = "Carcinoma\n(all samples)",
  "Carcinoma (paired samples)" = "Carcinoma\n(paired samples)"
)

# ── Nature-style heatmap ──────────────────────────────────────────────────────
p_heat <- ggplot(heat, aes(x = cell_type, y = fct_rev(group), fill = delta)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(aes(label = sig), size = 2.6, fontface = "bold",
            colour = "grey10", vjust = 0.8) +
  facet_grid(row_group ~ cell_group, scales = "free", space = "free") +
  scale_x_discrete(labels = disp) +
  scale_y_discrete(labels = y_labels) +
  scale_fill_gradient2(
    low = "#2166AC", mid = "white", high = "#B2182B",
    midpoint = 0, limits = c(-lim, lim),
    breaks = scales::breaks_pretty(3),
    name = expression(Delta*" z-score (BrM \u2212 extracranial relative abundance)")
  ) +
  coord_fixed() +
  labs(x = NULL, y = NULL) +
  theme_classic(base_size = 7, base_family = "Helvetica") +
  theme(
    axis.text.x        = element_text(angle = 45, hjust = 1,colour = "black"),
    axis.text.y = element_text(colour = "black", hjust = 0.5),
    axis.line          = element_blank(),
    axis.ticks         = element_blank(),
    strip.background.x = element_rect(fill = NA, colour = NA),
    strip.background.y = element_blank(),      # remove row strip boxes
    strip.text.x       = element_text(colour = "black",
                                      margin = margin(b = 2)),
    strip.text.y       = element_blank(),
    panel.spacing.x    = unit(2, "pt"),
    panel.spacing.y    = unit(6, "pt"),        # visual gap = melanoma vs carcinoma
    legend.position    = "bottom",
    legend.key.height  = unit(3, "mm"),
    legend.key.width   = unit(10, "mm"),
    legend.margin      = margin(t = 2),
    plot.margin        = margin(4, 6, 4, 4)
  ) +
  guides(fill = guide_colourbar(title.position = "top", title.hjust = 0.5,
                                ticks.colour = "grey30", frame.colour = "grey30",
                                frame.linewidth = 0.3))

p_heat

ggsave("BMCA/plots/fig4/heatmap_delta_z.pdf", p_heat,
       width = 180, height = 80, units = "mm", dpi = 300)
# -----------------------------------------------------------------------------

carcinoma_sites <- c("NSCLC", "breast", "colorectal")

# ── Helpers ───────────────────────────────────────────────────────────────────
ratio_of <- function(df, min_count = 5) {                 # ≥3 in each compartment
  df %>%
    mutate(lymphoid = T_NK + B_Plasma) %>%
    filter(Myeloid >= min_count, lymphoid >= min_count) %>%
    mutate(lr = log2((Myeloid / lymphoid) + 1))
}
center_by_tech <- function(df) {                           # mean-centre lr within seq_tech
  df %>% group_by(seq_tech) %>%
    mutate(lr_c = lr - mean(lr, na.rm = TRUE)) %>% ungroup()
}

# ── Row 1: Melanoma (all samples) ─────────────────────────────────────────────
mel <- celltype_patient_matrix %>%
  filter(primary_site == "melanoma") %>%
  ratio_of() %>% center_by_tech() %>%
  transmute(facet = "Melanoma", location, value = lr_c, Donor = NA_character_)

# ── Row 2: Carcinoma (all samples) ────────────────────────────────────────────
carc <- celltype_patient_matrix %>%
  filter(primary_site %in% carcinoma_sites) %>%
  ratio_of() %>% center_by_tech() %>%
  transmute(facet = "Carcinoma (all samples)", location, value = lr_c, Donor = NA_character_)

# ── Row 3: Carcinoma (paired samples) → donor×location collapse ───────────────
paired <- celltype_patient_matrix_donor %>%
  filter(primary_site %in% carcinoma_sites) %>%
  ratio_of() %>% center_by_tech() %>%
  group_by(Donor, location) %>%
  summarise(value = mean(lr_c, na.rm = TRUE), .groups = "drop") %>%
  group_by(Donor) %>% filter(n_distinct(location) == 2) %>% ungroup() %>%  # complete pairs
  mutate(facet = "Carcinoma (paired samples)")

# ── Combine + relabel BrM / Primary ───────────────────────────────────────────
box_df <- bind_rows(mel, carc, paired) %>%
  mutate(
    location = recode(location, Intracranial = "BrM", Extracranial = "Extracranial"),
    location = factor(location, levels = c("BrM", "Extracranial")),
    facet    = factor(facet, levels = c("Melanoma",
                                        "Carcinoma (all samples)",
                                        "Carcinoma (paired samples)"))
  )

# ── Stats: unpaired (rows 1–2) + paired (row 3), BH across the three ──────────
stat_unpaired <- box_df %>%
  filter(facet %in% c("Melanoma", "Carcinoma (all samples)")) %>%
  group_by(facet) %>% wilcox_test(value ~ location) %>% ungroup()

stat_paired <- box_df %>%
  filter(facet == "Carcinoma (paired samples)") %>%
  arrange(Donor, location) %>%
  wilcox_test(value ~ location, paired = TRUE) %>%
  mutate(facet = "Carcinoma (paired samples)")

rng <- diff(range(box_df$value, na.rm = TRUE))
stats_all <- bind_rows(stat_unpaired, stat_paired) %>%
  mutate(
    p.adj = p.adjust(p, method = "BH"),
    significance = case_when(p.adj < 0.001 ~ "***", p.adj < 0.01 ~ "**",
                             p.adj < 0.05 ~ "*", TRUE ~ "ns"),
    facet = factor(facet, levels = levels(box_df$facet))
  ) %>%
  left_join(box_df %>% group_by(facet) %>%
              summarise(ymax = max(value, na.rm = TRUE), .groups = "drop"),
            by = "facet") %>%
  mutate(y_br = ymax + 0.08 * rng, y_lab = ymax + 0.14 * rng)

pal <- c(BrM = "#B2182B", Extracranial = "#2166AC")

# ── Plot ──────────────────────────────────────────────────────────────────────
p_box <- ggplot(box_df, aes(x = location, y = value)) +
  geom_hline(yintercept = 0, linetype = "22", linewidth = 0.3, colour = "grey75") +
  geom_boxplot(aes(fill = location), width = 0.58, alpha = 0.55,
               outlier.shape = NA, linewidth = 0.4, colour = "grey25") +
  geom_line(data = filter(box_df, facet == "Carcinoma (paired samples)"),
            aes(group = Donor), colour = "grey70", linewidth = 0.3, alpha = 0.6) +
  geom_point(data = filter(box_df, facet == "Carcinoma (paired samples)"),
             aes(colour = location), size = 1, alpha = 0.9) +
  geom_point(data = filter(box_df, facet != "Carcinoma (paired samples)"),
             aes(colour = location),
             position = position_jitter(width = 0.13, height = 0, seed = 1),
             size = 1, alpha = 0.55) +
  geom_segment(data = stats_all, aes(x = 1, xend = 2, y = y_br, yend = y_br),
               inherit.aes = FALSE, linewidth = 0.3, colour = "grey25") +
  geom_text(data = stats_all, aes(x = 1.5, y = y_lab, label = significance),
            inherit.aes = FALSE, colour = "grey10") +
  facet_wrap(~ facet, nrow = 1) +
  scale_fill_manual(values = pal, guide = "none") +
  scale_colour_manual(values = pal, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0.03, 0.12))) +
  labs(x = NULL, y = expression("Myeloid / Lymphoid (log2, mean-centred)")) +
  theme_classic(base_size = 7, base_family = "Helvetica") +
  theme(
    axis.title.y     = element_text(margin = margin(r = 4)),
    axis.line        = element_line(linewidth = 0.4, colour = "black"),
    axis.ticks       = element_line(linewidth = 0.3, colour = "black"),
    strip.background = element_rect(fill = "grey95", colour = NA),
    strip.text       = element_text(face = "bold", colour = "grey10"),
    panel.spacing    = unit(6, "pt"),
    plot.margin      = margin(6, 8, 4, 6)
  )

p_box

ggsave("BMCA/plots/fig4/myeloid_lymph_ratio_box.pdf", p_box,
       width = 120, height = 60, units = "mm", dpi = 300)
# -----------------------------------------------------------------------------

# Get all sheet names
sheets <- excel_sheets("BMCA/csv/all_MPS.xlsx")

# Read each sheet into a named list of tibbles
tibble_list <- lapply(sheets, read_excel, path = "BMCA/csv/all_MPS.xlsx") |>
  setNames(sheets)

myelo_sigs <- tibble_list$Myeloid

tyler_sigs <- read_csv("BMCA/csv/tyler_myeloid_MP.csv")


# Function to calculate Jaccard similarity
jaccard <- function(set1, set2) {
  intersection <- length(intersect(set1, set2))
  union <- length(union(set1, set2))
  intersection / union
}

# Extract signatures as named lists
myelo_list <- as.list(myelo_sigs) %>% map(~ .x[!is.na(.x)])
tyler_list <- as.list(tyler_sigs) %>% map(~ .x[!is.na(.x)])
names(myelo_list)
# --- recode maps defined ONCE, reused for both the data and the level orders ---
myelo_map <- c(
  "sn_Microglia"                   = "Microglia (sn)",
  "sn_Macrophage_1"                = "Macrophage1 (sn)",
  "sn_Macrophage_2"                = "Macrophage2 (sn)",
  "sn_Stress"                      = "Stress (sn)",
  "sn_Scavenger_Immunosuppressive" = "Scavenger IS (sn)",
  "sn_Hypoxia"                     = "Hypoxia (sn)",
  "sn_CC"                          = "Cell Cycle (sn)",
  "sn_DC"                          = "DC (sn)",
  "sn_UA_2"                        = "UA2 (sn)",
  "sn_UA_3"                        = "UA3 (sn)",
  "sn_Monocyte_Neutrophile"        = "Mono/Neut (sn)",
  "sc_Microglia"                   = "Microglia (sc)",
  "sc_Macrophage_1"                = "Macrophage1 (sc)",
  "sc_Monocyte_Neutrophile"        = "Mono/Neut (sc)",
  "sc_Hypoxia"                     = "Hypoxia (sc)",
  "sc_DC"                          = "DC (sc)",
  "sc_INFR_1"                      = "INFR1 (sc)",
  "sc_INFR_2"                      = "INFR2 (sc)",
  "sc_CC"                          = "Cell Cycle (sc)",
  "sc_Scavenger_Immunosuppressive" = "Scavenger IS (sc)",
  "sc_UA_1"                        = "UA1 (sc)",
  "sc_MYC_Mitochondria"            = "MYC/Mito (sc)"
)

tyler_map <- c(
  "Microglia"                    = "Microglia",
  "Macrophage"                   = "Macrophage",
  "Monocyte"                     = "Monocyte",
  "cDC"                          = "DC",
  "Neutrophil"                   = "Neutrophil",
  "Systemic Inflammatory"        = "Systemic Inflam.",
  "Microglial Inflammatory"      = "Microglial Inflam.",
  "Complement Immunosuppressive" = "Complement IS",
  "Scavenger Immunosuppressive"  = "Scavenger IS",
  "Hypoxia"                      = "Hypoxia",
  "IFN Response"                 = "IFN",
  "HS-UPR"                       = "HS/UPR",
  "G2-M"                         = "G2/M",
  "G1-S"                         = "G1/S"
)
# apply recode to the data
jaccard_matrix <- jaccard_matrix %>%
  mutate(
    myelo_sig = recode(myelo_sig, !!!myelo_map),
    tyler_sig = recode(tyler_sig, !!!tyler_map)
  )
# translate the clustering orders through the SAME maps so levels match
row_order_new <- unname(myelo_map[row_order])
col_order_new <- unname(tyler_map[col_order])

jaccard_matrix <- jaccard_matrix %>%
  mutate(
    myelo_sig = factor(myelo_sig, levels = row_order_new),
    tyler_sig = factor(tyler_sig, levels = col_order_new),
    jaccard_capped = pmin(jaccard, 0.3)
  )


jaccard_matrix <- mutate(
  jaccard_matrix,
  myelo_sig = replace_na(myelo_sig, "UA2 (sn)"),
  myelo_sig = factor(myelo_sig)
)

# --- plot ---
p <- ggplot(jaccard_matrix, aes(x = tyler_sig, y = myelo_sig, fill = jaccard_capped)) +
  geom_tile(color = "white", linewidth = 0.5) +
  scale_fill_gradient(
    low = "white", high = "darkgreen",
    limits = c(0, 0.3), name = "Jaccard"
  ) +
  coord_fixed() +
  theme_minimal(base_size = 7) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title  = element_blank()
  )

p
ggsave("BMCA/plots/fig4/tyler_jaccard.pdf", p,
       width = 100, height = 100, units = "mm", dpi = 300)
# -----------------------------------------------------------------------------
# Get all sheet names
sheets <- excel_sheets("BMCA/csv/all_MPS.xlsx")

# Read each sheet into a named list of tibbles
tibble_list <- lapply(sheets, read_excel, path = "BMCA/csv/all_MPS.xlsx") |>
  setNames(sheets)

cd4_sigs <- tibble_list$CD4
cd8_sigs <- tibble_list$CD8

T_sigs <- cbind(cd4_sigs, cd8_sigs)
colnames(T_sigs) <- c(
  "CD4 N/M1 (sn)",
  "CD4 Treg1 (sn)",
  "CD4 Treg2 (sn)",
  "CD4 Cycle (sn)",
  "CD4 N/M2 (sn)",
  "CD4 Tfh (sn)",
  "CD4 UA1 (sn)",
  "CD4 N/M1 (sc)",
  "CD4 Treg (sc)",
  "CD4 IFN (sc)",
  "CD4 Cycle (sc)",
  "CD4 Tfh (sc)",
  "CD4 UA2 (sc)",
  "CD4 HSP (sc)",
  "CD8 Dys1 (sn)",
  "CD8 N/M1 (sn)",
  "CD8 Cycle (sn)",
  "CD8 IFN (sn)",
  "CD8 EM1 (sn)",
  "CD8 Dys2 (sn)",
  "CD8 N/M2 (sn)",
  "CD8 CT (sc)",
  "CD8 IFN (sc)",
  "CD8 UA4 (sc)",
  "CD8 Cycle (sc)",
  "CD8 NK (sc)",
  "CD8 CT1 (sc)",
  "CD8 HSP (sc)",
  "CD8 UA1 (sc)",
  "CD8 MYC (sc)",
  "CD8 CT2 (sc)",
  "CD8 EM2 (sc)",
  "CD8 UA2 (sc)"
)
joyce_sigs <- read_csv("BMCA/csv/joyce_mps.csv")
joyce_sigs_transformed <- joyce_sigs %>%
  filter(p_val_adj < 0.01) %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC), .by_group = TRUE) %>%
  slice_head(n = 50) %>%
  ungroup() %>%
  select(cluster, gene) %>%
  group_by(cluster) %>%
  mutate(row = row_number()) %>%
  pivot_wider(
    names_from = cluster,
    values_from = gene,
    values_fill = NA
  ) %>%
  select(-row) %>%
  ungroup()

# Function to calculate Jaccard similarity
jaccard <- function(set1, set2) {
  intersection <- length(intersect(set1, set2))
  union <- length(union(set1, set2))
  intersection / union
}

# Extract signatures as named lists
T_list <- as.list(T_sigs) %>% map(~ .x[!is.na(.x)])
joyce_list <- as.list(joyce_sigs_transformed) %>% map(~ .x[!is.na(.x)])

# Create matrix of Jaccard similarities
jaccard_matrix <- expand_grid(
  T_sig = names(T_list),
  joyce_sig = names(joyce_list)
) %>%
  mutate(
    jaccard = map2_dbl(
      T_sig, joyce_sig,
      ~ jaccard(T_list[[.x]], joyce_list[[.y]])
    )
  )

# Convert to wide matrix for clustering
jaccard_wide <- jaccard_matrix %>%
  pivot_wider(
    names_from = joyce_sig,
    values_from = jaccard,
    values_fill = 0
  ) %>%
  column_to_rownames("T_sig") %>%
  as.matrix()

# Hierarchical clustering on rows
row_dist <- dist(jaccard_wide, method = "euclidean")
row_hclust <- hclust(row_dist)
row_order <- row_hclust$labels[row_hclust$order]

# Hierarchical clustering on columns
col_dist <- dist(t(jaccard_wide), method = "euclidean")
col_hclust <- hclust(col_dist)
col_order <- col_hclust$labels[col_hclust$order]

# Plot with hierarchical clustering order
p <- jaccard_matrix %>%
  mutate(
    T_sig = factor(T_sig, levels = row_order),
    joyce_sig = factor(joyce_sig, levels = col_order),
    jaccard_capped = pmin(jaccard, 0.3)
  ) %>%
  ggplot(aes(x = joyce_sig, y = T_sig, fill = jaccard_capped)) +
  geom_tile(color = "white", linewidth = 0.5) +
  scale_fill_gradient(
    low = "white",
    high = "darkgreen",
    limits = c(0, 0.3),
    name = "Jaccard"
  ) +
  coord_fixed() +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    axis.title = element_blank()
  )

p

ggsave("BMCA/plots/fig4/joyce_jaccard.pdf", p,
       width = 100, height = 200, units = "mm", dpi = 300)

# -----------------------------------------------------------------------------
ecotypes_tib_proccesed_filtered <- readRDS("BMCA/data/ecotypes_tib_proccesed_filtered.RDS")
ecotypes_tib_proccesed_filtered <- ecotypes_tib_proccesed_filtered %>%
  rename(`Myeloid_Dendritic` = Myeloid_DC)

metadata_states_proccesed_filtered <- readRDS("BMCA/data/metadata_states_proccesed_filtered.RDS")
metadata_states_proccesed_filtered <- metadata_states_proccesed_filtered %>%
  mutate(cancer_type = if_else(primary_site == "melanoma", "Melanoma", "Carcinoma"))


# --- 1. Cell types to loop over ---
cell_types <- c("Myeloid", "B_Plasma", "CD4", "CD8", "GD", "NK",
                "Fibroblast", "Pericyte", "Endothelial", "Malignant")

# --- 2. Join metadata once ---
df_full <- ecotypes_tib_proccesed_filtered %>%
  left_join(metadata_states_proccesed_filtered %>%
              dplyr::select(patient, location, seq_tech, cancer_type),
            by = "patient") %>%
  filter(!is.na(location))

# --- 3. Effect size function (rank-biserial r) ---
compute_effects <- function(data, cols, prefix) {
  cols %>%
    map_dfr(function(col) {
      tmp <- data %>%
        dplyr::select(value = all_of(col), location) %>%
        filter(!is.na(value))

      if (length(unique(tmp$location)) < 2) return(NULL)
      if (sum(tmp$location == "Intracranial") < 3 ||
          sum(tmp$location == "Extracranial") < 3) return(NULL)

      test <- wilcox.test(value ~ location, data = tmp, exact = FALSE)
      n1 <- sum(tmp$location == "Intracranial")
      n2 <- sum(tmp$location == "Extracranial")
      r  <- 1 - (2 * test$statistic) / (n1 * n2)

      tibble(
        state       = str_remove(col, paste0("^", prefix, "_")),
        effect_size = as.numeric(r),
        p_value     = test$p.value
      )
    })
}

for (ct in cell_types) {

  print(ct)
  ct_cols <- names(df_full) %>% str_subset(paste0("^", ct, "_"))
  if (length(ct_cols) == 0) { message("No substate columns for: ", ct); next }

  # ── sc and sn separately, per cancer type ────────────────────────────────────
  per_modality <- list()
  for (ctype in c("Melanoma", "Carcinoma")) {
    df_ct <- df_full %>% filter(cancer_type == ctype)
    for (grp in c("sc", "sn")) {
      res <- tryCatch(
        compute_effects(df_ct %>% filter(seq_tech == grp), ct_cols, ct) %>%
          mutate(seq_tech = grp, cancer_type = ctype),
        error = function(e) NULL)
      if (!is.null(res) && nrow(res)) per_modality <- append(per_modality, list(res))
    }
  }
  if (!length(per_modality)) next

  long <- bind_rows(per_modality) %>%
    mutate(significant = p_value < 0.05)

  # ── Cross-cancer-type masking (per cell, no integration) ─────────────────────
  long <- long %>%
    mutate(
      mask = (cancer_type == "Carcinoma" &
                state %in% c("Neuronal_related_1", "Skin_Pigmentation")) |
        (cancer_type == "Melanoma" &
           state %in% c("secretory_CRC_enriched", "PDAC_related")),
      effect_size = if_else(mask, 0, effect_size),
      significant = if_else(mask, FALSE, significant)
    ) %>%
    select(-mask)

  # ── Combo on x-axis: seq_tech × cancer_type ──────────────────────────────────
  long <- long %>%
    mutate(combo = factor(
      paste(cancer_type, seq_tech, sep = " · "),
      levels = c("Carcinoma · sn", "Carcinoma · sc",
                 "Melanoma · sn",  "Melanoma · sc")
    ))

  # ── Order states by effect in "Carcinoma · sn" (high at top) ─────────────────
  state_order <- long %>%
    filter(combo == "Carcinoma · sn") %>%
    arrange(effect_size) %>%          # ascending -> highest ends up on top
    pull(state) %>%
    as.character()

  # states not present in that combo (dropped by the filter) go to the bottom
  missing_states <- setdiff(as.character(unique(long$state)), state_order)
  state_order <- c(missing_states, state_order)

  long <- long %>%
    mutate(state = factor(state, levels = state_order))

  # ── Plot: one tile per state × (seq_tech × cancer_type) ──────────────────────
  p <- ggplot(long, aes(x = combo, y = state, fill = effect_size)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = ifelse(significant, "*", "")), size = 2, vjust = 0.75) +
    scale_fill_gradient2(
      low = "#4575b4", mid = "white", high = "#d73027",
      midpoint = 0, limits = c(-1, 1),
      name = "BrM vs. Primary\n(rank-biserial r)"
    ) +
    labs(title = ct, x = NULL, y = NULL) +
    theme_minimal(base_size = 7) +
    theme(
      axis.text.x     = element_text(angle = 30, hjust = 1),
      panel.grid      = element_blank(),
      plot.title      = element_text(hjust = 0.5),
      legend.position = "right"
    )

  print(p)
  ggsave(paste0("BMCA/plots/fig4/supp_comper_", ct, ".pdf"), p,
         width = 80, height = 80, units = "mm")
}

# -----------------------------------------------------------------------------
ecotypes_tib_proccesed_filtered <- readRDS("BMCA/data/ecotypes_tib_proccesed_filtered.RDS")
ecotypes_tib_proccesed_filtered <- ecotypes_tib_proccesed_filtered %>%
  rename(`Myeloid_Dendritic` = Myeloid_DC)

metadata_states_proccesed_filtered <- readRDS("BMCA/data/metadata_states_proccesed_filtered.RDS")
metadata_states_proccesed_filtered <- metadata_states_proccesed_filtered %>%
  mutate(cancer_type = if_else(primary_site == "melanoma", "Melanoma", "Carcinoma"))


# --- 1. Cell types to loop over ---
cell_types <- c("Myeloid", "B_Plasma", "CD4", "CD8", "GD", "NK",
                "Fibroblast", "Pericyte", "Endothelial", "Malignant")

# --- 2. Join metadata once ---
df_full <- ecotypes_tib_proccesed_filtered %>%
  left_join(metadata_states_proccesed_filtered %>%
              dplyr::select(patient, location, seq_tech, cancer_type),
            by = "patient") %>%
  filter(!is.na(location))

# --- 3. Effect size function (rank-biserial r) ---
compute_effects <- function(data, cols, prefix) {
  cols %>%
    map_dfr(function(col) {
      tmp <- data %>%
        dplyr::select(value = all_of(col), location) %>%
        filter(!is.na(value))

      if (length(unique(tmp$location)) < 2) return(NULL)
      if (sum(tmp$location == "Intracranial") < 3 ||
          sum(tmp$location == "Extracranial") < 3) return(NULL)

      test <- wilcox.test(value ~ location, data = tmp, exact = FALSE)
      n1 <- sum(tmp$location == "Intracranial")
      n2 <- sum(tmp$location == "Extracranial")
      r  <- 1 - (2 * test$statistic) / (n1 * n2)

      tibble(
        state       = str_remove(col, paste0("^", prefix, "_")),
        effect_size = as.numeric(r),
        p_value     = test$p.value
      )
    })
}

# --- 4. Loop over cell types ---
for (ct in cell_types) {

  print(ct)
  ct_cols <- names(df_full) %>% str_subset(paste0("^", ct, "_"))
  if (length(ct_cols) == 0) { message("No substate columns for: ", ct); next }

  # ── sc and sn separately, per cancer type ────────────────────────────────────
  per_modality <- list()
  for (ctype in c("Melanoma", "Carcinoma")) {
    df_ct <- df_full %>% filter(cancer_type == ctype)
    for (grp in c("sc", "sn")) {
      res <- tryCatch(
        compute_effects(df_ct %>% filter(seq_tech == grp), ct_cols, ct) %>%
          mutate(seq_tech = grp, cancer_type = ctype),
        error = function(e) NULL)
      if (!is.null(res) && nrow(res)) per_modality <- append(per_modality, list(res))
    }
  }
  if (!length(per_modality)) next

  # ── Integrate from long form (robust to a missing modality) ──────────────────
  # mean effect across available modalities; star only if both modalities are
  # present, agree in sign, and at least one is individually significant (p<0.05)
  long <- bind_rows(per_modality)   # state, cancer_type, seq_tech, effect_size, p_value

  combined <- long %>%
    group_by(state, cancer_type) %>%
    summarise(
      n_modalities = n_distinct(seq_tech),
      concordant   = n_modalities == 2 &&
        n_distinct(sign(effect_size[is.finite(effect_size)])) == 1,
      sig_any      = any(p_value < 0.05, na.rm = TRUE),
      effect_size  = mean(effect_size, na.rm = TRUE),   # <- AFTER the sign check
      significant  = concordant & sig_any,
      .groups = "drop"
    )

  # ── Cross-cancer-type masking ────────────────────────────────────────────────
  combined <- combined %>%
    mutate(across(c(effect_size, significant),
                  ~ if_else(cancer_type == "Carcinoma" &
                              state %in% c("Neuronal_related_1", "Skin_Pigmentation"),
                            if (is.numeric(.x)) 0 else FALSE, .x))) %>%
    mutate(across(c(effect_size, significant),
                  ~ if_else(cancer_type == "Melanoma" &
                              state %in% c("secretory_CRC_enriched", "PDAC_related"),
                            if (is.numeric(.x)) 0 else FALSE, .x)))

  # ── Order states by mean combined effect across cancer types ─────────────────
  state_order <- combined %>%
    group_by(state) %>%
    summarise(m = mean(effect_size, na.rm = TRUE), .groups = "drop") %>%
    arrange(m) %>%
    pull(state)

  # ── Add NA row for Malignant Skin_Pigmentation in Carcinoma ──────────────────
  if (ct == "Malignant") {
    combined <- bind_rows(combined, tibble(
      state = "Skin_Pigmentation", cancer_type = "Carcinoma",
      effect_size = 0, significant = FALSE))
    state_order <- union(state_order, "Skin_Pigmentation")
  }

  combined <- combined %>%
    mutate(state       = factor(state, levels = state_order),
           cancer_type = factor(cancer_type, levels = c("Carcinoma", "Melanoma")))

  # ── Plot: one column per cancer type, no sc/sn split ─────────────────────────
  p <- ggplot(combined, aes(x = cancer_type, y = state, fill = effect_size)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = ifelse(significant, "*", "")), size = 5, vjust = 0.75) +
    scale_fill_gradient2(
      low      = "#4575b4",
      mid      = "white",
      high     = "#d73027",
      midpoint = 0,
      limits   = c(-1, 1),
      name     = "BrM vs. Primary\n(mean rank-biserial r)"
    ) +
    labs(title = ct, x = NULL, y = NULL) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x     = element_text(size = 11),
      axis.text.y     = element_text(size = 10),
      panel.grid      = element_blank(),
      strip.text      = element_text(face = "bold", size = 11),
      plot.title      = element_text(hjust = 0.5, face = "bold"),
      legend.position = "right"
    )

  print(p)
}

# -----------------------------------------------------------------------------

cell_types_facet <- c("Malignant", "Myeloid", "CD8", "CD4")

ecotypes_tib_proccesed_filtered_renamed <- ecotypes_tib_proccesed_filtered %>%
  rename(Malignant_MNR1 = Malignant_Neuronal_related_1,
         Malignant_MNR2 = Malignant_Neuronal_related_2,
         Malignant_MNR3 = Malignant_Neuronal_related_3,
         Malignant_MHCI = Malignant_MHC.I,
         `Malignant_MYC targets/respiration` = Malignant_MYC_targets.respiration,
         `Myeloid_MYC/Mitochondria`         = Myeloid_MYC_Mitochondria,
         `Myeloid_Monocyte/Neutrophile`     = Myeloid_Monocyte_Neutrophile,
         `CD8_Cytotoxic MHCII INFG+`        = CD8_CT_MHCII_INFG,
         CD8_Cytotoxic                      = CD8_CT,
         `CD4_Stress/HSP`                   = CD4_Stress_HSP,
         'Malignant_Chromatin Remodeling' = 'Malignant_Chromatin_Remodeling',
         'Malignant_Skin Pigmentation' = `Malignant_Skin_Pigmentation`,
         `Malignant_Cell Cycle` = 'Malignant_CC',
         `Malignant_INFR/MHC` = 'Malignant_INFR_MHC',
         `Malignant_Epithelial Lineage diff. 1` = 'Malignant_secretory_CRC_enriched',
         `Malignant_Epithelial Lineage diff. 2` = 'Malignant_PDAC_related',
         `Malignant_EpiSen/MHC` = 'Malignant_EpiSen_MHC',
         `Malignant_Proteasomal degradation` = 'Malignant_Proteasomal_degradation',
         `Myeloid_Scavenger Immunosuppressive` = 'Myeloid_Scavenger_Immunosuppressive',
         `Myeloid_Cell Cycle` = 'Myeloid_CC',
         `CD8_NK like` = 'CD8_NK_like',
         `CD8_MYC Targets` = 'CD8_MYC_Targets',
         `CD8_Effector-Memory` = 'CD8_Effector_Memory',
         `CD8_Naive-Memory` = 'CD8_Naive/Memory',
         `CD4_Naive-Memory` = 'CD4_Naive/Memory',
         `CD4_Follicular Helper` = 'CD4_Follicular_Helper'
         )

ecotypes_tib_proccesed_filtered_renamed$Malignant_Skin_Pigmentation <- NULL
df_full <- ecotypes_tib_proccesed_filtered_renamed %>%
  left_join(metadata_states_proccesed_filtered %>%
              dplyr::select(patient, location, seq_tech, cancer_type),
            by = "patient") %>%
  filter(!is.na(location))


compute_effects <- function(data, cols, prefix) {
  cols %>% map_dfr(function(col) {
    tmp <- data %>% dplyr::select(value = all_of(col), location) %>% filter(!is.na(value))
    if (length(unique(tmp$location)) < 2) return(NULL)
    if (sum(tmp$location == "Intracranial") < 3 ||
        sum(tmp$location == "Extracranial") < 3) return(NULL)
    test <- wilcox.test(value ~ location, data = tmp, exact = FALSE)
    n1 <- sum(tmp$location == "Intracranial"); n2 <- sum(tmp$location == "Extracranial")
    r  <- 1 - (2 * test$statistic) / (n1 * n2)
    tibble(state = str_remove(col, paste0("^", prefix, "_")),
           effect_size = as.numeric(r), p_value = test$p.value)
  })
}

build_combined <- function(ct) {
  ct_cols <- names(df_full) %>% str_subset(paste0("^", ct, "_"))
  if (length(ct_cols) == 0) return(NULL)

  per_modality <- list()
  for (ctype in c("Melanoma", "Carcinoma")) {
    df_ct <- df_full %>% filter(cancer_type == ctype)
    for (grp in c("sc", "sn")) {
      res <- tryCatch(
        compute_effects(df_ct %>% filter(seq_tech == grp), ct_cols, ct) %>%
          mutate(seq_tech = grp, cancer_type = ctype),
        error = function(e) NULL)
      if (!is.null(res) && nrow(res)) per_modality <- append(per_modality, list(res))
    }
  }
  if (!length(per_modality)) return(NULL)

  combined <- bind_rows(per_modality) %>%
    group_by(state, cancer_type) %>%
    summarise(
      n_modalities = n_distinct(seq_tech),
      concordant   = n_modalities == 2 &&
        n_distinct(sign(effect_size[is.finite(effect_size)])) == 1,
      sig_any      = any(p_value < 0.05, na.rm = TRUE),
      effect_size  = mean(effect_size, na.rm = TRUE),   # <- AFTER the sign check
      significant  = concordant & sig_any,
      .groups = "drop"
    ) %>%
    mutate(across(c(effect_size, significant),
                  ~ if_else(cancer_type == "Carcinoma" &
                              state %in% c("Neuronal_related_1", "Skin_Pigmentation"),
                            if (is.numeric(.x)) 0 else FALSE, .x))) %>%
    mutate(across(c(effect_size, significant),
                  ~ if_else(cancer_type == "Melanoma" &
                              state %in% c("Epithelial Lineage diff. 1",
                                           "Epithelial Lineage diff. 2",
                                           "EpiSen/MHC"),
                            if (is.numeric(.x)) 0 else FALSE, .x)))

  if (ct == "Malignant") {
    combined <- bind_rows(
      combined,
      tibble(
        state       = "Skin Pigmentation",
        cancer_type = "Carcinoma",
        effect_size = 0,
        significant = FALSE
      ),
      tibble(
        state       = c("Epithelial Lineage diff. 1",
                        "Epithelial Lineage diff. 2",
                        "EpiSen/MHC"),
        cancer_type = "Melanoma",
        effect_size = 0,
        significant = FALSE
      )
    )
  }

  combined %>% mutate(cell_type = ct)
}

make_panel <- function(ct) {
  d <- build_combined(ct)
  if (is.null(d)) return(NULL)

  state_order <- d %>%
    filter(cancer_type == "Carcinoma") %>%
    arrange(effect_size) %>%
    pull(state)

  # any states missing from Carcinoma (e.g. melanoma-only) go at the bottom
  state_order <- union(state_order, unique(d$state))

  d <- d %>%
    mutate(state       = factor(state, levels = state_order),
           cancer_type = factor(cancer_type, levels = c("Carcinoma", "Melanoma")))

  ggplot(d, aes(x = cancer_type, y = state, fill = effect_size)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = ifelse(significant, "*", "")), size = 4.5, vjust = 0.75) +
    scale_fill_gradient2(
      low = "#4575b4", mid = "white", high = "#d73027",
      midpoint = 0, limits = c(-1, 1),
      name = "BrM vs. Primary\n(mean rank-biserial r)"
    ) +
    labs(title = ct, x = NULL, y = NULL) +
    theme_minimal(base_size = 7) +
    theme(
      axis.text.x  = element_text(angle = 45, hjust = 1),
      panel.grid   = element_blank(),
      plot.title   = element_text(hjust = 0.5, face = "bold", size = 0),
      legend.position = "bottom"
    )
}

panels <- map(cell_types_facet, make_panel) %>% compact()

# ── Combine side by side, one shared legend at the bottom ─────────────────────
combined_plot <- wrap_plots(panels, nrow = 1) +
  plot_layout(guides = "collect") &
  theme(
    legend.position   = "none",
    legend.title      = element_text(size = 7),
    legend.text       = element_text(size = 7),
    legend.key.size   = unit(2, "mm")
  ) &
  guides(fill = guide_colourbar(title.position = "top", title.hjust = 0.5,
                                barwidth = unit(3, "mm"), barheight = unit(30, "mm")))
combined_plot

ggsave("BMCA/plots/fig4/cellstate_heatmaps_combined.pdf", combined_plot,
       width = 180, height = 120, units = "mm", dpi = 300)
# -----------------------------------------------------------------------------
# --- 1. Define the antigen presentation states per cell type ---
antigen_states <- list(
  Malignant   = c("Malignant_INFR_MHC", "Malignant_MHC.I", "Malignant_EpiSen_MHC"),
  Myeloid     = c("Myeloid_Dendritic"),
  B_Plasma    = c("B_Plasma_Antigen_presentation"),
  GD          = c("GD_MHCII_Dysfunctional"),
  CD8         = c("CD8_CT_MHCII_INFG"),
  Endothelial = c("Endothelial_HEV")
)

# --- 2. Join metadata ---
df_full <- ecotypes_tib_proccesed_filtered %>%
  left_join(metadata_states_proccesed_filtered %>% dplyr::select(patient, location),
            by = "patient") %>%
  filter(!is.na(location))

# --- 3. Compute mean difference + Wilcoxon p-value per cell type ---
results <- imap_dfr(antigen_states, function(cols, cell_type) {

  tmp <- df_full %>%
    dplyr::select(patient, location, all_of(cols)) %>%
    rowwise() %>%
    mutate(mean_val = mean(c_across(all_of(cols)), na.rm = TRUE)) %>%
    ungroup() %>%
    filter(!is.na(mean_val))

  if (length(unique(tmp$location)) < 2) return(NULL)

  # Mean difference
  means <- tmp %>%
    group_by(location) %>%
    summarise(mean_val = mean(mean_val, na.rm = TRUE), .groups = "drop")

  intra <- means %>% filter(location == "Intracranial") %>% pull(mean_val)
  extra <- means %>% filter(location == "Extracranial") %>% pull(mean_val)

  if (length(intra) == 0 || length(extra) == 0) return(NULL)

  # Wilcoxon test
  test <- wilcox.test(mean_val ~ location, data = tmp, exact = FALSE)

  tibble(
    cell_type = cell_type,
    delta     = intra - extra,
    p_value   = test$p.value
  )
}) %>%
  mutate(
    p_adj = p.adjust(p_value, method = "BH"),
    stars = case_when(
      p_adj < 0.001 ~ "***",
      p_adj < 0.01  ~ "**",
      p_adj < 0.05  ~ "*",
      TRUE          ~ ""
    ),
    # Order by delta low to high
    cell_type = fct_reorder(cell_type, delta)
  )

results$cell_type <- recode(results$cell_type,
                            "GD" = "Gamma Delta T Cells",
                            "CD8" = "CD8 T Cells")
# --- 4. Plot ---
ggplot(results, aes(x = cell_type, y = delta)) +
  geom_bar(stat = "identity", fill = "steelblue", color = "black", linewidth = 0.4) +
  geom_hline(yintercept = 0, linewidth = 0.5) +
  geom_text(
    aes(
      label = stars,
      y     = delta + ifelse(delta >= 0, 0.3, -0.3)
    ),
    size = 5, vjust = 1
  ) +
  labs(
    title = "",
    x     = "Cell Type",
    y     = ""
  ) +
  theme_minimal(base_size = 16) +
  theme(
    axis.text.x        = element_text(angle = 45, hjust = 1, size = 16),
    panel.grid.major.x = element_blank(),
    plot.title         = element_text(hjust = 0.5, face = "bold")
  )

# -----------------------------------------------------------------------------

# --- 1. Define states per theme ---
theme_states <- list(
  `Antigen Presentation` = list(
    Malignant   = c("Malignant_INFR_MHC", "Malignant_MHC.I", "Malignant_EpiSen_MHC"),
    Myeloid     = c("Myeloid_Dendritic"),
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
    x     = "Cell Type",
    y     = "Mean Intracranial − Mean Extracranial"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x        = element_text(angle = 45, hjust = 1),
    panel.grid.major.x = element_blank(),
    plot.title         = element_text(hjust = 0.5, face = "bold"),
    strip.text         = element_text(face = "bold", size = 12)
  )

p <- ggplot(results, aes(x = cell_type, y = delta, fill = theme)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.4) +
  geom_hline(yintercept = 0, linewidth = 0.5) +
  geom_text(
    aes(
      label = stars,
      y     = delta + ifelse(delta >= 0, 0.3, -1)
    ),
    size = 2, vjust = 0.5, inherit.aes = TRUE
  ) +
  facet_wrap(~ theme, scales = "free_x") +
  scale_fill_manual(values = c(
    "Antigen Presentation" = "steelblue",
    "Hypoxia"              = "#C64D40"   # firebrick, the red equivalent
  )) +
  labs(
    title = "",
    x     = "Transcriptional state per cell type",
    y     = "Relative abundance (BrM − Extracranial)"
  ) +
  theme_minimal(base_size = 7) +
  theme(
    axis.text.x        = element_text(angle = 45, hjust = 1),
    panel.grid.major.x = element_blank(),
    plot.title         = element_text(hjust = 0.5),
    legend.position    = "none"        # facet labels already show the theme
  )
p

ggsave("BMCA/plots/fig4/themes_enrich.pdf", p,
       width = 80, height = 60, units = "mm", dpi = 300)
# -----------------------------------------------------------------------------

md_tib <- readRDS("BMCA/data/md_tib_Myeloid.RDS")
score_tib <- readRDS("BMCA/data/score_tib_Myeloid.RDS")
relevant_states <- c("sn_Scavenger_Immunosuppressive", "sn_DC",
                     "sn_Macrophage_2",                "sn_Monocyte_Neutrophile",
                     "sc_Macrophage_1",                "sc_Monocyte_Neutrophile",
                     "sc_DC", "sc_Microglia", "sn_Macrophage_1",
                     "sn_Hypoxia", "sc_Hypoxia", "sc_Scavenger_Immunosuppressive",
                     "sn_Microglia")

md_tib_new <- md_tib %>%
  filter(annotation %in% relevant_states) %>%
  mutate(
    Hypoxia                     = pmax(sc_Hypoxia,                       sn_Hypoxia,                    na.rm = TRUE),
    DC                          = pmax(sc_DC,                            sn_DC,                         na.rm = TRUE),
    Macrophage                  = pmax(sc_Macrophage_1,                  sn_Macrophage_1, sn_Macrophage_2, na.rm = TRUE),
    Microglia                   = pmax(sc_Microglia,                     sn_Microglia,                  na.rm = TRUE),
    Monocyte_Neutrophile        = pmax(sc_Monocyte_Neutrophile,          sn_Monocyte_Neutrophile,       na.rm = TRUE),
    Scavenger_Immunosuppressive = pmax(sc_Scavenger_Immunosuppressive,   sn_Scavenger_Immunosuppressive, na.rm = TRUE)
  )  %>%
  dplyr::select(1:7, annotation, Hypoxia, DC, Macrophage, Microglia,
                Monocyte_Neutrophile, Scavenger_Immunosuppressive)  # keep metadata + annotation


md_tib_new <- md_tib_new %>%
  mutate(annotation = str_remove(annotation, "^sc_|^sn_"))
md_tib_new <- md_tib_new %>%
  mutate(annotation = str_remove(annotation, "^sc_|^sn_") %>%
           str_remove("_\\d+$"))

# -----------------------------------------------------------------------------

make_density_plot <- function(target, data = md_tib_new,
                              loc_filter = NULL, facet_loc = TRUE) {
  d <- data %>% filter(annotation %in% c(target, "Scavenger_Immunosuppressive"))
  if (!is.null(loc_filter)) d <- d %>% filter(location %in% loc_filter)

  facet <- if (facet_loc) facet_grid(seq_tech ~ location)
  else           facet_grid(rows = vars(seq_tech))

  ggplot(d, aes(x = Scavenger_Immunosuppressive, y = .data[[target]])) +
    geom_density_2d(color = "steelblue", linewidth = 0.5) +
    facet +
    theme_classic() +
    theme(
      strip.background = element_rect(fill = "grey90", color = NA),
      strip.text   = element_text(face = "bold", size = 12),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
      plot.title   = element_text(face = "bold", size = 13, hjust = 0.5)
    ) +
    labs(title = gsub("_", "/", target),
         x = "Scavenger Immunosuppressive", y = target)
}

# Hypoxia removed
shared_targets <- c("DC", "Macrophage", "Monocyte_Neutrophile")
plots <- lapply(shared_targets, make_density_plot)

# Microglia keeps its Intracranial-only, seq_tech-faceted layout
p_micro <- make_density_plot("Microglia", loc_filter = "Intracranial",
                             facet_loc = FALSE)

# All panels in one figure, side by side
combined <- wrap_plots(c(plots, list(p_micro)), nrow = 1,
                       widths = c(2, 2, 2, 1))   # Microglia is half-width (1 col vs 2)
combined
# -----------------------------------------------------------------------------

metadata_p <- metadata %>%
  dplyr::select(patient, primary_site, study, location, seq_tech) %>%
  distinct()

T_NK_pb <- readRDS("BMCA/data/UMI_pb_T_NK.RDS")


ifng_tbl <- imap_dfr(T_NK_pb, function(mat, dataset) {
  if (!"IFNG" %in% rownames(mat)) {
    warning(dataset, ": no IFNG row, skipping")
    return(NULL)
  }
  tibble(
    patient  = colnames(mat),
    IFNG_exp = mat["IFNG", ]
  )
})
ifng_tbl$IFNG_exp <- log2(ifng_tbl$IFNG_exp + 1)
ifng_tbl <- ifng_tbl %>% left_join(metadata_p, by = "patient")


plot_df <- ifng_tbl %>%
  filter(primary_site %in% c("breast", "colorectal", "NSCLC", "melanoma")) %>%
  mutate(
    cancer_type = if_else(primary_site == "melanoma", "Melanoma", "Carcinoma"),
    cancer_type = factor(cancer_type, levels = c("Carcinoma", "Melanoma")),
    location    = factor(location, levels = c("Extracranial", "Intracranial"))
  )

loc_cols <- c("Extracranial" = "#4575b4", "Intracranial" = "#d73027")

ggplot(plot_df, aes(x = location, y = IFNG_exp)) +
  geom_boxplot(aes(fill = location), outlier.shape = NA, width = 0.6,
               alpha = 0.9, colour = "grey25", linewidth = 0.4) +
  geom_jitter(aes(colour = location), width = 0.15, size = 0.7,
              alpha = 0.5, show.legend = FALSE) +
  stat_compare_means(method = "wilcox.test",
                     comparisons = list(c("Extracranial", "Intracranial")),
                     label = "p.signif", size = 3.2, tip.length = 0.01) +
  facet_nested(~ cancer_type + seq_tech, nest_line = element_line(colour = "black")) +
  scale_fill_manual(values = loc_cols) +
  scale_colour_manual(values = loc_cols) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.12))) +
  labs(x = NULL,
       y = expression(italic("IFNG")~"(T NK Pseudo-bulk expression, log2)"),
       fill = "Location") +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x      = element_text(angle = 45, hjust = 1, colour = "black"),
    axis.text.y      = element_text(colour = "black"),
    axis.line        = element_line(colour = "black", linewidth = 0.4),
    axis.ticks       = element_line(colour = "black"),
    strip.background  = element_rect(fill = "grey95", colour = NA),
    strip.text        = element_text(face = "bold"),
    legend.position   = "none",
    panel.spacing     = unit(0.4, "lines")
  )

plot_df <- plot_df %>%
  group_by(cancer_type, seq_tech, location) %>%        # literal: per box
  # group_by(cancer_type, seq_tech) %>%                # alt: per facet, keeps Extra vs Intra
  mutate(IFNG_centered = IFNG_exp - mean(IFNG_exp, na.rm = TRUE)) %>%
  ungroup()

plot_df$seq_tech[plot_df$seq_tech == "sn"] <- "Single-nucleus"
plot_df$seq_tech[plot_df$seq_tech == "sc"] <- "Single-cell"

ggplot(plot_df, aes(x = location, y = IFNG_centered)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60", linewidth = 0.3) +
  geom_boxplot(aes(fill = location), outlier.shape = NA, width = 0.6,
               alpha = 0.9, colour = "grey25", linewidth = 0.4) +
  geom_jitter(aes(colour = location), width = 0.15, size = 1,
              alpha = 0.5, show.legend = FALSE) +
  stat_compare_means(method = "wilcox.test",
                     comparisons = list(c("Extracranial", "Intracranial")),
                     label = "p.signif", size = 2, tip.length = 0.01) +
  facet_nested(~ cancer_type + seq_tech, nest_line = element_line(colour = "black")) +
  scale_fill_manual(values = loc_cols) +
  scale_colour_manual(values = loc_cols) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.12))) +
  labs(x = NULL,
       y = expression(italic("IFNG")~"(T/NK pseudo-bulk, log2, mean-centered)"),
       fill = "Location") +
  theme_classic(base_size = 7) +
  theme(
    axis.text.x      = element_text(angle = 45, hjust = 1, colour = "black"),
    axis.text.y      = element_text(colour = "black"),
    axis.line        = element_line(colour = "black", linewidth = 0.4),
    axis.ticks       = element_line(colour = "black"),
    strip.background  = element_rect(fill = "grey95", colour = NA),
    legend.position   = "none",
    panel.spacing     = unit(0.4, "lines")
  )
ggsave("BMCA/plots/fig4/IFNG_TNK_centered.pdf", width = 3, height = 3,
       dpi = 300)


# ── Integrate sc + sn: mean-centre within seq_tech (per cancer type) ──────────
plot_df_int <- plot_df %>%
  group_by(cancer_type, seq_tech) %>%
  mutate(IFNG_centered_tech = IFNG_exp - mean(IFNG_exp, na.rm = TRUE)) %>%  # remove tech baseline
  ungroup()

loc_cols <- c("Extracranial" = "#4575b4", "Intracranial" = "#d73027")

p_ifng_int <- ggplot(plot_df_int, aes(x = location, y = IFNG_centered_tech)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60", linewidth = 0.3) +
  geom_boxplot(aes(fill = location), outlier.shape = NA, width = 0.6,
               alpha = 0.9, colour = "grey25", linewidth = 0.4) +
  geom_jitter(aes(colour = location), width = 0.15, size = 0.7,
              alpha = 0.5, show.legend = FALSE) +
  stat_compare_means(method = "wilcox.test",
                     comparisons = list(c("Extracranial", "Intracranial")),
                     label = "p.signif", size = 3.2, tip.length = 0.01) +
  facet_wrap(~ cancer_type) +                       # ← seq_tech integrated out
  scale_fill_manual(values = loc_cols) +
  scale_colour_manual(values = loc_cols) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.12))) +
  labs(x = NULL,
       y = expression(atop(italic("IFNG")~"(T/NK pseudo-bulk, log2,",
                           "seq. tech. mean-centered)")),
       fill = "Location") +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x      = element_text(angle = 45, hjust = 1, colour = "black"),
    axis.text.y      = element_text(colour = "black"),
    axis.line        = element_line(colour = "black", linewidth = 0.4),
    axis.ticks       = element_line(colour = "black"),
    strip.background = element_rect(fill = "grey95", colour = NA),
    strip.text       = element_text(face = "bold"),
    legend.position  = "none",
    panel.spacing    = unit(0.4, "lines")
  )

p_ifng_int

ggsave("BMCA/plots/fig4/IFNG_TNK_integrated.pdf", p_ifng_int, width = 4.5, height = 4)

