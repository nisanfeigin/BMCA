# =============================================================================
# Figure 3 - Malignant meta-programs and malignant neural-related (MNR) programs
# =============================================================================
#
# Jaccard similarity between malignant meta-programs derived from single-cell
# and single-nucleus data, comparison against published 3CA meta-programs,
# expression heatmaps of the selected malignant MPs, MNR program scoring,
# and CNA-correlation analyses for the MNR programs.
#
# Part of: A brain metastasis cell atlas reveals multicellular adaptation to
#          the neural microenvironment (Feigin et al.)
#
# Run from the project root (see README.md). All paths below are relative to it.
#
# Inputs:
#   BMCA/csv/all_MPS.xlsx
#   BMCA/csv/meta_programs_2025.xlsx
#   BMCA/csv/mps_<cell_type>.csv   # path built at run time
#   BMCA/csv/mps_MNRs.csv
#   BMCA/data/centered_cpm_<cohort>.RDS   # path built at run time
#   BMCA/data/ecotypes_tib_proccesed_filtered.RDS
#   BMCA/data/list_of_cna_cor_sig_tibs_MNRs.RDS
#   BMCA/data/M_cor_all_for_violin.RDS
#   BMCA/data/Malignant_MP_sc_2.rds
#   BMCA/data/Malignant_MP_sn_2.rds
#   BMCA/data/md_tib_<cell_type>.RDS   # path built at run time
#   BMCA/data/metadata_all_studies.rds
#   BMCA/data/metadata_states_proccesed_filtered.RDS
#   BMCA/data/mps_selected_mal_sc.RDS
#   BMCA/data/mps_selected_mal_sn.RDS
#   BMCA/data/neuron_gene_list.RDS
#   BMCA/data/programs_all_for_fig_3.RDS
#   BMCA/data/score_tib_<cell_type>.RDS   # path built at run time
#   BMCA/data/score_tib_mnrs.RDS
#   BMCA/data/UMI_pb_<cell_type>.RDS   # path built at run time
#   code/inferGEP/data/hotmap.rds
#
# Outputs:
#   BMCA/data/  (6 files)
#   BMCA/plots/fig3/  (29 files)
#   BMCA/plots/figure_3/  (6 files)
#
# Original filename: Figure_3_2.R
# =============================================================================

# ---- Packages ---------------------------------------------------------------
pkgs <- c(
  "biomaRt", "cowplot", "data.table", "dendextend", "doParallel", "dplyr",
  "enrichR", "forcats", "foreach", "ggdendro", "ggplot2", "ggrastr",
  "ggrepel", "Matrix", "NMF", "parallel", "patchwork", "pheatmap",
  "purrr", "RColorBrewer", "readxl", "reshape2", "scales", "scalop",
  "stringr", "tibble", "tidyr", "tidyverse", "viridis"
)
invisible(lapply(pkgs, library, character.only = TRUE))

# ---- Project paths ----------------------------------------------------------
source("R/config.R")

# =============================================================================
metadata <- as_tibble(readRDS("BMCA/data/metadata_all_studies.rds"))


# -----------------------------------------------------------------------------
programs_sc <- readRDS("BMCA/data/Malignant_MP_sc_2.rds")

programs_sc <- as_tibble(programs_sc)
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Function to calculate Jaccard similarity
jaccard_similarity <- function(genes1, genes2) {
  length(intersect(genes1, genes2)) / length(union(genes1, genes2))
}

# Get the number of programs_sc (columns) in the tibble
num_programs_sc <- ncol(programs_sc)

# Initialize a matrix to store the Jaccard similarities
jaccard_matrix <- matrix(0, nrow = num_programs_sc, ncol = num_programs_sc)

# Calculate Jaccard similarities for each pair of programs_sc
for (i in 1:num_programs_sc) {
  for (j in 1:num_programs_sc) {
    # Get the sets of genes for program i and program j
    set1 <- programs_sc[[i]]  # Extract the gene set from the ith column
    set2 <- programs_sc[[j]]  # Extract the gene set from the jth column

    # Compute the Jaccard similarity between these two gene sets
    jaccard_matrix[i, j] <- jaccard_similarity(set1, set2)
  }
}

# jaccard_matrix now contains the Jaccard similarity for each pair of programs_sc
jaccard_matrix[1:5,1:5]

# Convert the matrix to a data frame for better visualization if needed
jaccard_df <- as.data.frame(jaccard_matrix)

# Print the Jaccard similarity matrix
#print(jaccard_df)
color_palette <- readRDS("code/inferGEP/data/hotmap.rds")

M_cor_sc <- jaccard_matrix
M_hc  <- hclust(as.dist(1-M_cor_sc), method="average")             ### create an object of class hclust which is a list that describes the tree produced by the clustering process
M_hc  <- reorder(as.dendrogram(M_hc), colMeans(M_cor_sc))          ### the leaves of the dendrogram are reordered so as to be in an order as consistent as possible with the weights (which are here given as colMeans(M_cor_sc))
M_cor_sc <- M_cor_sc[order.dendrogram(M_hc), order.dendrogram(M_hc)]

M_melt <- reshape2::melt(M_cor_sc)

# Define the values for the gradient from -0.4 to 0.4
col_co <- 0.2
gradient_values <- seq(-col_co, col_co, length.out = 11)

# Create the bins for the gradient values
breaks <- seq(-col_co, col_co, length.out = length(color_palette) + 1)

# Bin the values in the dataset
M_melt <- M_melt %>%
  mutate(value_to_plot = ifelse(value > col_co, col_co, value),
         binned_value = cut(value_to_plot, breaks = breaks, include.lowest = TRUE, labels = FALSE))

# Handle NA values in binned_value by assigning them to the closest bin
M_melt$binned_value[is.na(M_melt$binned_value) & M_melt$value_to_plot <= -col_co] <- 1
M_melt$binned_value[is.na(M_melt$binned_value) & M_melt$value_to_plot >= col_co] <- length(color_palette)

# Convert binned_value to factor for consistent coloring
M_melt$binned_value <- factor(M_melt$binned_value, levels = 1:length(color_palette))

# Reverse the levels for the legend
M_melt$binned_value <- factor(M_melt$binned_value, levels = rev(levels(M_melt$binned_value)))

# -----------------------------------------------------------------------------
# Create custom labels for the legend
legend_names <- c("20%", rep("", length(color_palette) - 2), "-0.40")

# Reverse the color palette for the legend
color_palette_reversed <- rev(color_palette)

# Create the plot with binned colors and customized legend labels
cor_plot <- ggplot(M_melt, aes(x = Var1, y = Var2, fill = binned_value, color = binned_value)) +
  geom_tile() +
  scale_fill_manual(values = setNames(color_palette_reversed, levels(M_melt$binned_value)),
                    labels = legend_names, name = "Jaccard\nsimilarity") +
  scale_color_manual(values = setNames(color_palette_reversed, levels(M_melt$binned_value)),
                     labels = legend_names, name = "Jaccard\nsimilarity") +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
  )

cor_plot

# -----------------------------------------------------------------------------

programs_sn <- readRDS("BMCA/data/Malignant_MP_sn_2.rds")
programs_sn <- as_tibble(programs_sn)

# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Function to calculate Jaccard similarity
jaccard_similarity <- function(genes1, genes2) {
  length(intersect(genes1, genes2)) / length(union(genes1, genes2))
}

# Get the number of programs_sn (columns) in the tibble
num_programs_sn <- ncol(programs_sn)

# Initialize a matrix to store the Jaccard similarities
jaccard_matrix <- matrix(0, nrow = num_programs_sn, ncol = num_programs_sn)

# Calculate Jaccard similarities for each pair of programs_sn
for (i in 1:num_programs_sn) {
  for (j in 1:num_programs_sn) {
    # Get the sets of genes for program i and program j
    set1 <- programs_sn[[i]]  # Extract the gene set from the ith column
    set2 <- programs_sn[[j]]  # Extract the gene set from the jth column

    # Compute the Jaccard similarity between these two gene sets
    jaccard_matrix[i, j] <- jaccard_similarity(set1, set2)
  }
}

# jaccard_matrix now contains the Jaccard similarity for each pair of programs_sn
jaccard_matrix[1:5,1:5]

# Convert the matrix to a data frame for better visualization if needed
jaccard_df <- as.data.frame(jaccard_matrix)

# Print the Jaccard similarity matrix
#print(jaccard_df)
color_palette <- readRDS("code/inferGEP/data/hotmap.rds")

M_cor_sn <- jaccard_matrix
M_hc  <- hclust(as.dist(1-M_cor_sn), method="average")             ### create an object of class hclust which is a list that describes the tree produced by the clustering process
M_hc  <- reorder(as.dendrogram(M_hc), colMeans(M_cor_sn))          ### the leaves of the dendrogram are reordered so as to be in an order as consistent as possible with the weights (which are here given as colMeans(M_cor_sn))
M_cor_sn <- M_cor_sn[order.dendrogram(M_hc), order.dendrogram(M_hc)]

M_melt <- reshape2::melt(M_cor_sn)

# Define the values for the gradient from -0.4 to 0.4
col_co <- 0.2
gradient_values <- seq(-col_co, col_co, length.out = 11)

# Create the bins for the gradient values
breaks <- seq(-col_co, col_co, length.out = length(color_palette) + 1)

# Bin the values in the dataset
M_melt <- M_melt %>%
  mutate(value_to_plot = ifelse(value > col_co, col_co, value),
         binned_value = cut(value_to_plot, breaks = breaks, include.lowest = TRUE, labels = FALSE))

# Handle NA values in binned_value by assigning them to the closest bin
M_melt$binned_value[is.na(M_melt$binned_value) & M_melt$value_to_plot <= -col_co] <- 1
M_melt$binned_value[is.na(M_melt$binned_value) & M_melt$value_to_plot >= col_co] <- length(color_palette)

# Convert binned_value to factor for consistent coloring
M_melt$binned_value <- factor(M_melt$binned_value, levels = 1:length(color_palette))

# Reverse the levels for the legend
M_melt$binned_value <- factor(M_melt$binned_value, levels = rev(levels(M_melt$binned_value)))

# -----------------------------------------------------------------------------
# Create custom labels for the legend
legend_names <- c("20%", rep("", length(color_palette) - 2), "-0.40")

# Reverse the color palette for the legend
color_palette_reversed <- rev(color_palette)

# Create the plot with binned colors and customized legend labels
cor_plot <- ggplot(M_melt, aes(x = Var1, y = Var2, fill = binned_value, color = binned_value)) +
  geom_tile() +
  scale_fill_manual(values = setNames(color_palette_reversed, levels(M_melt$binned_value)),
                    labels = legend_names, name = "Jaccard\nsimilarity") +
  scale_color_manual(values = setNames(color_palette_reversed, levels(M_melt$binned_value)),
                     labels = legend_names, name = "Jaccard\nsimilarity") +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
  )

cor_plot

# -----------------------------------------------------------------------------
programs_sc<- programs_sc[,colnames(M_cor_sc)]
programs_sn<- programs_sn[,colnames(M_cor_sn)]

programs_all <- as_tibble(cbind(programs_sc, programs_sn))
saveRDS(programs_all, "BMCA/data/programs_all_for_fig_3.RDS")
programs_all <- readRDS("BMCA/data/programs_all_for_fig_3.RDS")

programs <- programs_all

# Function to calculate Jaccard similarity
jaccard_similarity <- function(genes1, genes2) {
  length(intersect(genes1, genes2)) / length(union(genes1, genes2))
}

# Get the number of programs (columns) in the tibble
num_programs <- ncol(programs)

# Initialize a matrix to store the Jaccard similarities
jaccard_matrix <- matrix(0, nrow = num_programs, ncol = num_programs)

# Calculate Jaccard similarities for each pair of programs
for (i in 1:num_programs) {
  for (j in 1:num_programs) {
    # Get the sets of genes for program i and program j
    set1 <- programs[[i]]  # Extract the gene set from the ith column
    set2 <- programs[[j]]  # Extract the gene set from the jth column

    # Compute the Jaccard similarity between these two gene sets
    jaccard_matrix[i, j] <- jaccard_similarity(set1, set2)
  }
}

# jaccard_matrix now contains the Jaccard similarity for each pair of programs
jaccard_matrix[1:5,1:5]

# Convert the matrix to a data frame for better visualization if needed
jaccard_df <- as.data.frame(jaccard_matrix)

# Print the Jaccard similarity matrix
#print(jaccard_df)
color_palette <- readRDS("code/inferGEP/data/hotmap.rds")

M_cor_all <- jaccard_matrix
M_cor_all <- M_cor_all[colnames(programs), colnames(programs)]

M_melt <- reshape2::melt(M_cor_all)

# Define the values for the gradient from -0.4 to 0.4
col_co <- 0.2
gradient_values <- seq(-col_co, col_co, length.out = 11)

# Create the bins for the gradient values
breaks <- seq(-col_co, col_co, length.out = length(color_palette) + 1)

# Bin the values in the dataset
M_melt <- M_melt %>%
  mutate(value_to_plot = ifelse(value > col_co, col_co, value),
         binned_value = cut(value_to_plot, breaks = breaks, include.lowest = TRUE, labels = FALSE))

# Handle NA values in binned_value by assigning them to the closest bin
M_melt$binned_value[is.na(M_melt$binned_value) & M_melt$value_to_plot <= -col_co] <- 1
M_melt$binned_value[is.na(M_melt$binned_value) & M_melt$value_to_plot >= col_co] <- length(color_palette)

# Convert binned_value to factor for consistent coloring
M_melt$binned_value <- factor(M_melt$binned_value, levels = 1:length(color_palette))

# Reverse the levels for the legend
M_melt$binned_value <- factor(M_melt$binned_value, levels = rev(levels(M_melt$binned_value)))

# -----------------------------------------------------------------------------
# Create custom labels for the legend
legend_names <- c("20%", rep("", length(color_palette) - 2), "-0.40")

# Reverse the color palette for the legend
color_palette_reversed <- rev(color_palette)
# Create the plot with binned colors and customized legend labels
cor_plot <- ggplot(M_melt, aes(x = Var1, y = Var2, fill = binned_value, color = binned_value)) +
  geom_tile() +
  scale_fill_manual(values = setNames(color_palette_reversed, levels(M_melt$binned_value)),
                    labels = legend_names, name = "Jaccard\nsimilarity") +
  scale_color_manual(values = setNames(color_palette_reversed, levels(M_melt$binned_value)),
                     labels = legend_names, name = "Jaccard\nsimilarity") +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
  )


cor_plot <- ggplot(M_melt, aes(x = Var1, y = Var2, fill = binned_value, color = binned_value)) +
  rasterise(geom_tile(), dpi = 300) +
  scale_fill_manual(values = setNames(color_palette_reversed, levels(M_melt$binned_value)),
                    labels = legend_names, name = "Jaccard\nsimilarity") +
  scale_color_manual(values = setNames(color_palette_reversed, levels(M_melt$binned_value)),
                     labels = legend_names, name = "Jaccard\nsimilarity") +
  theme_minimal(base_size = 5) +
  theme(axis.text.x  = element_blank(),
        axis.text.y  = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "right")
cor_plot

ggsave(
  filename = "BMCA/plots/fig3/jaccarad_sn_sc_mal.pdf",
  plot = cor_plot,
  width = 100,
  height = 100,
  units = "mm",
  dpi = 300   # lower than the default 300
)


# -----------------------------------------------------------------------------
programs <- programs_all
# Function to calculate Jaccard similarity
jaccard_similarity <- function(genes1, genes2) {
  length(intersect(genes1, genes2)) / length(union(genes1, genes2))
}

# Get the number of programs (columns) in the tibble
num_programs <- ncol(programs)

# Initialize a matrix to store the Jaccard similarities
jaccard_matrix <- matrix(0, nrow = num_programs, ncol = num_programs)

# Calculate Jaccard similarities for each pair of programs
for (i in 1:num_programs) {
  for (j in 1:num_programs) {
    # Get the sets of genes for program i and program j
    set1 <- programs[[i]]  # Extract the gene set from the ith column
    set2 <- programs[[j]]  # Extract the gene set from the jth column

    # Compute the Jaccard similarity between these two gene sets
    jaccard_matrix[i, j] <- jaccard_similarity(set1, set2)
  }
}

# jaccard_matrix now contains the Jaccard similarity for each pair of programs
jaccard_matrix[1:5,1:5]

# Convert the matrix to a data frame for better visualization if needed
jaccard_df <- as.data.frame(jaccard_matrix)

# Print the Jaccard similarity matrix
#print(jaccard_df)
color_palette <- readRDS("code/inferGEP/data/hotmap.rds")

M_cor <- jaccard_matrix
M_hc  <- hclust(as.dist(1-M_cor), method="average")             ### create an object of class hclust which is a list that describes the tree produced by the clustering process
M_hc  <- reorder(as.dendrogram(M_hc), colMeans(M_cor))          ### the leaves of the dendrogram are reordered so as to be in an order as consistent as possible with the weights (which are here given as colMeans(M_cor))
M_cor <- M_cor[order.dendrogram(M_hc), order.dendrogram(M_hc)]

M_melt <- reshape2::melt(M_cor)

# Define the values for the gradient from -0.4 to 0.4
col_co <- 0.2
gradient_values <- seq(-col_co, col_co, length.out = 11)

# Create the bins for the gradient values
breaks <- seq(-col_co, col_co, length.out = length(color_palette) + 1)

# Bin the values in the dataset
M_melt <- M_melt %>%
  mutate(value_to_plot = ifelse(value > col_co, col_co, value),
         binned_value = cut(value_to_plot, breaks = breaks, include.lowest = TRUE, labels = FALSE))

# Handle NA values in binned_value by assigning them to the closest bin
M_melt$binned_value[is.na(M_melt$binned_value) & M_melt$value_to_plot <= -col_co] <- 1
M_melt$binned_value[is.na(M_melt$binned_value) & M_melt$value_to_plot >= col_co] <- length(color_palette)

# Convert binned_value to factor for consistent coloring
M_melt$binned_value <- factor(M_melt$binned_value, levels = 1:length(color_palette))

# Reverse the levels for the legend
M_melt$binned_value <- factor(M_melt$binned_value, levels = rev(levels(M_melt$binned_value)))

# -----------------------------------------------------------------------------
# Create custom labels for the legend
legend_names <- c("20%", rep("", length(color_palette) - 2), "-0.40")

# Reverse the color palette for the legend
color_palette_reversed <- rev(color_palette)

# Create the plot with binned colors and customized legend labels
cor_plot <- ggplot(M_melt, aes(x = Var1, y = Var2, fill = binned_value, color = binned_value)) +
  geom_tile() +
  scale_fill_manual(values = setNames(color_palette_reversed, levels(M_melt$binned_value)),
                    labels = legend_names, name = "Jaccard\nsimilarity") +
  scale_color_manual(values = setNames(color_palette_reversed, levels(M_melt$binned_value)),
                     labels = legend_names, name = "Jaccard\nsimilarity") +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
  )

cor_plot

# -----------------------------------------------------------------------------
ca <- 100

cut_avg <- cutree(M_hc, k = ca)

samples <- names(cut_avg)
group_by_ass <- tibble(prog = samples, group = cut_avg)
top_genes_tibble <- as_tibble(programs)

condition <- (table(group_by_ass$group) > 20) == TRUE

# Get the names of the columns that meet the condition
true_cols <- names(condition)[condition]
group_by_ass <- group_by_ass %>%
  filter(group %in% true_cols)
order <- colnames(M_cor)

# Use match to create an index that orders group_by_ass by the order in M_cor
group_by_ass_ordered <- group_by_ass[match(order, group_by_ass$prog), ]
# Assuming group_by_ass_ordered is your data frame
# Check the result
group_by_ass_ordered <- group_by_ass_ordered[!is.na(group_by_ass_ordered$group), ]

print(group_by_ass_ordered)

saveRDS(group_by_ass_ordered, "BMCA/data/mps_selected_mal_sn.RDS")

 str(M_melt)
# Keep the same order as the heatmap
annotation_df <- data.frame(
  Var1 = levels(M_melt$Var1),
  program_type = ifelse(grepl("_sc$", levels(M_melt$Var1)), "_sc", "_sn")
)


# Reorder factor to match heatmap x-axis
annotation_df$Var1 <- factor(annotation_df$Var1, levels = levels(M_melt$Var1))
annotation_df$program_type <- factor(annotation_df$program_type)

annotation_plot <- ggplot(annotation_df, aes(x = Var1, y = 1, fill = program_type)) +
  geom_tile() +
  scale_fill_manual(
    name = "Program type",
    values = c("_sc" = "#778873", "_sn" = "#8ca9ff"),
    labels = c("_sc" = "Single-cell", "_sn" = "Single-nucleus")
  ) +
  theme_minimal() +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "right",
    panel.grid = element_blank()
  ) +
  coord_cartesian(expand = FALSE)
annotation_plot
ggsave(
  filename = "BMCA/plots/figure_3/annotation_bar_sn_sc.pdf",
  plot = annotation_plot,
  width = 6,
  height = 2,
  units = "in"
)

group_by_ass_ordered_sn <- readRDS("BMCA/data/mps_selected_mal_sn.RDS") %>%
  pull(prog)

group_by_ass_ordered_sc <- readRDS("BMCA/data/mps_selected_mal_sc.RDS") %>%
  pull(prog)
selected <- group_by_ass_ordered$prog
# Keep the same order as the heatmap
# Remove the underscore and everything following it
#M_melt$Var1 <- sub("_[^_]*$", "", M_melt$Var1)
selected <- sub("_[^_]*$", "", selected)


# Create the dataframe using unique values
annotation_df <- data.frame(
  Var1 = unique(M_melt$Var1),
  selected = ifelse(unique(M_melt$Var1) %in% selected, "selected", "not_selected")
)
# Reorder factor to match heatmap x-axis
progs <- sub("_[^_]*$", "", colnames(programs_all))
progs <- factor(progs)

annotation_df$Var1 <- factor(annotation_df$Var1, levels = progs)
annotation_df$selected <- factor(annotation_df$selected)

annotation_plot <- ggplot(annotation_df, aes(x = Var1, y = 1, fill = selected)) +
  geom_tile() +
  scale_fill_manual(
    name = "selected into MP",
    values = c("selected" = "#000000", "not_selected" = "#FFFFFF")
  ) +
  theme_minimal() +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "right",
    panel.grid = element_blank()
  ) +
  coord_cartesian(expand = FALSE)
annotation_plot

ggplot(annotation_df, aes(x = Var1, fill = selected)) +
  geom_bar(width = 1) + # Each bar represents one observation
  scale_fill_manual(values = c("selected" = "black", "not_selected" = "white"))

ggsave(
  filename = "BMCA/plots/figure_3/annotation_bar_sn_sc_selected.pdf",
  plot = annotation_plot,
  width = 6,
  height = 2,
  units = "in"
)
# -----------------------------------------------------------------------------

saveRDS(M_cor_all, "BMCA/data/M_cor_all_for_violin.RDS")
M_cor_all <- readRDS("BMCA/data/M_cor_all_for_violin.RDS")

# Melt the correlation / Jaccard matrix
M_melt <- reshape2::melt(M_cor_all)

# Add program_type columns
M_melt <- M_melt %>%
  mutate(program_type_1 = ifelse(grepl("_sc$", Var1), "_sc", "_sn"),
         program_type_2 = ifelse(grepl("_sc$", Var2), "_sc", "_sn"))

# Function to get top N Jaccard values per Var1
get_top_jaccard <- function(df, row_type, col_type, top_n = 3) {
  df %>%
    filter(program_type_1 == row_type,
           program_type_2 == col_type,
           Var1 != Var2) %>%
    group_by(Var1) %>%
    # Use reframe instead of summarise
    reframe(top_jaccard = sort(value, decreasing = TRUE)[1:min(top_n, n())]) %>%
    pull(top_jaccard)
}

tn <- 5
# Get top tn Jaccard values for each comparison
sn_sn <- get_top_jaccard(M_melt, "_sn", "_sn", top_n = tn)
sc_sc <- get_top_jaccard(M_melt, "_sc", "_sc", top_n = tn)
sn_sc <- get_top_jaccard(M_melt, "_sn", "_sc", top_n = tn)
sc_sn <- get_top_jaccard(M_melt, "_sc", "_sn", top_n = tn)

# Pool cross values (_sn vs _sc and _sc vs _sn)
cross_pool <- c(sn_sc, sc_sn)

# Sample random values from the full matrix (excluding self-comparisons)
set.seed(67)  # reproducibility
random_vals <- M_melt %>%
  filter(Var1 != Var2) %>%
  group_by(Var1) %>%
  group_modify(~ {
    n_take <- min(tn, nrow(.x))
    .x %>% slice_sample(n = n_take)
  }) %>%
  ungroup() %>%
  pull(value)

# Combine all values into one data frame
boxplot_df <- data.frame(
  jaccard = c(sn_sn, sc_sc, cross_pool, random_vals),
  comparison = rep(c("SN vs. SN", "SC vs. SC", "SN vs. SC", "Randomly sampled similarities"),
                   times = c(length(sn_sn), length(sc_sc), length(cross_pool), length(random_vals)))
)
# Define colorblind-friendly palette
colors <- c("SN vs. SN" = "#8ca9ff",    # blue
            "SC vs. SC" = "#939f90",    # green
            "SN vs. SC" = "#7570b3",
            "Randomly sampled similarities" = "#E97F4A")       # orange

# Ensure desired order on x-axis
boxplot_df$comparison <- factor(
  boxplot_df$comparison,
  levels = c("SN vs. SN", "SN vs. SC", "SC vs. SC", "Randomly sampled similarities")
)

# Pairwise comparisons to display
comparisons <- list(
  c("SN vs. SN", "SN vs. SC"),
  c("SN vs. SN", "SC vs. SC"),
  c("SN vs. SC", "SC vs. SC"),
  c("SN vs. SN", "Randomly sampled similarities")
)

ggplot(boxplot_df, aes(x = comparison, y = jaccard, fill = comparison)) +
  geom_boxplot(
    alpha = 1,
    width = 0.6,
    color = "black"
  ) +

  stat_compare_means(
    comparisons = comparisons,
    method = "wilcox.test",
    label = "p.signif",   # *, **, ***
    size = 5,
    tip.length = 0.01
  ) +

  scale_fill_manual(values = colors) +

  labs(
    x = "",
    y = "Jaccard similarity",
    title = "",
    caption = ""
  ) +

  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    axis.text = element_text(size = 12),
    panel.grid.major.x = element_blank()
  )
# Define your desired order
target_order <- c("SC vs. SC", "SN vs. SN", "SN vs. SC", "Randomly sampled similarities")

# Reorder the factor in your dataframe
boxplot_df$comparison <- factor(boxplot_df$comparison, levels = target_order)


# --- optional but recommended: register a sans-serif font ---
# Nature wants Helvetica/Arial. Easiest reliable route:
# library(showtext); font_add("Arial", "/path/to/Arial.ttf"); showtext_auto()
# or just leave family unset and set it at export via the Cairo/PDF device.

# n per group for annotation
n_labels <- boxplot_df |>
  dplyr::group_by(comparison) |>
  dplyr::summarize(n = dplyr::n(), .groups = "drop")

p <- ggplot(boxplot_df, aes(x = comparison, y = jaccard, fill = comparison)) +
  geom_violin(
    alpha    = 0.7,
    trim     = TRUE,          # keep: Jaccard is bounded [0,1], don't let the KDE spill past the data
    scale    = "width",
    linewidth = 0.3,          # thin outline for print
    colour   = "grey20"
  ) +
  scale_fill_manual(values = colors) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.25),
    expand = expansion(mult = c(0.02, 0.05))
  ) +
  scale_x_discrete(labels = c(
    "SC vs. SC"                     = "Single-cell vs.\nsingle-cell",
    "SN vs. SN"                     = "Single-nuclei vs.\nsingle-nuclei",
    "SN vs. SC"                     = "Single-cell vs.\nsingle-nuclei",
    "Randomly sampled similarities" = "Random pairs"
  )) +
  geom_text(
    data = n_labels,
    aes(x = comparison, y = -0.02, label = paste0("n=", n)),
    inherit.aes = FALSE,
    size = 2.1, colour = "grey30", vjust = 1
  ) +
  labs(x = NULL, y = "Jaccard similarity") +
  theme_classic(base_size = 7) +          # 7pt base ≈ Nature body text
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    axis.line       = element_line(linewidth = 0.3),
    axis.ticks      = element_line(linewidth = 0.3)
  )
p

                    # embeds fonts, true vector
ggsave(
  filename = "BMCA/plots/fig3/jaccard_comparison_boxplot.pdf",
  plot = p,   # saves the most recent ggplot
  width = 88,
  height = 70,
  units = "mm",
  dpi = 300
)

# -----------------------------------------------------------------------------
glimpse(metadata)

metadata$cell_type[cell_type == "Oilgodendrocyte"] <- "Oligodendrocyte"
cts <- unique(metadata$cell_type)
cts <- c(cts, "CD8", "CD4", "GD")
cts <- cts[!cts %in% c(
  "Doublet", NA, "Mesenchimal", "T_NK", "Other",
  "Astrocyte", "Neuron", "Epithelial", "Malignant"
)]
list_of_jaccared_plots <- list()
list_of_scores_cor_plots <- list()
ct <- "Malignant"
for (ct in cts) {
  print(ct)
  sigs <- as_tibble(read.csv(paste0("BMCA/csv/mps_", ct, ".csv")))
  sigs$sn_Hypoxia_2 <- NULL
  sigs$sc_Respiration_2 <- NULL
  sigs$sc_Stress_2 <- NULL
  sigs$sc_lncRNA_MT <- NULL

  clean_names <- c(
    "sn_CC"                       = "sn_CC",
    "sn_Hypoxia"                  = "sn_Hypoxia",
    "sn_Neuronal_related_1"       = "sn_MNR1",
    "sn_INFR"                     = "sn_INFR",
    "sn_UA1"                      = "sn_UA",
    "sn_Neuronal_related_2"       = "sn_MNR2",
    "sn_Neuronal_related_3"       = "sn_MNR3",
    "sn_secretory_CRC_enriched"   = "sn_Epi. diff. 2",
    "sn_Stress"                   = "sn_Stress",
    "sn_Skin_Pigmentation"        = "sn_Skin Pigmentation",
    "sc_Stress"                   = "sc_Stress",
    "sc_MHC.I"                    = "sc_MHCI",
    "sc_CC"                       = "sc_CC",
    "sc_INFR_MHC"                 = "sc_INFR MHC",
    "sc_EpiSen_MHC"               = "sc_EpiSen/MHC",
    "sc_Hypoxia"                  = "sc_Hypoxia",
    "sc_Proteasomal_degradation"  = "sc_Proteasomal deg.",
    "sc_Skin_Pigmentation"        = "sc_Skin Pigmentation",
    "sc_MYC_targets.respiration"  = "sc_MYC targets/respiration",
    "sc_Chromatin_Remodeling"     = "sc_Chromatin remodeling",
    "sc_Respiration_1"            = "sc_Respiration",
    "sc_PDAC_related"             = "sc_Epi. diff. 1",
    "sc_secretory_CRC_enriched"   = "sc_Epi. diff. 2",
    "sc_Neuronal_related_2"       = "sc_MNR2"
  )

  colnames(sigs) <- clean_names[colnames(sigs)]
  # 1. Select only the relevant columns
  sn_cols <- sigs %>% dplyr::select(starts_with("sn_"))
  sc_cols <- sigs %>% dplyr::select(starts_with("sc_"))

  # Convert columns to gene-set lists
  sc_sets <- sc_cols %>%
    as.list() %>%
    map(~ unique(na.omit(.x)))

  sn_sets <- sn_cols %>%
    as.list() %>%
    map(~ unique(na.omit(.x)))

  # Jaccard function
  jaccard <- function(a, b) {
    length(intersect(a, b)) / length(union(a, b))
  }

  # Compute Jaccard matrix
  jaccard_mat <- map_dfr(
    sc_sets,
    function(sc) map_dbl(sn_sets, function(sn) jaccard(sc, sn)),
    .id = "sc_signature"
  ) %>%
    column_to_rownames("sc_signature") %>%
    as.matrix()

  # Remove prefixes for nicer labels
  rownames(jaccard_mat) <- sub("^sc_", "", rownames(jaccard_mat))
  colnames(jaccard_mat) <- sub("^sn_", "", colnames(jaccard_mat))

  # --- Order rows and columns by their max values ---
  row_order <- order(apply(jaccard_mat, 1, max), decreasing = TRUE)
  col_order <- order(apply(jaccard_mat, 2, max), decreasing = TRUE)

  jaccard_mat <- jaccard_mat[row_order, col_order]

  # White → dark green palette
  green_palette <- colorRampPalette(c("white", "darkgreen"))(100)

  # Format numbers (2 decimals)
  jaccard_labels <- matrix(
    sprintf("%.2f", jaccard_mat),
    nrow = nrow(jaccard_mat),
    ncol = ncol(jaccard_mat),
    dimnames = dimnames(jaccard_mat)
  )

  # Heatmap with numbers
  pheatmap(
    jaccard_mat,
    color = green_palette,
    display_numbers = jaccard_labels,
    number_color = "black",
    fontsize_number = 8,
    border_color = NA,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    main = "Jaccard index: sc vs sn signatures (ordered by max)",
    angle_col = 45
  )

  # --- Order rows and columns by their max values ---
  colnames(jaccard_mat)
  cols <- colnames(jaccard_mat)
  rownames(jaccard_mat)
  jaccard_mat <- jaccard_mat[, cols]


  # Heatmap with numbers
  pheatmap(
    jaccard_mat,
    color = green_palette,
    display_numbers = jaccard_labels,
    number_color = "black",
    fontsize_number = 8,
    border_color = NA,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    main = "",
    angle_col = 45
  )

  row_order <- rownames(jaccard_mat)
  col_order <- colnames(jaccard_mat)
  col_order <- c( "INFR", "Hypoxia",
                  "Epi. diff. 2",
                  "UA",       "CC",     "Stress",
                  "Skin Pigmentation",      "MNR2",
                  "MNR3",
                  "MNR1"  )
  col_order <- rev(col_order)


  # Convert matrix to tidy data frame
  heatmap_df <- as.data.frame(jaccard_mat) %>%
    rownames_to_column("Row") %>%
    pivot_longer(
      cols = -Row,
      names_to = "Column",
      values_to = "Jaccard"
    ) %>%
    mutate(
      Row = factor(Row, levels = row_order),
      Column = factor(Column, levels = col_order)
    )


  # Create heatmap
  p <- ggplot(heatmap_df, aes(x = Row, y = Column, fill = Jaccard)) +
    geom_tile(color = NA) +

    geom_text(
      aes(label = sprintf("%.2f", Jaccard)),
      size  = 0,
      color = "black"
    ) +

    scale_fill_gradientn(
      colors = green_palette,
      name   = "Jaccard\nsimilarity",
      guide  = guide_colorbar(
        frame.colour    = "black",
        frame.linewidth = 0.3,
        ticks.colour    = "black",
        ticks.linewidth = 0.3,
        title.position  = "top",
        barwidth        = unit(0.3, "cm"),
        barheight       = unit(2.5, "cm")
      )
    ) +

    coord_fixed() +

    labs(
      x     = "",
      y     = "",
      title = paste0(ct)
    ) +

    theme_minimal(base_size = 7) +
    theme(
      axis.text.x     = element_text(angle = 45, hjust = 1, vjust = 1),
      plot.title      = element_text(face = "bold", hjust = 0.5, size = 0),
      panel.grid      = element_blank(),
      legend.position = "right",

      # shrink the legend
      legend.text     = element_text(size = 5),
      legend.title    = element_text(size = 5),
      legend.spacing  = unit(0.2, "cm"),
      legend.margin   = margin(1, 1, 1, 1)
    )
  p <- p + theme(panel.border = element_rect(color = "black", fill = NA, size = 1))
  p

  ggsave(
    filename = "BMCA/plots/fig3/mal_mps_jacc.pdf",
    plot = p,
    width = 160,
    height = 80,
    units = "mm",
    dpi = 300
  )


  legend <- get_plot_component(
    p +
      theme(legend.position = "right") +
      guides(fill = guide_colourbar(
        title.position = "right",   # title below the bar
        title.hjust    = 0.5,        # centered
        barwidth       = unit(0.5, "cm"),
        barheight      = unit(4, "cm")
      )),
    "guide-box-right",
    return_all = TRUE
  )
  ggsave("BMCA/plots/fig3/jacc_heatmap_legend.pdf", legend,
         width = 20, height = 10, units = "mm",
         dpi = 300)

  list_of_jaccared_plots[[ct]] <- p
# -----------------------------------------------------------------------------

  score_tib <- readRDS(paste0("BMCA/data/score_tib_", ct, ".RDS"))
  md_tib <- readRDS(paste0("BMCA/data/md_tib_", ct, ".RDS"))

  # Keep only rows that **do NOT** contain "oubl"
  md_tib <- md_tib %>%
    filter(!grepl("oubl", annotation))# 1. Select only the relevant columns

  score_tib <- score_tib %>%
    dplyr::select(-contains("oubl"))
  colnames(score_tib)

  sn_cols <- score_tib %>% dplyr::select(starts_with("sn_"))
  sc_cols <- score_tib %>% dplyr::select(starts_with("sc_"))

  # 2. Compute correlation matrix
  cor_matrix <- cor(as.matrix(sn_cols), as.matrix(sc_cols), use = "pairwise.complete.obs", method = "pearson")

  cor_matrix[is.na(cor_matrix)] <- 0

  # --- Sort rows/columns by descending max similarity ---
  row_order <- order(apply(cor_matrix, 1, max), decreasing = TRUE)
  col_order <- order(apply(cor_matrix, 2, max), decreasing = TRUE)

  # Apply sort
  cor_matrix_sorted <- cor_matrix[row_order, col_order]

  # Remove "sn_" from rownames
  rownames(cor_matrix_sorted) <- sub("^sn_", "", rownames(cor_matrix_sorted))

  # Remove "sc_" from colnames
  colnames(cor_matrix_sorted) <- sub("^sc_", "", colnames(cor_matrix_sorted))


  # 1. Define your color palette
  custom_diverging <- colorRampPalette(c("#6BAED6", "white", "#F8766D"))
  #f58442
  custom_diverging <- colorRampPalette(c("#7A3E9D", "white", "#E68A00"))


  # 2. Define breaks such that 0 is at the center
  # Get min and max from your correlation matrix
  cor_min <- min(cor_matrix_sorted, na.rm = TRUE)
  cor_max <- max(cor_matrix_sorted, na.rm = TRUE)
  abs_max <- max(abs(cor_min), abs(cor_max))  # for symmetric scale around 0

  # Create breaks from -abs_max to +abs_max
  breaks <- seq(-abs_max, abs_max, length.out = 101)

  # 3. Generate color palette with same number of steps
  colors <- custom_diverging(length(breaks) - 1)

  # 4. Plot with pheatmap
  pheatmap(cor_matrix_sorted,
           cluster_rows = FALSE,
           cluster_cols = FALSE,
           display_numbers = TRUE,
           number_format = "%.2f",
           fontsize_number = 8,
           fontsize =12,
           color = colors,
           breaks = breaks,
           main = "Correlation between MPs scores based on cells scored for all MPs",
           border_color = NA,
           angle_col = 45)


  # 1️⃣ Convert matrix to long format for ggplot
  cor_df <- as.data.frame(cor_matrix_sorted) %>%
    rownames_to_column(var = "sn") %>%
    pivot_longer(-sn, names_to = "sc", values_to = "correlation")

  # 2️⃣ Optional: remove "sn_" and "sc_" if not already done
  cor_df$sn <- sub("^sn_", "", cor_df$sn)
  cor_df$sc <- sub("^sc_", "", cor_df$sc)

  # 3️⃣ Set factor levels to preserve your sorted order
  cor_df$sn <- factor(cor_df$sn, levels = rownames(cor_matrix_sorted))
  cor_df$sc <- factor(cor_df$sc, levels = colnames(cor_matrix_sorted))

  # 4️⃣ Define color palette (diverging, centered at 0)
  abs_max <- max(abs(cor_df$correlation))
  colors <- colorRampPalette(c("#7A3E9D", "white", "#E68A00"))(100)
  cor_df$sc <- factor(cor_df$sc, levels = rev(unique(cor_df$sc)))

  p <- ggplot(cor_df, aes(x = sn, y = sc, fill = correlation)) +  # swap x and y
    geom_tile(color = NA) +
    scale_fill_gradient2(
      low = "#7A3E9D",
      mid = "white",
      high = "#E68A00",
      midpoint = 0,
      limits = c(-abs_max, abs_max),
      name = "Pearson\nCorrelation"
    ) +
    theme_minimal(base_size = 5) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid = element_blank(),
      ,
      legend.position = "right"
    ) +
    labs(
      title = "",
      x = "",     # now x-axis
      y = ""   # now y-axis
    )

  p
  ggsave("BMCA/plots/fig3/exp_mal_mps_heatmap_.pdf", p,
         width = 100, height = 100, units = "mm",
         dpi = 300)
  p <- p + theme(panel.border = element_rect(color = "black", fill = "black", linewidth = 1))
  p

  list_of_scores_cor_plots[[ct]] <- p
}
p <- list_of_scores_cor_plots[[ct]]

combined_plot_scores <- wrap_plots(list_of_scores_cor_plots, ncol = 3)

combined_plot_scores

ggsave(
  filename = "BMCA/plots/figure_3/combined_plot_scores.pdf",
  plot = combined_plot_scores,
  width = 12,
  height = 12,
  units = "in"
)

combined_plot_Jacc <- wrap_plots(list_of_jaccared_plots, ncol = 4)

combined_plot_Jacc

ggsave(
  filename = "BMCA/plots/figure_3/combined_plot_Jacc.pdf",
  plot = combined_plot_Jacc,
  width = 12,
  height = 12,
  units = "in"
)

# -----------------------------------------------------------------------------
#panel C

mal_cells_assig <- metadata %>%
  filter(cell_type == "Malignant") %>%
  dplyr::select(CellID, seq_tech, state)

mal_cells_assig <- mal_cells_assig %>%
  mutate(
    state = ifelse(is.na(state), "Not Assigned", state),
    state = ifelse(state == "Malignant_Hypoxia_2", "Malignant_Hypoxia", state),
    state = ifelse(state == "Malignant_Respiration_2", "Malignant_Respiration", state),
    state = ifelse(state == "Malignant_Respiration_1", "Malignant_Respiration", state),
    state = ifelse(state == "Malignant_Stress_2", "Malignant_Stress", state),
    state = str_remove(state, "^Malignant_")
  )


# ── Palette ──────────────────────────────────────────────────────────────────
state_colors <- c(
  "CC"                     = "#4DBBD5",
  "MYC_targets.respiration"= "#8491B4",
  "Stress"                 = "#E64B35",
  "Hypoxia"                = "#F39B7F",
  "INFR_MHC"               = "#3C5488",
  "INFR"                   = "#5B84B1",
  "MHC.I"                  = "#85B7EB",
  "EpiSen_MHC"             = "#B5D4F4",
  "Respiration"            = "#00A087",
  "lncRNA_MT"              = "#58B09C",
  "Proteasomal_degradation"= "#A8D8CE",
  "Chromatin_Remodeling"   = "#7B4F9E",
  "EMT"                    = "#F2A900",
  "Skin_Pigmentation"      = "#D5A76B",
  "secretory_CRC_enriched" = "#C0A882",
  "PDAC_related"           = "#8B6D45",
  "Neuronal_related_1"     = "#D2006E",
  "Neuronal_related_2"     = "#FF4DAD",
  "Neuronal_related_3"     = "#FFB3DC",
  "Not Assigned"           = "#B2B2B2",
  "UA1"                    = "#D3D3D3"
)

# ── Data ─────────────────────────────────────────────────────────────────────
pie_data <- mal_cells_assig %>%
  count(seq_tech, state) %>%
  group_by(seq_tech) %>%
  mutate(
    pct   = n / sum(n),
    label = ifelse(pct >= 0.04, paste0(round(pct * 100, 1), "%"), ""),
    state = factor(state, levels = names(state_colors))
  ) %>%
  ungroup()

# ── Plot ──────────────────────────────────────────────────────────────────────
n_tech <- n_distinct(pie_data$seq_tech)

tech_labels <- c(
  "sn" = "snRNA-seq",
  "sc" = "scRNA-seq"
)

p <- ggplot(pie_data, aes(x = "", y = pct, fill = state)) +
  geom_col(width = 1, color = "white", linewidth = 0.25) +
  geom_text(
    aes(label = label),
    position = position_stack(vjust = 0.5),
    size     = 2.5,
    color    = "white",
    fontface = "bold"
  ) +
  coord_polar(theta = "y", clip = "off") +
  facet_wrap(~ seq_tech, nrow = 1, labeller = labeller(seq_tech = tech_labels)) +
  scale_fill_manual(
    values = state_colors,
    name   = NULL,
    guide  = guide_legend(
      ncol         = 1,
      byrow        = TRUE,
      keywidth     = unit(0.3, "cm"),
      keyheight    = unit(0.3, "cm"),
      override.aes = list(color = NA)
    )
  ) +
  theme_void(base_size = 5, base_family = "Helvetica") +
  theme(
    strip.text       = element_text(size = 5,
                                    margin = margin(b = 4)),
    legend.text      = element_text(size = 5),
    legend.position  = "bottom",
  )
p

ggsave(
  filename = "BMCA/plots/fig3/mps_dist.pdf",
  plot = p,
  width = 160,
  height = 300,
  units = "mm"
)
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------

# Get all sheet names
sheets <- excel_sheets("BMCA/csv/all_MPS.xlsx")

# Read each sheet into a named list of tibbles
tibble_list <- lapply(sheets, read_excel, path = "BMCA/csv/all_MPS.xlsx") |>
  setNames(sheets)


colnames_list <- lapply(tibble_list, colnames)
colnames_list_raw <- lapply(colnames_list, function(x) {
  x[!grepl("oubl|[Tt]ech", x)]
})
colnames_list <- lapply(colnames_list_raw, function(x) {
  gsub("(?<!UA)(?<!related)_\\d+$", "", x, perl = TRUE)
})
colnames_list <- lapply(colnames_list, unique)


mps_tib <- lapply(names(colnames_list), function(cell_type) {
  x <- colnames_list[[cell_type]]

  # Astrocyte and Neuron are all sn
  if (cell_type %in% c("Astrocyte", "Neuron")) {
    sn_only <- length(x)
    shared  <- 0
    sc_only <- 0
  } else {
    sn <- gsub("^sn_", "", x[grepl("^sn_", x)])
    sc <- gsub("^sc_", "", x[grepl("^sc_", x)])

    shared  <- length(intersect(sn, sc))
    sn_only <- length(sn) - shared
    sc_only <- length(sc) - shared
  }

  data.frame(cell_type = cell_type, sn_only = sn_only, shared = shared, sc_only = sc_only)
}) |> do.call(rbind, args = _)


mps_tib %>%
  mutate(total = sn_only + shared + sc_only) %>%
  arrange(total) %>%
  mutate(cell_type = factor(cell_type, levels = cell_type)) %>%
  pivot_longer(cols = c(sn_only, shared, sc_only), names_to = "category", values_to = "count") %>%
  mutate(category = factor(category, levels = c("Single Nucleui", "Shared", "Single Cell"))) %>%
  ggplot(aes(x = cell_type, y = count, fill = category)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("sn_only" = "#8ca9ff", "shared" = "#7570b3", "sc_only" = "#778873")) +
  coord_flip() +
  labs(x = NULL, y = "States", fill = NULL) +
  theme_minimal(base_size = 5) +
  theme(legend.position = "bottom")


plot_df <- mps_tib %>%
  mutate(total = sn_only + shared + sc_only) %>%
  arrange(total) %>%
  mutate(cell_type = factor(cell_type, levels = cell_type)) %>%
  pivot_longer(c(sn_only, shared, sc_only),
               names_to = "category", values_to = "count") %>%
  mutate(category = factor(category, levels = c("sc_only", "shared", "sn_only")))

cols <- c("sc_only" = "#778873", "shared" = "#7570b3", "sn_only" = "#8ca9ff")
labs <- c("sc_only" = "Single Cell", "shared" = "Shared", "sn_only" = "Single Nucleui")

p <- ggplot(plot_df, aes(x = cell_type, y = count, fill = category)) +
  geom_col(width = 0.72) +
  #geom_text(aes(label = ifelse(count > 0, count, "")),
  #         position = position_stack(vjust = 0.5),
  #        size = 2.6, colour = "white", fontface = "bold") +
  scale_fill_manual(values = cols, labels = labs, name = NULL,
                    guide = guide_legend(reverse = TRUE)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.04))) +
  coord_flip() +
  labs(x = NULL, y = "Number of states") +
  theme_minimal(base_size = 7) +
  theme(
    axis.text.y      = element_text(colour = "black"),
    axis.text.x      = element_text(colour = "black"),
    axis.title.x     = element_text(margin = margin(t = 6)),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_line(colour = "grey88", linewidth = 0.3),
    legend.position  = "bottom",
    legend.key.size  = unit(0.35, "cm"),
    plot.margin      = margin(6, 10, 6, 6)
  )

p

ggsave("BMCA/plots/fig3/states_sc_sn_shared.pdf", p,
       width = 120, height = 90, units = "mm", dpi = 300)
# -----------------------------------------------------------------------------

themes <- list(
  "Cell Cycle"            = "CC",
  "Stress"                = "Stress|stress|HSP|hsp",
  "Interferon Response"   = "INFR",
  "Memory"                = "Memory",
  "Naive"                 = "Naive",
  "Cytotoxicity"          = "CT",
  "Regulation"            = "Regul|regul",
  "Antigen Presentation"  = "HEV|MHC|DC|[Aa]ntigen",
  "Hypoxia"               = "Hypoxia",
  "Respiration"           = "Respiration",
  "NK like" = "_NK"
)

cell_groups <- list(
  "Malignant" = c("Malignant"),
  "Immune" = c("CD8", "CD4", "NK", "GD", "B_Plasma", "Myeloid"),
  "Stromal" = c("Fibroblast", "Pericyte", "Endothelial"),
  "Brain" = c("Neuron", "Astrocyte", "Oligodendrocyte")
)

# ── State groups ──────────────────────────────────────────────────────────────
state_groups <- list(
  "General"                  = c("Cell Cycle", "Stress", "Interferon Response",
                                 "Antigen Presentation", "Hypoxia", "Respiration"),
  "Immune Cell Type Specific" = c("Memory", "Naive", "Cytotoxicity",
                                  "Regulation", "NK like")
)

# reverse lookup: cell_type -> its cell group
ct_to_group <- stack(cell_groups) %>%
  transmute(cell_type = values, group = as.character(ind))

# for each cell type, count how many programs match each theme's regex
heatmap_df <- imap_dfr(colnames_list, function(progs, cell_type) {
  # strip sc_/sn_ prefixes so themes match the bare program names
  bare <- gsub("^s[cn]_", "", progs)
  map_dfr(names(themes), function(th) {
    pat <- themes[[th]]
    tibble(
      cell_type = cell_type,
      theme     = th,
      count     = sum(str_detect(bare, pat))
    )
  })
}) %>%
  left_join(ct_to_group, by = "cell_type") %>%
  mutate(count = ifelse(count == 0, NA_integer_, count))   # 0 -> NA so tiles show white/blank

# ── Add state_group to heatmap_df ─────────────────────────────────────────────
heatmap_df <- heatmap_df %>%
  mutate(
    state_group = case_when(
      theme %in% state_groups$General                   ~ "General",
      theme %in% state_groups$`Immune Cell Type Specific` ~ "Immune Cell Type Specific",
      TRUE ~ "Other"
    ),
    state_group = factor(state_group, levels = c("General", "Immune Cell Type Specific")),
    theme = factor(theme, levels = c(
      # General first
      "Cell Cycle", "Stress", "Interferon Response",
      "Antigen Presentation", "Hypoxia", "Respiration",
      # Immune specific second
      "Memory", "Naive", "Cytotoxicity", "Regulation", "NK like"
    ))
  )

# ── Plot ──────────────────────────────────────────────────────────────────────
p <- ggplot(heatmap_df, aes(x = theme, y = cell_type, fill = count)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = ifelse(is.na(count), "", count)), size = 2, color = "white") +
  scale_fill_gradient(low = "grey", high = "black", na.value = "white") +
  facet_grid(
    group ~ state_group,        # ← rows = cell group, cols = state group
    scales = "free",
    space  = "free"
  ) +
  labs(x = NULL, y = NULL, fill = "Count") +
  theme_minimal(base_size = 7) +
  theme(
    axis.text.x  = element_text(angle = 40, hjust = 1),
    panel.grid   = element_blank(),
    strip.text.y = element_text(angle = 0),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  )

p

ggsave("BMCA/plots/fig3/shared_mps.pdf", p,
       width = 150, height = 120, units = "mm", dpi = 300)
# -----------------------------------------------------------------------------
patients <- c("lin_sc_BM08T", "diaz_SF14720")
cohorts <- c("lin_sc", "diaz")

sample <- patients[1]
cohort <- cohorts[1]

ct <- "Malignant"

expmat_all <- readRDS(paste0("BMCA/data/centered_cpm_", cohort, ".RDS"))

expmat_all <- as.matrix(expmat_all)

metadata <- as_tibble(readRDS("BMCA/data/metadata_all_studies.rds"))

metadata_fil_all <- metadata %>%
  filter(study == cohort) %>%
  filter(cell_type == ct) %>%
  dplyr::select(c(colnames(metadata)[1:6], seq_tech))

unique_patients <- unique(metadata_fil_all$patient)
expmat_all[1:5, 1:5]


sheets <- excel_sheets("BMCA/csv/all_MPS.xlsx")

tibble_list <- set_names(
  map(sheets, ~ read_excel("BMCA/csv/all_MPS.xlsx", sheet = .x)),
  sheets
)

sigs <- tibble_list$Malignant
sig_list <- sigs %>%
  pivot_longer(cols = everything(), names_to = "signature", values_to = "gene") %>%
  filter(!is.na(gene)) %>%
  group_by(signature) %>%
  summarise(genes = list(unique(gene)), .groups = "drop") %>%
  deframe()

names(sig_list)

metadata_fil <- metadata_fil_all %>%
  filter(patient == sample) %>%
  mutate(CellID = gsub("\\.", "-", CellID))

cells <- unique(metadata_fil$CellID)

expmat <- expmat_all[, colnames(expmat_all) %in% cells]

score_mat <- sigScores(expmat, sig_list, conserved.genes = 0.25)

score_df <- score_mat %>%
  as.data.frame() %>%
  rownames_to_column("CellID")

score_df <- score_df %>%
  rowwise() %>%
  mutate(
    max_value = max(c_across(-CellID), na.rm = TRUE),
    annotation = ifelse(
      max_value > 0.7,
      colnames(score_df)[which.max(c_across(-CellID)) + 1],
      "not_assigned"
    )
  ) %>%
  ungroup()

metadata_fil <- metadata_fil %>%
  left_join(score_df, by = "CellID") %>%
  mutate(
    annotation = factor(annotation, levels = c(colnames(sigs), "not_assigned"), ordered = TRUE)
  ) %>%
  arrange(annotation, desc(max_value))
expmat <- expmat - rowMeans(expmat)
mps <- c( "sc_Hypoxia", "sc_CC", "sc_PDAC_related",
          "sc_Neuronal_related_2")
sigs <- sigs %>% dplyr::select(all_of(mps))

# Heatmap preparation
genes <- unique(na.omit(unlist(sigs)))
genes <- genes[genes %in% rownames(expmat)]
cells <- metadata_fil$CellID
cells <- cells[cells %in% colnames(expmat)]

expmat_to_plot <- expmat[genes, cells]

# per-gene % of values in [-0.5, 0.5]
pct_in_range <- rowMeans(expmat_to_plot >= -0.5 & expmat_to_plot <= 0.5) * 100

hist(pct_in_range,
     breaks = 30,
     main   = "Per-gene % of values in [-0.5, 0.5]",
     xlab   = "% of cells with value in [-0.5, 0.5]",
     col    = "grey80",
     border = "white")

relevant_genes <- names(pct_in_range)[pct_in_range < 50]
# 3. remove them from matrix
expmat_to_plot <- expmat_to_plot[relevant_genes, ]

df_heatmap <- as.data.frame(expmat_to_plot)
df_heatmap$Gene <- rownames(expmat_to_plot)
col_co <- 4

df_long <- melt(df_heatmap, id.vars = "Gene", variable.name = "Cell", value.name = "Expression") %>%
  distinct(Gene, Cell, .keep_all = TRUE)

color_palette <- readRDS("code/inferGEP/data/hotmap.rds")
breaks <- seq(-col_co, col_co, length.out = length(color_palette) + 1)

df_long <- df_long %>%
  mutate(Expression_to_plot = pmin(pmax(Expression, -col_co), col_co),
         binned_Expression = cut(Expression_to_plot, breaks = breaks, include.lowest = TRUE, labels = FALSE))

sigs_long <- sigs %>%
  pivot_longer(cols = everything(), names_to = "Category", values_to = "Gene") %>%
  group_by(Gene) %>%
  filter(n() == 1) %>%
  ungroup()

df_long <- df_long %>%
  left_join(sigs_long, by = "Gene") %>%
  dplyr::rename(Genes_MP = Category) %>%
  left_join(metadata_fil %>% dplyr::select(CellID, annotation), by = c("Cell" = "CellID")) %>%
  dplyr::rename(Cells_MP = annotation)

df_long$binned_Expression <- factor(df_long$binned_Expression, levels = rev(sort(unique(df_long$binned_Expression))))

df_long <- df_long %>%
  filter(complete.cases(.)) %>%
  mutate(
    Gene = factor(Gene, levels = genes),
    Cell = factor(Cell, levels = cells),
    Genes_MP = factor(Genes_MP, levels = levels(df_long$Cells_MP))
  )

legend_names <- c(paste0(col_co), rep("", length(color_palette) - 2), paste0("-",col_co))
color_palette_reversed <- rev(color_palette)
mps <- c( "sc_Hypoxia", "sc_CC", "sc_PDAC_related",
          "sc_Neuronal_related_2")
df_long <- df_long %>%
  filter(Cells_MP %in% mps & Genes_MP %in% mps)

df_long <-  df_long %>% filter((Gene %in% relevant_genes))

mps <- c( "sc_Hypoxia", "sc_CC", "sc_PDAC_related",
          "sc_Neuronal_related_2")
intersect(df_long$Gene, sigs$sc_Neuronal_related_2)
intersect(df_long$Gene, sigs$sc_CC)
intersect(df_long$Gene, sigs$sc_Hypoxia)
intersect(df_long$Gene, sigs$sc_PDAC_related)

mp_labels <- c(
  "sc_CC" = "Cell Cycle",
  "sc_Hypoxia" = "Hypoxia",
  "sc_PDAC_related" = "Epi. diff. 1",
  "sc_Neuronal_related_2" = "MNR2"
)

exp_plot <- ggplot(df_long, aes(x = Cell, y = Gene, fill = binned_Expression, color = binned_Expression)) +
  geom_tile() +
  facet_grid(rows = vars(Genes_MP), scales = "free", space = "free",
             labeller = as_labeller(mp_labels),
             switch = "y") +
  scale_fill_manual(
    values = setNames(color_palette_reversed, levels(df_long$binned_Expression)),
    labels = legend_names,
    name = "Expression"
  ) +
  scale_color_manual(
    values = setNames(color_palette_reversed, levels(df_long$binned_Expression)),
    labels = legend_names,
    name = "Expression"
  ) +
  labs(title = paste0("BM08T (Colorectal Carc. BrM)"), x = "", y = "") +
  theme_minimal(base_size = 7) +
  theme(
    axis.text.x     = element_blank(),
    axis.text.y     = element_blank(),
    legend.position = "none",
    legend.title    = element_text(size = 6),
    legend.text     = element_text(size = 5),
    legend.key.size = unit(0.2, "cm"),
    legend.key      = element_rect(color = "black", linewidth = 0.2, fill = NA),
    strip.text.x    = element_text(size = 7),
    strip.text.y.left = element_text(size = 7, angle = 0),
    strip.placement = "outside",
    panel.spacing   = unit(0.1, "lines"),
    panel.border    = element_rect(color = "black", fill = NA, linewidth = 0.3)
  )
exp_plot

ggsave("BMCA/plots/fig3/exp_plot_lin_sc.pdf", exp_plot,
       width = 70, height = 40, units = "mm", dpi = 300)


legend <- get_plot_component(
  exp_plot + theme(legend.position = "right"),
  "guide-box-right"
)

ggsave("BMCA/plots/fig3/exp_plot_legend.pdf", plot = legend,
       width = 1, height = 2, units = "in", dpi = 300)


list_of_cna_cor_sig_tibs <- readRDS("BMCA/data/list_of_cna_cor_sig_tibs_MNRs.RDS")
names(list_of_cna_cor_sig_tibs)
to_keep <- c("biermann_sn_MBM05_sn", "biermann_sn_MBM11_sn", "lin_sc_BM03T",
             "lin_sc_BM04T", "lin_sc_BM10T_rep", "lin_sn_B2255687",
             "lin_sn_B2349898", "lin_sn_B2421638")
list_of_cna_cor_sig_tibs <- list_of_cna_cor_sig_tibs[to_keep]

list_of_cna_cor_sig_tib <- bind_rows(list_of_cna_cor_sig_tibs, .id = "patient")

meta_to_join <- metadata %>%
  filter(CellID %in% list_of_cna_cor_sig_tib$CellID) %>%
  dplyr::select(CellID, cell_type)

list_of_cna_cor_sig_tib <- list_of_cna_cor_sig_tib %>%
  left_join(meta_to_join, by = "CellID")

list_of_cna_cor_sig_tib$cell_type_state <- list_of_cna_cor_sig_tib$cell_type

mnrs_cells <- metadata %>%
  filter(state %in% c("Malignant_Neuronal_related_1", "Malignant_Neuronal_related_2", "Malignant_Neuronal_related_3")) %>%
  pull(CellID)

list_of_cna_cor_sig_tib$cell_type_state[list_of_cna_cor_sig_tib$CellID %in% mnrs_cells] <-"Malignant_Neuronal_related"

list_of_cna_cor_sig_tib <- list_of_cna_cor_sig_tib %>%
  filter(cell_type %in% c("Malignant",       "Mesenchymal",     "Myeloid" ))


ggplot(list_of_cna_cor_sig_tib, aes(x = CNAsig, color = cell_type_state)) +
  geom_density(size = 1) +
  facet_wrap(~ patient, scales = "free")

p <- ggplot(list_of_cna_cor_sig_tib,
            aes(x = CNAsig, color = cell_type_state)) +
  geom_density(size = 1) +
  facet_wrap(~ patient, scales = "free") +
  scale_color_manual(
    values = c(
      "Malignant"                  = "#740A03",
      "Malignant_Neuronal_related" = "#C3110C",
      "Mesenchymal"                = "#09637E",
      "Myeloid"                    = "#377EB8"
    ),
    name = "Cell type / state"
  )

ggsave(
  "BMCA/plots/figure_3/CNA_density_low_quality.png",
  plot = p,
  width = 10,
  height = 6,
  units = "in",
  dpi = 300          # ← low quality
)

# -----------------------------------------------------------------------------
patients <- c("lin_sc_BM08T", "diaz_SF14720")
cohorts <- c("lin_sc", "diaz")

sample <- patients[2]
cohort <- cohorts[2]

ct <- "Malignant"

expmat_all <- readRDS(paste0("BMCA/data/centered_cpm_", cohort, ".RDS"))

expmat_all <- as.matrix(expmat_all)

metadata <- as_tibble(readRDS("BMCA/data/metadata_all_studies.rds"))

metadata_fil_all <- metadata %>%
  filter(study == cohort) %>%
  filter(cell_type == ct) %>%
  dplyr::select(c(colnames(metadata)[1:6], seq_tech))

unique_patients <- unique(metadata_fil_all$patient)
expmat_all[1:5, 1:5]


sheets <- excel_sheets("BMCA/csv/all_MPS.xlsx")

tibble_list <- set_names(
  map(sheets, ~ read_excel("BMCA/csv/all_MPS.xlsx", sheet = .x)),
  sheets
)

sigs <- tibble_list$Malignant
sig_list <- sigs %>%
  pivot_longer(cols = everything(), names_to = "signature", values_to = "gene") %>%
  filter(!is.na(gene)) %>%
  group_by(signature) %>%
  summarise(genes = list(unique(gene)), .groups = "drop") %>%
  deframe()

names(sig_list)

metadata_fil <- metadata_fil_all %>%
  filter(patient == sample) %>%
  mutate(CellID = gsub("\\.", "-", CellID))

cells <- unique(metadata_fil$CellID)

expmat <- expmat_all[, colnames(expmat_all) %in% cells]

score_mat <- sigScores(expmat, sig_list, conserved.genes = 0.25)

score_df <- score_mat %>%
  as.data.frame() %>%
  rownames_to_column("CellID")

score_df <- score_df %>%
  rowwise() %>%
  mutate(
    max_value = max(c_across(-CellID), na.rm = TRUE),
    annotation = ifelse(
      max_value > 0.7,
      colnames(score_df)[which.max(c_across(-CellID)) + 1],
      "not_assigned"
    )
  ) %>%
  ungroup()

metadata_fil <- metadata_fil %>%
  left_join(score_df, by = "CellID") %>%
  mutate(
    annotation = factor(annotation, levels = c(colnames(sigs), "not_assigned"), ordered = TRUE)
  ) %>%
  arrange(annotation, desc(max_value))
expmat <- expmat - rowMeans(expmat)

mps <- c( "sn_Hypoxia", "sn_CC", "sn_Skin_Pigmentation",
          "sn_Neuronal_related_1")
sigs <- sigs %>% dplyr::select(all_of(mps))
# Heatmap preparation
genes <- unique(na.omit(unlist(sigs)))
genes <- genes[genes %in% rownames(expmat)]
cells <- metadata_fil$CellID
cells <- cells[cells %in% colnames(expmat)]

expmat_to_plot <- expmat[genes, cells]

# per-gene % of values in [-0.5, 0.5]
pct_in_range <- rowMeans(expmat_to_plot >= -0.5 & expmat_to_plot <= 0.5) * 100

hist(pct_in_range,
     breaks = 30,
     main   = "Per-gene % of values in [-0.5, 0.5]",
     xlab   = "% of cells with value in [-0.5, 0.5]",
     col    = "grey80",
     border = "white")

relevant_genes <- names(pct_in_range)[pct_in_range < 50]
# 3. remove them from matrix
expmat_to_plot <- expmat_to_plot[relevant_genes, ]

df_heatmap <- as.data.frame(expmat_to_plot)
df_heatmap$Gene <- rownames(expmat_to_plot)
col_co <- 4

df_long <- melt(df_heatmap, id.vars = "Gene", variable.name = "Cell", value.name = "Expression") %>%
  distinct(Gene, Cell, .keep_all = TRUE)

color_palette <- readRDS("code/inferGEP/data/hotmap.rds")
breaks <- seq(-col_co, col_co, length.out = length(color_palette) + 1)

df_long <- df_long %>%
  mutate(Expression_to_plot = pmin(pmax(Expression, -col_co), col_co),
         binned_Expression = cut(Expression_to_plot, breaks = breaks, include.lowest = TRUE, labels = FALSE))

sigs_long <- sigs %>%
  pivot_longer(cols = everything(), names_to = "Category", values_to = "Gene") %>%
  group_by(Gene) %>%
  filter(n() == 1) %>%
  ungroup()

df_long <- df_long %>%
  left_join(sigs_long, by = "Gene") %>%
  dplyr::rename(Genes_MP = Category) %>%
  left_join(metadata_fil %>% dplyr::select(CellID, annotation), by = c("Cell" = "CellID")) %>%
  dplyr::rename(Cells_MP = annotation)

df_long$binned_Expression <- factor(df_long$binned_Expression, levels = rev(sort(unique(df_long$binned_Expression))))

df_long <- df_long %>%
  filter(complete.cases(.)) %>%
  mutate(
    Gene = factor(Gene, levels = genes),
    Cell = factor(Cell, levels = cells),
    Genes_MP = factor(Genes_MP, levels = levels(df_long$Cells_MP))
  )

legend_names <- c(paste0(col_co), rep("", length(color_palette) - 2), paste0("-",col_co))
color_palette_reversed <- rev(color_palette)

df_long <- df_long %>%
  filter(Cells_MP %in% mps & Genes_MP %in% mps)

df_long <-  df_long %>% filter((Gene %in% relevant_genes))

mps <- c( "sn_Hypoxia", "sn_CC", "sn_Skin_Pigmentation",
          "sn_Neuronal_related_1")

intersect(df_long$Gene, sigs$sn_Neuronal_related_1)
intersect(df_long$Gene, sigs$sn_CC)
intersect(df_long$Gene, sigs$sn_Hypoxia)
intersect(df_long$Gene, sigs$sn_Skin_Pigmentation)

mp_labels <- c(
  "sn_CC" = "Cell Cycle",
  "sn_Hypoxia" = "Hypoxia",
  "sn_Skin_Pigmentation" = "Skin Pig.",
  "sn_Neuronal_related_1" = "MNR1"
)


exp_plot <- ggplot(df_long, aes(x = Cell, y = Gene, fill = binned_Expression, color = binned_Expression)) +
  geom_tile() +
  facet_grid(rows = vars(Genes_MP), scales = "free", space = "free",
             labeller = as_labeller(mp_labels),
             switch = "y") +
  scale_fill_manual(
    values = setNames(color_palette_reversed, levels(df_long$binned_Expression)),
    labels = legend_names,
    name = "Expression"
  ) +
  scale_color_manual(
    values = setNames(color_palette_reversed, levels(df_long$binned_Expression)),
    labels = legend_names,
    name = "Expression"
  ) +
  labs(title = paste0("SF14720 (Melanoma BrM)"), x = "", y = "") +
  theme_minimal(base_size = 7) +
  theme(
    axis.text.x     = element_blank(),
    axis.text.y     = element_blank(),
    legend.position = "none",
    legend.title    = element_text(size = 6),
    legend.text     = element_text(size = 5),
    legend.key.size = unit(0.2, "cm"),
    legend.key      = element_rect(color = "black", linewidth = 0.2, fill = NA),
    strip.text.x    = element_text(size = 7),
    strip.text.y.left = element_text(size = 7, angle = 0),
    strip.placement = "outside",
    panel.spacing   = unit(0.1, "lines"),
    panel.border    = element_rect(color = "black", fill = NA, linewidth = 0.3)
  )
exp_plot

ggsave("BMCA/plots/fig3/exp_plot_diaz.pdf", exp_plot,
       width = 70, height = 40, units = "mm", dpi = 300)


# -----------------------------------------------------------------------------
sigs <- as_tibble(read.csv(paste0("BMCA/csv/mps_MNRs.csv")))[1:25,1:9]

# turn each column into a clean character vector (drop NA / blanks)
gene_sets <- lapply(sigs, function(x) {
  x <- x[!is.na(x) & x != ""]
  unique(x)
})

nf_sets <- gene_sets[grep("^NF_", names(gene_sets))]
p_sets  <- gene_sets[grep("^P_",  names(gene_sets))]

jaccard <- function(a, b) {
  length(intersect(a, b)) / length(union(a, b))
}

# rows = NF_, cols = P_
J <- outer(
  seq_along(nf_sets), seq_along(p_sets),
  Vectorize(function(i, j) jaccard(nf_sets[[i]], p_sets[[j]]))
)

pheatmap(
  J,
  display_numbers = TRUE,        # show Jaccard values in cells
  number_format = "%.2f",
  cluster_rows = FALSE,          # keep your row/col order; set TRUE to cluster
  cluster_cols = FALSE,
  color = colorRampPalette(c("white", "#2166AC"))(100),
  main = "Jaccard overlap: NF_ vs P_ signatures"
)


J_long <- as.data.frame(as.table(J))
names(J_long) <- c("NF", "P", "Jaccard")

p <- ggplot(J_long, aes(P, NF, fill = Jaccard)) +
  geom_tile(color = "grey90") +
  geom_text(aes(label = sprintf("%.2f", Jaccard)), size = 2) +
  scale_fill_gradient(low = "white", high = "#2166AC", limits = c(0, NA)) +
  theme_minimal(base_size = 6) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_blank()) +
  labs(title = "")
p
ggsave("BMCA/plots/fig3/MNRs_jaccard_heatmap.pdf", plot = p, width = 80, height = 80,
       units = "mm", dpi = 300)


# 1. Setup
setEnrichrSite("Enrichr")
target_databases <- c("GO_Biological_Process_2025", "KEGG_2026", "MSigDB_Hallmark_2020",
                      "GO_Cellular_Component_2025", "GO_Molecular_Function_2025",
                      "CellMarker_2024", "ProteomicsDB_2020", "Human_Gene_Atlas")
# 2. Run the loop to collect hits
raw_hits <- data.frame()

for (col in colnames(sigs)) {
  message("Analyzing: ", col)

  # Remove NAs and empty strings
  sig_genes <- sigs[[col]][!is.na(sigs[[col]]) & sigs[[col]] != ""]

  # Run enrichr
  enrich_res <- enrichr(sig_genes, target_databases)

  # Process and filter using Case-Insensitive Regex
  res_processed <- enrich_res %>%
    map_df(as.data.frame, .id = "source_db") %>%
    # Use | for OR. This search is not case sensitive because of ignore.case = TRUE
    filter(grepl("brain|cortex|glia|neuron|synap|astrocyt|oligodendro|microglia",
                 Term,
                 ignore.case = TRUE)) %>%
    mutate(Signature = col)

  raw_hits <- bind_rows(raw_hits, res_processed)
}

# 3. Create the "Universe" of all unique Neuro terms found across all signatures
brain_universe <- raw_hits %>%
  distinct(Term, source_db)

# 4. Expand grid to include every signature vs every found neuro term (even 0 overlaps)
template <- expand_grid(
  Signature = colnames(sigs),
  Term = brain_universe$Term
) %>%
  left_join(brain_universe, by = "Term")

# 5. Final Join and formatting the Tibble
final_complete_tibble <- template %>%
  left_join(
    raw_hits %>% dplyr::select(Signature, Term, pvalue = P.value, adj_pvalue = Adjusted.P.value, Genes),
    by = c("Signature", "Term")
  ) %>%
  mutate(
    # Set significance to 0 for missing terms
    `-log10(fdr)` = ifelse(is.na(adj_pvalue), 0, -log10(adj_pvalue)),
    pvalue = ifelse(is.na(pvalue), 1, pvalue),
    adj_pvalue = ifelse(is.na(adj_pvalue), 1, adj_pvalue),
    Genes = ifelse(is.na(Genes), "", Genes)
  ) %>%
  as_tibble() %>%
  arrange(Signature, desc(`-log10(fdr)`))

# View results
print(final_complete_tibble)


# 1. Identify the top terms to keep the plot readable
top_terms <- final_complete_tibble %>%
  dplyr::group_by(Term) %>%
  dplyr::summarize(max_sig = max(`-log10(fdr)`, na.rm = TRUE)) %>%
  dplyr::arrange(desc(max_sig)) %>%
  dplyr::slice_head(n = 12) %>%
  dplyr::pull(Term)
# 2. Filter the tibble for these top terms
plot_data <- final_complete_tibble %>%
  filter(Term %in% top_terms)

# 3. Create the faceted histogram
ggplot(plot_data, aes(x = `-log10(fdr)`)) +
  # Add a vertical line for the significance threshold (p=0.05 -> 1.3)
  geom_vline(xintercept = 1.301, linetype = "dashed", color = "red", alpha = 0.7) +
  geom_histogram(bins = 15, fill = "steelblue", color = "white") +
  # Facet by Term, allow the names to wrap
  facet_wrap(~Term, ncol = 3, labeller = label_wrap_gen(width = 30)) +
  theme_minimal() +
  labs(
    title = "FDR Distribution for Top Neuro-Related Terms",
    subtitle = "Red dashed line indicates significance (FDR < 0.05)",
    x = "-log10(FDR)",
    y = "Count of Signatures"
  ) +
  theme(
    strip.text = element_text(size = 8, face = "bold"),
    panel.spacing = unit(1, "lines")
  )
summary_plot_data <- final_complete_tibble %>%
  dplyr::group_by(Term, source_db) %>%
  dplyr::summarize(
    mean_score = mean(`-log10(fdr)`, na.rm = TRUE),
    sd_score   = sd(`-log10(fdr)`, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::filter(mean_score > 0) %>%
  dplyr::arrange(desc(mean_score)) %>%
  dplyr::slice_head(n = 30)

# 2. Create the Bar Plot
ggplot(summary_plot_data, aes(x = reorder(Term, mean_score), y = mean_score, fill = source_db)) +
  geom_col(color = "white", width = 0.8) +
  # Add a line for the standard significance threshold (FDR = 0.05)
  geom_hline(yintercept = 1.301, linetype = "dashed", color = "darkred", alpha = 0.5) +
  coord_flip() + # Flip for better readability of long GO terms
  theme_minimal() +
  scale_fill_brewer(palette = "Set2") +
  labs(
    title = "Average Neuro-Enrichment across All Signatures",
    subtitle = "Ranked by mean -log10(FDR). Dashed line at 1.3 represents p = 0.05.",
    x = "Enrichment Term",
    y = "Average -log10(FDR)",
    fill = "Database Source"
  ) +
  theme(
    axis.text.y = element_text(size = 9),
    legend.position = "bottom"
  )

final_complete_tibble <- final_complete_tibble %>%
  filter(Term %in% c("Generation of Neurons (GO:0048699)", "Neuron Brain Mouse",
                     "PrefrontalCortex"))

final_complete_tibble <- final_complete_tibble %>%
  group_by(Signature) %>%
  # Keep the row with the maximum -log10(fdr) value
  slice_min(`adj_pvalue`, n = 1, with_ties = FALSE) %>%
  ungroup()

final_complete_tibble_2 <- final_complete_tibble %>%
  # Create a 'Prefix' column by extracting everything up to the second underscore
  # And update 'Signature' by removing that same pattern
  mutate(
    Prefix = str_extract(Signature, "^[^_]+"),
    Signature = str_remove(Signature, "^[^_]+")
  ) %>%
  # Reorder columns to put Prefix first
  dplyr::select(Prefix, Signature, everything())


# Remove any leading underscores from the Signature column
final_complete_tibble_2 <- final_complete_tibble_2 %>%
  mutate(Signature = str_remove(Signature, "^_+"))

# View the cleaned result
print(final_complete_tibble_2)
# View the result
print(final_complete_tibble)


# 1. Clean data: Convert NAs to 1 and apply the p < 0.05 threshold
plot_data <- final_complete_tibble_2 %>%
  mutate(
    # Convert NA p-values to 1 (non-significant)
    pvalue = replace_na(pvalue, 1),

    # Calculate -log10 of the regular p-value
    neg_log10_p = -log10(pvalue),

    # If p >= 0.05 (neg_log10_p < 1.301), make it NA to trigger grey color
    plot_val = ifelse(neg_log10_p < 0.11, NA, neg_log10_p),

    # Shorten term names
    Term_Short = str_trunc(Term, 40)
  )

plot_data$Term_Short <- "Highest Enriched Neuronal Signature"
plot_data <- plot_data %>%
  mutate(Prefix = case_when(
    Prefix == "NF" ~ "This Study",
    Prefix == "P" ~ "Previously Published",
    TRUE ~ Prefix # Keeps others as they are
  ))
plot_data <- plot_data %>%
  mutate(Prefix = factor(Prefix, levels = c("This Study", "Previously Published")))
# 2. Create the Heatmap
ggplot(plot_data, aes(x = Term_Short, y = Signature, fill = plot_val)) +
  geom_tile(color = "white", size = 0.3) +

  # Gradient starts at p = 0.05 (1.301)
  scale_fill_gradient(
    low = "#e0f3f8", # Light blue for p ~ 0.05
    high = "#084594", # Dark blue for highly significant p-values
    na.value = "grey85",
    name = "-log10(p-value)"
  ) +

  facet_grid(Prefix ~ ., scales = "free_y", space = "free_y") +

  theme_minimal() +
  labs(
    title = " ",
    subtitle = " ",
    x = "Enrichment Term",
    y = ""
  ) +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    strip.background = element_rect(fill = "grey95", color = NA),
    strip.text = element_text(face = "bold"),
    panel.grid = element_blank()
  )

final_complete_tibble_2 <- final_complete_tibble_2 %>%
  dplyr::select(Prefix, Signature, Term, source_db, pvalue)

sig_enrich <- plot_data %>%
  dplyr::select(Signature, Prefix, neg_log10_p)

sig_enrich <- sig_enrich %>%
  mutate(significance = case_when(
    neg_log10_p > 3            ~ "***",
    neg_log10_p > 2            ~ "**",
    neg_log10_p > -log10(0.05) ~ "*",     # 1.301
    TRUE                       ~ "ns"
  ))
# -----------------------------------------------------------------------------

ct <- "Malignant"
study_pb_list_all <- readRDS(paste0("BMCA/data/UMI_pb_", ct, ".RDS"))


# 1. Get the union of all unique gene names across all matrices
all_genes <- unique(unlist(lapply(study_pb_list_all, rownames)))

# 2. Align each matrix to this master list of genes
aligned_list <- lapply(study_pb_list_all, function(mat) {
  # Identify missing genes
  missing_genes <- setdiff(all_genes, rownames(mat))

  # Create a sparse matrix of zeros for the missing genes
  # Note: Sparse matrices treat "NA" as 0 to maintain efficiency.
  if (length(missing_genes) > 0) {
    zero_mat <- Matrix(0,
                       nrow = length(missing_genes),
                       ncol = ncol(mat),
                       sparse = TRUE,
                       dimnames = list(missing_genes, colnames(mat)))

    # Combine and reorder to match the master gene list
    res <- rbind(mat, zero_mat)
    return(res[all_genes, , drop = FALSE])
  } else {
    return(mat[all_genes, , drop = FALSE])
  }
})

# 3. Use do.call with cbind to merge them all
combined_matrix <- do.call(cbind, aligned_list)

# Calculate the mean for each row (gene)
gene_means <- rowMeans(combined_matrix)

# Quick look at the distribution of averages

# Set the number of top genes you want to keep
n_top <- 3000

# Get the names of the genes with the highest averages
top_genes <- names(sort(gene_means, decreasing = TRUE))[1:n_top]

sigs_flat <- as.character(unlist(sigs))
genes <- unique(c(sigs_flat,top_genes))
# Subset the matrix
combined_matrix <- combined_matrix[rownames(combined_matrix) %in% genes, ]

combined_matrix[1:5, 1:5]
combined_matrix <- log2(combined_matrix + 1)

sigs_list <- as.list(sigs)

sig_scores_manual <- function(mat, sig_list, conserved.genes = 0.25) {
  # mat: genes x cells, centered log-expression; rownames = gene symbols
  detected <- rowMeans(mat != 0) >= conserved.genes
  vapply(sig_list, function(genes) {
    g <- intersect(genes, rownames(mat)[detected])
    if (length(g) == 0) return(rep(NA_real_, ncol(mat)))
    colMeans(mat[g, , drop = FALSE], na.rm = TRUE)
  }, numeric(ncol(mat)))
}
combined_matrix <- as.matrix(combined_matrix)
score_mat <- sig_scores_manual(combined_matrix, sigs_list, conserved.genes = 0.5)

metadata_p <- metadata %>%
  dplyr::select(patient, location, primary_site) %>%
  distinct()

score_mat <- as_tibble(rownames_to_column(as.data.frame(score_mat), "patient"))
# 1. Prepare the data
# Ensure patient IDs are in a column for joining
merged_df <- inner_join(score_mat, metadata_p, by = "patient")

# Define our special signatures
melanoma_sigs <- c("NF_sn_Neuronal_related_1", "P_Biermann_et_al_MP7")
all_sigs <- setdiff(colnames(score_mat), "patient")

# 2. Perform t-tests with conditional filtering
results <- lapply(all_sigs, function(sig) {

  # Apply filtering logic based on the signature
  if (sig %in% melanoma_sigs) {
    # Rule 1: Use ONLY melanoma samples
    sub_data <- merged_df %>% filter(primary_site == "melanoma")
  } else {
    # Rule 2: EXCLUDE melanoma samples
    sub_data <- merged_df %>% filter(primary_site != "melanoma")
  }

  # Split into groups for the t-test
  intra <- sub_data %>% filter(location == "Intracranial") %>% pull(!!sym(sig))
  extra <- sub_data %>% filter(location == "Extracranial") %>% pull(!!sym(sig))

  # Safety check: Ensure we have enough samples in both groups to test
  if(length(intra) < 3 | length(extra) < 3) {
    return(data.frame(Signature = sig, Delta = NA, PValue = NA, Significance = "N/A"))
  }

  # Conduct t-test
  t_res <- t.test(intra, extra)

  # Calculate delta: mean(Intra) - mean(Extra)
  delta <- mean(intra, na.rm = TRUE) - mean(extra, na.rm = TRUE)
  p_val <- t_res$p.value

  data.frame(
    Signature = sig,
    Delta = delta,
    PValue = p_val,
    Significance = ifelse(p_val < 0.05, "*", "")
  )
}) %>% bind_rows()
results$Signature <- gsub("^(NF_|P_)", "", results$Signature)

levs <- c("sn_Neuronal_related_3","sn_Neuronal_related_2", "sn_Neuronal_related_1",
          "sc_Neuronal_related_2", "sc_lncRNA_MT", "Xing_et_al_MP9",
          "Tagore_et_al_MP7", "Tagore_et_al_MP6", "Biermann_et_al_MP7")
levs <- rev(levs)
results <- results %>%
  mutate(Signature = factor(Signature, levels = levs))
# 3. Generate the refined Heatmap
ggplot(results, aes(x = "Intra vs Extra", y = Signature, fill = Delta)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "darkblue", mid = "white", high = "darkred",
                       midpoint = 0.5, na.value = "grey90") +
  labs(
    title = "Adjusted Signature Shift (Location-based)",
    subtitle = "Melanoma sigs: Melanoma only | Other sigs: Non-melanoma only\n* = p < 0.05",
    x = "",
    fill = "Delta Mean"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(),
    panel.grid = element_blank()
  )


# 1. Prepare the data: Join score_mat with metadata_p
# Assuming rownames(score_mat) match metadata_p$patient
merged_df <- inner_join(score_mat, metadata_p, by = "patient")

# Define the signatures to test (all columns except patient and location)
signatures <- setdiff(colnames(score_mat), "patient")

# 2. Perform t-tests and calculate deltas
results <- lapply(signatures, function(sig) {
  # Split data by location
  intra <- merged_df %>% filter(location == "Intracranial") %>% pull(!!sym(sig))
  extra <- merged_df %>% filter(location == "Extracranial") %>% pull(!!sym(sig))

  # Conduct t-test
  t_res <- t.test(intra, extra)

  # Calculate delta: mean(Intra) - mean(Extra)
  delta <- mean(intra, na.rm = TRUE) - mean(extra, na.rm = TRUE)
  p_val <- t_res$p.value

  data.frame(
    Signature = sig,
    Delta = delta,
    PValue = p_val,
    Significance = ifelse(p_val < 0.05, "*", "")
  )
}) %>% bind_rows()

# 3. Create the One-Column Heatmap
ggplot(results, aes(x = "Intra vs Extra", y = Signature, fill = Delta)) +
  geom_tile(color = "white") +
  # Add the significance star inside the tile
  geom_text(aes(label = Significance), color = "black", size = 8, vjust = 0.7) +
  # Define the color scale (Red for higher in Intra, Blue for higher in Extra)
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  labs(
    title = "Signature Shift: Intracranial vs Extracranial",
    subtitle = "Values = Delta Mean; * = p < 0.05 (t-test)",
    x = "",
    fill = "Delta Mean"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(face = "bold"),
    panel.grid = element_blank()
  )
results
# -----------------------------------------------------------------------------
# ============================================================
# NF_ (rows) vs P_ (cols) context-specific correlation heatmap
# Each P_ column computed on its OWN cell subset. Rows & cols
# ordered by hierarchical clustering (average linkage).
# ============================================================
score_tib <- readRDS("BMCA/data/score_tib_mnrs.RDS")
# ---- 0. Attach the metadata fields we filter on ---------------------------
score_meta <- score_tib %>%
  left_join(
    metadata %>% dplyr::select(CellID, patient, location, primary_site, seq_tech),
    by = "CellID"
  )

# ---- 0b. VERIFY filter values before trusting anything --------------------
message("locations:    ", paste(unique(metadata$location),     collapse = ", "))
message("primary_site: ", paste(unique(metadata$primary_site), collapse = ", "))
message("seq_tech:     ", paste(unique(metadata$seq_tech),     collapse = ", "))
message("join matched ", sum(!is.na(score_meta$location)), " of ",
        nrow(score_meta), " cells")

# ---- 1. Signature columns -------------------------------------------------
nf_cols <- grep("^NF_", names(score_tib), value = TRUE)

# ---- 2. Per-column cell filters (each P_ gets its OWN subset) -------------
p_filters <- list(
  P_Biermann_et_al_MP7 = quote(location == "Intracranial" & primary_site == "melanoma"),
  P_Tagore_et_al_MP6   = quote(location == "Intracranial" & primary_site == "NSCLC"),
  P_Tagore_et_al_MP7   = quote(location == "Intracranial" & primary_site == "NSCLC"),
  P_Xing_et_al_MP9     = quote(location == "Intracranial" & seq_tech == "sc")
)

min_cells <- 20   # columns with fewer usable cells become NA

cor_mat <- matrix(NA_real_, nrow = length(nf_cols), ncol = length(p_filters),
                  dimnames = list(nf_cols, names(p_filters)))
p_mat   <- cor_mat                       # same shape, holds raw p-values
n_cells <- setNames(integer(length(p_filters)), names(p_filters))
for (p in names(p_filters)) {
  keep <- with(score_meta, eval(p_filters[[p]]))
  keep[is.na(keep)] <- FALSE
  sub <- score_meta[keep, ]
  n_cells[p] <- nrow(sub)

  if (nrow(sub) < min_cells) {
    warning(sprintf("Column '%s' has only %d cells (< %d) -> set to NA",
                    p, nrow(sub), min_cells))
    next
  }

  for (nf in nf_cols) {
    ct <- suppressWarnings(
      cor.test(sub[[nf]], sub[[p]], method = "spearman", exact = FALSE)
    )
    cor_mat[nf, p] <- ct$estimate
    p_mat[nf, p]   <- ct$p.value
  }
}
print(n_cells)
print(round(cor_mat, 3))

# ---- 4. Hierarchical clustering (rows and cols SEPARATELY) ----------------
# NA -> 0 for clustering only (does not affect displayed values).
cm <- cor_mat
cm[is.na(cm)] <- 0

# rows: cluster NF signatures on their profile across P columns
row_order <- if (nrow(cm) > 2) {
  hc_r <- hclust(dist(cm), method = "average")
  hc_r$labels[hc_r$order]
} else rownames(cm)

# cols: cluster P signatures on their profile across NF rows
col_order <- if (ncol(cm) > 2) {
  hc_c <- hclust(dist(t(cm)), method = "average")
  hc_c$labels[hc_c$order]
} else colnames(cm)

p_long <- p_mat %>%
  as.data.frame() %>%
  rownames_to_column("NF") %>%
  pivot_longer(-NF, names_to = "P", values_to = "pval")

cor_long <- cor_mat %>%
  as.data.frame() %>%
  rownames_to_column("NF") %>%
  pivot_longer(-NF, names_to = "P", values_to = "r") %>%
  left_join(p_long, by = c("NF", "P")) %>%
  mutate(
    r_capped = pmax(pmin(r, 0.5), -0.5),
    padj     = p.adjust(pval, method = "BH"),   # <- FDR across all tiles
    star     = ifelse(!is.na(padj) & padj < 0.05, "*", ""),
    NF = factor(NF, levels = row_order),
    P  = factor(P,  levels = col_order)
  )
# x-axis labels annotated with n per column (respect clustered order)
p_labels <- setNames(
  sprintf("%s\n(n = %d)", names(n_cells), n_cells),
  names(n_cells)
)[col_order]

# ---- Print the current row/col names (in the order they appear) -----------
cat("ROWS (y-axis, top -> bottom as plotted):\n")
print(rev(row_order))     # rev() because y-axis plots bottom-up
cat("\nCOLS (x-axis, left -> right):\n")
print(col_order)

# ---- Manual display labels (edit the right-hand sides) --------------------
row_labels <- c(
  NF_sc_Neuronal_related_2 = "MNR 2 (sc)",
  NF_sc_lncRNA_MT          = "lncRNA-MT genes",
  NF_sn_Neuronal_related_1 = "MNR 1 (sn)",
  NF_sn_Neuronal_related_2 = "MNR 2 (sn)",
  NF_sn_Neuronal_related_3 = "MNR 3 (sn)"
)

col_labels <- c(
  P_Biermann_et_al_MP7 = "Biermann MP7",
  P_Tagore_et_al_MP6   = "Tagore MP6",
  P_Tagore_et_al_MP7   = "Tagore MP7",
  P_Xing_et_al_MP9     = "Xing MP9"
)
# ---- 6. Heatmap -----------------------------------------------------------
p_heat <- ggplot(cor_long, aes(P, NF, fill = r_capped)) +
  geom_tile(color = "white", linewidth = 0.4) +
  #geom_text(aes(label = star), size = 5, vjust = 0.75) +      # * instead of %.2f
  scale_fill_gradient2(low = "darkblue", mid = "white", high = "darkred",
                       midpoint = 0, limits = c(-0.5, 0.5),
                       name = "Spearman r") +
  scale_x_discrete(limits = col_order, labels = col_labels, position = "top") +
  scale_y_discrete(limits = row_order, labels = row_labels) +
  coord_fixed() +
  theme_minimal(base_size = 20) +
  theme(axis.text.x.top = element_text(angle = 45, hjust = 0, vjust = 0,
                                       margin = margin(b = -2)),
        axis.title = element_blank(),
        panel.grid = element_blank())

n_col <- length(col_order)
n_row <- length(row_order)

# ---- shared legend guide: horizontal bar, outline, title below ------------
leg_guide <- guides(fill = guide_colorbar(
  direction       = "horizontal",   # horizontal bar
  title.position  = "bottom",       # title under the bar
  title.hjust     = 0.5,            # centered
  frame.colour    = "black",        # outline around the bar
  frame.linewidth = 0.4,
  ticks.colour    = "black",
  ticks.linewidth = 0.4,
  barwidth        = unit(2, "cm"),
  barheight       = unit(0.3, "cm")
))

# ---- main heatmap ---------------------------------------------------------
p_heat <- ggplot(cor_long, aes(P, NF, fill = r_capped)) +
  geom_tile(color = "white", linewidth = 0.4) +
  #geom_text(aes(label = star), size = 5, vjust = 0.75) +
  annotate("rect", xmin = 0.5, xmax = n_col + 0.5, ymin = 0.5, ymax = n_row + 0.5,
           fill = NA, color = "black", linewidth = 0.4) +
  scale_fill_gradient2(low = "darkblue", mid = "white", high = "darkred",
                       midpoint = 0, limits = c(-0.5, 0.5),
                       name = "Spearman r") +
  scale_x_discrete(limits = col_order, labels = col_labels, position = "top") +
  scale_y_discrete(limits = row_order, labels = row_labels) +
  coord_fixed() +
  theme_minimal(base_size = 7) +
  theme(axis.text.x.top  = element_text(angle = 45, hjust = 0, vjust = 0,
                                        margin = margin(b = -2)),
        axis.text.y      = element_blank(),
        axis.ticks.y     = element_blank(),
        axis.title       = element_blank(),
        panel.grid       = element_blank(),
        legend.position  = "bottom",
        legend.direction = "horizontal") +
  leg_guide
p_heat

# ---- signature prefixing / joins (unchanged) ------------------------------
sig_enrich <- sig_enrich %>%
  mutate(Signature = paste0(
    c(rep("NF_", 5), rep("P_", 4)),
    Signature
  ))

sig_enrich <- sig_enrich %>%
  mutate(Signature = sub("^(NF_|P_)\\1", "\\1", Signature))

anno_data <- sig_enrich %>% left_join(results, by = "Signature")
colnames(anno_data) <- c("Signature", "Study",
                         "-log10(p-value) Highest Enriched Neuronal Signature",
                         "Significance Neuronal Enrichment",
                         "BrM - Primary Site Enrichment",
                         "PValue", "Significance (BrM - Primary Site Enrichment)")


# ---- 7. NF-only annotation data, aligned to heatmap rows ------------------
anno_nf <- anno_data %>%
  filter(Signature %in% row_order) %>%
  mutate(
    Signature = factor(Signature, levels = row_order),
    Neuronal  = `-log10(p-value) Highest Enriched Neuronal Signature`,
    BrM       = `BrM - Primary Site Enrichment`,
    BrM_star  = ifelse(grepl("\\*", `Significance (BrM - Primary Site Enrichment)`),
                       "*", "")
  )

# shared theme for the annotation strips — x labels on top to match the heatmap
anno_theme <- theme_minimal(base_size = 7) +
  theme(axis.title       = element_blank(),
        axis.text.x.top  = element_text(angle = 0,
                                        margin = margin(b = 0)),
        axis.text.x      = element_blank(),
        axis.text.y      = element_blank(),
        axis.ticks.y     = element_blank(),
        panel.grid       = element_blank(),
        legend.position  = "bottom",
        legend.direction = "horizontal")   # <- horizontal bar

# ---- Neuronal enrichment strip -------------------------------------------
n_sig <- nlevels(anno_nf$Signature)

p_neur <- ggplot(anno_nf, aes(x = 1, y = Signature, fill = Neuronal)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = ifelse(`Significance Neuronal Enrichment` == "ns", "", "*")),
            color = "black", size = 10, vjust = 0.75) +
  annotate("rect", xmin = 0.5, xmax = 1.5, ymin = 0.5, ymax = n_sig + 0.5,
           fill = NA, color = "black", linewidth = 0.4) +
  scale_fill_gradient(low = "white", high = "#08519c",
                      name = "Highest Neuronal Sig.\nEnrichment (-log10(p))") +
  scale_x_continuous(breaks = 1, labels = "Highest Neuronal Sig.\nEnrichment (-log10(p))",
                     position = "top") +
  coord_fixed() + anno_theme + leg_guide
p_neur

# ---- BrM enrichment strip -------------------------------------------------
p_brm <- ggplot(anno_nf, aes(x = 1, y = Signature, fill = BrM)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = BrM_star), size = 10, vjust = 0.75) +
  annotate("rect", xmin = 0.5, xmax = 1.5, ymin = 0.5, ymax = n_sig + 0.5,
           fill = NA, color = "black", linewidth = 0.4) +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b",
                       midpoint = 0, na.value = "grey90",
                       name = "BrM - Primary\nEnrichment") +
  scale_x_continuous(breaks = 1, labels = "BrM - Primary\nEnrichment",
                     position = "top") +
  coord_fixed() + anno_theme + leg_guide
p <- p_neur + plot_spacer() + p_brm + plot_spacer() + p_heat +
  plot_layout(widths = c(1, 0.1, 1, 0.01, 3.83)) &
  theme(legend.title = element_text(size = 5),   # legend title text
        legend.text  = element_text(size = 5))   # numeric tick labels on the bar
p
p

ggsave("BMCA/plots/fig3/mnr_heats.pdf", plot = p, width = 150, height = 100,
       unit = "mm", dpi = 300)
ggsave("BMCA/plots/fig3/p_heat.pdf", plot = p_heat, width = 45, height = 65,
       unit = "mm", dpi = 300)
ggsave("BMCA/plots/fig3/p_neur.pdf", plot = p_neur, width = 25,  height = 65,
       unit = "mm", dpi = 300)
ggsave("BMCA/plots/fig3/p_brm.pdf",  plot = p_brm, width = 20,  height = 65,
       unit = "mm", dpi = 300)

leg_guide <- guides(fill = guide_colorbar(title.position = "bottom",
                                          title.hjust    = 0.5,
                                          frame.colour    = "black",
                                          frame.linewidth = 0.4,
                                          ticks.colour    = "black"))
leg_theme <- theme(legend.position  = "bottom",
                   legend.direction = "horizontal",
                   legend.title = element_text(size = 6),
                   legend.text  = element_text(size = 5))

leg_heat <- cowplot::get_legend(p_heat + leg_theme + leg_guide)
leg_neur <- cowplot::get_legend(p_neur + leg_theme + leg_guide)
leg_brm  <- cowplot::get_legend(p_brm  + leg_theme + leg_guide)

ggsave("BMCA/plots/fig3/leg_heat.pdf", plot = leg_heat, width = 3, height = 1.5)
ggsave("BMCA/plots/fig3/leg_neur.pdf", plot = leg_neur, width = 3, height = 1.5)
ggsave("BMCA/plots/fig3/leg_brm.pdf",  plot = leg_brm,  width = 3, height = 1.5)


# -----------------------------------------------------------------------------


sheets <- excel_sheets("BMCA/csv/all_MPS.xlsx")

tibble_list <- set_names(
  map(sheets, ~ read_excel("BMCA/csv/all_MPS.xlsx", sheet = .x)),
  sheets
)

tibble_list$Astrocyte <- tibble_list$Astrocyte %>%
  rename_with(~ paste0("sn_", .x))

tibble_list$Neuron <- tibble_list$Neuron %>%
  rename_with(~ paste0("sn_", .x))

colname_list <- lapply(tibble_list, colnames)

imap_dfr(
  colname_list,
  ~ tibble(
    cell_type = .y,
    program   = .x
  )
) %>%
  mutate(
    prefix = if_else(str_starts(program, "sn_"), "sn", "sc"),
    base_program = str_remove(program, "^(sn_|sc_)")
  ) %>%
  distinct(cell_type, base_program, prefix) %>%
  group_by(cell_type, base_program) %>%
  summarise(
    category = case_when(
      all(c("sn", "sc") %in% prefix) ~ "sn_and_sc",
      prefix == "sn"                 ~ "sn_only",
      prefix == "sc"                 ~ "sc_only"
    ),
    .groups = "drop"
  ) %>%
  count(cell_type, category) %>%
  group_by(cell_type) %>%
  mutate(total = sum(n)) %>%
  ungroup() %>%
  mutate(
    cell_type = fct_reorder(cell_type, total, .desc = TRUE),
    category  = factor(category, levels = c("sn_only", "sn_and_sc", "sc_only"))
  ) %>%
  ggplot(aes(x = cell_type, y = n, fill = category)) +
  geom_col() +
  geom_text(
    aes(y = total, label = total),
    hjust = -0.1         # a bit outside the bar
  ) +
  coord_flip() +
  scale_fill_manual(
    values = c(
      sn_only   = "#8ca9ff",
      sn_and_sc = "#7570b3",
      sc_only   = "#778873"
    ),
    breaks = c("sn_only", "sn_and_sc", "sc_only"),
    labels = c("Single Nuclues only", "Both", "Single Cels only"),
    name   = "Modality"
  ) +
  labs(
    x = "Cell type",
    y = "Number of Meta Programs",
    title = " "
  ) +
  theme_classic(base_size = 20) +
  theme()

ggsave(
  filename = "BMCA/plots/figure_3/MPs_composition.pdf",  # name of your PDF file
  plot = last_plot(),                     # or store your plot in a variable and use it
  width = 12.5,                             # width in inches
  height = 8,                             # height in inches
  units = "in"                             # units
)

# -----------------------------------------------------------------------------
sheets <- excel_sheets("BMCA/csv/all_MPS.xlsx")

tibble_list_mine <- set_names(
  map(sheets, ~ read_excel("BMCA/csv/all_MPS.xlsx", sheet = .x)),
  sheets
)

sheets <- excel_sheets("BMCA/csv/meta_programs_2025.xlsx")

tibble_list_3ca <- set_names(
  map(sheets, ~ read_excel("BMCA/csv/meta_programs_2025.xlsx", sheet = .x)),
  sheets
)

names(tibble_list_mine)
tibble_list_3ca$Epithelial <- NULL
names(tibble_list_3ca) <- c("Malignant",   "B_Plasma",     "Endothelial",
                            "Fibroblast", "Myeloid", "CD4", "CD8")


ct <- "Malignant"
for (ct in names(tibble_list_3ca)) {
  print(ct)
  # --- Generate plots ---
  mal1 <- tibble_list_3ca[[ct]]
  mal2 <- tibble_list_mine[[ct]]
  mal2$sc_lncRNA_MT <- NULL
  # Compute pairwise Jaccard
  jaccard_df <- expand.grid(
    Program1 = colnames(mal1),
    Program2 = colnames(mal2),
    stringsAsFactors = FALSE
  ) %>%
    rowwise() %>%
    mutate(
      Jaccard = {
        genes1 <- mal1[[Program1]]
        genes2 <- mal2[[Program2]]
        length(intersect(genes1, genes2)) / length(union(genes1, genes2))
      }
    ) %>%
    ungroup()

  best_matches <- jaccard_df %>%
    group_by(Program2) %>%
    # Select the row with the highest Jaccard value for each Program2
    slice_max(order_by = Jaccard, n = 1, with_ties = FALSE) %>%
    # Select only the two columns you want
    dplyr::select(Program2, Program1) %>%
    ungroup()

  # View the result
  print(best_matches)  # Order rows/cols by maximum value
  jaccard_df <- jaccard_df %>%
    filter(Program1 %in% unique(best_matches$Program1))
  # Applying the transformation to your dataframe
  jaccard_df <- jaccard_df %>%
    mutate(Program2 = Program2 %>%
             # 1. Capture the prefix (sc/sn) and the rest of the string
             # Regex: ^(s[cn])_ (start with sc or sn then underscore)
             #        (.*)$     (capture everything else)
             str_replace("^(s[cn])_(.*)$", "\\2 (\\1)") %>%

             # 2. Replace all remaining underscores with spaces
             str_replace_all("_", " ") %>%

             # 3. Replace all dots with hyphens
             str_replace_all("\\.", "-")
    )

  row_order <- jaccard_df %>%
    group_by(Program1) %>%
    summarise(max_val = max(Jaccard)) %>%
    arrange(desc(max_val)) %>%
    pull(Program1)

  col_order <- jaccard_df %>%
    group_by(Program2) %>%
    summarise(max_val = max(Jaccard)) %>%
    arrange(desc(max_val)) %>%
    pull(Program2)

  col_order <- c(
    "INFR MHC (sc)",
    "INFR (sn)",
    "CC (sc)",
    "Proteasomal degradation (sc)",
    "Stress (sc)",
    "Stress (sn)",
    "MHC-I (sc)",
    "PDAC related (sc)",
    "Hypoxia (sc)",
    "Hypoxia (sn)",
    "MYC targets-respiration (sc)",
    "Stress 2 (sc)",
    "secretory CRC enriched (sc)",
    "secretory CRC enriched (sn)",
    "CC (sn)",
    "EMT (sc)",
    "Skin Pigmentation (sc)",
    "Skin Pigmentation (sn)",
    "Chromatin Remodeling (sc)",
    "Hypoxia 2 (sn)",
    "Respiration 2 (sc)",
    "Respiration 1 (sc)",
    "UA1 (sn)",
    "EpiSen MHC (sc)",
    "Neuronal related 1 (sn)",
    "Neuronal related 2 (sn)",
    "Neuronal related 3 (sn)",
    "Neuronal related 2 (sc)"
  )
  col_order <- as.factor(col_order)
  # Factor for plotting
  jaccard_df$Program1 <- factor(jaccard_df$Program1, levels = rev(row_order))
  jaccard_df$Program2 <- factor(jaccard_df$Program2, levels = col_order)
  unique(jaccard_df$Program2)


  # Verify the changes
  unique(jaccard_df$Program2)

  rename_map <- c(
    "CC (sn)"                      = "Cell cycle (sn)",
    "Hypoxia (sn)"                 = "Hypoxia (sn)",
    "Neuronal related 1 (sn)"      = "MNR1 (sn)",
    "INFR (sn)"                    = "INFR (sn)",
    "Hypoxia 2 (sn)"               = "Hypoxia2 (sn)",
    "UA1 (sn)"                     = "UA1 (sn)",
    "Neuronal related 2 (sn)"      = "MNR2 (sn)",
    "Neuronal related 3 (sn)"      = "MNR3 (sn)",
    "secretory CRC enriched (sn)"  = "Epi. diff. 2 (sn)",
    "Stress (sn)"                  = "Stress (sn)",
    "Skin Pigmentation (sn)"       = "Skin Pig. (sn)",
    "Stress (sc)"                  = "Stress (sc)",
    "MHC-I (sc)"                   = "MHCI (sc)",
    "CC (sc)"                      = "Cell cycle (sc)",
    "INFR MHC (sc)"                = "INFR/MHC (sc)",
    "EpiSen MHC (sc)"              = "EpiSen/MHC (sc)",
    "Hypoxia (sc)"                 = "Hypoxia (sc)",
    "Proteasomal degradation (sc)" = "Proteasomal deg. (sc)",
    "Skin Pigmentation (sc)"       = "Skin Pig. (sc)",
    "MYC targets-respiration (sc)" = "MYC targets/resp. (sc)",
    "Chromatin Remodeling (sc)"    = "Chromatin Rel. (sc)",
    "Respiration 1 (sc)"           = "Resp. 1 (sc)",
    "PDAC related (sc)"            = "Epi. diff. 1 (sc)",   # <-- inferred, please confirm
    "Respiration 2 (sc)"           = "Resp. 2 (sc)",
    "Stress 2 (sc)"                = "Stress2 (sc)",
    "secretory CRC enriched (sc)"  = "Epi. diff. 2 (sc)",   # <-- inferred, please confirm
    "Neuronal related 2 (sc)"      = "MNR2 (sc)",
    "EMT (sc)"                     = "EMT (sc)"
  )

  new_order <- c(
    "INFR/MHC (sc)","INFR (sn)","Cell cycle (sc)","Proteasomal deg. (sc)",
    "Stress (sc)","Stress (sn)","MHCI (sc)","Epi. diff. 1 (sc)","Hypoxia (sc)",
    "Hypoxia (sn)","MYC targets/resp. (sc)","Stress2 (sc)","Epi. diff. 2 (sc)",
    "Epi. diff. 2 (sn)","Cell cycle (sn)","EMT (sc)","Skin Pig. (sc)","Skin Pig. (sn)",
    "Chromatin Rel. (sc)","Hypoxia2 (sn)","Resp. 2 (sc)","Resp. 1 (sc)","UA1 (sn)",
    "EpiSen/MHC (sc)","MNR1 (sn)","MNR2 (sn)","MNR3 (sn)","MNR2 (sc)"
  )

  # safety check BEFORE overwriting: any level in the data not covered by the map?
  setdiff(as.character(unique(jaccard_df$Program2)), names(rename_map))  # should be character(0)

  jaccard_df$Program2 <- factor(
    rename_map[as.character(jaccard_df$Program2)],
    levels = new_order
  )

  stopifnot(!anyNA(jaccard_df$Program2))   # catches any missed mapping

  unique(jaccard_df$Program1)

  prog1_rename <- c(
    "MP1 Cell Cycle - G2/M"         = "Cell cycle",
    "MP4 Chromatin"                 = "Chromatin",
    "MP5 Cell cycle single-nucleus" = "Cell cycle (sn)",
    "MP6 Stress 1"                  = "Stress",
    "MP7 Hypoxia"                   = "Hypoxia",
    "MP10 Proteasomal degradation"  = "Proteasomal deg.",
    "MP12 Protein maturation"       = "Protein mat.",
    "MP15 EMT-I"                    = "EMT-I",
    "MP19 EMT-V"                    = "EMT-V",
    "MP22 Interferon/MHC-II (I)"    = "INFR/MHCII",
    "MP28 NRF2 targets"             = "NRF2 targets",
    "MP30 Respiration 1"            = "Resp. 1",
    "MP31 Respiration 2"            = "Resp. 2",
    "MP34 Secreted II"              = "Secretory",
    "MP39 Oligo Progenitor"         = "Oligo prog.",
    "MP42 Glioma single-nucleus"    = "Glioma (sn)",
    "MP45 Skin-pigmentation"        = "Skin Pig.",
    "MP46 CRC stemness"             = "CRC stemness",
    "MP47 Colon-related"            = "Colon-related",
    "MP59 PDAC-related 1"           = "PDAC-related"
  )

  # apply, keeping the existing factor order
  old_levels <- levels(jaccard_df$Program1)
  jaccard_df$Program1 <- factor(
    prog1_rename[as.character(jaccard_df$Program1)],
    levels = prog1_rename[old_levels]
  )

  # verify nothing became NA
  setdiff(as.character(unique(jaccard_df$Program1)), prog1_rename)  # want character(0)

  # --- 1. build the plot with the colourbar title BELOW the bar ---
  p <- ggplot(jaccard_df, aes(x = Program2, y = Program1, fill = Jaccard)) +
    geom_tile(color = "white") +
    scale_fill_gradient2(
      low = "#7A3E9D", mid = "white", high = "#E68A00", midpoint = 0,
      name = "Jaccard\nsimilarity",
      guide = guide_colorbar(
        frame.colour    = "black",
        frame.linewidth = 0.4,
        ticks.colour    = "black",
        ticks.linewidth = 0.4
      )
    ) +
    labs(
      x = "This study states",
      y = "Most Similar Gavish et al. States"
    ) +
    theme_minimal(base_size = 6) +
    theme(
      axis.title.x    = element_text(size = 7),
      axis.title.y    = element_text(size = 7),
      axis.text.x     = element_text(angle = 45, hjust = 1),
      axis.text.y     = element_text(hjust = 1),
      panel.border    = element_rect(color = "black", fill = NA, linewidth = 1),
      legend.position = "right"
    )
  p

  ggsave("BMCA/plots/fig3/3ca_jaccard.pdf", p,
         width = 110, height = 80, units = "mm", dpi = 300)
# --- 2. extract the legend (bottom slot matches legend.position = "bottom") ---
  legend <- get_plot_component(p, "guide-box-bottom", return_all = TRUE)

  ggsave("BMCA/plots/fig3/3ca_jaccard_legend.pdf", legend,
         width = 4, height = 1.2, units = "in", device = cairo_pdf)

  # --- 3. save the plot WITHOUT the legend ---
  p_nolegend <- p + theme(legend.position = "none")

  ggsave("BMCA/plots/fig3/3ca_jaccard_heatmap.pdf", p_nolegend,
         width = 6, height = 5, units = "in", device = cairo_pdf)
  # --- Save to PDF ---
}

# -----------------------------------------------------------------------------
sheets <- excel_sheets("BMCA/csv/all_MPS.xlsx")

tibble_list_mine <- set_names(
  map(sheets, ~ read_excel("BMCA/csv/all_MPS.xlsx", sheet = .x)),
  sheets
)


# Check available databases if you want to explore others
dbs <- listEnrichrDbs()
dbs$libraryName
databases <- c("GO_Biological_Process_2025", "KEGG_2026", "MSigDB_Hallmark_2020",
               "GO_Cellular_Component_2025", "GO_Molecular_Function_2025",
               "CellMarker_2024", "ProteomicsDB_2020", "Human_Gene_Atlas",
               "Mouse_Gene_Atlas")
ct <- "Malignant"

# 1. Initialize an empty tibble with the correct column names but 0 rows
tib_of_enrich_per_ct <- tibble(
  term = character(),
  `-log10(fdr)` = numeric(),
  MP = character()
)

print(paste("Processing Cell Type:", ct))
mps <- tibble_list_mine[[ct]]

for (col in colnames(mps)) {
  print(paste("  Processing Column:", col))

  # Extract gene list
  sig <- unlist(mps[[col]])

  # Run Enrichr
  enrich_sig <- enrichr(sig, databases)

  # Optional: Plotting
  p <- plotEnrich(enrich_sig[[1]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value") +
    ggtitle(paste(ct, col, "Upregulated Pathways"))

  # 2. Extract the top hits
  result_tibble <- enrich_sig %>%
    map_df(as.data.frame, .id = "database") %>%
    mutate(
      `-log10(fdr)` = -log10(Adjusted.P.value),
      MP = col  # Create the column here where it's safe
    ) %>%
    arrange(desc(`-log10(fdr)`)) %>%
    slice(1:10) %>%
    dplyr::select(
      term = Term,
      database,
      `-log10(fdr)`,
      MP
    )

  # 3. Append to the summary table
  tib_of_enrich_per_ct <- bind_rows(tib_of_enrich_per_ct, result_tibble)
}


terms <-c(3, 11, 30, 31, 43 ,51, 68, 74, 81, 91, 104, 111, 124, 135, 142, 152,
          161, 171, 187, 193, 201, 211, 221, 236, 241, 251, 261)
tib_of_enrich_per_ct_b <- tib_of_enrich_per_ct[terms, ]
tib_of_enrich_per_ct_fil <- tib_of_enrich_per_ct[terms, ] %>%
  dplyr::select(term, database)
tib_of_enrich_per_ct_fil

tib_of_enrich_per_ct_fil <- tib_of_enrich_per_ct_fil %>%
  add_row(
    term = "Chromatin Remodeling (GO:0006338)",
    database = "GO_Biological_Process_2025"
  )


# 1. Initialize result container
final_plot_data <- tibble()

# Use unique databases from your filter list to speed up enrichr
target_databases <- unique(tib_of_enrich_per_ct_fil$database)

for (col in colnames(mps)) {
  print(col)
  sig_genes <- mps[[col]][!is.na(mps[[col]])] # Remove NAs

  # Run enrichr for this MP
  enrich_res <- enrichr(sig_genes, target_databases)

  # Combine results and filter for only your 28 target terms
  res_processed <- enrich_res %>%
    map_df(as.data.frame, .id = "database") %>%
    mutate(`-log10(fdr)` = -log10(Adjusted.P.value)) %>%
    # INNER JOIN acts as the filter: only keeps Term + Database pairs in your list
    inner_join(tib_of_enrich_per_ct_fil, by = c("Term" = "term", "database")) %>%
    mutate(MP = col) %>%
    dplyr::select(Term, database, `-log10(fdr)`, MP)

  final_plot_data <- bind_rows(final_plot_data, res_processed)
}

# Create a color mapping for the databases
db_colors <- c(
  "MSigDB_Hallmark_2020"       = "#E41A1C", # Red
  "CellMarker_2024"            = "#377EB8", # Blue
  "GO_Biological_Process_2025" = "#4DAF4A", # Green
  "Mouse_Gene_Atlas"           = "#FF7F00", # Orange
  "ProteomicsDB_2020"          = "#A65628", # Brown
  "KEGG_2026"                  = "#F781BF", # Pink
  "GO_Cellular_Component_2025" = "#00CED1"  # Dark Turquoise
)

# Sort plot data so terms appear grouped/ordered
final_plot_data <- final_plot_data %>%
  arrange(database, Term)

# Create the color vector for the axis labels
axis_colors <- db_colors[tib_of_enrich_per_ct_fil$database]

ggplot(final_plot_data, aes(x = Term, y = MP, size = `-log10(fdr)`, color = `-log10(fdr)`)) +
  geom_point() +
  scale_color_viridis_c(option = "magma", direction = -1) +
  theme_minimal() +
  labs(title = "Enrichment Analysis of Meta-Programs",
       x = "Gene Signature (Colored by Database)",
       y = "Meta-Program (MP)") +
  theme_minimal(base_size = 6) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold",
                               colour = axis_colors) # This applies the database colors
  )


final_plot_data <- final_plot_data %>%
  filter(`-log10(fdr)`  > 1.3)
# 1. Cap the values at 5
plot_df <- final_plot_data %>%
  mutate(capped_fdr = pmin(`-log10(fdr)`, 5))

# 2. Determine order for Y-axis (MPs)
# Sorted by the highest single enrichment value in that MP
mp_order <- plot_df %>%
  group_by(MP) %>%
  summarize(max_val = max(capped_fdr, na.rm = TRUE)) %>%
  arrange(max_val) %>% # ggplot plots from bottom to top, so highest is last
  pull(MP)

# 3. Determine order for X-axis (Terms)
# Sorted by the highest single enrichment value for that Term
term_order_df <- plot_df %>%
  group_by(Term, database) %>%
  summarize(max_val = max(capped_fdr, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(max_val))

# 4. Apply levels to factors
plot_df$MP <- factor(plot_df$MP, levels = mp_order)
plot_df$Term <- factor(plot_df$Term, levels = term_order_df$Term)

# 5. Create color vector for axis labels based on NEW term order
axis_colors <- db_colors[term_order_df$database]

# 6. Generate the Plot
ggplot(plot_df, aes(x = Term, y = MP)) +
  # This hidden geom creates the Database Legend
  geom_point(aes(fill = database), alpha = 0) +
  guides(fill = guide_legend(title = "Database Source",
                             override.aes = list(alpha = 1, size = 4))) +
  scale_fill_manual(values = db_colors) +

  # The actual data points
  geom_point(aes(size = capped_fdr, color = capped_fdr)) +
  scale_color_viridis_c(option = "magma", direction = -1, name = "-log10(FDR)\n(capped at 5)") +
  scale_size_continuous(range = c(1, 6), name = "Significance") +
  theme_minimal(base_size = 6) +
  labs(
    title = "Meta-Program Enrichment (Ordered by Signal Intensity)",
    x = "Enriched Term",
    y = "Meta-Program"
  ) +

  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
      colour = axis_colors # Applies database colors to labels
    ),
    legend.box = "vertical"
  )


# 1. Prepare Wide Matrix for Clustering
matrix_data <- plot_df %>%
  # Ensure only one value per pair by taking the maximum capped_fdr
  group_by(MP, Term) %>%
  summarize(capped_fdr = max(capped_fdr, na.rm = TRUE), .groups = "drop") %>%
  # Now pivot
  pivot_wider(
    names_from = Term,
    values_from = capped_fdr,
    values_fill = 0      # This fills missing entries with numeric 0
  ) %>%
  as.data.frame()

# Set MP as rownames
rownames(matrix_data) <- matrix_data$MP
matrix_data <- as.matrix(matrix_data[,-1])


# 2. Perform Hierarchical Clustering
# Compute distance (Euclidean) and clustering (Complete linkage)
row_hc <- hclust(dist(matrix_data))
col_hc <- hclust(dist(t(matrix_data)))

# 3. Extract the Order
ordered_mps <- rownames(matrix_data)[row_hc$order]
ordered_terms <- colnames(matrix_data)[col_hc$order]

# 4. Re-apply orders to the Plotting Dataframe
plot_df$MP <- factor(plot_df$MP, levels = ordered_mps)
plot_df$Term <- factor(plot_df$Term, levels = ordered_terms)

# 5. Re-sync Axis Colors for the new clustered order
# We need to find the database for each term in the NEW order
label_info_clustered <- plot_df %>%
  distinct(Term, database) %>%
  slice(match(ordered_terms, Term))

# 1. Sync Axis Colors (same as before)
label_info_clustered <- plot_df %>%
  distinct(Term, database) %>%
  slice(match(ordered_terms, Term))

# 1. Ensure label_info_clustered is exactly synced with factor levels
label_info_clustered <- tibble(Term = levels(plot_df$Term)) %>%
  left_join(plot_df %>% distinct(Term, database), by = "Term")

axis_colors_clustered <- db_colors[label_info_clustered$database]
# 2. Generate the Plot
p <- ggplot(plot_df, aes(x = Term, y = MP)) +
  # This creates the legend. We map 'database' to 'fill'
  # and use a point shape that supports filling (shape 21)
  geom_point(aes(fill = database), color = "transparent", size = 0.01) +

  # The actual data points
  geom_point(aes(color = `-log10(fdr)`, size = `-log10(fdr)`)) +

  # Define the manual fill for the legend
  scale_fill_manual(
    values = db_colors,
    name = "Database Source"
  ) +

  scale_color_gradientn(
    colors = c("#F39F5A", "#AE445A", "#662549", "#451952"),
    name = "-log10(FDR)"
  ) +

  # Crucial: Override the legend aesthetics so they are visible
  guides(
    fill = guide_legend(override.aes = list(color = db_colors, size = 4, alpha = 1)),
    color = guide_colorbar(order = 1)
  ) +
  theme_minimal(base_size = 5) +
  labs(
    title = " ",
    x = "Enriched Term",
    y = "Meta-Program"
  ) +
  theme(
    axis.text.x = element_text(
      angle = 45, hjust = 1, vjust = 1,
      colour = axis_colors_clustered
    ),
    legend.position = "right"
  )
p
ggsave(
  filename = "BMCA/plots/fig3/metaprogram_heatmap_enrichment.pdf",
  plot = p,
  width = 140,       # Width in inches
  height = 150,        # Height in inches
  units = "mm",      # Can be "in", "cm", or "mm"
  dpi = 300          # High resolution
)

for (ct in names(tibble_list_mine)) {
  # 1. Initialize an empty tibble with the correct column names but 0 rows
  tib_of_enrich_per_ct <- tibble(
    term = character(),
    `-log10(fdr)` = numeric(),
    MP = character()
  )

  print(paste("Processing Cell Type:", ct))
  mps <- tibble_list_mine[[ct]]

  for (col in colnames(mps)) {
    print(paste("  Processing Column:", col))

    # Extract gene list
    sig <- unlist(mps[[col]])

    # Run Enrichr
    enrich_sig <- enrichr(sig, databases)

    # Optional: Plotting
    p <- plotEnrich(enrich_sig[[1]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value") +
      ggtitle(paste(ct, col, "Upregulated Pathways"))

    # 2. Extract the top hit across ALL databases
    result_tibble <- enrich_sig %>%
      map_df(as.data.frame) %>%
      mutate(`-log10(fdr)` = -log10(Adjusted.P.value)) %>%
      arrange(desc(`-log10(fdr)`)) %>%
      slice(1) %>%
      dplyr::select(term = Term, `-log10(fdr)`) %>%
      mutate(MP = col) # Add the MP column here

    # 3. Append to the summary table
    tib_of_enrich_per_ct <- bind_rows(tib_of_enrich_per_ct, result_tibble)
  }

  # Final summary for this cell type
  print(tib_of_enrich_per_ct)
}

# -----------------------------------------------------------------------------
ecotypes_tib <- readRDS("BMCA/data/ecotypes_tib_proccesed_filtered.RDS")
metadata_states <- readRDS("BMCA/data/metadata_states_proccesed_filtered.RDS")

metadata_states <- metadata_states %>%
  filter(primary_site %in% c("colorectal", "breast", "NSCLC", "melanoma") &
           location == "Intracranial") %>%
  mutate(
    primary_site = recode(primary_site,
                          "melanoma" = "Melanoma",
                          "breast" = "Breast",
                          "colorectal" = "Colorectal",
                          "NSCLC" = "NSCLC"
    ))

patients <- intersect(metadata_states$patient, ecotypes_tib$patient)

metadata_states <- metadata_states %>%
  filter(patient %in% patients) %>%
  dplyr::select(patient, primary_site, seq_tech)

ecotypes_tib <- ecotypes_tib %>%
  filter(patient %in% patients) %>%
  dplyr::select(patient,
                Malignant_Neuronal_related_1,
                Malignant_Neuronal_related_2,
                Malignant_Neuronal_related_3)

ecotypes_tib <- left_join(ecotypes_tib, metadata_states, by = "patient")

ecotypes_long <- ecotypes_tib %>%
  pivot_longer(
    cols = starts_with("Malignant_Neuronal"),
    names_to = "Signature",
    values_to = "Score"
  ) %>%
  # Remove NAs so they don't break the plot
  filter(!is.na(Score)) %>%
  mutate(
    primary_site = recode(primary_site,
                       "melanoma" = "Melanoma",
                       "breast" = "Breast",
                       "colorectal" = "Colorectal",
                       "NSCLC" = "NSCLC"
    ))


# ── Statistics: each cancer type vs. the other three pooled ──────────────────
sites <- unique(ecotypes_long$primary_site)
sigs  <- unique(ecotypes_long$Signature)

wilcox_results <- map_dfr(sigs, function(sig) {
  map_dfr(sites, function(site) {

    dat     <- ecotypes_long %>% filter(Signature == sig)
    group_a <- dat %>% filter(primary_site == site) %>% pull(Score)
    group_b <- dat %>% pull(Score)

    if (length(group_a) < 2 | length(group_b) < 2) return(NULL)

    wt <- wilcox.test(group_a, group_b, exact = FALSE, alternative = "greater")

    tibble(
      Signature    = sig,
      primary_site = site,
      p_value      = wt$p.value
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
    ),
    Signature = recode(Signature,
                       "Malignant_Neuronal_related_1" = "MNR 1",
                       "Malignant_Neuronal_related_2" = "MNR 2",
                       "Malignant_Neuronal_related_3" = "MNR 3"
    )
  )

# ── FDR + recode BEFORE join ──────────────────────────────────────────────────
wilcox_results <- wilcox_results %>%
  mutate(
    p_adj = p.adjust(p_value, method = "BH"),
    stars = case_when(
      p_adj < 0.001 ~ "**",
      p_adj < 0.01  ~ "*",
      p_adj < 0.05  ~ "",
      TRUE          ~ ""
    ),
    Signature = recode(Signature,
                       "Malignant_Neuronal_related_1" = "MNR 1",
                       "Malignant_Neuronal_related_2" = "MNR 2",
                       "Malignant_Neuronal_related_3" = "MNR 3"
    )
  )

# ── Heatmap data ──────────────────────────────────────────────────────────────
heatmap_data <- ecotypes_long %>%
  mutate(
    Signature = recode(Signature,
                       "Malignant_Neuronal_related_1" = "MNR 1",
                       "Malignant_Neuronal_related_2" = "MNR 2",
                       "Malignant_Neuronal_related_3" = "MNR 3"
    )
  ) %>%
  group_by(primary_site, Signature) %>%
  summarise(mean_frac = mean(Score, na.rm = TRUE), .groups = "drop") %>%
  left_join(
    wilcox_results %>% dplyr::select(primary_site, Signature, stars),
    by = c("primary_site", "Signature")
  ) %>%
  mutate(
    primary_site = factor(primary_site,
                          levels = c("Colorectal", "Breast", "NSCLC", "Melanoma")),
    tile_label = ifelse(
      stars == "",
      paste0(round(mean_frac, 1), "%"),
      paste0(round(mean_frac, 1), "% ", stars)
    )
  )

# ── Verify join ───────────────────────────────────────────────────────────────
print(heatmap_data %>% dplyr::select(primary_site, Signature, mean_frac, stars))

# ── Plot ──────────────────────────────────────────────────────────────────────
p <- ggplot(heatmap_data,
            aes(x = Signature, y = primary_site, fill = mean_frac)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(
    aes(
      label = tile_label,
      color = mean_frac > 5
    ),
    size       = 4,
    fontface   = "bold",
    lineheight = 0.9
  ) +
  scale_fill_gradientn(
    colors = c("#EED9B9", "#D53E0F", "#5E0006"),
    limits = c(0, NA),
    name   = "Mean % of\nmalignant cells"
  ) +
  scale_color_manual(
    values = c("TRUE" = "white", "FALSE" = "#5E0006"),
    guide  = "none"
  ) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 8, base_family = "Helvetica") +
  theme(
    axis.text.x      = element_text(size = 10),
    axis.text.y      = element_text(size = 10),
    panel.grid       = element_blank(),
    legend.title     = element_text(size = 10),
    legend.text      = element_text(size = 8),
    legend.key.width = unit(0.35, "cm"),
    plot.margin      = margin(4, 4, 4, 4)
  )
p

metadata_mnr <- metadata %>%
  filter(cell_type == "Malignant") %>%
  dplyr::select(CellID, patient, state)

# ── Step 1: count cells per patient × MNR state ───────────────────────────────
mnr_states <- c("Malignant_Neuronal_related_1",
                "Malignant_Neuronal_related_2",
                "Malignant_Neuronal_related_3")

cell_counts <- metadata_mnr %>%
  filter(state %in% mnr_states) %>%
  count(patient, state, name = "n_cells")

# ── Step 2: get all patients (so those with 0 cells are included) ─────────────
all_patients <- metadata_mnr %>%
  distinct(patient)

patient_meta <- ecotypes_tib %>%
  distinct(patient, primary_site)

# complete: every patient × every MNR state, fill missing with 0
cell_counts_complete <- all_patients %>%
  cross_join(tibble(state = mnr_states)) %>%
  left_join(cell_counts, by = c("patient", "state")) %>%
  mutate(n_cells = replace_na(n_cells, 0)) %>%
  left_join(patient_meta, by = "patient")

# ── Step 3: % of patients with > 20 cells per cancer type × state ────────────
heatmap_input <- cell_counts_complete %>%
  filter(!is.na(primary_site),
         primary_site %in% c("Colorectal", "Breast", "NSCLC", "Melanoma")) %>%
  group_by(primary_site, state) %>%
  summarise(
    pct_above20 = mean(n_cells > 20, na.rm = TRUE) * 100,
    .groups     = "drop"
  ) %>%
  mutate(
    Signature = recode(state,
                       "Malignant_Neuronal_related_1" = "MNR 1",
                       "Malignant_Neuronal_related_2" = "MNR 2",
                       "Malignant_Neuronal_related_3" = "MNR 3"
    )
  )

# ── Step 4: Wilcoxon — each cancer type vs. all pooled, one-sided ─────────────
sigs  <- unique(cell_counts_complete$state)
sites <- c("Colorectal", "Breast", "NSCLC", "Melanoma")

wilcox_results <- map_dfr(sigs, function(sig) {
  map_dfr(sites, function(site) {

    dat     <- cell_counts_complete %>%
      filter(state == sig,
             primary_site %in% sites)
    group_a <- dat %>% filter(primary_site == site) %>% pull(n_cells)
    group_b <- dat %>% pull(n_cells)

    if (length(group_a) < 2 | length(group_b) < 2) return(NULL)

    wt <- wilcox.test(group_a, group_b, exact = FALSE, alternative = "greater")

    tibble(
      state        = sig,
      primary_site = site,
      p_value      = wt$p.value
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
    ),
    Signature = recode(state,
                       "Malignant_Neuronal_related_1" = "MNR 1",
                       "Malignant_Neuronal_related_2" = "MNR 2",
                       "Malignant_Neuronal_related_3" = "MNR 3"
    )
  )

# ── Step 5: join stats into heatmap data ──────────────────────────────────────
heatmap_data <- heatmap_input %>%
  left_join(
    wilcox_results %>% dplyr::select(primary_site, Signature, stars),
    by = c("primary_site", "Signature")
  ) %>%
  mutate(
    primary_site = factor(primary_site,
                          levels = c("Colorectal", "Breast", "NSCLC", "Melanoma")),
    tile_label = ifelse(
      stars == "",
      paste0(round(pct_above20, 1), "%"),
      paste0(round(pct_above20, 1), "%", stars)
    )
  )

# ── Step 6: plot ──────────────────────────────────────────────────────────────
p <- ggplot(heatmap_data,
            aes(x = Signature, y = primary_site, fill = pct_above20)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(
    aes(
      label = tile_label,
      color = pct_above20 > 40
    ),
    size       = 4,
    fontface   = "bold",
    lineheight = 0.9
  ) +
  scale_fill_gradientn(
    colors = c("#E0D9D9", "#5A9690", "#2F5755"),
    limits = c(0, NA),
    name   = "% patients\n>20 cells in state"
  ) +
  scale_color_manual(
    values = c("TRUE" = "white", "FALSE" = "black"),
    guide  = "none"
  ) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 8, base_family = "Helvetica") +
  theme(
    axis.text.x      = element_text(size = 10),
    axis.text.y      = element_text(size = 10),
    panel.grid       = element_blank(),
    legend.title     = element_text(size = 10),
    legend.text      = element_text(size = 8),
    legend.key.width = unit(0.35, "cm"),
    plot.margin      = margin(4, 4, 4, 4)
  )
p

# -----------------------------------------------------------------------------

# ── seq_tech lookup per patient ───────────────────────────────────────────────
patient_meta <- ecotypes_tib %>%
  distinct(patient, primary_site, seq_tech)

# ═══════════════════════════════════════════════════════════════════════════════
# HEATMAP 1 — mean fraction score, faceted by seq_tech
# ═══════════════════════════════════════════════════════════════════════════════

# ── Statistics ────────────────────────────────────────────────────────────────
sites <- unique(ecotypes_long$primary_site)
sigs  <- unique(ecotypes_long$Signature)
techs <- c("sn", "sc")

wilcox_results_h1 <- map_dfr(techs, function(tech) {
  map_dfr(sigs, function(sig) {
    map_dfr(sites, function(site) {

      dat     <- ecotypes_long %>%
        filter(Signature == sig, seq_tech == tech)
      group_a <- dat %>% filter(primary_site == site) %>% pull(Score)
      group_b <- dat %>% pull(Score)

      if (length(group_a) < 2 | length(group_b) < 2) return(NULL)

      wt <- wilcox.test(group_a, group_b, exact = FALSE, alternative = "greater")

      tibble(
        seq_tech     = tech,
        Signature    = sig,
        primary_site = site,
        p_value      = wt$p.value
      )
    })
  })
}) %>%
  mutate(
    p_adj = p.adjust(p_value, method = "BH"),
    stars = case_when(
      p_adj < 0.001 ~ "***",
      p_adj < 0.01  ~ "**",
      p_adj < 0.05  ~ "*",
      TRUE          ~ ""
    ),
    Signature = recode(Signature,
                       "Malignant_Neuronal_related_1" = "MNR 1",
                       "Malignant_Neuronal_related_2" = "MNR 2",
                       "Malignant_Neuronal_related_3" = "MNR 3"
    )
  )

# ── Heatmap data ──────────────────────────────────────────────────────────────
heatmap_data_h1 <- ecotypes_long %>%
  mutate(
    Signature = recode(Signature,
                       "Malignant_Neuronal_related_1" = "MNR 1",
                       "Malignant_Neuronal_related_2" = "MNR 2",
                       "Malignant_Neuronal_related_3" = "MNR 3"
    )
  ) %>%
  group_by(primary_site, Signature, seq_tech) %>%
  summarise(mean_frac = mean(Score, na.rm = TRUE), .groups = "drop") %>%
  left_join(
    wilcox_results_h1 %>% dplyr::select(primary_site, Signature, seq_tech, stars),
    by = c("primary_site", "Signature", "seq_tech")
  ) %>%
  mutate(
    primary_site = factor(primary_site,
                          levels = c("Colorectal", "Breast", "NSCLC", "Melanoma")),
    seq_tech = recode(seq_tech, "sn" = "snRNA-seq", "sc" = "scRNA-seq"),
    tile_label = ifelse(
      stars == "",
      paste0(round(mean_frac, 1), "%"),
      paste0(round(mean_frac, 1), "% ", stars)
    )
  )

# ── Plot ──────────────────────────────────────────────────────────────────────
p_h1 <- ggplot(heatmap_data_h1,
               aes(x = Signature, y = primary_site, fill = mean_frac)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(
    aes(
      label = tile_label,
      color = mean_frac > 5
    ),
    size       = 4,
    fontface   = "bold",
    lineheight = 0.9
  ) +
  facet_wrap(~ seq_tech, ncol = 2) +
  scale_fill_gradientn(
    colors = c("#EED9B9", "#D53E0F", "#5E0006"),
    limits = c(0, NA),
    name   = "Mean % of\nmalignant cells"
  ) +
  scale_color_manual(
    values = c("TRUE" = "white", "FALSE" = "#5E0006"),
    guide  = "none"
  ) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 8, base_family = "Helvetica") +
  theme(
    strip.text       = element_text(size = 10, face = "bold"),
    axis.text.x      = element_text(size = 10),
    axis.text.y      = element_text(size = 10),
    panel.grid       = element_blank(),
    legend.title     = element_text(size = 10),
    legend.text      = element_text(size = 8),
    legend.key.width = unit(0.35, "cm"),
    plot.margin      = margin(4, 4, 4, 4)
  )
p_h1

# ── Two separate plots for Heatmap 1 ──────────────────────────────────────────

make_h1_plot <- function(data, tech_label) {
  ggplot(data,
         aes(x = Signature, y = primary_site, fill = mean_frac)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(
      aes(
        label = tile_label,
        color = mean_frac > 5
      ),
      size       = 4,
      lineheight = 0.9
    ) +
    scale_fill_gradientn(
      colors = c("#EED9B9", "#D53E0F", "#5E0006"),
      limits = c(0, NA),
      name   = "Mean % of\nmalignant cells"
    ) +
    scale_color_manual(
      values = c("TRUE" = "white", "FALSE" = "#5E0006"),
      guide  = "none"
    ) +
    labs(x = NULL, y = NULL, title = tech_label) +
    theme_minimal(base_size = 7, base_family = "Helvetica") +
    theme(
      panel.grid       = element_blank(),
      legend.key.width = unit(0.35, "cm"),
      plot.margin      = margin(4, 4, 4, 4)
    )
}

p_h1_sn <- make_h1_plot(
  heatmap_data_h1 %>% filter(seq_tech == "snRNA-seq"),
  ""
)

p_h1_sc <- make_h1_plot(
  heatmap_data_h1 %>% filter(seq_tech == "scRNA-seq"),
  ""
)

p_h1_sn
p_h1_sc

ggsave("BMCA/plots/fig3/heatmap1_sn.pdf", p_h1_sn, width = 5, height = 4, units = "in", dpi = 300)
ggsave("BMCA/plots/fig3/heatmap1_sc.pdf", p_h1_sc, width = 5, height = 4, units = "in", dpi = 300)
# ═══════════════════════════════════════════════════════════════════════════════
# HEATMAP 2 — % patients with >20 cells, faceted by seq_tech
# ═══════════════════════════════════════════════════════════════════════════════

# ── cell_counts_complete now includes seq_tech ────────────────────────────────
cell_counts <- metadata_mnr %>%
  filter(state %in% mnr_states) %>%
  count(patient, state, name = "n_cells")

cell_counts_complete <- metadata_mnr %>%
  distinct(patient) %>%
  cross_join(tibble(state = mnr_states)) %>%
  left_join(cell_counts, by = c("patient", "state")) %>%
  mutate(n_cells = replace_na(n_cells, 0)) %>%
  left_join(patient_meta, by = "patient")   # brings in primary_site + seq_tech

# ── Statistics ────────────────────────────────────────────────────────────────
wilcox_results_h2 <- map_dfr(techs, function(tech) {
  map_dfr(mnr_states, function(sig) {
    map_dfr(sites, function(site) {

      dat     <- cell_counts_complete %>%
        filter(state == sig, seq_tech == tech, primary_site %in% sites)
      group_a <- dat %>% filter(primary_site == site) %>% pull(n_cells)
      group_b <- dat %>% pull(n_cells)

      if (length(group_a) < 2 | length(group_b) < 2) return(NULL)

      wt <- wilcox.test(group_a, group_b, exact = FALSE, alternative = "greater")

      tibble(
        seq_tech     = tech,
        state        = sig,
        primary_site = site,
        p_value      = wt$p.value
      )
    })
  })
}) %>%
  mutate(
    p_adj = p.adjust(p_value, method = "BH"),
    stars = case_when(
      p_adj < 0.001 ~ "***",
      p_adj < 0.01  ~ "**",
      p_adj < 0.05  ~ "*",
      TRUE          ~ ""
    ),
    Signature = recode(state,
                       "Malignant_Neuronal_related_1" = "MNR 1",
                       "Malignant_Neuronal_related_2" = "MNR 2",
                       "Malignant_Neuronal_related_3" = "MNR 3"
    )
  )

# ── Heatmap data ──────────────────────────────────────────────────────────────
heatmap_data_h2 <- cell_counts_complete %>%
  filter(!is.na(primary_site), !is.na(seq_tech),
         primary_site %in% sites) %>%
  mutate(
    Signature = recode(state,
                       "Malignant_Neuronal_related_1" = "MNR 1",
                       "Malignant_Neuronal_related_2" = "MNR 2",
                       "Malignant_Neuronal_related_3" = "MNR 3"
    )
  ) %>%
  group_by(primary_site, Signature, seq_tech) %>%
  summarise(
    pct_above20 = mean(n_cells > 20, na.rm = TRUE) * 100,
    .groups     = "drop"
  ) %>%
  left_join(
    wilcox_results_h2 %>% dplyr::select(primary_site, Signature, seq_tech, stars),
    by = c("primary_site", "Signature", "seq_tech")
  ) %>%
  mutate(
    primary_site = factor(primary_site,
                          levels = c("Colorectal", "Breast", "NSCLC", "Melanoma")),
    seq_tech = recode(seq_tech, "sn" = "snRNA-seq", "sc" = "scRNA-seq"),
    tile_label = ifelse(
      stars == "",
      paste0(round(pct_above20, 1), "%"),
      paste0(round(pct_above20, 1), "% ", stars)
    )
  )

# ── Plot ──────────────────────────────────────────────────────────────────────
p_h2 <- ggplot(heatmap_data_h2,
               aes(x = Signature, y = primary_site, fill = pct_above20)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(
    aes(
      label = tile_label,
      color = pct_above20 > 40
    ),
    size       = 4,
    fontface   = "bold",
    lineheight = 0.9
  ) +
  facet_wrap(~ seq_tech, ncol = 2) +
  scale_fill_gradientn(
    colors = c("#E0D9D9", "#5A9690", "#2F5755"),
    limits = c(0, NA),
    name   = "% patients\n>20 cells in state"
  ) +
  scale_color_manual(
    values = c("TRUE" = "white", "FALSE" = "black"),
    guide  = "none"
  ) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 7, base_family = "Helvetica") +
  theme(
    panel.grid       = element_blank(),
    legend.key.width = unit(0.35, "cm"),
    plot.margin      = margin(4, 4, 4, 4)
  )

p_h2
ggsave("BMCA/plots/fig3/heatmap_cells_mnrs.pdf", p_h2, width = 160, height = 100,
       units = "mm",dpi = 300)

# -----------------------------------------------------------------------------


ecotypes_tib <- readRDS("BMCA/data/ecotypes_tib_proccesed_filtered.RDS")
metadata_states <- readRDS("BMCA/data/metadata_states_proccesed_filtered.RDS")

metadata_states <- metadata_states %>%
  filter(primary_site %in% c("colorectal", "breast", "NSCLC", "melanoma") &
           location == "Intracranial") %>%
  mutate(
    primary_site = recode(primary_site,
                          "melanoma" = "Melanoma",
                          "breast" = "Breast",
                          "colorectal" = "Colorectal",
                          "NSCLC" = "NSCLC"
    ))

patients <- intersect(metadata_states$patient, ecotypes_tib$patient)

metadata_states <- metadata_states %>%
  filter(patient %in% patients) %>%
  dplyr::select(patient, primary_site, seq_tech)

ecotypes_tib <- ecotypes_tib %>%
  filter(patient %in% patients) %>%
  dplyr::select(patient,
                Malignant_Neuronal_related_1,
                Malignant_Neuronal_related_2,
                Malignant_Neuronal_related_3)

ecotypes_tib <- left_join(ecotypes_tib, metadata_states, by = "patient")


dens_df <- bind_rows(
  ecotypes_tib %>%
    filter(primary_site == "Melanoma", seq_tech == "sn",
           !is.na(Malignant_Neuronal_related_1)) %>%
    transmute(value = Malignant_Neuronal_related_1, program = "MNR 1 — Melanoma"),
  ecotypes_tib %>%
    filter(primary_site == "Breast", seq_tech == "sn",
           !is.na(Malignant_Neuronal_related_2)) %>%
    transmute(value = Malignant_Neuronal_related_2, program = "MNR 2 — Breast")
)

p_overlay <- ggplot(dens_df, aes(x = value, colour = program, fill = program)) +
  geom_density(aes(y = after_stat(density)), linewidth = 0.8, alpha = 0.15) +
  scale_colour_manual(values = c("MNR 1 — Melanoma" = "#D53E0F",
                                 "MNR 2 — Breast"   = "#5A9690")) +
  scale_fill_manual(values   = c("MNR 1 — Melanoma" = "#D53E0F",
                                 "MNR 2 — Breast"   = "#5A9690")) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "MNR cells (% of malignant cells)", y = "Density",
       colour = NULL, fill = NULL, title = "") +
  theme_minimal(base_size = 7, base_family = "Helvetica") +
  theme(panel.grid.minor = element_blank(),
        legend.position = "bottom")
p_overlay

ggsave("BMCA/plots/fig3/mnrs_histo.pdf", p_overlay , width = 80, height = 80,
       units = "mm", dpi = 300)


# ============================================================
# Neuron gene expression: Intracranial vs Extracranial
# Per-patient pseudobulk differential analysis
#   - Separated version  : tested within each seq_tech (sn / sc)
#   - Integrated version : average of the sn & sc results
# ============================================================


# ------------------------------------------------------------
# Parameters
# ------------------------------------------------------------

MIN_CELLS_PER_PATIENT <- 50
MIN_PATIENTS_PER_SIDE <- 3

gene_groups <- list(
  neurotrophins  = c("NGF",    "BDNF",   "NTF3",   "NTF4",  "NTRK1", "NTRK2", "NTRK3"),
  cgrp_signaling = c("CALCA",  "CALCB",  "CALCRL", "RAMP1"),
  glutamate_ampa = c("GRIA1",  "GRIA2",  "GRIA3",  "GRIA4"),
  gaba_a         = c("GABRA1", "GABRB1", "GABRG2"),
  muscarinic     = c("CHRM1",  "CHRM2",  "CHRM3",  "CHRM5"),
  sensory_misc   = c("TRPV1",  "NLGN3",  "SEMA4F")
)

gene_cols <- unlist(gene_groups, use.names = FALSE)

# Integrated version excludes the sensory_misc group
gene_groups_integrated <- gene_groups[setdiff(names(gene_groups), "sensory_misc")]

# ------------------------------------------------------------
# Publication-ready facet / axis labels
# ------------------------------------------------------------

gene_group_labels <- c(
  neurotrophins  = "Neurotrophins",
  cgrp_signaling = "CGRP signaling",
  glutamate_ampa = "Glutamate (AMPA)",
  gaba_a         = "GABA",
  muscarinic     = "Muscarinic",
  sensory_misc   = "Sensory / misc"
)

# EDIT to match your actual seq_tech values (used for the separated-plot rows).
# Anything not listed here passes through unchanged.
seq_tech_labels <- c(
  sn = "Single-nucleus",
  sc = "Single-cell"
)

# Display names for cell types (applied at plot time only)
cell_type_relabel <- c(T_NK = "T/NK", B_Plasma = "B/Plasma")

# ============================================================
# Load and combine the per-dataset gene tables
# ============================================================

neuron_gene_list <- readRDS("BMCA/data/neuron_gene_list.RDS")

combined_neurons <- neuron_gene_list %>%
  rbindlist(fill = TRUE, use.names = TRUE) %>%   # faster than bind_rows
  as_tibble()

# ------------------------------------------------------------
# Clean metadata and join onto expression
# ------------------------------------------------------------

metadata <- metadata %>%
  filter(!cell_type %in% c("Astrocyte", "Oligodendrocyte",
                           "CD4", "CD8", "GD", "NK", "Mesenchimal",
                           "Other"))

cells <- intersect(metadata$CellID, combined_neurons$CellID)

metadata_n <- metadata %>%
  filter(CellID %in% cells) %>%
  dplyr::select(CellID, patient, primary_site, cell_type,
                location, study, seq_tech, Donor, state)

combined_neurons <- combined_neurons %>%
  filter(CellID %in% cells) %>%
  left_join(metadata_n, by = "CellID")


# ============================================================
# Optional exploratory block: intracranial fibroblasts (NTRK)
# ============================================================

combined_neurons_IC <- combined_neurons %>%
  filter(location == "Intracranial", cell_type == "Fibroblast") %>%
  dplyr::select(CellID, NTRK2, NTRK3, patient, study, seq_tech)

plot_ntrk_hist <- function(df, gene) {
  ggplot(df, aes(x = .data[[gene]], fill = seq_tech)) +
    geom_histogram(bins = 30, alpha = 0.7, color = "white") +
    scale_fill_brewer(palette = "Set2") +
    labs(title = paste0("Distribution of ", gene, " Expression"),
         x = paste0(gene, " Expression"), y = "Count", fill = "seq_tech") +
    theme_minimal()
}

plot_ntrk_hist(combined_neurons_IC, "NTRK3")
plot_ntrk_hist(combined_neurons_IC, "NTRK2")

result <- combined_neurons_IC %>%
  group_by(patient) %>%
  summarise(
    NTRK2_mean = log2(mean(2^NTRK2, na.rm = TRUE)),
    NTRK3_mean = log2(mean(2^NTRK3, na.rm = TRUE)),
    .groups = "drop"
  )

NTRK_CAFs <- combined_neurons_IC$CellID[
  combined_neurons_IC$NTRK2 > 4 | combined_neurons_IC$NTRK3 > 4
]

saveRDS(NTRK_CAFs, "BMCA/data/NTKR_CAFs.RDS")
saveRDS(result,    "BMCA/data/patient_means_fibro.RDS")

# ============================================================
# Center each gene within (cell_type x seq_tech)
#   Removes the per-platform baseline so sn and sc become
#   comparable. Uses ALL cells in the group (both locations).
# ============================================================

centered_neurons <- combined_neurons %>%
  group_by(cell_type, seq_tech) %>%
  mutate(across(all_of(gene_cols), ~ .x - mean(.x, na.rm = TRUE))) %>%
  ungroup()


# ============================================================
# Per-patient pseudobulk for one seq_tech / cell_type / location
#   Excludes patients with < MIN_CELLS_PER_PATIENT cells.
# ============================================================

build_pb <- function(df, st, ct, loc) {
  sub <- df %>% filter(seq_tech == st, cell_type == ct, location == loc)
  if (nrow(sub) == 0) return(NULL)

  valid_patients <- sub %>%
    count(patient, name = "n_cells") %>%
    filter(n_cells >= MIN_CELLS_PER_PATIENT) %>%
    pull(patient)

  if (length(valid_patients) == 0) return(NULL)

  sub %>%
    filter(patient %in% valid_patients) %>%
    group_by(patient) %>%
    summarise(across(all_of(gene_cols), ~ mean(.x, na.rm = TRUE)),
              .groups = "drop")
}

# ============================================================
# SEPARATED: Wilcoxon Intra vs Extra, within each seq_tech
# ============================================================

compute_stats_separated <- function(df, genes = gene_cols) {
  seq_techs  <- unique(df$seq_tech)
  cell_types <- unique(df$cell_type)
  results    <- list()

  for (st in seq_techs) {
    message("Processing seq_tech: ", st)

    for (ct in cell_types) {
      message("  cell_type: ", ct)

      g_intra <- build_pb(df, st, ct, "Intracranial")
      g_extra <- build_pb(df, st, ct, "Extracranial")

      if (is.null(g_intra) || is.null(g_extra)) next
      if (nrow(g_intra) < MIN_PATIENTS_PER_SIDE ||
          nrow(g_extra) < MIN_PATIENTS_PER_SIDE) next

      for (gene in genes) {
        vi <- g_intra[[gene]]
        ve <- g_extra[[gene]]

        if (all(vi == 0) && all(ve == 0)) {
          results[[length(results) + 1]] <- tibble(
            seq_tech = st, cell_type = ct, gene = gene,
            p_value = NA_real_, log2FC = 0
          )
          next
        }

        wt <- tryCatch(wilcox.test(vi, ve, exact = FALSE),
                       error = function(e) list(p.value = NA_real_))

        results[[length(results) + 1]] <- tibble(
          seq_tech = st, cell_type = ct, gene = gene,
          p_value = wt$p.value,
          log2FC  = mean(vi, na.rm = TRUE) - mean(ve, na.rm = TRUE)
        )
      }
    }
  }
  bind_rows(results)
}

# ============================================================
# Run separated analysis + FDR within each seq_tech
# ============================================================

set.seed(67)

stats_separated <- compute_stats_separated(centered_neurons) %>%
  group_by(seq_tech) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  ungroup() %>%
  mutate(neg_log10_p = pmin(-log10(p_adj), 10))

# ============================================================
# INTEGRATED: simply average the sn & sc results
#   log2FC = mean of the two seq_tech log2FCs
#   p_adj  = mean of the two per-seq_tech BH-corrected p-values
#   sensory_misc genes are dropped. (No further FDR correction,
#   since we are averaging values that are already adjusted.)
# ============================================================

stats_integrated <- stats_separated %>%
  filter(!gene %in% gene_groups$sensory_misc) %>%
  group_by(cell_type, gene) %>%
  summarise(
    log2FC = mean(log2FC, na.rm = TRUE),
    p_adj  = mean(p_adj,  na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(neg_log10_p = pmin(-log10(p_adj), 10))

# ============================================================
# Publication-ready dot plot
# ============================================================

theme_pub <- function(base_size = 11) {
  theme_bw(base_size = base_size) +
    theme(
      axis.text.x       = element_text(angle = 45, hjust = 1, face = "italic"),
      axis.text.y       = element_text(size = base_size),
      strip.background  = element_rect(fill = "grey92", color = NA),
      strip.text        = element_text(face = "bold", size = base_size - 1),
      panel.grid.minor  = element_blank(),
      panel.grid.major  = element_line(color = "grey93"),
      panel.spacing     = unit(0.4, "lines"),
      legend.position   = "right",
      legend.key.height = unit(0.8, "lines")
    )
}

make_dotplot <- function(stats_df, groups, facet_rows = TRUE) {
  group_map <- imap(groups, ~ tibble(gene = .x, Group = .y)) %>%
    bind_rows() %>%
    mutate(
      Group = factor(Group, levels = names(groups)),
      gene  = factor(gene,  levels = unlist(groups, use.names = FALSE))
    )

  plot_df <- stats_df %>%
    filter(!is.na(p_adj)) %>%
    inner_join(group_map, by = "gene") %>%   # drops genes not in `groups`
    mutate(
      cell_type = recode(as.character(cell_type), !!!cell_type_relabel),
      cell_type = factor(cell_type)
    )

  facet <- if (facet_rows) {
    facet_grid(
      seq_tech ~ Group, scales = "free_x", space = "free_x",
      labeller = labeller(
        seq_tech = as_labeller(seq_tech_labels, default = label_value),
        Group    = as_labeller(gene_group_labels)
      )
    )
  } else {
    facet_grid(
      . ~ Group, scales = "free_x", space = "free_x",
      labeller = labeller(Group = as_labeller(gene_group_labels))
    )
  }

  ggplot(plot_df, aes(x = gene, y = cell_type)) +
    geom_point(aes(color = log2FC, size = neg_log10_p), alpha = 0.9) +
    scale_color_gradient2(
      low = "#2166AC", mid = "white", high = "#B2182B",
      midpoint = 0, na.value = "grey90",
      name = "log2FC\n(BrMs - Primary)"
    ) +
    scale_size_continuous(
      range = c(1, 7),
      name  = expression(-log[10] ~ italic(p)[adj])
    ) +
    facet +
    labs(x = "Gene", y = "Cell type") +
    theme_pub()
}

# Separated plot (rows = sn / sc), all gene groups
p_separated <- make_dotplot(stats_separated, gene_groups, facet_rows = TRUE)

# Integrated plot (single row, sensory_misc removed)
p_integrated <- make_dotplot(stats_integrated, gene_groups_integrated, facet_rows = FALSE)

p_separated
p_integrated

ggsave("BMCA/plots/fig3/dotplot_separated.pdf",  p_separated,  width = 12, height = 4)
ggsave("BMCA/plots/fig3/dotplot_integrated.pdf", p_integrated, width = 10, height = 4)

