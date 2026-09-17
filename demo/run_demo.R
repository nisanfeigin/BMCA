# =============================================================================
# Demo analysis
# =============================================================================
#
# Runs the two core analysis patterns used throughout the manuscript on the
# small simulated dataset produced by demo/make_demo_data.R:
#
#   1. Cell-type composition, intracranial vs. extracranial
#      per-patient proportions -> Wilcoxon rank-sum -> Benjamini-Hochberg
#      (the method behind Figure 4 composition panels)
#
#   2. Pseudobulk differential expression, intracranial vs. extracranial
#      CPM -> log2 -> per-gene Wilcoxon + log2 fold change -> volcano
#      (the method behind the Figure 1 stromal DGE panels)
#
# Run from the project root:
#   Rscript demo/make_demo_data.R
#   Rscript demo/run_demo.R
#
# Expected run time: under 1 minute on a standard laptop.
# Expected output:   4 files in demo/output/ (see README.md)
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(ggplot2)
  library(readr)
})

t_start <- Sys.time()
set.seed(1)

dir.create("demo/output", recursive = TRUE, showWarnings = FALSE)

stopifnot(file.exists("demo/data/demo_metadata.rds"))
metadata <- readRDS("demo/data/demo_metadata.rds")
pb_list  <- readRDS("demo/data/demo_UMI_pb_Fibroblast.rds")
truth    <- readRDS("demo/data/demo_truth.rds")

# =============================================================================
# 1. Cell-type composition: intracranial vs. extracranial
# =============================================================================

comp <- metadata %>%
  count(patient, location, seq_tech, cell_type, name = "n") %>%
  group_by(patient) %>%
  mutate(frac = 100 * n / sum(n)) %>%
  ungroup()

# Brain-resident types are only recovered by single-nucleus sequencing, so they
# are tested within the sn cohort only, as in the manuscript.
brain_resident <- c("Astrocyte", "Oligodendrocyte", "Neuron")

comp_test <- comp %>%
  filter(!(cell_type %in% brain_resident) | seq_tech == "sn")

stats <- comp_test %>%
  group_by(cell_type) %>%
  filter(n_distinct(location) == 2) %>%
  summarise(
    n_ic    = sum(location == "Intracranial"),
    n_ec    = sum(location == "Extracranial"),
    med_ic  = median(frac[location == "Intracranial"]),
    med_ec  = median(frac[location == "Extracranial"]),
    p_value = wilcox.test(frac ~ location, exact = FALSE)$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    log2FC = log2((med_ic + 0.01) / (med_ec + 0.01)),
    p_adj  = p.adjust(p_value, method = "BH")
  ) %>%
  arrange(p_adj)

write_csv(stats, "demo/output/demo_celltype_composition_stats.csv")

p_comp <- comp_test %>%
  ggplot(aes(x = location, y = frac, fill = location)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, linewidth = 0.3) +
  geom_jitter(width = 0.15, size = 0.6, alpha = 0.6) +
  facet_wrap(~ cell_type, scales = "free_y", nrow = 2) +
  scale_fill_manual(values = c(Intracranial = "#3B6FB6", Extracranial = "#C9772E")) +
  labs(x = NULL, y = "Fraction of cells per patient (%)",
       title = "Demo: cell-type composition by tumour location",
       subtitle = "Simulated data - not patient-derived") +
  theme_classic(base_size = 9) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 30, hjust = 1),
        strip.background = element_blank())

ggsave("demo/output/demo_celltype_composition.pdf", p_comp,
       width = 9, height = 4.5, dpi = 300)

# =============================================================================
# 2. Pseudobulk differential expression (fibroblasts)
# =============================================================================

cts <- pb_list[["lin_sn"]]

# CPM normalise each sample, then log2 - the normalisation used throughout
cpm  <- t(t(cts) / colSums(cts)) * 1e6
lcpm <- log2(cpm + 1)

# keep expressed genes
keep <- rowMeans(lcpm) > 1
lcpm <- lcpm[keep, , drop = FALSE]

sample_loc <- metadata %>%
  distinct(patient, location) %>%
  arrange(match(patient, colnames(lcpm)))
stopifnot(identical(sample_loc$patient, colnames(lcpm)))

is_ic <- sample_loc$location == "Intracranial"

dge <- tibble(
  gene   = rownames(lcpm),
  mean_ic = rowMeans(lcpm[, is_ic, drop = FALSE]),
  mean_ec = rowMeans(lcpm[, !is_ic, drop = FALSE])
) %>%
  mutate(
    log2FC  = mean_ic - mean_ec,
    p_value = apply(lcpm, 1, function(g)
      wilcox.test(g[is_ic], g[!is_ic], exact = FALSE)$p.value),
    p_adj   = p.adjust(p_value, method = "BH"),
    direction = case_when(
      p_adj < 0.05 & log2FC >  1 ~ "Up intracranial",
      p_adj < 0.05 & log2FC < -1 ~ "Up extracranial",
      TRUE                       ~ "n.s."
    )
  ) %>%
  arrange(p_adj)

write_csv(dge, "demo/output/demo_pseudobulk_dge.csv")

p_volcano <- ggplot(dge, aes(log2FC, -log10(p_adj), colour = direction)) +
  geom_point(size = 0.7, alpha = 0.7) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", linewidth = 0.3) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.3) +
  scale_colour_manual(values = c("Up intracranial" = "#3B6FB6",
                                 "Up extracranial" = "#C9772E",
                                 "n.s."            = "grey75")) +
  labs(x = "log2 fold change (intracranial / extracranial)",
       y = "-log10 adjusted p",
       colour = NULL,
       title = "Demo: pseudobulk differential expression, fibroblasts",
       subtitle = "Simulated data - not patient-derived") +
  theme_classic(base_size = 9) +
  theme(strip.background = element_blank())

ggsave("demo/output/demo_volcano_fibroblast.pdf", p_volcano,
       width = 5.5, height = 4.5, dpi = 300)

# =============================================================================
# Summary and self-check
# =============================================================================

hits_ic <- dge$gene[dge$direction == "Up intracranial"]
hits_ec <- dge$gene[dge$direction == "Up extracranial"]

recall_ic <- mean(truth$up_intracranial %in% hits_ic)
recall_ec <- mean(truth$up_extracranial %in% hits_ec)

elapsed <- round(as.numeric(difftime(Sys.time(), t_start, units = "secs")), 1)

cat("\n=============================================================\n")
cat("Demo complete in", elapsed, "seconds\n")
cat("=============================================================\n\n")
cat("Cell-type composition:\n")
cat("  cell types tested:      ", nrow(stats), "\n")
cat("  significant (BH < 0.05):", sum(stats$p_adj < 0.05), "\n\n")
cat("Pseudobulk DGE (fibroblasts):\n")
cat("  genes tested:           ", nrow(dge), "\n")
cat("  up intracranial:        ", length(hits_ic), "\n")
cat("  up extracranial:        ", length(hits_ec), "\n")
cat("  recovery of planted intracranial genes:",
    sprintf("%.0f%%", 100 * recall_ic), "\n")
cat("  recovery of planted extracranial genes:",
    sprintf("%.0f%%", 100 * recall_ec), "\n\n")
cat("Output written to demo/output/:\n")
cat("  demo_celltype_composition.pdf\n")
cat("  demo_celltype_composition_stats.csv\n")
cat("  demo_volcano_fibroblast.pdf\n")
cat("  demo_pseudobulk_dge.csv\n")
