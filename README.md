# Brain Metastasis Cell Atlas (BMCA) — analysis code

Analysis and figure-generation code for:

> **A brain metastasis cell atlas reveals multicellular adaptation to the neural microenvironment**
> Feigin et al.

The atlas integrates single-cell and single-nucleus RNA-seq across brain
metastases and matched extracranial tumours, comparing the malignant and
tumour-microenvironment compartments between the two sites.

This repository contains the analysis scripts, a simulated demo dataset, and
instructions for reproducing the figures from the processed data objects.

---

## 1. System requirements

### Software dependencies

| | Version |
|---|---|
| R | ≥ 4.2.0 — the exact tested version is recorded in `sessionInfo.txt` |
| RStudio (optional) | ≥ 2022.07 |

R packages, all loaded at the top of each script:

**CRAN** — `cowplot`, `data.table`, `dendextend`, `doParallel`, `dplyr`,
`enrichR`, `forcats`, `foreach`, `ggdendro`, `ggh4x`, `ggnewscale`, `ggplot2`,
`ggpubr`, `ggraph`, `ggrastr`, `ggrepel`, `Hmisc`, `igraph`, `janitor`,
`lubridate`, `Matrix`, `NMF`, `openxlsx`, `patchwork`, `pheatmap`, `purrr`,
`RColorBrewer`, `readr`, `readxl`, `reshape2`, `rlang`, `rstatix`, `scales`,
`stringr`, `survival`, `survminer`, `tibble`, `tidyr`, `tidytext`, `tidyverse`,
`viridis`

**Bioconductor** — `biomaRt`

**GitHub** — `scalop` (Tirosh lab; `remotes::install_github("jlaffy/scalop")`)

Exact tested versions of every package are listed in `package_versions.csv` and
`sessionInfo.txt`. Regenerate them at any time with:

```bash
Rscript R/capture_session_info.R
```

### Operating systems

The code is platform independent and has been run on:

- macOS
- Linux, on the institutional cluster

The exact OS release the final figures were produced on is recorded in the
`Running under:` line of `sessionInfo.txt`. Windows is expected to work but has
not been tested.

> **Note on case sensitivity.** The scripts are run on a case-sensitive
> filesystem on Linux. Keep input filenames exactly as listed in each script
> header.

### Non-standard hardware

None. No GPU or specialised hardware is required.

- The **demo** runs in under 1 GB of RAM on any standard laptop.
- The **full analyses** load the complete atlas objects and pseudobulk
  matrices, and are memory-bound rather than CPU-bound; they were run on a
  standard compute node rather than a laptop. A few steps use
  `doParallel`/`foreach` and benefit from multiple cores but do not require
  them.

---

## 2. Installation guide

No compilation and no installation of the analysis code itself is required —
the scripts are run directly from the repository.

```bash
git clone https://github.com/nisanfeigin/BMCA.git
cd BMCA
```

Install the R dependencies:

```r
install.packages(c(
  "tidyverse", "data.table", "Matrix", "Hmisc", "janitor", "openxlsx",
  "readxl", "reshape2", "rlang", "RColorBrewer", "viridis", "scales",
  "ggplot2", "ggpubr", "ggrepel", "ggh4x", "ggnewscale", "ggrastr",
  "ggdendro", "dendextend", "cowplot", "patchwork", "pheatmap", "tidytext",
  "rstatix", "survival", "survminer", "igraph", "ggraph", "enrichR",
  "NMF", "doParallel", "foreach"
))

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("biomaRt")

if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
remotes::install_github("jlaffy/scalop")
```

### Typical install time

On a standard desktop with a normal broadband connection:

| | Time |
|---|---|
| macOS / Windows (pre-built binaries) | ~10–15 minutes |
| Linux (compiled from source) | ~30–45 minutes |

`NMF` and `survminer` account for most of the compile time on Linux.

---

## 3. Demo

A small **simulated** dataset is included so that the code can be run end to end
without access to the controlled-access patient data. It reproduces the
structure of the real objects (column names, cohort design, pseudobulk matrix
shape) and contains a planted differential-expression signal, so the demo
exercises the same statistical machinery as the manuscript analyses. It is not
derived from any patient sample.

### Instructions to run on the demo data

From the repository root:

```bash
Rscript demo/make_demo_data.R   # writes demo/data/  (~0.2 MB)
Rscript demo/run_demo.R         # writes demo/output/
```

The demo needs only `dplyr`, `tidyr`, `tibble`, `readr` and `ggplot2`, so it can
be run before the full dependency list is installed.

The demo reproduces the two core analysis patterns used throughout the paper:

1. **Cell-type composition, intracranial vs. extracranial** — per-patient
   proportions, Wilcoxon rank-sum test, Benjamini–Hochberg correction. This is
   the method behind the Figure 4 composition panels.
2. **Pseudobulk differential expression** — CPM normalisation, log2 transform,
   per-gene Wilcoxon test and log2 fold change, volcano plot. This is the
   method behind the Figure 1 stromal DGE panels.

### Expected output

Four files in `demo/output/`:

| File | Contents |
|---|---|
| `demo_celltype_composition.pdf` | Faceted boxplots of per-patient cell-type fractions by location |
| `demo_celltype_composition_stats.csv` | Per cell type: medians, log2FC, p, BH-adjusted p |
| `demo_volcano_fibroblast.pdf` | Volcano plot of the pseudobulk comparison |
| `demo_pseudobulk_dge.csv` | Per-gene log2FC, p, BH-adjusted p, direction call |

Reference copies of all four are in `demo/expected_output/`, together with the
console log in `demo/expected_output/console_output.txt`. Both demo scripts are
seeded, so the numbers below are reproduced exactly:

```
Cell-type composition:
  cell types tested:       7
  significant (BH < 0.05): 2          # Fibroblast and Myeloid

Pseudobulk DGE (fibroblasts):
  genes tested:            1997
  up intracranial:         59
  up extracranial:         61
  recovery of planted intracranial genes: 98%
  recovery of planted extracranial genes: 100%
```

The recovery percentages are a self-check: the demo reports how many of the
genes it planted as differentially expressed were recovered by the analysis.

### Expected run time

On a standard desktop: **~1 second** for `make_demo_data.R` and **~3 seconds**
for `run_demo.R`.

---

## 4. Instructions for use

### Repository structure

```
.
├── R/
│   ├── config.R                        # project root, directory layout, output dirs
│   ├── capture_session_info.R          # records package versions
│   ├── cell_state_tables.R             # shared upstream: abundance matrices
│   ├── fibroblast_pericyte_analysis.R  # Figure 1 — stromal compartment
│   ├── figure_3.R                      # Figure 3 — malignant MPs and MNR programs
│   ├── figure_4.R                      # Figure 4 — TME composition and cell states
│   └── figure_5.R                      # Figure 5 — ecotypes, networks, survival
├── demo/
│   ├── make_demo_data.R                # simulated dataset generator
│   ├── run_demo.R                      # end-to-end demo analysis
│   └── expected_output/                # reference output of the demo
├── BMCA/
│   ├── data/                           # .RDS inputs (not tracked — see Data)
│   ├── csv/                            # signatures and clinical tables (not tracked)
│   └── plots/                          # figure panels, created at run time
├── BMCA.Rproj
└── README.md
```

| Script | Figure panels | Original filename |
|---|---|---|
| `cell_state_tables.R` | upstream tables (no panels) | `states_tibs.R` |
| `fibroblast_pericyte_analysis.R` | `plots/fig1/` | `22_perycites_scores.R` |
| `figure_3.R` | `plots/fig3/` | `Figure_3_2.R` |
| `figure_4.R` | `plots/fig4/` | `figure_4_1.R` |
| `figure_5.R` | `plots/fig5/` | `Figure_5_4.R` |

### Data

The `.RDS` and `.csv` inputs named in each script header are not included here
because of size and patient-privacy restrictions. Processed objects are
available from **[add repository / Zenodo DOI]**; raw sequencing data are at
**[add GEO / EGA accessions]**. Place them so the layout matches `BMCA/data/`
and `BMCA/csv/` above.

Every script header lists the exact files it reads. Paths built at run time are
shown with a placeholder, e.g. `BMCA/data/UMI_pb_<cell_type>.RDS`.

### Running on your own data

The scripts expect the following object schemas:

- **`metadata_all_studies.rds`** — one row per cell, with columns `CellID`,
  `patient`, `study`, `primary_site`, `location` (`Intracranial` /
  `Extracranial`), `seq_tech` (`sc` / `sn`), `cell_type`, `state`.
- **`celltype_patient_matrix_proccesed.RDS`** — one row per patient, one column
  per cell type or cell state, values are cell counts.
- **`UMI_pb_<cell_type>.RDS`** — a named list of gene × sample count matrices,
  one element per cohort (e.g. `lin_sn`).

`demo/make_demo_data.R` is the shortest description of these schemas: it builds
each object from scratch in about 100 lines.

To run on your own cohort, reshape your data into those objects, place them
under `BMCA/data/`, and set the project root:

```bash
export BMCA_ROOT=/path/to/BMCA-code
```

or open `BMCA.Rproj` in RStudio. `R/config.R` resolves the root, checks the
input directories exist, and creates the plot output directories.

### Reproducing the manuscript figures

Run `R/cell_state_tables.R` first — it builds the abundance matrices the figure
scripts consume. The four figure scripts are independent of one another and can
be run in any order.

```r
source("R/cell_state_tables.R")
source("R/fibroblast_pericyte_analysis.R")
source("R/figure_3.R")
source("R/figure_4.R")
source("R/figure_5.R")
```

Each script is written to be read and run section by section in the order it
appears; individual panels are written to `BMCA/plots/` as PDFs and assembled
into the final figure layouts externally.

---

## Citation

If you use this code, please cite the manuscript above.

## License

**[Add a license before making the repository public — MIT or BSD-3-Clause are
the usual choices for analysis code accompanying a paper.]**

## Contact

Nisan Feigin — Tirosh Lab, Weizmann Institute of Science.
