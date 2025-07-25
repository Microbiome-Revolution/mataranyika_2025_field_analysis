# Mataranyika et al. 2025 - 16S rRNA Field Analysis

## Project Overview
This repository contains the analysis pipeline for 16S rRNA gene sequencing data from agricultural field samples, investigating the relationship between tillage methods, fertilizer use, and soil microbial communities in relation to take-all disease.

## Directory Structure
```
├── README.md                     # Project documentation
├── data/                        # Data files
│   ├── raw/                     # Raw input data (read-only)
│   │   ├── metadata.csv
│   │   ├── R1SeqDataTable.RDS
│   │   └── field_data.csv
│   ├── processed/               # Processed/cleaned data
│   └── results/                 # Analysis outputs
├── scripts/                     # Analysis scripts
│   ├── 01_setup_environment.R   # Package installation & loading
│   ├── 02_data_import.R         # Data loading and cleaning
│   ├── 03_quality_control.R     # Filtering and QC
│   ├── 04_alpha_diversity.R     # Alpha diversity analysis
│   ├── 05_beta_diversity.R      # Ordination and distance analysis
│   ├── 06_taxonomic_composition.R # Relative abundance plots
│   ├── 07_statistical_tests.R   # ANOVA, PERMANOVA, etc.
│   ├── 08_differential_abundance.R # DESeq2 analysis
│   └── 09_spatial_analysis.R    # Mantel tests, distance decay
├── functions/                   # Custom R functions
│   ├── plotting_functions.R
│   ├── analysis_functions.R
│   └── utility_functions.R
├── outputs/                     # Generated files
│   ├── figures/                 # Plots and visualizations
│   ├── tables/                  # Summary statistics
│   └── reports/                 # Automated reports
├── .gitignore                   # Git ignore file
└── renv.lock                    # R environment (if using renv)
```

## Quick Start
1. Clone this repository
2. Open R/RStudio in the project directory
3. Run scripts in numerical order (01-09)
4. Check `outputs/` for results

## Data Description
- **Study Design**: Agricultural fields with different tillage methods and fertilizer treatments
- **Response Variable**: Take-all disease presence/absence
- **Sequencing**: 16S rRNA gene amplicon sequencing
- **Samples**: [Add number] soil samples from [Add number] fields

## Key Analyses
- Alpha diversity (Shannon, Hill numbers)
- Beta diversity (NMDS, PERMANOVA)
- Taxonomic composition
- Differential abundance analysis
- Spatial autocorrelation (Mantel tests)

## Dependencies
See `01_setup_environment.R` for required packages.

## Citation
[Add citation when published]
