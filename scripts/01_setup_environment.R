# =============================================================================
# 01_setup_environment.R
# Package installation and loading for 16S microbiome analysis
# =============================================================================

# Install all required packages
required_packages <- c(
  "phyloseq", "ggplot2", "vegan", "dplyr", "tidyr",
  "SimpleMetaPackage", "metagMisc", "purrr", "scales",
  "igraph", "hillR", "tibble", "data.table", "geosphere", 
  "lubridate", "DESeq2", "ComplexHeatmap"
)

# Install missing packages
install.packages(setdiff(required_packages, rownames(installed.packages())))

# Install Bioconductor packages
if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install(c("phyloseq", "DESeq2", "ComplexHeatmap"))

# Install GitHub packages
if (!require("devtools", quietly = TRUE)) {
  install.packages("devtools")
}
devtools::install_github("J-Cos/SimpleMetaPackage")

if (!require("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
remotes::install_github("vmikk/metagMisc")

# Load Required Libraries
library(phyloseq)
library(ggplot2)
library(vegan)
library(dplyr)
library(tidyr)
library(SimpleMetaPackage)
library(metagMisc)
library(purrr)
library(scales)
library(igraph)
library(hillR)
library(tibble)
library(data.table)
library(geosphere)
library(lubridate)
library(DESeq2)
library(ComplexHeatmap)

# Set up project paths
PROJECT_ROOT <- here::here()
DATA_RAW <- file.path(PROJECT_ROOT, "data", "raw")
DATA_PROCESSED <- file.path(PROJECT_ROOT, "data", "processed")
OUTPUTS <- file.path(PROJECT_ROOT, "outputs")
FIGURES <- file.path(OUTPUTS, "figures")
TABLES <- file.path(OUTPUTS, "tables")

# Create output directories if they don't exist
dir.create(FIGURES, recursive = TRUE, showWarnings = FALSE)
dir.create(TABLES, recursive = TRUE, showWarnings = FALSE)

# Set ggplot theme
theme_set(theme_minimal())

cat("Environment setup complete!\n")
cat("Project root:", PROJECT_ROOT, "\n")
