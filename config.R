# =============================================================================
# config.R
# Project configuration and settings
# =============================================================================

# Project metadata
PROJECT_NAME <- "Mataranyika 2025 Field Analysis"
PROJECT_VERSION <- "1.0.0"
ANALYSIS_DATE <- Sys.Date()

# Analysis parameters
RAREFACTION_DEPTH <- 1250
RAREFACTION_ITERATIONS <- 50
MIN_SAMPLE_READS <- 1250
ALPHA_THRESHOLD <- 0.05

# File paths (relative to project root)
PATHS <- list(
  data_raw = "data/raw",
  data_processed = "data/processed", 
  outputs = "outputs",
  figures = "outputs/figures",
  tables = "outputs/tables",
  scripts = "scripts",
  functions = "functions"
)

# Expected input files
INPUT_FILES <- list(
  sequence_data = "R1SeqDataTable.RDS",
  metadata = "R1_field_data1.csv",  # Primary metadata file
  metadata_alt = "metadata.csv"     # Alternative metadata file
)

# Output file naming conventions
OUTPUT_PREFIX <- paste0("microbiome_analysis_", format(Sys.Date(), "%Y%m%d"))

# Color palettes for consistent plotting
COLOR_PALETTES <- list(
  treatment = c("#E69F00", "#56B4E9", "#009E73", "#F0E442"),
  takeall = c("#d62728", "#2ca02c"),
  tillage = c("#1f77b4", "#ff7f0e", "#2ca02c")
)

# Print configuration
cat("Project Configuration Loaded:\n")
cat("- Name:", PROJECT_NAME, "\n")
cat("- Version:", PROJECT_VERSION, "\n") 
cat("- Analysis Date:", as.character(ANALYSIS_DATE), "\n")
cat("- Rarefaction Depth:", RAREFACTION_DEPTH, "\n")
