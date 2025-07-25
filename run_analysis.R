# =============================================================================
# run_analysis.R
# Master script to run the complete 16S microbiome analysis pipeline
# =============================================================================

# Set working directory to project root
if (!require(here)) install.packages("here")
setwd(here::here())

cat("Starting 16S microbiome analysis pipeline...\n")
cat("=" , rep("=", 50), "\n", sep = "")

# Track execution time
start_time <- Sys.time()

# Define scripts to run in order
scripts <- c(
  "01_setup_environment.R",
  "02_data_import.R", 
  "03_quality_control.R",
  "04_alpha_diversity.R",
  "05_beta_diversity.R",
  "06_taxonomic_composition.R",
  "07_statistical_tests.R",
  "08_differential_abundance.R",
  "09_spatial_analysis.R"
)

# Run each script
for (script in scripts) {
  script_path <- file.path("scripts", script)
  
  if (file.exists(script_path)) {
    cat("\nRunning:", script, "\n")
    cat("-", rep("-", 40), "\n", sep = "")
    
    tryCatch({
      source(script_path)
      cat("✓ Completed:", script, "\n")
    }, error = function(e) {
      cat("✗ Error in", script, ":", e$message, "\n")
    })
  } else {
    cat("⚠ Script not found:", script, "\n")
  }
}

# Summary
end_time <- Sys.time()
total_time <- end_time - start_time

cat("\n", rep("=", 60), "\n", sep = "")
cat("Analysis pipeline completed!\n")
cat("Total execution time:", round(total_time, 2), attr(total_time, "units"), "\n")
cat("Check 'outputs/' directory for results\n")
cat(rep("=", 60), "\n", sep = "")
