# =============================================================================
# migrate_data.R
# Script to move existing data files to the new structure
# =============================================================================

cat("Migrating data files to new structure...\n")

# Create source and destination paths
old_files <- list(
  " metadata.csv" = "data/raw/metadata.csv",
  "sample_field_specific_environmental_data.xlsx" = "data/raw/field_environmental_data.xlsx",
  "sample_field_specific_environmental_data_FULL copy.csv" = "data/raw/field_environmental_data_full.csv"
)

# Move files
for (old_path in names(old_files)) {
  new_path <- old_files[[old_path]]
  
  if (file.exists(old_path)) {
    # Create directory if needed
    dir.create(dirname(new_path), recursive = TRUE, showWarnings = FALSE)
    
    # Move and rename file
    file.rename(old_path, new_path)
    cat("Moved:", old_path, "->", new_path, "\n")
  } else {
    cat("File not found:", old_path, "\n")
  }
}

# Note: You'll need to manually copy R1SeqDataTable.RDS to data/raw/ if it exists
cat("\nNote: Please copy R1SeqDataTable.RDS to data/raw/ if it exists in your system\n")
cat("Data migration complete!\n")
