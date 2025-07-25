# =============================================================================
# 02_data_import.R
# Data loading and initial cleaning
# =============================================================================

source("scripts/01_setup_environment.R")

# Convert SMP output to phyloseq object
cat("Loading sequence data...\n")
ps_raw <- SimpleMetaPackage::SeqDataTable2Phyloseq(
  file.path(DATA_RAW, "R1SeqDataTable.RDS")
)

# Load the metadata CSV
cat("Loading metadata...\n")
metadata_file <- file.path(DATA_RAW, "R1_field_data1.csv")
if (!file.exists(metadata_file)) {
  # Try alternative names if main file doesn't exist
  metadata_file <- list.files(DATA_RAW, pattern = "field.*data.*csv", full.names = TRUE)[1]
}

metadata <- read.csv(metadata_file, stringsAsFactors = FALSE)

# Clean metadata
metadata <- metadata %>% 
  filter(SampleID != "") %>%
  # Remove any trailing/leading whitespace
  mutate(across(where(is.character), str_trim))

# Ensure row names match sample IDs
rownames(metadata) <- metadata$SampleID
metadata$SampleID <- NULL  # Remove redundant column

# Replace the current sam_data with the new metadata
sample_data(ps_raw) <- sample_data(metadata)

# Subset data for analysis (adjust based on your project column)
ps_project <- subset_samples(ps_raw, project %in% c("sample", "farmer_kit"))
# Alternative: keep only main samples
# ps_project <- subset_samples(ps_raw, project == "sample")

# Transpose OTU table if needed
otu_table(ps_project) <- t(otu_table(ps_project))

# Save intermediate results
saveRDS(ps_project, file.path(DATA_PROCESSED, "phyloseq_raw.rds"))
fwrite(metadata, file.path(DATA_PROCESSED, "metadata_clean.csv"))

cat("Data import complete!\n")
cat("Samples:", nsamples(ps_project), "\n")
cat("Taxa:", ntaxa(ps_project), "\n")
