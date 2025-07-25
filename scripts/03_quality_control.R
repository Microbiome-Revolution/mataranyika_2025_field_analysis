# =============================================================================
# 03_quality_control.R
# Data filtering and quality control
# =============================================================================

source("scripts/02_data_import.R")

# Quality control filtering
cat("Applying quality filters...\n")

ps_clean <- ps_project %>%
  # Remove problematic taxa
  subset_taxa(order != "unclassified_Cyanobacteriia") %>% 
  subset_taxa(class != "unclassified_Cyanobacteria") %>%  
  subset_taxa(phylum != "unclassified_Bacteria") %>%
  subset_taxa(domain != "unclassified_Root") %>% 
  subset_taxa(domain != "Eukaryota") %>%
  subset_taxa(order != "Chloroplast") %>%
  subset_taxa(family != "Mitochondria")

cat("Taxa before filtering:", ntaxa(ps_project), "\n")
cat("Taxa after filtering:", ntaxa(ps_clean), "\n")

# Convert to matrix for rarefaction analysis
otu_matrix <- as(otu_table(ps_clean), "matrix")
otu_matrix <- otu_matrix[rowSums(otu_matrix) > 0, colSums(otu_matrix) > 0]

# Generate rarefaction curves
cat("Generating rarefaction curves...\n")
tryCatch({
  RCurve <- rarecurve(otu_matrix, step = 500, tidy = FALSE, label = FALSE)
  
  # Convert to tidy format
  RCurve_tidy <- imap_dfr(RCurve, ~ {
    data.frame(
      Sample = rownames(otu_matrix)[.y],
      Reads = attr(.x, "Subsample"),
      OTUs = as.vector(.x),
      stringsAsFactors = FALSE
    )
  })
  
  # Plot rarefaction curves
  p_rarefaction <- ggplot(RCurve_tidy, aes(x = Reads, y = OTUs)) +
    geom_line(aes(group = Sample), alpha = 0.3, linewidth = 0.5) +
    geom_vline(xintercept = 1250, linetype = "dashed", color = "red") +
    labs(
      x = "Sequencing Depth",
      y = "Observed OTUs",
      title = "Rarefaction Curves",
      subtitle = paste("Samples:", nrow(otu_matrix), "| OTUs:", ncol(otu_matrix))
    ) +
    scale_x_continuous(labels = scales::comma) +
    scale_y_continuous(labels = scales::comma)
  
  ggsave(file.path(FIGURES, "rarefaction_curves.png"), p_rarefaction, 
         width = 10, height = 6, dpi = 300)
  
}, error = function(e) {
  cat("Rarefaction curve generation failed:", e$message, "\n")
})

# Plot sequencing depth distribution
seq_depths <- sample_sums(ps_clean)
depth_df <- data.frame(Sample = names(seq_depths), Reads = seq_depths)

p_depth_dist <- ggplot(depth_df, aes(x = Reads)) +
  geom_histogram(binwidth = 100, fill = "blue", alpha = 0.5, color = "black") +
  geom_vline(xintercept = 750, linetype = "dashed", color = "red") +
  geom_vline(xintercept = 1000, linetype = "dashed", color = "orange") +
  geom_vline(xintercept = 1250, linetype = "dashed", color = "green") +
  labs(
    title = "Distribution of Sequencing Depths",
    x = "Sequencing Depth (Total Reads per Sample)",
    y = "Number of Samples"
  )

ggsave(file.path(FIGURES, "sequencing_depth_distribution.png"), p_depth_dist, 
       width = 10, height = 6, dpi = 300)

# Filter samples by minimum read count
min_reads <- 1250
ps_filtered <- prune_samples(sample_sums(ps_clean) >= min_reads, ps_clean)

cat("Samples before depth filtering:", nsamples(ps_clean), "\n")
cat("Samples after depth filtering:", nsamples(ps_filtered), "\n")

# Save filtered data
saveRDS(ps_filtered, file.path(DATA_PROCESSED, "phyloseq_filtered.rds"))

cat("Quality control complete!\n")
