# =============================================================================
# 05_beta_diversity.R
# Beta diversity analysis including NMDS and PERMANOVA
# =============================================================================

source("scripts/04_alpha_diversity.R")

cat("Starting beta diversity analysis...\n")

# Compute Bray-Curtis distance
bray_dist <- phyloseq::distance(ps_filtered, method = "bray")

# Run NMDS
cat("Running NMDS ordination...\n")
nmds_result <- metaMDS(bray_dist, k = 2, trymax = 100)
nmds_scores <- as.data.frame(scores(nmds_result, "sites"))

# Get sample names and metadata
dist_samples <- rownames(as.matrix(bray_dist))
metadata_nmds <- sample_data(ps_filtered) %>% 
  as("data.frame") %>% 
  filter(rownames(.) %in% dist_samples)

# Merge NMDS scores with metadata
nmds_scores <- cbind(nmds_scores, metadata_nmds[dist_samples, ])

# Create NMDS plots
if ("take_all_seen" %in% colnames(nmds_scores) && "tillage_method" %in% colnames(nmds_scores)) {
  p_nmds <- plot_nmds(nmds_scores, "take_all_seen", "tillage_method", nmds_result$stress)
  ggsave(file.path(FIGURES, "nmds_takeall_tillage.png"), p_nmds, 
         width = 10, height = 8, dpi = 300)
}

# Basic NMDS plot
p_nmds_basic <- ggplot(nmds_scores, aes(x = NMDS1, y = NMDS2)) +
  geom_point(size = 3) +
  labs(
    title = "NMDS Plot (Bray-Curtis)",
    subtitle = paste("Stress =", round(nmds_result$stress, 3))
  ) +
  theme_minimal()

ggsave(file.path(FIGURES, "nmds_basic.png"), p_nmds_basic, 
       width = 8, height = 6, dpi = 300)

# PERMANOVA tests
cat("Running PERMANOVA tests...\n")
permanova_results <- list()

# Test main effects and interactions
if (all(c("take_all_seen", "tillage_method") %in% colnames(metadata_nmds))) {
  permanova_results$main <- adonis2(
    bray_dist ~ take_all_seen * tillage_method, 
    data = metadata_nmds
  )
}

if ("fertiliser_use" %in% colnames(metadata_nmds)) {
  permanova_results$fertilizer <- adonis2(
    bray_dist ~ fertiliser_use, 
    data = metadata_nmds
  )
}

# ANOSIM tests
anosim_results <- list()

if ("take_all_seen" %in% colnames(metadata_nmds)) {
  anosim_results$takeall <- anosim(bray_dist, grouping = metadata_nmds$take_all_seen)
  
  # Plot ANOSIM results
  png(file.path(FIGURES, "anosim_takeall.png"), width = 800, height = 600, res = 150)
  plot(anosim_results$takeall,
       xlab = "Take-All Disease Presence", 
       ylab = "Rank Dissimilarity", 
       main = "ANOSIM: Microbial Community vs. Take-All Disease")
  dev.off()
}

if ("tillage_method" %in% colnames(metadata_nmds)) {
  anosim_results$tillage <- anosim(bray_dist, grouping = metadata_nmds$tillage_method)
  
  # Plot ANOSIM results
  png(file.path(FIGURES, "anosim_tillage.png"), width = 800, height = 600, res = 150)
  plot(anosim_results$tillage,
       xlab = "Tillage Method",
       ylab = "Rank Dissimilarity",
       main = "ANOSIM: Microbial Community by Tillage Method")
  dev.off()
}

# Save results
saveRDS(list(
  permanova = permanova_results,
  anosim = anosim_results,
  nmds = nmds_result,
  nmds_scores = nmds_scores
), file.path(DATA_PROCESSED, "beta_diversity_results.rds"))

# Print summary
cat("PERMANOVA Results:\n")
print(permanova_results)

cat("\nANOSIM Results:\n")
print(anosim_results)

cat("Beta diversity analysis complete!\n")
