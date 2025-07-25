# =============================================================================
# 04_alpha_diversity.R
# Alpha diversity analysis including Shannon and Hill numbers
# =============================================================================

source("scripts/03_quality_control.R")
source("functions/analysis_functions.R")
source("functions/plotting_functions.R")

cat("Starting alpha diversity analysis...\n")

# Rarefy the dataset to reduce randomness
ps_rarefied <- phyloseq_mult_raref(ps_filtered, SampSize = 1250, iter = 50)
reference_order <- rownames(sample_data(ps_filtered))

# Calculate Hill diversity metrics
cat("Calculating Hill diversity metrics...\n")
hill_metrics <- calculate_hill_diversity(ps_rarefied, reference_order)

# Update phyloseq metadata with diversity metrics
sample_data(ps_filtered)$Shannon <- hill_metrics$Shannon
sample_data(ps_filtered)$Hill_q0 <- hill_metrics$Hill_q0
sample_data(ps_filtered)$Hill_q1 <- hill_metrics$Hill_q1
sample_data(ps_filtered)$Hill_q2 <- hill_metrics$Hill_q2

# Extract metadata for plotting
metadata_div <- data.frame(sample_data(ps_filtered))
metadata_div <- metadata_div[!is.na(metadata_div$Shannon), ]

# Plot alpha diversity by different variables
variables_to_plot <- c("tillage_method", "fertiliser_use", "take_all_seen")
metrics_to_plot <- c("Shannon", "Hill_q0", "Hill_q1", "Hill_q2")

for (var in variables_to_plot) {
  if (var %in% colnames(metadata_div)) {
    for (metric in metrics_to_plot) {
      if (metric %in% colnames(metadata_div)) {
        p <- plot_alpha_diversity(
          metadata_div, 
          metric, 
          var, 
          paste(metric, "Diversity by", gsub("_", " ", var))
        )
        
        filename <- paste0("alpha_diversity_", metric, "_by_", var, ".png")
        ggsave(file.path(FIGURES, filename), p, width = 8, height = 6, dpi = 300)
      }
    }
  }
}

# Create faceted plot for take-all vs fertilizer use
if (all(c("take_all_seen", "fertiliser_use", "Shannon") %in% colnames(metadata_div))) {
  p_faceted <- ggplot(metadata_div, aes(x = take_all_seen, y = Shannon)) +
    geom_boxplot(aes(fill = take_all_seen), outlier.shape = NA, alpha = 0.7) +
    geom_jitter(aes(colour = take_all_seen), size = 1.5, width = 0.2) +
    stat_summary(fun = mean, geom = "point", shape = 20, size = 3, colour = "red") +
    labs(
      title = "Shannon Diversity by Take All observed vs Fertiliser Use",
      x = "Evidence of Take-All observed",
      y = "Shannon Diversity"
    ) +
    facet_wrap(~ fertiliser_use, scales = "free_y") +
    theme_minimal() +
    theme(legend.position = "none")
  
  ggsave(file.path(FIGURES, "shannon_takeall_fertilizer_faceted.png"), 
         p_faceted, width = 12, height = 8, dpi = 300)
}

# Save diversity data
fwrite(metadata_div, file.path(TABLES, "alpha_diversity_metrics.csv"))
saveRDS(ps_filtered, file.path(DATA_PROCESSED, "phyloseq_with_diversity.rds"))

cat("Alpha diversity analysis complete!\n")
