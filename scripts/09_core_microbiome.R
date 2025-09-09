# --------------------------------------------
# Core Microbiome Analysis
# --------------------------------------------


#########################################################

# --------------------------------------------
# Robust Core Microbiome Analysis
# --------------------------------------------

# Load needed packages
library(phyloseq)
library(microbiome)
library(ggplot2)
library(dplyr)
library(tidyr)
library(vegan)

# Use ps_filtered as the starting object

# ------------------------
# Safe Core Microbiome Analysis
# ------------------------

# Function to safely calculate core microbiome with error handling
safe_core <- function(ps, prevalence, detection = 0) {
  result <- tryCatch({
    core_obj <- core(ps, detection = detection, prevalence = prevalence)
    if (ntaxa(core_obj) == 0) {
      return(NULL)
    }
    return(core_obj)
  }, error = function(e) {
    message(paste("No core taxa found at", prevalence*100, "% prevalence"))
    return(NULL)
  })
  return(result)
}

# Define prevalence thresholds
prevalence_thresholds <- c(0.05, 0.20, 0.50, 0.80)
names(prevalence_thresholds) <- paste0("prev", prevalence_thresholds * 100)

# Calculate core microbiome for each threshold
# ------------------------
# Robust Core Microbiome Analysis
# ------------------------

# Function to safely calculate core microbiome with error handling
safe_core <- function(ps, prevalence, detection = 0) {
  result <- tryCatch({
    core_obj <- core(ps, detection = detection, prevalence = prevalence)
    if (ntaxa(core_obj) == 0) {
      return(NULL)
    }
    return(core_obj)
  }, error = function(e) {
    message(paste("No core taxa found at", prevalence*100, "% prevalence"))
    return(NULL)
  })
  return(result)
}

# Define prevalence thresholds
prevalence_thresholds <- c(0.05, 0.20, 0.50, 0.80)
core_list <- list()

# Calculate core microbiome for each threshold
for (prev in prevalence_thresholds) {
  threshold_name <- paste0("prev", prev * 100)
  core_obj <- safe_core(ps_filtered, prevalence = prev)
  if (!is.null(core_obj)) {
    core_list[[threshold_name]] <- core_obj
  }
}

# Print summary safely
cat("Core Microbiome Summary:\n")
cat("------------------------\n")
for (prev in prevalence_thresholds) {
  threshold_name <- paste0("prev", prev * 100)
  if (threshold_name %in% names(core_list)) {
    core_obj <- core_list[[threshold_name]]
    cat(sprintf("Core taxa at %.0f%% prevalence: %d taxa\n", 
                prev * 100, ntaxa(core_obj)))
  } else {
    cat(sprintf("No core taxa at %.0f%% prevalence\n", prev * 100))
  }
}

# Show what thresholds were successful
cat("\nSuccessful thresholds:", names(core_list), "\n")

# ------------------------
# Prevalence Analysis
# ------------------------

# Calculate prevalence for all taxa
prevalence_df <- microbiome::prevalence(ps_filtered) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Taxon") %>%
  rename(Prevalence = 2) %>%
  arrange(desc(Prevalence))

# Add taxonomy information
taxa_info <- as.data.frame(tax_table(ps_filtered)) %>%
  tibble::rownames_to_column("Taxon")

prevalence_df <- prevalence_df %>%
  left_join(taxa_info, by = "Taxon")

# Plot prevalence distribution
p1 <- ggplot(prevalence_df, aes(x = Prevalence)) +
  geom_histogram(binwidth = 0.05, fill = "steelblue", alpha = 0.7, 
                 color = "black", boundary = 0) +
  geom_vline(xintercept = prevalence_thresholds, linetype = "dashed", 
             color = "red", alpha = 0.7) +
  geom_text(data = data.frame(x = prevalence_thresholds, 
                              y = max(hist(prevalence_df$Prevalence, 
                                           breaks = seq(0, 1, 0.05), 
                                           plot = FALSE)$counts) * 0.9),
            aes(x = x, y = y, label = paste0(x*100, "%")), 
            vjust = -0.5, color = "red", size = 3) +
  labs(
    title = "Distribution of Taxa Prevalence",
    subtitle = "Red dashed lines show analysis thresholds",
    x = "Prevalence (Fraction of Samples)",
    y = "Number of Taxa"
  ) +
  theme_minimal() +
  theme(plot.subtitle = element_text(size = 10, color = "gray50"))

print(p1)

# ------------------------
# Heatmap Creation - Fixed Version
# ------------------------

# Check if we have data to work with
if (exists("ps_top_rel") && ntaxa(ps_top_rel) > 0 && nsamples(ps_top_rel) > 0) {
  
  cat("Creating heatmap for", ntaxa(ps_top_rel), "taxa across", nsamples(ps_top_rel), "samples\n")
  
  # Extract data step by step with validation
  abund_matrix <- as.matrix(otu_table(ps_top_rel))
  taxonomy_df <- as.data.frame(tax_table(ps_top_rel))
  sample_meta <- as.data.frame(sample_data(ps_top_rel))
  
  # Check orientation and transpose if needed
  if (nrow(abund_matrix) == nsamples(ps_top_rel) && ncol(abund_matrix) == ntaxa(ps_top_rel)) {
    cat("Transposing abundance matrix: samples as columns, taxa as rows\n")
    abund_matrix <- t(abund_matrix)
  }
  
  # Verify dimensions match
  if (nrow(abund_matrix) != nrow(taxonomy_df)) {
    cat("Warning: Abundance matrix and taxonomy dimensions don't match!\n")
    cat("Abundance taxa:", nrow(abund_matrix), "Taxonomy taxa:", nrow(taxonomy_df), "\n")
    
    # Use only common taxa
    common_taxa <- intersect(rownames(abund_matrix), rownames(taxonomy_df))
    if (length(common_taxa) > 0) {
      abund_matrix <- abund_matrix[common_taxa, , drop = FALSE]
      taxonomy_df <- taxonomy_df[common_taxa, , drop = FALSE]
      cat("Using", length(common_taxa), "common taxa\n")
    } else {
      stop("No common taxa between abundance and taxonomy tables!")
    }
  }
  
  # Create meaningful row names
  create_taxon_label <- function(tax_row) {
    taxon_name <- rownames(tax_row)
    if (!is.na(tax_row["genus"]) && tax_row["genus"] != "" && tax_row["genus"] != "NA") {
      return(paste0(tax_row["genus"], " (", taxon_name, ")"))
    } else if (!is.na(tax_row["family"]) && tax_row["family"] != "" && tax_row["family"] != "NA") {
      return(paste0("f__", tax_row["family"], " (", taxon_name, ")"))
    } else if (!is.na(tax_row["phylum"]) && tax_row["phylum"] != "" && tax_row["phylum"] != "NA") {
      return(paste0("p__", tax_row["phylum"], " (", taxon_name, ")"))
    } else {
      return(taxon_name)
    }
  }
  
  taxon_labels <- apply(taxonomy_df, 1, create_taxon_label)
  rownames(abund_matrix) <- taxon_labels
  
  # Handle sample ordering
  if ("field" %in% colnames(sample_meta) && ncol(abund_matrix) > 0) {
    tryCatch({
      if (ncol(abund_matrix) == nrow(sample_meta)) {
        sample_order <- order(sample_meta$field)
        abund_matrix <- abund_matrix[, sample_order, drop = FALSE]
        sample_meta <- sample_meta[sample_order, , drop = FALSE]
        cat("Samples ordered by field variable\n")
      }
    }, error = function(e) {
      cat("Could not order samples by field:", e$message, "\n")
    })
  }
  
  # Create heatmap with pheatmap
  tryCatch({
    # Prepare column annotations
    col_annot <- NULL
    if ("field" %in% colnames(sample_meta)) {
      col_annot <- data.frame(Field = as.factor(sample_meta$field))
      rownames(col_annot) <- colnames(abund_matrix)
    }
    
    # Add tillage method if available
    if ("tillage_method" %in% colnames(sample_meta)) {
      if (is.null(col_annot)) {
        col_annot <- data.frame(Tillage = as.factor(sample_meta$tillage_method))
        rownames(col_annot) <- colnames(abund_matrix)
      } else {
        col_annot$Tillage <- as.factor(sample_meta$tillage_method)
      }
    }
    
    # Apply log transformation
    transformed_matrix <- log1p(abund_matrix)
    
    # Create heatmap
    pheatmap::pheatmap(
      transformed_matrix,
      color = colorRampPalette(c("white", "lightblue", "blue", "darkblue"))(100),
      main = paste("Top", nrow(abund_matrix), "Most Prevalent Taxa"),
      cluster_rows = TRUE,
      cluster_cols = FALSE,
      show_colnames = ncol(abund_matrix) <= 30,
      show_rownames = TRUE,
      annotation_col = col_annot,
      fontsize_row = 8,
      fontsize_col = 7,
      angle_col = 45
    )
    
    # Save the plot
    png("top_prevalent_taxa_heatmap.png", width = 10, height = 8, units = "in", res = 300)
    pheatmap::pheatmap(
      transformed_matrix,
      color = colorRampPalette(c("white", "lightblue", "blue", "darkblue"))(100),
      main = paste("Top", nrow(abund_matrix), "Most Prevalent Taxa"),
      cluster_rows = TRUE,
      cluster_cols = FALSE,
      show_colnames = ncol(abund_matrix) <= 30,
      annotation_col = col_annot,
      fontsize_row = 8,
      fontsize_col = 7,
      angle_col = 45
    )
    dev.off()
    
    cat("Heatmap saved as 'top_prevalent_taxa_heatmap.png'\n")
    
  }, error = function(e) {
    cat("Error creating pheatmap:", e$message, "\n")
  })
  
} else {
  cat("Cannot create heatmap: ps_top_rel not found or has no data\n")
}

# ------------------------
# Abundance-Prevalence Plot
# ------------------------

# Calculate mean relative abundance
mean_abundance <- taxa_sums(ps_filtered) / nsamples(ps_filtered)
mean_abundance_df <- data.frame(
  Taxon = names(mean_abundance),
  MeanAbundance = as.numeric(mean_abundance)
)

# Combine with prevalence
abundance_prevalence_df <- prevalence_df %>%
  left_join(mean_abundance_df, by = "Taxon")

# Create abundance-prevalence plot
p3 <- ggplot(abundance_prevalence_df, aes(x = Prevalence, y = MeanAbundance)) +
  geom_point(alpha = 0.6, size = 2, color = "steelblue") +
  geom_smooth(method = "lm", color = "red", se = FALSE) +
  scale_y_log10() +
  labs(
    title = "Abundance vs Prevalence Relationship",
    x = "Prevalence (Fraction of Samples)",
    y = "Mean Abundance (log scale)"
  ) +
  theme_minimal()

print(p3)

# ------------------------
# Export Results
# ------------------------

# Save prevalence data
write.csv(prevalence_df, "taxa_prevalence_analysis.csv", row.names = FALSE)

# Save top prevalent taxa information
write.csv(top_prevalent, "top_prevalent_taxa.csv", row.names = FALSE)

# Save plots
ggsave("prevalence_distribution.png", p1, width = 10, height = 6)
ggsave("top_taxa_heatmap.png", p2, width = 12, height = 8)
ggsave("abundance_prevalence_plot.png", p3, width = 10, height = 6)

# ------------------------
# Session Info
# ------------------------

cat("\nAnalysis completed.\n")
cat("Most prevalent taxon:", top_prevalent$Taxon[1], 
    "(Prevalence:", round(top_prevalent$Prevalence[1], 3), ")\n")
cat("Number of taxa with prevalence > 0:", sum(prevalence_df$Prevalence > 0), "\n")
cat("Total number of taxa:", ntaxa(ps_filtered), "\n")


------------------------
  # Corrected Heatmap Creation
  # ------------------------

# Check if we have data to work with
if (exists("ps_top_rel") && ntaxa(ps_top_rel) > 0 && nsamples(ps_top_rel) > 0) {
  
  cat("Creating heatmap for", ntaxa(ps_top_rel), "taxa across", nsamples(ps_top_rel), "samples\n")
  
  # Extract data - careful with orientation
  abund_matrix <- as.matrix(otu_table(ps_top_rel))
  
  # Check orientation: phyloseq expects taxa as rows, samples as columns
  if (nrow(abund_matrix) == nsamples(ps_top_rel) && ncol(abund_matrix) == ntaxa(ps_top_rel)) {
    cat("Data is sample x taxa orientation - transposing to taxa x samples\n")
    abund_matrix <- t(abund_matrix)
  }
  
  taxonomy_df <- as.data.frame(tax_table(ps_top_rel))
  sample_meta <- as.data.frame(sample_data(ps_top_rel))
  
  # Verify we have the right orientation now
  cat("Final dimensions - Taxa:", nrow(abund_matrix), "Samples:", ncol(abund_matrix), "\n")
  cat("Taxonomy taxa:", nrow(taxonomy_df), "Sample metadata:", nrow(sample_meta), "\n")
  
  # Create meaningful row names (taxa labels)
  create_taxon_label <- function(tax_row) {
    if (!is.na(tax_row["genus"]) && tax_row["genus"] != "" && tax_row["genus"] != "NA") {
      return(paste0(tax_row["genus"], " (", rownames(tax_row), ")"))
    } else if (!is.na(tax_row["family"]) && tax_row["family"] != "" && tax_row["family"] != "NA") {
      return(paste0("f__", tax_row["family"], " (", rownames(tax_row), ")"))
    } else if (!is.na(tax_row["phylum"]) && tax_row["phylum"] != "" && tax_row["phylum"] != "NA") {
      return(paste0("p__", tax_row["phylum"], " (", rownames(tax_row), ")"))
    } else {
      return(rownames(tax_row))
    }
  }
  
  taxon_labels <- apply(taxonomy_df, 1, create_taxon_label)
  rownames(abund_matrix) <- taxon_labels
  
  # Handle sample ordering safely
  if ("field" %in% colnames(sample_meta) && ncol(abund_matrix) > 0) {
    tryCatch({
      sample_order <- order(sample_meta$field)
      if (length(sample_order) == ncol(abund_matrix)) {
        abund_matrix <- abund_matrix[, sample_order, drop = FALSE]
        sample_meta <- sample_meta[sample_order, , drop = FALSE]
        cat("Samples ordered by field variable\n")
      }
    }, error = function(e) {
      cat("Could not order samples by field:", e$message, "\n")
    })
  }
  
  # Create heatmap with pheatmap
  tryCatch({
    # Prepare column annotations if field exists
    col_annot <- NULL
    if ("field" %in% colnames(sample_meta)) {
      col_annot <- data.frame(Field = as.factor(sample_meta$field))
      rownames(col_annot) <- colnames(abund_matrix)
    }
    
    # Apply log transformation
    transformed_matrix <- log1p(abund_matrix)  # log(x + 1)
    
    # Create the heatmap
    heatmap_plot <- pheatmap::pheatmap(
      transformed_matrix,
      color = colorRampPalette(c("white", "lightblue", "blue", "darkblue"))(100),
      main = paste("Top", nrow(abund_matrix), "Most Prevalent Taxa"),
      cluster_rows = TRUE,
      cluster_cols = FALSE,
      show_colnames = ncol(abund_matrix) <= 50,  # Show sample names if not too many
      show_rownames = TRUE,
      annotation_col = col_annot,
      fontsize_row = 8,
      fontsize_col = 7,
      angle_col = 45,
      silent = TRUE  # Don't plot immediately
    )
    
    # Save the plot
    png("top_prevalent_taxa_heatmap.png", width = 12, height = 8, units = "in", res = 300)
    grid::grid.draw(heatmap_plot$gtable)
    dev.off()
    
    cat("Heatmap saved as 'top_prevalent_taxa_heatmap.png'\n")
    
    # Also display it
    grid::grid.draw(heatmap_plot$gtable)
    
  }, error = function(e) {
    cat("Error creating pheatmap:", e$message, "\n")
  })
  
} else {
  cat("Cannot create heatmap: ps_top_rel not found or has no data\n")
}

# ------------------------
# Simplified Alternative: ggplot2 Heatmap
# ------------------------

tryCatch({
  # Prepare data for ggplot heatmap
  heatmap_data <- as.data.frame(abund_matrix) %>%
    tibble::rownames_to_column("Taxon") %>%
    pivot_longer(cols = -Taxon, names_to = "Sample", values_to = "Abundance") %>%
    left_join(as.data.frame(sample_meta) %>% 
                tibble::rownames_to_column("Sample"), by = "Sample")
  
  # Check if the required columns exist
  if (!"field" %in% colnames(heatmap_data)) {
    stop("'field' column not found in sample metadata")
  }
  if (!"tillage_method" %in% colnames(heatmap_data)) {
    stop("'tillage_method' column not found in sample metadata")
  }
  
  # Create a combined field-tillage label for x-axis
  heatmap_data <- heatmap_data %>%
    mutate(Field_Tillage = paste(field, tillage_method, sep = " - "))
  
  # Calculate mean abundance for each taxon by field and tillage method
  summary_data <- heatmap_data %>%
    group_by(Taxon, field, tillage_method, Field_Tillage) %>%
    summarize(MeanAbundance = mean(Abundance, na.rm = TRUE),
              .groups = 'drop')
  
  # Order fields and tillage methods logically
  summary_data$field <- factor(summary_data$field)
  summary_data$tillage_method <- factor(summary_data$tillage_method)
  summary_data$Field_Tillage <- factor(summary_data$Field_Tillage,
                                       levels = unique(summary_data$Field_Tillage[order(summary_data$field, summary_data$tillage_method)]))
  
  # Create ggplot heatmap - Option 1: By Field with Tillage Colors
  p_heatmap1 <- ggplot(summary_data, aes(x = field, y = Taxon, fill = log1p(MeanAbundance))) +
    geom_tile() +
    scale_fill_gradientn(colors = c("white", "lightblue", "blue", "darkblue"),
                         name = "log(Mean Abundance + 1)") +
    labs(title = "Top Prevalent Taxa Heatmap - By Field",
         x = "Field", y = "Taxa") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
          axis.text.y = element_text(size = 8))
  
  print(p_heatmap1)
  ggsave("top_taxa_by_field.png", p_heatmap1, width = 12, height = 10)
  
  # Option 2: By Tillage Method with Field Colors
  p_heatmap2 <- ggplot(summary_data, aes(x = tillage_method, y = Taxon, fill = log1p(MeanAbundance))) +
    geom_tile() +
    scale_fill_gradientn(colors = c("white", "lightblue", "blue", "darkblue"),
                         name = "log(Mean Abundance + 1)") +
    labs(title = "Top Prevalent Taxa Heatmap - By Tillage Method",
         x = "Tillage Method", y = "Taxa") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
          axis.text.y = element_text(size = 8))
  
  print(p_heatmap2)
  ggsave("top_taxa_by_tillage.png", p_heatmap2, width = 10, height = 10)
  
  # Option 3: Faceted by Field and Tillage
  p_heatmap3 <- ggplot(summary_data, aes(x = tillage_method, y = Taxon, fill = log1p(MeanAbundance))) +
    geom_tile() +
    scale_fill_gradientn(colors = c("white", "lightblue", "blue", "darkblue"),
                         name = "log(Mean Abundance + 1)") +
    facet_grid(. ~ field, scales = "free_x", space = "free_x") +
    labs(title = "Top Prevalent Taxa Heatmap - By Field and Tillage Method",
         x = "Tillage Method", y = "Taxa") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
          axis.text.y = element_text(size = 8),
          strip.text = element_text(size = 9))
  
  print(p_heatmap3)
  ggsave("top_taxa_faceted_field_tillage.png", p_heatmap3, width = 16, height = 10)
  
  # Option 4: Combined Field-Tillage on x-axis
  p_heatmap4 <- ggplot(summary_data, aes(x = Field_Tillage, y = Taxon, fill = log1p(MeanAbundance))) +
    geom_tile() +
    scale_fill_gradientn(colors = c("white", "lightblue", "blue", "darkblue"),
                         name = "log(Mean Abundance + 1)") +
    labs(title = "Top Prevalent Taxa Heatmap - Field x Tillage Method",
         x = "Field - Tillage Method", y = "Taxa") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
          axis.text.y = element_text(size = 8))
  
  print(p_heatmap4)
  ggsave("top_taxa_combined_field_tillage.png", p_heatmap4, width = 14, height = 10)
  
}, error = function(e) {
  cat("Error creating ggplot heatmap:", e$message, "\n")
})

# ------------------------
# Additional Analysis: Tillage Method Effects
# ------------------------

tryCatch({
  # Statistical analysis of tillage effects
  if (exists("summary_data")) {
    # ANOVA or Kruskal-Wallis test for each taxon
    tillage_effects <- heatmap_data %>%
      group_by(Taxon) %>%
      summarize(
        p_value_kruskal = tryCatch({
          kruskal.test(Abundance ~ tillage_method)$p.value
        }, error = function(e) NA_real_),
        p_value_anova = tryCatch({
          summary(aov(Abundance ~ tillage_method))[[1]]$"Pr(>F)"[1]
        }, error = function(e) NA_real_),
        .groups = 'drop'
      ) %>%
      arrange(p_value_kruskal)
    
    cat("Top taxa with significant tillage method effects (Kruskal-Wallis):\n")
    print(head(tillage_effects, 10))
    
    write.csv(tillage_effects, "tillage_effects_analysis.csv", row.names = FALSE)
  }
  
}, error = function(e) {
  cat("Error in tillage effects analysis:", e$message, "\n")
})

# ------------------------
# Bar Plot by Tillage Method
# ------------------------

tryCatch({
  # Prepare data for bar plot by tillage method
  barplot_tillage <- heatmap_data %>%
    group_by(Taxon, tillage_method) %>%
    summarize(MeanAbundance = mean(Abundance, na.rm = TRUE),
              SE = sd(Abundance, na.rm = TRUE) / sqrt(n()),
              .groups = 'drop')
  
  p_bar_tillage <- ggplot(barplot_tillage, aes(x = tillage_method, y = MeanAbundance, fill = Taxon)) +
    geom_bar(stat = "identity", position = "fill") +
    labs(title = "Relative Abundance by Tillage Method",
         x = "Tillage Method", y = "Relative Abundance") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
          legend.position = "right")
  
  print(p_bar_tillage)
  ggsave("abundance_by_tillage.png", p_bar_tillage, width = 12, height = 8)
  
  # Individual taxon abundance by tillage
  p_point_tillage <- ggplot(barplot_tillage, aes(x = tillage_method, y = MeanAbundance, color = Taxon)) +
    geom_point(position = position_dodge(width = 0.5), size = 3) +
    geom_errorbar(aes(ymin = MeanAbundance - SE, ymax = MeanAbundance + SE),
                  position = position_dodge(width = 0.5), width = 0.2) +
    scale_y_log10() +
    labs(title = "Taxon Abundance by Tillage Method",
         x = "Tillage Method", y = "Mean Abundance (log scale)") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10))
  
  print(p_point_tillage)
  ggsave("taxon_abundance_by_tillage.png", p_point_tillage, width = 14, height = 8)
  
}, error = function(e) {
  cat("Error creating tillage bar plots:", e$message, "\n")
})

# ------------------------
# Bar Plot Alternative
# ------------------------

tryCatch({
  # Prepare data for bar plot
  barplot_data <- as.data.frame(abund_matrix) %>%
    tibble::rownames_to_column("Taxon") %>%
    pivot_longer(cols = -Taxon, names_to = "Sample", values_to = "Abundance") %>%
    group_by(Sample) %>%
    mutate(RelativeAbundance = Abundance / sum(Abundance)) %>%
    ungroup() %>%
    left_join(as.data.frame(sample_meta) %>% 
                tibble::rownames_to_column("Sample"), by = "Sample")
  
  # Create bar plot
  p_bar <- ggplot(barplot_data, aes(x = Sample, y = RelativeAbundance, fill = Taxon)) +
    geom_bar(stat = "identity", position = "stack") +
    labs(title = "Relative Abundance of Top Prevalent Taxa",
         y = "Relative Abundance", x = "Sample") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
          legend.position = "right")
  
  # Add field faceting if available
  if ("field" %in% colnames(barplot_data)) {
    p_bar <- p_bar + facet_grid(. ~ field, scales = "free_x", space = "free_x")
  }
  
  print(p_bar)
  ggsave("top_taxa_barplot.png", p_bar, width = 18, height = 8)
  
}, error = function(e) {
  cat("Error creating bar plot:", e$message, "\n")
})


tryCatch({
  # Prepare data for bar plot with field info
  barplot_data <- as.data.frame(abund_matrix) %>%
    tibble::rownames_to_column("Taxon") %>%
    pivot_longer(cols = -Taxon, names_to = "Sample", values_to = "Abundance") %>%
    group_by(Sample) %>%
    mutate(RelativeAbundance = Abundance / sum(Abundance)) %>%
    ungroup() %>%
    # join metadata (which must have Sample rownames and a column called 'field')
    left_join(
      as.data.frame(sample_meta) %>%
        tibble::rownames_to_column("Sample"),
      by = "Sample"
    ) %>%
    # now group by field (instead of just Sample) if you want summaries by field
    group_by(field, Taxon) %>%
    summarise(MeanRelativeAbundance = mean(RelativeAbundance), .groups = "drop")
}, error = function(e) {
  message("Error in preparing barplot data: ", e$message)
})

tryCatch({
  # Prepare data aggregated by field
  barplot_data <- as.data.frame(abund_matrix) %>%
    tibble::rownames_to_column("Taxon") %>%
    pivot_longer(cols = -Taxon, names_to = "Sample", values_to = "Abundance") %>%
    # relative abundance per sample
    group_by(Sample) %>%
    mutate(RelativeAbundance = Abundance / sum(Abundance)) %>%
    ungroup() %>%
    # join metadata (must include 'field' column)
    left_join(
      as.data.frame(sample_meta) %>%
        tibble::rownames_to_column("Sample"),
      by = "Sample"
    ) %>%
    # aggregate by field
    group_by(field, Taxon) %>%
    summarise(RelativeAbundance = mean(RelativeAbundance), .groups = "drop")
  
  # Create bar plot aggregated by field
  p_bar <- ggplot(barplot_data, aes(x = field, y = RelativeAbundance, fill = Taxon)) +
    geom_bar(stat = "identity", position = "stack") +
    labs(
      title = "Relative Abundance of Top Prevalent Taxa by Field",
      y = "Relative Abundance",
      x = "Field"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      legend.position = "right"
    )
  
  print(p_bar)
}, error = function(e) {
  message("Error in preparing or plotting barplot: ", e$message)
})

