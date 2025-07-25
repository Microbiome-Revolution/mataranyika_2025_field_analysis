# =============================================================================
# plotting_functions.R
# Custom plotting functions for microbiome analysis
# =============================================================================

#' Plot alpha diversity metrics
#' @param metadata_df Data frame with diversity metrics
#' @param metric_col Column name for diversity metric
#' @param group_col Column name for grouping variable
#' @param title Plot title
plot_alpha_diversity <- function(metadata_df, metric_col, group_col, title) {
  ggplot(metadata_df, aes_string(x = group_col, y = metric_col, fill = group_col)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    geom_jitter(aes_string(colour = group_col), size = 1.5, width = 0.2) +
    stat_summary(fun = mean, geom = "point", shape = 20, size = 3, colour = "red") +
    labs(title = title, x = group_col, y = metric_col) +
    theme_minimal() +
    theme(legend.position = "none")
}

#' Create relative abundance plot
#' @param ps_rel Phyloseq object with relative abundances
#' @param tax_level Taxonomic level to plot
#' @param top_n Number of top taxa to show
#' @param facet_var Variable to facet by
plot_relative_abundance <- function(ps_rel, tax_level, top_n = 25, facet_var = NULL) {
  # Agglomerate at specified taxonomic level
  ps_agg <- tax_glom(ps_rel, taxrank = tax_level)
  ps_melt <- psmelt(ps_agg)
  
  # Get top taxa
  top_taxa <- ps_melt %>%
    filter(!is.na(!!sym(tax_level)), !!sym(tax_level) != "Other") %>%
    group_by(!!sym(tax_level)) %>%
    summarise(Total = sum(Abundance)) %>%
    arrange(desc(Total)) %>%
    slice(1:top_n) %>%
    pull(!!sym(tax_level))
  
  # Prepare data for plotting
  ps_melt_filtered <- ps_melt %>%
    filter(!!sym(tax_level) %in% top_taxa) %>%
    mutate(!!tax_level := factor(!!sym(tax_level), levels = top_taxa))
  
  # Create color palette
  n_taxa <- length(top_taxa)
  taxa_colors <- scales::hue_pal()(n_taxa)
  names(taxa_colors) <- top_taxa
  
  p <- ggplot(ps_melt_filtered, aes(x = project, y = Abundance, fill = !!sym(tax_level))) +
    geom_bar(stat = "identity", position = position_fill(reverse = TRUE)) +
    scale_fill_manual(name = tax_level, values = taxa_colors, drop = FALSE) +
    labs(
      title = paste("Relative Abundance of Top", top_n, tax_level),
      x = "Site",
      y = "Relative Abundance (%)"
    ) +
    theme_minimal() +
    theme(legend.position = "right", legend.text = element_text(size = 6.5)) +
    guides(fill = guide_legend(ncol = 1)) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1))
  
  if (!is.null(facet_var)) {
    p <- p + facet_wrap(as.formula(paste("~", facet_var)), ncol = 3)
  }
  
  return(p)
}

#' Create NMDS plot
#' @param nmds_scores NMDS coordinates with metadata
#' @param color_var Variable for coloring points
#' @param shape_var Variable for point shapes
#' @param stress NMDS stress value
plot_nmds <- function(nmds_scores, color_var, shape_var = NULL, stress) {
  p <- ggplot(nmds_scores, aes_string(x = "NMDS1", y = "NMDS2", color = color_var))
  
  if (!is.null(shape_var)) {
    p <- p + aes_string(shape = shape_var)
  }
  
  p <- p +
    geom_point(size = 3) +
    stat_ellipse() +
    labs(
      title = "NMDS Plot (Bray-Curtis)",
      subtitle = paste("Stress =", round(stress, 3))
    ) +
    theme_minimal()
  
  return(p)
}

#' Create distance decay plot
#' @param decay_df Data frame with geographic and community distances
#' @param group_var Grouping variable for coloring
plot_distance_decay <- function(decay_df, group_var = NULL) {
  p <- ggplot(decay_df, aes(x = GeoDist, y = BrayCurtis))
  
  if (!is.null(group_var)) {
    p <- p + aes_string(color = group_var) +
      geom_smooth(method = "lm", se = FALSE)
  } else {
    p <- p + geom_smooth(method = "lm", color = "blue")
  }
  
  p <- p +
    geom_point(alpha = 0.4) +
    labs(
      title = "Distance Decay Relationship",
      x = "Geographic Distance (km)",
      y = "Bray-Curtis Dissimilarity"
    ) +
    theme_minimal()
  
  return(p)
}
