# =============================================================================
# analysis_functions.R
# Custom analysis functions for microbiome data
# =============================================================================

#' Calculate Hill diversity numbers from phyloseq object
#' @param ps_rarefied List of rarefied phyloseq objects
#' @param reference_order Sample order for alignment
calculate_hill_diversity <- function(ps_rarefied, reference_order) {
  all_shannon <- list()
  all_hill_q0 <- list()
  all_hill_q1 <- list()
  all_hill_q2 <- list()
  
  for (i in seq_along(ps_rarefied)) {
    phy <- ps_rarefied[[i]]
    
    # Extract OTU table (ensure samples are rows)
    otu <- as(otu_table(phy), "matrix")
    if (taxa_are_rows(phy)) { otu <- t(otu) }
    
    # Compute diversity metrics
    shannon <- vegan::diversity(otu, "shannon")
    hill_q0 <- hillR::hill_taxa(otu, q = 0)
    hill_q1 <- hillR::hill_taxa(otu, q = 1)
    hill_q2 <- hillR::hill_taxa(otu, q = 2)
    
    # Align with reference order
    all_shannon[[i]] <- shannon[reference_order]
    all_hill_q0[[i]] <- hill_q0[reference_order]
    all_hill_q1[[i]] <- hill_q1[reference_order]
    all_hill_q2[[i]] <- hill_q2[reference_order]
  }
  
  # Calculate averages
  results <- data.frame(
    Shannon = rowMeans(do.call(cbind, all_shannon), na.rm = TRUE),
    Hill_q0 = rowMeans(do.call(cbind, all_hill_q0), na.rm = TRUE),
    Hill_q1 = rowMeans(do.call(cbind, all_hill_q1), na.rm = TRUE),
    Hill_q2 = rowMeans(do.call(cbind, all_hill_q2), na.rm = TRUE)
  )
  
  return(results)
}

#' Prepare distance decay data
#' @param coords Matrix of coordinates (lat, long)
#' @param otu_dist Community distance matrix
#' @param metadata Metadata data frame
prepare_distance_decay <- function(coords, otu_dist, metadata) {
  # Compute geographic distance (in km)
  geo_dist <- distm(as.matrix(coords), fun = distHaversine) / 1000
  
  # Get sample IDs and pairwise combinations
  sample_ids <- rownames(coords)
  pairwise_idx <- combn(seq_along(sample_ids), 2)
  
  # Extract pairwise values
  decay_df <- data.frame(
    Sample1 = sample_ids[pairwise_idx[1,]],
    Sample2 = sample_ids[pairwise_idx[2,]],
    GeoDist = geo_dist[lower.tri(geo_dist)],
    BrayCurtis = as.vector(as.matrix(otu_dist)[lower.tri(otu_dist)])
  )
  
  # Add metadata for both samples
  metadata$SampleID <- rownames(metadata)
  
  decay_df <- decay_df %>%
    left_join(metadata, by = c("Sample1" = "SampleID")) %>%
    rename_with(~paste0(.x, "1"), -c(Sample1, Sample2, GeoDist, BrayCurtis)) %>%
    left_join(metadata, by = c("Sample2" = "SampleID")) %>%
    rename_with(~paste0(.x, "2"), -c(Sample1, Sample2, GeoDist, BrayCurtis, ends_with("1")))
  
  return(decay_df)
}

#' Run comprehensive statistical tests
#' @param metadata_df Metadata with diversity metrics
#' @param response_vars Vector of response variable names
#' @param predictor_vars Vector of predictor variable names
run_diversity_stats <- function(metadata_df, response_vars, predictor_vars) {
  results <- list()
  
  for (response in response_vars) {
    for (predictor in predictor_vars) {
      formula_str <- paste(response, "~", predictor)
      
      # ANOVA
      anova_result <- aov(as.formula(formula_str), data = metadata_df)
      
      # Store results
      results[[paste(response, predictor, sep = "_vs_")]] <- list(
        anova = summary(anova_result),
        tukey = TukeyHSD(anova_result)
      )
    }
  }
  
  return(results)
}

#' Extract significant results from statistical tests
#' @param stats_results Results from run_diversity_stats
#' @param alpha Significance threshold
extract_significant_results <- function(stats_results, alpha = 0.05) {
  sig_results <- list()
  
  for (test_name in names(stats_results)) {
    tukey_res <- stats_results[[test_name]]$tukey
    
    if (length(tukey_res) > 0) {
      for (factor_name in names(tukey_res)) {
        tukey_df <- as.data.frame(tukey_res[[factor_name]]) %>%
          rownames_to_column("comparison") %>%
          rename(p.adj = `p adj`) %>%
          mutate(
            test = test_name,
            factor = factor_name,
            comparison = gsub("-", " vs. ", comparison)
          ) %>%
          filter(p.adj < alpha)
        
        if (nrow(tukey_df) > 0) {
          sig_results[[paste(test_name, factor_name)]] <- tukey_df
        }
      }
    }
  }
  
  if (length(sig_results) == 0) {
    return(NULL)
  } else {
    return(bind_rows(sig_results))
  }
}
