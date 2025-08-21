#------------------------------------------------------------------------------
# Relative and Absolute Abundance Analysis
# ---------------------------------------

# Function to aggregate and calculate abundances at a specific taxonomic level
calculate_abundances <- function(physeq, level) {
  # Aggregate taxa at the specified level
  physeq_agg <- tax_glom(physeq, taxrank = level, NArm = FALSE)
  
  # Calculate relative abundances (percentages)
  physeq_rel <- transform_sample_counts(physeq_agg, function(x) x / sum(x) * 100)
  
  # Melt the data for plotting
  rel_melt <- psmelt(physeq_rel)
  abs_melt <- psmelt(physeq_agg)
  
  # Return both relative and absolute abundances
  list(relative = rel_melt, absolute = abs_melt)
}

# Calculate abundances at genus and phylum levels
genus_abundances <- calculate_abundances(ps_filtered, "genus")
phylum_abundances <- calculate_abundances(ps_filtered, "phylum")

# Plot relative abundance at phylum level
phylum_abundances$relative %>%
  group_by(field, phylum) %>%
  summarize(Abundance = sum(Abundance)) %>%
  ggplot(aes(x = field, y = Abundance, fill = phylum)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(title = "Absolute Abundance at Phylum Level",
       y = "Absolute Abundance",
       x = "Sample") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


# Plot relative abundance at phylum level
phylum_abundances$relative %>%
  group_by(field, phylum) %>%
  summarize(Abundance = sum(Abundance)) %>%
  ggplot(aes(x = field, y = Abundance, fill = phylum)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(title = "Absolute Abundance at Phylum Level",
       y = "Absolute Abundance (%)",
       x = "Field") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))



# Plot relative abundance at genus level (top 20 genera)
top_genera <- genus_abundances$relative %>%
  group_by(genus) %>%
  summarize(total = sum(Abundance)) %>%
  arrange(desc(total)) %>%
  head(20) %>%
  pull(genus)

genus_abundances$relative %>%
  filter(genus %in% top_genera) %>%
  group_by(Sample, genus) %>%
  summarize(Abundance = sum(Abundance)) %>%
  ggplot(aes(x = Sample, y = Abundance, fill = genus)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(title = "Relative Abundance of Top 20 Genera",
       y = "Relative Abundance (%)",
       x = "Sample") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


# Plot relative abundance at genus level (top 20 genera)
top_genera <- genus_abundances$relative %>%
  group_by(genus) %>%
  summarize(total = sum(Abundance)) %>%
  arrange(desc(total)) %>%
  head(20) %>%
  pull(genus)

genus_abundances$relative %>%
  filter(genus != "NA") %>%  # Explicitly exclude "NA" if present
  filter(genus %in% top_genera) %>%
  group_by(field, genus) %>%
  summarize(Abundance = sum(Abundance)) %>%
  ggplot(aes(x = field, y = Abundance, fill = genus)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(title = "Absolute Abundance of Top 20 Genera",
       y = "Absolute Abundance (%)",
       x = "Sample") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

### Heatmap calculations

# 1. Aggregate at genus level and merge by field
ps_field <- tax_glom(ps_filtered, "genus", NArm = FALSE)
ps_field <- merge_samples(ps_field, "field")
sample_data(ps_field)$field <- levels(sample_data(ps_filtered)$field)

# 2. Get relative abundance matrix
otu_mat <- as.data.frame(otu_table(ps_field))
if(!taxa_are_rows(ps_field)) otu_mat <- t(otu_mat)
otu_rel <- otu_mat / colSums(otu_mat)

# 3. Get genus names for labels
tax_info <- tax_table(ps_field)
genus_names <- tax_info[rownames(otu_rel), "genus"]
genus_names[is.na(genus_names)] <- "Unknown"

# 4. Keep top 50 genera
top_taxa <- names(sort(rowSums(otu_rel), decreasing = TRUE))[1:50]
otu_rel <- otu_rel[top_taxa, ]
genus_names <- genus_names[top_taxa]

# 5. Create heatmap
pheatmap(otu_rel,
         scale = "row",
         labels_row = genus_names,
         main = "Top Genera by Field",
         color = colorRampPalette(c("blue", "white", "red"))(100),
         fontsize_row = 9,
         fontsize_col = 10,
         angle_col = 90)
##Phyla 
# 1. Aggregate at phylum level and merge by field
ps_field <- tax_glom(ps_filtered, "phylum", NArm = FALSE)
ps_field <- merge_samples(ps_field, "field")
sample_data(ps_field)$field <- levels(sample_data(ps_filtered)$field)

# 2. Get relative abundance matrix
otu_mat <- as.data.frame(otu_table(ps_field))
if(!taxa_are_rows(ps_field)) otu_mat <- t(otu_mat)
otu_rel <- otu_mat / colSums(otu_mat)

# 3. Get phylum names for labels
tax_info <- tax_table(ps_field)
phylum_names <- tax_info[rownames(otu_rel), "phylum"]
phylum_names[is.na(phylum_names)] <- "Unknown"

# 4. Keep top 20 phyla
top_taxa <- names(sort(rowSums(otu_rel), decreasing = TRUE))[1:20]
otu_rel <- otu_rel[top_taxa, ]
phylum_names <- phylum_names[top_taxa]

# 5. Create heatmap
pheatmap(otu_rel,
         scale = "row",
         labels_row = phylum_names,
         main = "Top Phyla by Field",
         color = colorRampPalette(c("blue", "white", "red"))(100),
         fontsize_row = 9,
         fontsize_col = 10,
         angle_col = 90)

