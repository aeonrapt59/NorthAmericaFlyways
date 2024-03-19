library("data.table"); packageVersion("data.table")
library("phyloseq"); packageVersion("phyloseq") 
library("ggplot2"); packageVersion("ggplot2")       # graphics
library("tibble"); packageVersion("tibble")        # Needed for converting column to row names

# Load data
otu_mat_ebd <- fread("F:/ebd_relDec-2023/DV/final/ebd_USnCA.txt", header = TRUE, sep = "\t", fill = TRUE)
tax_mat_ebd <- fread("F:/ebd_relDec-2023/DV/final/taxanomy.txt", header = TRUE, sep = "\t", fill = TRUE)
meta_df_ebd <- fread("F:/ebd_relDec-2023/DV/final/meta_USnCA.txt", header = TRUE, sep = "\t", fill = TRUE)

# Convert first column to row names
otu_mat_ebd <- otu_mat_ebd %>% column_to_rownames("otu")
tax_mat_ebd <- tax_mat_ebd %>% column_to_rownames("otu")
meta_df_ebd <- meta_df_ebd %>% column_to_rownames("sample")

# Convert to matrix
otu_mat_ebd <- as.matrix(otu_mat_ebd)
tax_mat_ebd <- as.matrix(tax_mat_ebd)

# Create phyloseq objects
OTU_ebd <- otu_table(otu_mat_ebd, taxa_are_rows = TRUE)
TAX_ebd <- tax_table(tax_mat_ebd)
META_ebd <- sample_data(meta_df_ebd)

# Create phyloseq object
ps_ebd_USnCA <- phyloseq(OTU_ebd, TAX_ebd, META_ebd)
ps_ebd_USnCA

rm(otu_mat_ebd, meta_df_ebd, tax_mat_ebd, OTU_ebd, TAX_ebd, META_ebd)

# Perform garbage collection
gc()

# Step 1: Filter samples based on the 'Year' criteria (1900 or higher)
# Access the sample metadata
sample_metadata <- sample_data(ps_ebd_USnCA)

# Filter samples based on the 'Year' criteria (1910 or higher)
samples_to_keep <- sample_names(sample_metadata)[sample_metadata$Year >= 1910]

# Prune the phyloseq object to keep only the filtered samples
ps_ebd_USnCA_filtered_1910_or_higher <- prune_samples(samples_to_keep, ps_ebd_USnCA)

# Step 2: Filter taxa with non-zero abundance
# Calculate the total abundance of each taxon in the filtered dataset
taxa_abundance <- taxa_sums(ps_ebd_USnCA_filtered_1910_or_higher)

# Identify taxa with non-zero abundance
nonzero_taxa <- names(taxa_abundance[taxa_abundance > 0])

# Prune taxa with non-zero abundance from the filtered phyloseq object
ps_ebd_USnCA_1910 <- prune_taxa(nonzero_taxa, ps_ebd_USnCA_filtered_1910_or_higher)

# ps_ebd_USnCA_filtered_nonzero_taxa now contains only samples from 1910 or later and taxa with non-zero abundance

rm(samples_to_keep, ps_ebd_USnCA_filtered_1910_or_higher, taxa_abundance, nonzero_taxa)

# Perform garbage collection
gc()

library(microViz)

ps.Species_1910 <- tax_fix(ps_ebd_USnCA_1910)

ps.Species_1910 <- phyloseq_validate(ps.Species_1910, remove_undetected = TRUE)

p3.Species.aitchison_1910 <- ps.Species_1910 %>% 
  tax_transform("identity", rank = "Species") %>% 
  dist_calc(dist = "aitchison") %>% 
  ord_calc("PCoA") %>% 
  ord_plot(color = "State", size = 4, alpha = 0.5) +
  theme_bw() +
  ggside::geom_xsideboxplot(aes(fill = State), orientation = "y") +
  ggside::geom_ysideboxplot(aes(fill = State), orientation = "x") +
  ggside::scale_xsidey_discrete(labels = NULL) +
  ggside::scale_ysidex_discrete(labels = NULL) +
  ggside::theme_ggside_void() + 
  theme(panel.background = element_rect(fill = "white", color = 'black', linetype = 'solid', linewidth=1.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  geom_hline(yintercept = 0, linetype = 3, linewidth = 1)+ 
  geom_vline(xintercept = 0, linetype = 3, linewidth = 1)+
  theme(axis.title = element_text(size = 20, face = "bold"))+
  theme(axis.text.x = element_text(size = 17, color = "black"), 
        axis.text.y = element_text(size=17, color = "black")) +
  theme(aspect.ratio = 1)+
  theme(legend.title = element_text(size = 20, face = "bold")) +
  theme(legend.text = element_text(size = 17)) + 
  theme(plot.title = element_text(size = 24, face="bold"), 
        plot.caption = element_text(size = 15)) + 
  ggtitle("A) 1910-2023") + 
  geom_text(aes(label = State, colour = State, vjust=1.3), check_overlap = F)

p3.Species.aitchison_1910

rm(ps_ebd_USnCA_1910, ps.Species_1910)

# Perform garbage collection
gc()


# Step 1: Filter samples based on the 'Year' criteria (1920 or higher)

# Filter samples based on the 'Year' criteria (1920 or higher)
samples_to_keep <- sample_names(sample_metadata)[sample_metadata$Year >= 1920]

# Prune the phyloseq object to keep only the filtered samples
ps_ebd_USnCA_filtered_1920_or_higher <- prune_samples(samples_to_keep, ps_ebd_USnCA)

# Step 2: Filter taxa with non-zero abundance
# Calculate the total abundance of each taxon in the filtered dataset
taxa_abundance <- taxa_sums(ps_ebd_USnCA_filtered_1920_or_higher)

# Identify taxa with non-zero abundance
nonzero_taxa <- names(taxa_abundance[taxa_abundance > 0])

# Prune taxa with non-zero abundance from the filtered phyloseq object
ps_ebd_USnCA_1920 <- prune_taxa(nonzero_taxa, ps_ebd_USnCA_filtered_1920_or_higher)

# ps_ebd_USnCA_filtered_nonzero_taxa now contains only samples from 1920 or later and taxa with non-zero abundance

rm(samples_to_keep, ps_ebd_USnCA_filtered_1920_or_higher, taxa_abundance, nonzero_taxa)

# Perform garbage collection
gc()

ps.Species_1920 <- tax_fix(ps_ebd_USnCA_1920)

ps.Species_1920 <- phyloseq_validate(ps.Species_1920, remove_undetected = TRUE)

p3.Species.aitchison_1920 <- ps.Species_1920 %>% 
  tax_transform("identity", rank = "Species") %>% 
  dist_calc(dist = "aitchison") %>% 
  ord_calc("PCoA") %>% 
  ord_plot(color = "State", size = 4, alpha = 0.5) +
  theme_bw() +
  ggside::geom_xsideboxplot(aes(fill = State), orientation = "y") +
  ggside::geom_ysideboxplot(aes(fill = State), orientation = "x") +
  ggside::scale_xsidey_discrete(labels = NULL) +
  ggside::scale_ysidex_discrete(labels = NULL) +
  ggside::theme_ggside_void() + 
  theme(panel.background = element_rect(fill = "white", color = 'black', linetype = 'solid', linewidth=1.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  geom_hline(yintercept = 0, linetype = 3, linewidth = 1)+ 
  geom_vline(xintercept = 0, linetype = 3, linewidth = 1)+
  theme(axis.title = element_text(size = 20, face = "bold"))+
  theme(axis.text.x = element_text(size = 17, color = "black"), 
        axis.text.y = element_text(size=17, color = "black")) +
  theme(aspect.ratio = 1)+
  theme(legend.title = element_text(size = 20, face = "bold")) +
  theme(legend.text = element_text(size = 17)) + 
  theme(plot.title = element_text(size = 24, face="bold"), 
        plot.caption = element_text(size = 15)) + 
  ggtitle("B) 1920-2023") + 
  geom_text(aes(label = State, colour = State, vjust=1.3), check_overlap = F)

p3.Species.aitchison_1920

rm(ps_ebd_USnCA_1920, ps.Species_1920)

# Perform garbage collection
gc()


# Step 1: Filter samples based on the 'Year' criteria (1930 or higher)

# Filter samples based on the 'Year' criteria (1930 or higher)
samples_to_keep <- sample_names(sample_metadata)[sample_metadata$Year >= 1930]

# Prune the phyloseq object to keep only the filtered samples
ps_ebd_USnCA_filtered_1930_or_higher <- prune_samples(samples_to_keep, ps_ebd_USnCA)

# Step 2: Filter taxa with non-zero abundance
# Calculate the total abundance of each taxon in the filtered dataset
taxa_abundance <- taxa_sums(ps_ebd_USnCA_filtered_1930_or_higher)

# Identify taxa with non-zero abundance
nonzero_taxa <- names(taxa_abundance[taxa_abundance > 0])

# Prune taxa with non-zero abundance from the filtered phyloseq object
ps_ebd_USnCA_1930 <- prune_taxa(nonzero_taxa, ps_ebd_USnCA_filtered_1930_or_higher)

# ps_ebd_USnCA_filtered_nonzero_taxa now contains only samples from 1930 or later and taxa with non-zero abundance

rm(samples_to_keep, ps_ebd_USnCA_filtered_1930_or_higher, taxa_abundance, nonzero_taxa)

# Perform garbage collection
gc()

ps.Species_1930 <- tax_fix(ps_ebd_USnCA_1930)

ps.Species_1930 <- phyloseq_validate(ps.Species_1930, remove_undetected = TRUE)

p3.Species.aitchison_1930 <- ps.Species_1930 %>% 
  tax_transform("identity", rank = "Species") %>%
  dist_calc(dist = "aitchison") %>% 
  ord_calc("PCoA") %>% 
  ord_plot(color = "State", size = 4, alpha = 0.5) +
  theme_bw() +
  ggside::geom_xsideboxplot(aes(fill = State), orientation = "y") +
  ggside::geom_ysideboxplot(aes(fill = State), orientation = "x") +
  ggside::scale_xsidey_discrete(labels = NULL) +
  ggside::scale_ysidex_discrete(labels = NULL) +
  ggside::theme_ggside_void() + 
  theme(panel.background = element_rect(fill = "white", color = 'black', linetype = 'solid', linewidth=1.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  geom_hline(yintercept = 0, linetype = 3, linewidth = 1)+ 
  geom_vline(xintercept = 0, linetype = 3, linewidth = 1)+
  theme(axis.title = element_text(size = 20, face = "bold"))+
  theme(axis.text.x = element_text(size = 17, color = "black"), 
        axis.text.y = element_text(size=17, color = "black")) +
  theme(aspect.ratio = 1)+
  theme(legend.title = element_text(size = 20, face = "bold")) +
  theme(legend.text = element_text(size = 17)) + 
  theme(plot.title = element_text(size = 24, face="bold"), 
        plot.caption = element_text(size = 15)) + 
  ggtitle("C) 1930-2023") + 
  geom_text(aes(label = State, colour = State, vjust=1.3), check_overlap = F)

p3.Species.aitchison_1930

rm(ps_ebd_USnCA_1930, ps.Species_1930)

# Perform garbage collection
gc()


# Step 1: Filter samples based on the 'Year' criteria (1940 or higher)

# Filter samples based on the 'Year' criteria (1940 or higher)
samples_to_keep <- sample_names(sample_metadata)[sample_metadata$Year >= 1940]

# Prune the phyloseq object to keep only the filtered samples
ps_ebd_USnCA_filtered_1940_or_higher <- prune_samples(samples_to_keep, ps_ebd_USnCA)

# Step 2: Filter taxa with non-zero abundance
# Calculate the total abundance of each taxon in the filtered dataset
taxa_abundance <- taxa_sums(ps_ebd_USnCA_filtered_1940_or_higher)

# Identify taxa with non-zero abundance
nonzero_taxa <- names(taxa_abundance[taxa_abundance > 0])

# Prune taxa with non-zero abundance from the filtered phyloseq object
ps_ebd_USnCA_1940 <- prune_taxa(nonzero_taxa, ps_ebd_USnCA_filtered_1940_or_higher)

# ps_ebd_USnCA_filtered_nonzero_taxa now contains only samples from 1940 or later and taxa with non-zero abundance

rm(samples_to_keep, ps_ebd_USnCA_filtered_1940_or_higher, taxa_abundance, nonzero_taxa)

# Perform garbage collection
gc()

ps.Species_1940 <- tax_fix(ps_ebd_USnCA_1940)

ps.Species_1940 <- phyloseq_validate(ps.Species_1940, remove_undetected = TRUE)

p3.Species.aitchison_1940 <- ps.Species_1940 %>% 
  tax_transform("identity", rank = "Species") %>%
  dist_calc(dist = "aitchison") %>% 
  ord_calc("PCoA") %>% 
  ord_plot(color = "State", size = 4, alpha = 0.5) +
  theme_bw() +
  ggside::geom_xsideboxplot(aes(fill = State), orientation = "y") +
  ggside::geom_ysideboxplot(aes(fill = State), orientation = "x") +
  ggside::scale_xsidey_discrete(labels = NULL) +
  ggside::scale_ysidex_discrete(labels = NULL) +
  ggside::theme_ggside_void() + 
  theme(panel.background = element_rect(fill = "white", color = 'black', linetype = 'solid', linewidth=1.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  geom_hline(yintercept = 0, linetype = 3, linewidth = 1)+ 
  geom_vline(xintercept = 0, linetype = 3, linewidth = 1)+
  theme(axis.title = element_text(size = 20, face = "bold"))+
  theme(axis.text.x = element_text(size = 17, color = "black"), 
        axis.text.y = element_text(size=17, color = "black")) +
  theme(aspect.ratio = 1)+
  theme(legend.title = element_text(size = 20, face = "bold")) +
  theme(legend.text = element_text(size = 17)) + 
  theme(plot.title = element_text(size = 24, face="bold"), 
        plot.caption = element_text(size = 15)) + 
  ggtitle("D) 1940-2023") + 
  geom_text(aes(label = State, colour = State, vjust=1.3), check_overlap = F)

p3.Species.aitchison_1940

rm(ps_ebd_USnCA_1940, ps.Species_1940)

# Perform garbage collection
gc()


# Step 1: Filter samples based on the 'Year' criteria (1950 or higher)

# Filter samples based on the 'Year' criteria (1950 or higher)
samples_to_keep <- sample_names(sample_metadata)[sample_metadata$Year >= 1950]

# Prune the phyloseq object to keep only the filtered samples
ps_ebd_USnCA_filtered_1950_or_higher <- prune_samples(samples_to_keep, ps_ebd_USnCA)

# Step 2: Filter taxa with non-zero abundance
# Calculate the total abundance of each taxon in the filtered dataset
taxa_abundance <- taxa_sums(ps_ebd_USnCA_filtered_1950_or_higher)

# Identify taxa with non-zero abundance
nonzero_taxa <- names(taxa_abundance[taxa_abundance > 0])

# Prune taxa with non-zero abundance from the filtered phyloseq object
ps_ebd_USnCA_1950 <- prune_taxa(nonzero_taxa, ps_ebd_USnCA_filtered_1950_or_higher)

# ps_ebd_USnCA_filtered_nonzero_taxa now contains only samples from 1950 or later and taxa with non-zero abundance

rm(samples_to_keep, ps_ebd_USnCA_filtered_1950_or_higher, taxa_abundance, nonzero_taxa)

# Perform garbage collection
gc()

ps.Species_1950 <- tax_fix(ps_ebd_USnCA_1950)

ps.Species_1950 <- phyloseq_validate(ps.Species_1950, remove_undetected = TRUE)

p3.Species.aitchison_1950 <- ps.Species_1950 %>% 
  tax_transform("identity", rank = "Species") %>% 
  dist_calc(dist = "aitchison") %>% 
  ord_calc("PCoA") %>% 
  ord_plot(color = "State", size = 4, alpha = 0.5) +
  theme_bw() +
  ggside::geom_xsideboxplot(aes(fill = State), orientation = "y") +
  ggside::geom_ysideboxplot(aes(fill = State), orientation = "x") +
  ggside::scale_xsidey_discrete(labels = NULL) +
  ggside::scale_ysidex_discrete(labels = NULL) +
  ggside::theme_ggside_void() + 
  theme(panel.background = element_rect(fill = "white", color = 'black', linetype = 'solid', linewidth=1.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  geom_hline(yintercept = 0, linetype = 3, linewidth = 1)+ 
  geom_vline(xintercept = 0, linetype = 3, linewidth = 1)+
  theme(axis.title = element_text(size = 20, face = "bold"))+
  theme(axis.text.x = element_text(size = 17, color = "black"), 
        axis.text.y = element_text(size=17, color = "black")) +
  theme(aspect.ratio = 1)+
  theme(legend.title = element_text(size = 20, face = "bold")) +
  theme(legend.text = element_text(size = 17)) + 
  theme(plot.title = element_text(size = 24, face="bold"), 
        plot.caption = element_text(size = 15)) + 
  ggtitle("E) 1950-2023") + 
  geom_text(aes(label = State, colour = State, vjust=1.3), check_overlap = F)

p3.Species.aitchison_1950

rm(ps_ebd_USnCA_1950, ps.Species_1950)

# Perform garbage collection
gc()


# Step 1: Filter samples based on the 'Year' criteria (1960 or higher)

# Filter samples based on the 'Year' criteria (1960 or higher)
samples_to_keep <- sample_names(sample_metadata)[sample_metadata$Year >= 1960]

# Prune the phyloseq object to keep only the filtered samples
ps_ebd_USnCA_filtered_1960_or_higher <- prune_samples(samples_to_keep, ps_ebd_USnCA)

# Step 2: Filter taxa with non-zero abundance
# Calculate the total abundance of each taxon in the filtered dataset
taxa_abundance <- taxa_sums(ps_ebd_USnCA_filtered_1960_or_higher)

# Identify taxa with non-zero abundance
nonzero_taxa <- names(taxa_abundance[taxa_abundance > 0])

# Prune taxa with non-zero abundance from the filtered phyloseq object
ps_ebd_USnCA_1960 <- prune_taxa(nonzero_taxa, ps_ebd_USnCA_filtered_1960_or_higher)

# ps_ebd_USnCA_filtered_nonzero_taxa now contains only samples from 1960 or later and taxa with non-zero abundance

rm(samples_to_keep, ps_ebd_USnCA_filtered_1960_or_higher, taxa_abundance, nonzero_taxa)

# Perform garbage collection
gc()

ps.Species_1960 <- tax_fix(ps_ebd_USnCA_1960)

ps.Species_1960 <- phyloseq_validate(ps.Species_1960, remove_undetected = TRUE)

p3.Species.aitchison_1960 <- ps.Species_1960 %>% 
  tax_transform("identity", rank = "Species") %>%
  dist_calc(dist = "aitchison") %>% 
  ord_calc("PCoA") %>% 
  ord_plot(color = "State", size = 4, alpha = 0.5) +
  theme_bw() +
  ggside::geom_xsideboxplot(aes(fill = State), orientation = "y") +
  ggside::geom_ysideboxplot(aes(fill = State), orientation = "x") +
  ggside::scale_xsidey_discrete(labels = NULL) +
  ggside::scale_ysidex_discrete(labels = NULL) +
  ggside::theme_ggside_void() + 
  theme(panel.background = element_rect(fill = "white", color = 'black', linetype = 'solid', linewidth=1.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  geom_hline(yintercept = 0, linetype = 3, linewidth = 1)+ 
  geom_vline(xintercept = 0, linetype = 3, linewidth = 1)+
  theme(axis.title = element_text(size = 20, face = "bold"))+
  theme(axis.text.x = element_text(size = 17, color = "black"), 
        axis.text.y = element_text(size=17, color = "black")) +
  theme(aspect.ratio = 1)+
  theme(legend.title = element_text(size = 20, face = "bold")) +
  theme(legend.text = element_text(size = 17)) + 
  theme(plot.title = element_text(size = 24, face="bold"), 
        plot.caption = element_text(size = 15)) + 
  ggtitle("F) 1960-2023") + 
  geom_text(aes(label = State, colour = State, vjust=1.3), check_overlap = F)

p3.Species.aitchison_1960

rm(ps_ebd_USnCA_1960, ps.Species_1960)

# Perform garbage collection
gc()


# Step 1: Filter samples based on the 'Year' criteria (1970 or higher)

# Filter samples based on the 'Year' criteria (1970 or higher)
samples_to_keep <- sample_names(sample_metadata)[sample_metadata$Year >= 1970]

# Prune the phyloseq object to keep only the filtered samples
ps_ebd_USnCA_filtered_1970_or_higher <- prune_samples(samples_to_keep, ps_ebd_USnCA)

# Step 2: Filter taxa with non-zero abundance
# Calculate the total abundance of each taxon in the filtered dataset
taxa_abundance <- taxa_sums(ps_ebd_USnCA_filtered_1970_or_higher)

# Identify taxa with non-zero abundance
nonzero_taxa <- names(taxa_abundance[taxa_abundance > 0])

# Prune taxa with non-zero abundance from the filtered phyloseq object
ps_ebd_USnCA_1970 <- prune_taxa(nonzero_taxa, ps_ebd_USnCA_filtered_1970_or_higher)

# ps_ebd_USnCA_filtered_nonzero_taxa now contains only samples from 1970 or later and taxa with non-zero abundance

rm(samples_to_keep, ps_ebd_USnCA_filtered_1970_or_higher, taxa_abundance, nonzero_taxa)

# Perform garbage collection
gc()

ps.Species_1970 <- tax_fix(ps_ebd_USnCA_1970)

ps.Species_1970 <- phyloseq_validate(ps.Species_1970, remove_undetected = TRUE)

p3.Species.aitchison_1970 <- ps.Species_1970 %>% 
  tax_transform("identity", rank = "Species") %>%
  dist_calc(dist = "aitchison") %>% 
  ord_calc("PCoA") %>% 
  ord_plot(color = "State", size = 4, alpha = 0.5) +
  theme_bw() +
  ggside::geom_xsideboxplot(aes(fill = State), orientation = "y") +
  ggside::geom_ysideboxplot(aes(fill = State), orientation = "x") +
  ggside::scale_xsidey_discrete(labels = NULL) +
  ggside::scale_ysidex_discrete(labels = NULL) +
  ggside::theme_ggside_void() + 
  theme(panel.background = element_rect(fill = "white", color = 'black', linetype = 'solid', linewidth=1.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  geom_hline(yintercept = 0, linetype = 3, linewidth = 1)+ 
  geom_vline(xintercept = 0, linetype = 3, linewidth = 1)+
  theme(axis.title = element_text(size = 20, face = "bold"))+
  theme(axis.text.x = element_text(size = 17, color = "black"), 
        axis.text.y = element_text(size=17, color = "black")) +
  theme(aspect.ratio = 1)+
  theme(legend.title = element_text(size = 20, face = "bold")) +
  theme(legend.text = element_text(size = 17)) + 
  theme(plot.title = element_text(size = 24, face="bold"), 
        plot.caption = element_text(size = 15)) + 
  ggtitle("G) 1970-2023") + 
  geom_text(aes(label = State, colour = State, vjust=1.3), check_overlap = F)

p3.Species.aitchison_1970

rm(ps_ebd_USnCA_1970, ps.Species_1970)

# Perform garbage collection
gc()


# Step 1: Filter samples based on the 'Year' criteria (1980 or higher)

# Filter samples based on the 'Year' criteria (1980 or higher)
samples_to_keep <- sample_names(sample_metadata)[sample_metadata$Year >= 1980]

# Prune the phyloseq object to keep only the filtered samples
ps_ebd_USnCA_filtered_1980_or_higher <- prune_samples(samples_to_keep, ps_ebd_USnCA)

# Step 2: Filter taxa with non-zero abundance
# Calculate the total abundance of each taxon in the filtered dataset
taxa_abundance <- taxa_sums(ps_ebd_USnCA_filtered_1980_or_higher)

# Identify taxa with non-zero abundance
nonzero_taxa <- names(taxa_abundance[taxa_abundance > 0])

# Prune taxa with non-zero abundance from the filtered phyloseq object
ps_ebd_USnCA_1980 <- prune_taxa(nonzero_taxa, ps_ebd_USnCA_filtered_1980_or_higher)

# ps_ebd_USnCA_filtered_nonzero_taxa now contains only samples from 1980 or later and taxa with non-zero abundance

rm(samples_to_keep, ps_ebd_USnCA_filtered_1980_or_higher, taxa_abundance, nonzero_taxa)

# Perform garbage collection
gc()

ps.Species_1980 <- tax_fix(ps_ebd_USnCA_1980)

ps.Species_1980 <- phyloseq_validate(ps.Species_1980, remove_undetected = TRUE)

p3.Species.aitchison_1980 <- ps.Species_1980 %>% 
  tax_transform("identity", rank = "Species") %>%
  dist_calc(dist = "aitchison") %>% 
  ord_calc("PCoA") %>% 
  ord_plot(color = "State", size = 4, alpha = 0.5) +
  theme_bw() +
  ggside::geom_xsideboxplot(aes(fill = State), orientation = "y") +
  ggside::geom_ysideboxplot(aes(fill = State), orientation = "x") +
  ggside::scale_xsidey_discrete(labels = NULL) +
  ggside::scale_ysidex_discrete(labels = NULL) +
  ggside::theme_ggside_void() + 
  theme(panel.background = element_rect(fill = "white", color = 'black', linetype = 'solid', linewidth=1.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  geom_hline(yintercept = 0, linetype = 3, linewidth = 1)+ 
  geom_vline(xintercept = 0, linetype = 3, linewidth = 1)+
  theme(axis.title = element_text(size = 20, face = "bold"))+
  theme(axis.text.x = element_text(size = 17, color = "black"), 
        axis.text.y = element_text(size=17, color = "black")) +
  theme(aspect.ratio = 1)+
  theme(legend.title = element_text(size = 20, face = "bold")) +
  theme(legend.text = element_text(size = 17)) + 
  theme(plot.title = element_text(size = 24, face="bold"), 
        plot.caption = element_text(size = 15)) + 
  ggtitle("H) 1980-2023") + 
  geom_text(aes(label = State, colour = State, vjust=1.3), check_overlap = F)

p3.Species.aitchison_1980

rm(ps_ebd_USnCA_1980, ps.Species_1980)

# Perform garbage collection
gc()


# Step 1: Filter samples based on the 'Year' criteria (1990 or higher)

# Filter samples based on the 'Year' criteria (1990 or higher)
samples_to_keep <- sample_names(sample_metadata)[sample_metadata$Year >= 1990]

# Prune the phyloseq object to keep only the filtered samples
ps_ebd_USnCA_filtered_1990_or_higher <- prune_samples(samples_to_keep, ps_ebd_USnCA)

# Step 2: Filter taxa with non-zero abundance
# Calculate the total abundance of each taxon in the filtered dataset
taxa_abundance <- taxa_sums(ps_ebd_USnCA_filtered_1990_or_higher)

# Identify taxa with non-zero abundance
nonzero_taxa <- names(taxa_abundance[taxa_abundance > 0])

# Prune taxa with non-zero abundance from the filtered phyloseq object
ps_ebd_USnCA_1990 <- prune_taxa(nonzero_taxa, ps_ebd_USnCA_filtered_1990_or_higher)

# ps_ebd_USnCA_filtered_nonzero_taxa now contains only samples from 1990 or later and taxa with non-zero abundance

rm(samples_to_keep, ps_ebd_USnCA_filtered_1990_or_higher, taxa_abundance, nonzero_taxa)

# Perform garbage collection
gc()

ps.Species_1990 <- tax_fix(ps_ebd_USnCA_1990)

ps.Species_1990 <- phyloseq_validate(ps.Species_1990, remove_undetected = TRUE)

p3.Species.aitchison_1990 <- ps.Species_1990 %>% 
  tax_transform("identity", rank = "Species") %>%
  dist_calc(dist = "aitchison") %>% 
  ord_calc("PCoA") %>% 
  ord_plot(color = "State", size = 4, alpha = 0.5) +
  theme_bw() +
  ggside::geom_xsideboxplot(aes(fill = State), orientation = "y") +
  ggside::geom_ysideboxplot(aes(fill = State), orientation = "x") +
  ggside::scale_xsidey_discrete(labels = NULL) +
  ggside::scale_ysidex_discrete(labels = NULL) +
  ggside::theme_ggside_void() + 
  theme(panel.background = element_rect(fill = "white", color = 'black', linetype = 'solid', linewidth=1.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  geom_hline(yintercept = 0, linetype = 3, linewidth = 1)+ 
  geom_vline(xintercept = 0, linetype = 3, linewidth = 1)+
  theme(axis.title = element_text(size = 20, face = "bold"))+
  theme(axis.text.x = element_text(size = 17, color = "black"), 
        axis.text.y = element_text(size=17, color = "black")) +
  theme(aspect.ratio = 1)+
  theme(legend.title = element_text(size = 20, face = "bold")) +
  theme(legend.text = element_text(size = 17)) + 
  theme(plot.title = element_text(size = 24, face="bold"), 
        plot.caption = element_text(size = 15)) + 
  ggtitle("I) 1990-2023") + 
  geom_text(aes(label = State, colour = State, vjust=1.3), check_overlap = F)

p3.Species.aitchison_1990

rm(ps_ebd_USnCA_1990, ps.Species_1990)

# Perform garbage collection
gc()


# Load data
otu_mat_ebd <- fread("F:/ebd_relDec-2023/DV/final/ebd_USnCA.txt", header = TRUE, sep = "\t", fill = TRUE)
tax_mat_ebd <- fread("F:/ebd_relDec-2023/DV/final/taxanomy.txt", header = TRUE, sep = "\t", fill = TRUE)
meta_df_ebd <- fread("F:/ebd_relDec-2023/DV/final/meta_USnCA_flyways.txt", header = TRUE, sep = "\t", fill = TRUE)

# Convert first column to row names
otu_mat_ebd <- otu_mat_ebd %>% column_to_rownames("otu")
tax_mat_ebd <- tax_mat_ebd %>% column_to_rownames("otu")
meta_df_ebd <- meta_df_ebd %>% column_to_rownames("sample")

# Convert to matrix
otu_mat_ebd <- as.matrix(otu_mat_ebd)
tax_mat_ebd <- as.matrix(tax_mat_ebd)

# Create phyloseq objects
OTU_ebd <- otu_table(otu_mat_ebd, taxa_are_rows = TRUE)
TAX_ebd <- tax_table(tax_mat_ebd)
META_ebd <- sample_data(meta_df_ebd)

# Create phyloseq object
ps_ebd_USnCA_fly <- phyloseq(OTU_ebd, TAX_ebd, META_ebd)
ps_ebd_USnCA_fly

rm(otu_mat_ebd, meta_df_ebd, tax_mat_ebd, OTU_ebd, TAX_ebd, META_ebd)

# Perform garbage collection
gc()

# Step 1: Filter samples based on the 'Year' criteria (1990 or higher)

# Filter samples based on the 'Year' criteria (1990 or higher)
samples_to_keep <- sample_names(sample_metadata)[sample_metadata$Year >= 1990]

# Prune the phyloseq object to keep only the filtered samples
ps_ebd_USnCA_filtered_1990_or_higher <- prune_samples(samples_to_keep, ps_ebd_USnCA_fly)

# Step 2: Filter taxa with non-zero abundance
# Calculate the total abundance of each taxon in the filtered dataset
taxa_abundance <- taxa_sums(ps_ebd_USnCA_filtered_1990_or_higher)

# Identify taxa with non-zero abundance
nonzero_taxa <- names(taxa_abundance[taxa_abundance > 0])

# Prune taxa with non-zero abundance from the filtered phyloseq object
ps_ebd_USnCA_1990_fly <- prune_taxa(nonzero_taxa, ps_ebd_USnCA_filtered_1990_or_higher)

# ps_ebd_USnCA_filtered_nonzero_taxa now contains only samples from 1990 or later and taxa with non-zero abundance

rm(samples_to_keep, ps_ebd_USnCA_filtered_1990_or_higher, taxa_abundance, nonzero_taxa)

# Perform garbage collection
gc()

ps.Species_1990_fly <- tax_fix(ps_ebd_USnCA_1990_fly)

ps.Species_1990_fly <- phyloseq_validate(ps.Species_1990_fly, remove_undetected = TRUE)

p3.Species.aitchison_1990_fly <- ps.Species_1990_fly %>% 
  tax_transform("identity", rank = "Species") %>%
  dist_calc(dist = "aitchison") %>% 
  ord_calc("PCoA") %>% 
  ord_plot(color = "Flyway", size = 4, alpha = 0.5) +
  theme_bw() +
  ggside::geom_xsideboxplot(aes(fill = Flyway), orientation = "y") +
  ggside::geom_ysideboxplot(aes(fill = Flyway), orientation = "x") +
  ggside::scale_xsidey_discrete(labels = NULL) +
  ggside::scale_ysidex_discrete(labels = NULL) +
  ggside::theme_ggside_void() + 
  theme(panel.background = element_rect(fill = "white", color = 'black', linetype = 'solid', linewidth=1.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  geom_hline(yintercept = 0, linetype = 3, linewidth = 1)+ 
  geom_vline(xintercept = 0, linetype = 3, linewidth = 1)+
  theme(axis.title = element_text(size = 20, face = "bold"))+
  theme(axis.text.x = element_text(size = 17, color = "black"), 
        axis.text.y = element_text(size=17, color = "black")) +
  theme(aspect.ratio = 1)+
  theme(legend.title = element_text(size = 20, face = "bold")) +
  theme(legend.text = element_text(size = 17)) + 
  theme(plot.title = element_text(size = 24, face="bold"), 
        plot.caption = element_text(size = 15)) + 
  ggtitle("Bird flyways from beta-diversity")

p3.Species.aitchison_1990_fly


ps.Species_1990_fly %>%
  ps_mutate(
  ) %>%
  tax_transform("clr", rank = "Genus") %>%
  ord_calc(
    constraints = c(),
    # method = "RDA", # Note: you can specify RDA explicitly, and it is good practice to do so, but microViz can guess automatically that you want an RDA here (helpful if you don't remember the name?)
    scale_cc = FALSE # doesn't make a difference
  ) %>%
  ord_plot(
    colour = "Flyway", size = 5, alpha = 0.8,
    plot_taxa = 1:20, tax_lab_style = tax_lab_style(size = 6, alpha = 0.5),
    constraint_lab_style = constraint_lab_style(
      alpha = 0.8, size = 6, perpendicular = TRUE
    )
  ) + 
  theme(panel.background = element_rect(fill = "white", color = 'black', linetype = 'solid', size=1.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  geom_hline(yintercept = 0, linetype = 3, linewidth = 1)+ 
  geom_vline(xintercept = 0, linetype = 3, linewidth = 1)+
  theme(axis.title = element_text(size = 20, face = "bold"))+
  theme(axis.text.x = element_text(size = 17, color = "black"), 
        axis.text.y = element_text(size=17, color = "black")) +
  theme(aspect.ratio = 1)+
  theme(legend.title = element_text(size = 20, face = "bold")) +
  theme(legend.text = element_text(size = 17)) + 
  theme(plot.title = element_text(size = 24, face="bold"), 
        plot.caption = element_text(size = 15)) + 
  ggtitle("Redundancy Analysis") + 
  #geom_text(aes(label = Samples, colour = site), check_overlap = TRUE) + 
  stat_ellipse(geom = "polygon", aes(color = Flyway, fill = Flyway), alpha = 0.25, size = 0.7, linetype = 1)

