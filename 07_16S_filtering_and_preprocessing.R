library(tidyverse)
library(readxl)
library(phyloseq)
library(vegan)

# Loading data ----

metadata <- read_xlsx('Raw_data/Metadata_all.xlsx', sheet = 4) %>% 
  column_to_rownames(var = 'SampleID')

asv_data <- read_tsv('Raw_data/feature_frequency_filtered_table_16S.tsv', skip = 1) %>% 
  column_to_rownames(var = '#OTU ID')

taxonomy <- read_tsv('Raw_data/taxonomy.tsv') %>% 
  separate(Taxon, 
           into = c('Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species'), 
           sep = '; ') %>% 
  select(-Confidence) %>% 
  mutate(across(everything(), ~str_remove(.x, '[a-z]__'))) %>% 
  filter(`Feature ID` %in% rownames(asv_data)) %>% 
  column_to_rownames(var = 'Feature ID') %>% 
  as.matrix()

taxonomy[taxonomy == ""] <- NA  


ps <- phyloseq(tax_table(taxonomy), sample_data(metadata), otu_table(asv_data, taxa_are_rows = TRUE))


# Number original taxa
paste0('Number of initial taxa: ', nrow(otu_table(ps)))


ps0 <- subset_taxa(ps, !is.na(Phylum))

# Number taxa after Phylum filtering

paste0('Number of taxa no phylum: ', nrow(otu_table(ps0)))

# Compute prevalence of each feature, store as data.frame
prevdf = apply(X = otu_table(ps0),
               MARGIN = ifelse(taxa_are_rows(ps0), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this data.frame
prevdf = data.frame(Prevalence = prevdf,
                    TotalAbundance = taxa_sums(ps0),
                    tax_table(ps0))

# Checking prevalence
plyr::ddply(prevdf, "Phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})


# Subset to the remaining phyla
prevdf1 = subset(prevdf, Phylum %in% get_taxa_unique(ps0, "Phylum"))
ggplot(prevdf1, aes(TotalAbundance, Prevalence / nsamples(ps0),color=Phylum)) +
  # Include a guess for parameter
  geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) + geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Phylum) + theme(legend.position="none")


#  Define prevalence threshold manually

prevalenceThreshold <- 2


# Execute prevalence filter, using `prune_taxa()` function
keepTaxa = rownames(prevdf1)[(prevdf1$Prevalence >= prevalenceThreshold)]
ps1 = prune_taxa(keepTaxa, ps0)

paste0('Number of taxa after prevalence filter: ', nrow(otu_table(ps1)))

# Agglomerating data based on genera (for stacked columns)

ps_agglomerated <- tax_glom(ps1, "Genus", NArm = TRUE)

tax_agglomerated_table <- as.data.frame(tax_table(ps_agglomerated)) %>% 
  rownames_to_column(var = 'FeatureID')

otu_agglomerated_table <- as.data.frame(otu_table(ps_agglomerated)) %>% 
  rownames_to_column(var = 'FeatureID') %>% 
  pivot_longer(!FeatureID, names_to = 'SampleID', values_to = 'abundance') %>% 
  inner_join(tax_agglomerated_table, by = 'FeatureID') 

otu_genus <- otu_agglomerated_table %>% 
  group_by(SampleID, Genus) %>% 
  summarise(abundance = sum(abundance)) %>% 
  filter(abundance > 0) %>% 
  group_by(SampleID) %>% 
  mutate(perc_abundance = abundance / sum(abundance)) %>% 
  ungroup() %>% 
  mutate(Genus_lump = fct_lump_n(Genus, w = perc_abundance, n = 9))

otu_genus_plot <- otu_genus %>% 
  filter(str_detect(SampleID, 'S$')) %>%
  inner_join(rownames_to_column(metadata, var = 'SampleID'), by = 'SampleID') %>% 
  filter(Treatment != 'Glucose_13C') %>% 
  ggplot() +
  geom_col(aes(x = SampleID,
               y = perc_abundance,
               fill = Genus_lump)) +
  theme_bw() +
  facet_grid(~ TimePoint + Location, scales = 'free') +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

otu_genus_plot

# Agglomerating data based on class (for stacked columns)

ps_agglomerated_class <- tax_glom(ps1, "Class", NArm = TRUE)

tax_agglomerated_table_class <- as.data.frame(tax_table(ps_agglomerated_class)) %>% 
  rownames_to_column(var = 'FeatureID')

otu_agglomerated_table_class <- as.data.frame(otu_table(ps_agglomerated_class)) %>% 
  rownames_to_column(var = 'FeatureID') %>% 
  pivot_longer(!FeatureID, names_to = 'SampleID', values_to = 'abundance') %>% 
  inner_join(tax_agglomerated_table_class, by = 'FeatureID') 

otu_class <- otu_agglomerated_table_class %>% 
  group_by(SampleID, Class) %>% 
  summarise(abundance = sum(abundance)) %>% 
  filter(abundance > 0) %>% 
  group_by(SampleID) %>% 
  mutate(perc_abundance = abundance / sum(abundance)) %>% 
  ungroup() %>% 
  mutate(Class_lump = fct_lump_n(Class, w = perc_abundance, n = 9))

otu_class_plot <- otu_class %>% 
  filter(str_detect(SampleID, 'S$')) %>%
  inner_join(rownames_to_column(metadata, var = 'SampleID'), by = 'SampleID') %>% 
  filter(Treatment != 'Glucose_13C') %>% 
  ggplot() +
  geom_col(aes(x = SampleID,
               y = perc_abundance,
               fill = Class_lump)) +
  theme_bw() +
  facet_grid(~ TimePoint + Location, scales = 'free') +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

otu_class_plot

ggsave('unfractionated_class_stacked.svg', dpi = 300, height = 6, width = 6)  

# Extracting data no agglomerated (w/wo Chord transformation)

tax_no_agglomerated_table <- as.data.frame(tax_table(ps1)) %>% 
  rownames_to_column(var = 'FeatureID')

otu_no_agglomerated_table <- as.data.frame(otu_table(ps1))

otu_chord <- otu_no_agglomerated_table

for(i in 1:ncol(otu_chord)){
  otu_chord[,i] <- otu_no_agglomerated_table[, i] / sqrt(sum(otu_no_agglomerated_table[,i] ^ 2))
}

otu_untransformed_tax <- otu_no_agglomerated_table %>% 
  rownames_to_column(var = 'FeatureID') %>% 
  inner_join(tax_no_agglomerated_table, by = 'FeatureID')

write_csv(otu_untransformed_tax, '16S/otu_untransformed_tax.csv')

otu_chord_tax <- otu_chord %>% 
  rownames_to_column(var = 'FeatureID') %>% 
  inner_join(tax_no_agglomerated_table, by = 'FeatureID')

write_csv(otu_chord_tax, '16S/otu_chord_transformed_tax.csv')

saveRDS(ps1, 'Preprocess_data/ps1.rds')

# diversity calculations

unfractionated_only_for_diversity <- subset(otu_untransformed_tax, select = -c(Kingdom,Phylum,Class,Order,Family,Genus,Species))%>%
  column_to_rownames(var = 'FeatureID')%>% 
  select(contains('S'))%>% 
  select(-contains('13C'))

Diversity_input <- merge(unfractionated_only_for_diversity, taxonomy, by=0, all=F)%>%
  remove_rownames %>% column_to_rownames(var="Row.names")

Shannon_12C_1000 <- diversity(Diversity_input$`12C-1000frS`, index = "shannon", equalize.groups = FALSE,
                              MARGIN = 1, base = exp(1))
Shannon_12C_1070 <- diversity(Diversity_input$`12C-10070frS`, index = "shannon", equalize.groups = FALSE,
                              MARGIN = 1, base = exp(1))
Shannon_12C_200 <- diversity(Diversity_input$`12C-200frS`, index = "shannon", equalize.groups = FALSE,
                              MARGIN = 1, base = exp(1))
Shannon_12C_2070 <- diversity(Diversity_input$`12C-2070frS`, index = "shannon", equalize.groups = FALSE,
                             MARGIN = 1, base = exp(1))

Shannon_Una_1000 <- diversity(Diversity_input$`Una-1000frS`, index = "shannon", equalize.groups = FALSE,
                              MARGIN = 1, base = exp(1))
Shannon_Una_1070 <- diversity(Diversity_input$`Una-10070frS`, index = "shannon", equalize.groups = FALSE,
                              MARGIN = 1, base = exp(1))
Shannon_Una_200 <- diversity(Diversity_input$`Una-200frS`, index = "shannon", equalize.groups = FALSE,
                             MARGIN = 1, base = exp(1))
Shannon_Una_2070 <- diversity(Diversity_input$`Una-2070frS`, index = "shannon", equalize.groups = FALSE,
                              MARGIN = 1, base = exp(1))

