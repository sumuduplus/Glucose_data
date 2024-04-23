library(tidyverse)
library(readxl)
library(phyloseq)
library(HTSSIP)


# Loading data ----

ps1 <- readRDS('Preprocess_data/ps1.rds')

colnames(sample_data(ps1))[5] <- 'Buoyant_density' 

tax_df <- as.data.frame(tax_table(ps1)) %>% 
  rownames_to_column(var = 'ASV_ID') %>% 
  unite(taxonomy, Phylum, Class, Order, Family, Genus, sep = '_', remove = FALSE, na.rm = TRUE)

# Filter samples

m0 <- data.frame(sample_data(ps1)) %>% 
  filter(Location == 'Deep',
         Treatment != 'Control',
         Days == 'day 0') %>% 
  drop_na(Buoyant_density)

m70 <- data.frame(sample_data(ps1)) %>% 
  filter(Location == 'Deep',
         Treatment != 'Control',
         Days == 'day 70') %>% 
  drop_na(Buoyant_density)

phy0 <- prune_samples(rownames(m0), ps1)
phy70 <- prune_samples(rownames(m70), ps1)

# Calculate enrichment

df_l2fc_0 <- HRSIP(phy0,
                   design = ~ Treatment,
                   padj_cutoff = 0.1,
                   density_windows = data.frame(density_min = 1.73,
                                                density_max =  1.85))

df_l2fc_0_sig <- df_l2fc_0 %>% 
  filter(padj < 0.1,
         log2FoldChange > 0)

df_l2fc_70 <- HRSIP(phy70,
                    design = ~ Treatment,
                    padj_cutoff = 0.1,
                    density_windows = data.frame(density_min = 1.73,
                                                 density_max = 1.85))

df_l2fc_70_sig <- df_l2fc_70 %>% 
  filter(padj < 0.1,
         log2FoldChange > 0)

# Plotting abundances by density T0

otu_table_0 <- as.data.frame(otu_table(phy0)) %>% 
  rownames_to_column(var = 'ASV_ID') %>% 
  pivot_longer(!ASV_ID, names_to = 'SampleID', values_to = 'abundance') %>% 
  group_by(SampleID) %>% 
  mutate(rel_abundance = abundance / sum(abundance)) %>% 
  inner_join(rownames_to_column(m0, var = 'SampleID'), by = 'SampleID') %>% 
  ungroup() %>% 
  select(ASV_ID, Treatment, Buoyant_density, rel_abundance) %>% 
  inner_join(tax_df, by = 'ASV_ID')

tax_names_0 <- set_names(otu_table_0$taxonomy,
                       nm = otu_table_0$ASV_ID)

otu_check_heavy_0 <- otu_table_0 %>% 
  filter(ASV_ID %in% df_l2fc_0_sig$OTU) %>%
  mutate(fraction = ifelse(Buoyant_density >= 1.73, 'heavy', 'light')) %>% 
  select(-Buoyant_density) %>% 
  filter(Treatment == 'Glucose_13C') %>% 
  group_by(ASV_ID, fraction) %>%
  summarise(mean_rel_abundance = mean(rel_abundance)) %>% 
  pivot_wider(names_from = 'fraction', values_from = 'mean_rel_abundance') %>% 
  filter(heavy > light)

line_plot_0 <- otu_table_0 %>% 
  filter(ASV_ID %in% otu_check_heavy_0$ASV_ID) %>%
  ggplot(aes(x = Buoyant_density, 
             y = rel_abundance,
             color = Treatment, 
             linetype = Treatment)) +
  geom_point() +
  geom_line() +
  geom_vline(xintercept = 1.73,
             linetype = 2) +
  labs(x = "Buoyant density (g/ml) ", 
       y = "Relative abundance",
       color = "Substrate", linetype="Substrate") +
  scale_color_manual(values=c("blue", "red")) +
  scale_linetype_manual(values=c(2, 1)) +
  facet_wrap(~ ASV_ID, scales = 'free_y',
             labeller = labeller(ASV_ID = tax_names_0)) +
  theme_bw()

line_plot_0

ggsave('16S/T0_DNA_SIP_taxa.svg', dpi = 300, height = 8, width = 16)  

otu_table_70 <- as.data.frame(otu_table(phy70)) %>% 
  rownames_to_column(var = 'ASV_ID') %>% 
  pivot_longer(!ASV_ID, names_to = 'SampleID', values_to = 'abundance') %>% 
  group_by(SampleID) %>% 
  mutate(rel_abundance = abundance / sum(abundance)) %>% 
  inner_join(rownames_to_column(m70, var = 'SampleID'), by = 'SampleID') %>% 
  ungroup() %>% 
  select(ASV_ID, Treatment, Buoyant_density, rel_abundance) %>% 
  inner_join(tax_df, by = 'ASV_ID')

tax_names_70 <- set_names(otu_table_70$taxonomy,
                       nm = otu_table_70$ASV_ID)

otu_check_heavy_70 <- otu_table_70 %>% 
  filter(ASV_ID %in% df_l2fc_70_sig$OTU) %>%
  mutate(fraction = ifelse(Buoyant_density >= 1.73, 'heavy', 'light')) %>% 
  select(-Buoyant_density) %>% 
  filter(Treatment == 'Glucose_13C') %>% 
  group_by(ASV_ID, fraction) %>%
  summarise(mean_rel_abundance = mean(rel_abundance)) %>% 
  pivot_wider(names_from = 'fraction', values_from = 'mean_rel_abundance') %>% 
  filter(heavy > light)

line_plot_70 <- otu_table_70 %>% 
  filter(ASV_ID %in% otu_check_heavy_70$ASV_ID) %>%
  ggplot(aes(x = Buoyant_density, 
             y = rel_abundance,
             color = Treatment, 
             linetype = Treatment)) +
  geom_point() +
  geom_line() +
  geom_vline(xintercept = 1.73,
             linetype = 2) +
  labs(x = "Buoyant density (g/ml) ", 
       y = "Relative abundance",
       color = "Substrate", linetype="Substrate") +
  scale_color_manual(values=c("blue", "red")) +
  scale_linetype_manual(values=c(2, 1)) +
  facet_wrap(~ ASV_ID, scales = 'free_y',
             labeller = labeller(ASV_ID = tax_names_70)) +
  theme_bw()

line_plot_70

ggsave('16S/T70_DNA_SIP_taxa.svg', dpi = 300, height = 8, width = 20)  

t0_heavy_taxa <- otu_table_0 %>% 
  filter(ASV_ID %in% otu_check_heavy_0$ASV_ID)%>%
  subset(select = -c(Treatment,Buoyant_density,rel_abundance, Species))%>%
  distinct(ASV_ID, .keep_all = TRUE)

t70_heavy_taxa <- otu_table_70 %>% 
  filter(ASV_ID %in% otu_check_heavy_70$ASV_ID)%>%
  subset(select = -c(Treatment,Buoyant_density,rel_abundance, Species))%>%
  distinct(ASV_ID, .keep_all = TRUE)

write_csv(t0_heavy_taxa, '16S/t0_heavy_taxa.csv')
write_csv(t70_heavy_taxa, '16S/t70_heavy_taxa.csv')

