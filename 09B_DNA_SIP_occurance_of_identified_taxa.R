library(tidyverse)
library(readxl)
library(patchwork)
library(phyloseq)
library(HTSSIP)

# Loading data ----
ps1 <- readRDS('Preprocess_data/ps1.rds')

T0_SIP_taxa <- read_csv("16S/t0_heavy_taxa.csv")
T70_SIP_taxa <- read_csv("16S/t70_heavy_taxa.csv")

# Calculating Rel abundances ----

tax_df <- as.data.frame(tax_table(ps1)) %>% 
  rownames_to_column(var = 'ASV_ID') %>% 
  unite(taxonomy, Phylum, Class, Order, Family, Genus, sep = '_', remove = FALSE, na.rm = TRUE)

mDeep0 <- data.frame(sample_data(ps1)) %>% 
  filter(Location == 'Deep',
         Days == 'day 0')

phy_deep0 <- prune_samples(rownames(mDeep0), ps1)

otu_table_deep_0<- as.data.frame(otu_table(phy_deep0)) %>% 
  rownames_to_column(var = 'ASV_ID') %>% 
  pivot_longer(!ASV_ID, names_to = 'SampleID', values_to = 'abundance') %>% 
  group_by(SampleID) %>% 
  mutate(rel_abundance = abundance / sum(abundance)) %>% 
  inner_join(rownames_to_column(mDeep0, var = 'SampleID'), by = 'SampleID') %>% 
  ungroup() %>% 
  select(ASV_ID, Treatment, rel_abundance, SIPFraction) %>% 
  inner_join(tax_df, by = 'ASV_ID')



mDeep70 <- data.frame(sample_data(ps1)) %>% 
  filter(Location == 'Deep',
         Days == 'day 70')

phy_deep70 <- prune_samples(rownames(mDeep70), ps1)

otu_table_deep_70<- as.data.frame(otu_table(phy_deep70)) %>% 
  rownames_to_column(var = 'ASV_ID') %>% 
  pivot_longer(!ASV_ID, names_to = 'SampleID', values_to = 'abundance') %>% 
  group_by(SampleID) %>% 
  mutate(rel_abundance = abundance / sum(abundance)) %>% 
  inner_join(rownames_to_column(mDeep70, var = 'SampleID'), by = 'SampleID') %>% 
  ungroup() %>% 
  select(ASV_ID, Treatment, rel_abundance, SIPFraction) %>% 
  inner_join(tax_df, by = 'ASV_ID')



# filtering only DNA-SIP identified taxa ----

T0_filtered_taxa <- otu_table_deep_0 %>% 
  filter(otu_table_deep_0$ASV_ID %in% T0_SIP_taxa$ASV_ID) %>% 
  mutate(SIP= case_when(SIPFraction == '1' ~ 'F_12',
                        SIPFraction == '2' ~ 'F_11',
                        SIPFraction == '3' ~ 'F_10',
                        SIPFraction == '4' ~ 'F_09',
                        SIPFraction == '5' ~ 'F_08',
                        SIPFraction == '6' ~ 'F_07',
                        SIPFraction == '7' ~ 'F_06',
                        SIPFraction == '8' ~ 'F_05',
                        SIPFraction == '9' ~ 'F_04',
                        SIPFraction == '10' ~ 'F_03',
                        SIPFraction == '11' ~ 'F_02',
                        SIPFraction == '12' ~ 'F_01',
                        SIPFraction == 'S' ~ 'UN'))

T0_filtered_taxa_atT70 <- otu_table_deep_70 %>% 
  filter(otu_table_deep_70$ASV_ID %in% T0_SIP_taxa$ASV_ID)%>% 
  mutate(SIP= case_when(SIPFraction == '1' ~ 'F_12',
                        SIPFraction == '2' ~ 'F_11',
                        SIPFraction == '3' ~ 'F_10',
                        SIPFraction == '4' ~ 'F_09',
                        SIPFraction == '5' ~ 'F_08',
                        SIPFraction == '6' ~ 'F_07',
                        SIPFraction == '7' ~ 'F_06',
                        SIPFraction == '8' ~ 'F_05',
                        SIPFraction == '9' ~ 'F_04',
                        SIPFraction == '10' ~ 'F_03',
                        SIPFraction == '11' ~ 'F_02',
                        SIPFraction == '12' ~ 'F_01',
                        SIPFraction == 'S' ~ 'UN'))


T70_filtered_taxa <- otu_table_deep_70 %>% 
  filter(otu_table_deep_70$ASV_ID %in% T70_SIP_taxa$ASV_ID)%>% 
  mutate(SIP= case_when(SIPFraction == '1' ~ 'F_12',
                        SIPFraction == '2' ~ 'F_11',
                        SIPFraction == '3' ~ 'F_10',
                        SIPFraction == '4' ~ 'F_09',
                        SIPFraction == '5' ~ 'F_08',
                        SIPFraction == '6' ~ 'F_07',
                        SIPFraction == '7' ~ 'F_06',
                        SIPFraction == '8' ~ 'F_05',
                        SIPFraction == '9' ~ 'F_04',
                        SIPFraction == '10' ~ 'F_03',
                        SIPFraction == '11' ~ 'F_02',
                        SIPFraction == '12' ~ 'F_01',
                        SIPFraction == 'S' ~ 'UN'))

T70_filtered_taxa_atT0 <- otu_table_deep_0 %>% 
  filter(otu_table_deep_0$ASV_ID %in% T70_SIP_taxa$ASV_ID) %>% 
  mutate(SIP= case_when(SIPFraction == '1' ~ 'F_12',
                        SIPFraction == '2' ~ 'F_11',
                        SIPFraction == '3' ~ 'F_10',
                        SIPFraction == '4' ~ 'F_09',
                        SIPFraction == '5' ~ 'F_08',
                        SIPFraction == '6' ~ 'F_07',
                        SIPFraction == '7' ~ 'F_06',
                        SIPFraction == '8' ~ 'F_05',
                        SIPFraction == '9' ~ 'F_04',
                        SIPFraction == '10' ~ 'F_03',
                        SIPFraction == '11' ~ 'F_02',
                        SIPFraction == '12' ~ 'F_01',
                        SIPFraction == 'S' ~ 'UN'))


# Occurance in fractionated samples across ----

tax_names_0 <- set_names(otu_table_deep_0$taxonomy,
                         nm = otu_table_deep_0$ASV_ID)

tax_names_70 <- set_names(otu_table_deep_70$taxonomy,
                          nm = otu_table_deep_70$ASV_ID)


p1 <- ggplot(T0_filtered_taxa, aes(x = SIP, y = rel_abundance)) + 
  geom_line(aes(color = Treatment, linetype = Treatment, group = Treatment)) + 
  scale_color_manual(values = c("black", "blue","green"))+
  ggtitle("T0 heavy taxa at T0") +
  facet_wrap(~ASV_ID, scales = 'free_y',labeller = labeller(ASV_ID = tax_names_0), ncol = 1) +
  theme_bw() +theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))

p2 <- ggplot(T0_filtered_taxa_atT70, aes(x = SIP, y = rel_abundance)) + 
  geom_line(aes(color = Treatment, linetype = Treatment, group = Treatment)) + 
  ggtitle("T0 heavy taxa at T70") +
  scale_color_manual(values = c("black", "blue","green"))+
  facet_wrap(~ASV_ID, scales = 'free_y',labeller = labeller(ASV_ID = tax_names_70), ncol = 1) +
  theme_bw() +theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))

pT0 <- (p1 | p2) +
  plot_layout(guides = 'collect')

pT0

ggsave('T0_DNA_SIP_taxa.svg', dpi = 300, height = 7, width = 12)  


P3 <- ggplot(T70_filtered_taxa, aes(x = SIP, y = rel_abundance)) + 
  geom_line(aes(color = Treatment, linetype = Treatment, group = Treatment)) + 
  ggtitle("T70 heavy taxa at T70") +
  scale_color_manual(values = c("black", "blue","green"))+
  facet_wrap(~ASV_ID, scales = 'free_y',labeller = labeller(ASV_ID = tax_names_70),ncol = 1) +
  theme_bw() + theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))

P4 <- ggplot(T70_filtered_taxa_atT0, aes(x = SIP, y = rel_abundance)) + 
  geom_line(aes(color = Treatment, linetype = Treatment, group = Treatment)) + 
  ggtitle("T70 heavy taxa at T0") +
  scale_color_manual(values = c("black", "blue","green"))+
  facet_wrap(~ASV_ID, scales = 'free_y',labeller = labeller(ASV_ID = tax_names_70), ncol = 1) +
  theme_bw() +theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))

pT70 <- (P3 | P4) +
  plot_layout(guides = 'collect')

pT70

ggsave('T70_DNA_SIP_taxa.svg', dpi = 300, height = 16, width = 12)  


