# Libraries ----

library(tidyverse)
library(readxl)
library(ggplot2)
library(factoextra)
library(vegan)
library(dplyr)
library(ape)

# load data -----

nmr_metadata <- read_xlsx("Raw_data/Metadata_all.xlsx", sheet = 3) %>% 
  unite(uniqueID, Location, Treatment, Days,sep = "_", remove = FALSE) %>% 
  mutate(Days = factor(Days, levels = c('7', '14', '28', '42', '56', '70')))

nmr_data <- read_csv("Preprocess_data/NMR_normalized.csv") %>% 
  left_join(select(nmr_metadata, SampleID, uniqueID), by = 'SampleID') %>% 
  select(-SampleID)

metadata_adonis <- nmr_metadata %>% 
  column_to_rownames(var = 'uniqueID')

metadata_adonis_surface <- nmr_metadata %>% 
  filter(!str_detect(uniqueID, 'Deep')) %>% 
  column_to_rownames(var = 'uniqueID')

metadata_adonis_deep <- nmr_metadata %>% 
  filter(!str_detect(uniqueID, 'Surface')) %>% 
  column_to_rownames(var = 'uniqueID')

# PCoA LCMS all samples ----

nmr_matrix_all <- nmr_data %>% 
  column_to_rownames(var = 'uniqueID')

nmr_dist_all <- vegdist(nmr_matrix_all, method = 'bray')

PCoA_nmr_all <- pcoa(nmr_dist_all, correction="cailliez", rn=NULL)

pcoa_nmr_scores_all <- as.data.frame(PCoA_nmr_all$vectors) %>% 
  rownames_to_column(var = 'uniqueID') %>% 
  left_join(nmr_metadata, by = 'uniqueID')

pcoa1_nmr_all <- paste0('PCo1 ', round(PCoA_nmr_all$values$Rel_corr_eig[1] * 100, 2), '%')
pcoa2_nmr_all <- paste0('PCo2 ', round(PCoA_nmr_all$values$Rel_corr_eig[2]* 100, 2), '%')

pcoa_nmr_plot_all <- pcoa_nmr_scores_all %>% 
  ggplot() +
  geom_point(aes(x = Axis.1,
                 y = Axis.2,
                 color = Location,
                 shape = Days),
             size = 5) +
  labs(title = 'NMR all data (PCoA)',
       x = pcoa1_nmr_all,
       y = pcoa2_nmr_all) +
  scale_color_manual(values = c(Deep = '#977E59', Surface = '#4EB967')) +
  theme_bw() +
  theme(plot.title = element_text(face = 'bold', hjust = 0.5))

pcoa_nmr_plot_all


nmr_permanova_all <- adonis2(nmr_dist_all ~ Location * Treatment * Time, 
                              data = metadata_adonis,
                              method = 'bray',
                              permutations = 999) 
nmr_permanova_all

ggsave('nmr_permanova_all.svg', dpi = 300, height = 5, width = 6)  

# PCoA nmr surface only ----

nmr_matrix_surface <- nmr_data %>% 
  filter(!str_detect(uniqueID, 'Deep')) %>% 
  column_to_rownames(var = 'uniqueID')

nmr_dist_surface <- vegdist(nmr_matrix_surface, method = 'bray')

PCoA_nmr_surface <- pcoa(nmr_dist_surface, correction="cailliez", rn=NULL)

pcoa_nmr_scores_surface <- as.data.frame(PCoA_nmr_surface$vectors) %>% 
  rownames_to_column(var = 'uniqueID') %>% 
  left_join(nmr_metadata, by = 'uniqueID')

pcoa1_nmr_s <- paste0('PCo1 ', round(PCoA_nmr_surface$values$Rel_corr_eig[1] * 100, 2), '%')
pcoa2_nmr_s <- paste0('PCo2 ', round(PCoA_nmr_surface$values$Rel_corr_eig[2] * 100, 2), '%')

pcoa_nmr_plot_surface <- pcoa_nmr_scores_surface %>% 
  ggplot() +
  geom_point(aes(x = Axis.1,
                 y = Axis.2,
                 color = Treatment,
                 shape = Days),
             size = 5) +
  labs(title = 'NMR surface (PCoA)',
       x = pcoa1_nmr_s,
       y = pcoa2_nmr_s) +
  scale_color_manual(values = c(Control = '#EC9EC5', Glucose = '#14BCD8')) +
  theme_bw() +
  theme(plot.title = element_text(face = 'bold', hjust = 0.5))

pcoa_nmr_plot_surface


nmr_permanova_surface <- adonis2(nmr_dist_surface ~ Treatment * Time, 
                                  data = metadata_adonis_surface,
                                  method = 'bray',
                                  permutations = 999) 
nmr_permanova_surface

ggsave('pcoa_nmr_plot_surface.svg', dpi = 300, height = 5, width = 6)  

#Time split differently across early and late??

# PCoA nmr deep only ----

nmr_matrix_deep <- nmr_data %>% 
  filter(!str_detect(uniqueID, 'Surface')) %>% 
  column_to_rownames(var = 'uniqueID')

nmr_dist_deep <- vegdist(nmr_matrix_deep, method = 'bray')

PCoA_nmr_deep <- pcoa(nmr_dist_deep, correction="cailliez", rn=NULL)

pcoa_nmr_scores_deep <- as.data.frame(PCoA_nmr_deep$vectors) %>% 
  rownames_to_column(var = 'uniqueID') %>% 
  left_join(nmr_metadata, by = 'uniqueID')

pcoa1_nmr_d <- paste0('PCo1 ', round(PCoA_nmr_deep$values$Rel_corr_eig[1] * 100, 2), '%')
pcoa2_nmr_d <- paste0('PCo2 ', round(PCoA_nmr_deep$values$Rel_corr_eig[2] * 100, 2), '%')

pcoa_nmr_plot_deep <- pcoa_nmr_scores_deep %>% 
  ggplot() +
  geom_point(aes(x = Axis.1,
                 y = Axis.2,
                 color = Treatment,
                 shape = Days),
             size = 5) +
  labs(title = 'NMR deep (PCoA)',
       x = pcoa1_nmr_d,
       y = pcoa2_nmr_d) +
  scale_color_manual(values = c(Control = '#EC9EC5', Glucose = '#14BCD8')) +
  theme_bw() +
  theme(plot.title = element_text(face = 'bold', hjust = 0.5))

pcoa_nmr_plot_deep


nmr_permanova_deep <- adonis2(nmr_dist_deep ~ Treatment * Time, 
                               data = metadata_adonis_deep,
                               method = 'bray',
                               permutations = 999) 
nmr_permanova_deep

ggsave('pcoa_nmr_plot_deep.svg', dpi = 300, height = 5, width = 6)  
