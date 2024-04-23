# Libraries ----

library(tidyverse)
library(readxl)
library(ggplot2)
library(factoextra)
library(vegan)
library(dplyr)
library(ape)

# load data -----

lcms_metadata <- read_xlsx("Raw_data/Metadata_all.xlsx", sheet = 2) %>% 
  unite(uniqueID, Location, Treatment, Days,sep = "_", remove = FALSE) %>% 
  mutate(Days = factor(Days, levels = c('7', '14', '28', '42', '56', '70')))

lcms_data <-read_csv("Preprocess_data/LCMS_normalized.csv") %>% 
  left_join(select(lcms_metadata, SampleID, uniqueID), by = 'SampleID') %>% 
  select(-SampleID)

metadata_adonis <- lcms_metadata %>% 
  column_to_rownames(var = 'uniqueID')

metadata_adonis_surface <- lcms_metadata %>% 
  filter(!str_detect(uniqueID, 'Deep')) %>% 
  column_to_rownames(var = 'uniqueID')

metadata_adonis_deep <- lcms_metadata %>% 
  filter(!str_detect(uniqueID, 'Surface')) %>% 
  column_to_rownames(var = 'uniqueID')


# PCoA LCMS all samples ----

lcms_matrix_all <- lcms_data %>% 
  column_to_rownames(var = 'uniqueID')

lcms_dist_all <- vegdist(lcms_matrix_all, method = 'bray')

PCoA_LCMS_all <- pcoa(lcms_dist_all, correction="cailliez", rn=NULL)

pcoa_lcms_scores_all <- as.data.frame(PCoA_LCMS_all$vectors) %>% 
  rownames_to_column(var = 'uniqueID') %>% 
  left_join(lcms_metadata, by = 'uniqueID')

pcoa1_lc_all <- paste0('PCo1 ', round(PCoA_LCMS_all$values$Relative_eig[1] * 100, 2), '%')
pcoa2_lc_all <- paste0('PCo2 ', round(PCoA_LCMS_all$values$Relative_eig[2] * 100, 2), '%')

pcoa_lcms_plot_all <- pcoa_lcms_scores_all %>% 
  ggplot() +
  geom_point(aes(x = Axis.1,
                 y = Axis.2,
                 color = Location,
                 shape = Days),
             size = 5) +
  labs(title = 'LCMS all data (PCoA)',
       x = pcoa1_lc_all,
       y = pcoa2_lc_all) +
  scale_color_manual(values = c(Deep = '#977E59', Surface = '#4EB967')) +
  theme_bw() +
  theme(plot.title = element_text(face = 'bold', hjust = 0.5))

pcoa_lcms_plot_all


lcms_permanova_all <- adonis2(lcms_dist_all ~ Location * Treatment * Time, 
                              data = metadata_adonis,
                              method = 'bray',
                              permutations = 999) 
lcms_permanova_all

ggsave('lcms_permanova_all.svg', dpi = 300, height = 5, width = 6)  


# PCoA LCMS surface only ----

lcms_matrix_surface <- lcms_data %>% 
  filter(!str_detect(uniqueID, 'Deep')) %>% 
  column_to_rownames(var = 'uniqueID')

lcms_dist_surface <- vegdist(lcms_matrix_surface, method = 'bray')

PCoA_LCMS_surface <- pcoa(lcms_dist_surface, correction="cailliez", rn=NULL)

pcoa_lcms_scores_surface <- as.data.frame(PCoA_LCMS_surface$vectors) %>% 
  rownames_to_column(var = 'uniqueID') %>% 
  left_join(lcms_metadata, by = 'uniqueID')

pcoa1_lc_s <- paste0('PCo1 ', round(PCoA_LCMS_surface$values$Relative_eig[1] * 100, 2), '%')
pcoa2_lc_s <- paste0('PCo2 ', round(PCoA_LCMS_surface$values$Relative_eig[2] * 100, 2), '%')

pcoa_lcms_plot_surface <- pcoa_lcms_scores_surface %>% 
  ggplot() +
  geom_point(aes(x = Axis.1,
                 y = Axis.2,
                 color = Treatment,
                 shape = Days),
             size = 5) +
  labs(title = 'LCMS surface (PCoA)',
       x = pcoa1_lc_s,
       y = pcoa2_lc_s) +
  scale_color_manual(values = c(Control = '#EC9EC5', Glucose = '#14BCD8')) +
  theme_bw() +
  theme(plot.title = element_text(face = 'bold', hjust = 0.5))

pcoa_lcms_plot_surface


lcms_permanova_surface <- adonis2(lcms_dist_surface ~ Treatment * Time, 
                              data = metadata_adonis_surface,
                              method = 'bray',
                              permutations = 999) 
lcms_permanova_surface

ggsave('pcoa_lcms_plot_surface.svg', dpi = 300, height = 5, width = 6)  

# PCoA LCMS deep only ----

lcms_matrix_deep <- lcms_data %>% 
  filter(!str_detect(uniqueID, 'Surface')) %>% 
  column_to_rownames(var = 'uniqueID')

lcms_dist_deep <- vegdist(lcms_matrix_deep, method = 'bray')

PCoA_LCMS_deep <- pcoa(lcms_dist_deep, correction="cailliez", rn=NULL)

pcoa_lcms_scores_deep <- as.data.frame(PCoA_LCMS_deep$vectors) %>% 
  rownames_to_column(var = 'uniqueID') %>% 
  left_join(lcms_metadata, by = 'uniqueID')

pcoa1_lc_d <- paste0('PCo1 ', round(PCoA_LCMS_deep$values$Relative_eig[1] * 100, 2), '%')
pcoa2_lc_d <- paste0('PCo2 ', round(PCoA_LCMS_deep$values$Relative_eig[2] * 100, 2), '%')

pcoa_lcms_plot_deep <- pcoa_lcms_scores_deep %>% 
  ggplot() +
  geom_point(aes(x = Axis.1,
                 y = Axis.2,
                 color = Treatment,
                 shape = Days),
             size = 5) +
  labs(title = 'LCMS deep (PCoA)',
       x = pcoa1_lc_d,
       y = pcoa2_lc_d) +
  scale_color_manual(values = c(Control = '#EC9EC5', Glucose = '#14BCD8')) +
  theme_bw() +
  theme(plot.title = element_text(face = 'bold', hjust = 0.5))

pcoa_lcms_plot_deep


lcms_permanova_deep <- adonis2(lcms_dist_deep ~ Treatment * Time, 
                                  data = metadata_adonis_deep,
                                  method = 'bray',
                                  permutations = 999) 
lcms_permanova_deep

ggsave('pcoa_lcms_plot_deep.svg', dpi = 300, height = 5, width = 6)  

