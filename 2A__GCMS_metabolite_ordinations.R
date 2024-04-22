# Libraries ----

library(tidyverse)
library(readxl)
library(ggplot2)
library(factoextra)
library(vegan)
library(dplyr)
library(ape)

# load data -----

gcms_metadata <- read_xlsx("Raw_data/Metadata_all.xlsx", sheet = 1) %>% 
  unite(uniqueID, Location, Treatment, Days,sep = "_", remove = FALSE) %>% 
  mutate(Days = factor(Days, levels = c('7', '14', '28', '42', '56', '70')))

gcms_data <- read_csv("Preprocess_data/GCMS_normalized.csv") %>% 
  left_join(select(gcms_metadata, SampleID, uniqueID), by = 'SampleID') %>% 
  select(-SampleID)

metadata_adonis <- gcms_metadata %>% 
  filter(!str_detect(Location, 'Mid')) %>% 
  column_to_rownames(var = 'uniqueID')

metadata_adonis_surface <- gcms_metadata %>% 
  filter(!str_detect(Location, 'Mid')) %>% 
  filter(!str_detect(uniqueID, 'Deep')) %>% 
  column_to_rownames(var = 'uniqueID')

metadata_adonis_deep <- gcms_metadata %>% 
  filter(!str_detect(Location, 'Mid')) %>% 
  filter(!str_detect(uniqueID, 'Surface')) %>% 
  column_to_rownames(var = 'uniqueID')


# PCoA GCMS all samples ----

gcms_matrix_all <- gcms_data %>% 
  column_to_rownames(var = 'uniqueID')

gcms_dist_all <- vegdist(gcms_matrix_all, method = 'bray')

PCoA_gcMS_all <- pcoa(gcms_dist_all, correction="cailliez", rn=NULL)

pcoa_gcms_scores_all <- as.data.frame(PCoA_gcMS_all$vectors) %>% 
  rownames_to_column(var = 'uniqueID') %>% 
  left_join(gcms_metadata, by = 'uniqueID')

pcoa1_gc_all <- paste0('PCo1 ', round(PCoA_gcMS_all$values$Rel_corr_eig[1] * 100, 2), '%')
pcoa2_gc_all <- paste0('PCo2 ', round(PCoA_gcMS_all$values$Rel_corr_eig[2] * 100, 2), '%')

pcoa_gcms_plot_all <- pcoa_gcms_scores_all %>% 
  ggplot() +
  geom_point(aes(x = Axis.1,
                 y = Axis.2,
                 color = Location,
                 shape = Days),
             size = 5) +
  labs(title = 'GCMS all data (PCoA)',
       x = pcoa1_gc_all,
       y = pcoa2_gc_all) +
  scale_color_manual(values = c(Deep = '#977E59', Surface = '#4EB967')) +
  theme_bw() +
  theme(plot.title = element_text(face = 'bold', hjust = 0.5))

pcoa_gcms_plot_all


gcms_permanova_all <- adonis2(gcms_dist_all ~ Location * Treatment * Time, 
                              data = metadata_adonis,
                              method = 'bray',
                              permutations = 999) 
gcms_permanova_all

ggsave('gcms_permanova_all.svg', dpi = 300, height = 5, width = 6)  

# PCoA gcMS surface only ----

gcms_matrix_surface <- gcms_data %>% 
  filter(!str_detect(uniqueID, 'Deep')) %>% 
  column_to_rownames(var = 'uniqueID')

gcms_dist_surface <- vegdist(gcms_matrix_surface, method = 'bray')

PCoA_gcMS_surface <- pcoa(gcms_dist_surface, correction="cailliez", rn=NULL)

pcoa_gcms_scores_surface <- as.data.frame(PCoA_gcMS_surface$vectors) %>% 
  rownames_to_column(var = 'uniqueID') %>% 
  left_join(gcms_metadata, by = 'uniqueID')

pcoa1_gc_s <- paste0('PCo1 ', round(PCoA_gcMS_surface$values$Relative_eig[1] * 100, 2), '%')
pcoa2_gc_s <- paste0('PCo2 ', round(PCoA_gcMS_surface$values$Relative_eig[2] * 100, 2), '%')

pcoa_gcms_plot_surface <- pcoa_gcms_scores_surface %>% 
  ggplot() +
  geom_point(aes(x = Axis.1,
                 y = Axis.2,
                 color = Treatment,
                 shape = Days),
             size = 5) +
  labs(title = 'GCMS surface (PCoA)',
       x = pcoa1_gc_s,
       y = pcoa2_gc_s) +
  scale_color_manual(values = c(Control = '#EC9EC5', Glucose = '#14BCD8')) +
  theme_bw() +
  theme(plot.title = element_text(face = 'bold', hjust = 0.5))

pcoa_gcms_plot_surface


gcms_permanova_surface <- adonis2(gcms_dist_surface ~ Treatment * Time, 
                              data = metadata_adonis_surface,
                              method = 'bray',
                              permutations = 999) 
gcms_permanova_surface

ggsave('pcoa_gcms_plot_surface.svg', dpi = 300, height = 5, width = 6)  

# PCoA gcMS deep only ----

gcms_matrix_deep <- gcms_data %>% 
  filter(!str_detect(uniqueID, 'Surface')) %>% 
  column_to_rownames(var = 'uniqueID')

gcms_dist_deep <- vegdist(gcms_matrix_deep, method = 'bray')

PCoA_gcMS_deep <- pcoa(gcms_dist_deep, correction="cailliez", rn=NULL)

pcoa_gcms_scores_deep <- as.data.frame(PCoA_gcMS_deep$vectors) %>% 
  rownames_to_column(var = 'uniqueID') %>% 
  left_join(gcms_metadata, by = 'uniqueID')

pcoa1_gc_d <- paste0('PCo1 ', round(PCoA_gcMS_deep$values$Relative_eig[1] * 100, 2), '%')
pcoa2_gc_d <- paste0('PCo2 ', round(PCoA_gcMS_deep$values$Relative_eig[2] * 100, 2), '%')

pcoa_gcms_plot_deep <- pcoa_gcms_scores_deep %>% 
  ggplot() +
  geom_point(aes(x = Axis.1,
                 y = Axis.2,
                 color = Treatment,
                 shape = Days),
             size = 5) +
  labs(title = 'GCMS deep (PCoA)',
       x = pcoa1_gc_d,
       y = pcoa2_gc_d) +
  scale_color_manual(values = c(Control = '#EC9EC5', Glucose = '#14BCD8')) +
  theme_bw() +
  theme(plot.title = element_text(face = 'bold', hjust = 0.5))

pcoa_gcms_plot_deep


gcms_permanova_deep <- adonis2(gcms_dist_deep ~ Treatment * Time, 
                                  data = metadata_adonis_deep,
                                  method = 'bray',
                                  permutations = 999) 
gcms_permanova_deep


ggsave('pcoa_gcms_plot_deep.svg', dpi = 300, height = 5, width = 6)
