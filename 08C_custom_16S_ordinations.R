library(tidyverse)
library(readxl)
library(ggplot2)
library(ggrepel)
library(factoextra)
library(ape)
library(vegan)


microbial_metadata <- read_xlsx("Raw_data/Metadata_all.xlsx", sheet = 4) %>% 
  unite(uniqueID, Location, Treatment, TimePoint,SIPFraction, sep = "_", remove = FALSE)%>%
  column_to_rownames(var = 'SampleID')

microbial_data <- read_csv("16S/otu_chord_transformed_tax.csv")

microbial_tax <- microbial_data %>% 
  select(FeatureID, Kingdom, Phylum, Class, Order, Family, Genus, Species) %>% 
  mutate(name_tax_plot = case_when(
    !is.na(Genus) ~ Genus,
    !is.na(Family) ~ paste0('Un_tax f_', Family),
    !is.na(Order) ~ paste0('Un_tax o_', Order),
    !is.na(Class) ~ paste0('Un_tax c_', Class),
    !is.na(Phylum) ~ paste0('Un_tax p_', Phylum)
  ))

# Filtering 13C taxa out
Microbial_data_values <-  subset(microbial_data, select = -c(Kingdom,Phylum,Class,Order,Family,Genus,Species))%>%
column_to_rownames(var = 'FeatureID')%>% 
  select(-contains('13C'))


Microbial_data_values_t <- t(Microbial_data_values)

microbial_metadata_updated <- subset(microbial_metadata, select = -c(Depth,Location,TimePoint,SIPFraction,density,Days,Treatment))


All_w_out_13C <-  merge(Microbial_data_values_t, microbial_metadata_updated, by=0, all=F)%>%
  remove_rownames %>% column_to_rownames(var="uniqueID")

PCA_input <- All_w_out_13C[,!names(All_w_out_13C) %in% c("Row.names")]


# control vs glucose (no 13C) ----

pca_all <- prcomp(PCA_input, center = TRUE)
eigen_all <- get_eigenvalue(pca_all)
var_info_all <- get_pca_var(pca_all)

top_features_all <- as.data.frame(var_info_all$contrib) %>% 
  rownames_to_column(var = 'FeatureID') %>% 
  slice_max(Dim.1, n = 10)

top_features_coord_all <- as.data.frame(var_info_all$coord) %>% 
  rownames_to_column(var = 'FeatureID') %>%
  filter(FeatureID %in% top_features_all$FeatureID) %>% 
  inner_join(microbial_tax, by = 'FeatureID')

pca_coordinates_all <- as_tibble(pca_all$x)
pca_coordinates_all$uniqueID <- rownames(pca_all$x)

pca_coordinates_all_UP <- pca_coordinates_all %>% 
  left_join(microbial_metadata, by = 'uniqueID')%>% 
  mutate(SIPFraction_new = case_when(
    between(as.numeric(SIPFraction), 1, 6) ~ 'b_heavy (1-6)',
    between(as.numeric(SIPFraction), 7, 12) ~ 'a_light(7-12)',
    SIPFraction == 'S' ~ 'Unfractionated'
  ))


# labeling and formatting --

pc1_all <- paste0('PC1 (', round(eigen_all$variance.percent[1], digits = 1), '%)')
pc2_all <- paste0('PC1 (', round(eigen_all$variance.percent[2], digits = 1), '%)')

pca_plot_all <- pca_coordinates_all_UP %>% 
  ggplot() +
  geom_point(aes(x = PC1,
                 y = PC2,
                 color = Location,
                 shape = interaction(Treatment, Days),
             size = SIPFraction_new)) +
  scale_shape_manual(values = c(1, 2, 16, 17)) +
  scale_color_manual(values = c(Deep = '#977E59', Surface = '#4EB967')) +
  ggnewscale::new_scale_color() +
  geom_segment(data = top_features_coord_all,
               aes(xend = Dim.1,
                   yend = Dim.2,
                   group = FeatureID,
                   color = FeatureID),
               x = 0,
               y = 0,
               arrow = arrow(),
               show.legend = FALSE) +
  geom_text_repel(data = top_features_coord_all,
            aes(x = Dim.1 * 1.1,
                y = Dim.2 * 1.1,
                label = name_tax_plot,
                color = FeatureID),
            position = position_jitter(),
            show.legend = FALSE) +
  ggtitle("PCA control vs treatment") +
  labs(x = pc1_all,
       y = pc2_all) +
  theme_bw()

pca_plot_all  

ggsave('PCA all samples control vs treatment.svg', dpi = 300, height = 6, width = 8)  

#all samples (12C, 13C and control) ----

Microbial_data_values_all <-  subset(microbial_data, select = -c(Kingdom,Phylum,Class,Order,Family,Genus,Species))%>%
  column_to_rownames(var = 'FeatureID')

Microbial_data_values_t_all <- t(Microbial_data_values_all)

microbial_metadata_updated_all <- subset(microbial_metadata, select = -c(Depth,Location,TimePoint,SIPFraction,density,Days,Treatment))


All_everything <-  merge(Microbial_data_values_t_all, microbial_metadata_updated_all, by=0, all=F)%>%
  remove_rownames %>% column_to_rownames(var="uniqueID")

PCA_input_everything <- All_everything[,!names(All_everything) %in% c("Row.names")]

# PCA all samples (12C, 13C and control) ----

pca_everything <- prcomp(PCA_input_everything, center = TRUE)
eigen_everything <- get_eigenvalue(pca_everything)

var_info_everything <- get_pca_var(pca_everything)

top_features_everything <- as.data.frame(var_info_everything$contrib) %>% 
  rownames_to_column(var = 'FeatureID') %>% 
  slice_max(Dim.1, n = 10)

top_features_coord_everything <- as.data.frame(var_info_everything$coord) %>% 
  rownames_to_column(var = 'FeatureID') %>%
  filter(FeatureID %in% top_features_everything$FeatureID) %>% 
  inner_join(microbial_tax, by = 'FeatureID')

pca_coordinates_everything <- as_tibble(pca_everything$x)
pca_coordinates_everything$uniqueID <- rownames(pca_everything$x)

pca_coordinates_everything <- pca_coordinates_everything %>% 
  left_join(microbial_metadata, by = 'uniqueID')%>% 
  mutate(SIPFraction_new = case_when(
    between(as.numeric(SIPFraction), 1, 6) ~ 'b_heavy (1-6)',
    between(as.numeric(SIPFraction), 7, 12) ~ 'a_light(7-12)',
    SIPFraction == 'S' ~ 'Unfractionated'
  ))

pc1_everything <- paste0('PC1 (', round(eigen_everything$variance.percent[1], digits = 1), '%)')
pc2_everything <- paste0('PC1 (', round(eigen_everything$variance.percent[2], digits = 1), '%)')

pca_plot_everything <- pca_coordinates_everything %>% 
  ggplot() +
  geom_point(aes(x = PC1,
                 y = PC2,
                 color = Location,
                 shape = interaction(Treatment, Days),
             size = SIPFraction_new)) +
  scale_shape_manual(values = c(1, 2, 9, 10, 17,18)) +
  scale_color_manual(values = c(Deep = '#977E59', Surface = '#4EB967')) +
  ggnewscale::new_scale_color() +
  geom_segment(data = top_features_coord_everything,
               aes(xend = Dim.1,
                   yend = Dim.2,
                   group = FeatureID,
                   color = FeatureID),
               x = 0,
               y = 0,
               arrow = arrow(),
               show.legend = FALSE) +
  geom_text_repel(data = top_features_coord_everything,
                  aes(x = Dim.1 * 1.1,
                      y = Dim.2 * 1.1,
                      label = name_tax_plot,
                      color = FeatureID),
                  position = position_jitter(),
                  show.legend = FALSE) +
  ggtitle("PCA all microbial samples") +
  labs(x = pc1_everything,
       y = pc2_everything) +
  theme_bw()

pca_plot_everything  

ggsave('PCA all microbial samples.svg', dpi = 300, height = 6, width = 8) 


# 12c vs 13C (glucose only)----


Microbial_data_values_both_gl <-  subset(microbial_data, select = -c(Kingdom,Phylum,Class,Order,Family,Genus,Species))%>%
  column_to_rownames(var = 'FeatureID')%>% 
  select(-contains('Una'))


Microbial_data_values_t_both_gl <- t(Microbial_data_values_both_gl)

microbial_metadata_updated_both_gl <- subset(microbial_metadata, select = -c(Depth,Location,TimePoint,SIPFraction,density,Days,Treatment))


Both_glu <-  merge(Microbial_data_values_t_both_gl, microbial_metadata_updated_both_gl, by=0, all=F)%>%
  remove_rownames %>% column_to_rownames(var="uniqueID")

PCA_input_both_glu <- Both_glu[,!names(Both_glu) %in% c("Row.names")]


# 12C vs 13C surface (both glucose) only ----

surface_both_Glu__only <- PCA_input_both_glu %>%
  filter(grepl("Surface", rownames(PCA_input_both_glu)))

pca_surface_both_Glu <- prcomp(surface_both_Glu__only, center = TRUE)
eigen_surface_both_Glu <- get_eigenvalue(pca_surface_both_Glu)

var_info_surface_both_Glu <- get_pca_var(pca_surface_both_Glu)

top_features_surface_both_Glu <- as.data.frame(var_info_surface_both_Glu$contrib) %>% 
  rownames_to_column(var = 'FeatureID') %>% 
  slice_max(Dim.1, n = 10)


top_features_coord_surface_both_Glu <- as.data.frame(var_info_surface_both_Glu$coord) %>% 
  rownames_to_column(var = 'FeatureID') %>%
  filter(FeatureID %in% top_features_surface_both_Glu$FeatureID) %>% 
  inner_join(microbial_tax, by = 'FeatureID')


pca_coordinates_surface_both_Glu <- as_tibble(pca_surface_both_Glu$x)
pca_coordinates_surface_both_Glu$uniqueID <- rownames(pca_surface_both_Glu$x)

pca_coordinates_surface_both_Glu <- pca_coordinates_surface_both_Glu %>% 
  left_join(microbial_metadata, by = 'uniqueID')%>% 
mutate(SIPFraction_new = case_when(
  between(as.numeric(SIPFraction), 1, 6) ~ 'b_heavy (1-6)',
  between(as.numeric(SIPFraction), 7, 12) ~ 'a_light(7-12)',
  SIPFraction == 'S' ~ 'Unfractionated'
))

# labeling and formatting --

pc1_surface_both_Glu <- paste0('PC1 (', round(eigen_surface_both_Glu$variance.percent[1], digits = 1), '%)')
pc2_surface_both_Glu <- paste0('PC1 (', round(eigen_surface_both_Glu$variance.percent[2], digits = 1), '%)')

pca_plot_surface_both_Glu <- pca_coordinates_surface_both_Glu %>% 
  ggplot() +
  geom_point(aes(x = PC1,
                 y = PC2,
                 color = Treatment,
                 shape = Days,
             size = SIPFraction_new)) +
  scale_color_manual(values = c(Glucose_13C = '#ADD8E6', Glucose = '#14BCD8')) +
  ggnewscale::new_scale_color() +
  geom_segment(data = top_features_coord_surface_both_Glu,
               aes(xend = Dim.1,
                   yend = Dim.2,
                   group = FeatureID,
                   color = FeatureID),
               x = 0,
               y = 0,
               arrow = arrow(),
               show.legend = FALSE) +
  geom_text_repel(data = top_features_coord_surface_both_Glu,
                  aes(x = Dim.1 * 1.1,
                      y = Dim.2 * 1.1,
                      label = name_tax_plot,
                      color = FeatureID),
                  position = position_jitter(),
                  show.legend = FALSE) +
  ggtitle("PCA surface 12C vs 13C") +
  labs(x = pc1_surface_both_Glu,
       y = pc2_surface_both_Glu) +
  theme_bw()

pca_plot_surface_both_Glu

ggsave('12C vs 13C surface PCA.svg', dpi = 300, height = 6, width = 8) 

# Deep samples 12C vs 13C (glucose only) ----

deep_both_Glu__only <- PCA_input_both_glu %>%
  filter(grepl("Deep", rownames(PCA_input_both_glu)))

pca_deep_both_Glu <- prcomp(deep_both_Glu__only, center = TRUE)
eigen_deep_both_Glu <- get_eigenvalue(pca_deep_both_Glu)

var_info_deep_both_Glu <- get_pca_var(pca_deep_both_Glu)

top_features_deep_both_Glu <- as.data.frame(var_info_deep_both_Glu$contrib) %>% 
  rownames_to_column(var = 'FeatureID') %>% 
  slice_max(Dim.1, n = 10)


top_features_coord_deep_both_Glu <- as.data.frame(var_info_deep_both_Glu$coord) %>% 
  rownames_to_column(var = 'FeatureID') %>%
  filter(FeatureID %in% top_features_deep_both_Glu$FeatureID) %>% 
  inner_join(microbial_tax, by = 'FeatureID')


pca_coordinates_deep_both_Glu <- as_tibble(pca_deep_both_Glu$x)
pca_coordinates_deep_both_Glu$uniqueID <- rownames(pca_deep_both_Glu$x)

pca_coordinates_deep_both_Glu <- pca_coordinates_deep_both_Glu %>% 
  left_join(microbial_metadata, by = 'uniqueID')%>% 
  mutate(SIPFraction_new = case_when(
    between(as.numeric(SIPFraction), 1, 6) ~ 'b_heavy (1-6)',
    between(as.numeric(SIPFraction), 7, 12) ~ 'a_light(7-12)',
    SIPFraction == 'S' ~ 'Unfractionated'
  ))

# labeling and formatting --

pc1_deep_both_Glu <- paste0('PC1 (', round(eigen_deep_both_Glu$variance.percent[1], digits = 1), '%)')
pc2_deep_both_Glu <- paste0('PC1 (', round(eigen_deep_both_Glu$variance.percent[2], digits = 1), '%)')

pca_plot_deep_both_Glu <- pca_coordinates_deep_both_Glu %>% 
  ggplot() +
  geom_point(aes(x = PC1,
                 y = PC2,
                 color = Treatment,
                 shape = Days,
                 size = SIPFraction_new)) +
  scale_color_manual(values = c(Glucose_13C = '#ADD8E6', Glucose = '#14BCD8')) +
  ggnewscale::new_scale_color() +
  geom_segment(data = top_features_coord_deep_both_Glu,
               aes(xend = Dim.1,
                   yend = Dim.2,
                   group = FeatureID,
                   color = FeatureID),
               x = 0,
               y = 0,
               arrow = arrow(),
               show.legend = FALSE) +
  geom_text_repel(data = top_features_coord_deep_both_Glu,
                  aes(x = Dim.1 * 1.1,
                      y = Dim.2 * 1.1,
                      label = name_tax_plot,
                      color = FeatureID),
                  position = position_jitter(),
                  show.legend = FALSE) +
  ggtitle("PCA deep 12C vs 13C") +
  labs(x = pc1_deep_both_Glu,
       y = pc2_deep_both_Glu) +
  theme_bw()

pca_plot_deep_both_Glu

ggsave('12C vs 13C deep PCA.svg', dpi = 300, height = 6, width = 8)

write.csv(Microbial_data_values, 'OPLS_da/16S_fractionated_for_OPLS.csv')


