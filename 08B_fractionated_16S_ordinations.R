library(tidyverse)
library(readxl)
library(ggplot2)
library(ggrepel)
library(factoextra)


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


Microbial_data_values_fractionated <-  subset(microbial_data, select = -c(Kingdom,Phylum,Class,Order,Family,Genus,Species))%>%
  column_to_rownames(var = 'FeatureID')%>% 
  select(-contains('S'))

Microbial_data_values_t_frac <- t(Microbial_data_values_fractionated)

microbial_metadata_updated_frac <- subset(microbial_metadata, select = -c(Depth,Location,TimePoint,SIPFraction,density,Days,Treatment))

All_fractionated <-  merge(Microbial_data_values_t_frac, microbial_metadata_updated_frac, by=0, all=F)%>%
  remove_rownames %>% column_to_rownames(var="uniqueID")

PCA_input_frac <- All_fractionated[,!names(All_fractionated) %in% c("Row.names")]


# Filtering fractionated data

All_fractionated <-  merge(Microbial_data_values_t_frac, microbial_metadata_updated_frac, by=0, all=F)%>%
  remove_rownames %>% column_to_rownames(var="uniqueID")

PCA_input <- All_fractionated[,!names(All_fractionated) %in% c("Row.names")]



# all samples ----

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

pca_coordinates_all_final <- pca_coordinates_all%>% 
  left_join(microbial_metadata, by = 'uniqueID')

# labeling and formatting --

pc1_all <- paste0('PC1 (', round(eigen_all$variance.percent[1], digits = 1), '%)')
pc2_all <- paste0('PC1 (', round(eigen_all$variance.percent[2], digits = 1), '%)')

pca_plot_all <- pca_coordinates_all_final %>% 
  ggplot() +
  geom_point(aes(x = PC1,
                 y = PC2,
                 color = Location,
                 shape = Treatment, size = SIPFraction)) +
  scale_color_manual(values = c(Deep = '#977E59', Surface = '#4EB967')) +
  ggtitle("PCA All fractionated") +
  labs(x = pc1_all,
       y = pc2_all) +
  theme_bw()

pca_plot_all  

# surface samples ----

Surface_only_frac <- PCA_input %>%
  filter(grepl("Surface", rownames(All_fractionated)))

pca_Surface <- prcomp(Surface_only_frac, center = TRUE)
eigen_Surface <- get_eigenvalue(pca_Surface)

var_info_surface <- get_pca_var(pca_Surface)

top_features_surface <- as.data.frame(var_info_surface$contrib) %>% 
  rownames_to_column(var = 'FeatureID') %>% 
  slice_max(Dim.1, n = 10)


top_features_coord_surface <- as.data.frame(var_info_surface$coord) %>% 
  rownames_to_column(var = 'FeatureID') %>%
  filter(FeatureID %in% top_features_surface$FeatureID) %>% 
  inner_join(microbial_tax, by = 'FeatureID')


pca_coordinates_Surface <- as_tibble(pca_Surface$x)
pca_coordinates_Surface$uniqueID <- rownames(pca_Surface$x)

pca_coordinates_Surface <- pca_coordinates_Surface %>% 
  left_join(microbial_metadata, by = 'uniqueID')%>% 
  mutate(SIPFraction = as.numeric(SIPFraction))

# labeling and formatting --

pc1_Surface <- paste0('PC1 (', round(eigen_Surface$variance.percent[1], digits = 1), '%)')
pc2_Surface <- paste0('PC1 (', round(eigen_Surface$variance.percent[2], digits = 1), '%)')

pca_plot_Surface <- pca_coordinates_Surface %>% 
  ggplot() +
  geom_point(aes(x = PC1,
                 y = PC2,
                 color = Treatment,
                 shape = Days,
             size = SIPFraction)) +
  scale_color_manual(values = c(Control = '#EC9EC5', Glucose = '#14BCD8', Glucose_13C = '#1C4966')) +
  ggnewscale::new_scale_color() +
  geom_segment(data = top_features_coord_surface,
               aes(xend = Dim.1,
                   yend = Dim.2,
                   group = FeatureID,
                   color = FeatureID),
               x = 0,
               y = 0,
               arrow = arrow(),
               show.legend = FALSE) +
  geom_text_repel(data = top_features_coord_surface,
                  aes(x = Dim.1 * 1.1,
                      y = Dim.2 * 1.1,
                      label = name_tax_plot,
                      color = FeatureID),
                  position = position_jitter(),
                  show.legend = FALSE) +
  ggtitle("PCA Surface fractionated") +
  labs(x = pc1_Surface,
       y = pc2_Surface) +
  theme_bw()

pca_plot_Surface 

# Deep samples ----

deep_only_frac <- PCA_input %>%
  filter(grepl("Deep", rownames(All_fractionated)))

pca_deep <- prcomp(deep_only_frac, center = TRUE)
eigen_deep <- get_eigenvalue(pca_deep)

var_info_deep <- get_pca_var(pca_deep)

top_features_deep <- as.data.frame(var_info_deep$contrib) %>% 
  rownames_to_column(var = 'FeatureID') %>% 
  slice_max(Dim.1, n = 10)


top_features_coord_deep <- as.data.frame(var_info_deep$coord) %>% 
  rownames_to_column(var = 'FeatureID') %>%
  filter(FeatureID %in% top_features_deep$FeatureID) %>% 
  inner_join(microbial_tax, by = 'FeatureID')


pca_coordinates_deep <- as_tibble(pca_deep$x)
pca_coordinates_deep$uniqueID <- rownames(pca_deep$x)

pca_coordinates_deep <- pca_coordinates_deep %>% 
  left_join(microbial_metadata, by = 'uniqueID')%>% 
  mutate(SIPFraction = as.numeric(SIPFraction))

# labeling and formatting --

pc1_deep <- paste0('PC1 (', round(eigen_deep$variance.percent[1], digits = 1), '%)')
pc2_deep <- paste0('PC1 (', round(eigen_deep$variance.percent[2], digits = 1), '%)')

pca_plot_deep <- pca_coordinates_deep %>% 
  ggplot() +
  geom_point(aes(x = PC1,
                 y = PC2,
                 color = Treatment,
                 shape = Days,
             size = SIPFraction)) +
  scale_color_manual(values = c(Control = '#EC9EC5', Glucose = '#14BCD8', Glucose_13C = '#1C4966')) +
  ggnewscale::new_scale_color() +
  geom_segment(data = top_features_coord_deep,
               aes(xend = Dim.1,
                   yend = Dim.2,
                   group = FeatureID,
                   color = FeatureID),
               x = 0,
               y = 0,
               arrow = arrow(),
               show.legend = FALSE) +
  geom_text_repel(data = top_features_coord_deep,
                  aes(x = Dim.1 * 1.1,
                      y = Dim.2 * 1.1,
                      label = name_tax_plot,
                      color = FeatureID),
                  position = position_jitter(),
                  show.legend = FALSE) +
  ggtitle("PCA deep fractionated") +
  labs(x = pc1_deep,
       y = pc2_deep) +
  theme_bw()

pca_plot_deep 

