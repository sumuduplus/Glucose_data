# libraries ----

library(readxl)
library(tidyverse)

# loading data ----

metadata <- read_xlsx('Raw_data/Metadata_all.xlsx', sheet = 2)
LC_compound_annotations <- read_xlsx('Raw_data/Metadata_all.xlsx', sheet = 5)


lcms_raw_data <- read_csv("Raw_data/Raw_data_LCMS.csv") 

# Filtering ----

lcms_long <- lcms_raw_data %>% 
  pivot_longer(!SampleID, names_to = 'FeatureID', 
               values_to = 'abundance') %>% 
  mutate(abundance = ifelse(is.na(abundance), 0, abundance)) %>% 
  inner_join(metadata, by = 'SampleID')

# Counting features for filtering ----

lcms_counts <- lcms_long %>% 
  filter(abundance > 0) %>% 
  group_by(Location) %>% 
  count(FeatureID) %>% 
  pivot_wider(names_from = 'Location', 
              values_from = 'n',
              values_fill = 0)


# filtering cutoff across all samples (0.5 used)

features_for_filtering <- lcms_counts %>% 
  mutate(across(2:3, ~ .x / 12)) %>% 
  filter(Deep >= 0.5 | Surface >= 0.5) %>% 
  pull(FeatureID)


lcms_filtered <- lcms_long %>% 
  filter(FeatureID %in% features_for_filtering)

# Normalization ----

# Pseudo count added cube root transformation

pseudo <- min(lcms_filtered$abundance[lcms_filtered$abundance > 0])

lcms_normalized <- lcms_filtered %>%  
  mutate(normalized_abundance = (abundance + pseudo)^(1/3))

# Mean centering
lcms_scaled <- lcms_normalized %>% 
  select(FeatureID, SampleID, normalized_abundance) %>% 
  group_by(SampleID) %>% 
  mutate(median = median(normalized_abundance)) %>% 
  ungroup() %>% 
  mutate(scaled_normalized_abundance = 
           normalized_abundance  / median * mean(median)) %>% 
  select(-median, -normalized_abundance)


lcms_boxplot <- lcms_scaled %>% 
  ggplot() +
  geom_boxplot(aes(x = SampleID,
                   y = scaled_normalized_abundance)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

lcms_boxplot


# Saving data ----

lcms_scaled_wider <- lcms_scaled %>% 
  pivot_wider(names_from = 'FeatureID', 
              values_from = 'scaled_normalized_abundance')

write_csv(lcms_scaled_wider, 'Preprocess_data/LCMS_normalized.csv')
