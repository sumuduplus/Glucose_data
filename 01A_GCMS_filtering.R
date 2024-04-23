# libraries ----

library(readxl)
library(tidyverse)

# loading data ----

metadata <- read_xlsx('Raw_data/Metadata_all.xlsx', sheet = 1)

gcms_raw_data <- read_csv("Raw_data/Raw_data_GCMS.csv") %>% 
  select(-contains('Unknown'))

# Filtering ----

gcms_long <- gcms_raw_data %>% 
  pivot_longer(!SampleID, names_to = 'FeatureID', 
               values_to = 'abundance') %>% 
  mutate(abundance = ifelse(is.na(abundance), 0, abundance)) %>% 
  inner_join(metadata, by = 'SampleID')

gcms_filtered <- gcms_long %>%
  filter(Location == "Surface" | Location == "Deep")

# Counting features for filtering

gcms_counts <- gcms_filtered %>% 
  filter(abundance > 0) %>% 
  group_by(Location) %>% 
  count(FeatureID) %>% 
  pivot_wider(names_from = 'Location', 
              values_from = 'n',
              values_fill = 0)


# filtering cutoff within each depth (0.25 used)

features_for_filtering <- gcms_counts %>% 
  mutate(across(2:3, ~ .x / 12)) %>% 
  filter(Deep >= 0.25 | Surface >= 0.25) %>% 
  pull(FeatureID)

gcms_filtered_updated <- gcms_filtered %>% 
  filter(FeatureID %in% features_for_filtering)

# Normalization ----

# Pseudo count added cube root transformation

pseudo <- min(gcms_filtered_updated$abundance[gcms_filtered_updated$abundance > 0])

gcms_normalized <- gcms_filtered_updated %>%  
  mutate(normalized_abundance = (abundance + pseudo)^(1/3))

# scaling

gcms_scaled <- gcms_normalized %>% 
  select(FeatureID, SampleID, normalized_abundance) %>% 
  group_by(SampleID) %>% 
  mutate(median = median(normalized_abundance)) %>% 
  ungroup() %>% 
  mutate(scaled_normalized_abundance = 
           normalized_abundance  / median * mean(median)) %>% 
  select(-median, -normalized_abundance)

gcms_boxplot <- gcms_scaled %>% 
  ggplot() +
  geom_boxplot(aes(x = SampleID,
                   y = scaled_normalized_abundance)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

gcms_boxplot


# Saving data ----

gcms_scaled_wider <- gcms_scaled %>% 
  pivot_wider(names_from = 'FeatureID', 
              values_from = 'scaled_normalized_abundance')

write_csv(gcms_scaled_wider, 'Preprocess_data/GCMS_normalized.csv')
