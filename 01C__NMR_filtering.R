# libraries ----

library(readxl)
library(tidyverse)

# loading data ----

metadata <- read_xlsx('Raw_data/Metadata_all.xlsx', sheet = 3)
nmr_raw_data <- read_csv("Raw_data/Raw_data_NMR.csv")

# Normalization ----

nmr_long <- nmr_raw_data %>% 
  pivot_longer(!SampleID, names_to = 'FeatureID', 
               values_to = 'abundance') 

# Pseudo count added cube root transformation

pseudo <- min(nmr_long$abundance[nmr_long$abundance > 0])

nmr_normalized <- nmr_long %>%  
  mutate(normalized_abundance = (abundance + pseudo)^(1/3))

# Scaling

nmr_scaled <- nmr_normalized %>% 
  select(FeatureID, SampleID, normalized_abundance) %>% 
  group_by(SampleID) %>% 
  mutate(median = median(normalized_abundance)) %>% 
  ungroup() %>% 
  mutate(scaled_normalized_abundance = 
           normalized_abundance  / median * mean(median)) %>% 
  select(-median, -normalized_abundance)

nmr_boxplot <- nmr_scaled %>% 
  ggplot() +
  geom_boxplot(aes(x = SampleID,
                   y = scaled_normalized_abundance)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

nmr_boxplot

# Saving data ----

nmr_scaled_wider <- nmr_scaled %>% 
  pivot_wider(names_from = 'FeatureID', 
              values_from = 'scaled_normalized_abundance')

write_csv(nmr_scaled_wider, 'Preprocess_data/NMR_normalized.csv')
