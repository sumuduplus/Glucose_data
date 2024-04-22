# Libraries ----

library(readxl)
library(ropls)
library(tidyverse)
library(ggpubr)
library(dplyr)
library(patchwork)


# load data -----

gcms_metadata <- read_xlsx("Raw_data/Metadata_all.xlsx", sheet = 1) %>% 
  unite(uniqueID, Location, Treatment, Days,sep = "_", remove = FALSE)
lcms_metadata <- read_xlsx("Raw_data/Metadata_all.xlsx", sheet = 2) %>% 
  unite(uniqueID, Location, Treatment, Days,sep = "_", remove = FALSE)
nmr_metadata <- read_xlsx("Raw_data/Metadata_all.xlsx", sheet = 3) %>% 
  unite(uniqueID, Location, Treatment, Days,sep = "_", remove = FALSE) 

gcms_data <- read_csv("Preprocess_data/GCMS_normalized.csv") %>% 
  left_join(select(gcms_metadata, SampleID, uniqueID), by = 'SampleID') %>% 
  select(-SampleID)

lcms_data <-read_csv("Preprocess_data/LCMS_normalized.csv") %>% 
  left_join(select(lcms_metadata, SampleID, uniqueID), by = 'SampleID') %>% 
  select(-SampleID)

nmr_data <- read_csv("Preprocess_data/NMR_normalized.csv") %>% 
  left_join(select(nmr_metadata, SampleID, uniqueID), by = 'SampleID') %>% 
  select(-SampleID)


# Divenn_input upload ----

Control_data <- read_csv("OPLS_da/control_diVenn_mets.csv")
glucose_data <- read_csv("OPLS_da/glucose_diVenn_mets.csv")

control_unique <- Control_data %>% 
  filter(!(FeatureID %in% glucose_data$FeatureID)) %>% 
  pull(FeatureID)

glucose_unique <- glucose_data %>% 
  filter(!(FeatureID %in% Control_data$FeatureID)) %>% 
  pull(FeatureID)

shared <- rbind(Control_data, glucose_data) %>% 
  select(FeatureID) %>% 
  distinct() %>% 
  filter(!FeatureID %in% c(glucose_unique, control_unique)) %>% 
  pull(FeatureID)

all_data <- gcms_data %>% 
  left_join(lcms_data, by = 'uniqueID') %>% 
  left_join(nmr_data, by = 'uniqueID') %>% 
  pivot_longer(!uniqueID, names_to = 'FeatureID', values_to = 'value') %>% 
  left_join(gcms_metadata, by = 'uniqueID') %>% 
  mutate(Days = factor(Days, levels = c('7', '14', '28', '42', '56', '70')),
         Location = factor(Location, levels = c('Surface', 'Deep')))

# Control_unique heatmap ----

control_abundances <- all_data %>% 
  filter(FeatureID %in% control_unique,
         Treatment == 'Control') %>% 
  group_by(FeatureID, Location, Days) %>% 
  summarise(mean_val = mean(value)) %>% 
  group_by(FeatureID) %>% 
  mutate(zscore = (mean_val - mean(mean_val)) / sd(mean_val)) 


control_order <- control_abundances %>% 
  filter(Location == 'Surface') %>% 
  group_by(FeatureID) %>% 
  summarise(zmean = mean(zscore)) %>% 
  arrange(zmean) %>% 
  mutate(order = n():1) %>% 
  select(FeatureID, order)

control_hmp <- control_abundances %>% 
  left_join(control_order, by = "FeatureID") %>% 
  ggplot() + ggtitle("unique metabolites of control treatment: across depth") +
  geom_tile(aes(x = fct_reorder(FeatureID, order),
                y = Days,
                fill = zscore),
            color = 'white') +
  facet_grid(rows = vars(Location)) +
  scale_fill_distiller(palette = 'RdBu') +

  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        axis.title.x = element_blank())

control_hmp  

# NO LC Control_unique heatmap ----

subset_control_nolc <- control_abundances %>% 
  filter(!str_detect(FeatureID, 'FT[0-9][0-9]'))


control_order_nolc <- subset_control_nolc %>% 
  filter(Location == 'Surface') %>% 
  group_by(FeatureID) %>% 
  summarise(zmean = mean(zscore)) %>% 
  arrange(zmean) %>% 
  mutate(order = n():1) %>% 
  select(FeatureID, order)

control_hmp_nolC <- subset_control_nolc %>% 
  left_join(control_order, by = "FeatureID") %>% 
  ggplot() + ggtitle("unique metabolites of control treatment: across depth no LC") +
  geom_tile(aes(y = fct_reorder(FeatureID, order),
                x = Days,
                fill = zscore),
            color = 'white') +
  facet_grid(cols = vars(Location)) +
  scale_fill_distiller(palette = 'RdBu') +
  
  theme_bw() 

control_hmp_nolC

ggsave('control_hmp_nolc.svg', dpi = 300, height = 6, width = 6)  

# Glucose_unique heatmap ----

glucose_abundances <- all_data %>% 
  filter(FeatureID %in% glucose_unique,
         Treatment == 'Glucose') %>% 
  group_by(FeatureID, Location, Days) %>% 
  summarise(mean_val = mean(value)) %>% 
  group_by(FeatureID) %>% 
  mutate(zscore = (mean_val - mean(mean_val)) / sd(mean_val)) 

subset_glucose_nolc <- glucose_abundances %>% 
  filter(!str_detect(FeatureID, 'FT[0-9][0-9]'))

glucose_order <- glucose_abundances %>% 
  filter(Location == 'Surface') %>% 
  group_by(FeatureID) %>% 
  summarise(zmean = mean(zscore)) %>% 
  arrange(zmean) %>% 
  mutate(order = n():1) %>% 
  select(FeatureID, order)

glucose_hmp <- glucose_abundances %>% 
  left_join(glucose_order, by = "FeatureID") %>% 
  ggplot() + ggtitle("unique metabolites of glucose treatment: across depth") +
  geom_tile(aes(x = fct_reorder(FeatureID, order),
                y = Days,
                fill = zscore),  color = 'white') +
  facet_grid(rows = vars(Location)) +
  scale_fill_distiller(palette = 'RdBu') +
  
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        axis.title.x = element_blank())

glucose_hmp 

# NO LC Glucose_unique heatmap ----

subset_glucose_nolc <- glucose_abundances %>% 
  filter(!str_detect(FeatureID, 'FT[0-9][0-9]'))


glucose_order_nolc <- subset_glucose_nolc %>% 
  filter(Location == 'Surface') %>% 
  group_by(FeatureID) %>% 
  summarise(zmean = mean(zscore)) %>% 
  arrange(zmean) %>% 
  mutate(order = n():1) %>% 
  select(FeatureID, order)

glucose_hmp_nolC <- subset_glucose_nolc %>% 
  left_join(glucose_order, by = "FeatureID") %>% 
  ggplot() + ggtitle("unique metabolites of glucose treatment: across depth no LC") +
  geom_tile(aes(y = fct_reorder(FeatureID, order),
                x = Days,
                fill = zscore),
            color = 'white') +
  facet_grid(cols = vars(Location)) +
  scale_fill_distiller(palette = 'RdBu') +
  
  theme_bw() 

glucose_hmp_nolC

ggsave('glucose_hmp_nolc.svg', dpi = 300, height = 6, width = 6)  

# Shared features heatmap ----

shared_abundances <- all_data %>% 
  filter(FeatureID %in% shared,
         Treatment == 'Control'|Treatment == 'Glucose') %>% 
  group_by(FeatureID, Treatment, Location, Days) %>% 
  summarise(mean_val = mean(value)) %>% 
  group_by(FeatureID) %>% 
  mutate(zscore = (mean_val - mean(mean_val)) / sd(mean_val)) 

shared_order <- shared_abundances %>% 
  filter(Location == 'Surface') %>% 
  group_by(FeatureID) %>% 
  summarise(zmean = mean(zscore)) %>% 
  arrange(zmean) %>% 
  mutate(order = n():1) %>% 
  select(FeatureID, order)

shared_hmp_control <- shared_abundances %>% 
  filter(Treatment == 'Control') %>% 
  left_join(shared_order, by = "FeatureID") %>% 
  ggplot() +
  labs(title = 'Control') +
  geom_tile(aes(y = fct_reorder(FeatureID, order),
                x = Days,
                fill = zscore),
            color = 'white') +
  facet_grid(rows = vars(Location)) +
  scale_fill_distiller(palette = 'RdBu', limits = c(-3.7, 3.7)) +
  
  theme_bw() +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_text(size = 3.5))

shared_hmp_control

shared_hmp_glucose <- shared_abundances %>% 
  filter(Treatment == 'Glucose') %>% 
  left_join(shared_order, by = "FeatureID") %>% 
  ggplot() +
  labs(title = 'Glucose') +
  geom_tile(aes(y = fct_reorder(FeatureID, order),
                x = Days,
                fill = zscore),
            color = 'white') +
  facet_grid(rows = vars(Location)) +
  scale_fill_distiller(palette = 'RdBu', limits = c(-3.7, 3.7)) +
  
  theme_bw() +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank())

shared_hmp_glucose


shared_hmp <- (shared_hmp_control | shared_hmp_glucose) +
  plot_layout(guides = 'collect')

shared_hmp

ggsave('shared_hmp.png', dpi = 300, height = 14, width = 10)  
