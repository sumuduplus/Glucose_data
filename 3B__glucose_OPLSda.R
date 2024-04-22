# Libraries ----

library(readxl)
library(ropls)
library(tidyverse)
library(dplyr)


# load data -----

gcms_metadata <- read_xlsx("Raw_data/Metadata_all.xlsx", sheet = 1) %>% 
  unite(uniqueID, Location, Treatment, Days,sep = "_", remove = FALSE)
lcms_metadata <- read_xlsx("Raw_data/Metadata_all.xlsx", sheet = 2) %>% 
  unite(uniqueID, Location, Treatment, Days,sep = "_", remove = FALSE)
nmr_metadata <- read_xlsx("Raw_data/Metadata_all.xlsx", sheet = 3) %>% 
  unite(uniqueID, Location, Treatment, Days,sep = "_", remove = FALSE) 

gcms_data <- read_csv("Preprocess_data/GCMS_normalized.csv")

lcms_data <-read_csv("Preprocess_data/LCMS_normalized.csv")

nmr_data <- read_csv("Preprocess_data/NMR_normalized.csv")

# GCMS_oplsda ----

gcms_all <- gcms_data %>% 
  inner_join(gcms_metadata, by = 'SampleID') 

gcms_glucose <- gcms_all %>%
  filter(Treatment == "Glucose")

gcms_glucose_class <- as.factor (gcms_glucose$Location) 
gcms_glucose_data <- gcms_glucose[, !colnames(gcms_glucose) %in% colnames(gcms_metadata)]

gcms_opls <- opls (gcms_glucose_data, gcms_glucose_class, predI = 1, crossvalI = 11, orthoI = NA, permI = 1000) 
gcms_VIPS <- as.data.frame(getVipVn(gcms_opls,orthoL = FALSE)) %>% 
  rownames_to_column(var = 'FeatureID') %>% 
  rename(VIP = 2)

gcms_glucose_t <- gcms_data %>% 
  filter(!str_detect(SampleID, '_C_')) %>% 
  column_to_rownames(var = 'SampleID') %>% 
  t()

gc_glucose_deep <- gcms_glucose_t[, 1:6]
gc_glucose_surface <- gcms_glucose_t[, 7:12]

gcms_da_res <- tibble(FeatureID = NA, comp_abundance = NA, .rows = nrow(gc_glucose_deep))


for(i in 1:nrow(gcms_da_res)){
  
  gcms_da_res$FeatureID[i] <- rownames(gc_glucose_deep)[i]
  
  if(mean(gc_glucose_deep[i,])> mean(gc_glucose_surface[i,])){
    gcms_da_res$comp_abundance[i] <- 'Deep_abundant'
  } else {
    gcms_da_res$comp_abundance[i] <- 'Surface_abundant'
  }
}

# gcms opls-da_results filtering ----

gcms_final_res <- gcms_da_res %>% 
  left_join(gcms_VIPS, by = 'FeatureID')

gcms_final_res_filtered <- gcms_final_res %>% 
  filter(VIP > 1)

write_csv(gcms_final_res_filtered, 'OPLS_da/GCMS_opls_results_glucose.csv')


# LCMS_oplsda ----

lcms_all <- lcms_data %>% 
  inner_join(lcms_metadata, by = 'SampleID') 

lcms_glucose <- lcms_all %>%
  filter(Treatment == "Glucose")

lcms_glucose_class <- as.factor (lcms_glucose$Location) 
lcms_glucose_data <- lcms_glucose[, !colnames(lcms_glucose) %in% colnames(lcms_metadata)]

lcms_opls <- opls (lcms_glucose_data, lcms_glucose_class, predI = 1, crossvalI = 11, orthoI = NA, permI = 1000) 
lcms_VIPS <- as.data.frame(getVipVn(lcms_opls,orthoL = FALSE)) %>% 
  rownames_to_column(var = 'FeatureID') %>% 
  rename(VIP = 2)

lcms_glucose_t <- lcms_data %>% 
  filter(!str_detect(SampleID, '_C')) %>% 
  column_to_rownames(var = 'SampleID') %>% 
  t()

lc_glucose_surface <- lcms_glucose_t[, 1:6]
lc_glucose_deep <- lcms_glucose_t[, 7:12]

lcms_da_res <- tibble(FeatureID = NA, comp_abundance = NA, .rows = nrow(lc_glucose_deep))

for(i in 1:nrow(lcms_da_res)){
  
  lcms_da_res$FeatureID[i] <- rownames(lc_glucose_deep)[i]
  
  if(mean(lc_glucose_deep[i,])> mean(lc_glucose_surface[i,])){
    lcms_da_res$comp_abundance[i] <- 'Deep_abundant'
  } else {
    lcms_da_res$comp_abundance[i] <- 'Surface_abundant'
  }
}

# lcms opls-da_results filtering ----

lcms_final_res <- lcms_da_res %>% 
  left_join(lcms_VIPS, by = 'FeatureID')

lcms_final_res_filtered <- lcms_final_res %>% 
  filter(VIP > 1)

write_csv(lcms_final_res_filtered, 'OPLS_da/LCMS_opls_results_glucose.csv')


# NMR_oplsda ----

nmr_all <- nmr_data %>% 
  inner_join(nmr_metadata, by = 'SampleID') 

nmr_glucose <- nmr_all %>%
  filter(Treatment == "Glucose")

nmr_glucose_class <- as.factor (nmr_glucose$Location) 
nmr_glucose_data <- nmr_glucose[, !colnames(nmr_glucose) %in% colnames(nmr_metadata)]

nmr_opls <- opls (nmr_glucose_data, nmr_glucose_class, predI = 1, crossvalI = 11, orthoI = NA, permI = 1000) 
nmr_VIPS <- as.data.frame(getVipVn(nmr_opls,orthoL = FALSE)) %>% 
  rownames_to_column(var = 'FeatureID') %>% 
  rename(VIP = 2)

nmr_glucose_t <- nmr_data %>% 
  filter(!str_detect(SampleID, '_CT_')) %>% 
  column_to_rownames(var = 'SampleID') %>% 
  t()

nmr_glucose_surface <- nmr_glucose_t[, 1:6]
nmr_glucose_deep <- nmr_glucose_t[, 7:12]

nmr_da_res <- tibble(FeatureID = NA, comp_abundance = NA, .rows = nrow(nmr_glucose_deep))

for(i in 1:nrow(nmr_da_res)){
  
  nmr_da_res$FeatureID[i] <- rownames(nmr_glucose_deep)[i]
  
  if(mean(nmr_glucose_deep[i,])> mean(nmr_glucose_surface[i,])){
    nmr_da_res$comp_abundance[i] <- 'Deep_abundant'
  } else {
    nmr_da_res$comp_abundance[i] <- 'Surface_abundant'
  }
}

# nmr opls-da_results filtering ----

nmr_final_res <- nmr_da_res %>% 
  left_join(nmr_VIPS, by = 'FeatureID')

nmr_final_res_filtered <- nmr_final_res %>% 
  filter(VIP > 1)

write_csv(nmr_final_res_filtered, 'OPLS_da/NMR_opls_results_glucose.csv')

# Di_venn input ----

Glucose_divenn_mets <- rbind(gcms_final_res_filtered,lcms_final_res_filtered, nmr_final_res_filtered)
write_csv(Glucose_divenn_mets, 'OPLS_da/glucose_diVenn_mets.csv')
