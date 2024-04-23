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

gcms_control <- gcms_all %>%
  filter(Treatment == "Control")

gcms_control_class <- as.factor (gcms_control$Location) 
gcms_control_data <- gcms_control[, !colnames(gcms_control) %in% colnames(gcms_metadata)]

gcms_opls <- opls (gcms_control_data, gcms_control_class, predI = 1, crossvalI = 11, orthoI = NA, permI = 1000) 
gcms_VIPS <- as.data.frame(getVipVn(gcms_opls,orthoL = FALSE)) %>% 
  rownames_to_column(var = 'FeatureID') %>% 
  rename(VIP = 2)

gcms_control_t <- gcms_data %>% 
  filter(!str_detect(SampleID, '_U_')) %>% 
  column_to_rownames(var = 'SampleID') %>% 
  t()

gc_control_deep <- gcms_control_t[, 1:6]
gc_control_surface <- gcms_control_t[, 7:12]

gcms_da_res <- tibble(FeatureID = NA, comp_abundance = NA, .rows = nrow(gc_control_deep))


for(i in 1:nrow(gcms_da_res)){
  
  gcms_da_res$FeatureID[i] <- rownames(gc_control_deep)[i]
  
  if(mean(gc_control_deep[i,])> mean(gc_control_surface[i,])){
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

write_csv(gcms_final_res_filtered, 'OPLS_da/GCMS_opls_results_control.csv')

# LCMS_oplsda ----

lcms_all <- lcms_data %>% 
  inner_join(lcms_metadata, by = 'SampleID') 

lcms_control <- lcms_all %>%
  filter(Treatment == "Control")

lcms_control_class <- as.factor (lcms_control$Location) 
lcms_control_data <- lcms_control[, !colnames(lcms_control) %in% colnames(lcms_metadata)]

lcms_opls <- opls (lcms_control_data, lcms_control_class, predI = 1, crossvalI = 11, orthoI = NA, permI = 1000) 
lcms_VIPS <- as.data.frame(getVipVn(lcms_opls,orthoL = FALSE)) %>% 
  rownames_to_column(var = 'FeatureID') %>% 
  rename(VIP = 2)

lcms_control_t <- lcms_data %>% 
  filter(!str_detect(SampleID, '_U')) %>% 
  column_to_rownames(var = 'SampleID') %>% 
  t()

lc_control_surface <- lcms_control_t[, 1:6]
lc_control_deep <- lcms_control_t[, 7:12]

lcms_da_res <- tibble(FeatureID = NA, comp_abundance = NA, .rows = nrow(lc_control_deep))

for(i in 1:nrow(lcms_da_res)){
  
  lcms_da_res$FeatureID[i] <- rownames(lc_control_deep)[i]
  
  if(mean(lc_control_deep[i,])> mean(lc_control_surface[i,])){
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

write_csv(lcms_final_res_filtered, 'OPLS_da/LCMS_opls_results_control.csv')

# NMR_oplsda ----

nmr_all <- nmr_data %>% 
  inner_join(nmr_metadata, by = 'SampleID') 

nmr_control <- nmr_all %>%
  filter(Treatment == "Control")

nmr_control_class <- as.factor (nmr_control$Location) 
nmr_control_data <- nmr_control[, !colnames(nmr_control) %in% colnames(nmr_metadata)]

nmr_opls <- opls (nmr_control_data, nmr_control_class, predI = 1, crossvalI = 11, orthoI = NA, permI = 1000) 
nmr_VIPS <- as.data.frame(getVipVn(nmr_opls,orthoL = FALSE)) %>% 
  rownames_to_column(var = 'FeatureID') %>% 
  rename(VIP = 2)

nmr_control_t <- nmr_data %>% 
  filter(!str_detect(SampleID, '_Gl_')) %>% 
  column_to_rownames(var = 'SampleID') %>% 
  t()

nmr_control_surface <- nmr_control_t[, 1:6]
nmr_control_deep <- nmr_control_t[, 7:12]


nmr_da_res <- tibble(FeatureID = NA, comp_abundance = NA, .rows = nrow(nmr_control_deep))

for(i in 1:nrow(nmr_da_res)){
  
  nmr_da_res$FeatureID[i] <- rownames(nmr_control_deep)[i]
  
  if(mean(nmr_control_deep[i,])> mean(nmr_control_surface[i,])){
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

write_csv(nmr_final_res_filtered, 'OPLS_da/NMR_opls_results_control.csv')


# Di_venn input ----

Control_divenn_mets <- rbind(gcms_final_res_filtered,lcms_final_res_filtered, nmr_final_res_filtered)
write_csv(Control_divenn_mets, 'OPLS_da/control_diVenn_mets.csv')
