library(readxl)
library(ropls)
library(tidyverse)
library(dplyr)

# input data ----

microbial_data <- read_csv("16S/otu_chord_transformed_tax.csv")%>%
  column_to_rownames(var="FeatureID") %>%
  select(-contains('frs'))%>% 
  select(-contains('Una'))
  

Microbial_data_values <- subset(microbial_data, select = -c(Kingdom,Phylum,Class,Order,Family,Genus,Species))%>%
  t()

metadata <- read_xlsx('Raw_data/Metadata_all.xlsx', sheet = 4) %>% 
  column_to_rownames(var = 'SampleID')

All_combined <- merge(Microbial_data_values,metadata, by=0, all.x=TRUE)%>% 
  column_to_rownames(var = 'Row.names')


microbial_tax <- microbial_data %>% 
  select(Kingdom, Phylum, Class, Order, Family, Genus, Species)


# DEEP only OPLS-da

Deep_microbes <- All_combined %>%
  filter(Location == "Deep")

Deep_calss <- as.factor (Deep_microbes$Treatment) 
Deep_op_data <- Deep_microbes[, !colnames(Deep_microbes) %in% colnames(metadata)]

Deep_OPLS <- opls (Deep_op_data, Deep_calss, predI = 1, crossvalI = 47, orthoI = NA, permI = 1000) 

Deep_da_VIPS <- as.data.frame(getVipVn(Deep_OPLS,orthoL = FALSE)) %>% 
  rename(VIP = 1)%>%
  filter(VIP > 1)

Deep_da_ASVs <-  merge(Deep_da_VIPS,microbial_tax, by=0, all.x=TRUE)
  
write_csv(Deep_da_ASVs, 'OPLS_da/Deep_differential_members_isotopic_frac.csv')

# Surface only OPLS-da

Surface_microbes <- All_combined %>%
  filter(Location == "Surface")

Surface_calss <- as.factor (Surface_microbes$Treatment) 
Surface_op_data <- Surface_microbes[, !colnames(Surface_microbes) %in% colnames(metadata)]

Surface_OPLS <- opls (Surface_op_data, Surface_calss, predI = 1, crossvalI = 47, orthoI = NA, permI = 1000) 

Surface_da_VIPS <- as.data.frame(getVipVn(Surface_OPLS,orthoL = FALSE)) %>% 
  rename(VIP = 1)%>%
  filter(VIP > 1)

Surface_da_ASVs <-  merge(Surface_da_VIPS,microbial_tax, by=0, all.x=TRUE)

write_csv(Surface_da_ASVs, 'OPLS_da/Surface_differential_members_isotopic_frac.csv')

