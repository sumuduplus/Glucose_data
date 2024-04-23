# Libraries ----

library(tidyverse)
library(readxl)
library(KEGGREST)

# load data -----
gcms_metadata <- read_xlsx("Raw_data/Metadata_all.xlsx", sheet = 1) %>% 
  unite(uniqueID, Location, Treatment, Days,sep = "_", remove = FALSE) %>% 
  mutate(Days = factor(Days, levels = c('7', '14', '28', '42', '56', '70')))

nmr_metadata <- read_xlsx("Raw_data/Metadata_all.xlsx", sheet = 3) %>% 
  unite(uniqueID, Location, Treatment, Days,sep = "_", remove = FALSE) %>% 
  mutate(Days = factor(Days, levels = c('7', '14', '28', '42', '56', '70')))

gcms_data <- read_csv("Preprocess_data/GCMS_normalized.csv") %>% 
  left_join(select(gcms_metadata, SampleID, uniqueID), by = 'SampleID') %>% 
  select(-SampleID)

nmr_data <- read_csv("Preprocess_data/NMR_normalized.csv") %>% 
  left_join(select(nmr_metadata, SampleID, uniqueID), by = 'SampleID') %>% 
  select(-SampleID)

# Filtering data ----

gcms_surface <- gcms_data %>% 
  filter(str_detect(uniqueID, 'Surface_')) %>% 
  column_to_rownames(var = 'uniqueID') %>% 
  t()

gcms_deep <- gcms_data %>% 
  filter(str_detect(uniqueID, 'Deep_'))  %>% 
  column_to_rownames(var = 'uniqueID') %>% 
  t()

nmr_surface <- nmr_data %>% 
  filter(str_detect(uniqueID, 'Surface_')) %>% 
  column_to_rownames(var = 'uniqueID') %>% 
  t()

nmr_deep <- nmr_data %>% 
  filter(str_detect(uniqueID, 'Deep_'))  %>% 
  column_to_rownames(var = 'uniqueID') %>% 
  t()


# Function for differential abundance  control vs glucose----

calculate_da <- function(df){
  df_untransformed <- df ^ 3
  
  control_samples <- colnames(df_untransformed)[str_detect(colnames(df_untransformed), 'Control')]
  glucose_samples <- colnames(df_untransformed)[str_detect(colnames(df_untransformed), 'Glucose')]
  
  control_initial <- control_samples[str_detect(control_samples, '7$|14|28')]
  control_final <- control_samples[str_detect(control_samples, '42|56|70')]
  
  glucose_initial <- glucose_samples[str_detect(glucose_samples, '7$|14|28')]
  glucose_final <- glucose_samples[str_detect(glucose_samples, '42|56|70')]
  
  
  df_control_initial <- df_untransformed[,control_initial]
  df_control_final <- df_untransformed[,control_final]
  df_glucose_initial <- df_untransformed[,glucose_initial]
  df_glucose_final <- df_untransformed[,glucose_final]
  
  
  df_results_initial <- tibble(FeatureID = NA,
                               Log2FC = NA,
                               time = 'initial',
                               .rows = nrow(df_control_initial))
  
  df_results_final <- tibble(FeatureID = NA,
                             Log2FC = NA,
                             time = 'final',
                             .rows = nrow(df_control_final))
  
  for(i in 1:nrow(df_control_initial)){
    
    mean_control_i <- mean(df_control_initial[i,])
    mean_glucose_i <- mean(df_glucose_initial[i,])
    
    log2FC_i <- log2(mean_glucose_i/mean_control_i)
    
    df_results_initial$FeatureID[i] <- rownames(df_control_initial)[i]
    df_results_initial$Log2FC[i] <- log2FC_i
    
    mean_control_f <- mean(df_control_final[i,])
    mean_glucose_f <- mean(df_glucose_final[i,])
    
    log2FC_f <- log2(mean_control_f/mean_glucose_f)
    
    df_results_final$FeatureID[i] <- rownames(df_control_final)[i]
    df_results_final$Log2FC[i] <- log2FC_f
    
  }
  
  df_results <- rbind(df_results_initial, df_results_final)
  
  return(df_results)
  
}

# Function for mapping into KEGG

get_kegg_by_name <- function(compound_names){
  
  kegg_info <- map(compound_names, function(name){
    
    print(paste0('Working on ', name))
    
    cpd_ids <- names(keggFind('compound', name))
    
    if(!is.null(cpd_ids)){
      searches <- map(cpd_ids, function(k_id){
        cpd_info <- keggGet(k_id)
        
        df <- tibble(
          query = name,
          KEGG_id = k_id,
          KEGG_name = paste(cpd_info[[1]]$NAME, collapse = ''),
          KEGG_formula = ifelse(!is.null(cpd_info[[1]]$FORMULA), paste(cpd_info[[1]]$FORMULA, collapse = ';'), NA),
          KEGG_pathway = ifelse(!is.null(cpd_info[[1]]$PATHWAY), paste(cpd_info[[1]]$PATHWAY, collapse = ';'), NA),
          KEGG_module = ifelse(!is.null(cpd_info[[1]]$MODULE), paste(cpd_info[[1]]$MODULE, collapse = ';'), NA),
          KEGG_brite = ifelse(!is.null(cpd_info[[1]]$BRITE), paste(cpd_info[[1]]$BRITE, collapse = ';'), NA),
          KEGG_enzyme = ifelse(!is.null(cpd_info[[1]]$ENZYME), paste(cpd_info[[1]]$ENZYME, collapse = ';'), NA),
          KEGG_reaction = ifelse(!is.null(cpd_info[[1]]$REACTION), paste(cpd_info[[1]]$REACTION, collapse = ';'), NA)
        )
      })
      
      cpd_df <- reduce(searches, rbind)
      
      Sys.sleep(30)
      
    } else {
      
      cpd_df <- tibble(
        query = name,
        KEGG_id = NA,
        KEGG_name = NA,
        KEGG_formula = NA,
        KEGG_pathway = NA,
        KEGG_module = NA,
        KEGG_brite = NA,
        KEGG_enzyme = NA,
        KEGG_reaction = NA
      )
    }
    
    return(cpd_df)
  })
  
  final <- reduce(kegg_info, rbind)
}

# Calculating Diff Abundance across timepoints ----

#gcms
gcms_surface_da <- calculate_da(gcms_surface) %>% 
  filter(abs(Log2FC) > 1) %>% 
  mutate(comparison = 'gcms_surface')

gcms_deep_da <- calculate_da(gcms_deep)  %>% 
  filter(abs(Log2FC) > 1) %>% 
  mutate(comparison = 'gcms_deep')

gcms_deep_da_initial <- gcms_deep_da %>% 
  filter(time == 'initial')

gcms_deep_da_final <- gcms_deep_da %>% 
  filter(time == 'final')

gcms_surface_da_initial <- gcms_surface_da %>% 
  filter(time == 'initial')

gcms_surface_da_final <- gcms_surface_da %>% 
  filter(time == 'final')

#nmr

nmr_surface_da <- calculate_da(nmr_surface) %>% 
  filter(abs(Log2FC) > 1) %>% 
  mutate(comparison = 'NMR_surface')

nmr_deep_da <- calculate_da(nmr_deep)  %>% 
  filter(abs(Log2FC) > 1) %>% 
  mutate(comparison = 'NMR_deep')

nmr_deep_da_initial <- nmr_deep_da %>% 
  filter(time == 'initial')

nmr_deep_da_final <- nmr_deep_da %>% 
  filter(time == 'final')

nmr_surface_da_initial <- nmr_surface_da %>% 
  filter(time == 'initial')

nmr_surface_da_final <- nmr_surface_da %>% 
  filter(time == 'final')

#binding

da_all <- rbind(gcms_surface_da_initial,
                nmr_surface_da_initial,
                gcms_surface_da_final,
                nmr_surface_da_final,
                gcms_deep_da_initial,
                gcms_deep_da_final,
                nmr_deep_da_initial,
                nmr_deep_da_final
                )

# Saving differential abundant metabolites ----

write_csv(da_all, 'Pathways/differential_temporal_mets.csv')


# Mapping to KEGG ----

all_differ_mets <- unique(da_all$FeatureID)
all_with_kegg <- get_kegg_by_name(all_differ_mets)

# modify and upload filtered differential abundant metabolites
write_csv(all_with_kegg, 'Pathways/all_kegg_hits.csv')


all_with_kegg_updated <- read_csv("Pathways/all_kegg_hits_updated.csv")


# filtering

da_surface_initial <- rbind(gcms_surface_da_initial,
                nmr_surface_da_initial)

da_surface_final <- rbind(gcms_surface_da_final,
                            nmr_surface_da_final)

da_deep_initial <- rbind(gcms_deep_da_initial,
                            nmr_deep_da_initial)

da_deep_final <- rbind(gcms_deep_da_final,
                            nmr_deep_da_final)


surface_kegg_initial <- all_with_kegg_updated %>% 
  filter(query %in% da_surface_initial$FeatureID)

surface_kegg_final <- all_with_kegg_updated %>% 
  filter(query %in% da_surface_final$FeatureID)

deep_kegg_initial <- all_with_kegg_updated %>% 
  filter(query %in% da_deep_initial$FeatureID)

deep_kegg_final <- all_with_kegg_updated %>% 
  filter(query %in% da_deep_final$FeatureID)

write_csv(surface_kegg_initial, 'Pathways/surface_kegg_initial.csv')
write_csv(surface_kegg_final, 'Pathways/surface_kegg_final.csv')
write_csv(deep_kegg_initial, 'Pathways/deep_kegg_initial.csv')
write_csv(deep_kegg_final, 'Pathways/deep_kegg_final.csv')
