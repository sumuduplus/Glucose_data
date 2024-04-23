# Libraries ----

library(tidyverse)

# Load data ----

da_metabolites <- read_csv('Pathways/differential_temporal_mets.csv') %>% 
  mutate(enriched = ifelse(Log2FC > 0, 'Glucose abundant', 'Control abundant'))

da_sorted <- da_metabolites %>% 
  mutate(location = gsub(".*_","",comparison))

surface_paths_initial <- read_csv('Pathways/surface_kegg_initial.csv') %>% 
  mutate(KEGG_id = str_remove(KEGG_id, 'cpd:')) %>% 
  select(query, KEGG_id, KEGG_pathway) %>% 
  drop_na(KEGG_pathway) %>%
  separate_rows(KEGG_pathway, sep = ';') %>% 
  mutate(group = 'surface_initial')

surface_paths_final <- read_csv('Pathways/surface_kegg_final.csv') %>% 
  mutate(KEGG_id = str_remove(KEGG_id, 'cpd:')) %>% 
  select(query, KEGG_id, KEGG_pathway) %>% 
  drop_na(KEGG_pathway) %>%
  separate_rows(KEGG_pathway, sep = ';') %>% 
  mutate(group = 'surface_final')

deep_paths_initial <- read_csv('Pathways/deep_kegg_initial.csv') %>% 
  mutate(KEGG_id = str_remove(KEGG_id, 'cpd:')) %>% 
  select(query, KEGG_id, KEGG_pathway) %>% 
  drop_na(KEGG_pathway) %>%
  separate_rows(KEGG_pathway, sep = ';') %>% 
  mutate(group = 'deep_initial')

deep_paths_final <- read_csv('Pathways/deep_kegg_final.csv') %>% 
  mutate(KEGG_id = str_remove(KEGG_id, 'cpd:')) %>% 
  select(query, KEGG_id, KEGG_pathway) %>% 
  drop_na(KEGG_pathway) %>%
  separate_rows(KEGG_pathway, sep = ';') %>% 
  mutate(group = 'deep_final')


# Process data ----

da_initial_surface <- da_sorted %>% 
  filter(location == 'surface',
         time == 'initial') %>% 
  left_join(surface_paths_initial, by = c('FeatureID' = 'query')) %>% 
  drop_na(KEGG_pathway) %>%  
  group_by(time,location, enriched) %>% 
  count(KEGG_pathway)


da_final_surface <- da_sorted %>% 
  filter(location == 'surface',
         time == 'final') %>% 
  left_join(surface_paths_final, by = c('FeatureID' = 'query')) %>% 
  drop_na(KEGG_pathway) %>%  
  group_by(time,location, enriched)%>% 
  count(KEGG_pathway)

da_initial_deep <- da_sorted %>% 
  filter(location == 'deep',
         time == 'initial') %>% 
  left_join(deep_paths_initial, by = c('FeatureID' = 'query')) %>% 
  drop_na(KEGG_pathway) %>%  
  group_by(time,location, enriched) %>% 
  count(KEGG_pathway)

da_final_deep <- da_sorted %>% 
  filter(location == 'deep',
         time == 'final') %>% 
  left_join(deep_paths_final, by = c('FeatureID' = 'query')) %>% 
  drop_na(KEGG_pathway) %>%  
  group_by(time,location, enriched)%>% 
  count(KEGG_pathway)


pathway_counts <- rbind(da_initial_surface,
                        da_final_surface,
                        da_initial_deep,
                        da_final_deep) %>% 
  group_by(KEGG_pathway, location, enriched, time) %>% 
  summarise(n = sum(n)) %>% 
  unite(treatment_enrich, location, enriched,  time, sep = '_') %>% 
  pivot_wider(names_from = 'treatment_enrich', values_from = 'n', values_fill = 0)

write_csv(pathway_counts, 'Pathways/all_pathways_compared.csv')

selected_pathways <- read_csv('Pathways/all_pathways_compared_filtered.csv')

selected_pathways_plots <- selected_pathways %>%
  ungroup() %>% 
  pull(KEGG_pathway)


da_pathways_plot <- rbind(da_initial_surface,
                          da_final_surface,
                          da_initial_deep,
                          da_final_deep) %>%
  #separate(comparison, into = c('dataset', 'treatment'), sep = '_') %>% 
  filter(KEGG_pathway %in% selected_pathways_plots) %>% 
  group_by(KEGG_pathway, time, location, enriched) %>% 
  summarise(n = sum(n)) %>% 
  mutate(counts = ifelse(enriched == 'Control abundant', -n, n),
         time = factor(time, levels = c('initial', 'final')),
         location = factor(location, levels = c('surface', 'deep'))) %>% 
  ggplot() +
  geom_col(aes(x = counts,
               y = KEGG_pathway,
               fill = enriched),
           color = 'black') +
  facet_grid(rows = vars(location), cols = vars(time)) +
  scale_x_continuous(labels = abs) +
  theme_bw()

da_pathways_plot  

ggsave('da_pathways_plot.svg', dpi = 300, height = 14, width = 10)

# Extracting metabolites per pathway ----

metabolites_group <- da_sorted %>% 
  unite(group, location, time, enriched, sep = '_') %>% 
  select(-comparison, -Log2FC) 

pathways_group <- rbind(surface_paths_initial,
                        surface_paths_final,
                        deep_paths_initial,
                        deep_paths_final) %>% 
  select(-group) %>% 
  distinct() %>% 
  inner_join(metabolites_group, by = c('query' = 'FeatureID'))

# Function to check ids ----
extract_ids <- function(path){
  
  df <- pathways_group %>% 
    mutate(presence = TRUE) %>% 
    arrange(group) %>% 
    pivot_wider(names_from = group, values_from = 'presence', 
                values_fill = FALSE) %>% 
    filter(KEGG_pathway == path)
    
}

# Individual pathway extraction

aromatic_compounds_degradation <- extract_ids('Degradation of aromatic compounds')
write_csv(aromatic_compounds_degradation, 'Pathways/aromatic_compounds_degradation_hits.csv')

Carbon_metabolism <- extract_ids('Carbon metabolism')
write_csv(Carbon_metabolism, 'Pathways/Carbon_metabolism_hits.csv')

Tyrosine_metabolism <- extract_ids('Tyrosine metabolism')
write_csv(Tyrosine_metabolism, 'Pathways/Tyrosine_metabolism_hits.csv')
