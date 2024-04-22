# Libraries ----

library(tidyverse)
library(readxl)
library(ggplot2)
library(factoextra)
library(vegan)
library(dplyr)
library(ape)
library(ade4)
library(made4)
library(patchwork)
library(reshape2)

# Load data ----


lcms_metadata <- read_xlsx("Raw_data/Metadata_all.xlsx", sheet = 2) %>% 
  unite(uniqueID, Location, Treatment, Days,sep = "_", remove = FALSE) %>% 
  mutate(Days = factor(Days, levels = c('7', '14', '28', '42', '56', '70')))

lcms_data <-read_csv("Preprocess_data/LCMS_normalized.csv") %>% 
  left_join(select(lcms_metadata, SampleID, uniqueID), by = 'SampleID') %>% 
  select(-SampleID)

gcms_metadata <- read_xlsx("Raw_data/Metadata_all.xlsx", sheet = 1) %>% 
  unite(uniqueID, Location, Treatment, Days,sep = "_", remove = FALSE) %>% 
  mutate(Days = factor(Days, levels = c('7', '14', '28', '42', '56', '70')))

gcms_data <- read_csv("Preprocess_data/GCMS_normalized.csv") %>% 
  left_join(select(gcms_metadata, SampleID, uniqueID), by = 'SampleID') %>% 
  select(-SampleID)

nmr_metadata <- read_xlsx("Raw_data/Metadata_all.xlsx", sheet = 3) %>% 
  unite(uniqueID, Location, Treatment, Days,sep = "_", remove = FALSE) %>% 
  mutate(Days = factor(Days, levels = c('7', '14', '28', '42', '56', '70')))

nmr_data <- read_csv("Preprocess_data/NMR_normalized.csv") %>% 
  left_join(select(nmr_metadata, SampleID, uniqueID), by = 'SampleID') %>% 
  select(-SampleID)

microbial_metadata <- read_xlsx("Raw_data/Metadata_all.xlsx", sheet = 4) %>% 
  unite(uniqueID, Location, Treatment, TimePoint,SIPFraction, sep = "_", remove = FALSE)%>%
  column_to_rownames(var = 'SampleID')

microbial_metadata$Days <- sub('day 0','day 7',microbial_metadata$Days)

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


# sample sorting ----

LCms_matrix <- lcms_data %>% 
  filter(grepl("7|70", lcms_data$uniqueID)) 

LCms_matrix_surface <- LCms_matrix%>% 
  filter(grepl("Surface", LCms_matrix$uniqueID))%>%
  column_to_rownames(var = 'uniqueID')
LCms_matrix_deep <- LCms_matrix%>% 
  filter(grepl("Deep", LCms_matrix$uniqueID))%>%
  column_to_rownames(var = 'uniqueID')

gcms_matrix<- gcms_data %>% 
  filter(grepl("7|70", gcms_data$uniqueID)) 

gcms_matrix_surface <- gcms_matrix%>% 
  filter(grepl("Surface", gcms_matrix$uniqueID))%>%
  column_to_rownames(var = 'uniqueID')
gcms_matrix_deep <- gcms_matrix%>% 
  filter(grepl("Deep", gcms_matrix$uniqueID))%>%
  column_to_rownames(var = 'uniqueID')

nmr_matrix <- nmr_data %>% 
  filter(grepl("7|70", nmr_data$uniqueID))

nmr_matrix_surface <- nmr_matrix%>% 
  filter(grepl("Surface", nmr_matrix$uniqueID))%>%
  column_to_rownames(var = 'uniqueID')
nmr_matrix_deep <- nmr_matrix%>% 
  filter(grepl("Deep", nmr_matrix$uniqueID))%>%
  column_to_rownames(var = 'uniqueID')

Microbial_data_values <-  subset(microbial_data, select = -c(Kingdom,Phylum,Class,Order,Family,Genus,Species))%>%
  column_to_rownames(var = 'FeatureID')%>% 
  select(contains('S'))%>% 
  select(-contains('13C'))

Microbial_data_values_t <- t(Microbial_data_values)

microbial_metadata_updated <- subset(microbial_metadata, select = -c(Depth,Location,TimePoint,SIPFraction,density,Days,Treatment))

All_unfractionated <-  merge(Microbial_data_values_t, microbial_metadata_updated, by=0, all=F)%>%
  remove_rownames 

Surface_unfractionated <- All_unfractionated %>% 
  filter(grepl("Surface", All_unfractionated$uniqueID))%>%
  column_to_rownames(var="uniqueID") %>% 
  select(-contains('Row.names'))

deep_unfractionated <- All_unfractionated %>% 
  filter(grepl("Deep", All_unfractionated$uniqueID))%>%
  column_to_rownames(var="uniqueID") %>% 
  select(-contains('Row.names'))


# Surface samples ----

lcms_dist_surface <- vegdist(LCms_matrix_surface, method = 'bray')
PCoA_LCMS_surface <- pcoa(lcms_dist_surface, correction="cailliez", rn=NULL)

gcms_dist_surface <- vegdist(gcms_matrix_surface, method = 'bray')
PCoA_gcMS_surface <- pcoa(gcms_dist_surface, correction="cailliez", rn=NULL)

nmr_dist_surface <- vegdist(nmr_matrix_surface, method = 'bray')
PCoA_nmr_surface <- pcoa(nmr_dist_surface, correction="cailliez", rn=NULL)

pca_surface <- prcomp(Surface_unfractionated, center = TRUE)
eigen_surface <- get_eigenvalue(pca_surface)


# CIA surface samples ----

LCMS_surface <- PCoA_LCMS_surface$vectors
GCMS_surface <- PCoA_gcMS_surface$vectors
NMR_surface <- PCoA_nmr_surface$vectors

ASV_components_S <- pca_surface$x
rownames(ASV_components_S) = gsub("_Tinitial_S", "_7", rownames(ASV_components_S))
rownames(ASV_components_S) = gsub("_Tfinal_S", "_70", rownames(ASV_components_S))

ASV_surface <- ASV_components_S[rownames(GCMS_surface), ]

CIA_LC_Surface <- cia(t(ASV_surface), t(LCMS_surface))
attach(CIA_LC_Surface )
summary(CIA_LC_Surface)
summary(CIA_LC_Surface$coinertia)

plot(CIA_LC_Surface)
plot(CIA_LC_Surface, classvec=NCI60$classes[,2], clab=0, cpoint=3)

CIA_GC_Surface <- cia(t(ASV_surface), t(GCMS_surface))
attach(CIA_GC_Surface)
summary(CIA_GC_Surface )
summary(CIA_GC_Surface$coinertia)

plot(CIA_GC_Surface)
plot(CIA_GC_Surface, classvec=NCI60$classes[,2], clab=0, cpoint=3)

CIA_NMR_Surface <- cia(t(ASV_surface), t(NMR_surface))
attach(CIA_NMR_Surface )
summary(CIA_NMR_Surface)
summary(CIA_NMR_Surface$coinertia)

plot(CIA_NMR_Surface)
plot(CIA_NMR_Surface, classvec=NCI60$classes[,2], clab=0, cpoint=3)

# CIA deep samples ----

lcms_dist_deep <- vegdist(LCms_matrix_deep, method = 'bray')
PCoA_LCMS_deep <- pcoa(lcms_dist_deep, correction="cailliez", rn=NULL)

gcms_dist_deep <- vegdist(gcms_matrix_deep, method = 'bray')
PCoA_gcMS_deep <- pcoa(gcms_dist_deep, correction="cailliez", rn=NULL)

nmr_dist_deep <- vegdist(nmr_matrix_deep, method = 'bray')
PCoA_nmr_deep <- pcoa(nmr_dist_deep, correction="cailliez", rn=NULL)

pca_deep <- prcomp(deep_unfractionated, center = TRUE)
eigen_deep <- get_eigenvalue(pca_deep)

LCMS_deep <- PCoA_LCMS_deep$vectors
GCMS_deep <- PCoA_gcMS_deep$vectors
NMR_deep <- PCoA_nmr_deep$vectors

ASV_components_D <- pca_deep$x
rownames(ASV_components_D) = gsub("_Tinitial_S", "_7", rownames(ASV_components_D))
rownames(ASV_components_D) = gsub("_Tfinal_S", "_70", rownames(ASV_components_D))

ASV_deep <- ASV_components_D[rownames(GCMS_deep), ]

CIA_LC_deep <- cia(t(ASV_deep), t(LCMS_deep))
attach(CIA_LC_deep )
summary(CIA_LC_deep)
summary(CIA_LC_deep$coinertia)

plot(CIA_LC_deep)
plot(CIA_LC_deep, classvec=NCI60$classes[,2], clab=0, cpoint=3)

CIA_GC_deep <- cia(t(ASV_deep), t(GCMS_deep),cia.nf=2 )
attach(CIA_GC_deep)
summary(CIA_GC_deep )
summary(CIA_GC_deep$coinertia)

plot(CIA_GC_deep)
plot(CIA_GC_deep, classvec=NCI60$classes[,2], clab=0, cpoint=3)

CIA_NMR_deep <- cia(t(ASV_deep), t(NMR_deep))
attach(CIA_NMR_deep )
summary(CIA_NMR_deep)
summary(CIA_NMR_deep$coinertia)

plot(CIA_NMR_deep)
plot(CIA_NMR_deep, classvec=NCI60$classes[,2], clab=0, cpoint=3)


# Monte Carlo simulations ----

ASV_surface <- as.data.frame(ASV_surface)
LCMS_surface <- as.data.frame(LCMS_surface)
GCMS_surface <- as.data.frame(GCMS_surface)
NMR_surface <- as.data.frame(NMR_surface)


RV_Surface_LC <- RV.rtest (ASV_surface, LCMS_surface, nrepet = 999) 
RV_Surface_GC <- RV.rtest (ASV_surface, GCMS_surface, nrepet = 999) 
RV_Surface_NMR <- RV.rtest (ASV_surface, NMR_surface, nrepet = 999) 

ASV_deep <- as.data.frame(ASV_deep)
LCMS_deep <- as.data.frame(LCMS_deep)
GCMS_deep <- as.data.frame(GCMS_deep)
NMR_deep <- as.data.frame(NMR_deep)

RV_deep_LC <- RV.rtest (ASV_deep, LCMS_deep, nrepet = 999) 
RV_deep_GC <- RV.rtest (ASV_deep, GCMS_deep, nrepet = 999) 
RV_deep_NMR <- RV.rtest (ASV_deep, NMR_deep, nrepet = 999) 


# plot cooordinates for biplots ----

# surface plots ----

# Surface LC vs ASV
Surface_LC_ASV <- CIA_LC_Surface$coinertia$mX %>% 
  mutate(Origin=c("ASV", "ASV", "ASV","ASV")) %>%
  rownames_to_column(var="Combined") %>% 
  unite(col=UniqueID,c("Origin", "Combined"), sep="", remove = FALSE)

Surface_LC_LC <- CIA_LC_Surface$coinertia$mY%>%
  mutate(Origin=c("LC", "LC", "LC","LC"))%>%
  rownames_to_column(var="Combined") %>% 
  unite(col=UniqueID,c("Origin", "Combined"), sep="", remove = FALSE)

Surface_CIA_LC <- rbind(Surface_LC_ASV,Surface_LC_LC ) %>% 
  mutate(Days = str_replace(Combined,"Surface_Control_|Surface_Glucose_", "D")) %>% 
  mutate(Treatment = str_replace(Combined,"Surface_Control_|Surface_Glucose_", "D"))


CIA_surface_LC <- Surface_CIA_LC %>% 
  ggplot () +
  geom_line(aes(x=NorS1, 
                y=NorS2,
                group = Combined)) +
  xlim(-2.0, 2.0) +
  ylim(-1.6,1.6) +
  geom_point(aes(x=NorS1, 
                 y=NorS2,
                 shape = Days,
                 color = Origin),
             size = 5) +
  labs(title = 'CIA: Microbial PCA and LC PCoA') +
  scale_color_manual(values = c(ASV = '#B026FF', LC = '#40E0D0'))

CIA_surface_LC

# Surface GC vs ASV

Surface_GC_ASV <- CIA_GC_Surface$coinertia$mX %>% 
  mutate(Origin=c("ASV", "ASV", "ASV","ASV")) %>%
  rownames_to_column(var="Combined") %>% 
  unite(col=UniqueID,c("Origin", "Combined"), sep="", remove = FALSE)

Surface_GC_GC <- CIA_GC_Surface$coinertia$mY%>%
  mutate(Origin=c("GC", "GC", "GC","GC"))%>%
  rownames_to_column(var="Combined") %>% 
  unite(col=UniqueID,c("Origin", "Combined"), sep="", remove = FALSE)

Surface_CIA_GC <- rbind(Surface_GC_ASV,Surface_GC_GC ) %>% 
  mutate(Days = str_replace(Combined,"Surface_Control_|Surface_Glucose_", "D")) %>% 
  mutate(Treatment = str_replace(Combined,"Surface_Control_|Surface_Glucose_", "D"))


CIA_surface_GC <- Surface_CIA_GC %>% 
  ggplot () +
  geom_line(aes(x=NorS1, 
                y=NorS2,
                group = Combined)) +
  xlim(-2.0, 2.0) +
  ylim(-1.6,1.6) +
  geom_point(aes(x=NorS1, 
                 y=NorS2,
                 shape = Days,
                 color = Origin),
             size = 5)  +
  labs(title = 'CIA: Microbial PCA and GC PCoA') +
  scale_color_manual(values = c(ASV = '#B026FF', GC = '#007FFF'))

CIA_surface_GC


# Surface NMR vs ASV

Surface_NMR_ASV <- CIA_NMR_Surface$coinertia$mX %>% 
  mutate(Origin=c("ASV", "ASV", "ASV","ASV")) %>%
  rownames_to_column(var="Combined") %>% 
  unite(col=UniqueID,c("Origin", "Combined"), sep="", remove = FALSE)

Surface_NMR_NMR <- CIA_NMR_Surface$coinertia$mY%>%
  mutate(Origin=c("NMR", "NMR", "NMR","NMR"))%>%
  rownames_to_column(var="Combined") %>% 
  unite(col=UniqueID,c("Origin", "Combined"), sep="", remove = FALSE)

Surface_CIA_NMR <- rbind(Surface_NMR_ASV,Surface_NMR_NMR ) %>% 
  mutate(Days = str_replace(Combined,"Surface_Control_|Surface_Glucose_", "D")) %>% 
  mutate(Treatment = str_replace(Combined,"Surface_Control_|Surface_Glucose_", "D"))


CIA_surface_NMR <- Surface_CIA_NMR %>% 
  ggplot () +
  geom_line(aes(x=NorS1, 
                y=NorS2,
                group = Combined)) +
  xlim(-2.0, 2.0) +
  ylim(-1.6,1.6) +
  geom_point(aes(x=NorS1, 
                 y=NorS2,
                 shape = Days,
                 color = Origin),
             size = 5) +
  labs(title = 'CIA: Microbial PCA and NMR PCoA') +
  scale_color_manual(values = c(ASV = '#B026FF', NMR = '#00308F'))

CIA_surface_NMR

# deep plots ----

# deep LC vs ASV
deep_LC_ASV <- CIA_LC_deep$coinertia$mX %>% 
  mutate(Origin=c("ASV", "ASV", "ASV","ASV")) %>%
  rownames_to_column(var="Combined") %>% 
  unite(col=UniqueID,c("Origin", "Combined"), sep="", remove = FALSE)

deep_LC_LC <- CIA_LC_deep$coinertia$mY%>%
  mutate(Origin=c("LC", "LC", "LC","LC"))%>%
  rownames_to_column(var="Combined") %>% 
  unite(col=UniqueID,c("Origin", "Combined"), sep="", remove = FALSE)

deep_CIA_LC <- rbind(deep_LC_ASV,deep_LC_LC ) %>% 
  mutate(Days = str_replace(Combined,"Deep_Control_|Deep_Glucose_", "D")) %>% 
  mutate(Treatment = str_replace(Combined,"Deep_Control_|Deep_Glucose_", "D"))


CIA_deep_LC <- deep_CIA_LC %>% 
  ggplot () +
  geom_line(aes(x=NorS1, 
                y=NorS2,
                group = Combined)) +
  xlim(-2.0, 2.0) +
  ylim(-1.6,1.6) +
  geom_point(aes(x=NorS1, 
                 y=NorS2,
                 shape = Days,
                 color = Origin),
             size = 5) +
  labs(title = 'CIA: Microbial PCA and LC PCoA') +
  scale_color_manual(values = c(ASV = '#B026FF', LC = '#40E0D0'))

CIA_deep_LC

# deep GC vs ASV

deep_GC_ASV <- CIA_GC_deep$coinertia$mX %>% 
  mutate(Origin=c("ASV", "ASV", "ASV","ASV")) %>%
  rownames_to_column(var="Combined") %>% 
  unite(col=UniqueID,c("Origin", "Combined"), sep="", remove = FALSE)

deep_GC_GC <- CIA_GC_deep$coinertia$mY%>%
  mutate(Origin=c("GC", "GC", "GC","GC"))%>%
  rownames_to_column(var="Combined") %>% 
  unite(col=UniqueID,c("Origin", "Combined"), sep="", remove = FALSE)

deep_CIA_GC <- rbind(deep_GC_ASV,deep_GC_GC ) %>% 
  mutate(Days = str_replace(Combined,"Deep_Control_|Deep_Glucose_", "D")) %>% 
  mutate(Treatment = str_replace(Combined,"Deep_Control_|Deep_Glucose_", "D"))


CIA_deep_GC <- deep_CIA_GC %>% 
  ggplot () +
  geom_line(aes(x=NorS1, 
                y=NorS2,
                group = Combined)) +
  xlim(-2.0, 2.0) +
  ylim(-1.6,1.6) +
  geom_point(aes(x=NorS1, 
                 y=NorS2,
                 shape = Days,
                 color = Origin),
             size = 5)  +
  labs(title = 'CIA: Microbial PCA and GC PCoA') +
  scale_color_manual(values = c(ASV = '#B026FF', GC = '#007FFF'))

CIA_deep_GC


# deep NMR vs ASV

deep_NMR_ASV <- CIA_NMR_deep$coinertia$mX %>% 
  mutate(Origin=c("ASV", "ASV", "ASV","ASV")) %>%
  rownames_to_column(var="Combined") %>% 
  unite(col=UniqueID,c("Origin", "Combined"), sep="", remove = FALSE)

deep_NMR_NMR <- CIA_NMR_deep$coinertia$mY%>%
  mutate(Origin=c("NMR", "NMR", "NMR","NMR"))%>%
  rownames_to_column(var="Combined") %>% 
  unite(col=UniqueID,c("Origin", "Combined"), sep="", remove = FALSE)

deep_CIA_NMR <- rbind(deep_NMR_ASV,deep_NMR_NMR ) %>% 
  mutate(Days = str_replace(Combined,"Deep_Control_|Deep_Glucose_", "D")) %>% 
  mutate(Treatment = str_replace(Combined,"deep_Control_|deep_Glucose_", "D"))


CIA_deep_NMR <- deep_CIA_NMR %>% 
  ggplot () +
  geom_line(aes(x=NorS1, 
                y=NorS2,
                group = Combined)) +
  xlim(-2.0, 2.0) +
  ylim(-1.6,1.6) +
  geom_point(aes(x=NorS1, 
                 y=NorS2,
                 shape = Days,
                 color = Origin),
             size = 5) +
  labs(title = 'CIA: Microbial PCA and NMR PCoA') +
  scale_color_manual(values = c(ASV = '#B026FF', NMR = '#00308F'))

CIA_deep_NMR



# Combined plots ----

Surface_CIA_plots <- (CIA_surface_LC | CIA_surface_GC | CIA_surface_NMR ) +
  plot_layout(guides = 'collect')

Surface_CIA_plots

ggsave('Surface_CIA_plots.svg', dpi = 300, height = 5, width = 8)  


Deep_CIA_plots <- (CIA_deep_LC | CIA_deep_GC | CIA_deep_NMR ) +
  plot_layout(guides = 'collect')

Deep_CIA_plots

ggsave('Deep_CIA_plots.svg', dpi = 300, height = 5, width = 8)  


