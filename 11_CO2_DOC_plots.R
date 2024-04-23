### CO2 and DOC Data 


# Load Libraries & Functions ----------------------------------------------------------
library(tidyverse)
library(ggplot2)
library(purrr)
library(magrittr)
library(cowplot)
library(gridExtra)



# Function to perform unpaired t-test given mean, standard deviation, and sample size
perform_t_test <- function(mean1, mean2, sd1, sd2, n1, n2) {
  se <- sqrt((sd1^2 / n1) + (sd2^2 / n2))  # Calculate standard error
  t_value <- (mean1 - mean2) / se  # Calculate t-value
  p_value <- 2 * pt(-abs(t_value), df = n1 + n2 - 2)  # Calculate two-tailed p-value
  return(p_value)
}


# Load Data ---------------------------------------------------------------
data <- read.csv("Raw_data/CO2_data.csv")


# 1. CO2 Significance  ----------------------------------------------------
co2_conc <- data %>%
  select(Treatment, Depth, Timepoints, CO2_produced, CO2_produced_sd, Obs_numb)%>%
  filter(Timepoints == c(7, 70))


## 1.1 Statistics --------------------------------------------------------------

# Perform comparisons for all combinations of treatments, depths, and time points
co2_stats_test <- data.frame()

for (treatment1 in unique(co2_conc$Treatment)) {
  for (treatment2 in unique(co2_conc$Treatment)) {
    for (depth1 in unique(co2_conc$Depth)) {
      for (depth2 in unique(co2_conc$Depth)) {
        for (timepoint1 in unique(co2_conc$Timepoint)) {
          for (timepoint2 in unique(co2_conc$Timepoint)) {
            # Ensure treatment1, depth1, and timepoint1 are less than or equal to treatment2, depth2, and timepoint2
            if (treatment1 <= treatment2 & depth1 <= depth2 & timepoint1 <= timepoint2) {
              p_value <- perform_t_test(
                mean1 = mean(co2_conc$CO2_produced[co2_conc$Treatment == treatment1 & co2_conc$Depth == depth1 & co2_conc$Timepoint == timepoint1]),
                mean2 = mean(co2_conc$CO2_produced[co2_conc$Treatment == treatment2 & co2_conc$Depth == depth2 & co2_conc$Timepoint == timepoint2]),
                sd1 = mean(co2_conc$CO2_produced_sd[co2_conc$Treatment == treatment1 & co2_conc$Depth == depth1 & co2_conc$Timepoint == timepoint1]),
                sd2 = mean(co2_conc$CO2_produced_sd[co2_conc$Treatment == treatment2 & co2_conc$Depth == depth2 & co2_conc$Timepoint == timepoint2]),
                n1 = sum(co2_conc$Obs_numb[co2_conc$Treatment == treatment1 & co2_conc$Depth == depth1 & co2_conc$Timepoint == timepoint1]),
                n2 = sum(co2_conc$Obs_numb[co2_conc$Treatment == treatment2 & co2_conc$Depth == depth2 & co2_conc$Timepoint == timepoint2])
              )
              comparison <- data.frame(
                treatment1 = treatment1,
                depth1 = depth1,
                timepoint1 = timepoint1,
                treatment2 = treatment2,
                depth2 = depth2,
                timepoint2 = timepoint2,
                p_value = p_value
              )
              co2_stats_test <- rbind(co2_stats_test, comparison)
              
            }
          }
        }
      }
    }
  }
}


# Print summary table
print(co2_stats_test)

#write.csv(co2_stats_test, "co2_stats_test.csv", row.names = F)


## 1.2. Plotting -----------------------------------------------------------
CO2_plot_data <- co2_stats_test %>% 
  filter(treatment1 == 'Control',
         treatment2 == 'Glucose',
         timepoint1 == timepoint2,
         depth1 == depth2) %>% 
  mutate(signif = case_when(p_value < 0.05 ~ '*',
                            p_value < 0.01 ~ '**',
                            p_value < 0.001 ~ '***',
                            TRUE ~ 'NS'),
         y_pos = 14500) %>% 
  rename(Depth = depth1,
         Timepoints = timepoint1)

Depth_order <- c("Surface", "Deep")
co2_conc$Depth <- factor(co2_conc$Depth, levels = Depth_order)
CO2_plot_data$Depth <- factor(CO2_plot_data$Depth, levels = Depth_order)

# Plotting
CO2_plot <- ggplot(co2_conc, aes(x=Treatment, y=CO2_produced, fill=Treatment)) +
  geom_bar(stat="identity", position="dodge") +
  geom_errorbar(aes(ymin=CO2_produced-CO2_produced_sd, ymax=CO2_produced+CO2_produced_sd), position="dodge", width=0.2) +

  facet_grid(~Depth * Timepoints) +
  labs(title="Cumulative CO"[2]~"Concentration", y="CO"[2]~"Concentration in ppm") +
  
  theme_classic() +
  geom_text(data = CO2_plot_data,
            aes(x = 1.5,
                y = y_pos,
                label = signif),
            size = 3,
            inherit.aes = FALSE) +
  geom_segment(x = 1,
               xend = 2,
               y = 13750,
               yend = 13750) +
  geom_segment(x = 1,
               xend = 1,
               y = 13750,
               yend = 13500) +
  geom_segment(x = 2,
               xend = 2,
               y = 13750,
               yend = 13500)+
  theme(plot.title = element_text(hjust = 0.5, size = 11),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 9),
        axis.ticks.x = element_blank())


CO2_plot

ggsave("CO2_plot.svg", CO2_plot, width = 6, height = 4, device = "svg")


# 2. DOC Significance -----------------------------------------------------

DOC_conc <- data %>%
  select(Treatment, Depth, Timepoints, DOC_conc, DOC_conc_sd, Obs_numb)%>%
  filter(Timepoints == c(7, 70))


## 2.1. Statistics ---------------------------------------------------------

# Perform comparisons for all combinations of treatments, depths, and time points
DOC_stats_test <- data.frame()

for (treatment1 in unique(DOC_conc$Treatment)) {
  for (treatment2 in unique(DOC_conc$Treatment)) {
    for (depth1 in unique(DOC_conc$Depth)) {
      for (depth2 in unique(DOC_conc$Depth)) {
        for (timepoint1 in unique(DOC_conc$Timepoint)) {
          for (timepoint2 in unique(DOC_conc$Timepoint)) {
            # Ensure treatment1, depth1, and timepoint1 are less than or equal to treatment2, depth2, and timepoint2
            if (treatment1 <= treatment2 & depth1 <= depth2 & timepoint1 <= timepoint2) {
              p_value <- perform_t_test(
                mean1 = mean(DOC_conc$DOC_conc[DOC_conc$Treatment == treatment1 & DOC_conc$Depth == depth1 & DOC_conc$Timepoint == timepoint1]),
                mean2 = mean(DOC_conc$DOC_conc[DOC_conc$Treatment == treatment2 & DOC_conc$Depth == depth2 & DOC_conc$Timepoint == timepoint2]),
                sd1 = mean(DOC_conc$DOC_conc_sd[DOC_conc$Treatment == treatment1 & DOC_conc$Depth == depth1 & DOC_conc$Timepoint == timepoint1]),
                sd2 = mean(DOC_conc$DOC_conc_sd[DOC_conc$Treatment == treatment2 & DOC_conc$Depth == depth2 & DOC_conc$Timepoint == timepoint2]),
                n1 = sum(DOC_conc$Obs_numb[DOC_conc$Treatment == treatment1 & DOC_conc$Depth == depth1 & DOC_conc$Timepoint == timepoint1]),
                n2 = sum(DOC_conc$Obs_numb[DOC_conc$Treatment == treatment2 & DOC_conc$Depth == depth2 & DOC_conc$Timepoint == timepoint2])
              )
              comparison <- data.frame(
                treatment1 = treatment1,
                depth1 = depth1,
                timepoint1 = timepoint1,
                treatment2 = treatment2,
                depth2 = depth2,
                timepoint2 = timepoint2,
                p_value = p_value
              )
              DOC_stats_test <- rbind(DOC_stats_test, comparison)
              
            }
          }
        }
      }
    }
  }
}


# Print summary table
print(DOC_stats_test)

#write.csv(DOC_stats_test, "DOC_Stats_test.csv", row.names = F)


## 2.2. Plotting -----------------------------------------------------------

DOC_plot_data <- DOC_stats_test %>% 
  filter(treatment1 == 'Control',
         treatment2 == 'Glucose',
         timepoint1 == timepoint2,
         depth1 == depth2) %>% 
  mutate(signif = case_when(p_value < 0.05 ~ '*',
                            p_value < 0.01 ~ '**',
                            p_value < 0.001 ~ '***',
                            TRUE ~ 'NS'),
         y_pos = 600) %>% 
  rename(Depth = depth1,
         Timepoints = timepoint1)

Depth_order <- c("Surface", "Deep")
DOC_conc$Depth <- factor(DOC_conc$Depth, levels = Depth_order)
DOC_plot_data$Depth <- factor(DOC_plot_data$Depth, levels = Depth_order)

# Plotting
DOC_plot <- ggplot(DOC_conc, aes(x=Treatment, y=DOC_conc, fill=Treatment)) +
  geom_bar(stat="identity", position="dodge") +
  geom_errorbar(aes(ymin=DOC_conc-DOC_conc_sd, ymax=DOC_conc+DOC_conc_sd), position="dodge", width=0.2) +
  facet_grid(~Depth * Timepoints) +
  labs(title="Porewater DOC Concentration", y="DOC Concentration in mg/L") +
  theme_classic() +
  geom_text(data = DOC_plot_data,
            aes(x = 1.5,
                y = y_pos,
                label = signif),
            size = 3,
            inherit.aes = FALSE) +
  geom_segment(x = 1,
               xend = 2,
               y = 570,
               yend = 570) +
  geom_segment(x = 1,
               xend = 1,
               y = 570,
               yend = 555) +
  geom_segment(x = 2,
               xend = 2,
               y = 570,
               yend = 555)+
  theme(plot.title = element_text(hjust = 0.5, size = 11),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 9),
        axis.ticks.x =  element_blank())


DOC_plot
ggsave("DOC_plot.svg", DOC_plot, width = 6, height = 4, device = "svg")


# 3. CO2 Production potential - rate --------------------------------------
CO2_Conc <- data %>%
  select(Treatment, Depth, Timepoints, CO2_produced, CO2_produced_sd, Obs_numb)


## 3.1. Linear regression --------------------------------------------------
CO2_s1 <- CO2_Conc %>%
  filter(Timepoints != c(56, 70)) %>%
  group_by(Treatment, Depth) %>%
  nest() %>%
  mutate(mod = purrr::map(data, function(x) lm(CO2_produced ~ Timepoints, data = x))) %>%
  mutate(coef = purrr::map_dbl(mod, function(x) coefficients(x)[2])) %>%
  dplyr::select(Treatment, Depth, coef)

CO2_s2 <- CO2_Conc %>%
  filter(Timepoints != c(70)) %>%
  group_by(Treatment, Depth) %>%
  nest() %>%
  mutate(mod = purrr::map(data, function(x) lm(CO2_produced ~ Timepoints, data = x))) %>%
  mutate(coef = purrr::map_dbl(mod, function(x) coefficients(x)[2])) %>%
  dplyr::select(Treatment, Depth, coef)

CO2_s3 <- CO2_Conc  %>%
  group_by(Treatment, Depth) %>%
  nest() %>%
  mutate(mod = purrr::map(data, function(x) lm(CO2_produced ~ Timepoints, data = x))) %>%
  mutate(coef = purrr::map_dbl(mod, function(x) coefficients(x)[2])) %>%
  dplyr::select(Treatment, Depth, coef)

all_mod <- rbind(CO2_s1, CO2_s2) %>%
  rbind(CO2_s3) %>%
  magrittr::inset('vol', value = 22.19755827) %>%
  mutate(ppm_mlww_d = coef/vol) %>%
  magrittr::inset('head', value = 97.80244173) %>%
  mutate(mmol_mlww_d = head/(0.08205*298)*1/1000*ppm_mlww_d*1/1000) %>%
  group_by(Treatment, Depth) %>%
  reframe(mu = mean(mmol_mlww_d), 
          sig = sd(mmol_mlww_d)) %>%
  modify_at('Depth', factor, levels = c('Surface', 'Deep'))



## 3.2. Plotting -----------------------------------------------------------

Rate_plot <- ggplot(all_mod, aes(x = Treatment, y = mu, fill = Treatment)) +
  facet_wrap(~Depth)+
  geom_bar(stat="identity", position="dodge") +
  geom_errorbar(aes(ymin=mu-sig, ymax=mu+sig), position="dodge", width=0.2) +
  facet_grid(~Depth)+
  labs(title="CO"[2]~"Production Rate", y= "CO"[2]~"Production in mmoles cm"^"-3"~"day"^"-1") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 11),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 9),
        axis.ticks.x = element_blank())

Rate_plot
ggsave("Rate_plot.svg", Rate_plot, width = 6, height = 4, device = "svg")


# 4. CO2 production source - glucose treated samples ----------------------

CO2_data <- data %>%
  select(Treatment, Depth, Timepoints, CO2_produced, delta13C_total, delta13C_soil)

filtered_CO2_source <- CO2_data %>%
  filter(Treatment == "Glucose") %>%
  filter(Timepoints == c("7", "70"))


## 4.1. Calculations -------------------------------------------------------
Production_source <- filtered_CO2_source %>%
  mutate(RSOM = CO2_produced*((delta13C_total - 31.6547701058918)/(delta13C_soil - 31.6547701058918))) %>%
  mutate(Rg = CO2_produced - RSOM) %>%
  mutate(SOM = (RSOM/CO2_produced) * 100) %>%
  mutate(Glucose = 100- SOM)


## 4.2. Plotting -----------------------------------------------------------
plot_data <- Production_source %>%
  select(Treatment, Depth, Timepoints, SOM, Glucose)%>%
  pivot_longer(cols = c(SOM, Glucose),
               names_to = "Source",
               values_to = "Percent") %>%
  modify_at('Depth', factor, levels = c('Surface', 'Deep'))
plot_data$Timepoints <- as.factor(plot_data$Timepoints)


CO2_source_plot <- ggplot(plot_data, aes(x = Timepoints, y = Percent, fill = Source)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~ Depth)+
  scale_fill_manual(values = c("SOM" = "yellow3", "Glucose" = "blue3"))+
  labs(title="CO"[2]~"Production Source", y= "Percent of CO"[2]~"Produced") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 11),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 9))+
  scale_y_continuous(labels = function(x) paste0(x, "%"))
CO2_source_plot

ggsave("CO2_source_plot.svg", CO2_source_plot, width = 6, height = 4, device = "svg")


# 5. Priming of Old SOM - glucose induced ---------------------------------

## 5.1. Calculations -------------------------------------------------------
control_data <- CO2_data %>%
  filter(Treatment == "Control") %>%
  filter(Timepoints == c("7", "70"))%>%
  select(Depth, Timepoints, CO2_produced, delta13C_total)%>%
  rename(co2_produced_control = CO2_produced)%>%
  rename(delta13C_total_control = delta13C_total)

Priming <- Production_source%>%
  rename(CO2_produced_glucose = CO2_produced)%>%
  rename(delta13C_total_glucose = delta13C_total)%>%
  select(-Treatment)%>%
  left_join(control_data, by=c("Depth", "Timepoints"))%>%
  mutate(Rp = RSOM - co2_produced_control)%>%
  mutate(Priming = (Rp/CO2_produced_glucose) * 100)


## 5.2. Plotting -----------------------------------------------------------
priming_plot_data <- Priming %>%
  select(Depth, Timepoints, Priming)%>%
  modify_at('Depth', factor, levels = c('Surface', 'Deep'))

priming_plot <- ggplot(priming_plot_data, aes(x = Depth, y = Priming, fill = Depth)) +
  facet_wrap(~ Timepoints)+
  geom_bar(stat="identity", position="dodge")+
  scale_fill_manual(values = c("Surface" = "#4FB766", "Deep" = "#977E59"))+
  labs(title="Glucose Induced CO"[2]~"Production", y= "Percent of CO"[2]~"Produced") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 11),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 9))+
  scale_y_continuous(labels = function(x) paste0(x, "%"))
priming_plot

ggsave("priming_plot.svg", priming_plot, width = 6, height = 4, device = "svg")


# Merging overall plots ---------------------------------------------------

# Extract individual legends from each plot
legend_CO2 <- cowplot::get_legend(CO2_plot)
legend_CO2_source <- cowplot::get_legend(CO2_source_plot)
legend_priming <- cowplot::get_legend(priming_plot)

# Combine the individual legends into one
combined_legend <- grid.arrange(legend_CO2, legend_CO2_source, legend_priming,
                                ncol = 1, nrow = 3)

#overall figure


overall_figure <- plot_grid(CO2_plot + theme(legend.position = "none"),
          DOC_plot + theme(legend.position = "none"), 
          Rate_plot + theme(legend.position = "none"), 
          CO2_source_plot + theme(legend.position = "none"), 
          priming_plot + theme(legend.position = "none"), 
          combined_legend,
          labels= c("A.", "B.", "C.", "D.", "E."),
          label_size = 11,
          ncol = 3, nrow = 2)

overall_figure
ggsave("overall_figure_with_combined_legend.svg", overall_figure, width = 10, height = 6, device = "svg")
