library(tidyverse)
library(readxl)
library(ggh4x)


root_comp <- read_xlsx('Raw_data/Metadata_all.xlsx', sheet = 6) %>% 
  column_to_rownames(var = 'metabolite')
#root_comp is a dataframe with the weights (could be m/z or neutrualized weight) of all features by each method, with method as columns


#count the total features for include that in legend.

root_comp_long <- tidyr::gather(root_comp, key = "Method", value = "Exact_mass", na.rm = TRUE)


counts <- root_comp_long %>% 
  group_by(Method) %>% 
  summarize(n = n())


# Merge counts with root_comp_long
root_comp_long <- left_join(root_comp_long, counts, by = "Method")

root_comp_long$Method <- paste0(root_comp_long$Method, " (n = ", root_comp_long$n, ")")
legend_order <- c("NMR (n = 34)", "GC (n = 110)","LC (n = 514)")

# Specify colors, RGB works good for overlapped color, Specify labels for secondary axis.

plot_root_comp <- ggplot(root_comp_long, aes(x = Exact_mass, fill = Method)) +
  geom_histogram(data = filter(root_comp_long, Method %in% c("LC (n = 514)", "NMR (n = 34)","GC (n = 110)")),
                 position = "identity", binwidth = 25, alpha = 0.5, aes(y = ..count..)) +
  scale_fill_manual(values = c("LC (n = 514)" = "blue", "NMR (n = 34)" = "green", "GC (n = 110)" = "red"), limits = legend_order) + 
  labs(title = "Histogram of all filtered metabolites",
       x = "Molecular Weight",
       y = "Frequency") +
  cowplot::theme_minimal_grid() +
  theme(legend.position = "bottom", legend.direction = "horizontal")+
  scale_y_continuous(
    name = "Frequency")  

plot_root_comp + ggh4x::stat_theodensity(
  aes(color = Method,
      y = after_stat(count * 25)),
  linewidth = 1
) +
  scale_color_manual(values = c("LC (n = 514)" = "blue", "NMR (n = 34)" = "green", "GC (n = 110)" = "red"), limits = legend_order)

ggsave('plot_root_comp.svg', dpi = 300, height = 5, width = 12)  
