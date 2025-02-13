############################################################################
### Import libraries
library (ggplot2)
library (reshape2)
library (tidyverse)
library (raster)
library (scales)
library (factoextra)
library (vegan)

############################################################################
### Import files
# Set categories
regions = c("FG-dDMR_FG-siren","FG-dDMR_only","FG-siren_only")
samples = c("FG1","FG2","FG3","FG4","FG5","FG6","FL")
motifs = c("FRE_0","FRE_1","FRE_2","FRE_3")

# Import & process data
total = read.csv("24nt_total.csv", header = T)

# data_exp (24nt-siRNA expression)
data_exp = data.frame (matrix(ncol = 12, nrow = 0))

for (region in regions){
  temp_exp = read.csv (paste0("count.",region,".txt"), sep = "\t")
  temp_exp$Region = region
  temp_exp <- temp_exp %>%
    dplyr::select (Chr, Start, End, Region,
                   FG1, FG2, FG3, FG4, FG5, FG6, FL) %>%
    mutate (Chr = as.integer(gsub("Chr", "", Chr))) %>%
    arrange (Chr, Start)
  temp_exp$id = row.names(temp_exp)
  temp_exp <- temp_exp %>%
    mutate (id = paste0(Region,"_",id))
  data_exp <- rbind(data_exp, temp_exp)
  rm(temp_exp)
}

colnames(data_exp) = c("Chr","Start","End","Region",
                       "FG1","FG2","FG3","FG4","FG5","FG6","FL","ID")
data_exp <- data_exp %>%
  melt (id.vars = c("ID","Chr","Start","End","Region"),
        variable.name = "Sample", value.name = "Count") %>%
  left_join(total, by = "Sample") %>%
  mutate (RPM = Count/Total_read * (10^6),
          logRPM = log10(RPM + 1)) 

data_exp <- data_exp %>%
  dplyr::select (ID, Region, Sample, logRPM) %>%
  pivot_wider (names_from = Sample, values_from = logRPM) %>%
  arrange (FG6)
ID_order = unique(data_exp$ID)

#data_count (FREs)
data_count = data.frame (matrix(ncol = 12, nrow = 0))
for (region in regions){
  temp_count = read.csv (paste0("FRE.",region,".bed"), sep = "\t", header = F)
  temp_count$V6 = region
  colnames (temp_count) = c("Chr","Start","End","Motif","Count","Region")
  temp_count <- temp_count %>%
    pivot_wider(names_from = Motif, values_from = Count) 
  temp_count$ID = row.names(temp_count)
  temp_count <- temp_count %>%
    mutate (ID = paste0(Region,"_",ID))
  data_count <- rbind(data_count, temp_count)
  rm(temp_count)
}
data_count <- data_count %>% 
  dplyr::select(ID, FRE_0, FRE_1, FRE_2, FRE_3) %>%
  arrange(FRE_0, FRE_1, FRE_2, FRE_3)
data_count$ID = factor (data_count$ID, levels = ID_order)

# Merge, then split data
data <- left_join (data_exp, data_count, by = "ID")
data$Region = factor(data$Region, levels = regions)
data$ID = factor(data$ID, levels = ID_order)

data_plot_exp <- data %>%
  dplyr::select(ID, Region, FG1, FG2, FG3, FG4, FG5, FG6, FL) %>%
  melt (id.vars = c("ID", "Region"),
        variable.name = "Sample", value.name = "logRPM") 
data_plot_exp$Sample = factor (data_plot_exp$Sample, levels = samples)

data_plot_count <- data %>%
  dplyr::select(ID, Region, FRE_0, FRE_1, FRE_2, FRE_3) %>%
  melt (id.vars = c("ID", "Region"),
        variable.name = "Motif", value.name = "Count") 
data_plot_count$Motif = factor (data_plot_count$Motif, levels = motifs)

############################################################################
###Plotting
exp_high <- "#A42B2A"
exp_low <- "#f3f4cc"
count_high <- "#488f31"
count_low <- "#f1f1f1"

# Heatmap - expression
p_exp <- ggplot (data_plot_exp, aes(x = Sample, y = ID, fill = logRPM)) + 
  geom_tile() + 
  facet_grid(Region~., scales = "free_y", space = "free_y") + 
  scale_fill_gradient2(limits = c(0, +4.0),
                       low = NA, high = exp_high, mid = exp_low,oob=squish) + 
  theme_classic() + 
  theme(text = element_text(color="black", size = 5),
        axis.title = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(color="black", size = 5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA, color = "black",linewidth = 0.3),
        panel.spacing = unit (0.2,"mm"),
        strip.background = element_blank(),
        strip.text = element_text(color="black", size = 5),
        legend.background = element_blank(),
        legend.title = element_text(color="black", size = 6, 
                                    angle = -90, hjust = 0.5),
        legend.title.position = "right",
        legend.key.size = unit (2,"mm"),
        legend.text = element_text(color="black", size = 5, margin = margin(l = 2))) + 
  labs(fill = "log10(1+RPM)")

pdf(file = "heatmap_24nt_region.pdf", width = 1.8, height = 3)
print(p_exp)
dev.off()

# Heatmap - FRE_count
p_count <- ggplot (data_plot_count, aes(x = Motif, y = ID, fill = Count)) + 
  geom_tile() + 
  facet_grid(Region~., scales = "free_y", space = "free_y") + 
  scale_fill_gradient2(limits = c(0, 1), breaks = c(0, 1),
                       low = NA, high = count_high, mid = count_low ,oob=squish) + 
  theme_classic() + 
  theme(text = element_text(color="black", size = 5),
        axis.title = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(color="black", size = 5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA, color = "black",linewidth = 0.3),
        panel.spacing = unit (0.2,"mm"),
        strip.background = element_blank(),
        strip.text = element_text(color="black", size = 5),
        legend.background = element_blank(),
        legend.title = element_text(color="black", size = 6, 
                                    angle = -90, hjust = 0.5),
        legend.title.position = "right",
        legend.key.size = unit (2,"mm"),
        legend.text = element_text(color="black", size = 5, margin = margin(l = 2))) + 
  labs(fill = "Count")

pdf(file = paste0("heatmap_FRE_region.pdf"), width = 1.7, height = 3)
print(p_count)
dev.off()

