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
# Import data
file = "motif_enrich_TE.csv"
data = read.csv(file, header = T)
unique(data$Region)
regions = c("LTR/Ty1-copia", "LTR/Ty3-gypsy", "LTR/TRIM", "LTR/unknown",
            "nonLTR/LINE", "nonLTR/SINE", "nonLTR/unknown",
            "TIR/CACTA", "TIR/Mutator", "TIR/PIF_Harbinger", "TIR/Tc1_Mariner",
            "TIR/hAT", "TIR/unknown", "nonTIR/helitron", "Centromeric", "Repeat", "Satellite")
motifs = unique (data$Motif)
offsets = unique (data$Offset)

data$Region = factor(data$Region, levels = rev(regions))
data$Motif = factor(data$Motif, levels = motifs)
data$Offset = factor(data$Offset, levels = offsets)

############################################################################
## Data preparation
# process data
data <- data %>%
  filter (Offset == "±0.0kb") %>%
  filter (!Region %in% c("Centromeric", "Repeat", "Satellite")) %>%
  mutate (Perc = Count/Total*100)

############################################################################
###Plotting
high_col <- "#A42B2A"
mid_col <- "#f3f4cc"  
low_col <- "#f1f1f1"

p <- ggplot (data, aes(x = Motif, y = Region, fill = Perc)) + 
  geom_tile() + 
  coord_fixed() +
  scale_fill_gradientn(colors = c(low_col, mid_col, high_col),
                       limits = c(0, 60), breaks = c(0,20,40,60),
                       values = rescale(c(0,10,60))) + 
  theme_classic() + 
  theme(text = element_text(color="black", size = 5),
        axis.title.x = element_blank(),
        axis.title.y = element_text(color="black", size = 6),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_text(color="black", size = 5),
        axis.text.x = element_text(color="black", size = 5, 
                                   angle = 90, hjust = 0, vjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
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
  labs(fill = "%Enrichment", y = "TE families")

pdf(file = "heatmap_motif_TE.pdf", width = 5, height = 2)
print(p)
dev.off()

