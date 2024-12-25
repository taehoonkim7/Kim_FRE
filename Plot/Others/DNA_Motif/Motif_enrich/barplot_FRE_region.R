############################################################################
### Import libraries
library(ggplot2)
library(ggsignif)
library(dplyr)

############################################################################
### Data preparation
# Import data
file = "motif_enrich_region.csv"
data = read.csv(file, header = T)

# Process data
data <- data %>%
  filter (Motif == "FG-siren_01") %>%
  mutate (Perc = Count/Total*100)

regions = c("FG-dDMR_FG-siren","FG-dDMR_only","FG-siren_only","Shuffle")
offsets = unique (data$Offset)

data$Region <- factor (data$Region, levels = regions)
data$Offset <- factor (data$Offset, levels = offsets)


############################################################################
### Plotting

palette <- c("#f8a66c", "#f17272", "#ffd966", "#cccccc")
  
p <- ggplot(data, aes(x = Offset, y = Perc, fill = Region)) + 
  geom_bar (stat="identity", position=position_dodge(),
            color = "black", linewidth = 0.2) + 
  scale_fill_manual(values = palette) +
  theme_linedraw() +
  labs(y="%Occurrence") +
  theme(axis.line = element_line(color="black", linewidth = 0.2),
        text = element_text(color="black", size = 5),
        axis.title.x = element_blank(),
        axis.title.y = element_text(color="black", size = 6),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(color="black", size = 5),
        axis.text.x = element_text(color="black", size = 5),
        axis.ticks.y = element_line(color="black", linewidth = 0.2),
        axis.ticks.x = element_blank(),
        panel.border = element_blank(),
        legend.position = c(.3, .9),
        legend.background = element_blank(),
        legend.title = element_blank(),
        legend.key.size = unit (2,"mm"),
        legend.text = element_text(color="black", size = 5, margin = margin(l = 1))) +
  scale_y_continuous(limits = c(0, 125), breaks = c(0,25,50,75,100),
                     expand = c(0,1))
  
pdf(file = "barplot_FRE_region.pdf", width = 2, height = 2)
print(p)
dev.off()
  