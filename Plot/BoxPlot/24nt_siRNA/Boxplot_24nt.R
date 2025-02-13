############################################################################
### Import libraries
library (ggplot2)
library (tidyverse)
library (ggthemes)
library (ggpubr)
library (cowplot)
library (rstatix)
library (forcats)       
library (reshape2)

############################################################################
### Import files
#Import data
data <- read.csv("RPM.24nt.csv", header = T)

############################################################################
### Pre-process files
#Set categories
regions <- c("FG-dDMR_FG-siren","FG-dDMR_only","FG-siren_only","Shuffle")
samples <- c("FG1","FG2","FG3","FG4","FG5","FG6","FL")

#Process data
data <- data %>%
  select (-c(Chr,Start,End))

data_melted <- data %>%
  melt (id.vars = "Region", variable.name = "Sample", value.name = "RPM") %>%
  mutate(log10RPM = log10(1+RPM))

data_melted$Sample <- factor(data_melted$Sample, levels = samples)
data_melted$Region <- factor(data_melted$Region, levels = regions)

############################################################################
###Boxplot of body methylation

fill_pal <- c("#f17272", "#f8c578", "#d1e47b", "#aad9bf",
              "#6b89d7", "#ae88c8", "#666666")

#24nt-siRNA
p <- data_melted %>%
  ggplot(aes(x = Sample, y = log10RPM, fill = Sample))+
  stat_boxplot(geom="errorbar", color = "black", linewidth = 0.2)+
  geom_boxplot(outliers = FALSE, color="black", linewidth = 0.2)+
  theme_classic()+
  scale_fill_manual(values = fill_pal)+
  facet_grid(Region~.) +
  coord_cartesian(clip = 'off') + 
  labs(y = "log10(1+RPM)", title = "24nt-siRNA") +
  theme(axis.line = element_line(color="black", linewidth = 0.2),
        axis.title.x = element_blank(),
        axis.title.y = element_text(color="black", size = 6),
        axis.text.y = element_text(color="black", size = 5),
        axis.text.x = element_text(color="black", size = 5,
                                   angle = 90, vjust = 0.5),
        axis.ticks.y = element_line(color="black", linewidth = 0.2),
        axis.ticks.x = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_text(color="black", size = 5),
        plot.title = element_blank(),
  )+
  scale_y_continuous(expand = c(0, 0.5))

pdf(file = "./boxplot_sRNA.pdf", width = 1.4, height = 2.2)
print(p)
dev.off()
