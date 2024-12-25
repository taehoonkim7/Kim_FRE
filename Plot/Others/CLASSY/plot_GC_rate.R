############################################################################
### Import libraries
library (ggplot2)
library (tidyverse)
library (ggthemes)
library (ggpubr)
library (cowplot)

############################################################################
### Import files
## Import data
data <- read.csv("OsCLSY3_GC.csv", header = T)
colnames (data) <- c("Gene","Position","GC_rate")

## Data processing
data <- data %>%
  mutate (GC_rate = GC_rate * 100) %>%
  arrange(Position)

############################################################################
### Metaplot

p <- data %>%
  ggplot(aes(x=Position, y = GC_rate)) +
  geom_line(linewidth = .2, color = "#de425b")+
  theme_classic()+
  labs(x = "Nucleotide number", y = "GC rate(%)") +
  theme(axis.line = element_line(color="black", linewidth = .3),
        axis.text = element_text(color="black", size = 5),
        axis.title = element_text(color="black", size = 6),
        axis.ticks = element_line(color="black", linewidth = .3),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none"
  )+
  scale_x_continuous(breaks = c(0,1000,2000,3000,4000), expand = c(0,100))+
  scale_y_continuous(breaks = c(30,50,70,90), expand = c(0, 1))

pdf(file = "plot_GC.pdf", width = 3, height = 1)
print(p)
dev.off()

