############################################################################
### Import libraries
library (ggplot2)
library (tidyverse)
library (ggthemes)
library (ggpubr)
library (cowplot)

############################################################################
### Import files
#Import data
data <- read.csv("confidence.csv", header = T)
colnames (data) <- c("Chain","Res_Number","Res_name","pLDDT")

# Data processing
data <- data %>%
  group_by(Chain, Res_Number, Res_name) %>%
  summarize(avg_pLDDT = mean(pLDDT)) %>%
  ungroup() %>%
  arrange(Res_Number)

############################################################################
### Lineplot

p <- data %>%
  ggplot(aes(x=Res_Number, y = avg_pLDDT)) +
  geom_vline(xintercept = 870,linetype="dashed", linewidth = .2, color = "#666666")+
  geom_line(linewidth = .2, color = "#488f31")+
  theme_classic()+
  labs(x = "Residue number", y = "pLDDT") +
  theme(axis.line = element_line(color="black", linewidth = .3),
        axis.text = element_text(color="black", size = 5),
        axis.title = element_text(color="black", size = 6),
        axis.ticks = element_line(color="black", linewidth = .3),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none"
  )+
  scale_x_continuous(breaks = c(0,400,800,1200))+
  scale_y_continuous(expand = c(0, 1))

pdf(file = "plot_pLDDT.pdf", width = 3, height = 1)
print(p)
dev.off()

