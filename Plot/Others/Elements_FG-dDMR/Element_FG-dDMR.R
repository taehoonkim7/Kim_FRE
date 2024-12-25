############################################################################
### Import libraries
library(stats)
library(tidyverse)
library(ggplot2)
library(dplyr)
############################################################################
### Import files
data <- read.csv("Element_FG-dDMR.csv")

############################################################################
### Plotting

data$Element = factor (data$Element,
                       levels = c("2kb-upstream","5'UTR","CDS","3'UTR","Intronic","2kb-downstream","Intergenic"))
data$if_TE = factor (data$if_TE,
                     levels = c("TE","Non-TE"))

p <- ggplot(data, aes(x = Element, y = Count, fill = if_TE)) +
  geom_bar(stat = "identity", position = "stack",
           color = "black", linewidth = 0.2, width=0.75)+
  theme_classic()+
  scale_fill_manual(values = c("#f8c578","#6b89d7"))+
  labs(y = "Counts",title = "")+
  theme(axis.line = element_line(color="black", linewidth = 0.2),
        axis.title.y = element_text(color = "black", size = 7),
        axis.title.x = element_blank(),
        axis.text.y = element_text(color = "black", size = 5),
        axis.text.x = element_text(color = "black", size = 5,
                                   angle=45, vjust=1, hjust=1),
        axis.ticks.y = element_line(color="black", linewidth = 0.2),
        axis.ticks.x = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.title = element_blank(),
        legend.background = element_blank(),
        legend.text = element_text(color = "black", size = 5),
        legend.key.size = unit(2.5,"mm"),
        legend.margin = margin(2,2,2,2,"mm"),
        legend.position = "inside",
        legend.position.inside = c(0.27,0.95))+
  scale_y_continuous(expand = c(0, 0), limits = c(0,300))

pdf(file = "Element_FG-dDMR.pdf", width = 1.4, height = 1.7)
p
dev.off()

