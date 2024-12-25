############################################################################
### Import libraries
library (ggplot2)
library (tidyverse)
library (reshape2)

############################################################################
### Import files
# Import & process data
data = read.csv("PCA_ChIP_osclsy3.tab",sep = "\t", header = T, skip = 1)

# Set categories
components <- c(1,2)
samples <- c("GFP1","GFP2","MYC1","MYC2","Input1","Input2")

# Process data
data <- data %>%
  filter (Component %in% components) %>%
  select (-c(Component,Eigenvalue))

data <- data.frame(t(data))
colnames(data) <- c("PC1","PC2")
data$Sample <- row.names(data)
data$Sample <- factor(data$Sample, 
                      levels = samples)

############################################################################
###Plotting
palette <- c("#488f31","#83af70","#de425b","#ec838a","#004c6d","#5886a5")
shapes <- rep(c(1,2),3)

#PCA plot
p <- ggplot (data, aes(x = PC1, y = PC2, col = Sample, shape = Sample)) + 
  geom_point() + 
  scale_color_manual(values = palette) + 
  scale_shape_manual(values = shapes) +
  theme_linedraw() + 
  labs(x = "PC1", y = "PC2") +
  theme(line = element_line(color="black", linewidth = 0.2),
        axis.title = element_text(color="black", size = 6),
        axis.text = element_text(color="black", size = 5),
        axis.ticks = element_line(color="black", linewidth = 0.2),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "right",
        legend.background = element_blank(),
        legend.title = element_blank(),
        legend.key.size = unit(3,'mm'),
        legend.text = element_text(color="black", size = 5)
  ) + 
  scale_y_continuous(expand = c(0,0.1)) + 
  scale_x_continuous(limits = c(0,0.45)) 

pdf(file = "PCA_ChIP_osclsy3.pdf", width = 2.5, height = 1.7)
print(p)
dev.off()

