############################################################################
### Import libraries
library (ggplot2)
library (dplyr)
library (reshape2)
library (pheatmap)
library (vegan)
library (scales)

############################################################################
### Import files
# Import data
data_FG <- read.csv("FG_OsCLSY.csv", header = TRUE)
data_Evo <- read.csv("EvoRepro_OsCLSY.csv", header = TRUE)
############################################################################
# Process data
data_FG <- data_FG %>%
  mutate (FG1 = rowMeans(data_FG[,c("FG1.1","FG1.2","FG1.3")]),
          FG2 = rowMeans(data_FG[,c("FG2.1","FG2.2","FG2.3")]),
          FG3 = rowMeans(data_FG[,c("FG3.1","FG3.2","FG3.3")]),
          FG4 = rowMeans(data_FG[,c("FG4.1","FG4.2","FG4.3")]),
          FG5 = rowMeans(data_FG[,c("FG5.1","FG5.2","FG5.3")]),
          FG6 = rowMeans(data_FG[,c("FG6.1","FG6.2","FG6.3")]),
          FL  = rowMeans(data_FG[,c("FL.1", "FL.2", "FL.3")])) %>%
  select (Gene,FG1,FG2,FG3,FG4,FG5,FG6,FL) 

sample_FG = colnames(data_FG[,-1])
sample_Evo = colnames(data_Evo[,-1])

data <- left_join(data_FG, data_Evo, by = join_by(Gene))
data$Blank = NA

data <- data %>%
  melt(id.vars = "Gene",
       variable.name = "Sample", value.name = "TPM")

data$Gene = factor(data$Gene,
                       levels = rev(unique(data$Gene)))
data$Sample = factor(data$Sample,
                       levels = c(sample_FG, "Blank", sample_Evo))

data <- data %>%
  mutate (log10TPM = log10(1+TPM))

high_col <- "#A42B2A"
mid_col <- "#F3F4CC"  
low_col <- NA

############################################################################
#ALL 

p <- data %>%
  ggplot(aes(y = Gene, x = Sample, fill = log10TPM)) + 
  geom_tile() + 
  geom_rect(aes(xmin = 0.5, xmax = 7.5, ymin = 0.5, ymax = 3.5),
            color = "black", fill = NA, linewidth = 0.2)+
  geom_rect(aes(xmin = 8.5, xmax = 21.5, ymin = 0.5, ymax = 3.5),
            color = "black", fill = NA, linewidth = 0.2)+
  scale_fill_gradient2(limits = c(0, 1.2), breaks = c(0.0, 0.6, 1.2),
                       low = low_col, high = high_col, mid = mid_col, oob=squish, na.value = NA) +
  scale_x_discrete(labels = function(x) ifelse(x == "Blank", "", x))+
  theme_void() +
  labs(title = "", x = "", y = "", fill = "log10(1+TPM)") +
  theme(axis.text.x = element_text(size = 6, color = "black", angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 6, color = "black", hjust = 1),
        legend.title = element_text(size = 6, color = "black"),
        legend.text = element_text(size = 5, color = "black"),
        legend.key.size = unit(2.5,'mm'),
        legend.position = "top") 
  
pdf ("expression_OsCLSY.pdf", width = 3, height = 1.2)
p
dev.off()

