############################################################################
### Import libraries
library (ggplot2)
library (tidyverse)
library (ggthemes)
library (ggpubr)
library (cowplot)
library (lemon)
############################################################################
### Import files
# Set directory
# Import all csv file in the directory
lst_gene   <- list.files (path = "./Gene", pattern = ".txt")
lst_TE     <- list.files(path = "./TE_split", pattern = "txt")

Desc_file <- read.csv ("Sample_Description.csv")

############################################################################
### Pre-process files
#Import data
rownames(Desc_file) <- Desc_file[,1]
data <- data.frame (matrix(ncol = 7, nrow = 0))

#Merge data
for (i in lst_gene){
  temp_data <- read.table(paste(getwd(),"/Gene/",i,sep=""))
  temp_data["Sample"]  <- Desc_file[i,"Sample"]
  temp_data["Context"] <- Desc_file[i, "Context"]
  temp_data["Family"]  <- Desc_file[i, "Family"]
  
  data <- rbind(data, temp_data)
}

for (i in lst_TE){
  temp_data <- read.table(paste(getwd(),"/TE_split/",i,sep=""))
  temp_data["Sample"]  <- Desc_file[i,"Sample"]
  temp_data["Context"] <- Desc_file[i, "Context"]
  temp_data["Family"]  <- Desc_file[i, "Family"]
  
  data <- rbind(data, temp_data)
}

colnames (data) <- c("Position","mC","C","Percentage",
                     "Sample","Context","Family")
############################################################################
### Metaplot
  
## Data processing
data <- data %>%
  filter(!Sample %in% c("FG1-2_S2")) %>%
  mutate(Sample = case_when(Sample %in% c("FG1-1_S1")         ~ "FG1",
                            Sample %in% c("FG2-1_S3","FG2-2_S4") ~ "FG2",
                            Sample %in% c("FG3-1_S5","FG3-2_S6") ~ "FG3",
                            Sample %in% c("FG4-1_S7","FG4-2_S8") ~ "FG4",
                            Sample %in% c("FG5-1_S9","FG5-2_S10") ~ "FG5",
                            Sample %in% c("FG6-1_S11","FG6-2_S12") ~ "FG6",
                            Sample %in% c("FL-1_S13","FL-2_S14") ~ "FL",
                            .default = NA)) %>%
  mutate(GeneTE = case_when(Family %in% c("Gene") ~ "Gene",
                            .default = "TE")) %>%
  select(Position, mC, C, Sample, Context, Family, GeneTE) 

data <- na.omit(data)

data$Sample <- factor(data$Sample,
                      levels = c("FL","FG6","FG5","FG4","FG3","FG2","FG1"))
data$Context <- factor(data$Context,
                       levels = c("CpG","CHG","CHH"))
data$Family <- factor(data$Family, 
                      levels = c("Gene","Copia","Gypsy","TRIM","LTRun",
                                 "CACTA","Mut","PIF","Tc1","hAT","TIRun",
                                 "Cen","LINE","SINE","Nonun","HEL","Rep","Sat"))
data$GeneTE <- factor(data$GeneTE,
                      levels = c("Gene","TE"))

data_GeneTE <- data %>%
  group_by(Position, Sample, Context, GeneTE) %>%
  summarise(mC = sum(mC), C = sum(C)) %>%
  mutate(Percentage = mC/(mC + C)*100) 

data_TE <- data %>%
  filter(!GeneTE %in% c("gene")) %>%
  group_by(Position, Sample, Context, Family) %>%
  summarise(mC = sum(mC), C = sum(C)) %>%
  mutate(Percentage = mC/(mC + C)*100) 

##Line graphs
col_pal <- c("black","#754696","#2A499B","#72BF94",
                          "#AECA2A","#F49F1E","#E41817")

p_gene <- data_GeneTE %>%
  filter(GeneTE == "Gene") %>%
  ggplot(aes(x = Position, y = Percentage, color = Sample))+
  geom_vline(xintercept = c(20,40),linetype="dashed", 
             linewidth = .3, color = "#666666")+
  geom_line(linewidth = .5)+
  theme_classic()+
  scale_color_manual(values = col_pal)+
  facet_grid(Context~., scales = "free") +
  coord_cartesian(ylim = c(0,NA), clip = 'off') + 
  labs(x = "", y = "%Methylation") +
  theme(axis.line = element_line(color="black", linewidth = .3),
        axis.text.y = element_text(color="black", size = 5),
        axis.text.x = element_text(color="black", size = 5,
                                   angle = 90, vjust = 0.5),
        axis.title.y = element_text(color="black", size = 7),
        axis.ticks = element_line(color="black", linewidth = .3),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_blank()
  )+
  scale_x_continuous(breaks = c(0,20,40,60),
                     labels = c("-2kb", "", "", "+2kb"))+
  scale_y_continuous(expand = c(0, 1))

p_TE <- data_GeneTE %>%
  filter(GeneTE == "TE") %>%
  ggplot(aes(x = Position, y = Percentage, color = Sample))+
  geom_vline(xintercept = c(20,40),linetype="dashed", 
             linewidth = .3, color = "#666666")+
  geom_line(linewidth = .5)+
  theme_classic()+
  scale_color_manual(values = col_pal)+
  facet_grid(Context~., scales = "free") +
  coord_cartesian(clip = 'off') + 
  labs(x = "", y = "") +
  theme(axis.line = element_line(color="black", linewidth = 0.3),
        axis.text.y = element_text(color="black", size = 5),
        axis.text.x = element_text(color="black", size = 5,
                                   angle = 90, vjust = 0.5),
        axis.ticks = element_line(color="black", linewidth = 0.3),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_blank()
  )+
  scale_x_continuous(breaks = c(0,20,40,60),
                     labels = c("-2kb", "", "", "+2kb"))+
  scale_y_continuous(expand = c(0, 2))

p_combined <- ggarrange(p_gene, p_TE,
                        ncol = 2, nrow = 1,
               widths = c(3, 3))

pdf(file = paste("metaplot_geneTE.pdf",sep=""), width = 3.5, height = 2.5)
print(p_combined)
dev.off()
