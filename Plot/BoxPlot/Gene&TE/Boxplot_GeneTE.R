############################################################################
### Import libraries
library (ggplot2)
library (tidyverse)
library (ggthemes)
library (ggpubr)
library (cowplot)

############################################################################
### Import files
# Import all csv file in the directory
lst_gene   <- list.files (path = "./Gene", pattern = ".txt")
lst_TE     <- list.files(path = "./TE_split", pattern = "txt")
Desc_file <- read.csv ("Sample_Description.csv")

############################################################################
### Pre-process files
##Import data
rownames(Desc_file) <- Desc_file[,1]
data <- data.frame (matrix(ncol = 7, nrow = 0))

##Merge data
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

colnames (data) <- c("ID","mC","C","Percentage","Sample","Context","Family")

## Remove temporary files
if (exists("Desc_file")) {rm (Desc_file)}
if (exists("temp_data")) {rm (temp_data)} 
if (exists("i"))         {rm (i)} 
if (exists("lst_gene"))  {rm (lst_gene)} 
if (exists("lst_TE"))    {rm (lst_TE)} 


## Process data
# Merge replicates of the same samples
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
  select(ID, mC, C, Sample, Context, Family, GeneTE)  %>%
  group_by(ID, Sample, Context, Family) %>%
  summarise(mC = sum(mC), C = sum(C)) %>%
  mutate (Percentage = mC/(mC+C)*100) %>%
  select(ID, mC, C, Percentage, Sample, Context, GeneTE)

data$Sample <- factor(data$Sample,
                      levels = c("FG1","FG2","FG3","FG4","FG5","FG6","FL"))
data$Context <- factor(data$Context,
                       levels = c("CpG","CHG","CHH"))

# Remove data where mC and C are both 0
# total (7,446,138) - na_rows (257,203) = remaining (7,188,935)
data <- na.omit(data) 

############################################################################
###Boxplot - GeneTE
fill_pal <- c("#f17272", "#f8c578", "#d1e47b", "#aad9bf",
                    "#6b89d7", "#ae88c8", "#666666")
context <- c("CG","CHG","CHH")
names(context) <- c("CpG","CHG","CHH")

p_gene <- data %>%
  filter(GeneTE == "Gene") %>%
  ggplot(aes(x = Sample, y = Percentage, fill = Sample))+
  stat_boxplot(geom="errorbar", color = "black", linewidth = 0.3)+
  geom_boxplot(outliers = FALSE, color = "black", linewidth = 0.3)+
  theme_classic()+
  scale_fill_manual(values = fill_pal)+
  facet_grid(Context~., scales = "free") +
  coord_cartesian(clip = 'off') + 
  labs(y = "%Methylation", title = "Gene") +
  theme(plot.title = element_text(color="black", size = 7),
        axis.line = element_line(color="black", linewidth = 0.3),
        axis.title.x = element_blank(),
        axis.title.y = element_text(color="black", size = 6),
        axis.text.y = element_text(color="black", size = 5),
        axis.text.x = element_text(color="black", size = 6, angle = 90, vjust = 0.5),
        axis.ticks.y = element_line(color="black", linewidth = 0.3),
        axis.ticks.x = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_blank()
  )+
  scale_y_continuous(expand = c(0, 1))

p_TE <- data %>%
  filter(GeneTE == "TE") %>%
  ggplot(aes(x = Sample, y = Percentage, fill = Sample))+
  stat_boxplot(geom="errorbar", color = "black", linewidth = 0.3)+
  geom_boxplot(outliers = FALSE, color = "black", linewidth = 0.3)+
  theme_classic()+
  scale_fill_manual(values = fill_pal)+
  facet_grid(Context~., scales = "free",
             labeller = labeller(Context = context)) +
  coord_cartesian(clip = 'off') + 
  labs(title = "TE") +
  theme(plot.title = element_text(color="black", size = 7),
        axis.line = element_line(color="black", linewidth = 0.3),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(color="black", size = 5),
        axis.text.x = element_text(color="black", size = 6, angle = 90, vjust = 0.5),
        axis.ticks.y = element_line(color="black", linewidth = 0.3),
        axis.ticks.x = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_text(color="black", size = 6)
  )+
  scale_y_continuous(expand = c(0, 1))

p_combined <- ggarrange(p_gene, p_TE, 
                        ncol = 2, nrow = 1, widths = c(2, 2))

pdf(file = paste("boxplot_geneTE.pdf",sep=""), width = 2.7, height = 2.2)
print(p_combined)
dev.off()
