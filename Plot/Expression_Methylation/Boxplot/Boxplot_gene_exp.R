############################################################################
### Import libraries
library (ggplot2)       # if failed to load, try: install.packages("ggplot2")
library (tidyverse)     # if failed to load, try: install.packages("tidyverse")
library (ggthemes)      # if failed to load, try: install.packages("ggthemes")
library (ggpubr)        # if failed to load, try: install.packages("ggpubr")
library (cowplot)        # if failed to load, try: install.packages("cowplot")
library (reshape2)

############################################################################
### Import files
lst_gene   <- list.files (path = "./Gene", pattern = ".txt")
Desc_data <- read.csv ("Sample_Description_GeneTE.csv")
Exp_data <- read.csv ("./Plots_exp/FG_TPM_decile.csv")

############################################################################
### Pre-process files
##Import data
rownames(Desc_data) <- Desc_data[,1]
data <- data.frame (matrix(ncol = 7, nrow = 0))

##Merge data
for (i in lst_gene){
  temp_data <- read.table(paste(getwd(),"/Gene/",i,sep=""))
  temp_data["Sample"]  <- Desc_data[i,"Sample"]
  temp_data["Context"] <- Desc_data[i, "Context"]
  temp_data["Family"]  <- Desc_data[i, "Family"]
  
  data <- rbind(data, temp_data)
}

colnames (data) <- c("ID","mC","C","Percentage","Sample","Context","Family")
head(data)
## Remove temporary files
if (exists("Desc_data")) {rm (Desc_data)}
if (exists("temp_data")) {rm (temp_data)} 
if (exists("i"))         {rm (i)} 
if (exists("lst_gene"))  {rm (lst_gene)} 

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
  select(ID, mC, C, Sample, Context)  %>%
  group_by(ID, Sample, Context) %>%
  summarise(mC = sum(mC), C = sum(C)) %>%
  mutate (Percentage = mC/(mC+C)*100) %>%
  select(ID, mC, C, Percentage, Sample, Context)

data <- na.omit(data) 

# Assign expression
Exp_data <- melt(Exp_data, 
                 id.vars = "gene", variable.name = "sample", value.name = "Expression")
data <- left_join(data, Exp_data,
                  by = c("ID" = "gene",
                         "Sample" = "sample"))


# Set categories
samples = c("FG1","FG2","FG3","FG4","FG5","FG6","FL")

contexts = c("CpG","CHG","CHH")
contexts_label = c("mCG","mCHG","mCHH")
contexts_mapping = setNames(contexts_label, contexts)

expressions = c("NoExp","2nd","4th","6th","8th","10th")
expressions_label = c("No expression", "2nd decile", "4th decile", 
                      "6th decile", "8th decile", "10th decile")
expressions_mapping = setNames(expressions_label, expressions)

exp_palette = c("#e86473","#f4a979","#e4d189",
                "#9dc489","#65a14f","#666666")

# Apply categories
data <- data %>%
  mutate (Context = contexts_mapping[Context],
          Expression = expressions_mapping[Expression])

data$Sample <- factor(data$Sample, levels = samples)
data$Context <- factor(data$Context, levels = contexts_label)
data$Expression <- factor (data$Expression, levels = expressions_label)

# Export table
write.csv (data, file = "data.csv")

############################################################################
###Boxplot - Genes by expression
p <- data %>%
  ggplot(aes(x = Sample, y = Percentage, fill = Expression))+
  stat_boxplot(geom="errorbar", linewidth = 0.3,
               position = position_dodge(width = 1))+
  geom_boxplot(color = "black", linewidth = 0.3, outliers = FALSE,
               position = position_dodge(width = 1))+
  theme_classic()+
  scale_fill_manual(values = rev(exp_palette))+
  facet_grid(Context~Sample, scales = "free") +
  coord_cartesian(clip = 'off') + 
  labs(y = "%Methylation") +
  theme(axis.line = element_line(color="black", linewidth = 0.3),
        axis.title.x = element_blank(),
        axis.title.y = element_text(color="black", size = 6),
        axis.text.y = element_text(color="black", size = 5),
        axis.text.x = element_blank(),
        axis.ticks.y = element_line(color="black", linewidth = 0.3),
        axis.ticks.x = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        panel.spacing = unit(1,"mm"),
        strip.background = element_blank(),
        strip.text = element_text(color="black", size = 6),
        legend.background = element_blank(),
        legend.spacing = unit(0,"mm"),
        legend.title = element_blank(),
        legend.key.size = unit(3,"mm"),
        legend.text = element_text(color="black", size = 5)
  )+
  scale_y_continuous(expand = c(0, 1)) + 
  scale_x_discrete(expand = c(0.3, 0.25))

pdf(file = paste("boxplot_gene_exp.pdf",sep=""), width = 7, height = 2)
print(p)
dev.off()

