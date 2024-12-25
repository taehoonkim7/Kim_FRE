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
lst_input   <- list.files (path = "./Input", pattern = ".txt")
Desc_file <- read.csv ("Sample_Description.csv")

############################################################################
### Pre-process files
#Import data
rownames(Desc_file) <- Desc_file[,1]
data <- data.frame (matrix(ncol = 7, nrow = 0))

#Merge data
for (i in lst_input){
  temp_data <- read.table(paste0("./Input/",i))
  temp_data["Sample"]  <- Desc_file[i,"Sample"]
  temp_data["Context"] <- Desc_file[i, "Context"]
  temp_data["Expression"]  <- Desc_file[i, "Expression"]
  
  data <- rbind(data, temp_data)
}
colnames (data) <- c("Position","mC","C","Percentage",
                     "Sample","Context","Expression")

# Set categories
samples = c("FG1","FG2","FG3","FG4","FG5","FG6","FL")

contexts = c("CpG","CHG","CHH")
contexts_label = c("mCG","mCHG","mCHH")
contexts_mapping = setNames(contexts_label, contexts)

expressions = c("NoExp","2nd","4th","6th","8th","10th")
expressions_label = c("No expression", "2nd decile", "4th decile", 
                      "6th decile", "8th decile", "10th decile")
expressions_mapping = setNames(expressions_label, expressions)
exp_palette = c("#de425b","#e8894a","#dbc667",
                "#79ab62","#488f31","#333333")
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
  select(Position, mC, C, Sample, Context, Expression) %>%
  mutate (Context = contexts_mapping[Context],
          Expression = expressions_mapping[Expression])

data <- na.omit(data)

data <- data %>%
  group_by(Position, Sample, Context, Expression) %>%
  summarise(mC = sum(mC), C = sum(C)) %>%
  mutate(Percentage = mC/(mC + C)*100) %>%
  select(Position, Sample, Context, Expression, Percentage)

data$Sample <- factor(data$Sample, levels = samples)
data$Context <- factor(data$Context, levels = contexts_label)
data$Expression <- factor(data$Expression, levels = expressions_label)

##Plotting
p_sample <- data %>%
  ggplot (aes (x = Position, y = Percentage, color = Expression)) + 
  geom_vline(xintercept = c(20,40),linetype="dashed", linewidth = 0.3, color = "#666666")+
  geom_line (linewidth = 0.4) + 
  theme_classic()+
  scale_color_manual(values = rev(exp_palette))+
  facet_grid(Context~Sample, scales = "free_y") +
  scale_y_continuous(limits = c(0,NA))+
  scale_x_continuous(breaks = c(0,20,40,60),
                     labels = c("-2kb", "TSS", "TES", "+2kb"))+
  labs(x = "", y = "DNA methylation (%)") + 
  theme (axis.text.y = element_text(color = "black", size = 5),
         axis.text.x = element_text(color = "black", size = 5, 
                                    angle = 90, vjust = 0.5, hjust = 0.5),
         axis.line = element_line(color = "black", linewidth = 0.3),
         axis.ticks = element_line(color = "black", linewidth = 0.3),
         axis.title.y = element_text(color = "black", size = 6),
         axis.title.x = element_blank(),
         strip.background = element_blank(),
         strip.text = element_text(color="black", size = 6),
         panel.background = element_blank(),
         panel.spacing = unit(0.5,'mm'),
         legend.background = element_blank(),
         legend.spacing = unit(0,"mm"),
         legend.title = element_blank(),
         legend.key.size = unit(3,"mm"),
         legend.text = element_text(color="black", size = 5))

pdf(file = "Metaplot_exp_by_sample.pdf", width = 7, height = 2)
print(p_sample)
dev.off()
