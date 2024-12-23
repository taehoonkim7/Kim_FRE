############################################################################
### Import libraries
library (ggplot2)
library (tidyverse)
library (ggthemes)
library (ggpubr)
library (lemon)
############################################################################
### Import files
# Set directory
# Import all csv file in the directory
lst_file <- list.files (path = "./Input", pattern = ".txt")
description <- read.csv ("Sample_Description.csv")
centromere <- read.csv("Centromere.csv")

############################################################################
### Pre-process files
#Import data
rownames(description) <- description[,1]
data <- data.frame (matrix(ncol = 7, nrow = 0))

#Merge data
for (i in lst_file){
  temp_data <- read.table(paste(getwd(),"/Input/",i,sep=""))
  temp_data["Sample"] <- description[i,"Sample"]
  temp_data["Context"] <- description[i, "Context"]
  
  data <- rbind(data, temp_data)
}

colnames (data) <- c("Chr","Position","mC","C","Percentage",
                     "Sample","Context")

data <- left_join(data, centromere, by = c("Chr" = "chr"))
data <- data %>%
  mutate(Centromere = ifelse(Position >= min & Position < max, "Y", "N"))

data <- data %>%
  mutate(Position = Position - 250000) %>%
  filter(Sample != "FG1-2") %>%
  mutate(Sample = case_when(Sample %in% c("FG1-1") ~ "FG1",
                               Sample %in% c("FG2-1","FG2-2") ~ "FG2",
                               Sample %in% c("FG3-1","FG3-2") ~ "FG3",
                               Sample %in% c("FG4-1","FG4-2") ~ "FG4",
                               Sample %in% c("FG5-1","FG5-2") ~ "FG5",
                               Sample %in% c("FG6-1","FG6-2") ~ "FG6",
                               Sample %in% c("FL-1","FL-2") ~ "FL")) %>%
  mutate(panel = case_when(
    Chr %in% c("Chr1","Chr3","Chr5","Chr7","Chr9","Chr11") ~ "1",
    Chr %in% c("Chr2","Chr4","Chr6","Chr8","Chr10","Chr12") ~ "2")) %>%
  select(Chr, Position, mC, C, Sample, Context, panel, Centromere) %>%
  group_by(Chr, Position, Sample, Context, panel, Centromere) %>%
  summarise(mC = sum(mC), C = sum(C)) %>%
  mutate(Percentage = mC/(mC + C)*100) 

data$Chr <- factor(data$Chr,
                   levels = c("Chr1","Chr2","Chr3","Chr4","Chr5","Chr6",
                              "Chr7","Chr8","Chr9","Chr10","Chr11","Chr12"))
data$Sample <- factor(data$Sample,
                       levels = c("FL","FG6","FG5","FG4","FG3","FG2","FG1"))
data$Context <- factor(data$Context,
                        levels = c("CpG","CHG","CHH"))
data$panel <- factor(data$panel,
                       levels = c("1","2"))

data <- data %>%
  group_by (Chr, Position) %>%
  mutate (Chr_Pos = cur_group_id()) %>%
  ungroup() 

rect_data <- data %>%
  arrange(Chr_Pos) %>%
  group_by(Chr) %>%
  summarise(LastUniqueID = last(Chr_Pos)) %>%
  mutate(ymin = 0, ymax = Inf, xmin = lag(LastUniqueID, default = 1), xmax = LastUniqueID) %>%
  select(Chr, xmin, xmax, ymin, ymax) %>%
  mutate(fill = rep(c("1","2"),6))

centromere = data %>%
  filter(Centromere == "Y") %>%
  select(Chr, Chr_Pos) %>%
  group_by(Chr) %>%
  summarise(Chr_Pos = as.integer(mean(Chr_Pos)))

chr_annotation <- data %>%
  select(Chr, Chr_Pos) %>%
  group_by(Chr) %>%
  summarise(Chr_Pos = as.integer(mean(Chr_Pos)))

##Plotting
col_pal <- c("#E41817","#F49F1E","#AECA2A","#72BF94",
             "#2A499B","#754696","black")
fill_pal <- c("#eeeeee", "#cccccc")
context <- c("CG","CHG","CHH")
names(context) <- c("CpG","CHG","CHH")

#Line graph
p <- ggplot ()+
  geom_rect(data = rect_data, 
            aes(xmin=xmin, xmax=xmax,ymin=ymin, ymax=ymax, fill = fill))+
  geom_vline(xintercept = centromere$Chr_Pos,
             linetype = "dotted", linewidth = .2, color = "#666666")+
  geom_line(data = data, linewidth = .2, 
            aes(x = Chr_Pos, y = Percentage, color = Sample))+
  theme_classic()+
  scale_color_manual(values = rev(col_pal))+
  scale_fill_manual(values = fill_pal)+
  facet_grid(Context~., scales = "free",
             labeller = labeller(Context = context)) +
  geom_text(data = chr_annotation, aes(x = Chr_Pos, y = 0, label = Chr),
            vjust=1.5, size = 2)+
  labs(x = "", y = "%Methylation") +
  theme(panel.spacing = unit (0, "null"),
        axis.title.y = element_text(color="black", size = 7),
        axis.line.x = element_blank(),
        axis.line.y = element_line(color="black", linewidth = 0.2),
        axis.text.x = element_blank(),
        axis.text.y = element_text(color="black", size = 5),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(color="black", linewidth = 0.2),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        strip.background = element_blank(),
        strip.text.y = element_text(color = "black",size = 7)
        )+
  coord_cartesian(clip = 'off')+
  scale_x_continuous(expand = c(0, 0))+
  scale_y_continuous(expand = c(0, 0))

pdf(file = paste("metaplot_chromosome.pdf",sep=""), width = 4.9, height = 2)
print(p)
dev.off()
