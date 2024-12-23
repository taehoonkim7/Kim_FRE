############################################################################
### Import libraries
library(ggplot2)
library(tidyverse)
library(ggthemes)
library(ggpubr)
library(cowplot)
library(kneedle)

############################################################################
### Import files
# Import all csv file in the directory
lst_input   <- list.files (path = "./input", pattern = ".bed")
Desc_file <- read.csv ("Sample_Description.csv")

############################################################################
### Pre-process files
#Import data
rownames(Desc_file) <- Desc_file[,1]
data <- data.frame (matrix(ncol = 5, nrow = 0))

#Merge data
for (i in lst_input){
  temp_data <- read.table(paste(getwd(),"/input/",i,sep=""), header = FALSE)
  temp_data["Sample"]  <- Desc_file[i,"Sample"]
  
  data <- rbind(data, temp_data)
}

colnames (data) <- c("chr","start","end","count","Sample")
############################################################################
### CumSumPlot
## Data processing

data <- data %>%
  filter(count!=0) %>%
  arrange(desc(count)) %>%
  group_by(Sample) %>%
  mutate(CumulativeSum = cumsum(count),
         Total = sum(count),
         CumSumPerc = CumulativeSum/Total*100,
         CumRegionPerc = row_number()/n()*100)  %>%
  ungroup()

data$Sample <- factor(data$Sample,
                      levels = c("FG1","FG2","FG3","FG4","FG5","FG6","FL"))
Knee <- data %>%
  group_by(Sample) %>%
  summarise(x = kneedle(CumRegionPerc, CumSumPerc, sensitivity = 2)[1],
            y = kneedle(CumRegionPerc, CumSumPerc, sensitivity = 2)[2])

## Export clusters based on Knee points

for (sample in c("FG1","FG2","FG3","FG4","FG5","FG6","FL")){
  cutoff <- Knee[Knee$Sample==sample,]$y
  
  data_export <- data %>%
    filter(Sample == sample,
           CumSumPerc <= cutoff) %>%
    select(chr,start,end,count,CumSumPerc) %>%
    write.csv(file = paste0("CumSum_",sample,".csv"),row.names = FALSE)
}

##Line graphs
col_pal <- c("#E41817","#F49F1E","#AECA2A","#72BF94",
             "#2A499B","#754696","black")

p <- data %>%
  ggplot(aes(x = CumRegionPerc, y = CumSumPerc,
             group = Sample, color = Sample)) +
  #geom_hline(data = Knee, aes(yintercept = y, color = Sample),
  #           linetype="dashed", linewidth = .3, show.legend = FALSE)+
  geom_segment(data = Knee, 
               aes(x = 0, xend = x, y = y, yend = y, color = Sample),
               linetype="dashed", linewidth = .2, show.legend = FALSE)+
  geom_segment(data = Knee, 
               aes(x = x, xend = x, y = 0, yend = y, color = Sample), 
               linetype="dashed", linewidth = .2, show.legend = FALSE)+
  geom_line(linewidth = .4)+
  geom_point(data = Knee,
             aes(x = x, y = y, color = Sample), size = 0.3)+
  theme_classic()+
  scale_color_manual(values = col_pal)+
  coord_cartesian(clip = 'off') + 
  labs(x = "%24nt-siRNA clusters", y = "%CumulativeSum") +
  theme(axis.line = element_line(color="black", linewidth = 0.2),
        axis.title = element_text(color="black", size = 7),
        axis.text = element_text(color="black", size = 5),
        axis.ticks = element_line(color="black", linewidth = 0.2),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = c(0.7,0.3),
        legend.background = element_blank(),
        legend.title = element_blank(),
        legend.key.size = unit(3,'mm'),
        legend.text = element_text(color="black", size = 5)
  )+
  guides (color = guide_legend(nrow = 4, ncol = 2)) +
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0,1))

pdf(file = paste("CumSumPlot.pdf",sep=""), width = 1.7, height = 1.7)
print(p)
dev.off()

