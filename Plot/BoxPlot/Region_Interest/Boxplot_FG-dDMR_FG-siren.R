############################################################################
### Import libraries
library (ggplot2)
library (tidyverse)
library (ggthemes)
library (ggpubr)
library (cowplot)
library (rstatix)
library (forcats)

############################################################################
### Import files
# Import all csv file in the directory
lst_file    <- list.files (path = "./FG-dDMR_FG-siren", pattern = ".txt")
Sample_desc <- read.csv ("Sample_Description_FG-dDMR_FG-siren.csv")

############################################################################
### Pre-process files
#Import data
rownames(Sample_desc) <- Sample_desc[,1]
data <- data.frame (matrix(ncol = 9, nrow = 0))

#Merge data
for (i in lst_file){
  temp_data <- read.table(paste("./FG-dDMR_FG-siren/",i,sep=""))
  temp_data["V8"]  <- Sample_desc[i,"Sample"]
  temp_data["V9"]  <- Sample_desc[i,"Region"]
  data <- rbind(data, temp_data)
}

colnames (data) <- c("chr","start","end","mC","C","mC_rate","Context","Sample","regionType")

data <- data %>%
  filter(Sample != "FG1-2_S2") %>%
  mutate(Sample = case_when(Sample %in% c("FG1-1_S1")         ~ "FG1",
                            Sample %in% c("FG2-1_S3","FG2-2_S4") ~ "FG2",
                            Sample %in% c("FG3-1_S5","FG3-2_S6") ~ "FG3",
                            Sample %in% c("FG4-1_S7","FG4-2_S8") ~ "FG4",
                            Sample %in% c("FG5-1_S9","FG5-2_S10") ~ "FG5",
                            Sample %in% c("FG6-1_S11","FG6-2_S12") ~ "FG6",
                            Sample %in% c("FL-1_S13","FL-2_S14")   ~ "FL",
                            TRUE ~ Sample
                            )) %>%
  select(chr, start, end, mC, C, Context, Sample, regionType) %>%
  group_by(chr, start, end, Context, Sample, regionType) %>%
  summarise(mC = sum(mC), C = sum(C)) %>%
  filter (mC+C >= 20) %>% #20 reads minimum
  mutate (mC_rate = mC/(mC+C)*100) %>%
  select(chr, start, end, mC, C, mC_rate, Context, Sample, regionType)

data$Sample <- factor(data$Sample,
                      levels = c("FG1","FG2","FG3","FG4","FG5","FG6","FL"))
data$Context <- factor(data$Context,
                       levels = c("CG","CHG","CHH"))
data$regionType <- factor(data$regionType,
                          levels = c("DMR_siren","DMR","siren","shuffle"))

#Filtering NA values
data <- na.omit(data)

############################################################################
###Boxplot of body methylation

fill_pal <- c("#f17272", "#f8c578", "#d1e47b", "#aad9bf",
              "#6b89d7", "#ae88c8", "#666666")

#CG
p_CG <- data %>%
  filter(Context == "CG") %>%
  ggplot(aes(x = Sample, y = mC_rate, fill = Sample))+
  stat_boxplot(geom="errorbar", color="black", linewidth = 0.2)+
  geom_boxplot(outliers = FALSE, color="black", linewidth = 0.2)+
  theme_classic()+
  scale_fill_manual(values = fill_pal)+
  facet_grid(regionType~.) +
  coord_cartesian(clip = 'off') + 
  labs(y = "%Methylation", title = "CG") +
  theme(axis.line = element_line(color="black", linewidth = 0.2),
        axis.title.x = element_blank(),
        axis.title.y = element_text(color="black", size = 6),
        axis.text.y = element_text(color="black", size = 5),
        axis.text.x = element_text(color="black", size = 5,
                                   angle = 90, vjust = 0.5),
        axis.ticks.y = element_line(color="black", linewidth = 0.2),
        axis.ticks.x = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_text(color="black", size = 5),
        plot.title = element_blank(),
  )+
  scale_y_continuous(expand = c(0, 1))

pdf(file = "./boxplot_CG.pdf", width = 1.4, height = 2.2)
print(p_CG)
dev.off()

#CHG
p_CHG <- data %>%
  filter(Context == "CHG") %>%
  ggplot(aes(x = Sample, y = mC_rate, fill = Sample))+
  stat_boxplot(geom="errorbar", color = "black", linewidth = 0.2)+
  geom_boxplot(outliers = FALSE, color="black", linewidth = 0.2)+
  theme_classic()+
  scale_fill_manual(values = fill_pal)+
  facet_grid(regionType~.) +
  coord_cartesian(clip = 'off') + 
  labs(y = "%Methylation", title = "CHG") +
  theme(axis.line = element_line(color="black", linewidth = 0.2),
        axis.title.x = element_blank(),
        axis.title.y = element_text(color="black", size = 6),
        axis.text.y = element_text(color="black", size = 5),
        axis.text.x = element_text(color="black", size = 5,
                                   angle = 90, vjust = 0.5),
        axis.ticks.y = element_line(color="black", linewidth = 0.2),
        axis.ticks.x = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_text(color="black", size = 5),
        plot.title = element_blank(),
  )+
  scale_y_continuous(expand = c(0, 1))

pdf(file = "./boxplot_CHG.pdf", width = 1.4, height = 2.2)
print(p_CHG)
dev.off()

#CHH
p_CHH <- data %>%
  filter(Context == "CHH") %>%
  ggplot(aes(x = Sample, y = mC_rate, fill = Sample))+
  stat_boxplot(geom="errorbar", color = "black", linewidth = 0.2)+
  geom_boxplot(outliers = FALSE, color="black", linewidth = 0.2)+
  theme_classic()+
  scale_fill_manual(values = fill_pal)+
  facet_grid(regionType~.) +
  coord_cartesian(clip = 'off') + 
  labs(y = "%Methylation", title = "CHH") +
  theme(axis.line = element_line(color="black", linewidth = 0.2),
        axis.title.x = element_blank(),
        axis.title.y = element_text(color="black", size = 6),
        axis.text.y = element_text(color="black", size = 5),
        axis.text.x = element_text(color="black", size = 5,
                                   angle = 90, vjust = 0.5),
        axis.ticks.y = element_line(color="black", linewidth = 0.2),
        axis.ticks.x = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_text(color="black", size = 5),
        plot.title = element_blank(),
  )+
  scale_y_continuous(expand = c(0, 1))

pdf(file = "./boxplot_CHH.pdf", width = 1.4, height = 2.2)
print(p_CHH)
dev.off()
