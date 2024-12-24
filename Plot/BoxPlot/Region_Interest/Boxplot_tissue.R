############################################################################
### Import libraries
library (ggplot2)
library (tidyverse)
library (dplyr)
library (FSA)
library (multcompView)
library (ggpubr)

############################################################################
### Import files
# Import all csv file in the directory
lst_file    <- list.files (path = "./tissue", pattern = ".txt")
Sample_desc <- read.csv ("Sample_Description_tissue.csv")

############################################################################
### Pre-process files
#Groups
contexts <- c("CG","CHG","CHH")
samples <- c("Seedling","FL","Panicle","Stamen","FG1","FG6","Pistil_Nip",
             "Pistil_rdr2_3","Pistil_rdr2_6","EC","CC", "Embryo", "Endosperm")
pal_sample <- c("#888888","#666666","#ae88c8","#6b89d7","#f17272","#a34c4c","#7d0000",
                "#d1a897","#e4c2b3","#98b357","#7f9d3a","#c4e1b2","#a9d19b")
regions <- c("DMR_siren","DMR","Siren","Shuffle")

#Import data
rownames(Sample_desc) <- Sample_desc[,1]
data <- data.frame (matrix(ncol = 9, nrow = 0))

#Merge data
for (i in lst_file){
  temp_data <- read.table(paste("./tissue/",i,sep=""))
  temp_data["V8"]  <- Sample_desc[i,"Sample"]
  temp_data["V9"]  <- Sample_desc[i,"Region"]
  data <- rbind(data, temp_data)  
}
rm (temp_data)

colnames (data) <- c("chr","start","end","mC","C","mC_rate",
                     "Context","Sample","regionType")

data <- data %>%
  mutate(Sample = case_when(Sample %in% c("FG1-1")         ~ "FG1",
                            Sample %in% c("FG6-1","FG6-2") ~ "FG6",
                            Sample %in% c("FL-1","FL-2")   ~ "FL",
                            Sample %in% c("Panicle-1","Panicle-2") ~ "Panicle",
                            Sample %in% c("Seedling-1","Seedling-2") ~ "Seedling",
                            Sample %in% c("Pistil") ~ "Pistil_Nip",
                            Sample %in% c("rdr2-3") ~ "Pistil_rdr2_3",
                            Sample %in% c("rdr2-6") ~ "Pistil_rdr2_6",
                            TRUE ~ Sample
                            )) %>%
  filter(!Sample %in% c("FG2","FG3","FG4","FG5")) %>%
  dplyr::select(chr,start,end, mC, C, Context, Sample, regionType) %>%
  group_by(chr,start,end, Context, Sample, regionType) %>%
  summarise(mC = sum(mC), C = sum(C)) %>%
  mutate (mC_rate = mC/(mC+C)*100) %>%
  dplyr::select(chr,start,end, mC, C, mC_rate, Context, Sample, regionType) %>%
  mutate (mC_rate = ifelse(mC+C >= 20, mC_rate, NA)) %>% #20 reads minimum 
  ungroup()
  
data$Sample <- factor(data$Sample, levels = samples)
data$Context <- factor(data$Context, levels = contexts)
data$regionType <- factor(data$regionType, levels = regions)

#Filtering NA values
data <- na.omit(data)

############################################################################
### Kruskal-Wallis / Dunn's test - comparing samples 
data_dunn <- data.frame (matrix(ncol = 6, nrow = 0))
data_dunn_letter <- data.frame (matrix(ncol = 4, nrow = 0))
data_kruskal <- data.frame (matrix(ncol = 7, nrow = 0))

for (context in contexts){
  for (region in regions){
    data_to_test <- data %>%
      filter (Context == context, regionType == region) %>%
      dplyr::select (Sample, mC_rate)
    
    temp_kruskal <- kruskal.test (mC_rate ~ Sample, data = data_to_test)
    temp_kruskal$Context = context
    temp_kruskal$Region = region
    data_kruskal <- rbind (data_kruskal, temp_kruskal)
    
    if (temp_kruskal$p.value < 0.00001){
      temp_dunn <- dunnTest (mC_rate ~ Sample, data = data_to_test, method = "bh")$res
      temp_dunn$Context = context
      temp_dunn$Region = region
      temp_dunn <- temp_dunn[,c("Context","Region","Comparison","Z","P.unadj","P.adj")]
      data_dunn <- rbind(data_dunn, temp_dunn)
      
      dunn_p <- temp_dunn$P.adj
      names(dunn_p) <- gsub(" ", "", temp_dunn$Comparison)
      
      temp_dunn_letter <- data.frame(multcompLetters(dunn_p)["Letters"])
      temp_dunn_letter$Sample <- rownames(temp_dunn_letter)
      rownames(temp_dunn_letter) <- NULL
      temp_dunn_letter$Context = context
      temp_dunn_letter$Region = region
      temp_dunn_letter <- temp_dunn_letter[,c("Context","Region","Sample","Letters")]
      data_dunn_letter <- rbind(data_dunn_letter, temp_dunn_letter)
    }
  }
}

write.csv(data_kruskal, file = "Kruskal.csv")
write.csv(data_dunn, file = "Dunn_Individual.csv")
write.csv(data_dunn_letter, file = "Dunn_letter.csv")

############################################################################
###Boxplot of body methylation

p_CG_legend <- data %>%
  filter(Context == "CG") %>%
  ggplot (aes(x = regionType, y = mC_rate, fill = Sample))+
  stat_boxplot(geom="errorbar")+
  geom_boxplot(outlier.shape = NA)+
  theme_classic()+
  scale_fill_manual(values = pal_sample)+
  labs(title = "CpG", x = "", y = "DNA methylation (%)")+
  ylim (c(0,100))
p_CG <- p_CG_legend+
  theme(legend.position = "none")
p_legend <- cowplot::get_legend(p_CG_legend)

p_CHG <- data %>%
  filter(Context == "CHG") %>%
  ggplot (aes(x = regionType, y = mC_rate, fill = Sample))+
  stat_boxplot(geom="errorbar")+
  geom_boxplot(outlier.shape = NA)+
  theme_classic()+
  scale_fill_manual(values = pal_sample)+
  labs(title = "CHG",x = "", y = "DNA methylation (%)")+
  ylim (c(0,100))+
  theme(legend.position = "none")

p_CHH <- data %>%
  filter(Context == "CHH") %>%
  ggplot (aes(x = regionType, y = mC_rate, fill = Sample))+
  stat_boxplot(geom="errorbar")+
  geom_boxplot(outlier.shape = NA)+
  theme_classic()+
  scale_fill_manual(values = pal_sample)+
  labs(title = "CHH",x = "", y = "DNA methylation (%)")+
  ylim (c(0,100))+
  theme(legend.position = "none")

p <- ggarrange(p_CG  , NULL,
               p_CHG , p_legend,
               p_CHH , NULL,
               ncol = 2, nrow = 3,
               widths = c(10, 2))

pdf(file = "./boxplot_tissue.pdf", width = 12, height = 6)
print(p)
dev.off()
  
