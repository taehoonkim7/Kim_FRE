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

## Export table
table_GeneTE <- data_GeneTE %>%
  select (-c(mC,C)) %>%
  mutate (Position = case_when(
    Position < 19 ~ paste0("Upstream (-",as.character(2000-(Position-0.5)*100),
                           " - -",as.character(2000-(Position+0.5)*100),")"),
    Position == 19.5 ~ "Upstream (-100 - TSS)",
    Position == 40.5 ~ "Downstream (TES - +100)",
    Position > 41 ~ paste0("Downstream (+",as.character((Position-40.5)*100),
                           " - +",as.character((Position-39.5)*100),")"),
    .default = paste0("Body_Proportional_Bin_",as.integer(Position-19.5)),
  )) %>%
  arrange(desc(Sample)) %>%
  arrange(Context) %>%
  arrange(GeneTE) %>%
  mutate (Context = case_when(
    Context == "CpG" ~ "mCG",
    Context == "CHG" ~ "mCHG",
    Context == "CHH" ~ "mCHH"
  )) %>%
  setNames(c("Position","Sample","Context","Element","%mC")) %>%
  select("Element","Context","Sample","Position","%mC")
write_csv(table_GeneTE, file = 'data_GeneTE.csv')


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
        #legend.title = element_blank(),
        #legend.position = "right",
        #legend.text = element_text(color="black",size=8),
        #legend.key.size = unit(1, "lines"),
        legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_blank()
  )+
  scale_x_continuous(breaks = c(0,20,40,60),
                     labels = c("-2kb", "", "", "+2kb"))+
  scale_y_continuous(expand = c(0, 2))

#legend <- cowplot::get_legend(p_TE)
#p_TE <- p_TE + 
#  theme (legend.position = "none")
p_combined <- ggarrange(p_gene, p_TE,
                        ncol = 2, nrow = 1,
               widths = c(3, 3))

pdf(file = paste("metaplot_geneTE.pdf",sep=""), width = 3.5, height = 2.5)
print(p_combined)
dev.off()

##############################################################################
##GeneTE
CpG_Gene_legend <- data_GeneTE %>%
  filter(Context == "CpG") %>%
  filter(GeneTE == "Gene") %>%
  ggplot (aes(x = Position, y = Percentage, color = Sample))+
  geom_line(linewidth = 1)+
  geom_vline(xintercept = c(20,40),linetype="dashed")+
  theme_classic()+
  scale_color_manual(values = vivid_palette)+
  scale_x_continuous(breaks = c(0,20,40,60),
                     labels = c("-2kb", "TSS", "TES", "+2kb"))+
  labs(x = "", y = "")+
  ylim (c(0,60))
CpG_Gene <- CpG_Gene_legend+
  theme(legend.position = "none")
legend <- cowplot::get_legend(CpG_Gene_legend)

CHG_Gene <- data_GeneTE %>%
  filter(Context == "CHG") %>%
  filter(GeneTE == "Gene") %>%
  ggplot (aes(x = Position, y = Percentage, color = Sample))+
  geom_line(linewidth = 1)+
  geom_vline(xintercept = c(20,40),linetype="dashed")+
  theme_classic()+
  scale_color_manual(values = vivid_palette)+
  scale_x_continuous(breaks = c(0,20,40,60),
                     labels = c("-2kb", "TSS", "TES", "+2kb"))+
  labs(x = "", y = "")+
  ylim (c(0,40))+
  theme(legend.position = "none")

CHH_Gene <- data_GeneTE %>%
  filter(Context == "CHH") %>%
  filter(GeneTE == "Gene") %>%
  ggplot (aes(x = Position, y = Percentage, color = Sample))+
  geom_line(linewidth = 1)+
  geom_vline(xintercept = c(20,40),linetype="dashed")+
  theme_classic()+
  scale_color_manual(values = vivid_palette)+
  scale_x_continuous(breaks = c(0,20,40,60),
                     labels = c("-2kb", "TSS", "TES", "+2kb"))+
  labs(x = "", y = "")+
  ylim (c(0,15))+
  theme(legend.position = "none")

CpG_TE <- data_GeneTE %>%
  filter(Context == "CpG") %>%
  filter(GeneTE == "TE") %>%
  ggplot (aes(x = Position, y = Percentage, color = Sample))+
  geom_line(linewidth = 1)+
  geom_vline(xintercept = c(20,40),linetype="dashed")+
  theme_classic()+
  scale_color_manual(values = vivid_palette)+
  scale_x_continuous(breaks = c(0,20,40,60),
                     labels = c("-2kb", "Start", "End", "+2kb"))+
  labs(x = "", y = "")+
  ylim (c(40,90))+
  theme(legend.position = "none")

CHG_TE <- data_GeneTE %>%
  filter(Context == "CHG") %>%
  filter(GeneTE == "TE") %>%
  ggplot (aes(x = Position, y = Percentage, color = Sample))+
  geom_line(linewidth = 1)+
  geom_vline(xintercept = c(20,40),linetype="dashed")+
  theme_classic()+
  scale_color_manual(values = vivid_palette)+
  scale_x_continuous(breaks = c(0,20,40,60),
                     labels = c("-2kb", "Start", "End", "+2kb"))+
  labs(x = "", y = "")+
  ylim (c(20,70))+
  theme(legend.position = "none")

CHH_TE <- data_GeneTE %>%
  filter(Context == "CHH") %>%
  filter(GeneTE == "TE") %>%
  ggplot (aes(x = Position, y = Percentage, color = Sample))+
  geom_line(linewidth = 1)+
  geom_vline(xintercept = c(20,40),linetype="dashed")+
  theme_classic()+
  scale_color_manual(values = vivid_palette)+
  scale_x_continuous(breaks = c(0,20,40,60),
                     labels = c("-2kb", "Start", "End", "+2kb"))+
  labs(x = "", y = "")+
  scale_y_continuous(breaks = c(6,8,10,12,14), limits = c(6,15))+ 
  theme(legend.position = "none")

p <- ggarrange(CpG_Gene, CHG_Gene, CHH_Gene, NULL,
          CpG_TE  , CHG_TE  , CHH_TE  , legend,
          ncol = 4, nrow = 2,
          widths = c(2, 2, 2, 1)) #5.75x9.0

pdf(file = paste("metaplot_geneTE.pdf",sep=""), width = 7, height = 2)
print(p)
dev.off()

##TE families
title_Copia   = "LTR/Ty1-copia"
title_Gypsy   = "LTR/Ty3-gypsy"
title_TRIM    = "LTR/TRIM"
title_LTRun   = "LTR/unknown"
title_CACTA   = "TIR/CACTA"
title_Mut     = "TIR/Mutator"
title_PIF     = "TIR/PIF_Harbinger"
title_Tc1     = "TIR/Tc1_Mariner"
title_hAT     = "TIR/hAT"
title_TIRun   = "TIR/unknown"
title_Cen     = "Centromeric repeat"
title_LINE    = "nonLTR/LINE element"
title_SINE    = "nonLTR/SINE element"
title_Nonun   = "nonLTR/unknown"
title_HEL     = "nonTIR/helitron"
title_Rep     = "Repeat region"
title_Sat     = "Satellite DNA"

n_Copia = "n=11547"
n_Gypsy = "n=41002"
n_TRIM  = "n=2829"
n_LTRun = "n=1359"
n_CACTA = "n=17827"
n_Mut   = "n=47744"
n_PIF   = "n=42068"
n_Tc1   = "n=46521"
n_hAT   = "n=14413"
n_TIRun = "n=6522"
n_Cen   = "n=454"
n_LINE  = "n=7292"
n_SINE  = "n=6308"
n_Nonun = "n=275"
n_HEL   = "n=47752"
n_Rep   = "n=1279"
n_Sat   = "n=31"

#CpG methylation
CpG_Copia_legend <- data_TE %>%
  filter(Context == "CpG") %>%
  filter(Family == "Copia") %>%
  ggplot (aes(x = Position, y = Percentage, color = Sample))+
  geom_line(linewidth = 1)+
  geom_vline(xintercept = c(20,40),linetype="dashed")+
  theme_classic()+
  theme(plot.subtitle = element_text(color = "#999999"))+
  scale_color_manual(values = vivid_palette)+
  scale_x_continuous(breaks = c(0,20,40,60),
                     labels = c("-2kb", "TSS", "TES", "+2kb"))+
  labs(title = title_Copia, subtitle = n_Copia, x = "", y = "DNA methylation (%)")+
  ylim (c(50,100))
CpG_Copia <- CpG_Copia_legend+
  theme(legend.position = "none")
TE_legend <- cowplot::get_legend(CpG_Copia_legend)

CpG_Gypsy <- data_TE %>%
  filter(Context == "CpG") %>%
  filter(Family == "Gypsy") %>%
  ggplot (aes(x = Position, y = Percentage, color = Sample))+
  geom_line(linewidth = 1)+
  geom_vline(xintercept = c(20,40),linetype="dashed")+
  theme_classic()+
  theme(plot.subtitle = element_text(color = "#999999"))+
  scale_color_manual(values = vivid_palette)+
  scale_x_continuous(breaks = c(0,20,40,60),
                     labels = c("-2kb", "TSS", "TES", "+2kb"))+
  labs(title = title_Gypsy, subtitle = n_Gypsy, x = "", y = "")+
  ylim (c(60,100))+ theme(legend.position = "none")

CpG_TRIM <- data_TE %>%
  filter(Context == "CpG") %>%
  filter(Family == "TRIM") %>%
  ggplot (aes(x = Position, y = Percentage, color = Sample))+
  geom_line(linewidth = 1)+
  geom_vline(xintercept = c(20,40),linetype="dashed")+
  theme_classic()+
  theme(plot.subtitle = element_text(color = "#999999"))+
  scale_color_manual(values = vivid_palette)+
  scale_x_continuous(breaks = c(0,20,40,60),
                     labels = c("-2kb", "TSS", "TES", "+2kb"))+
  labs(title = title_TRIM, subtitle = n_TRIM, x = "", y = "")+
  ylim (c(50,100))+ theme(legend.position = "none")

CpG_LTRun <- data_TE %>%
  filter(Context == "CpG") %>%
  filter(Family == "LTRun") %>%
  ggplot (aes(x = Position, y = Percentage, color = Sample))+
  geom_line(linewidth = 1)+
  geom_vline(xintercept = c(20,40),linetype="dashed")+
  theme_classic()+
  theme(plot.subtitle = element_text(color = "#999999"))+
  scale_color_manual(values = vivid_palette)+
  scale_x_continuous(breaks = c(0,20,40,60),
                     labels = c("-2kb", "TSS", "TES", "+2kb"))+
  labs(title = title_LTRun, subtitle = n_LTRun, x = "", y = "")+
  ylim (c(50,90))+ theme(legend.position = "none")

CpG_CACTA <- data_TE %>%
  filter(Context == "CpG") %>%
  filter(Family == "CACTA") %>%
  ggplot (aes(x = Position, y = Percentage, color = Sample))+
  geom_line(linewidth = 1)+
  geom_vline(xintercept = c(20,40),linetype="dashed")+
  theme_classic()+
  theme(plot.subtitle = element_text(color = "#999999"))+
  scale_color_manual(values = vivid_palette)+
  scale_x_continuous(breaks = c(0,20,40,60),
                     labels = c("-2kb", "TSS", "TES", "+2kb"))+
  labs(title = title_CACTA, subtitle = n_CACTA, x = "", y = "")+
  ylim (c(50,100))+ theme(legend.position = "none")

CpG_Mut <- data_TE %>%
  filter(Context == "CpG") %>%
  filter(Family == "Mut") %>%
  ggplot (aes(x = Position, y = Percentage, color = Sample))+
  geom_line(linewidth = 1)+
  geom_vline(xintercept = c(20,40),linetype="dashed")+
  theme_classic()+
  theme(plot.subtitle = element_text(color = "#999999"))+
  scale_color_manual(values = vivid_palette)+
  scale_x_continuous(breaks = c(0,20,40,60),
                     labels = c("-2kb", "TSS", "TES", "+2kb"))+
  labs(title = title_Mut, subtitle = n_Mut, x = "", y = "")+
  ylim (c(30,90))+ theme(legend.position = "none")

CpG_PIF <- data_TE %>%
  filter(Context == "CpG") %>%
  filter(Family == "PIF") %>%
  ggplot (aes(x = Position, y = Percentage, color = Sample))+
  geom_line(linewidth = 1)+
  geom_vline(xintercept = c(20,40),linetype="dashed")+
  theme_classic()+
  theme(plot.subtitle = element_text(color = "#999999"))+
  scale_color_manual(values = vivid_palette)+
  scale_x_continuous(breaks = c(0,20,40,60),
                     labels = c("-2kb", "TSS", "TES", "+2kb"))+
  labs(title = title_PIF, subtitle = n_PIF, x = "", y = "DNA methylation (%)")+
  ylim (c(30,100))+ theme(legend.position = "none")

CpG_Tc1 <- data_TE %>%
  filter(Context == "CpG") %>%
  filter(Family == "Tc1") %>%
  ggplot (aes(x = Position, y = Percentage, color = Sample))+
  geom_line(linewidth = 1)+
  geom_vline(xintercept = c(20,40),linetype="dashed")+
  theme_classic()+
  theme(plot.subtitle = element_text(color = "#999999"))+
  scale_color_manual(values = vivid_palette)+
  scale_x_continuous(breaks = c(0,20,40,60),
                     labels = c("-2kb", "TSS", "TES", "+2kb"))+
  labs(title = title_Tc1, subtitle = n_Tc1, x = "", y = "")+
  ylim (c(30,90))+ theme(legend.position = "none")

CpG_hAT <- data_TE %>%
  filter(Context == "CpG") %>%
  filter(Family == "hAT") %>%
  ggplot (aes(x = Position, y = Percentage, color = Sample))+
  geom_line(linewidth = 1)+
  geom_vline(xintercept = c(20,40),linetype="dashed")+
  theme_classic()+
  theme(plot.subtitle = element_text(color = "#999999"))+
  scale_color_manual(values = vivid_palette)+
  scale_x_continuous(breaks = c(0,20,40,60),
                     labels = c("-2kb", "TSS", "TES", "+2kb"))+
  labs(title = title_hAT, subtitle = n_hAT, x = "", y = "")+
  ylim (c(40,80))+ theme(legend.position = "none")

CpG_TIRun <- data_TE %>%
  filter(Context == "CpG") %>%
  filter(Family == "TIRun") %>%
  ggplot (aes(x = Position, y = Percentage, color = Sample))+
  geom_line(linewidth = 1)+
  geom_vline(xintercept = c(20,40),linetype="dashed")+
  theme_classic()+
  theme(plot.subtitle = element_text(color = "#999999"))+
  scale_color_manual(values = vivid_palette)+
  scale_x_continuous(breaks = c(0,20,40,60),
                     labels = c("-2kb", "TSS", "TES", "+2kb"))+
  labs(title = title_TIRun, subtitle = n_TIRun, x = "", y = "")+
  ylim (c(30,90))+ theme(legend.position = "none")

CpG_Cen <- data_TE %>%
  filter(Context == "CpG") %>%
  filter(Family == "Cen") %>%
  ggplot (aes(x = Position, y = Percentage, color = Sample))+
  geom_line(linewidth = 1)+
  geom_vline(xintercept = c(20,40),linetype="dashed")+
  theme_classic()+
  theme(plot.subtitle = element_text(color = "#999999"))+
  scale_color_manual(values = vivid_palette)+
  scale_x_continuous(breaks = c(0,20,40,60),
                     labels = c("-2kb", "TSS", "TES", "+2kb"))+
  labs(title = title_Cen, subtitle = n_Cen, x = "", y = "")+
  ylim (c(85,100))+ theme(legend.position = "none")

CpG_LINE <- data_TE %>%
  filter(Context == "CpG") %>%
  filter(Family == "LINE") %>%
  ggplot (aes(x = Position, y = Percentage, color = Sample))+
  geom_line(linewidth = 1)+
  geom_vline(xintercept = c(20,40),linetype="dashed")+
  theme_classic()+
  theme(plot.subtitle = element_text(color = "#999999"))+
  scale_color_manual(values = vivid_palette)+
  scale_x_continuous(breaks = c(0,20,40,60),
                     labels = c("-2kb", "TSS", "TES", "+2kb"))+
  labs(title = title_LINE, subtitle = n_LINE, x = "", y = "")+
  ylim (c(40,90))+ theme(legend.position = "none")

CpG_SINE <- data_TE %>%
  filter(Context == "CpG") %>%
  filter(Family == "SINE") %>%
  ggplot (aes(x = Position, y = Percentage, color = Sample))+
  geom_line(linewidth = 1)+
  geom_vline(xintercept = c(20,40),linetype="dashed")+
  theme_classic()+
  theme(plot.subtitle = element_text(color = "#999999"))+
  scale_color_manual(values = vivid_palette)+
  scale_x_continuous(breaks = c(0,20,40,60),
                     labels = c("-2kb", "TSS", "TES", "+2kb"))+
  labs(title = title_SINE, subtitle = n_SINE, x = "", y = "")+
  ylim (c(30,100))+ theme(legend.position = "none")

CpG_Nonun <- data_TE %>%
  filter(Context == "CpG") %>%
  filter(Family == "Nonun") %>%
  ggplot (aes(x = Position, y = Percentage, color = Sample))+
  geom_line(linewidth = 1)+
  geom_vline(xintercept = c(20,40),linetype="dashed")+
  theme_classic()+
  theme(plot.subtitle = element_text(color = "#999999"))+
  scale_color_manual(values = vivid_palette)+
  scale_x_continuous(breaks = c(0,20,40,60),
                     labels = c("-2kb", "TSS", "TES", "+2kb"))+
  labs(title = title_Nonun, subtitle = n_Nonun, x = "", y = "DNA methylation (%)")+
  ylim (c(50,100))+ theme(legend.position = "none")

CpG_HEL <- data_TE %>%
  filter(Context == "CpG") %>%
  filter(Family == "HEL") %>%
  ggplot (aes(x = Position, y = Percentage, color = Sample))+
  geom_line(linewidth = 1)+
  geom_vline(xintercept = c(20,40),linetype="dashed")+
  theme_classic()+
  theme(plot.subtitle = element_text(color = "#999999"))+
  scale_color_manual(values = vivid_palette)+
  scale_x_continuous(breaks = c(0,20,40,60),
                     labels = c("-2kb", "TSS", "TES", "+2kb"))+
  labs(title = title_HEL, subtitle = n_HEL, x = "", y = "")+
  ylim (c(35,75))+ theme(legend.position = "none")

CpG_Rep <- data_TE %>%
  filter(Context == "CpG") %>%
  filter(Family == "Rep") %>%
  ggplot (aes(x = Position, y = Percentage, color = Sample))+
  geom_line(linewidth = 1)+
  geom_vline(xintercept = c(20,40),linetype="dashed")+
  theme_classic()+
  theme(plot.subtitle = element_text(color = "#999999"))+
  scale_color_manual(values = vivid_palette)+
  scale_x_continuous(breaks = c(0,20,40,60),
                     labels = c("-2kb", "TSS", "TES", "+2kb"))+
  labs(title = title_Rep, subtitle = n_Rep, x = "", y = "")+
  ylim (c(20,60))+ theme(legend.position = "none")

CpG_Sat <- data_TE %>%
  filter(Context == "CpG") %>%
  filter(Family == "Sat") %>%
  ggplot (aes(x = Position, y = Percentage, color = Sample))+
  geom_line(linewidth = 1)+
  geom_vline(xintercept = c(20,40),linetype="dashed")+
  theme_classic()+
  theme(plot.subtitle = element_text(color = "#999999"))+
  scale_color_manual(values = vivid_palette)+
  scale_x_continuous(breaks = c(0,20,40,60),
                     labels = c("-2kb", "TSS", "TES", "+2kb"))+
  labs(title = title_Sat, subtitle = n_Sat, x = "", y = "")+
  ylim (c(60,100))+ theme(legend.position = "none")

ggarrange(CpG_Copia , CpG_Gypsy , CpG_TRIM  , CpG_LTRun , CpG_CACTA, CpG_Mut,
          CpG_PIF   , CpG_Tc1   , CpG_hAT   , CpG_TIRun , CpG_LINE, CpG_SINE,
          CpG_Nonun , CpG_HEL   , CpG_Cen   , CpG_Rep   , CpG_Sat , TE_legend,
          ncol = 6, nrow = 3,
          widths  = c(2, 2, 2, 2, 2, 2)) #8.0x15.0

#CHG methylation
CHG_Copia <- data_TE %>%
  filter(Context == "CHG") %>%
  filter(Family == "Copia") %>%
  ggplot (aes(x = Position, y = Percentage, color = Sample))+
  geom_line(linewidth = 1)+
  geom_vline(xintercept = c(20,40),linetype="dashed")+
  theme_classic()+
  theme(plot.subtitle = element_text(color = "#999999"))+
  scale_color_manual(values = vivid_palette)+
  scale_x_continuous(breaks = c(0,20,40,60),
                     labels = c("-2kb", "TSS", "TES", "+2kb"))+
  labs(title = title_Copia, subtitle = n_Copia, x = "", y = "DNA methylation (%)")+
  ylim (c(30,80))+ theme(legend.position = "none")

CHG_Gypsy <- data_TE %>%
  filter(Context == "CHG") %>%
  filter(Family == "Gypsy") %>%
  ggplot (aes(x = Position, y = Percentage, color = Sample))+
  geom_line(linewidth = 1)+
  geom_vline(xintercept = c(20,40),linetype="dashed")+
  theme_classic()+
  theme(plot.subtitle = element_text(color = "#999999"))+
  scale_color_manual(values = vivid_palette)+
  scale_x_continuous(breaks = c(0,20,40,60),
                     labels = c("-2kb", "TSS", "TES", "+2kb"))+
  labs(title = title_Gypsy, subtitle = n_Gypsy, x = "", y = "")+
  ylim (c(40,90))+ theme(legend.position = "none")

CHG_TRIM <- data_TE %>%
  filter(Context == "CHG") %>%
  filter(Family == "TRIM") %>%
  ggplot (aes(x = Position, y = Percentage, color = Sample))+
  geom_line(linewidth = 1)+
  geom_vline(xintercept = c(20,40),linetype="dashed")+
  theme_classic()+
  theme(plot.subtitle = element_text(color = "#999999"))+
  scale_color_manual(values = vivid_palette)+
  scale_x_continuous(breaks = c(0,20,40,60),
                     labels = c("-2kb", "TSS", "TES", "+2kb"))+
  labs(title = title_TRIM, subtitle = n_TRIM, x = "", y = "")+
  ylim (c(25,80))+ theme(legend.position = "none")

CHG_LTRun <- data_TE %>%
  filter(Context == "CHG") %>%
  filter(Family == "LTRun") %>%
  ggplot (aes(x = Position, y = Percentage, color = Sample))+
  geom_line(linewidth = 1)+
  geom_vline(xintercept = c(20,40),linetype="dashed")+
  theme_classic()+
  theme(plot.subtitle = element_text(color = "#999999"))+
  scale_color_manual(values = vivid_palette)+
  scale_x_continuous(breaks = c(0,20,40,60),
                     labels = c("-2kb", "TSS", "TES", "+2kb"))+
  labs(title = title_LTRun, subtitle = n_LTRun, x = "", y = "")+
  ylim (c(25,65))+ theme(legend.position = "none")

CHG_CACTA <- data_TE %>%
  filter(Context == "CHG") %>%
  filter(Family == "CACTA") %>%
  ggplot (aes(x = Position, y = Percentage, color = Sample))+
  geom_line(linewidth = 1)+
  geom_vline(xintercept = c(20,40),linetype="dashed")+
  theme_classic()+
  theme(plot.subtitle = element_text(color = "#999999"))+
  scale_color_manual(values = vivid_palette)+
  scale_x_continuous(breaks = c(0,20,40,60),
                     labels = c("-2kb", "TSS", "TES", "+2kb"))+
  labs(title = title_CACTA, subtitle = n_CACTA, x = "", y = "")+
  ylim (c(30,80))+ theme(legend.position = "none")

CHG_Mut <- data_TE %>%
  filter(Context == "CHG") %>%
  filter(Family == "Mut") %>%
  ggplot (aes(x = Position, y = Percentage, color = Sample))+
  geom_line(linewidth = 1)+
  geom_vline(xintercept = c(20,40),linetype="dashed")+
  theme_classic()+
  theme(plot.subtitle = element_text(color = "#999999"))+
  scale_color_manual(values = vivid_palette)+
  scale_x_continuous(breaks = c(0,20,40,60),
                     labels = c("-2kb", "TSS", "TES", "+2kb"))+
  labs(title = title_Mut, subtitle = n_Mut, x = "", y = "")+
  ylim (c(20,70))+ theme(legend.position = "none")

CHG_PIF <- data_TE %>%
  filter(Context == "CHG") %>%
  filter(Family == "PIF") %>%
  ggplot (aes(x = Position, y = Percentage, color = Sample))+
  geom_line(linewidth = 1)+
  geom_vline(xintercept = c(20,40),linetype="dashed")+
  theme_classic()+
  theme(plot.subtitle = element_text(color = "#999999"))+
  scale_color_manual(values = vivid_palette)+
  scale_x_continuous(breaks = c(0,20,40,60),
                     labels = c("-2kb", "TSS", "TES", "+2kb"))+
  labs(title = title_PIF, subtitle = n_PIF, x = "", y = "DNA methylation (%)")+
  ylim (c(10,80))+ theme(legend.position = "none")

CHG_Tc1 <- data_TE %>%
  filter(Context == "CHG") %>%
  filter(Family == "Tc1") %>%
  ggplot (aes(x = Position, y = Percentage, color = Sample))+
  geom_line(linewidth = 1)+
  geom_vline(xintercept = c(20,40),linetype="dashed")+
  theme_classic()+
  theme(plot.subtitle = element_text(color = "#999999"))+
  scale_color_manual(values = vivid_palette)+
  scale_x_continuous(breaks = c(0,20,40,60),
                     labels = c("-2kb", "TSS", "TES", "+2kb"))+
  labs(title = title_Tc1, subtitle = n_Tc1, x = "", y = "")+
  ylim (c(10,60))+ theme(legend.position = "none")

CHG_hAT <- data_TE %>%
  filter(Context == "CHG") %>%
  filter(Family == "hAT") %>%
  ggplot (aes(x = Position, y = Percentage, color = Sample))+
  geom_line(linewidth = 1)+
  geom_vline(xintercept = c(20,40),linetype="dashed")+
  theme_classic()+
  theme(plot.subtitle = element_text(color = "#999999"))+
  scale_color_manual(values = vivid_palette)+
  scale_x_continuous(breaks = c(0,20,40,60),
                     labels = c("-2kb", "TSS", "TES", "+2kb"))+
  labs(title = title_hAT, subtitle = n_hAT, x = "", y = "")+
  ylim (c(20,55))+ theme(legend.position = "none")

CHG_TIRun <- data_TE %>%
  filter(Context == "CHG") %>%
  filter(Family == "TIRun") %>%
  ggplot (aes(x = Position, y = Percentage, color = Sample))+
  geom_line(linewidth = 1)+
  geom_vline(xintercept = c(20,40),linetype="dashed")+
  theme_classic()+
  theme(plot.subtitle = element_text(color = "#999999"))+
  scale_color_manual(values = vivid_palette)+
  scale_x_continuous(breaks = c(0,20,40,60),
                     labels = c("-2kb", "TSS", "TES", "+2kb"))+
  labs(title = title_TIRun, subtitle = n_TIRun, x = "", y = "")+
  ylim (c(15,60))+ theme(legend.position = "none")

CHG_Cen <- data_TE %>%
  filter(Context == "CHG") %>%
  filter(Family == "Cen") %>%
  ggplot (aes(x = Position, y = Percentage, color = Sample))+
  geom_line(linewidth = 1)+
  geom_vline(xintercept = c(20,40),linetype="dashed")+
  theme_classic()+
  theme(plot.subtitle = element_text(color = "#999999"))+
  scale_color_manual(values = vivid_palette)+
  scale_x_continuous(breaks = c(0,20,40,60),
                     labels = c("-2kb", "TSS", "TES", "+2kb"))+
  labs(title = title_Cen, subtitle = n_Cen, x = "", y = "")+
  ylim (c(40,70))+ theme(legend.position = "none")

CHG_LINE <- data_TE %>%
  filter(Context == "CHG") %>%
  filter(Family == "LINE") %>%
  ggplot (aes(x = Position, y = Percentage, color = Sample))+
  geom_line(linewidth = 1)+
  geom_vline(xintercept = c(20,40),linetype="dashed")+
  theme_classic()+
  theme(plot.subtitle = element_text(color = "#999999"))+
  scale_color_manual(values = vivid_palette)+
  scale_x_continuous(breaks = c(0,20,40,60),
                     labels = c("-2kb", "TSS", "TES", "+2kb"))+
  labs(title = title_LINE, subtitle = n_LINE, x = "", y = "")+
  ylim (c(15,70))+ theme(legend.position = "none")

CHG_SINE <- data_TE %>%
  filter(Context == "CHG") %>%
  filter(Family == "SINE") %>%
  ggplot (aes(x = Position, y = Percentage, color = Sample))+
  geom_line(linewidth = 1)+
  geom_vline(xintercept = c(20,40),linetype="dashed")+
  theme_classic()+
  theme(plot.subtitle = element_text(color = "#999999"))+
  scale_color_manual(values = vivid_palette)+
  scale_x_continuous(breaks = c(0,20,40,60),
                     labels = c("-2kb", "TSS", "TES", "+2kb"))+
  labs(title = title_SINE, subtitle = n_SINE, x = "", y = "")+
  ylim (c(10,70))+ theme(legend.position = "none")

CHG_Nonun <- data_TE %>%
  filter(Context == "CHG") %>%
  filter(Family == "Nonun") %>%
  ggplot (aes(x = Position, y = Percentage, color = Sample))+
  geom_line(linewidth = 1)+
  geom_vline(xintercept = c(20,40),linetype="dashed")+
  theme_classic()+
  theme(plot.subtitle = element_text(color = "#999999"))+
  scale_color_manual(values = vivid_palette)+
  scale_x_continuous(breaks = c(0,20,40,60),
                     labels = c("-2kb", "TSS", "TES", "+2kb"))+
  labs(title = title_Nonun, subtitle = n_Nonun, x = "", y = "DNA methylation (%)")+
  ylim (c(40,100))+ theme(legend.position = "none")

CHG_HEL <- data_TE %>%
  filter(Context == "CHG") %>%
  filter(Family == "HEL") %>%
  ggplot (aes(x = Position, y = Percentage, color = Sample))+
  geom_line(linewidth = 1)+
  geom_vline(xintercept = c(20,40),linetype="dashed")+
  theme_classic()+
  theme(plot.subtitle = element_text(color = "#999999"))+
  scale_color_manual(values = vivid_palette)+
  scale_x_continuous(breaks = c(0,20,40,60),
                     labels = c("-2kb", "TSS", "TES", "+2kb"))+
  labs(title = title_HEL, subtitle = n_HEL, x = "", y = "")+
  ylim (c(15,40))+ theme(legend.position = "none")

CHG_Rep <- data_TE %>%
  filter(Context == "CHG") %>%
  filter(Family == "Rep") %>%
  ggplot (aes(x = Position, y = Percentage, color = Sample))+
  geom_line(linewidth = 1)+
  geom_vline(xintercept = c(20,40),linetype="dashed")+
  theme_classic()+
  theme(plot.subtitle = element_text(color = "#999999"))+
  scale_color_manual(values = vivid_palette)+
  scale_x_continuous(breaks = c(0,20,40,60),
                     labels = c("-2kb", "TSS", "TES", "+2kb"))+
  labs(title = title_Rep, subtitle = n_Rep, x = "", y = "")+
  ylim (c(15,40))+ theme(legend.position = "none")

CHG_Sat <- data_TE %>%
  filter(Context == "CHG") %>%
  filter(Family == "Sat") %>%
  ggplot (aes(x = Position, y = Percentage, color = Sample))+
  geom_line(linewidth = 1)+
  geom_vline(xintercept = c(20,40),linetype="dashed")+
  theme_classic()+
  theme(plot.subtitle = element_text(color = "#999999"))+
  scale_color_manual(values = vivid_palette)+
  scale_x_continuous(breaks = c(0,20,40,60),
                     labels = c("-2kb", "TSS", "TES", "+2kb"))+
  labs(title = title_Sat, subtitle = n_Sat, x = "", y = "")+
  ylim (c(40,90))+ theme(legend.position = "none")

ggarrange(CHG_Copia , CHG_Gypsy , CHG_TRIM  , CHG_LTRun , CHG_CACTA, CHG_Mut,
          CHG_PIF   , CHG_Tc1   , CHG_hAT   , CHG_TIRun , CHG_LINE, CHG_SINE,
          CHG_Nonun , CHG_HEL   , CHG_Cen   , CHG_Rep   , CHG_Sat , TE_legend,
          ncol = 6, nrow = 3,
          widths  = c(2, 2, 2, 2, 2, 2)) #8.0x15.0

#CHH methylation
CHH_Copia <- data_TE %>%
  filter(Context == "CHH") %>%
  filter(Family == "Copia") %>%
  ggplot (aes(x = Position, y = Percentage, color = Sample))+
  geom_line(linewidth = 1)+
  geom_vline(xintercept = c(20,40),linetype="dashed")+
  theme_classic()+
  theme(plot.subtitle = element_text(color = "#999999"))+
  scale_color_manual(values = vivid_palette)+
  scale_x_continuous(breaks = c(0,20,40,60),
                     labels = c("-2kb", "TSS", "TES", "+2kb"))+
  labs(title = title_Copia, subtitle = n_Copia, x = "", y = "DNA methylation (%)")+
  scale_y_continuous(breaks = c(6,8,10,12), limits = c(5,12.5))+ 
  theme(legend.position = "none")

CHH_Gypsy <- data_TE %>%
  filter(Context == "CHH") %>%
  filter(Family == "Gypsy") %>%
  ggplot (aes(x = Position, y = Percentage, color = Sample))+
  geom_line(linewidth = 1)+
  geom_vline(xintercept = c(20,40),linetype="dashed")+
  theme_classic()+
  theme(plot.subtitle = element_text(color = "#999999"))+
  scale_color_manual(values = vivid_palette)+
  scale_x_continuous(breaks = c(0,20,40,60),
                     labels = c("-2kb", "TSS", "TES", "+2kb"))+
  labs(title = title_Gypsy, subtitle = n_Gypsy, x = "", y = "")+
  scale_y_continuous(breaks = c(6,8,10,12), limits = c(5,12.5))+ 
  theme(legend.position = "none")

CHH_TRIM <- data_TE %>%
  filter(Context == "CHH") %>%
  filter(Family == "TRIM") %>%
  ggplot (aes(x = Position, y = Percentage, color = Sample))+
  geom_line(linewidth = 1)+
  geom_vline(xintercept = c(20,40),linetype="dashed")+
  theme_classic()+
  theme(plot.subtitle = element_text(color = "#999999"))+
  scale_color_manual(values = vivid_palette)+
  scale_x_continuous(breaks = c(0,20,40,60),
                     labels = c("-2kb", "TSS", "TES", "+2kb"))+
  labs(title = title_TRIM, subtitle = n_TRIM, x = "", y = "")+
  ylim (c(5,25))+ theme(legend.position = "none")

CHH_LTRun <- data_TE %>%
  filter(Context == "CHH") %>%
  filter(Family == "LTRun") %>%
  ggplot (aes(x = Position, y = Percentage, color = Sample))+
  geom_line(linewidth = 1)+
  geom_vline(xintercept = c(20,40),linetype="dashed")+
  theme_classic()+
  theme(plot.subtitle = element_text(color = "#999999"))+
  scale_color_manual(values = vivid_palette)+
  scale_x_continuous(breaks = c(0,20,40,60),
                     labels = c("-2kb", "TSS", "TES", "+2kb"))+
  labs(title = title_LTRun, subtitle = n_LTRun, x = "", y = "")+
  scale_y_continuous(breaks = c(4,8,12), limits = c(2.5,12.5))+ 
  theme(legend.position = "none")

CHH_CACTA <- data_TE %>%
  filter(Context == "CHH") %>%
  filter(Family == "CACTA") %>%
  ggplot (aes(x = Position, y = Percentage, color = Sample))+
  geom_line(linewidth = 1)+
  geom_vline(xintercept = c(20,40),linetype="dashed")+
  theme_classic()+
  theme(plot.subtitle = element_text(color = "#999999"))+
  scale_color_manual(values = vivid_palette)+
  scale_x_continuous(breaks = c(0,20,40,60),
                     labels = c("-2kb", "TSS", "TES", "+2kb"))+
  labs(title = title_CACTA, subtitle = n_CACTA, x = "", y = "")+
  ylim (c(5,20))+ theme(legend.position = "none")

CHH_Mut <- data_TE %>%
  filter(Context == "CHH") %>%
  filter(Family == "Mut") %>%
  ggplot (aes(x = Position, y = Percentage, color = Sample))+
  geom_line(linewidth = 1)+
  geom_vline(xintercept = c(20,40),linetype="dashed")+
  theme_classic()+
  theme(plot.subtitle = element_text(color = "#999999"))+
  scale_color_manual(values = vivid_palette)+
  scale_x_continuous(breaks = c(0,20,40,60),
                     labels = c("-2kb", "TSS", "TES", "+2kb"))+
  labs(title = title_Mut, subtitle = n_Mut, x = "", y = "")+
  ylim (c(5,20))+ theme(legend.position = "none")

CHH_PIF <- data_TE %>%
  filter(Context == "CHH") %>%
  filter(Family == "PIF") %>%
  ggplot (aes(x = Position, y = Percentage, color = Sample))+
  geom_line(linewidth = 1)+
  geom_vline(xintercept = c(20,40),linetype="dashed")+
  theme_classic()+
  theme(plot.subtitle = element_text(color = "#999999"))+
  scale_color_manual(values = vivid_palette)+
  scale_x_continuous(breaks = c(0,20,40,60),
                     labels = c("-2kb", "TSS", "TES", "+2kb"))+
  labs(title = title_PIF, subtitle = n_PIF, x = "", y = "DNA methylation (%)")+
  ylim (c(0,50))+ theme(legend.position = "none")

CHH_Tc1 <- data_TE %>%
  filter(Context == "CHH") %>%
  filter(Family == "Tc1") %>%
  ggplot (aes(x = Position, y = Percentage, color = Sample))+
  geom_line(linewidth = 1)+
  geom_vline(xintercept = c(20,40),linetype="dashed")+
  theme_classic()+
  theme(plot.subtitle = element_text(color = "#999999"))+
  scale_color_manual(values = vivid_palette)+
  scale_x_continuous(breaks = c(0,20,40,60),
                     labels = c("-2kb", "TSS", "TES", "+2kb"))+
  labs(title = title_Tc1, subtitle = n_Tc1, x = "", y = "")+
  ylim (c(0,35))+ theme(legend.position = "none")

CHH_hAT <- data_TE %>%
  filter(Context == "CHH") %>%
  filter(Family == "hAT") %>%
  ggplot (aes(x = Position, y = Percentage, color = Sample))+
  geom_line(linewidth = 1)+
  geom_vline(xintercept = c(20,40),linetype="dashed")+
  theme_classic()+
  theme(plot.subtitle = element_text(color = "#999999"))+
  scale_color_manual(values = vivid_palette)+
  scale_x_continuous(breaks = c(0,20,40,60),
                     labels = c("-2kb", "TSS", "TES", "+2kb"))+
  labs(title = title_hAT, subtitle = n_hAT, x = "", y = "")+
  scale_y_continuous(breaks = c(5,10,15), limits = c(5,17.5))+ 
  theme(legend.position = "none")

CHH_TIRun <- data_TE %>%
  filter(Context == "CHH") %>%
  filter(Family == "TIRun") %>%
  ggplot (aes(x = Position, y = Percentage, color = Sample))+
  geom_line(linewidth = 1)+
  geom_vline(xintercept = c(20,40),linetype="dashed")+
  theme_classic()+
  theme(plot.subtitle = element_text(color = "#999999"))+
  scale_color_manual(values = vivid_palette)+
  scale_x_continuous(breaks = c(0,20,40,60),
                     labels = c("-2kb", "TSS", "TES", "+2kb"))+
  labs(title = title_TIRun, subtitle = n_TIRun, x = "", y = "")+
  scale_y_continuous(breaks = c(5,10,15), limits = c(5,17.5))+ 
  theme(legend.position = "none")

CHH_Cen <- data_TE %>%
  filter(Context == "CHH") %>%
  filter(Family == "Cen") %>%
  ggplot (aes(x = Position, y = Percentage, color = Sample))+
  geom_line(linewidth = 1)+
  geom_vline(xintercept = c(20,40),linetype="dashed")+
  theme_classic()+
  theme(plot.subtitle = element_text(color = "#999999"))+
  scale_color_manual(values = vivid_palette)+
  scale_x_continuous(breaks = c(0,20,40,60),
                     labels = c("-2kb", "TSS", "TES", "+2kb"))+
  labs(title = title_Cen, subtitle = n_Cen, x = "", y = "")+
  ylim (c(0,15))+ theme(legend.position = "none")

CHH_LINE <- data_TE %>%
  filter(Context == "CHH") %>%
  filter(Family == "LINE") %>%
  ggplot (aes(x = Position, y = Percentage, color = Sample))+
  geom_line(linewidth = 1)+
  geom_vline(xintercept = c(20,40),linetype="dashed")+
  theme_classic()+
  theme(plot.subtitle = element_text(color = "#999999"))+
  scale_color_manual(values = vivid_palette)+
  scale_x_continuous(breaks = c(0,20,40,60),
                     labels = c("-2kb", "TSS", "TES", "+2kb"))+
  labs(title = title_LINE, subtitle = n_LINE, x = "", y = "")+
  scale_y_continuous(breaks = c(5,10,15), limits = c(4,15))+ 
  theme(legend.position = "none")

CHH_SINE <- data_TE %>%
  filter(Context == "CHH") %>%
  filter(Family == "SINE") %>%
  ggplot (aes(x = Position, y = Percentage, color = Sample))+
  geom_line(linewidth = 1)+
  geom_vline(xintercept = c(20,40),linetype="dashed")+
  theme_classic()+
  theme(plot.subtitle = element_text(color = "#999999"))+
  scale_color_manual(values = vivid_palette)+
  scale_x_continuous(breaks = c(0,20,40,60),
                     labels = c("-2kb", "TSS", "TES", "+2kb"))+
  labs(title = title_SINE, subtitle = n_SINE, x = "", y = "")+
  ylim (c(0,35))+ theme(legend.position = "none")

CHH_Nonun <- data_TE %>%
  filter(Context == "CHH") %>%
  filter(Family == "Nonun") %>%
  ggplot (aes(x = Position, y = Percentage, color = Sample))+
  geom_line(linewidth = 1)+
  geom_vline(xintercept = c(20,40),linetype="dashed")+
  theme_classic()+
  theme(plot.subtitle = element_text(color = "#999999"))+
  scale_color_manual(values = vivid_palette)+
  scale_x_continuous(breaks = c(0,20,40,60),
                     labels = c("-2kb", "TSS", "TES", "+2kb"))+
  labs(title = title_Nonun, subtitle = n_Nonun, x = "", y = "DNA methylation (%)")+
  ylim (c(0,15))+ theme(legend.position = "none")

CHH_HEL <- data_TE %>%
  filter(Context == "CHH") %>%
  filter(Family == "HEL") %>%
  ggplot (aes(x = Position, y = Percentage, color = Sample))+
  geom_line(linewidth = 1)+
  geom_vline(xintercept = c(20,40),linetype="dashed")+
  theme_classic()+
  theme(plot.subtitle = element_text(color = "#999999"))+
  scale_color_manual(values = vivid_palette)+
  scale_x_continuous(breaks = c(0,20,40,60),
                     labels = c("-2kb", "TSS", "TES", "+2kb"))+
  labs(title = title_HEL, subtitle = n_HEL, x = "", y = "")+
  ylim (c(0,15))+ theme(legend.position = "none")

CHH_Rep <- data_TE %>%
  filter(Context == "CHH") %>%
  filter(Family == "Rep") %>%
  ggplot (aes(x = Position, y = Percentage, color = Sample))+
  geom_line(linewidth = 1)+
  geom_vline(xintercept = c(20,40),linetype="dashed")+
  theme_classic()+
  theme(plot.subtitle = element_text(color = "#999999"))+
  scale_color_manual(values = vivid_palette)+
  scale_x_continuous(breaks = c(0,20,40,60),
                     labels = c("-2kb", "TSS", "TES", "+2kb"))+
  labs(title = title_Rep, subtitle = n_Rep, x = "", y = "")+
  scale_y_continuous(breaks = c(0,5,10,15), limits = c(0,17.5))+ 
  theme(legend.position = "none")

CHH_Sat <- data_TE %>%
  filter(Context == "CHH") %>%
  filter(Family == "Sat") %>%
  ggplot (aes(x = Position, y = Percentage, color = Sample))+
  geom_line(linewidth = 1)+
  geom_vline(xintercept = c(20,40),linetype="dashed")+
  theme_classic()+
  theme(plot.subtitle = element_text(color = "#999999"))+
  scale_color_manual(values = vivid_palette)+
  scale_x_continuous(breaks = c(0,20,40,60),
                     labels = c("-2kb", "TSS", "TES", "+2kb"))+
  labs(title = title_Sat, subtitle = n_Sat, x = "", y = "")+
  ylim (c(0,15))+ theme(legend.position = "none")

ggarrange(CHH_Copia , CHH_Gypsy , CHH_TRIM  , CHH_LTRun , CHH_CACTA, CHH_Mut,
          CHH_PIF   , CHH_Tc1   , CHH_hAT   , CHH_TIRun , CHH_LINE, CHH_SINE,
          CHH_Nonun , CHH_HEL   , CHH_Cen   , CHH_Rep   , CHH_Sat , TE_legend,
          ncol = 6, nrow = 3,
          widths  = c(2, 2, 2, 2, 2, 2)) #8.0x15.0
