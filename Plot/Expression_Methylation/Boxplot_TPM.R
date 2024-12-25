############################################################################
### Import libraries
library (tidyverse)
library (ggplot2) 
library (reshape2)
library (scales)
library (ggpubr)
library (cowplot)
############################################################################
### Import files
file <- "FG.TPM.txt"
data <- read.csv(file, header = T, sep = "\t")

############################################################################
### Pre-process files
exp_decile <- function(x){
  cat = ifelse(x<0.5, 0, cut(x, breaks = quantile(x[x>=0.5],
                                                  prob = c(0,0.2,0.4,0.6,0.8,1)),
                             include.lowest = TRUE))
  cat = case_when(
    cat == 0 ~ "NoExp",
    cat == 1 ~ "2nd",
    cat == 2 ~ "4th",
    cat == 3 ~ "6th",
    cat == 4 ~ "8th",
    cat == 5 ~ "10th"
  )
  return (cat)
}

data <- data %>%
  mutate (GENE = gsub("\\.MSUv7\\.0", "", GENE)) %>%
  mutate(FG1 = rowMeans(select(data, starts_with("FG1")),na.rm = TRUE),
         FG2 = rowMeans(select(data, starts_with("FG2")),na.rm = TRUE),
         FG3 = rowMeans(select(data, starts_with("FG3")),na.rm = TRUE),
         FG4 = rowMeans(select(data, starts_with("FG4")),na.rm = TRUE),
         FG5 = rowMeans(select(data, starts_with("FG5")),na.rm = TRUE),
         FG6 = rowMeans(select(data, starts_with("FG6")),na.rm = TRUE),
         FL  = rowMeans(select(data, starts_with("FL")), na.rm = TRUE)) %>%
  mutate(FG1_cat = exp_decile(FG1),
         FG2_cat = exp_decile(FG2),
         FG3_cat = exp_decile(FG3),
         FG4_cat = exp_decile(FG4),
         FG5_cat = exp_decile(FG5),
         FG6_cat = exp_decile(FG6),
         FL_cat  = exp_decile(FL)) %>%
  select (GENE,FG1,FG2,FG3,FG4,FG5,FG6,FL,
          FG1_cat,FG2_cat,FG3_cat,FG4_cat,FG5_cat,FG6_cat,FL_cat)

############################################################################

## Export table
data_FG1 <- data %>% select (FG1_cat, GENE) %>% 
  write.table(file = "FG1_exp.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
data_FG2 <- data %>% select (FG2_cat, GENE) %>% 
  write.table(file = "FG2_exp.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
data_FG3 <- data %>% select (FG3_cat, GENE) %>% 
  write.table(file = "FG3_exp.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
data_FG4 <- data %>% select (FG4_cat, GENE) %>% 
  write.table(file = "FG4_exp.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
data_FG5 <- data %>% select (FG5_cat, GENE) %>% 
  write.table(file = "FG5_exp.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
data_FG6 <- data %>% select (FG6_cat, GENE) %>% 
  write.table(file = "FG6_exp.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
data_FL <- data %>% select (FL_cat, GENE) %>% 
  write.table(file = "FL_exp.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

data_export <- data %>% select (GENE, FG1_cat, FG2_cat, FG3_cat, FG4_cat, FG5_cat, FG6_cat, FL_cat)
colnames(data_export) <- c("gene","FG1","FG2","FG3","FG4","FG5","FG6","FL")
write.csv(data_export, file = "gene_group_exp.csv", row.names = FALSE, quote = FALSE)

############################################################################

###Boxplot - TPM expression
exp_palette = c("#666666","#65a14f","#9dc489","#e4d189","#f4a979","#e86473")
cat = c("NoExp","2nd","4th","6th","8th","10th")
cat_label = c("No expression","Very low", "Low", "Moderate", "High", "Very high")
cat_mapping = setNames(cat_label, cat)

data_box <- melt (data, id.vars = c("GENE","FG1_cat","FG2_cat","FG3_cat",
                                           "FG4_cat","FG5_cat","FG6_cat","FL_cat"),
                  variable.name = "sample", value.name = "exp")
data_box <- data_box %>%
  mutate (cat = case_when(
    sample == "FG1" ~ FG1_cat,
    sample == "FG2" ~ FG2_cat,
    sample == "FG3" ~ FG3_cat,
    sample == "FG4" ~ FG4_cat,
    sample == "FG5" ~ FG5_cat,
    sample == "FG6" ~ FG6_cat,
    sample == "FL"  ~ FL_cat)) %>%
  mutate (cat = cat_mapping[cat]) %>%
  mutate (exp = log10(exp+1)) %>%
  select (GENE, sample, exp, cat)
data_box$cat <- factor (data_box$cat, levels = cat_label)

p <- data_box %>%
  ggplot(aes(y = exp, fill = cat))+
  stat_boxplot(geom="errorbar", linewidth = 0.3,
               position = position_dodge(width = 1))+
  geom_boxplot(color = "black", linewidth = 0.3, outliers = FALSE,
               position = position_dodge(width = 1))+
  theme_classic()+
  scale_fill_manual(values = exp_palette)+
  facet_grid(.~sample, scales = "fixed") + 
  coord_cartesian(clip = 'off') + 
  labs(y = "log10(1+TPM)") +
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
  scale_y_continuous(breaks = c(0.0, 1.0, 2.0), expand = c(0, 0.1),
                     labels = function(x) sprintf("%.1f",x))

pdf(file = "boxplot_TPM.pdf", width = 7, height = 1)
print(p)
dev.off()
