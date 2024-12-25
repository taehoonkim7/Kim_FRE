############################################################################
### Import libraries
library(stats)
library(ggplot2)
library(dplyr)

############################################################################
### Import files
data <- read.csv("TE_FG-dDMR.csv", header = TRUE)
TE_total = sum(data[,3])
DMR_total = sum(data[,2])
NonDMR_total = TE_total - DMR_total

############################################################################
### Process data & Fisher's exact test
TEs = c(data$TE_name)
res = data.frame(row.names = TEs)

for (i in TEs){
  total_class = data[data$TE_name==i,3]
  DMR_class = data[data$TE_name==i,2]
  NonDMR_class = total_class - DMR_class
  
  df <- data.frame(
    "DMR" = c(DMR_class, DMR_total-DMR_class),
    "NonDMR" = c(NonDMR_class, NonDMR_total-NonDMR_class),
    row.names = c("Class", "NonClass")
  )
  stat <- fisher.test(df)
  data[data$TE_name==i,"O_E"] = stat$estimate
  data[data$TE_name==i,"O_E_log2"] = log2(data[data$TE_name==i,"O_E"])
  data[data$TE_name==i,"Fisher.p"] = stat$p.value
}

write.csv(data, file = "TE_contingency_fisher.csv")

############################################################################
### Plotting
data_p <- data %>%
  filter(observed != 0) %>%
  filter(TE_name != "Repeat") %>%
  mutate(Highlight = case_when((Fisher.p < 0.05) & (O_E_log2 > 0) ~ "Y",
                               .default = "N"),
         Annotation = case_when((Fisher.p < 0.05) & (Fisher.p >= 0.01) ~ "*",
                                (Fisher.p < 0.01) & (Fisher.p >= 0.001) ~ "**",
                                Fisher.p < 0.001 ~ "***",
                                .default = ""),
         Annotation.Position = case_when (O_E_log2 > 0 ~ -0.4,
                                          O_E_log2 < 0 ~ 0.2)) %>%
  arrange(desc(O_E_log2))
data_p$TE_name = factor(data_p$TE_name, levels = data_p$TE_name)
  
p <- ggplot(data_p, aes(x = TE_name, y = O_E_log2, fill = Highlight)) +
  geom_bar(stat = "identity", position = "identity",
           color = "black", linewidth = 0.2, width = 0.75)+
  geom_hline(yintercept = 0, linewidth = 0.2)+
  geom_text(aes(y = Annotation.Position, label = Annotation),
            position = position_identity(), vjust = 0.5, size = 1.75)+
  theme_classic()+
  scale_fill_manual(values = c("Y" = "#7db6d9", "N" = "black"))+
  labs(x = element_blank(), y = "O/E (log2)")+
  theme(text = element_text(color = "black", size = 5),
        axis.line = element_line(color="black", linewidth = 0.2),
        axis.title.y = element_text(color = "black", size = 7),
        axis.text.y = element_text(color = "black", size = 5),
        axis.text.x = element_text(color = "black", size = 5,
                                   angle=45, vjust=1, hjust=1),
        axis.ticks = element_line(color="black", linewidth = 0.2),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "None")

pdf(file = "TE_FG-dDMR.pdf", width = 2, height = 1.4)
p
dev.off()

