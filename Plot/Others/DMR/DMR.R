############################################################################
### Import libraries
library (ggplot2)
library (tidyverse)
library (ggthemes)
library (ggpubr)
library (cowplot)

############################################################################
### Import file
#Import file
data <- read.csv("DMR.csv")
colnames (data) <- c("Context","Sample","Hypo","Hyper")

# Process data
data <- data %>%
  mutate (Hypo = -Hypo) %>%
  pivot_longer(cols = c("Hypo","Hyper"),
               names_to = "Type", values_to = "Value") %>%
  mutate (Count = abs(Value),
          Offset = case_when((Context == "CG") & (Type == "Hypo") & (Count > 500) ~ Value+30,
                             (Context == "CG") & (Type == "Hypo") & (Count <= 500) ~ Value-30,
                             (Context == "CG") & (Type == "Hyper") ~ Value+30,
                             (Context == "CHG") & (Type == "Hypo") ~ Value-10,
                             (Context == "CHG") & (Type == "Hyper")~ Value+10,
                             (Context == "CHH") & (Type == "Hypo") & (Count > 3000) ~ Value+500,
                             (Context == "CHH") & (Type == "Hypo") & (Count <= 3000) ~ Value-500,
                             (Context == "CHH") & (Type == "Hyper") & (Count > 3000) ~ Value-500,
                             (Context == "CHH") & (Type == "Hyper") & (Count <= 3000) ~ Value+500))

data$Sample <- factor(data$Sample,
                      levels = c("FG1","FG2","FG3","FG4","FG5","FG6"))
data$Context <- factor(data$Context,
                       levels = c("CG","CHG","CHH"))

############################################################################
### Plotting
col_hypo  <- "#ffd966"
col_hyper <- "#a9d18e"

#CG
p_CG <- data %>%
  filter (Context == "CG") %>%
  ggplot (aes(x = Sample, y = Value, fill = Type)) + 
  geom_bar(stat = "identity", position = "identity",
           color = "black", linewidth = 0.2, width = .8) +
  geom_text(aes(label = Count, y = Offset), position = position_identity(), 
            vjust = 0.5, color = "black", size = 1.6) +
  geom_hline(aes(yintercept = 0), color = "black", linewidth = 0.2)+
  scale_fill_manual(values = c(col_hyper, col_hypo)) +
  theme_classic() + 
  theme(axis.line = element_blank(),
        axis.title = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(color="black", size = 6),
        axis.ticks = element_blank(),
        legend.title = element_blank(),
        legend.key.size = unit(3,"mm"),
        legend.text = element_text(color = "black", size = 5),
        legend.position = c(.8,.3),
        legend.background = element_blank()
  )+
  scale_y_continuous(expand = c(0, 30))

pdf(file = paste("DMR_CG.pdf",sep=""), width = 1.4, height = 1.2)
print(p_CG)
dev.off()

#CHG
p_CHG <- data %>%
  filter (Context == "CHG") %>%
  ggplot (aes(x = Sample, y = Value, fill = Type)) + 
  geom_bar(stat = "identity", position = "identity",
           color = "black", linewidth = 0.2, width = .8) +
  geom_text(aes(label = Count, y = Offset), position = position_identity(), 
            vjust = 0.5, color = "black", size = 1.6) +
  geom_hline(aes(yintercept = 0), color = "black", linewidth = 0.2)+
  scale_fill_manual(values = c(col_hyper, col_hypo)) +
  theme_classic() + 
  theme(axis.line = element_blank(),
        axis.title = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(color="black", size = 6),
        axis.ticks = element_blank(),
        legend.title = element_blank(),
        legend.key.size = unit(3,"mm"),
        legend.text = element_text(color = "black", size = 5),
        legend.position = c(.375,.85),
        legend.background = element_blank()
  )+
  guides(fill = guide_legend(nrow=1,byrow=TRUE))+
  scale_y_continuous(expand = c(0, 10))

pdf(file = paste("DMR_CHG.pdf",sep=""), width = 1.4, height = 1.2)
print(p_CHG)
dev.off()

#CHH
p_CHH <- data %>%
  filter (Context == "CHH") %>%
  ggplot (aes(x = Sample, y = Value, fill = Type)) + 
  geom_bar(stat = "identity", position = "identity",
           color = "black", linewidth = 0.2, width = .8) +
  geom_text(aes(label = Count, y = Offset), position = position_identity(), 
            vjust = 0.5, color = "black", size = 1.6) +
  geom_hline(aes(yintercept = 0), color = "black", linewidth = 0.2)+
  scale_fill_manual(values = c(col_hyper, col_hypo)) +
  theme_classic() + 
  theme(axis.line = element_blank(),
        axis.title = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(color="black", size = 6),
        axis.ticks = element_blank(),
        legend.title = element_blank(),
        legend.key.size = unit(3,"mm"),
        legend.text = element_text(color = "black", size = 5),
        legend.position = c(.25,0.8),
        legend.background = element_blank()
  )+
  scale_y_continuous(expand = c(0, 100))

pdf(file = paste("DMR_CHH.pdf",sep=""), width = 1.4, height = 1.2)
print(p_CHH)
dev.off()



