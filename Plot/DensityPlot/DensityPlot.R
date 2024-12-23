############################################################################
### Import libraries
library(ggplot2)
library(tidyverse)
library(ggthemes)
library(ggpubr)
library(cowplot)
library(rstatix)
library(forcats)

############################################################################
##global parameters
min_C <- 20 #Minimal number of Cs (both samples)
min_ratio_CG <- 0.5 #Minimal ratio of CG methylation (either sample)
min_ratio_CHG <- 0.2 #Minimal ratio of CHG methylation (either sample)
min_ratio_CHH <- 0.1 #Minimal ratio of CHH methylation (either sample)
samples <- c("FG1", "FG2", "FG3", "FG4", "FG5", "FG6")

### Import files
## Set directory
dir_input  <- "./50-bp_bed"
dir_diff <- "./diff_csv_FG"

## Local parameters
context <- "CHH" #One of the followings: "CG", "CHG", "CHH"

if (context == "CG") { min_ratio <- min_ratio_CG
}else if (context == "CHG") { min_ratio <- min_ratio_CHG
}else if (context == "CHH") { min_ratio <- min_ratio_CHH}

### Process data
## Calculate fractional differences

for (i in 1:length(samples)) {
  # Sample assignment
  sample_target  <- samples[i]  # one FG stage
  sample_control <- samples[-i] # all other FG stages

  # File name assignment
  file_target <- paste0(sample_target, "_", context, "_w50.bed")
  files_control <- c()
  for (sample in sample_control){
    files_control <- c(files_control, paste0(sample, "_", context, "_w50.bed"))
  }
  file_diff <- paste0("diff_", sample_target, "_others_", context, ".csv")

  # Calculate fractional difference
  data_target  <- read.table(paste0(dir_input,"/",file_target))
  data_control <- data.frame (matrix(ncol = 6, nrow = 0))
  for (file_control in files_control){
    temp_data <- read.table(paste0(dir_input,"/",file_control))
    data_control <- rbind(data_control, temp_data)
    rm(temp_data)
  }
  colnames(data_target) <- c("chr","start","end","mC","C","ratio")
  colnames(data_control) <- c("chr","start","end","mC","C","ratio")
  
  data_target <- data_target %>%
    select (-c(end,ratio)) %>%
    filter (mC+C >= min_C) %>%
    mutate(ratio = mC/(mC+C)) %>%
    select(-c(mC,C))
  colnames(data_target) <- c("chr", "position", "ratio_target")
  
  data_control <- data_control %>%
    select(-c(end,ratio)) %>%
    group_by(chr, start) %>%
    summarise(mC = sum(mC), C = sum(C)) %>%
    filter (mC+C >= min_C) %>%
    mutate (ratio = mC/(mC+C)) %>%
    select(-c(mC,C))
  colnames(data_control) <- c("chr", "position", "ratio_control")
  
  data_diff <- inner_join(data_target, data_control, by = c("chr","position")) %>%
    filter ((ratio_target >= min_ratio ) | (ratio_control >= min_ratio)) %>%
    mutate (ratio_diff = ratio_target - ratio_control) %>%
    select (-c(ratio_target, ratio_control))
  
  #Export data
  write.csv(data_diff, file = paste0(dir_diff,"/",file_diff), row.names = FALSE)
}

################################################################################
###Plotting
## Plotting FG1-others (CHH)
ratio_filter <- -0.1
col_filter <- "#E41817"
file_filter <- "diff_FG1_others_CHH.csv"
data_filter <- read.csv(paste0(dir_diff,"/",file_filter))
density_filter <- density (data_filter$ratio_diff)
data_shade <- tibble (x = density_filter$x, y = density_filter$y) %>%
  mutate (shade = case_when(
    x < ratio_filter ~ "Y",
    .default = NA
  ))
p <- data_shade %>%
  ggplot(aes(x = x, y = y)) + 
  geom_vline(aes(xintercept = 0),linetype="dashed",color = "#666666", linewidth = .3) +
  geom_line(color = "black", linewidth = .4) + 
  geom_area(data = filter (data_shade, shade == "Y"), 
            color = col_filter, fill = col_filter, alpha = 0.5, linewidth = .5)+
  theme_classic()+
  annotate("segment", x=-0.1, xend=-0.1, y=0, yend = 1.7,
           color = col_filter, linewidth = .5)+
  annotate("text", x = -0.2, y = 2.5, label = "-0.1", color = col_filter, size = 2.5)+
  scale_x_continuous(limits = c(-0.4,0.4), breaks = c(-0.4,0.0,0.4), expand = c(0, 0))+
  scale_y_continuous(limits = function(limits) c(0, ceiling(limits[2])),
                     breaks = function(limits) c(0, ceiling(limits[2])),
                     expand = c(0, 0))+
  labs(title = "FG1")+
  theme(axis.line = element_line(color="black", linewidth = 0.4),
        axis.text = element_text(color="black", size=5),
        axis.ticks = element_line(color="black", linewidth = 0.4),
        axis.title = element_blank(),
        panel.border = element_blank(),
        strip.background = element_blank(),
        panel.background = element_blank(),
        strip.text = element_text(color="black"),
        plot.title = element_text(color="black",size=7)
  )
p_FG1_others_CHH <- p

## Plotting others
data_filter <- data_filter %>% 
  filter (ratio_diff < ratio_filter)

write.csv(data_filter, file = "./highlight.csv")

#Plotting & saving one-by-one
file_plot <- "diff_FG1_others_CG.csv" #Change this 
title_plot <- "FG6" #Change this 
data_plot <- read.csv(paste0(dir_diff,"/",file_plot,sep =""), header = TRUE)

data_plot_filter <- inner_join(data_plot, data_filter, by = c("chr","position")) %>%
  select(chr, position, ratio_diff.x)
colnames(data_plot_filter) <- c("chr", "position", "ratio_diff")

p <- data_plot %>%
  ggplot(aes(x=ratio_diff)) + 
  geom_vline(aes(xintercept = 0),linetype="dashed",color = "#666666", linewidth = .3) +
  geom_density(color = "black", linewidth = .4) + 
  geom_density(data = data_plot_filter, aes(x=ratio_diff), 
               fill = col_filter, alpha = 0.5, linewidth = 0) +
  theme_classic()+
  scale_x_continuous(limits = c(-0.4,0.4), breaks = c(-0.4,0,0.4), expand = c(0, 0))+
  # CG - 0.2
  # CHG, CHH - 0.4
  scale_y_continuous(limits = function(limits) c(0, ceiling(limits[2])),
                     breaks = function(limits) c(0, ceiling(limits[2])),
                     expand = c(0, 0))+
  labs(title = title_plot)+
  theme(axis.line = element_line(color="black", linewidth = 0.4),
        axis.text = element_text(color="black", size=5),
        axis.title = element_blank(),
        axis.ticks = element_line(color="black", linewidth = 0.4),
        panel.border = element_blank(),
        strip.background = element_blank(),
        panel.background = element_blank(),
        strip.text = element_text(color="black"),
        plot.title = element_text(color="black", size=7)
  )
p_FG6_others_CHH <- p #Change this

# Plot everything
p_combined <- ggarrange(p_FG1_others_CHH,p_FG2_others_CHH,p_FG3_others_CHH,
                        p_FG4_others_CHH,p_FG5_others_CHH,p_FG6_others_CHH,
                        nrow = 1, widths = c(2,2,2,2,2,2))

pdf(file = paste("densityplot_FG_CHH.pdf",sep=""), width = 6, height = 1.2)
print(p_combined)
dev.off()
