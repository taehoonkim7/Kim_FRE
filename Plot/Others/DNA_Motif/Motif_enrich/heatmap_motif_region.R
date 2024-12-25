############################################################################
### Import libraries
library (ggplot2)
library (reshape2)
library (tidyverse)
library (raster)
library (scales)
library (factoextra)
library (vegan)
library (stats)

############################################################################
### Import files
# Import data
file = "motif_enrich_region.csv"
data = read.csv(file, header = T)

regions = c("FG-dDMR_FG-siren","FG-dDMR_only","FG-siren_only","Shuffle")
motifs = unique (data$Motif)
offsets = unique (data$Offset)

data$Region = factor(data$Region, levels = rev(regions))
data$Motif = factor(data$Motif, levels = motifs)
data$Offset = factor(data$Offset, levels = offsets)

############################################################################
## Data preparation
# Fisher's exact test
df_pair <- as.data.frame(t(combn(regions, 2)), stringsAsFactors = FALSE)
colnames(df_pair) <- c("Region_A", "Region_B")
df_comb <- expand.grid(Motif = motifs,Offset = offsets,stringsAsFactors = FALSE)
data_fisher <- merge (df_pair, df_comb, by = NULL)
data_fisher$pval = NA

for (i in 1:nrow(data_fisher)){
  A <- data %>% 
    filter (Region == data_fisher[i,]$Region_A,
            Motif  == data_fisher[i,]$Motif,
            Offset == data_fisher[i,]$Offset)
  B <- data %>% 
    filter (Region == data_fisher[i,]$Region_B,
            Motif  == data_fisher[i,]$Motif,
            Offset == data_fisher[i,]$Offset)
  df_test <- data.frame(
    "A" = c (A$Count, A$Total-A$Count),
    "B" = c (B$Count, B$Total-B$Count),
    row.names = c("True", "False")
  )
  stat <- fisher.test(df_test)
  data_fisher[i,]$pval = stat$p.value
}

rm (A, B, df_comb, df_pair, df_test, stat, file, i)

# further processing
data <- data %>%
  mutate (Perc = Count/Total*100)

data_fisher <- data_fisher %>% 
  mutate (Comparison = paste0(Region_A,"_vs_",Region_B)) %>%
  mutate (pval = ifelse(pval > 1, 1, pval)) %>%
  mutate (neg_log10_pval = -log10(pval))
data_fisher$Comparison = factor (data_fisher$Comparison, 
                                 levels = unique(rev(data_fisher$Comparison)))

############################################################################
###Plotting
## Enrichment
high_col <- "#A42B2A"
mid_col <- "#f3f4cc"  
low_col <- "#f1f1f1"

p_enrich <- ggplot (data, aes(x = Motif, y = Region, fill = Perc)) + 
  geom_tile() + 
  coord_fixed() +
  facet_grid(Offset~.) + 
  scale_fill_gradientn(colors = c(low_col, mid_col, high_col),
                       limits = c(0, 100),
                       breaks = c(0,25,50,75,100),
                       values = rescale(c(0,40,100))) + 
  theme_classic() + 
  theme(text = element_text(color="black", size = 5),
        axis.title.x = element_blank(),
        axis.title.y = element_text(color="black", size = 6),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_text(color="black", size = 5),
        axis.text.x = element_text(color="black", size = 5, 
                                   angle = 90, hjust = 0, vjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill = NA, color = "black",linewidth = 0.3),
        panel.spacing = unit (0.2,"mm"),
        strip.background = element_blank(),
        strip.text = element_text(color="black", size = 5),
        legend.background = element_blank(),
        legend.title = element_text(color="black", size = 6, 
                                    angle = -90, hjust = 0.5),
        legend.title.position = "right",
        legend.key.size = unit (2,"mm"),
        legend.text = element_text(color="black", size = 5, margin = margin(l = 2))) + 
  labs(fill = "%Enrichment", y = "Genomic region category")

pdf(file = "heatmap_motif_enrich.pdf", width = 5, height = 2)
print(p_enrich)
dev.off()

## Fisher
p_fisher <- ggplot (data_fisher, aes(x = Motif, y = Comparison, fill = neg_log10_pval)) + 
  geom_tile() + 
  coord_fixed() +
  facet_grid(Offset~.) + 
  scale_fill_gradientn(colors = c(low_col, mid_col, high_col),
                       limits = c(0, 80),breaks = c(0,20,40,60,80),
                       values = rescale(c(0,1,80))) + 
  theme_classic() + 
  theme(text = element_text(color="black", size = 5),
        axis.title.x = element_blank(),
        axis.title.y = element_text(color="black", size = 6),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_text(color="black", size = 5),
        axis.text.x = element_text(color="black", size = 5, 
                                   angle = 90, hjust = 0, vjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill = NA, color = "black",linewidth = 0.3),
        panel.spacing = unit (0.2,"mm"),
        strip.background = element_blank(),
        strip.text = element_text(color="black", size = 5),
        legend.background = element_blank(),
        legend.title = element_text(color="black", size = 6, 
                                    angle = -90, hjust = 0.5),
        legend.title.position = "right",
        legend.key.size = unit (2,"mm"),
        legend.text = element_text(color="black", size = 5, margin = margin(l = 2))) + 
  labs(fill = "-log10(Pvalue)", y = "Comparisons")

pdf(file = "heatmap_motif_fisher.pdf", width = 5, height = 2.5)
print(p_fisher)
dev.off()
