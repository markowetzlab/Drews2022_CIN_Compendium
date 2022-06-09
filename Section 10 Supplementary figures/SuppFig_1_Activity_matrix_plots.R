## Make activity heatmaps (normalised)

rm(list=ls(all=TRUE))

library(data.table)
library(reshape2)
library(ggplot2)
library(ggthemes)
library(ggpubr)
library(lemon)
library(cowplot)
library(patchwork)
library(ggrastr)
library(this.path)

theme_set(theme_tufte(base_size = 6, base_family = "ArialMT"))
theme_update(text = element_text(size = 6),
             axis.text = element_text(size = 6),
             axis.line = element_line(size = 0.5), 
             axis.ticks = element_line(size = 0.5),
             axis.ticks.length = unit(.1, "cm"))

#PATHS
BASE=dirname(this.path())
INP=file.path(BASE, "input")
OUT=file.path(BASE, "output")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

## Plot 1: Normalised signatures
exp = readRDS(file.path(INP, "Signature_Compendium_v5_Cosine-0.74_Activities_THRESH95_NAMESAPRIL21.rds"))

hExp = hclust(dist(exp, method = "manhattan"), method = "ward.D2")
dtExp = data.table(melt(exp))
dtExp$Var1 = factor(dtExp$Var1, levels = rownames(exp)[ hExp$order ])
dtExp$Var2 = factor(as.character(dtExp$Var2), levels = paste0("CX", 17:1))

p1A = ggplot(dtExp, aes(y = Var2, x = Var1, fill = value)) + rasterise(geom_tile(), dpi = 600) + 
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) + 
  labs(y = "Signatures", x = "TCGA samples", fill = "Activity") +
  scale_fill_gradientn(colours = c("white", "#b3cde3", "#8c96c6", "#88419d"),
                       values = c(0.0e+00, 2.5e-05, 2.5e-02, 1),
                       breaks = c(0, 0.10, 0.5, 1),
                       labels = c(0, 0.10, 0.5, 1)) +
  coord_capped_cart(left = "both", bottom = "both")


# Panel A cancer types
meta = readRDS(file.path(INP, "Metadata_TCGA_ASCAT_penalty70.rds"))
cols = read.csv(file.path(INP, "TCGA_colour_scheme.txt"), header = FALSE, quote = "|", sep = "\t")
vCols = as.character(cols$V2)
names(vCols) = cols$V1

dtExp$Cancer = meta$cancer_type[ match(as.character(dtExp$Var1), meta$name) ]
dtExp$Cancer = factor(dtExp$Cancer, levels = sort(unique(meta$cancer_type), decreasing = TRUE))

p1ACancer = ggplot(dtExp[ dtExp$Var2 == "CX1", ], aes(y = 1, x = Var1, fill = Cancer)) + rasterise(geom_tile(), dpi = 600) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), 
        axis.text.x = element_blank(), axis.ticks.x =  element_blank()) + 
  labs(x = element_blank(), y = "Cancer") + scale_fill_manual(values = vCols) +
  guides(fill = guide_legend(nrow = 3, byrow = FALSE))

# Extract legend
p1ALegend = as_ggplot(get_legend(p1ACancer))
p1ACancer = p1ACancer + theme(legend.position = "none")




# Panel B: Sorted by cancer type
# Sort by cancer type and CX1 activity
dtCX1 = dtExp[ dtExp$Var2 == "CX1", ]
newOrder = as.character(dtCX1$Var1)[ order(c(dtCX1$Cancer, dtCX1$value), decreasing = TRUE) ]
dtExp$Var1 = factor(as.character(dtExp$Var1), levels = newOrder)

p1B = ggplot(dtExp, aes(y = Var2, x = Var1, fill = value)) + rasterise(geom_tile(), dpi = 600) + 
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) + 
  labs(y = "Signatures", x = "TCGA samples", fill = "Activity") +
  scale_fill_gradientn(colours = c("white", "#b3cde3", "#8c96c6", "#88419d"),
                       values = c(0.0e+00, 2.5e-05, 2.5e-02, 1),
                       breaks = c(0, 0.10, 0.5, 1),
                       labels = c(0, 0.10, 0.5, 1)) +
  coord_capped_cart(left = "both", bottom = "both")

p1BCancer = ggplot(dtExp[ dtExp$Var2 == "CX1", ], aes(y = 1, x = Var1, fill = Cancer)) + rasterise(geom_tile(), dpi = 600) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), 
        axis.text.x = element_blank(), axis.ticks.x =  element_blank(), legend.position = "none") + 
  labs(x = element_blank(), y = "Cancer") + scale_fill_manual(values = vCols) +
  guides(fill = guide_legend(nrow = 3, byrow = FALSE))


# On Windows this command saves a smaller svg than specified?! Manually editing the size in Affinity.
ggsave(file.path(OUT, "SuppFig_1_Activity_Panel_A.svg"), p1A, width = 160, height = 80, units = "mm")
# cairo_pdf(file.path(OUT, "SuppFig_Activity_Panel_A.pdf"), width = 160/25.4, height = 80/25.4); print(p1A); dev.off()
ggsave(file.path(OUT, "SuppFig_1_Activity_Panel_A_Cancer.svg"), p1ACancer, width = 160, height = 5, units = "mm")
ggsave(file.path(OUT, "SuppFig_1_Activity_Panel_A_Legend.svg"), p1ALegend, width = 160, height = 15, units = "mm")

ggsave(file.path(OUT, "SuppFig_1_Activity_Panel_B.svg"), p1B, width = 160, height = 80, units = "mm")
ggsave(file.path(OUT, "SuppFig_1_Activity_Panel_B_Cancer.svg"), p1BCancer, width = 160, height = 5, units = "mm")

