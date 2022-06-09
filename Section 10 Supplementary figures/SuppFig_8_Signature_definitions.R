## Create summary of biological covariates

rm(list=ls(all=TRUE))

library(this.path)
library(data.table)
library(reshape2)
library(ggplot2)
library(ggthemes)
library(lemon)
library(RColorBrewer)

theme_set(theme_tufte(base_size = 6, base_family = "ArialMT"))
theme_update(text = element_text(size = 6),
             axis.text = element_text(size = 6),
             axis.line = element_line(size = 0.5), 
             axis.ticks = element_line(size = 0.5),
             axis.ticks.length = unit(.1, "cm"))

BASE=dirname(this.path())
OUTFIGURES=file.path(BASE, "output")
SIG=file.path(BASE, "input/Signature_Compendium_v5_Cosine-0.74_Signatures_NAMESAPRIL21.rds")

## Load and normalise
sig = readRDS(SIG)
sigCols = apply(sig, 2, function(x) x/sum(x))

dtSigDef = data.table(melt(sigCols))
dtSigDef$Var1 = factor(as.character(dtSigDef$Var1), levels = rev(levels(dtSigDef$Var1)))

## Discretise weights
breaks = c(0, 10e-10, 0.01, 0.1, 0.2, 0.5, 0.75, 1)
# Names of categories
tags = c("0", "> 0", "> 0.01", "> 0.1", "> 0.2", "> 0.5", "> 0.75")
# Split p-vals into somewhat meaningful categories
dtSigDef$Disc = cut(dtSigDef$value, breaks = breaks, include.lowest=TRUE, right=FALSE, labels = tags)

discCols = c("0" = "#ffffff", "> 0" = "#ffffcc", "> 0.01" =  "#ffeda0", "> 0.1" = "#feb24c", 
             "> 0.2" = "#fd8d3c", "> 0.5" = "#f03b20", "> 0.75" = "#bd0026")


### Tilted vertical and elongated for A4 page with mixture components
dtSigDef$Var1 = factor(as.character(dtSigDef$Var1), levels = rev(levels(dtSigDef$Var1)))
dtSigDef$Var2 = factor(as.character(dtSigDef$Var2), levels = rev(levels(dtSigDef$Var2)))

pB = ggplot(dtSigDef, aes(y = Var2, x = Var1, fill = Disc)) + 
  geom_tile(aes(width = 0.94, height = 0.94)) + 
  theme_tufte(base_family = "", base_size = 18) + 
  theme(legend.position = "bottom", axis.line = element_line(size = 0.5), axis.ticks = element_line(size = 0.5), 
        axis.ticks.length = unit(.1, "cm"), plot.margin = unit(c(0, 0, 0, 0), "null"), aspect.ratio = 43/17) + 
  labs(x = "CIN signature", y = "Feature components") + 
  scale_x_discrete(position = "top", guide = guide_axis(angle = -90)) + labs(fill = "Weight") + 
  scale_fill_manual(values = discCols) + guides(fill = guide_legend(nrow = 1)) +
  coord_capped_cart(top = "both", left = "both")

## Save output
cairo_pdf(file.path(OUTFIGURES, "SuppFig_8_Signature_definition_vertical.pdf"), height = 230/25.4, width = 161/25.4)
print(pB); dev.off()

ggsave(file.path(OUTFIGURES, "SuppFig_8_Signature_definition_vertical.svg"), pB, height = 230/25.4, width = 161/25.4)

