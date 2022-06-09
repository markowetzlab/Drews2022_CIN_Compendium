## Plot estimated number of CNAs vs signatures and cancers

rm(list=ls(all=TRUE))

library(this.path)
library(data.table)
library(reshape2)
library(ggplot2)
library(ggthemes)
library(lemon)
library(RColorBrewer)
library(dplyr)
library(ggforce) # for 'geom_arc_bar'
library(gglorenz)

theme_set(theme_tufte(base_size = 6, base_family = "ArialMT"))
theme_update(text = element_text(size = 6),
             axis.text = element_text(size = 6),
             axis.line = element_line(size = 0.5), 
             axis.ticks = element_line(size = 0.5),
             axis.ticks.length = unit(.1, "cm"))

## Basics
BASE=dirname(this.path())
OUTPLOTS=file.path(BASE, "output")
OUTRESULTS=file.path(BASE, "output")

## Signature files
RAW=file.path(BASE, "input/Signature_Compendium_v5_Cosine-0.74_Activities_raw_THRESH95_NAMESAPRIL21.rds")
META=file.path(BASE, "input/Metadata_TCGA_ASCAT_penalty70.rds")


## Load data
raw = readRDS(RAW)
meta = readRDS(META)


## Estimate number of CNAs
dfRaw = as.data.frame(raw)
dfRaw$CNAs = meta$CNAs[ match(rownames(dfRaw), meta$name) ]

lmCNAs = lm(CNAs ~ . + 0, dfRaw)

saveRDS(lmCNAs, file.path(OUTRESULTS, "Estimated_CNAs_per_signature_and_sample_MODEL.rds"))

est = raw %*% diag(coefficients(lmCNAs))
colnames(est) = colnames(raw)

# Save as output file
estOut = signif(est, 4)
saveRDS(est, file.path(OUTRESULTS, "SuppMat_Estimated_CNAs_per_signature_and_sample.rds"))
write.table(estOut, file.path(OUTRESULTS, "SuppMat_Estimated_CNAs_per_signature_and_sample.txt"), 
            row.names = TRUE, col.names = TRUE, sep = "\t", quote = FALSE)

# Melt and transfer cancer type
dtEst = data.table(melt(est))
dtEst$Cancer = meta$cancer_type[ match(dtEst$Var1, meta$name) ]

dtSum = aggregate(value ~ Cancer, dtEst, sum)
dtSumSig = aggregate(value ~ Cancer + Var2, dtEst, sum)

dtSumSig$all = dtSum$value[ match(dtSumSig$Cancer, dtSum$Cancer) ]
dtSumSig$Prop = dtSumSig$value / dtSumSig$all

# Discretise 
BREAKS = c(0, 0.01, 0.025, 0.05, 0.1, 0.2, 0.3, 1)
LABELS = c("0%", ">1%", ">2.5%", ">5%", ">10%", ">20%", ">30%" )
dtSumSig$DiscExp = cut(dtSumSig$Prop, breaks = BREAKS, labels = LABELS, include.lowest = TRUE)

dtSumSig$Var2 = factor(as.character(dtSumSig$Var2), levels = rev(levels(dtSumSig$Var2)))

pA = ggplot(dtSumSig, aes(x = Cancer, y = Var2, fill = DiscExp)) + 
  geom_tile(aes(width = 0.94, height = 0.94)) + 
  labs(x = "TCGA cancer type", y = "CIN signatures", fill = "Estimated number\nof CNAs") +
  theme(legend.position = "bottom") +
  scale_fill_manual(values = c("white", "#ffffcc", "#fed976", "#fd8d3c", 
                               "#fc4e2a", "#b10026", "#67000d")) +
  theme(aspect.ratio = 17/33, legend.key.width = unit(1, "line"), legend.key.height = unit(0.4, "line"),
        axis.line = element_line(size = 1), axis.text.x = element_text(angle = 270, vjust = 0.5, hjust = 1)) +
  scale_x_discrete(position = "top") + coord_capped_cart(top = "both", left = "both") 

cairo_pdf(file.path(OUTPLOTS, "SuppFig_38_A_CNA_heatmap.pdf"), width = 90/25.4, height = 75/25.4)
print(pA); dev.off()

ggsave(file.path(OUTPLOTS, "SuppFig_38_A_CNA_heatmap.svg"), pA, width = 90/25.4, height = 75/25.4)


## Panel B
# How many CNAs across TCGA cohort?
dfEstPie = data.frame(colSums(est))
dfEstPie$Sig = factor(rownames(dfEstPie), levels = rownames(dfEstPie)[ order(dfEstPie$colSums.est.) ])
colnames(dfEstPie) = c("Val", "Sig")

dfEstPie = dfEstPie[ levels(dfEstPie$Sig), ]
dfEstPie$Label = paste0(dfEstPie$Sig, " ", signif(dfEstPie$Val/sum(dfEstPie$Val)*100, 2), "%")

## From: https://stackoverflow.com/questions/48184645/how-can-i-put-the-labels-outside-of-piechart
dfEstPie <- dfEstPie %>% 
  mutate(end = 2 * pi * cumsum(Val)/sum(Val),
         start = lag(end, default = 0),
         middle = 0.5 * (start + end),
         hjust = ifelse(middle > pi, 1, 0),
         vjust = ifelse(middle < pi/2 | middle > 3 * pi/2, 0, 1))

## Add colour
BREAKS = c(0, 0.01, 0.025, 0.05, 0.1, 0.2, 0.3, 1)
LABELS = c("0%", ">1%", ">2.5%", ">5%", ">10%", ">20%", ">30%" )
dfEstPie$DiscExp = cut(dfEstPie$Val/sum(dfEstPie$Val), breaks = BREAKS, labels = LABELS, 
                       include.lowest = TRUE)

pB = ggplot(dfEstPie, aes(fill = DiscExp)) + 
  geom_arc_bar(aes(x0 = 0, y0 = 0, r0 = 0, r = 1,
                   start = start, end = end, group = Sig)) +
  geom_text(aes(x = 1 * sin(middle), y = 1 * cos(middle), label = Label,
                hjust = hjust, vjust = vjust)) +
  coord_fixed() + theme(legend.position = "bottom") +
  scale_x_continuous(limits = c(-1.5, 1.5),  # Adjust so labels are not cut off
                     name = "", breaks = NULL, labels = NULL) +
  scale_y_continuous(limits = c(-1, 1.1),      # Adjust so labels are not cut off
                     name = "", breaks = NULL, labels = NULL) +
  scale_fill_manual(values = c("white", "#ffffcc", "#fed976", "#fd8d3c", 
                             "#fc4e2a", "#b10026", "#67000d"))

cairo_pdf(file.path(OUTPLOTS, "SuppFig_38_B_CNA_Pie.pdf"), width = 90/25.4, height = 75/25.4)
print(pB); dev.off()

ggsave(file.path(OUTPLOTS, "SuppFig_38_B_CNA_Pie.svg"), pB, width = 90/25.4, height = 75/25.4)



## Panel C
# Distribution of CNAs across samples with this signature
# Define the number of colors you want
numSigs = 17
myCols = colorRampPalette(brewer.pal(8, "Set2"))(numSigs)

myCols = c("#ff8f8d","#34b331","#de0082","#00ad7c","#d3005c","#c7cc76","#a591ff","#f7b810","#0283dd",
            "#ff724d","#0172a0","#7e7a00","#ff9af3","#714a22","#94bfff","#982a35","#665889")

pC = ggplot(dtEst, aes(colour = Var2)) + stat_lorenz(aes(value)) + geom_abline(color = "grey") +
  scale_colour_manual(values = myCols) + scale_x_continuous(labels = scales::percent) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "Proportion of samples", y = "Proportion of CNA", colour = "Signature") +
  theme(legend.position = c(0.1, 0.7)) +
  theme(legend.key.height = unit(0.4, "line"),
        axis.line = element_line(size = 1)) +
  coord_capped_cart(bottom = "both", left = "both") 

cairo_pdf(file.path(OUTPLOTS, "SuppFig_38_C_Lorenz.pdf"), width = 90/25.4, height = 75/25.4)
print(pC); dev.off()

ggsave(file.path(OUTPLOTS, "SuppFig_38_C_Lorenz.svg"), pC, width = 90/25.4, height = 75/25.4)

