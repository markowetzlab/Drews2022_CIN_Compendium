#This code is copyright (c) 2022, University of Cambridge and Spanish National Cancer Research Centre (CNIO).
#This code is published and distributed under the GAP Available Source License v1.0 (ASL). 
#This code is distributed in the hope that it will be useful for non-commercial academic research, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the ASL for more details. 

# Creating Extended Data Fig. 1, panels c-e

library(this.path)
library(data.table)
library(reshape2)
library(ggplot2)
library(ggthemes)
library(patchwork)
library(lemon)
library(RColorBrewer)

theme_set(theme_tufte(base_size = 6.5, base_family = "ArialMT"))
theme_update(text = element_text(size = 6.5),
             axis.text = element_text(size = 6.5),
             axis.line = element_line(size = 0.5), 
             axis.ticks = element_line(size = 0.5), axis.ticks.length = unit(.1, "cm"), 
             plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "null"))
  

## Paths
BASE=dirname(this.path())
METAFULL=file.path(BASE, "input/Metadata_TCGA_ASCAT_penalty70.rds")
COLS = file.path(BASE, "input/TCGA_colour_scheme.txt")
OUT=file.path(BASE, "output")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

WEAVERCLEVELAND2006=file.path(BASE, "input/external/WeaverCleveland2006_Table1.txt")
TAYLOR2018=file.path(BASE, "input/external/Taylor2018_S2.csv")
DAVOLI2017=file.path(BASE, "input/external/Davoli2017_Supplement_from_author.txt")

## Load data
meta = readRDS(METAFULL)
wc = fread(WEAVERCLEVELAND2006)
taylor = fread(TAYLOR2018)
davoli = fread(DAVOLI2017)

# Prepare colours
cols = read.table(COLS, comment.char = "'")
vCols = as.character(cols$V2)
names(vCols) = cols$V1


#### Figure 1: dCIN vs Weaver/Cleveland 2006 (Mitelman database and origin of the 90% aneuploid statement)
meta = meta[ ! is.na(meta$dCIN), ]
meta$cancer_type[ grepl("BRCA", meta$cancer_type) ] = "BRCA"
dtDCIN = aggregate(dCIN ~ cancer_type, meta, mean)
dtDCIN$Weaver = wc$dCIN[ match(dtDCIN$cancer_type, wc$TCGA) ]
dtDCIN = dtDCIN[ ! is.na(dtDCIN$Weaver),]
dtDCIN$cancer_type = factor(dtDCIN$cancer_type)

p0 = ggplot(dtDCIN, aes(y = Weaver, x = dCIN, colour = cancer_type)) + geom_abline(intercept = 0, slope = 1) + 
  geom_point(size = 2) + labs(y = "Proportion of CIN according to Weaver and Cleveland (2006, Mitelman database)",
                              x = "Proportion of patients with CIN based on detectable CIN (dCIN)",
                              colour = "Best matched\nTCGA code") +
  scale_colour_manual(values = vCols) + scale_y_continuous(labels = scales::percent) + 
  geom_smooth(method='lm', colour = "red", linetype = "dashed", se = FALSE) +
  scale_x_continuous(labels = scales::percent) + 
  coord_capped_cart(bottom = "both", left = "both", xlim = c(0.2, 1), ylim = c(0.2, 1)) +
  theme(legend.key.size = unit(0.25, "cm"), legend.position = c(0.2, 0.6))

cairo_pdf(file.path(OUT, "EDF_1_C_Weaver_Cleveland_vs_dCIN.pdf"), width = 82.6/25.4, height = 82.6/25.4)
print(p0); dev.off()
ggsave(file.path(OUT, "EDF_1_C_Weaver_Cleveland_vs_dCIN.svg"), p0, width = 82.6, height = 82.6, units = "mm")



#### Figure 2: dCIN vs Taylor et al. 2018
meta$AS = taylor$`AneuploidyScore(AS)`[ match(substr(meta$name,1,12), substr(taylor$Sample, 1, 12) ) ]

t.test(meta$AS[meta$dCIN], meta$AS[! meta$dCIN], var.equal = FALSE)

p1 = ggplot(meta, aes(x = dCIN, y = AS)) + geom_boxplot(outlier.colour = NA) + geom_jitter(size = 0.3 , height = 0, width = 0.2, alpha = 0.1) +
  labs(x = "Detectable CIN (dCIN)", y = "Aneuploidy Score (AS) from Taylor et al. (2018)") +
  coord_capped_cart(bottom = "both", left = "both")


cairo_pdf(file.path(OUT, "EDF_1_D_Taylor_AS_vs_dCIN.pdf"), width = 2, height = 3)
print(p1); dev.off()
ggsave(file.path(OUT, "EDF_1_DTaylor_AS_vs_dCIN.svg"), p1, width = 2, height = 3)



#### Figure 3: dCIN vs Davoli et al. 2017
meta$DavChrom = davoli$Unweighted.chrom.value.scale[ match(substr(meta$name,1,12), substr(davoli$TumorSample, 1, 12) ) ]
meta$DavArm = davoli$Unweighted.arm.value.scale[ match(substr(meta$name,1,12), substr(davoli$TumorSample, 1, 12) ) ]
meta$DavFocal = davoli$Unweighted.foc.value.scale[ match(substr(meta$name,1,12), substr(davoli$TumorSample, 1, 12) ) ]
meta$DavSCNA = davoli$Unweighted.SCNA.value.scale[ match(substr(meta$name,1,12), substr(davoli$TumorSample, 1, 12) ) ]

dtMetaDav = meta[ ! is.na(meta$DavSCNA), c("name", "dCIN", "DavChrom", "DavArm", "DavFocal", "DavSCNA")]
mMetaDav = melt(dtMetaDav, id.vars = c("name", "dCIN"), measure.vars = c("DavChrom", "DavArm", "DavFocal", "DavSCNA"))
mMetaDav$variable = factor(mMetaDav$variable, levels = c("DavSCNA", "DavChrom", "DavArm", "DavFocal"), 
                           labels = c("SCNA score", "Whole Chrom. score", "Arm score", "Focal score"))

t.test(mMetaDav$value[ mMetaDav$dCIN & mMetaDav$variable == "SCNA score" ], mMetaDav$value[ ! mMetaDav$dCIN & mMetaDav$variable == "SCNA score" ], var.equal = FALSE)
t.test(mMetaDav$value[ mMetaDav$dCIN & mMetaDav$variable == "Whole Chrom. score" ], mMetaDav$value[ ! mMetaDav$dCIN & mMetaDav$variable == "Whole Chrom. score" ], var.equal = FALSE)
t.test(mMetaDav$value[ mMetaDav$dCIN & mMetaDav$variable == "Arm score" ], mMetaDav$value[ ! mMetaDav$dCIN & mMetaDav$variable == "Arm score" ], var.equal = FALSE)
t.test(mMetaDav$value[ mMetaDav$dCIN & mMetaDav$variable == "Focal score" ], mMetaDav$value[ ! mMetaDav$dCIN & mMetaDav$variable == "Focal score" ], var.equal = FALSE)

p2 = ggplot(mMetaDav, aes(x = dCIN, y = value)) + geom_boxplot(outlier.size = 0.5) + facet_grid(. ~ variable) + 
  geom_hline(yintercept = 0, colour = "grey80", linetype = "dashed") + coord_capped_cart(bottom = "both", left = "both", ylim = c(-5, 8)) +
  labs(x = "Detectable CIN (dCIN)", y = "SCNA level from Davoli et al. (2017)")

cairo_pdf(file.path(OUT, "EDF_1_E_Davoli_vs_dCIN.pdf"), width = 137/25.4, height = 3)
print(p2); dev.off()
ggsave(file.path(OUT, "EDF_1_E_Davoli_vs_dCIN.svg"), p2, width = 137/25.4, height = 3)
