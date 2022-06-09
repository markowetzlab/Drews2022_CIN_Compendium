#This code is copyright (c) 2022, University of Cambridge and Spanish National Cancer Research Centre (CNIO).
#This code is published and distributed under the GAP Available Source License v1.0 (ASL). 
#This code is distributed in the hope that it will be useful for non-commercial academic research, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the ASL for more details. 

## Identify the threshold of when detectable CIN is not detectable anymore

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


## Use normalise_df_per_dim function from YAPSA package to use the same code as in Macintyre et al.
repeat_df = function (in_value, in_rows, in_cols)
{
  return(as.data.frame(matrix(rep(in_value, in_cols * in_rows),
                              ncol = in_cols)))
}

normalize_df_per_dim = function (in_df, in_dimension)
{
  out_df <- repeat_df(in_value = 0, in_rows = dim(in_df)[1],
                      in_cols = dim(in_df)[2])
  if (in_dimension == 1) {
    choice_ind <- which(rowSums(in_df) > 0)
    out_df[choice_ind, ] <- in_df[choice_ind, ]/rowSums(in_df)[choice_ind]
  }
  else if (in_dimension == 2) {
    t_df <- t(in_df)
    choice_ind <- which(rowSums(t_df) > 0)
    temp_df <- repeat_df(in_value = 0, in_rows = dim(in_df)[2],
                         in_cols = dim(in_df)[1])
    temp_df[choice_ind, ] <- t_df[choice_ind, ]/rowSums(t_df)[choice_ind]
    out_df <- as.data.frame(t(temp_df))
  }
  else {
    return(NULL)
  }
  colnames(out_df) <- colnames(in_df)
  rownames(out_df) <- rownames(in_df)
  return(out_df)
}


## Paths
BASE=dirname(this.path())
ABSEXP=file.path(BASE, "input/Export-matrix_OV_Sigs_on_TCGA-OV_rawExposures_12112019.rds")
NORMEXP=file.path(BASE, "input/Export-matrix_OV_Sigs_on_TCGA-OV_12112019.rds")
OVSIGS=file.path(BASE, "input/external/Macintyre2018_OV_feat_sig_mat.rds")
DETECTIONLIMIT=0.05
META=file.path( BASE, "input/Metadata_TCGA_ASCAT_penalty70.rds" )
CLINICALDATA=file.path( BASE, "input/external/metadata_CNA_12K.RDS" )
OUTFIGURE=file.path(BASE, "output")
OUTTABLE=file.path(BASE, "output")
dir.create(OUTFIGURE, showWarnings = FALSE, recursive = TRUE)
dir.create(OUTTABLE, showWarnings = FALSE, recursive = TRUE)

DIAGNOSIS="Serous cystadenocarcinoma, NOS"
PROJECT="TCGA-OV"

## Get samples with distinct HGSOC diagnosis
clinMetadata = readRDS(CLINICALDATA)
sampleIDs = clinMetadata$submitter_id[ clinMetadata$project_id == PROJECT & clinMetadata$primary_diagnosis == DIAGNOSIS ]

meta = readRDS(META)
metaFilt = meta[ ! is.na(meta$CNAs), ]
finalNames = metaFilt$name[ metaFilt$patient %in% sampleIDs ]


## Load exposures
absHgsoc = readRDS(ABSEXP)
absHgsoc = absHgsoc[ rownames(absHgsoc) %in% finalNames, ]

relHgsoc = readRDS(NORMEXP)
relHgsoc = relHgsoc[ rownames(relHgsoc) %in% finalNames, ]

ovSigs = readRDS(OVSIGS)

# Set everything we can detect (ie over the detection limit) to 0 and then calculate back the SxC matrix
absExpLim = absHgsoc
absExpLim[ relHgsoc > DETECTIONLIMIT ] = 0
notDetectSxC = absExpLim %*% t(normalize_df_per_dim(ovSigs, 2))


# Add elements up per feature
dfNotDetect = melt(notDetectSxC)
dfNotDetect$Var2 = gsub('[[:digit:]]+', '', dfNotDetect$Var2)
dtSummedNotDetect = aggregate(value ~ Var1 + Var2, dfNotDetect, sum)

dtSummedNotDetect = dtSummedNotDetect[ dtSummedNotDetect$value != 0, ]
p1 = ggplot(dtSummedNotDetect, aes( x = value)) + geom_histogram() + 
  facet_wrap( . ~ Var2, scales = "free") +
  coord_capped_cart(bottom = "both", left = "both") +
  labs(y = "Number of samples", x = "Estimated number of unassigned CNAs")

# 95% quantile of feature distribution is our threshold
lThresh = lapply(unique(dtSummedNotDetect$Var2), function(feature) {
    THRESHOLD = as.numeric( ceiling( quantile( dtSummedNotDetect$value[ dtSummedNotDetect$Var2 == feature ], 0.95 ) ) )    
    return(paste0(feature, ": ", THRESHOLD))
} )
dtThresh = do.call(rbind, lThresh)

# Save outputs
cairo_pdf(file.path(OUTFIGURE, "EDF_1_B_Not_detectable_CNA_signal.pdf"), width = 3, height = 2)
print(p1); dev.off()

ggsave(filename = file.path(OUTFIGURE, "EDF_1_B_Not_detectable_CNA_signal.svg"), p1, width = 3, height = 2)

write.table(dtThresh, file.path(OUTTABLE, paste0("EDF_1_B_Threshold_detectable_CIN_Limit_", DETECTIONLIMIT, ".txt")), 
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

