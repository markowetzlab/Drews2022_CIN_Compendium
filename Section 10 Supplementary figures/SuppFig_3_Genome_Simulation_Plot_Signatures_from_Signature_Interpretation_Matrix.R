# Plot copy number profiles created with CIN Genome Simulation 
## (see CINGenomeSimulation repo for more details on how to create the profiles)

rm(list=ls(all=TRUE))

## Libraries
library(this.path)
library(data.table)
library(reshape2)
library(ggplot2)
library(ggthemes)
library(lemon)

## Paths
BASE=dirname(this.path())
OUT=file.path(BASE, "output")

# Load relevant data and functions
source(file.path(BASE, "input/0_Simulation_Cancer_Genomes_Functions.R"))
lDFStep3=readRDS(file.path(BASE, "input/Out_2_Step_3_Segmentation_tables_20each_N240.rds"))
fileChrSizes = file.path(BASE, "input/hg19.chrom.sizes.txt")

dfChrSizes = loadChromSizes(fileChrSizes, rmChr = TRUE)

## Plot sample profiles for all combinations
lStep3Plots = plotAllSamples(lDFStep3, dfChrSizes)
examplePlot3 = wrap_plots(lStep3Plots[c(1, 101, 122, 231)], ncol = 2)
ggsave(file.path(OUT, "SuppFig_3_Example_plots_Genomic_Simulation.png"), examplePlot3, scale = 1.5, units = "cm", width = 16, height = 12,
       bg = "white")
ggsave(file.path(OUT, "SuppFig_3_Example_plots_Genomic_Simulation.svg"), examplePlot3, scale = 1.5, units = "cm", width = 16, height = 12,
       bg = "white")
