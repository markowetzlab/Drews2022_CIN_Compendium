library(this.path)
library(data.table)
library(reshape2)
library(ggplot2)
library(ggthemes)
library(lemon)
library(RColorBrewer)

theme_set(theme_tufte(base_size = 6, base_family = "ArialMT"))
theme_update(text = element_text(size = 6),
             axis.text = element_text(size = 6))

BASE=dirname(this.path())
threshs = fread(file.path(BASE, "input/3_Boxplots_per_sig_fullTCGA_1000sims_10pGaussian_10pSamplePoisson.txt"))
OUTPUT=file.path(BASE, "output")

## Plot
pOut = ggplot(threshs, aes(x = Sig, y = Thresh_ZeroGMM_0.05)) + 
  geom_crossbar(aes(ymin = Thresh_ZeroGMM_0.05, ymax = Thresh_ZeroGMM_0.05), width = 0.6) + 
  scale_y_continuous(labels = scales::percent) + labs(x = "Signatures", y = "Activity threshold") +
  theme(axis.line = element_line(size = 1), axis.ticks = element_line(size = 1),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
  coord_capped_cart(bottom = "both", left = "both")

cairo_pdf(file.path(OUTPUT, "MC_Simulation_thresholds.pdf"), width = 88/25.4, height = 45/25.4)
print(pOut); dev.off()
ggsave(file.path(OUTPUT, "MC_Simulation_thresholds.svg"), plot = pOut, 
       width = 88/25.4, height = 45/25.4)
