# info --------------------------------------------------------------------
# purpose: plot density plots for flow cytometry data
# author: Robert Gnuegge (robert.gnuegge@gmail.com)
# date: 02/25/22
# version: 1.0


# preamble ----------------------------------------------------------------

# set wd to this file's location
wd_path <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(wd_path)

# load helper files and functions
source(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Analysis_scripts/Misc/Misc_helper_functions.R")


# read extracted fluorescence data --------------------------------------------
df <- read.table(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Data/Flow_cytometry/LSY4753-3B_4754-3A/Extracted_Citrine_data.txt", header = TRUE, stringsAsFactors = FALSE)

# plot density distributions for all samples ----------------------------------
plot_dir <- "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Figures/19_Flow_cytometry/LSY4753-3B_4754-3A/Density_Plots"

for (gal in c(0.01, 2)){

  for (strain in c("4436", "4518", "4753", "4754")){

    for (t in c(0, 2, 6, 4, 8)){
      
      pdf(file = paste0(plot_dir, "/Gal_", sub(pattern = "\\.", replacement = "_", x = gal), "_",strain, "_", t, "h.pdf"), width = 2, height = 2)
        par(mar = rep(0,4))
        d <- density(x = log(x = df$Citrine[df$gal == gal & df$strain == strain & df$time == t], base = 10))
        plot(d, xlim = c(1.25, 4.75), axes = FALSE, xlab = NA, ylab = NA, main = NA, col = "transparent")
        polygon(d, col = "gray")
        points(d, type = "l")
      dev.off()
      
    }
  
  }
    
}