# info --------------------------------------------------------------------
# purpose: plot averages
# author: Robert Gnuegge (robert.gnuegge@gmail.com)
# date: 02/22/22
# version: 1.0


# preamble ----------------------------------------------------------------

# set working directory to this file's location
wd.path <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(wd.path)

# directories
plot_dir <- "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Figures/15_qPCR_Resection_assays/LSY4741-8C_4815-25B/Average"
data_dir <- "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Data/qPCR_Resection_assays/LSY4741-8C_4815-25B/Average"

# read in helper files and functions
source(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Analysis_scripts/Misc/Misc_helper_functions.R")
source(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Analysis_scripts/Misc/JFly_colors.R")

# read data
Cut_modes <- read.table(file = paste0(data_dir, "/Cut_fraction.txt"), header = TRUE)
Resect_modes <- read.table(file = paste0(data_dir, "/Resected_fraction.txt"), header = TRUE)


# Plot cutting ============================================================
strains <- c("4741-8C", "4815-25B")
MyColors <- JFly_colors[1:length(strains)]

tmp <- Cut_modes

# print to PDF
pdf(file = "tmp.pdf", width=2.75, height=2.5)
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(2.1, 1.0, 4.0, 2.0), tcl = -0.25, mgp = c(2.25, 0.5, 0), las = 1)

# start empty plot
# y.range <- range(c(tmp$mean + tmp$sd, tmp$mean - tmp$sd, 0, 1))
y.range <- c(-0.39, 1)
plot(x = NA, y = NA, ylim = y.range, xlim = range(tmp$time), xlab = NA, ylab = "Cut Fraction")
title(xlab =  "Time [h]", line = 2)

# iterate through all samples
for (n in 1:length(strains)){
  strain <- strains[n]
  
  # add HOcs cutting
  tmp <- Cut_modes[Cut_modes$primers == "oRG46_oRG47", ]
  
  # time values with jitter
  t <- jitter(x = tmp$time[tmp$strain == strain], factor = 0)
  
  # line plot of means
  points(x = t, 
         y = tmp$mean[tmp$strain == strain],
         pch = 20, col = MyColors[n], type = "o"
  )
  
  # add error bars
  arrows(x0 = t,
         y0 = tmp$mean[tmp$strain == strain] - tmp$sd[tmp$strain == strain],
         x1 = t,
         y1 = tmp$mean[tmp$strain == strain] + tmp$sd[tmp$strain == strain],
         length = 0.03, # length of arrow head
         angle = 90, # angle of arrow head
         code = 3, # to draw arrow head on both ends
         col = MyColors[n])
  
  # add gSrfIcs cutting
  tmp <- Cut_modes[Cut_modes$primers == "oRG341_oRG342", ]
  
  # time values with jitter
  t <- jitter(x = tmp$time[tmp$strain == strain], factor = 0)
  
  # line plot of means
  points(x = t, 
         y = tmp$mean[tmp$strain == strain],
         pch = 21, col = MyColors[n], type = "o", lty = "dashed"
  )
  
  # add error bars
  arrows(x0 = t,
         y0 = tmp$mean[tmp$strain == strain] - tmp$sd[tmp$strain == strain],
         x1 = t,
         y1 = tmp$mean[tmp$strain == strain] + tmp$sd[tmp$strain == strain],
         length = 0.03, # length of arrow head
         angle = 90, # angle of arrow head
         code = 3, # to draw arrow head on both ends
         col = MyColors[n])

}

dev.off()

GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/Cut_fraction.pdf"))


# plot legend =================================================================
LegTxt <- c(expression(italic("lexO-HO")), expression(italic("lexO-HO lexO-SrfI")))

pdf(file = "tmp.pdf", width=1.85, height=0.65)
par(cex = 1, mar = rep(0, 4))
plot(1, type="n", axes=FALSE, xlab="", ylab="")
legend(1, 1, xjust=0.5, yjust=0.5, legend = LegTxt, 
       col = MyColors, pch = 20)
dev.off()
GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/Legend.pdf"))


# plot resection ==========================================================
strains <- c("4741-8C", "4815-25B")
MyColors <- JFly_colors[1:length(strains)]


# +98 bp ------------------------------------------------------------------
tmp <- Resect_modes[Resect_modes$primers == "oRG50_oRG51", ]

# print to PDF
pdf(file = "tmp.pdf", width=2.75, height=2.5)
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(2.1, 1.0, 4.0, 2.0), tcl = -0.25, mgp = c(2.25, 0.5, 0), las = 1)

# start empty plot
# y.range <- range(c(tmp$mean + tmp$sd, tmp$mean - tmp$sd, 0, 1))
y.range <- c(-0.39, 1)
plot(x = NA, y = NA, ylim = y.range, xlim = range(tmp$time), xlab = NA, ylab = "ssDNA Fraction")
title(xlab =  "Time [h]", line = 2)

# iterate through all samples
for (n in 1:length(strains)){
  
  strain <- strains[n]
  
  # time values with jitter
  t <- jitter(x = tmp$time[tmp$strain == strain], factor = 0.2)
  
  # line plot of means
  points(x = t, 
         y = tmp$mean[tmp$strain == strain],
         pch = 20, col = MyColors[n], type = "o")
  
  # add error bars
  arrows(x0 = t,
         y0 = tmp$mean[tmp$strain == strain] - tmp$sd[tmp$strain == strain],
         x1 = t,
         y1 = tmp$mean[tmp$strain == strain] + tmp$sd[tmp$strain == strain],
         length = 0.03, # length of arrow head
         angle = 90, # angle of arrow head
         code = 3, # to draw arrow head on both ends
         col = MyColors[n])
}

dev.off()

GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/Resection_98bp.pdf"))


# +640 bp -----------------------------------------------------------------
tmp <- Resect_modes[Resect_modes$primers == "oRG52_oRG53", ]

# print to PDF
pdf(file = "tmp.pdf", width=2.75, height=2.5)
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(2.1, 1.0, 4.0, 2.0), tcl = -0.25, mgp = c(2.25, 0.5, 0), las = 1)

# start empty plot
# y.range <- range(c(tmp$mean + tmp$sd, tmp$mean - tmp$sd, 0, 1))
y.range <- c(-0.39, 1)
plot(x = NA, y = NA, ylim = y.range, xlim = range(tmp$time), xlab = NA, ylab = "ssDNA Fraction")
title(xlab =  "Time [h]", line = 2)

# iterate through all samples
for (n in 1:length(strains)){
  
  strain <- strains[n]
  
  # time values with jitter
  t <- jitter(x = tmp$time[tmp$strain == strain], factor = 0.2)
  
  # line plot of means
  points(x = t, 
         y = tmp$mean[tmp$strain == strain],
         pch = 20, col = MyColors[n], type = "o")
  
  # add error bars
  arrows(x0 = t,
         y0 = tmp$mean[tmp$strain == strain] - tmp$sd[tmp$strain == strain],
         x1 = t,
         y1 = tmp$mean[tmp$strain == strain] + tmp$sd[tmp$strain == strain],
         length = 0.03, # length of arrow head
         angle = 90, # angle of arrow head
         code = 3, # to draw arrow head on both ends
         col = MyColors[n])
}

dev.off()

GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/Resection_640bp.pdf"))