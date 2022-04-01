# info --------------------------------------------------------------------
# purpose: plot averages
# author: Robert Gnuegge (robert.gnuegge@gmail.com)
# date: 02/08/22
# version: 1.0


# preamble ----------------------------------------------------------------

# set working directory to this file's location
wd.path <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(wd.path)

# directories
plot_dir <- "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Figures/15_qPCR_Resection_assays/LSY5023-98C_5038-9C/Average"
data_dir <- "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Data/qPCR_Resection_assays/LSY5023-98C_5038-9C/Average"

# read in helper files and functions
source(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Analysis_scripts/Misc/Misc_helper_functions.R")
source(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Analysis_scripts/Misc/JFly_colors.R")

# read data
PHO5.cut.modes <- read.table(file = paste0(data_dir, "/Cutting_PHO5.txt"), header = TRUE)
PHO5.resect.modes <- read.table(file = paste0(data_dir, "/Resection_PHO5.txt"), header = TRUE)


# Plot PHO5 cutting -------------------------------------------------------
strains <- c("5038-9C", "5023-98C")
MyColors <- gray(level = c(0.75, 0))

tmp <- PHO5.cut.modes

# print to PDF
pdf(file = "tmp.pdf", width=2.75, height=2.5)
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(2.1, 1.0, 4.0, 2.0), tcl = -0.25, mgp = c(2.25, 0.5, 0), las = 1)

# start empty plot
y.range <- range(c(tmp$mean + tmp$sd, tmp$mean - tmp$sd, 0, 1))
# y.range <- c(-0.5, 1.5)
plot(x = NA, y = NA, ylim = y.range, xlim = range(tmp$time), xlab = NA, ylab = "Cut Fraction")
title(xlab =  "Time [h]", line = 2)

# iterate through all samples
for (n in 1:length(strains)){
  
  strain <- strains[n]
  
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
}

legend(x = "bottomright", legend = c(expression(italic("pho4"*Delta)), expression(italic("pho4-SA1234PA6"))), 
       col = MyColors, pch = 20, inset = 0.0, bty = "n")

dev.off()

GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/Cut_fraction_PHO5.pdf"))


# plot PHO5 resection -----------------------------------------------------
strains <- c("5038-9C", "5023-98C")
MyColors <- gray(level = c(0.75, 0))

# print to PDF
pdf(file = "tmp.pdf", width=2.75, height=2.5)
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(2.1, 1.0, 4.0, 2.0), tcl = -0.25, mgp = c(2.25, 0.5, 0), las = 1)

# -223 bp
tmp <- PHO5.resect.modes[PHO5.resect.modes$primers == "oRG802_oRG803", ]

# start empty plot
y.range <- range(c(tmp$mean + tmp$sd, tmp$mean - tmp$sd), na.rm = TRUE)
y.range <- c(0, 0.8)
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
         pch = 21, col = MyColors[n], type = "o", lty = "dashed")
  
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

# -538 bp
tmp <- PHO5.resect.modes[PHO5.resect.modes$primers == "oRG806_oRG807", ]

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

legend(x = "topleft", 
       legend = c(expression(italic("pho4-SA1234PA6")), expression(italic("pho4"*Delta))), 
       col = MyColors[2:1], pch = 20, bty = "n", pt.cex = 0.75)

legend(x = "topleft", 
       legend = c(NA, NA), 
       col = MyColors[2:1], pch = 21, bty = "n", pt.cex = 1.25)

dev.off()

GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/Resection_PHO5.pdf"))