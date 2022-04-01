# info --------------------------------------------------------------------
# purpose: evaluate meltability profile of preferred nt sequence
# author: Robert Gnuegge (robert.gnuegge@gmail.com)
# date: 01/14/22
# version: 1.0

# preamble ----------------------------------------------------------------

# set wd to this file's location
wd_path <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(wd_path)

# libraries
library(rmelting)

# load helper files and functions
source(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Analysis_scripts/Misc/Misc_helper_functions.R")
source(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Analysis_scripts/Misc/JFly_colors.R")

plot_dir <- "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Figures/09_Meltability_impact/LSY4377-12B_4377-15A"

# function definitions ----------------------------------------------------

# calculate melting profile along a nt sequence
# argument: character string (nt sequence)
# result: numeric vector
melting_profile <- function(sequence, width){
  out <- rep(NA, nchar(sequence))
  n <- 3
  for(n in ceiling(width / 2):(nchar(sequence) - floor(width / 2))){
    nt <- toupper(substr(x = sequence, start = n - ceiling(width / 2) + 1, stop = n + floor(width / 2)))
    tmp <- melting(sequence = nt, nucleic.acid.conc = 1e-9, Na.conc = 0.001, K.conc = 0.001, Mg.conc = 0.005, hybridisation.type = "dnadna")
    tmp <- unlist(tmp$Results)[3:5]
    out[n] <- tmp["Melting temperature (C)"]
  }
  names(out) <- strsplit(sequence, "")[[1]]
  return(out)
}


# calculate meltability profile of preferred and bg sequences -------------
nicks <- "TTTTTTTTTTTTTTTTATTTTTTAATTTTTTTTTTAAAATTTTAAATATACTATATTTAAAAATTTAATTATTTTTTTTATAAATTTATAATTTTTTTTT"
bg <- "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"

profile_nicks_5bp <- melting_profile(sequence = nicks, width = 5)
profile_nicks_10bp <- melting_profile(sequence = nicks, width = 10)
profile_nicks_15bp <- melting_profile(sequence = nicks, width = 15)

profile_bg_5bp <- melting_profile(sequence = bg, width = 5)
profile_bg_10bp <- melting_profile(sequence = bg, width = 10)
profile_bg_15bp <- melting_profile(sequence = bg, width = 15)


# plot profiles -----------------------------------------------------------

# plotting function
plot_melting_profile <- function(nicks, bg){
  plot(x = 1:100, y = nicks, type = "l", ylim = range(c(nicks, bg), na.rm = TRUE),
       xlab = NA, ylab = expression("Melting Temperature ["*degree*"C"*"]"), xaxt = "n")
  points(x = 1:100, y = bg, type = "l", lty = "dashed")
  title(xlab = "Distance from Nick Site [nt]", line = 2.5)
  
  axis(side = 1, at = 1:100, labels = FALSE, tcl = -0.2)
  axis(side = 1, at = c(0:5 * 10, 5:9 * 10 + 1), labels = c(NA, -4:-1 * 10, NA, NA, 1:4 * 10), tcl = -0.3)
  axis(side = 1, at = 50.5, labels = 0, tick = FALSE)
  
  abline(v = 50.5, col = "gray")
}

# zoomed in plotting function
plot_melting_profile_zoom <- function(nicks, bg){
  plot(x = 1:41, y = nicks[31:71], type = "l", ylim = range(c(nicks[31:71], bg[31:71]), na.rm = TRUE),
       xlab = NA, ylab = expression("Melting Temperature ["*degree*"C"*"]"), xaxt = "n")
  points(x = 1:41, y = bg[31:71], type = "l", lty = "dashed")
  title(xlab = "Distance from Nick Site [nt]", line = 2.5)
  
  axis(side = 1, at = 1:41, labels = FALSE, tcl = -0.2)
  axis(side = 1, at = c(0:2 * 10, 2:4 * 10 + 1), labels = c(-2:-1 * 10, NA, NA, 1:2 * 10), tcl = -0.3)
  axis(side = 1, at = 20.5, labels = 0, tick = FALSE)
  
  abline(v = 20.5, col = "gray")
}

# 15 bp
pdf(file = "tmp.pdf", width = 5.5, height = 3)
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(1.5, 0.1, 4, 1.7), tcl = -0.3, mgp = c(2.75, 0.6, 0), las = 1)
  plot_melting_profile(nicks = profile_nicks_15bp, bg = profile_bg_15bp)
dev.off()
GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/Preferred_Nt_sequence_Melting_profile_15bp.pdf"))

pdf(file = "tmp.pdf", width = 3.5, height = 3)
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(1.5, 0.1, 4, 1.7), tcl = -0.3, mgp = c(2.75, 0.6, 0), las = 1)
  plot_melting_profile_zoom(nicks = profile_nicks_15bp, bg = profile_bg_15bp)
dev.off()
GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/Preferred_Nt_sequence_Melting_profile_15bp_zoom.pdf"))


# 10 bp
pdf(file = "tmp.pdf", width = 5.5, height = 3)
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(1.5, 0.1, 4, 1.7), tcl = -0.3, mgp = c(2.75, 0.6, 0), las = 1)
  plot_melting_profile(nicks = profile_nicks_10bp, bg = profile_bg_10bp)
dev.off()
GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/Preferred_Nt_sequence_Melting_profile_10bp.pdf"))

pdf(file = "tmp.pdf", width = 3.5, height = 3)
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(1.5, 0.1, 4, 1.7), tcl = -0.3, mgp = c(2.75, 0.6, 0), las = 1)
  plot_melting_profile_zoom(nicks = profile_nicks_10bp, bg = profile_bg_10bp)
dev.off()
GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/Preferred_Nt_sequence_Melting_profile_10bp_zoom.pdf"))


# 5 bp
pdf(file = "tmp.pdf", width = 5.5, height = 3)
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(1.5, 0.1, 4, 1.7), tcl = -0.3, mgp = c(2.75, 0.6, 0), las = 1)
  plot_melting_profile(nicks = profile_nicks_5bp, bg = profile_bg_5bp)
dev.off()
GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/Preferred_Nt_sequence_Melting_profile_5bp.pdf"))

pdf(file = "tmp.pdf", width = 3.5, height = 3)
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(1.5, 0.1, 4, 1.7), tcl = -0.3, mgp = c(2.75, 0.6, 0), las = 1)
  plot_melting_profile_zoom(nicks = profile_nicks_5bp, bg = profile_bg_5bp)
dev.off()
GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/Preferred_Nt_sequence_Melting_profile_5bp_zoom.pdf"))