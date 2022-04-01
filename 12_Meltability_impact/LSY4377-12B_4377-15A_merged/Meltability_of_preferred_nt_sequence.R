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

plot_dir <- "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Figures/09_Meltability_impact/LSY4377-12B_4377-15A_merged"

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


# preferred and bg sequences ----------------------------------------------
preferred_seq <- "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTAAAATTTTAAATATACTATATTTAAAAATTTAATTATTTTTTTTATAAATTTATAATTTTTTTTT"
bg_seq <- "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"

profile_nicks_5bp <- melting_profile(sequence = nicks, width = 5)
profile_nicks_10bp <- melting_profile(sequence = nicks, width = 10)
profile_nicks_15bp <- melting_profile(sequence = nicks, width = 15)

profile_bg_5bp <- melting_profile(sequence = bg, width = 5)
profile_bg_10bp <- melting_profile(sequence = bg, width = 10)
profile_bg_15bp <- melting_profile(sequence = bg, width = 15)




# plotting ================================================================

# plotting function -------------------------------------------------------
plot_melting_profile <- function(nicks, bg){
  # plot data
  idx <- which(nicks$dist > -36 & nicks$dist < 36)
  plot(x = 1:70, y = nicks$Tm[idx], type = "l", ylim = range(c(nicks$Tm[idx], bg$Tm[idx])),
       xlab = NA, ylab = expression("Melting Temperature ["*degree*"C"*"]"), xaxt = "n")
  points(x = 1:70, y = bg$Tm[idx], type = "l", lty = "dashed")
  title(xlab = "Distance from Nick Site [nt]", line = 2)
  abline(v = 35.5, col = "gray", lty = "dashed")
  
  # add customized x axis
  x <- nicks$dist[idx]
  at <-  which(x %% 5 == 0)
  labels <- x[at]
  axis(side = 1, at = 1:70, labels = FALSE, tcl = -0.15)
  axis(side = 1, at = c(at, 35, 36), labels = NA, tcl = -0.3)
  axis(side = 1, at = at, labels = labels, tick = FALSE)
  
  # add legend
  legend(x = "bottomleft", legend = c("Expected*", "Observed"), lty = c("dashed", "solid"), bty = "n", inset = 0.03)
}


# 5 bp --------------------------------------------------------------------
nicks <- data.frame(dist = c(-50:-1, 1:50), Tm = melting_profile(sequence = preferred_seq, width = 5))
bg <- data.frame(dist = c(-50:-1, 1:50), Tm = melting_profile(sequence = bg_seq, width = 5))

pdf(file = "tmp.pdf", width = 8, height = 3)
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(2, 0.7, 4, 2), tcl = -0.25, mgp = c(2.25, 0.5, 0), las = 1)
  plot_melting_profile(nicks = nicks, bg = bg)
dev.off()
GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/Preferred_Nt_sequence_Melting_profile_5bp.pdf"))


# 10 bp -------------------------------------------------------------------
nicks <- data.frame(dist = c(-50:-1, 1:50), Tm = melting_profile(sequence = preferred_seq, width = 10))
bg <- data.frame(dist = c(-50:-1, 1:50), Tm = melting_profile(sequence = bg_seq, width = 10))

pdf(file = "tmp.pdf", width = 8, height = 3)
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(2, 0.7, 4, 2), tcl = -0.25, mgp = c(2.25, 0.5, 0), las = 1)
plot_melting_profile(nicks = nicks, bg = bg)
dev.off()
GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/Preferred_Nt_sequence_Melting_profile_10bp.pdf"))


# 15 bp -------------------------------------------------------------------
nicks <- data.frame(dist = c(-50:-1, 1:50), Tm = melting_profile(sequence = preferred_seq, width = 15))
bg <- data.frame(dist = c(-50:-1, 1:50), Tm = melting_profile(sequence = bg_seq, width = 15))

pdf(file = "tmp.pdf", width = 8, height = 3)
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(2, 0.7, 4, 2), tcl = -0.25, mgp = c(2.25, 0.5, 0), las = 1)
plot_melting_profile(nicks = nicks, bg = bg)
dev.off()
GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/Preferred_Nt_sequence_Melting_profile_15bp.pdf"))