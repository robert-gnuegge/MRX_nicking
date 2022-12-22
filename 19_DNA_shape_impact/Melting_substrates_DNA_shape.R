# info --------------------------------------------------------------------
# purpose: plot DNA shape features of in vitro substrates
# author: Robert Gnuegge (robert.gnuegge@gmail.com)
# date: 04/26/22
# version: 2.0

# preamble ----------------------------------------------------------------

# set wd to this file's location
wd_path <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(wd_path)

# libraries
library(DNAshapeR)
library(Biostrings)

# load helper files and functions
source(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Analysis_scripts/Misc/Misc_helper_functions.R")
source(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Analysis_scripts/Misc/JFly_colors.R")

plot_dir <- "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Figures/09_Shape_impact/Melting_substrates"


# our substrates ----------------------------------------------------------
low_flat <- "GTAAGTGCCGCGTTTACATTAGGATAGTTTACATTAGGATAGTTTACATTAGGATAGCACCTCATGCATC"
low_profile <- "GTAAGTGCCGCGTGGACATTAATATTGTGGACATTAATATTGTGGACATTAATATTGCACCTCATGCATC"
mid_flat <- "GTAAGTGCCGCGATGACAAGAGATGAGATGACAAGAGATGAGATGACAAGAGATGAGCACCTCATGCATC"
mid_profile <- "GTAAGTGCCGCGGGGACAAGATAAGATGGGACAAGATAAGATGGGACAAGATAAGATCACCTCATGCATC"
high_flat <- "GTAAGTGCCGCGAGGACAGGTGAGGTGAGGACAGGTGAGGTGAGGACAGGTGAGGTGCACCTCATGCATC"
high_profile <- "GTAAGTGCCGCGTGGACATAGGAGGGGTGGACATAGGAGGGGTGGACATAGGAGGGGCACCTCATGCATC"

# save as fasta files (required for using getShape command) 
sequences <- c(low_flat, low_profile, mid_flat, mid_profile, high_flat, high_profile)
names(sequences) <- c("low_flat", "low_profile", "mid_flat", "mid_profile", "high_flat", "high_profile")  # add names (will be fasta entry IDs)
sequences <- DNAStringSet(sequences)
file_path <- "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Data/Shape/Sequence_substrate_sequences.fasta"
writeXStringSet(x = sequences, filepath = file_path)

# calculate DNA shape features for each sequence -------------------------
res <- getShape(filename = file_path)
res <- lapply(res, function(x){rownames(x) <- names(sequences); return(x)})  # add row names

MGW <- res$MGW

# plot MGW low Tm -------------------------------------------------------
pdf(file = "tmp.pdf", width=8, height=2.5)
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(2, 0.5, 4, 2), tcl = -0.3, mgp = c(2.25, 0.6, 0), las = 1)

plot(x = 1:70, y = MGW["low_flat", ], type = "l", col = JFly_colors[1], ylim = range(MGW, na.rm = TRUE),
     ylab = expression("Minor groove width ["*ring(A)*"]"), xlab = NA, xaxt = "n")
axis(side = 1, at = c(1, 1:7 * 10))
title(xlab = "Position in Substrate [nt]", line = 2)
points(x = 1:70, y = MGW["low_profile", ], type = "l", col = JFly_colors[2])

abline(v = c(17, 32, 47), col = "gray", lty = "dashed")

dev.off()
GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/MGW_low.pdf"))


# plot legend -------------------------------------------------------------
Legend_txt <- c("- Profile", "+ Profile")

pdf(file = "tmp.pdf", width=1.35, height=0.65)
par(cex = 1, mar = rep(0, 4))
plot(1, type="n", axes=FALSE, xlab="", ylab="")
legend(1, 1, xjust=0.5, yjust=0.5, legend = Legend_txt, col = JFly_colors[1:2], lty = "solid")
dev.off()
GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/Legend.pdf"))


# plot MGW mid Tm ---------------------------------------------------------
pdf(file = "tmp.pdf", width=8, height=2.5)
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(2, 0.5, 4, 2), tcl = -0.3, mgp = c(2.25, 0.6, 0), las = 1)

plot(x = 1:70, y = MGW["mid_flat", ], type = "l", col = JFly_colors[1], ylim = range(MGW, na.rm = TRUE),
     ylab = expression("Minor groove width ["*ring(A)*"]"), xlab = NA, xaxt = "n")
axis(side = 1, at = c(1, 1:7 * 10))
title(xlab = "Position in Substrate [nt]", line = 2)
points(x = 1:70, y = MGW["mid_profile", ], type = "l", col = JFly_colors[2])

abline(v = c(17, 32, 47), col = "gray", lty = "dashed")

dev.off()
GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/MGW_mid.pdf"))


# plot MGW high Tm ---------------------------------------------------------
pdf(file = "tmp.pdf", width=8, height=2.5)
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(2, 0.5, 4, 2), tcl = -0.3, mgp = c(2.25, 0.6, 0), las = 1)

plot(x = 1:70, y = MGW["high_flat", ], type = "l", col = JFly_colors[1], ylim = range(MGW, na.rm = TRUE),
     ylab = expression("Minor groove width ["*ring(A)*"]"), xlab = NA, xaxt = "n")
axis(side = 1, at = c(1, 1:7 * 10))
title(xlab = "Position in Substrate [nt]", line = 2)
points(x = 1:70, y = MGW["high_profile", ], type = "l", col = JFly_colors[2])

abline(v = c(17, 32, 47), col = "gray", lty = "dashed")

dev.off()
GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/MGW_high.pdf"))
