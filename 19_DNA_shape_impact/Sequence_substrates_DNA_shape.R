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

plot_dir <- "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Figures/09_Shape_impact/Sequence_substrates"


# our substrates ----------------------------------------------------------
C <- "GTAAGTGCCGCGGAAATACTATATATAAACATTATTTATACTTTATTAAAACAATATCACCTCATGCATC"
G <- "GTAAGTGCCGCGGAAATAGTATATATAAAGATTATTTATAGTTTATTAAAAGAATATCACCTCATGCATC"
noC <- "GTAAGTGCCGCGGAAATAATATATATAAATATTATTTATAATTTATTAAAATAATATCACCTCATGCATC"

# save as fasta files (required for using getShape command) 
sequences <- c(C, G, noC)
names(sequences) <- c("C", "G", "noC")  # add names (will be fasta entry IDs)
sequences <- DNAStringSet(sequences)
file_path <- "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Data/Shape/Sequence_substrate_sequences.fasta"
writeXStringSet(x = sequences, filepath = file_path)

# calculate DNA shape features for each sequence -------------------------
res <- getShape(filename = file_path)
res <- lapply(res, function(x){rownames(x) <- c("C", "G", "noC"); return(x)})  # add row names


# plot MGW ----------------------------------------------------------------
MGW <- res$MGW

pdf(file = "tmp.pdf", width=8, height=2)
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(2, 0.5, 4, 2), tcl = -0.3, mgp = c(2.25, 0.6, 0), las = 1)

plot(x = 1:70, y = MGW["noC", ], type = "l", ylim = range(MGW, na.rm = TRUE),
     ylab = expression("MGW ["*ring(A)*"]"), xlab = NA, xaxt = "n")
axis(side = 1, at = c(1, 1:7 * 10))
title(xlab = "Position in Substrate [nt]", line = 2)
points(x = 1:70, y = MGW["C", ], type = "l", col = JFly_colors[2])
points(x = 1:70, y = MGW["G", ], type = "l", col = JFly_colors[3])
abline(v = c(19, 30, 41, 52), col = "gray", lty = "dashed")
# abline(v = 1, col = "red", xpd = TRUE)

dev.off()
GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/MGW.pdf"))


# rev strand labeled =========================================================

# our substrates ----------------------------------------------------------
noC <- "GATGCATGAGGTGATATTATTTTAATAAATTATAAATAATATTTATATATATTATTTCCGCGGCACTTAC"
G <-   "GATGCATGAGGTGATATTGTTTTAATAAAGTATAAATAATGTTTATATATAGTATTTCCGCGGCACTTAC"
C <-   "GATGCATGAGGTGATATTCTTTTAATAAACTATAAATAATCTTTATATATACTATTTCCGCGGCACTTAC"

# save as fasta files (required for using getShape command) 
sequences <- c(C, G, noC)
names(sequences) <- c("C", "G", "noC")  # add names (will be fasta entry IDs)
sequences <- DNAStringSet(sequences)
file_path <- "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Data/Shape/Sequence_substrate_sequences.fasta"
writeXStringSet(x = sequences, filepath = file_path)

# calculate DNA shape features for each sequence -------------------------
res <- getShape(filename = file_path)
res <- lapply(res, function(x){rownames(x) <- c("C", "G", "noC"); return(x)})  # add row names


# plot MGW ----------------------------------------------------------------
MGW <- res$MGW

pdf(file = "tmp.pdf", width=8, height=2)
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(2, 0.5, 4, 2), tcl = -0.3, mgp = c(2.25, 0.6, 0), las = 1)

plot(x = 1:70, y = MGW["noC", ], type = "l", ylim = range(MGW, na.rm = TRUE),
     ylab = expression("MGW ["*ring(A)*"]"), xlab = NA, xaxt = "n")
axis(side = 1, at = c(1, 1:7 * 10))
title(xlab = "Position in Substrate [nt]", line = 2)
points(x = 1:70, y = MGW["C", ], type = "l", col = JFly_colors[2])
points(x = 1:70, y = MGW["G", ], type = "l", col = JFly_colors[3])
abline(v = c(19, 30, 41, 52), col = "gray", lty = "dashed")
# abline(v = 1, col = "red", xpd = TRUE)

dev.off()
GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/MGW_rev.pdf"))