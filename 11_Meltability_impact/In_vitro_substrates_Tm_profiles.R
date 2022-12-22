# info --------------------------------------------------------------------
# purpose: plot meltability profiles of in vitro substrates
# author: Robert Gnuegge (robert.gnuegge@gmail.com)
# date: 01/13/21
# version: 2.0

# preamble ----------------------------------------------------------------

# set wd to this file's location
wd_path <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(wd_path)

# libraries
library(rmelting)

# load helper files and functions
source(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Analysis_scripts/Misc/Misc_helper_functions.R")
source(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Analysis_scripts/Misc/JFly_colors.R")

plot_dir <- "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Figures/09_Meltability_impact/In_vitro_substrates"

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
  return(out)
}


# in vitro substrate inserts ----------------------------------------------
# derived by random sequence (with positioned C) analysis and manual refinement

end5 <- "GTAAGTGCCGCG"  # substrate ends (to improve annealing)
end3 <- "CACCTCATGCATC"

insert <- "atgACAagagatgag"
medium_melt_sequence <- paste0(c(end5, rep(insert, 3), end3), collapse = "")

insert <- "gggACAagataagat"
medium_melt_sequence_w_profile <- paste0(c(end5, rep(insert, 3), end3), collapse = "")

insert <- "aggACAggtgaggtg"
high_melt_sequence <- paste0(c(end5, rep(insert, 3), end3), collapse = "")

insert <- "tggACAtaggagggg"
high_melt_sequence_w_profile <- paste0(c(end5, rep(insert, 3), end3), collapse = "")

insert <- "tttACAttaggatag"
low_melt_sequence <- paste0(c(end5, rep(insert, 3), end3), collapse = "")

insert <- "tggACAttaatattg"
low_melt_sequence_w_profile <- paste0(c(end5, rep(insert, 3), end3), collapse = "")


# plotting function -------------------------------------------------------

plot_all_melting_profiles <- function(low, medium, high, ylim = NULL){
  if(is.null(ylim)){
    ylim <- range(c(low, medium, high), na.rm = TRUE)
  }
  plot(x = NA, y = NA, xlab = "Position in Substrate [nt]", ylab = NA, xlim = c(1, 70), ylim = ylim)
  title(ylab = expression("Melting Temperature ["*degree*"C"*"]"), line = 2.5)
  
  abline(v = c(17, 32, 47), col = "gray", lty = "dashed")
  
  points(x = 1:70, y = low, type = "l", col = JFly_colors[2])
  points(x = 1:70, y = medium, type = "l", col = JFly_colors[3])
  points(x = 1:70, y = high, type = "l", col = JFly_colors[4])
}

# 5 bp w/ profile
low <- melting_profile(sequence = low_melt_sequence_w_profile, width = 5)
medium <- melting_profile(sequence = medium_melt_sequence_w_profile, width = 5)
high <- melting_profile(sequence = high_melt_sequence_w_profile, width = 5)

pdf(file = "tmp.pdf", width=3.5, height=3)
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(1.7, 0.3, 4, 2), tcl = -0.3, mgp = c(2.25, 0.6, 0), las = 1)
  plot_all_melting_profiles(low = low, medium = medium, high = high, ylim = c(-67, -8))
dev.off()
GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/Substrates_with_profile_5bp.pdf"))

# 5 bp w/o profile
low <- melting_profile(sequence = low_melt_sequence, width = 5)
medium <- melting_profile(sequence = medium_melt_sequence, width = 5)
high <- melting_profile(sequence = high_melt_sequence, width = 5)

pdf(file = "tmp.pdf", width=3.5, height=3)
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(1.7, 0.3, 4, 2), tcl = -0.3, mgp = c(2.25, 0.6, 0), las = 1)
  plot_all_melting_profiles(low = low, medium = medium, high = high, ylim = c(-67, -8))
dev.off()
GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/Substrates_5bp.pdf"))

# 10 bp w/ profile
low <- melting_profile(sequence = low_melt_sequence_w_profile, width = 10)
medium <- melting_profile(sequence = medium_melt_sequence_w_profile, width = 10)
high <- melting_profile(sequence = high_melt_sequence_w_profile, width = 10)

pdf(file = "tmp.pdf", width=3.5, height=3)
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(1.7, 0.3, 4, 2), tcl = -0.3, mgp = c(2.25, 0.6, 0), las = 1)
  plot_all_melting_profiles(low = low, medium = medium, high = high, ylim = c(2, 43))
dev.off()
GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/Substrates_with_profile_10bp.pdf"))

# 10 bp w/o profile
low <- melting_profile(sequence = low_melt_sequence, width = 10)
medium <- melting_profile(sequence = medium_melt_sequence, width = 10)
high <- melting_profile(sequence = high_melt_sequence, width = 10)

pdf(file = "tmp.pdf", width=3.5, height=3)
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(1.7, 0.3, 4, 2), tcl = -0.3, mgp = c(2.25, 0.6, 0), las = 1)
plot_all_melting_profiles(low = low, medium = medium, high = high, ylim = c(2, 43))
dev.off()
GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/Substrates_10bp.pdf"))

# 15 bp w/ profile
low <- melting_profile(sequence = low_melt_sequence_w_profile, width = 15)
medium <- melting_profile(sequence = medium_melt_sequence_w_profile, width = 15)
high <- melting_profile(sequence = high_melt_sequence_w_profile, width = 15)

pdf(file = "tmp.pdf", width=3.5, height=3)
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(1.7, 0.3, 4, 2), tcl = -0.3, mgp = c(2.25, 0.6, 0), las = 1)
  plot_all_melting_profiles(low = low, medium = medium, high = high, ylim = c(29, 56))
dev.off()
GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/Substrates_with_profile_15bp.pdf"))

# 15 bp w/o profile
low <- melting_profile(sequence = low_melt_sequence, width = 15)
medium <- melting_profile(sequence = medium_melt_sequence, width = 15)
high <- melting_profile(sequence = high_melt_sequence, width = 15)

pdf(file = "tmp.pdf", width=3.5, height=3)
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(1.7, 0.3, 4, 2), tcl = -0.3, mgp = c(2.25, 0.6, 0), las = 1)
plot_all_melting_profiles(low = low, medium = medium, high = high, ylim = c(29, 56))
dev.off()
GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/Substrates_15bp.pdf"))


# plot legend -----------------------------------------------------------------
pdf(file = "tmp.pdf", width=2.45, height=1.05)
par(cex = 1, mar = rep(0, 4))

LegTxt <- c(expression("''High T"["m"]*"'' substrate"), 
            expression("''Medium T"["m"]*"'' substrate"), 
            expression("''Low T"["m"]*"'' substrate"), 
            "Positioned C")
plot(1, type="n", axes=FALSE, xlab="", ylab="")
legend(1, 1, xjust=0.5, yjust=0.5, legend = LegTxt,
       lty = c(rep("solid", 3), "dashed"), col = c(JFly_colors[4:2], "gray"))

dev.off()
GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/Legend.pdf"))
