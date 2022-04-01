# info --------------------------------------------------------------------
# purpose: evaluate correlation between biological replicates
# author: Robert Gnuegge (robert.gnuegge@gmail.com)
# date: 01/25/21
# version: 2.0

# preamble ----------------------------------------------------------------

# set wd to this file's location
wd_path <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(wd_path)

# libraries
library(GenomicAlignments)

# load helper files and functions
source(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Analysis_scripts/Misc/Genomic_helper_functions.R")
source(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Analysis_scripts/Misc/Misc_helper_functions.R")
source(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Analysis_scripts/Misc/S_cerevisiae_SrfI_cut_sites.R")

plot_dir <- "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Figures/02_Replicate_correlation/MNase-seq"



# read and process all files ----------------------------------------------

load(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Data/MNase-seq/LSY4377-12B_4377-15A/LSY4377-12B_LSY4377-15A_trimmed.RData")
load(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Data/MNase-seq/LSY4377-12B_4377-15A_rep2/LSY4377-12B_LSY4377-15A_rep2_trimmed.RData")

roi <- DSB_regions(DSBs = SrfIcs, region_width = 4000)

LSY4377_12B_0_rep1 <- sort(x = as_nt_resolved_GRanges(GRanges = subsetByIntersect(subject = LSY4377_12B_0_MNase_seq_trimmed, query = roi)), ignore.strand = FALSE)
LSY4377_12B_1_rep1 <- sort(x = as_nt_resolved_GRanges(GRanges = subsetByIntersect(subject = LSY4377_12B_1_MNase_seq_trimmed, query = roi)), ignore.strand = FALSE)
LSY4377_12B_2_rep1 <- sort(x = as_nt_resolved_GRanges(GRanges = subsetByIntersect(subject = LSY4377_12B_2_MNase_seq_trimmed, query = roi)), ignore.strand = FALSE)
LSY4377_12B_4_rep1 <- sort(x = as_nt_resolved_GRanges(GRanges = subsetByIntersect(subject = LSY4377_12B_4_MNase_seq_trimmed, query = roi)), ignore.strand = FALSE)
LSY4377_12B_0_rep2 <- sort(x = as_nt_resolved_GRanges(GRanges = subsetByIntersect(subject = LSY4377_12B_0_rep2_MNase_seq_trimmed, query = roi)), ignore.strand = FALSE)
LSY4377_12B_1_rep2 <- sort(x = as_nt_resolved_GRanges(GRanges = subsetByIntersect(subject = LSY4377_12B_1_rep2_MNase_seq_trimmed, query = roi)), ignore.strand = FALSE)
LSY4377_12B_2_rep2 <- sort(x = as_nt_resolved_GRanges(GRanges = subsetByIntersect(subject = LSY4377_12B_2_rep2_MNase_seq_trimmed, query = roi)), ignore.strand = FALSE)
LSY4377_12B_4_rep2 <- sort(x = as_nt_resolved_GRanges(GRanges = subsetByIntersect(subject = LSY4377_12B_4_rep2_MNase_seq_trimmed, query = roi)), ignore.strand = FALSE)

# make sure all GRanges have same order in both replicate files
all(granges(LSY4377_12B_0_rep1) == granges(LSY4377_12B_0_rep2))
all(granges(LSY4377_12B_1_rep1) == granges(LSY4377_12B_1_rep2))
all(granges(LSY4377_12B_2_rep1) == granges(LSY4377_12B_2_rep2))
all(granges(LSY4377_12B_4_rep1) == granges(LSY4377_12B_4_rep2))


# plot correlation --------------------------------------------------------

# plotting function
plot_correlation <- function(rep1, rep2){
  xy_range <- range(c(rep1, rep2))  # to use same range on both axes
  plot(x = NA, y = NA, xlim = xy_range, ylim = xy_range, xaxt = "n", yaxt = "n",
       xlab = expression("log"[10]*"(Replicate 1 [RPM] + 1)"), ylab = expression("log"[10]*"(Replicate 2 [RPM] + 1)"))
  axis(side = 1, at = pretty(x = xy_range, n = 3)) 
  axis(side = 2, at = pretty(x = xy_range, n = 3))
  symbols(x = rep1, y = rep2, circles = rep(1, length(rep1)), inches = 0.01, add = TRUE, fg = NA, bg = gray(level = 0, alpha = 0.25))  # plot data as semitransparent dots
  abline(lm(rep2 ~ rep1), col = "red")  # add linear regression line
  # add Pearson's correlation coefficient
  r <- round(cor(x = rep1, y = rep2, method = "pearson"), 2)
  text(x = par("usr")[1] + 0.1 * (par("usr")[2] - par("usr")[1]), y = par("usr")[4] - 0.1 * (par("usr")[4]-par("usr")[3]), adj = c(0, 1),
       labels = substitute(italic("r")~"="~r))
}


# t = 0 -------------------------------------------------------------------

file_name <- paste0(plot_dir, "/LSY4377_12B_0.png")
png(filename = file_name, width = 3, height = 3, units = "in", res = 900)
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(1.7, 0.5, 3.8, 1.8), tcl = -0.3, mgp = c(2.5, 0.6, 0), las = 1)
  plot_correlation(rep1 = log10(x = LSY4377_12B_0_rep1$score + 1), rep2 = log10(x = LSY4377_12B_0_rep2$score + 1))
dev.off()

pdf(file = "tmp.pdf", width = 3, height = 3)
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(1.7, 0.5, 3.8, 1.8), tcl = -0.3, mgp = c(2.5, 0.6, 0), las = 1)
  plot_correlation(rep1 = log10(x = LSY4377_12B_0_rep1$score + 1), rep2 = log10(x = LSY4377_12B_0_rep2$score + 1))
dev.off()
GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/LSY4377_12B_0.pdf"))


# t = 1 -------------------------------------------------------------------

file_name <- paste0(plot_dir, "/LSY4377_12B_1.png")
png(filename = file_name, width = 3, height = 3, units = "in", res = 900)
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(1.7, 0.5, 3.8, 1.8), tcl = -0.3, mgp = c(2.5, 0.6, 0), las = 1)
plot_correlation(rep1 = log10(x = LSY4377_12B_1_rep1$score + 1), rep2 = log10(x = LSY4377_12B_1_rep2$score + 1))
dev.off()

pdf(file = "tmp.pdf", width = 3, height = 3)
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(1.7, 0.5, 3.8, 1.8), tcl = -0.3, mgp = c(2.5, 0.6, 0), las = 1)
plot_correlation(rep1 = log10(x = LSY4377_12B_1_rep1$score + 1), rep2 = log10(x = LSY4377_12B_1_rep2$score + 1))
dev.off()
GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/LSY4377_12B_1.pdf"))


# t = 2 -------------------------------------------------------------------

file_name <- paste0(plot_dir, "/LSY4377_12B_2.png")
png(filename = file_name, width = 3, height = 3, units = "in", res = 900)
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(1.7, 0.5, 3.8, 1.8), tcl = -0.3, mgp = c(2.5, 0.6, 0), las = 1)
plot_correlation(rep1 = log10(x = LSY4377_12B_2_rep1$score + 1), rep2 = log10(x = LSY4377_12B_2_rep2$score + 1))
dev.off()

pdf(file = "tmp.pdf", width = 3, height = 3)
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(1.7, 0.5, 3.8, 1.8), tcl = -0.3, mgp = c(2.5, 0.6, 0), las = 1)
plot_correlation(rep1 = log10(x = LSY4377_12B_2_rep1$score + 1), rep2 = log10(x = LSY4377_12B_2_rep2$score + 1))
dev.off()
GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/LSY4377_12B_2.pdf"))


# t = 4 -------------------------------------------------------------------

file_name <- paste0(plot_dir, "/LSY4377_12B_4.png")
png(filename = file_name, width = 3, height = 3, units = "in", res = 900)
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(1.7, 0.5, 3.8, 1.8), tcl = -0.3, mgp = c(2.5, 0.6, 0), las = 1)
plot_correlation(rep1 = log10(x = LSY4377_12B_4_rep1$score + 1), rep2 = log10(x = LSY4377_12B_4_rep2$score + 1))
dev.off()

pdf(file = "tmp.pdf", width = 3, height = 3)
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(1.7, 0.5, 3.8, 1.8), tcl = -0.3, mgp = c(2.5, 0.6, 0), las = 1)
plot_correlation(rep1 = log10(x = LSY4377_12B_4_rep1$score + 1), rep2 = log10(x = LSY4377_12B_4_rep2$score + 1))
dev.off()
GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/LSY4377_12B_4.pdf"))