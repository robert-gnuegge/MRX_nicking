# info --------------------------------------------------------------------
# purpose: analyze nt sequence impact on MRX nicking
# author: Robert Gnuegge (robert.gnuegge@gmail.com)
# date: 01/26/22
# version: 2.0


# preamble ----------------------------------------------------------------

# set wd to this file's location
wd_path <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(wd_path)

# libraries
library(GenomicRanges)
library(BSgenome)

# read modified S. cerevisiae genome
genome <- import(con = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Data/S_cerevisiae_reference_genomes/S288C_R64-2-1_W303_SNPs_MATa_hocs2SrfIcs_hml_hmr.fasta")

# load helper files and functions
source(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Analysis_scripts/Misc/Genomic_helper_functions.R")
source(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Analysis_scripts/Misc/S_cerevisiae_SrfI_cut_sites.R")
source(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Analysis_scripts/Misc/Misc_helper_functions.R")
source(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Analysis_scripts/Misc/JFly_colors.R")

plot_dir <- "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Figures/08_Sequence_impact/LSY4377-12B_4377-15A_merged"


# function definitions ----------------------------------------------------

# calculate position weight matrices
# argument: GRanges, even integer, genome (DNAStringSet)
# result: matrix
calc_PWM <- function(GRanges, width, genome, calc_beakground_PWM = FALSE){
  if(calc_beakground_PWM){
    nick_pos <- GRanges
  }else{
    nick_pos <- GRanges[rep(1:length(GRanges), GRanges$score)]
    # score counts how often a nick occurred at a specific genome position
    # let's repeat score times to retrieve the sequence context score times (see below)
  }
  if(width %% 2 == 1){
    width <- width + 1  
    warning("width must be an even integer. Changing width to ", width)
  }
  nick_context <- flank(x = nick_pos, width = 0.5 * width, both = TRUE)
  nt_sequences <- getSeq(x = genome, names = nick_context)
  PWM <- consensusMatrix(x = nt_sequences, as.prob = TRUE)[1:4, ]
  colnames(PWM) <- c(-(0.5 * width - 1):-0, 0:(0.5 * width - 1))
  return(PWM)
}

# derive sequence with highest frequency base at each position from PWM
# argument: matrix (PWM)
# result: character string
get_highest_freq_sequence <- function(PWM){
  paste0(rownames(PWM)[apply(X = PWM, MARGIN = 2, FUN = which.max)], collapse = "")
}


# process data for plotting ===============================================

# read S1-seq coverage
load(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Data/S1-seq/LSY4377-12B_4377-15A_merged/LSY4377-12B_LSY4377-15A_unnormalized_merged.RData")
# use unnormalized coverage data (see explanation below)!
load(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Data/S1-seq/LSY4377-12B_4377-15A_merged/LSY4377-12B_LSY4377-15A_merged_resection_extend.RData")

# keep only regions with correct orientation w.r.t DSBs
LSY4377_12B_1 <- subsetByIntersect(subject = LSY4377_12B_1_merged_S1_seq_unnormalized, query = DSB_regions(DSBs = SrfIcs[-c(9, 17)], region_width = 4000, up_rev_down_fw = TRUE))
LSY4377_12B_2 <- subsetByIntersect(subject = LSY4377_12B_2_merged_S1_seq_unnormalized, query = DSB_regions(DSBs = SrfIcs[-c(9, 17)], region_width = 4000, up_rev_down_fw = TRUE))
LSY4377_12B_4 <- subsetByIntersect(subject = LSY4377_12B_4_merged_S1_seq_unnormalized, query = DSB_regions(DSBs = SrfIcs[-c(9, 17)], region_width = 4000, up_rev_down_fw = TRUE))

# retain GRanges in regions around DSBs where resection was detected
# other regions are excluded to prevent pattern "dilution" by background noise
LSY4377_12B_1 <- subsetByIntersect(subject = LSY4377_12B_1, query = resection_extend_LSY4377_12B_1_merged)
LSY4377_12B_2 <- subsetByIntersect(subject = LSY4377_12B_2, query = resection_extend_LSY4377_12B_2_merged)
LSY4377_12B_4 <- subsetByIntersect(subject = LSY4377_12B_4, query = resection_extend_LSY4377_12B_4_merged)

# make nt-resolved GRanges
LSY4377_12B_1 <- as_nt_resolved_GRanges(GRanges = LSY4377_12B_1)
LSY4377_12B_2 <- as_nt_resolved_GRanges(GRanges = LSY4377_12B_2)
LSY4377_12B_4 <- as_nt_resolved_GRanges(GRanges = LSY4377_12B_4)

# remove SrfIcs sequence regions
# they give a lot of S1-seq signal, but derive from unprocessed DSBs and not from MRX nicking
DSBs <- DSB_regions(DSBs = SrfIcs, region_width = 8)
unique(getSeq(x = genome, names = DSBs))  # all SrfIcs sequence

LSY4377_12B_1 <- LSY4377_12B_1[-subjectHits(findOverlaps(query = DSBs, subject = LSY4377_12B_1))]
LSY4377_12B_2 <- LSY4377_12B_2[-subjectHits(findOverlaps(query = DSBs, subject = LSY4377_12B_2))]
LSY4377_12B_4 <- LSY4377_12B_4[-subjectHits(findOverlaps(query = DSBs, subject = LSY4377_12B_4))]

# calculate PWM and background PWM
PWM_1 <- calc_PWM(GRanges = LSY4377_12B_1, width = 100, genome = genome)
PWM_1_bg <- calc_PWM(GRanges = LSY4377_12B_1, width = 100, genome = genome, calc_beakground_PWM = TRUE)

PWM_2 <- calc_PWM(GRanges = LSY4377_12B_2, width = 100, genome = genome)
PWM_2_bg <- calc_PWM(GRanges = LSY4377_12B_2, width = 100, genome = genome, calc_beakground_PWM = TRUE)

PWM_4 <- calc_PWM(GRanges = LSY4377_12B_4, width = 100, genome = genome)
PWM_4_bg <- calc_PWM(GRanges = LSY4377_12B_4, width = 100, genome = genome, calc_beakground_PWM = TRUE)

PWM_all <- calc_PWM(GRanges = c(LSY4377_12B_1, LSY4377_12B_2, LSY4377_12B_4), width = 100, genome = genome)
PWM_all_bg <- calc_PWM(GRanges = c(LSY4377_12B_1, LSY4377_12B_2, LSY4377_12B_4), width = 100, genome = genome, calc_beakground_PWM = TRUE)

# save highest frequency sequences (for melting analysis)
sink(file = "Highest_freq_sequences_in_vivo.txt")
cat("Nicks: ", get_highest_freq_sequence(PWM = PWM_all))
cat("\n\nBackground: ", get_highest_freq_sequence(PWM = PWM_all_bg))
sink()


# plot nt fractions -------------------------------------------------------

# plotting function
plot_nt_fractions <- function(PWM, PWM_bg){
  plot(x = 1:dim(PWM)[2], y = PWM["A", ], type = "l", col = JFly_colors[1], ylim = range(c(PWM, PWM_bg)),
       xlab = "Distance from Nick Site [nt]", ylab = "Fraction", xaxt = "n")
  points(x = 1:dim(PWM_bg)[2], y = PWM_bg["A", ], type = "l", col = JFly_colors[1], lty = "dashed")
  
  axis(side = 1, at = 1:100, labels = FALSE, tcl = -0.2)
  axis(side = 1, at = c(0:5 * 10, 5:9 * 10 + 1), labels = c(NA, -4:-1 * 10, NA, NA, 1:4 * 10), tcl = -0.3)
  axis(side = 1, at = 50.5, labels = 0, tick = FALSE)
  
  points(x = 1:dim(PWM)[2], y = PWM["T", ], type = "l", col = JFly_colors[2])
  points(x = 1:dim(PWM_bg)[2], y = PWM_bg["T", ], type = "l", col = JFly_colors[2], lty = "dashed")
  
  points(x = 1:dim(PWM)[2], y = PWM["G", ], type = "l", col = JFly_colors[3])
  points(x = 1:dim(PWM_bg)[2], y = PWM_bg["G", ], type = "l", col = JFly_colors[3], lty = "dashed")
  
  points(x = 1:dim(PWM)[2], y = PWM["C", ], type = "l", col = JFly_colors[4])
  points(x = 1:dim(PWM_bg)[2], y = PWM_bg["C", ], type = "l", col = JFly_colors[4], lty = "dashed")
  
  abline(v = 50.5, col = "gray")
  
  legend(x = "top", legend = c("A", "T", "G", "C"), col = JFly_colors[1:4], inset = -0.275, lty = "solid", ncol = 4, xpd = TRUE)
}

pdf(file = "tmp.pdf", width = 6, height = 3.5)
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(1.5, 0.5, 1, 1.7), tcl = -0.3, mgp = c(2.5, 0.6, 0), las = 1)
 plot_nt_fractions(PWM = PWM_1, PWM_bg = PWM_1_bg)
dev.off()
GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/Sequence_preference_1h.pdf"))

pdf(file = "tmp.pdf", width = 6, height = 3.5)
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(1.5, 0.5, 1, 1.7), tcl = -0.3, mgp = c(2.5, 0.6, 0), las = 1)
plot_nt_fractions(PWM = PWM_2, PWM_bg = PWM_2_bg)
dev.off()
GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/Sequence_preference_2h.pdf"))

pdf(file = "tmp.pdf", width = 6, height = 3.5)
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(1.5, 0.5, 1, 1.7), tcl = -0.3, mgp = c(2.5, 0.6, 0), las = 1)
plot_nt_fractions(PWM = PWM_4, PWM_bg = PWM_4_bg)
dev.off()
GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/Sequence_preference_4h.pdf"))

pdf(file = "tmp.pdf", width = 6, height = 3.5)
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(1.5, 0.5, 1, 1.7), tcl = -0.3, mgp = c(2.5, 0.6, 0), las = 1)
plot_nt_fractions(PWM = PWM_all, PWM_bg = PWM_all_bg)
dev.off()
GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/Sequence_preference.pdf"))


# compare bg PWM to S. cerevisiae genome-wide nt frequencies --------------

# S. cerevisiae genome-wide N content
tmp <- alphabetFrequency(x = genome)
tmp <- tmp[c(1:16), c("A", "T", "G", "C")]
apply(X = tmp, MARGIN = 2, FUN = sum) / sum(tmp)
#         A         T         G         C 
# 0.3089952 0.3079593 0.1913379 0.1917076

# N content around SrfIcs
apply(X = PWM_all_bg, MARGIN = 1, FUN = mean)
#         A        T         G         C
# 0.2927857 0.3085557 0.1937455 0.2049131


# final figure ============================================================

PWM <- calc_PWM(GRanges = c(LSY4377_12B_1, LSY4377_12B_2, LSY4377_12B_4), width = 70, genome = genome)
PWM_bg <- calc_PWM(GRanges = c(LSY4377_12B_1, LSY4377_12B_2, LSY4377_12B_4), width = 70, genome = genome, calc_beakground_PWM = TRUE)

my_colors <- c(A = JFly_colors[1], T = JFly_colors[2], G = JFly_colors[3], C = JFly_colors[8])

pdf(file = "tmp.pdf", width = 8, height = 3)
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(2, 0.9, 4, 2), tcl = -0.25, mgp = c(2.25, 0.5, 0), las = 1)

# A
plot(x = 1:dim(PWM)[2], y = PWM["A", ], type = "l", col = my_colors["A"], ylim = range(c(PWM, PWM_bg, 0.555)),
     xlab = NA, ylab = "Fraction", xaxt = "n")
title(xlab = "Distance from Nick Site [nt]", line = 2)
points(x = 1:dim(PWM_bg)[2], y = PWM_bg["A", ], type = "l", col = my_colors["A"], lty = "dashed")
abline(v = 35.5, col = "gray", lty = "dashed")

# add customized x axis
x <- c(-35:-1, 1:35)
at <-  which(x %% 5 == 0)
labels <- x[at]
axis(side = 1, at = 1:70, labels = FALSE, tcl = -0.15)
axis(side = 1, at = c(at, 35, 36), labels = NA, tcl = -0.3)
axis(side = 1, at = at, labels = labels, tick = FALSE)


# add other nts
points(x = 1:dim(PWM)[2], y = PWM["T", ], type = "l", col = my_colors["T"])
points(x = 1:dim(PWM_bg)[2], y = PWM_bg["T", ], type = "l", col = my_colors["T"], lty = "dashed")
points(x = 1:dim(PWM)[2], y = PWM["G", ], type = "l", col = my_colors["G"])
points(x = 1:dim(PWM_bg)[2], y = PWM_bg["G", ], type = "l", col = my_colors["G"], lty = "dashed")
points(x = 1:dim(PWM)[2], y = PWM["C", ], type = "l", col = my_colors["C"])
points(x = 1:dim(PWM_bg)[2], y = PWM_bg["C", ], type = "l", col = my_colors["C"], lty = "dashed")

# add legend
text(x = 71, y = c(0.305, 0.305, 0.1925, 0.1925), pos = c(3, 1, 1, 3), offset = 0.175,
     labels = c("T", "A", "G", "C"), col = JFly_colors[c(2, 1, 3, 8)], cex = 1)

# add most frequent nts
most_freq <- rownames(PWM)[apply(X = PWM, MARGIN = 2, FUN = which.max)]
idx <- 1:70
text(x = idx, y = 0.5475, labels = most_freq[idx], col = my_colors[most_freq[idx]], cex = 1)
# underline
y <- 0.5275
x_expand <- 0.4
segments(x0 = c(19, 25, 37, 49) - x_expand, y0 = y, x1 = c(23, 35, 47, 53) + x_expand, y1 = y, col = "gray")
# # rectangles
# rect(xleft = c(9, 58) - x_expand, ybottom = 0.255, xright = c(14, 63) + x_expand, ytop = 0.38, border = "gray", lty = "dotted")

dev.off()
GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/Sequence_preference_final.pdf"))
