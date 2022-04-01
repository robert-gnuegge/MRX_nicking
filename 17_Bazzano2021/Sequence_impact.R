# info --------------------------------------------------------------------
# purpose: analyze nt sequence impact on MRX nicking in Bazzano et al., 2021 data
# author: Robert Gnuegge (robert.gnuegge@gmail.com)
# date: 02/13/22
# version: 1.0


# preamble ----------------------------------------------------------------

# set wd to this file's location
wd_path <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(wd_path)

# libraries
library(GenomicRanges)
library(BSgenome)

# read modified S. cerevisiae genome
genome <- import(con = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Resources/Bazzano2021/Reference/ILV1-L.fasta")

# load helper files and functions
source(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Analysis_scripts/Misc/Genomic_helper_functions.R")
source(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Analysis_scripts/Misc/Misc_helper_functions.R")
source(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Analysis_scripts/Misc/JFly_colors.R")

plot_dir <- "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Figures/16_Bazzano2021/Sequence_impact"


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


# process data for plotting ===============================================

# read S1-seq coverage
load(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Data/Bazzano2021/Coverage/exo1_sgs1_coverage_unnormalized.RData")
# use unnormalized coverage data (see explanation below)!

# make nt-resolved GRanges
exo1_sgs1_T35 <- as_nt_resolved_GRanges(GRanges = exo1_sgs1_T35_unnormalized)
exo1_sgs1_T60 <- as_nt_resolved_GRanges(GRanges = exo1_sgs1_T60_unnormalized)
exo1_sgs1_T75 <- as_nt_resolved_GRanges(GRanges = exo1_sgs1_T75_unnormalized)
exo1_sgs1_T90 <- as_nt_resolved_GRanges(GRanges = exo1_sgs1_T90_unnormalized)

# exclude HOcs site and HOcs-distal region, as only HOcs-proximal region will be nicked
# HOcs is at 530, let's include region until 520, to avoid unprocessed DSB coverage
# let's exclude first 101 nts, as they were added to reference sequence (couldn't find out why)
roi <- GRanges(seqnames = "dsb", ranges = IRanges(start = 101, end = 520))
T35 <- subsetByOverlaps(x = exo1_sgs1_T35, ranges = roi)
T60 <- subsetByOverlaps(x = exo1_sgs1_T60, ranges = roi)
T75 <- subsetByOverlaps(x = exo1_sgs1_T75, ranges = roi)
T90 <- subsetByOverlaps(x = exo1_sgs1_T90, ranges = roi)

# calculate PWM and background PWM
PWM <- calc_PWM(GRanges = c(T35, T60, T75, T90), width = 70, genome = genome)
PWM_bg <- calc_PWM(GRanges = c(T35, T60, T75, T90), width = 70, genome = genome, calc_beakground_PWM = TRUE)


# plotting ----------------------------------------------------------------

my_colors <- c(A = JFly_colors[1], T = JFly_colors[2], G = JFly_colors[3], C = JFly_colors[8])

pdf(file = "tmp.pdf", width = 8, height = 3)
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(2, 0.9, 4, 2), tcl = -0.25, mgp = c(2.25, 0.5, 0), las = 1)

# A
plot(x = 1:dim(PWM)[2], y = PWM["A", ], type = "l", col = my_colors["A"], ylim = range(c(PWM, PWM_bg)),
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
text(x = 70, y = PWM_bg[c("C", "G"), 70], adj = c(-0.33, 0.5),
     labels = names(PWM_bg[c("C", "G"), 70]), col = my_colors[c("C", "G")], cex = 1)
text(x = 70, y = PWM_bg["A", 70], adj = c(-0.33, 0),
     labels = "A", col = my_colors["A"], cex = 1)
text(x = 70, y = PWM_bg["T", 70], adj = c(-0.33, 1),
     labels = "T", col = my_colors["T"], cex = 1)

# # add most frequent nts
# most_freq <- rownames(PWM)[apply(X = PWM, MARGIN = 2, FUN = which.max)]
# idx <- 1:70
# text(x = idx, y = 0.5575, labels = most_freq[idx], col = my_colors[most_freq[idx]], cex = 1)

# # underline nts that are equal to S1-seq-derived motif
# my_most_freq <- strsplit(x = "TTTTTTTTTTTTTTTTTTTTAAAATTTTAAATATACTATATTTAAAAATTTAATTATTTTTTTTATAAAT", split = "")[[1]]
# which(most_freq == my_most_freq)
# y <- 0.5375
# x_expand <- 0.4
# segments(x0 = c(2, 5, 14, 18, 21, 23, 27, 33, 35, 38, 41, 44, 48, 50, 59, 62, 65, 68) - x_expand, y0 = y, 
#          x1 = c(2, 5, 16, 18, 21, 24, 31, 33, 36, 39, 42, 46, 48, 53, 60, 62, 65, 69) + x_expand, y1 = y, col = "gray")

dev.off()
GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/Sequence_preference.pdf"))
