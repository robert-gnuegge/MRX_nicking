# info --------------------------------------------------------------------
# purpose: analyze DNA shape impact on MRX nicking
# author: Robert Gnuegge (robert.gnuegge@gmail.com)
# date: 04/25/22
# version: 1.0


# preamble ----------------------------------------------------------------

# set wd to this file's location
wd_path <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(wd_path)

# libraries
library(GenomicRanges)


# load helper files and functions
source(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Analysis_scripts/Misc/Genomic_helper_functions.R")
source(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Analysis_scripts/Misc/S_cerevisiae_SrfI_cut_sites.R")
source(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Analysis_scripts/Misc/Misc_helper_functions.R")
source(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Analysis_scripts/Misc/JFly_colors.R")

plot_dir <- "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Figures/09_Shape_impact/LSY4377-12B_4377-15A_merged"


# function definitions ----------------------------------------------------

# calculate average meltability profiles
# argument: GRanges, logical, character
# result: numeric vector
average_shape <- function(GRanges_DNA_shape, GRanges_S1, bg = FALSE, shape = "MGW"){
  hits <- findOverlaps(query = GRanges_S1, subject = GRanges_DNA_shape)
  S1 <- GRanges_S1[queryHits(hits)]
  DNA_shape <- GRanges_DNA_shape[subjectHits(hits)]
  tmp <- mcols(DNA_shape)[shape]
  if(bg){
    out <- apply(X = as.matrix(tmp), MARGIN = 2, FUN = mean)
  }else{
    out <- apply(X = as.matrix(tmp) * S1$score, MARGIN = 2, FUN = sum) / sum(S1$score)
  }
  return(out)
}


# process data for plotting ===============================================

load(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Data/S1-seq/LSY4377-12B_4377-15A_merged/LSY4377-12B_LSY4377-15A_unnormalized_merged.RData")
# use unnormalized coverage data (for easy scaling with S1-seq score)
load(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Data/S1-seq/LSY4377-12B_4377-15A_merged/LSY4377-12B_LSY4377-15A_merged_resection_extend.RData")
load(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Data/Shape/DNA_shape_around_SrfIcs.RData")

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

# calculate average shapes
all_S1_seq_data <- c(LSY4377_12B_1_merged_S1_seq_unnormalized, LSY4377_12B_2_merged_S1_seq_unnormalized, LSY4377_12B_4_merged_S1_seq_unnormalized)
names(mcols(DNA_shape))
nicks <- data.frame(MGW = average_shape(GRanges_DNA_shape = DNA_shape, GRanges_S1 = all_S1_seq_data, shape = "MGW"))

bg <- data.frame(MGW = average_shape(GRanges_DNA_shape = DNA_shape, GRanges_S1 = all_S1_seq_data, shape = "MGW", bg = TRUE))

# plotting ================================================================

# plotting helper function ------------------------------------------------

# add custom x axis
add_custom_x_axis <- function(){
  at <- c(0:6 * 5 + 1, 8:14 * 5)
  labels <- c(-7:-1, 1:7) * 5
  axis(side = 1, at = 1:70, labels = FALSE, tcl = -0.15)
  axis(side = 1, at = c(at, 35, 36), labels = NA, tcl = -0.3)
  axis(side = 1, at = at, labels = labels, tick = FALSE)
}


# MGW --------------------------------------------------------------------
pdf(file = "tmp.pdf", width = 8, height = 2.5)
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(2, 0.5, 3.7, 2), tcl = -0.25, mgp = c(2.25, 0.5, 0), las = 1)

plot(x = 1:70, y = nicks$MGW, type = "l", ylim = range(c(nicks$MGW, bg$MGW), na.rm = TRUE),
     xlab = NA, ylab = expression("MGW ["*ring(A)*"]"), xaxt = "n")
points(x = 1:70, y = bg$MGW, type = "l", lty = "dashed")
title(xlab = "Distance from Nick Site [nt]", line = 2)
abline(v = 35.5, col = "gray", lty = "solid")
add_custom_x_axis()
legend(x = "topleft", legend = c("Background", "Observed"), lty = c("dashed", "solid"), bty = "n", inset = 0.025)

dev.off()
GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/MGW.pdf"))