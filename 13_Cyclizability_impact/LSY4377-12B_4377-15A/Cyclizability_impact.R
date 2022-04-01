# info --------------------------------------------------------------------
# purpose: analyze cyclizability impact on MRX nicking
# author: Robert Gnuegge (robert.gnuegge@gmail.com)
# date: 01/18/22
# version: 2.0


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

plot_dir <- "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Figures/10_Cyclizability_impact/LSY4377-12B_4377-15A"


# function definitions ----------------------------------------------------

# calculate average cyclizability profiles
# argument: list with GRanges and matrix (from get_cyclizability_profiles function)
# result: numeric vector
average_profile <- function(GRanges_and_profiles, GRanges, bg = FALSE){
  hits <- findOverlaps(query = GRanges, subject = GRanges_and_profiles$GRanges)
  if(bg){
    out <- apply(X = GRanges_and_profiles$cyclizability_profiles[subjectHits(hits), ], MARGIN = 2, FUN = mean, na.rm = TRUE)
  }else{
    out <- apply(X = GRanges_and_profiles$cyclizability_profiles[subjectHits(hits), ] * GRanges$score[queryHits(hits)], MARGIN = 2, FUN = sum, na.rm = TRUE) / sum(GRanges$score[queryHits(hits)], na.rm = TRUE)
  }
}


# process data for plotting ===============================================

load(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Data/S1-seq/LSY4377-12B_4377-15A/LSY4377-12B_LSY4377-15A_unnormalized.RData")
# use unnormalized coverage data (for easy scaling with S1-seq score)
load(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Data/S1-seq/LSY4377-12B_4377-15A/LSY4377-12B_LSY4377-15A_resection_extend.RData")
load(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Data/Cyclizability/Extracted_cyclizability_profiles.RData")

# keep only regions with correct orientation w.r.t DSBs
LSY4377_12B_1 <- subsetByIntersect(subject = LSY4377_12B_1_S1_seq_unnormalized, query = DSB_regions(DSBs = SrfIcs[4:5], region_width = 4000, up_rev_down_fw = TRUE))
LSY4377_12B_2 <- subsetByIntersect(subject = LSY4377_12B_2_S1_seq_unnormalized, query = DSB_regions(DSBs = SrfIcs[4:5], region_width = 4000, up_rev_down_fw = TRUE))
LSY4377_12B_4 <- subsetByIntersect(subject = LSY4377_12B_4_S1_seq_unnormalized, query = DSB_regions(DSBs = SrfIcs[4:5], region_width = 4000, up_rev_down_fw = TRUE))

# retain GRanges in regions around DSBs where resection was detected
# other regions are excluded to prevent pattern "dilution" by background noise
LSY4377_12B_1 <- subsetByIntersect(subject = LSY4377_12B_1, query = resection_extend_LSY4377_12B_1)
LSY4377_12B_2 <- subsetByIntersect(subject = LSY4377_12B_2, query = resection_extend_LSY4377_12B_2)
LSY4377_12B_4 <- subsetByIntersect(subject = LSY4377_12B_4, query = resection_extend_LSY4377_12B_4)

# make nt-resolved GRanges
LSY4377_12B_1 <- as_nt_resolved_GRanges(GRanges = LSY4377_12B_1)
LSY4377_12B_2 <- as_nt_resolved_GRanges(GRanges = LSY4377_12B_2)
LSY4377_12B_4 <- as_nt_resolved_GRanges(GRanges = LSY4377_12B_4)


# plot cyclizability profiles ---------------------------------------------------

# 200 bp
nicks <- average_profile(GRanges_and_profiles = LSY4377_12B_4_cyclizability_profiles_200bp_interpolated, GRanges = c(LSY4377_12B_1, LSY4377_12B_2, LSY4377_12B_4))
bg <- average_profile(GRanges_and_profiles = LSY4377_12B_4_cyclizability_profiles_200bp_interpolated, GRanges = c(LSY4377_12B_1, LSY4377_12B_2, LSY4377_12B_4), bg = TRUE)
nicks <- moving_average(x = nicks, k = 51, keep = 0)
bg <- moving_average(x = bg, k = 51, keep = 0)

pdf(file = "tmp.pdf", width = 4, height = 3)
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(1.5, -0.5, 4, 1.7), tcl = -0.3, mgp = c(3.5, 0.6, 0), las = 1)
  plot(x = 1:200, y = nicks, type = "l", ylim = range(c(nicks, bg)),
       xlab = NA, ylab = expression("Cyclizability [AU]"), xaxt = "n")
  points(x = 1:200, y = bg, type = "l", lty = "dashed")
  title(xlab = "Distance from Nick Site [nt]", line = 2.5)
  
  axis(side = 1, at = c(1:10 * 10, 10:20 * 10 + 1), labels = FALSE, tcl = -0.2)
  axis(side = 1, at = c(0:2 * 50, 2:4 * 50 + 1), labels = c(-2:-1 * 50, NA, NA, 1:2 * 50), tcl = -0.3)
  axis(side = 1, at = 100.5, labels = 0, tick = FALSE)
  
  abline(v = 100.5, col = "gray")
dev.off()
GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/Cyclizability_profile_200bp.pdf"))

# 500 bp
nicks <- average_profile(GRanges_and_profiles = LSY4377_12B_4_cyclizability_profiles_500bp_interpolated, GRanges = c(LSY4377_12B_1, LSY4377_12B_2, LSY4377_12B_4))
bg <- average_profile(GRanges_and_profiles = LSY4377_12B_4_cyclizability_profiles_500bp_interpolated, GRanges = c(LSY4377_12B_1, LSY4377_12B_2, LSY4377_12B_4), bg = TRUE)
nicks <- moving_average(x = nicks, k = 51, keep = 0)
bg <- moving_average(x = bg, k = 51, keep = 0)

pdf(file = "tmp.pdf", width = 4, height = 3)
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(1.5, -0.5, 3.8, 2), tcl = -0.3, mgp = c(3.5, 0.6, 0), las = 1)
  plot(x = 1:500, y = nicks, type = "l", ylim = range(c(nicks, bg)),
       xlab = NA, ylab = expression("Cyclizability [AU]"), xaxt = "n")
  points(x = 1:500, y = bg, type = "l", lty = "dashed")
  title(xlab = "Distance from Nick Site [nt]", line = 2.5)
  
  axis(side = 1, at = c(0:5 * 50, 5:10 * 50 + 1), labels = FALSE, tcl = -0.2)
  axis(side = 1, at = c(0:2 * 100 + 50, 2:4 * 100 + 51), labels = c(-2:-1 * 100, NA, NA, 1:2 * 100), tcl = -0.3)
  axis(side = 1, at = 250.5, labels = 0, tick = FALSE)
  
  abline(v = 250.5, col = "gray")
dev.off()
GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/Cyclizability_profile_500bp.pdf"))
