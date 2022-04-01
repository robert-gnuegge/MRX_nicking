# info --------------------------------------------------------------------
# purpose: analyze meltability impact on MRX nicking
# author: Robert Gnuegge (robert.gnuegge@gmail.com)
# date: 01/12/22
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

plot_dir <- "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Figures/09_Meltability_impact/LSY4377-12B_4377-15A"


# function definitions ----------------------------------------------------

# calculate average meltability profiles
# argument: list with GRanges and matrix (from get_meltability_profiles function)
# result: numeric vector
average_profile <- function(GRanges_and_profiles, GRanges, bg = FALSE){
  hits <- findOverlaps(query = GRanges, subject = GRanges_and_profiles$GRanges)
  if(bg){
    out <- apply(X = GRanges_and_profiles$meltability_profiles[subjectHits(hits), ], MARGIN = 2, FUN = mean)
  }else{
    out <- apply(X = GRanges_and_profiles$meltability_profiles[subjectHits(hits), ] * GRanges$score[queryHits(hits)], MARGIN = 2, FUN = sum) / sum(GRanges$score[queryHits(hits)])
  }
}


# process data for plotting ===============================================

load(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Data/S1-seq/LSY4377-12B_4377-15A/LSY4377-12B_LSY4377-15A_unnormalized.RData")
# use unnormalized coverage data (for easy scaling with S1-seq score)
load(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Data/S1-seq/LSY4377-12B_4377-15A/LSY4377-12B_LSY4377-15A_resection_extend.RData")

load(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Data/Melting/Extracted_melting_profiles.RData")

# keep only regions with correct orientation w.r.t DSBs
LSY4377_12B_1 <- subsetByIntersect(subject = LSY4377_12B_1_S1_seq_unnormalized, query = DSB_regions(DSBs = SrfIcs[-c(9, 17)], region_width = 4000, up_rev_down_fw = TRUE))
LSY4377_12B_2 <- subsetByIntersect(subject = LSY4377_12B_2_S1_seq_unnormalized, query = DSB_regions(DSBs = SrfIcs[-c(9, 17)], region_width = 4000, up_rev_down_fw = TRUE))
LSY4377_12B_4 <- subsetByIntersect(subject = LSY4377_12B_4_S1_seq_unnormalized, query = DSB_regions(DSBs = SrfIcs[-c(9, 17)], region_width = 4000, up_rev_down_fw = TRUE))

# retain GRanges in regions around DSBs where resection was detected
# other regions are excluded to prevent pattern "dilution" by background noise
LSY4377_12B_1 <- subsetByIntersect(subject = LSY4377_12B_1, query = resection_extend_LSY4377_12B_1)
LSY4377_12B_2 <- subsetByIntersect(subject = LSY4377_12B_2, query = resection_extend_LSY4377_12B_2)
LSY4377_12B_4 <- subsetByIntersect(subject = LSY4377_12B_4, query = resection_extend_LSY4377_12B_4)

# make nt-resolved GRanges
LSY4377_12B_1 <- as_nt_resolved_GRanges(GRanges = LSY4377_12B_1)
LSY4377_12B_2 <- as_nt_resolved_GRanges(GRanges = LSY4377_12B_2)
LSY4377_12B_4 <- as_nt_resolved_GRanges(GRanges = LSY4377_12B_4)


# plot melting profiles ---------------------------------------------------

# plotting function
plot_melting_profile <- function(nicks, bg){
  plot(x = 1:100, y = nicks, type = "l", ylim = range(c(nicks, bg)),
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
  plot(x = 1:41, y = nicks[31:71], type = "l", ylim = range(c(nicks[31:71], bg[31:71])),
       xlab = NA, ylab = expression("Melting Temperature ["*degree*"C"*"]"), xaxt = "n")
  points(x = 1:41, y = bg[31:71], type = "l", lty = "dashed")
  title(xlab = "Distance from Nick Site [nt]", line = 2.5)
  
  axis(side = 1, at = 1:41, labels = FALSE, tcl = -0.2)
  axis(side = 1, at = c(0:2 * 10, 2:4 * 10 + 1), labels = c(-2:-1 * 10, NA, NA, 1:2 * 10), tcl = -0.3)
  axis(side = 1, at = 20.5, labels = 0, tick = FALSE)
  
  abline(v = 20.5, col = "gray")
}

# 15 bp
nicks <- average_profile(GRanges_and_profiles = LSY4377_12B_4_15bp_melting_profiles, GRanges = c(LSY4377_12B_1, LSY4377_12B_2, LSY4377_12B_4))
bg <- average_profile(GRanges_and_profiles = LSY4377_12B_4_15bp_melting_profiles, GRanges = c(LSY4377_12B_1, LSY4377_12B_2, LSY4377_12B_4), bg = TRUE)
avg_profiles_15bp <- data.frame(dist = c(-49:0, 0:49), nicks = nicks, bg = bg)

pdf(file = "tmp.pdf", width = 5.5, height = 3)
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(1.5, 0.1, 4, 1.7), tcl = -0.3, mgp = c(2.75, 0.6, 0), las = 1)
  plot_melting_profile(nicks = nicks, bg = bg)
dev.off()
GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/Melting_profile_15bp.pdf"))

pdf(file = "tmp.pdf", width = 3.5, height = 3)
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(1.5, 0.1, 4, 1.7), tcl = -0.3, mgp = c(2.75, 0.6, 0), las = 1)
  plot_melting_profile_zoom(nicks = nicks, bg = bg)
dev.off()
GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/Melting_profile_15bp_zoom.pdf"))

# 10 bp
nicks <- average_profile(GRanges_and_profiles = LSY4377_12B_4_10bp_melting_profiles, GRanges = c(LSY4377_12B_1, LSY4377_12B_2, LSY4377_12B_4))
bg <- average_profile(GRanges_and_profiles = LSY4377_12B_4_10bp_melting_profiles, GRanges = c(LSY4377_12B_1, LSY4377_12B_2, LSY4377_12B_4), bg = TRUE)
avg_profiles_10bp <- data.frame(dist = c(-49:0, 0:49), nicks = nicks, bg = bg)

pdf(file = "tmp.pdf", width = 5.5, height = 3)
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(1.5, 0.1, 4, 1.7), tcl = -0.3, mgp = c(2.75, 0.6, 0), las = 1)
  plot_melting_profile(nicks = nicks, bg = bg)
dev.off()
GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/Melting_profile_10bp.pdf"))

pdf(file = "tmp.pdf", width = 3.5, height = 3)
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(1.5, 0.1, 4, 1.7), tcl = -0.3, mgp = c(2.75, 0.6, 0), las = 1)
  plot_melting_profile_zoom(nicks = nicks, bg = bg)
dev.off()
GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/Melting_profile_10bp_zoom.pdf"))

# 5 bp
nicks <- average_profile(GRanges_and_profiles = LSY4377_12B_4_5bp_melting_profiles, GRanges = c(LSY4377_12B_1, LSY4377_12B_2, LSY4377_12B_4))
bg <- average_profile(GRanges_and_profiles = LSY4377_12B_4_5bp_melting_profiles, GRanges = c(LSY4377_12B_1, LSY4377_12B_2, LSY4377_12B_4), bg = TRUE)
avg_profiles_5bp <- data.frame(dist = c(-49:0, 0:49), nicks = nicks, bg = bg)

pdf(file = "tmp.pdf", width = 5.5, height = 3)
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(1.5, 0.1, 4, 1.7), tcl = -0.3, mgp = c(2.75, 0.6, 0), las = 1)
  plot_melting_profile(nicks = nicks, bg = bg)
dev.off()
GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/Melting_profile_5bp.pdf"))

pdf(file = "tmp.pdf", width = 3.5, height = 3)
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(1.5, 0.1, 4, 1.7), tcl = -0.3, mgp = c(2.75, 0.6, 0), las = 1)
  plot_melting_profile_zoom(nicks = nicks, bg = bg)
dev.off()
GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/Melting_profile_5bp_zoom.pdf"))