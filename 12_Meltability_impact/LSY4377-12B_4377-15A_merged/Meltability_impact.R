# info --------------------------------------------------------------------
# purpose: analyze meltability impact on MRX nicking
# author: Robert Gnuegge (robert.gnuegge@gmail.com)
# date: 01/26/22
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

plot_dir <- "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Figures/09_Meltability_impact/LSY4377-12B_4377-15A_merged"


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

load(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Data/S1-seq/LSY4377-12B_4377-15A_merged/LSY4377-12B_LSY4377-15A_unnormalized_merged.RData")
# use unnormalized coverage data (for easy scaling with S1-seq score)
load(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Data/S1-seq/LSY4377-12B_4377-15A_merged/LSY4377-12B_LSY4377-15A_merged_resection_extend.RData")

load(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Data/Melting/Extracted_melting_profiles.RData")

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
  legend(x = "bottomleft", legend = c("Background", "Observed"), lty = c("dashed", "solid"), bty = "n", inset = 0.03)
}


# 5 bp --------------------------------------------------------------------
nicks <- data.frame(dist = c(-50:-1, 1:50), Tm = average_profile(GRanges_and_profiles = LSY4377_12B_4_5bp_melting_profiles, GRanges = c(LSY4377_12B_1, LSY4377_12B_2, LSY4377_12B_4)))
bg <- data.frame(dist = c(-50:-1, 1:50), Tm = average_profile(GRanges_and_profiles = LSY4377_12B_4_5bp_melting_profiles, GRanges = c(LSY4377_12B_1, LSY4377_12B_2, LSY4377_12B_4), bg = TRUE))

pdf(file = "tmp.pdf", width = 8, height = 3)
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(2, 0.7, 4, 2), tcl = -0.25, mgp = c(2.25, 0.5, 0), las = 1)

plot_melting_profile(nicks = nicks, bg = bg)

dev.off()
GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/Melting_preference_5bp.pdf"))


# 10 bp -------------------------------------------------------------------
nicks <- data.frame(dist = c(-50:-1, 1:50), Tm = average_profile(GRanges_and_profiles = LSY4377_12B_4_10bp_melting_profiles, GRanges = c(LSY4377_12B_1, LSY4377_12B_2, LSY4377_12B_4)))
bg <- data.frame(dist = c(-50:-1, 1:50), Tm = average_profile(GRanges_and_profiles = LSY4377_12B_4_10bp_melting_profiles, GRanges = c(LSY4377_12B_1, LSY4377_12B_2, LSY4377_12B_4), bg = TRUE))

pdf(file = "tmp.pdf", width = 8, height = 3)
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(2, 0.7, 4, 2), tcl = -0.25, mgp = c(2.25, 0.5, 0), las = 1)

plot_melting_profile(nicks = nicks, bg = bg)

dev.off()
GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/Melting_preference_10bp.pdf"))


# 15 bp -------------------------------------------------------------------
nicks <- data.frame(dist = c(-50:-1, 1:50), Tm = average_profile(GRanges_and_profiles = LSY4377_12B_4_15bp_melting_profiles, GRanges = c(LSY4377_12B_1, LSY4377_12B_2, LSY4377_12B_4)))
bg <- data.frame(dist = c(-50:-1, 1:50), Tm = average_profile(GRanges_and_profiles = LSY4377_12B_4_15bp_melting_profiles, GRanges = c(LSY4377_12B_1, LSY4377_12B_2, LSY4377_12B_4), bg = TRUE))

pdf(file = "tmp.pdf", width = 8, height = 3)
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(2, 0.7, 4, 2), tcl = -0.25, mgp = c(2.25, 0.5, 0), las = 1)

plot_melting_profile(nicks = nicks, bg = bg)

dev.off()
GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/Melting_preference_15bp.pdf"))

