# info --------------------------------------------------------------------
# purpose: plot meltability impact on MRX nicking
# author: Robert Gnuegge (robert.gnuegge@gmail.com)
# date: 05/02/22
# version: 3.0


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

plot_dir <- "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Figures/09_Meltability_impact/LSY4377-12B_4377-15A/"


# function definitions ----------------------------------------------------

# calculate average meltability profiles
# argument: GRanges, logical, character
# result: numeric vector
average_Tm_profile <- function(GRanges_Tm, GRanges_S1, bg = FALSE){
  hits <- findOverlaps(query = GRanges_S1, subject = GRanges_Tm)
  S1 <- GRanges_S1[queryHits(hits)]
  Tm <- GRanges_Tm[subjectHits(hits)]
  tmp <- Tm$T_m_profile
  if(bg){
    out <- apply(X = as.matrix(tmp), MARGIN = 2, FUN = mean)
  }else{
    out <- apply(X = as.matrix(tmp) * S1$score, MARGIN = 2, FUN = sum) / sum(S1$score)
  }
  return(out)
}


# process data for plotting ===============================================

load(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Data/S1-seq/LSY4377-12B_4377-15A/LSY4377-12B_LSY4377-15A_unnormalized_merged.RData")
# use unnormalized coverage data (for easy scaling with S1-seq score)
load(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Data/S1-seq/LSY4377-12B_4377-15A/LSY4377-12B_4377-15A_resection_extend.RData")
load(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Data/Melting/Melting_profiles.RData")

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

# calculate average Tm profiles
all_S1_seq_data <- c(LSY4377_12B_1_merged_S1_seq_unnormalized, LSY4377_12B_2_merged_S1_seq_unnormalized, LSY4377_12B_4_merged_S1_seq_unnormalized)
nicks <- data.frame(Tm_5bp = average_Tm_profile(GRanges_Tm = Melting_profiles_5bp, GRanges_S1 = all_S1_seq_data),
                    Tm_10bp = average_Tm_profile(GRanges_Tm = Melting_profiles_10bp, GRanges_S1 = all_S1_seq_data),
                    Tm_15bp = average_Tm_profile(GRanges_Tm = Melting_profiles_15bp, GRanges_S1 = all_S1_seq_data))

bg <- data.frame(Tm_5bp = average_Tm_profile(GRanges_Tm = Melting_profiles_5bp, GRanges_S1 = all_S1_seq_data, bg = TRUE),
                 Tm_10bp = average_Tm_profile(GRanges_Tm = Melting_profiles_10bp, GRanges_S1 = all_S1_seq_data, bg = TRUE),
                 Tm_15bp = average_Tm_profile(GRanges_Tm = Melting_profiles_15bp, GRanges_S1 = all_S1_seq_data, bg = TRUE))

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


# 10 bp -------------------------------------------------------------------
pdf(file = "tmp.pdf", width = 8, height = 2.75)
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(2.3, 0.5, 4, 2), tcl = -0.25, mgp = c(2.5, 0.5, 0), las = 1)

plot(x = 1:70, y = rep(NA, 70), ylim = range(c(nicks$Tm_10bp, bg$Tm_10bp)),
     xlab = NA, ylab = expression("Melting Temperature ["*degree*"C"*"]"), xaxt = "n")
add_custom_x_axis()
title(xlab = "Distance from Nick Site [nt]", line = 1.75)
abline(v = 35.5, col = "gray")

points(x = 1:70, y = nicks$Tm_10bp, type = "l")
points(x = 1:70, y = bg$Tm_10bp, type = "l", lty = "dashed")

legend(x = "bottomleft", legend = c("Background", "Observed"), lty = c("dashed", "solid"), bty = "n", inset = 0.03)

dev.off()
GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/Melting_preference_10bp.pdf"))


# 15 bp -------------------------------------------------------------------
pdf(file = "tmp.pdf", width = 8, height = 2.75)
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(2.3, 0.5, 4, 2), tcl = -0.25, mgp = c(2.5, 0.5, 0), las = 1)

plot(x = 1:70, y = rep(NA, 70), ylim = range(c(nicks$Tm_15bp, bg$Tm_15bp)),
     xlab = NA, ylab = expression("Melting Temperature ["*degree*"C"*"]"), xaxt = "n")
add_custom_x_axis()
title(xlab = "Distance from Nick Site [nt]", line = 1.75)
abline(v = 35.5, col = "gray")

points(x = 1:70, y = nicks$Tm_15bp, type = "l")
points(x = 1:70, y = bg$Tm_15bp, type = "l", lty = "dashed")

legend(x = "bottomleft", legend = c("Background", "Observed"), lty = c("dashed", "solid"), bty = "n", inset = 0.03)

dev.off()
GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/Melting_preference_15bp.pdf"))