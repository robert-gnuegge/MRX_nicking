# info --------------------------------------------------------------------
# purpose: analyze (a)symmetry of S1-seq up- and downstream of DSBs
# author: Robert Gnuegge (robert.gnuegge@gmail.com)
# date: 01/27/22
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

plot_dir <- "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Figures/06_Asymmetry/LSY4602-20C_1"


# process data for plotting ===============================================

# read S1-seq coverage
load(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Data/S1-seq/LSY4377-12B_4377-15A_merged/LSY4377-12B_LSY4377-15A_merged_around_DSBs.RData")
load(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Data/S1-seq/LSY4377-12B_4377-15A_merged/LSY4377-12B_LSY4377-15A_merged_resection_extend.RData")
load(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Data/S1-seq/LSY4602-20C_1/LSY4602_20C_1_around_DSBs.RData")
load(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Data/S1-seq/LSY4602-20C_1/LSY4602_20C_1_resection_extend.RData")

# retain GRanges in roi
LSY4377_12B_1 <- subsetByIntersect(subject = LSY4377_12B_1_merged_S1_seq, query = resection_extend_LSY4377_12B_1_merged)
LSY4602_20C_1 <- subsetByIntersect(subject = LSY4602_20C_1_S1_seq, query = resection_extend_LSY4602_20C_1)


# calculate asymmetry -----------------------------------------------------
# as sum(S1-seq scores upstream DSB) / sum(S1-seq scores downstream DSB)

calc_up_down_ratio <- function(x){ sum(x[strand(x) == "-"]$score) / sum(x[strand(x) == "+"]$score) }

tmp <- sapply(X = split(LSY4377_12B_1, ~ DSB_id), FUN = calc_up_down_ratio)
up_down_ratio <- data.frame(DSB_id = names(tmp), sample = "LSY4377_12B_1", score = tmp)

tmp <- sapply(X = split(LSY4602_20C_1, ~ DSB_id), FUN = calc_up_down_ratio)
up_down_ratio <- rbind(up_down_ratio, data.frame(DSB_id = names(tmp), sample = "LSY4602_20C_1", score = tmp))


# plot asymmetry ----------------------------------------------------------

mean_score <- sapply(X = split(up_down_ratio, ~ DSB_id), FUN = function(x){ mean(x$score) })
idx <- order(mean_score, decreasing = TRUE)
# idx <- order(up_down_ratio$score[up_down_ratio$sample == "LSY4377_12B_1"], decreasing = TRUE)

DSBs <- unique(LSY4377_12B_1$DSB_id)

pdf(file = "tmp.pdf", width = 7, height = 4)
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(0.3, -1.5, 0.7, 2), tcl = -0.3, mgp = c(3.5, 0.6, 0), las = 1)
  shift <- 0.05
  plot(x = 1:length(DSBs) - 1 * shift, y = up_down_ratio$score[up_down_ratio$sample == "LSY4377_12B_1"][idx], ylim = range(up_down_ratio$score),
       pch = 21, col = NA, bg = adjustcolor(col = JFly_colors[1], alpha.f = 0.5), 
       xlab = NA, xaxt = "n", ylab = "Ratio of Summed S1-seq\nUp-/Downstream of DSB")
  
  text(x = 1:length(DSBs), y = par("usr")[3] - 0.025 * (par("usr")[4] - par("usr")[3]), srt = 45, 
       labels = up_down_ratio$DSB_id[up_down_ratio$sample == "LSY4377_12B_1"][idx], xpd = TRUE, adj = c(1,1))
  
  points(x = 1:length(DSBs) + 0 * shift, y = up_down_ratio$score[up_down_ratio$sample == "LSY4602_20C_1"][idx],
         pch = 21, col = NA, bg = adjustcolor(col = JFly_colors[2], alpha.f = 0.5))
  
  abline(h = 1, lty = "dashed", col = "gray")
  
  legend(x = "top", inset = -0.27, ncol = 3, xpd = TRUE,
         pch = 21, col = NA, pt.bg = adjustcolor(col = JFly_colors[1:2], alpha.f = 0.5), 
         legend = c("WT", expression(italic("ku70"*Delta))))
  
dev.off()
GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/Asymmetry.pdf"))