# info --------------------------------------------------------------------
# purpose: plot average S1-seq spreading over time
# author: Robert Gnuegge (robert.gnuegge@gmail.com)
# date: 01/27/22
# version: 2.0


# preamble ----------------------------------------------------------------

# set wd to this file's location
wd_path <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(wd_path)

# libraries
library(GenomicRanges)
library(plotrix)

# load helper files and functions
source(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Analysis_scripts/Misc/Genomic_helper_functions.R")
source(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Analysis_scripts/Misc/S_cerevisiae_SrfI_cut_sites.R")
source(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Analysis_scripts/Misc/Misc_helper_functions.R")
source(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Analysis_scripts/Misc/JFly_colors.R")

plot_dir <- "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Figures/07_Average_S1-seq_spreading/LSY4602-20C_1"


# function definitions ----------------------------------------------------

# find x for which AUC is half of total AUC
# argument: numeric vectors
# result: double
find_x_where_half_AUC <- function(x, y){
  AUC <- sum(diff(x) * (head(y,-1)+tail(y,-1)))/2  # area under curve (source: https://stackoverflow.com/a/30280873/11705274)
  # find x where half AUC is reached  
  AUC.x <- 0
  z <- 1
  while (AUC.x <= 0.5 * AUC) {
    z <- z + 1
    AUC.x <- sum(diff(x[1:z]) * (head(y[1:z],-1)+tail(y[1:z],-1)))/2
  }
  return(x[z])
}

# transform data for axis break with scale change
# argument: numeric vector
# result: numeric vector
trans <- function(x, threshold, threshold_trans, factor){
  pmin(x, threshold) + factor * pmax(x - threshold_trans, 0)
}

# process data for plotting ===============================================

# read S1-seq coverage
load(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Data/S1-seq/LSY4377-12B_4377-15A_merged/LSY4377-12B_LSY4377-15A_merged_around_DSBs.RData")
load(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Data/S1-seq/LSY4602-20C_1/LSY4602_20C_1_around_DSBs.RData")

# keep only regions with correct orientation w.r.t DSBs
LSY4377_12B_1 <- subsetByIntersect(subject = LSY4377_12B_1_merged_S1_seq, query = DSB_regions(DSBs = SrfIcs[-c(9, 17)], region_width = 5000, up_rev_down_fw = TRUE))
LSY4602_20C_1 <- subsetByIntersect(subject = LSY4602_20C_1_S1_seq, query = DSB_regions(DSBs = SrfIcs[-c(9, 17)], region_width = 5000, up_rev_down_fw = TRUE))

# add slow/middle/fast DSB kinetics mcol
add_DSB_kinetics_category <- function(GRanges){
  mapping <- list(category = c("fast", "middle", "slow"), ranks = list(1:9, 10:16, 17:19))
  GRanges$DSB_kinetics <- NA
  for(n in 1:length(mapping$category)){
    GRanges$DSB_kinetics[GRanges$DSB_kinetics_rank %in% mapping$ranks[[n]]] <- mapping$category[n]
  }
  return(GRanges)
}

LSY4377_12B_1 <- add_DSB_kinetics_category(GRanges = LSY4377_12B_1)
LSY4602_20C_1 <- add_DSB_kinetics_category(GRanges = LSY4602_20C_1)


# calculate average (mean) as function of distance from DSB and grouped by DSB kinetics
agg_WT <- aggregate(score ~ distance_to_DSB + DSB_kinetics, data = LSY4377_12B_1, FUN = mean)
agg_ku70 <- aggregate(score ~ distance_to_DSB + DSB_kinetics, data = LSY4602_20C_1, FUN = mean)

# plotting ================================================================
k <- 51
keep <- 2
my_colors <- gray(level = c(0, 0.75))


# plot spreading from fast-cut DSBs ---------------------------------------
DSB_kinetics <- "fast"
y_shift <- 6

pdf(file = "tmp.pdf", width=3.5, height=3)
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(1.5, 0.5, 3.7, 1.7), tcl = -0.3, mgp = c(2.5, 0.6, 0), las = 1)

  x <- agg_WT$distance_to_DSB[agg_WT$DSB_kinetics == DSB_kinetics]
  y <- moving_average(x = agg_WT$score[agg_WT$DSB_kinetics == DSB_kinetics], k = k, keep = keep)
  max(y)
  max(y[-1:-10])
  th <- 65
  th_t <- 295
  f <- 1
  y <- trans(x = y, threshold = th, threshold_trans = th_t, factor = f)
  plot(x = x, y = y, ylim = c(0, trans(x = 310, threshold = th, threshold_trans = th_t, factor = f)), 
       xlim = c(0, 1800), ylab = "Average S1-seq [RPM]", xlab = "Distance from DSB [nt]", 
       type = "l", lwd = 1.5, col = my_colors[1], yaxt = "n")
  
  # y <- y[x > keep]
  # x <- x[x > keep]
  # i <- find_x_where_half_AUC(x = x, y = y)
  # yy <- par("usr")[4] - 0.2 * (par("usr")[4] - par("usr")[3]) + y_shift
  # lines(x = c(x[i], x[i]), y = c(0, yy), col = my_colors[1], lty = "dashed")
  # text(x = x[i], y = yy, labels = x[i], adj = c(0.5, -0.3), col = my_colors[1])
  
  x <- agg_ku70$distance_to_DSB[agg_WT$DSB_kinetics == DSB_kinetics]
  y <- moving_average(x = agg_ku70$score[agg_ku70$DSB_kinetics == DSB_kinetics], k = k, keep = keep)
  max(y)
  max(y[-1])
  y <- trans(x = y, threshold = th, threshold_trans = th_t, factor = f)
  points(x = x, y = y, type = "l", lwd = 1.5, col = my_colors[2])
  
  # y <- y[x > keep]
  # x <- x[x > keep]
  # i <- find_x_where_half_AUC(x = x, y = y)
  # yy <- par("usr")[4] - 0.2 * (par("usr")[4] - par("usr")[3])
  # lines(x = c(x[i], x[i]), y = c(0, yy), col = my_colors[2], lty = "dashed")
  # text(x = x[i], y = yy, labels = x[i], adj = c(1, -0.3), col = my_colors[2])
  
  below_break <- 0:6 * 10
  above_break <- c(300, 310)
  axis(side = 2, at = c(below_break, trans(x = above_break, threshold = th, threshold_trans = th_t, factor = f)), labels = c(below_break, above_break))
  
  axis.break(axis = 2, breakpos = th, style = "slash", bgcol = "white")
  rect(xleft = -10, ybottom = th - 0.5, xright = 10, ytop =  th + 0.5, col = "white", border = "white")
  
  legend(x = "topright", col = my_colors[1:2], lwd = 1.5, legend = c("WT", expression(italic("yku70"*Delta))), inset = c(0.05, 0.05))

dev.off()
GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/Avg_S1-seq_spreading_", DSB_kinetics, ".pdf"))

# zoom --------------------------------------------------------------------
DSB_kinetics <- "fast"

pdf(file = "tmp.pdf", width=3.5, height=3)
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(1.5, 0.5, 4, 2), tcl = -0.3, mgp = c(2.5, 0.6, 0), las = 1)

  x <- agg_WT$distance_to_DSB[agg_WT$DSB_kinetics == DSB_kinetics]
  y <- moving_average(x = agg_WT$score[agg_WT$DSB_kinetics == DSB_kinetics], k = k, keep = keep)
  y[1]
  max(y[10:70])
  th <- 25
  th_t <- 26
  f <- 1 / 40
  y <- trans(x = y, threshold = th, threshold_trans = th_t, factor = f)
  plot(x = x, y = y, ylim = c(0, trans(x = 308, threshold = th, threshold_trans = th_t, factor = f)), 
       xlim = c(0, 70), ylab = "Average S1-seq [RPM]", xlab = "Distance from DSB [nt]", 
       type = "l", lwd = 1.5, col = my_colors[1], yaxt = "n")
  
  x <- agg_ku70$distance_to_DSB[agg_WT$DSB_kinetics == DSB_kinetics]
  y <- moving_average(x = agg_ku70$score[agg_ku70$DSB_kinetics == DSB_kinetics], k = k, keep = keep)
  y[1]
  max(y[10:70])
  y <- trans(x = y, threshold = th, threshold_trans = th_t, factor = f)
  points(x = x, y = y, type = "l", lwd = 1.5, col = my_colors[2])
  
  below_break <- 0:4 * 5
  above_break <- 1:3 * 100
  axis(side = 2, at = c(below_break, trans(x = above_break, threshold = th, threshold_trans = th_t, factor = f)), labels = c(below_break, above_break))
  
  axis.break(axis = 2, breakpos = th, style = "slash", bgcol = "white")
  rect(xleft = -0.5, ybottom = th - 0.25, xright = 5, ytop =  th + 0.25, col = "white", border = "white")
  
  legend(x = "top", col = my_colors[1:2], lwd = 1.5, legend = c("WT", expression(italic("yku70"*Delta))), inset = 0.03)

dev.off()
GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/Avg_S1-seq_spreading_", DSB_kinetics, "_zoom.pdf"))


# plot spreading from middle-fast cut DSBs --------------------------------
DSB_kinetics <- "middle"
y_shift <- 5

pdf(file = "tmp.pdf", width=3.5, height=3)
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(1.5, 0.5, 3.7, 1.7), tcl = -0.3, mgp = c(2.5, 0.6, 0), las = 1)

  x <- agg_WT$distance_to_DSB[agg_WT$DSB_kinetics == DSB_kinetics]
  y <- moving_average(x = agg_WT$score[agg_WT$DSB_kinetics == DSB_kinetics], k = k, keep = keep)
  max(y)
  max(y[-1:-5])
  th <- 35
  th_t <- 150
  f <- 1/20
  y <- trans(x = y, threshold = th, threshold_trans = th_t, factor = f)
  plot(x = x, y = y, ylim = c(0, trans(x = 435, threshold = th, threshold_trans = th_t, factor = f)), 
       xlim = c(0, 1800), ylab = "Average S1-seq [RPM]", xlab = "Distance from DSB [nt]", 
       type = "l", lwd = 1.5, col = my_colors[1], yaxt = "n")
  
  # y <- y[x > keep]
  # x <- x[x > keep]
  # i <- find_x_where_half_AUC(x = x, y = y)
  # yy <- par("usr")[4] - 0.2 * (par("usr")[4] - par("usr")[3]) + y_shift
  # lines(x = c(x[i], x[i]), y = c(0, yy), col = my_colors[1], lty = "dashed")
  # text(x = x[i], y = yy, labels = x[i], adj = c(0.5, -0.3), col = my_colors[1])
  
  x <- agg_ku70$distance_to_DSB[agg_WT$DSB_kinetics == DSB_kinetics]
  y <- moving_average(x = agg_ku70$score[agg_ku70$DSB_kinetics == DSB_kinetics], k = k, keep = keep)
  max(y)
  max(y[-1:-5])
  y <- trans(x = y, threshold = th, threshold_trans = th_t, factor = f)
  points(x = x, y = y, type = "l", lwd = 1.5, col = my_colors[2])
  
  # y <- y[x > keep]
  # x <- x[x > keep]
  # i <- find_x_where_half_AUC(x = x, y = y)
  # yy <- par("usr")[4] - 0.2 * (par("usr")[4] - par("usr")[3])
  # lines(x = c(x[i], x[i]), y = c(0, yy), col = my_colors[2], lty = "dashed")
  # text(x = x[i], y = yy, labels = x[i], adj = c(0.5, -0.3), col = my_colors[2])
  
  below_break <- 0:3 * 10
  above_break <- c(200, 300, 400)
  axis(side = 2, at = c(below_break, trans(x = above_break, threshold = th, threshold_trans = th_t, factor = f)), labels = c(below_break, above_break))
  
  axis.break(axis = 2, breakpos = th, style = "slash", bgcol = "white")
  rect(xleft = -10, ybottom = th - 0.5, xright = 10, ytop =  th + 0.5, col = "white", border = "white")
  
  legend(x = "topright", col = my_colors[1:2], lwd = 1.5, legend = c("WT", expression(italic("yku70"*Delta))), inset = c(0.05, 0.05))

dev.off()
GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/Avg_S1-seq_spreading_", DSB_kinetics, ".pdf"))

# zoom --------------------------------------------------------------------
DSB_kinetics <- "middle"

pdf(file = "tmp.pdf", width=3.5, height=3)
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(1.5, 0.5, 4, 2), tcl = -0.3, mgp = c(2.5, 0.6, 0), las = 1)

  x <- agg_WT$distance_to_DSB[agg_WT$DSB_kinetics == DSB_kinetics]
  y <- moving_average(x = agg_WT$score[agg_WT$DSB_kinetics == DSB_kinetics], k = k, keep = keep)
  y[1]
  max(y[10:70])
  th <- 13
  th_t <- 140
  f <- 1/80
  y <- trans(x = y, threshold = th, threshold_trans = th_t, factor = f)
  plot(x = x, y = y, ylim = c(0, trans(x = 435, threshold = th, threshold_trans = th_t, factor = f)), 
       xlim = c(0, 70), ylab = "Average S1-seq [RPM]", xlab = "Distance from DSB [nt]", 
       type = "l", lwd = 1.5, col = my_colors[1], yaxt = "n")
  
  x <- agg_ku70$distance_to_DSB[agg_WT$DSB_kinetics == DSB_kinetics]
  y <- moving_average(x = agg_ku70$score[agg_ku70$DSB_kinetics == DSB_kinetics], k = k, keep = keep)
  y[1]
  max(y[10:70])
  y <- trans(x = y, threshold = th, threshold_trans = th_t, factor = f)
  points(x = x, y = y, type = "l", lwd = 1.5, col = my_colors[2])
  
  below_break <- 0:6 * 2
  above_break <- c(200, 300, 400)
  axis(side = 2, at = c(below_break, trans(x = above_break, threshold = th, threshold_trans = th_t, factor = f)), labels = c(below_break, above_break))
  
  axis.break(axis = 2, breakpos = th, style = "slash", bgcol = "white")
  rect(xleft = -1, ybottom = th - 0.5, xright = 5, ytop =  th + 0.5, col = "white", border = "white")
  
  legend(x = "top", col = my_colors[1:2], lwd = 1.5, legend = c("WT", expression(italic("yku70"*Delta))), inset = 0.03)

dev.off()
GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/Avg_S1-seq_spreading_", DSB_kinetics, "_zoom.pdf"))


# plot spreading from slowly cut DSBs -------------------------------------
DSB_kinetics <- "slow"
y_shift <- 3.5

pdf(file = "tmp.pdf", width=3.5, height=3)
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(1.5, 0.5, 3.7, 1.7), tcl = -0.3, mgp = c(2.5, 0.6, 0), las = 1)

  x <- agg_WT$distance_to_DSB[agg_WT$DSB_kinetics == DSB_kinetics]
  xx <- x
  xx[1:5] <-  xx[1:5] - 10  # to prevent covering by other lines
  y <- moving_average(x = agg_WT$score[agg_WT$DSB_kinetics == DSB_kinetics], k = k, keep = keep)
  max(y)
  max(y[-1:-5])
  th <- 4.75
  th_t <- 70
  f <- 1/40
  y <- trans(x = y, threshold = th, threshold_trans = th_t, factor = f)
  plot(x = xx, y = y, ylim = c(0, trans(x = 125, threshold = th, threshold_trans = th_t, factor = f)), 
       xlim = c(0, 1800), ylab = "Average S1-seq [RPM]", xlab = "Distance from DSB [nt]", 
       type = "l", lwd = 1.5, col = my_colors[1], yaxt = "n")
  
  # y <- y[x > keep]
  # x <- x[x > keep]
  # i <- find_x_where_half_AUC(x = x, y = y)
  # yy <- par("usr")[4] - 0.2 * (par("usr")[4] - par("usr")[3]) + y_shift
  # lines(x = c(x[i], x[i]), y = c(0, yy), col = my_colors[1], lty = "dashed")
  # text(x = x[i], y = yy, labels = x[i], adj = c(0.5, -0.3), col = my_colors[1])
  
  x <- agg_ku70$distance_to_DSB[agg_WT$DSB_kinetics == DSB_kinetics]
  y <- moving_average(x = agg_ku70$score[agg_ku70$DSB_kinetics == DSB_kinetics], k = k, keep = keep)
  max(y)
  max(y[-1:-5])
  y <- trans(x = y, threshold = th, threshold_trans = th_t, factor = f)
  points(x = x, y = y, type = "l", lwd = 1.5, col = my_colors[2])
  
  # y <- y[x > keep]
  # x <- x[x > keep]
  # i <- find_x_where_half_AUC(x = x, y = y)
  # yy <- par("usr")[4] - 0.2 * (par("usr")[4] - par("usr")[3])
  # lines(x = c(x[i], x[i]), y = c(0, yy), col = my_colors[2], lty = "dashed")
  # text(x = x[i], y = yy, labels = x[i], adj = c(0.25, -0.3), col = my_colors[2])
  
  below_break <- 0:4
  above_break <- c(80, 100, 120)
  axis(side = 2, at = c(below_break, trans(x = above_break, threshold = th, threshold_trans = th_t, factor = f)), labels = c(below_break, above_break))
  
  axis.break(axis = 2, breakpos = th, style = "slash", bgcol = "white")
  rect(xleft = -15, ybottom = th - 0.075, xright = 10, ytop =  th + 0.075, col = "white", border = "white")
  
  legend(x = "topright", col = my_colors[1:2], lwd = 1.5, legend = c("WT", expression(italic("yku70"*Delta))), inset = c(0.05, 0.05))

dev.off()
GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/Avg_S1-seq_spreading_", DSB_kinetics, ".pdf"))

# zoom --------------------------------------------------------------------
DSB_kinetics <- "slow"

pdf(file = "tmp.pdf", width=3.5, height=3)
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(1.5, 0.5, 4, 2), tcl = -0.3, mgp = c(2.5, 0.6, 0), las = 1)

  x <- agg_WT$distance_to_DSB[agg_WT$DSB_kinetics == DSB_kinetics]
  xx <- x
  xx[1:5] <-  xx[1:5] - 0.25  # to prevent covering by other lines
  y <- moving_average(x = agg_WT$score[agg_WT$DSB_kinetics == DSB_kinetics], k = k, keep = keep)
  y[1]
  max(y[10:70])
  th <- 2.75
  th_t <- 65
  f <- 1/50
  y <- trans(x = y, threshold = th, threshold_trans = th_t, factor = f)
  plot(x = xx, y = y, ylim = c(0, trans(x = 122, threshold = th, threshold_trans = th_t, factor = f)), 
       xlim = c(0, 70), ylab = "Average S1-seq [RPM]", xlab = "Distance from DSB [nt]", 
       type = "l", lwd = 1.5, col = my_colors[1], yaxt = "n")
  
  x <- agg_ku70$distance_to_DSB[agg_WT$DSB_kinetics == DSB_kinetics]
  y <- moving_average(x = agg_ku70$score[agg_ku70$DSB_kinetics == DSB_kinetics], k = k, keep = keep)
  y[1]
  max(y[10:70])
  y <- trans(x = y, threshold = th, threshold_trans = th_t, factor = f)
  points(x = x, y = y, type = "l", lwd = 1.5, col = my_colors[2])
  
  below_break <- 0:5 * 0.5
  above_break <- c(80, 100, 120)
  axis(side = 2, at = c(below_break, trans(x = above_break, threshold = th, threshold_trans = th_t, factor = f)), labels = c(below_break, above_break))
  
  axis.break(axis = 2, breakpos = th, style = "slash", bgcol = "white")
  rect(xleft = -0.5, ybottom = th - 0.05, xright = 5, ytop =  th + 0.05, col = "white", border = "white")
  
  legend(x = "top", col = my_colors[1:2], lwd = 1.5, legend = c("WT", expression(italic("yku70"*Delta))), inset = 0.03)

dev.off()
GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/Avg_S1-seq_spreading_", DSB_kinetics, "_zoom.pdf"))


# plot spreading from all DSBs zoom =======================================
keep <- 2
k <- 51
y_shift <- 3.5

agg_WT <- aggregate(score ~ distance_to_DSB, data = LSY4377_12B_1, FUN = mean)
agg_ku70 <- aggregate(score ~ distance_to_DSB, data = LSY4602_20C_1, FUN = mean)

pdf(file = "tmp.pdf", width=3.5, height=3)
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(1.5, 0.5, 4, 2), tcl = -0.3, mgp = c(2.5, 0.6, 0), las = 1)

  x <- agg_WT$distance_to_DSB
  y <- moving_average(x = agg_WT$score, k = k, keep = keep)
  y[1]
  max(y[10:70])
  th <- 17
  th_t <- 30
  f <- 1/65
  y <- trans(x = y, threshold = th, threshold_trans = th_t, factor = f)
  plot(x = x, y = y, ylim = c(0, trans(x = 320, threshold = th, threshold_trans = th_t, factor = f)), 
       xlim = c(0, 70), ylab = "Average S1-seq [RPM]", xlab = "Distance from DSB [nt]", 
       type = "l", lwd = 1.5, col = my_colors[1], yaxt = "n")
  
  x <- agg_ku70$distance_to_DSB
  y <- moving_average(x = agg_ku70$score, k = k, keep = keep)
  y[1]
  max(y[10:70])
  y <- trans(x = y, threshold = th, threshold_trans = th_t, factor = f)
  points(x = x, y = y, type = "l", lwd = 1.5, col = my_colors[2])
  
  below_break <- 0:3 * 5
  above_break <- c(100, 200, 300)
  axis(side = 2, at = c(below_break, trans(x = above_break, threshold = th, threshold_trans = th_t, factor = f)), labels = c(below_break, above_break))
  
  axis.break(axis = 2, breakpos = th, style = "slash", bgcol = "white")
  rect(xleft = -1, ybottom = th - 0.125, xright = 5, ytop =  th + 0.125, col = "white", border = "white")
  
  legend(x = "top", col = my_colors[1:2], lwd = 1.5, legend = c("WT", expression(italic("yku70"*Delta))), inset = 0.03, bty = "n")

dev.off()
GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/Avg_S1-seq_spreading_all_zoom.pdf"))

# plot spreading from all DSBs ============================================
pdf(file = "tmp.pdf", width=3.5, height=3)
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(1.5, 0.5, 4, 1.6), tcl = -0.3, mgp = c(2.5, 0.6, 0), las = 1)

x <- agg_WT$distance_to_DSB
y <- moving_average(x = agg_WT$score, k = k, keep = keep)
y[1]
max(y[-1:-10])
th <- 42.5
th_t <- 50
f <- 1/25
y <- trans(x = y, threshold = th, threshold_trans = th_t, factor = f)
plot(x = x, y = y, ylim = c(0, trans(x = 320, threshold = th, threshold_trans = th_t, factor = f)), 
     ylab = "Average S1-seq [RPM]", xlab = "Distance from DSB [nt]", xlim = c(0, 1000),
     type = "l", lwd = 1.5, col = my_colors[1], yaxt = "n")

y <- y[x > keep]
x <- x[x > keep]
i <- find_x_where_half_AUC(x = x, y = y)
yy <- par("usr")[4] - 0.2 * (par("usr")[4] - par("usr")[3]) - 0.5 * y_shift
lines(x = c(x[i], x[i]), y = c(0, yy), col = my_colors[1], lty = "dashed")
text(x = x[i], y = yy, labels = x[i], adj = c(0, -0.1), col = my_colors[1])

x <- agg_ku70$distance_to_DSB
y <- moving_average(x = agg_ku70$score, k = k, keep = keep)
y[1]
max(y[-1:-10])
y <- trans(x = y, threshold = th, threshold_trans = th_t, factor = f)
points(x = x, y = y, type = "l", lwd = 1.5, col = my_colors[2])

y <- y[x > keep]
x <- x[x > keep]
i <- find_x_where_half_AUC(x = x, y = y)
yy <- par("usr")[4] - 0.2 * (par("usr")[4] - par("usr")[3]) + 0.5 * y_shift
lines(x = c(x[i], x[i]), y = c(0, yy), col = my_colors[2], lty = "dashed")
text(x = x[i], y = yy, labels = x[i], adj = c(0.8, -0.3), col = my_colors[2])

below_break <- 0:8 * 5
above_break <- c(100, 200, 300)
axis(side = 2, at = c(below_break, trans(x = above_break, threshold = th, threshold_trans = th_t, factor = f)), labels = c(below_break, above_break))

axis.break(axis = 2, breakpos = th, style = "slash", bgcol = "white")
rect(xleft = -1, ybottom = th - 0.125, xright = 5, ytop =  th + 0.125, col = "white", border = "white")

legend(x = "topright", col = my_colors[1:2], lwd = 1.5, legend = c("WT", expression(italic("yku70"*Delta))), inset = 0.03)

dev.off()
GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/Avg_S1-seq_spreading_all.pdf"))
