# info --------------------------------------------------------------------
# purpose: plot average S1-seq spreading over time
# author: Robert Gnuegge (robert.gnuegge@gmail.com)
# date: 01/26/22
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

plot_dir <- "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Figures/07_Average_S1-seq_spreading/LSY4377-12B_4377-15A_merged"


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

# keep only regions with correct orientation w.r.t DSBs
LSY4377_12B_0 <- subsetByIntersect(subject = LSY4377_12B_0_merged_S1_seq, query = DSB_regions(DSBs = SrfIcs[-c(9, 17)], region_width = 5000, up_rev_down_fw = TRUE))
LSY4377_12B_1 <- subsetByIntersect(subject = LSY4377_12B_1_merged_S1_seq, query = DSB_regions(DSBs = SrfIcs[-c(9, 17)], region_width = 5000, up_rev_down_fw = TRUE))
LSY4377_12B_2 <- subsetByIntersect(subject = LSY4377_12B_2_merged_S1_seq, query = DSB_regions(DSBs = SrfIcs[-c(9, 17)], region_width = 5000, up_rev_down_fw = TRUE))
LSY4377_12B_4 <- subsetByIntersect(subject = LSY4377_12B_4_merged_S1_seq, query = DSB_regions(DSBs = SrfIcs[-c(9, 17)], region_width = 5000, up_rev_down_fw = TRUE))
LSY4377_15A_4 <- subsetByIntersect(subject = LSY4377_15A_4_merged_S1_seq, query = DSB_regions(DSBs = SrfIcs[-c(9, 17)], region_width = 5000, up_rev_down_fw = TRUE))

# add slow/middle/fast DSB kinetics mcol
add_DSB_kinetics_category <- function(GRanges){
  mapping <- list(category = c("fast", "middle", "slow"), ranks = list(1:9, 10:16, 17:19))
  GRanges$DSB_kinetics <- NA
  for(n in 1:length(mapping$category)){
    GRanges$DSB_kinetics[GRanges$DSB_kinetics_rank %in% mapping$ranks[[n]]] <- mapping$category[n]
  }
  return(GRanges)
}

LSY4377_12B_0 <- add_DSB_kinetics_category(GRanges = LSY4377_12B_0)
LSY4377_12B_1 <- add_DSB_kinetics_category(GRanges = LSY4377_12B_1)
LSY4377_12B_2 <- add_DSB_kinetics_category(GRanges = LSY4377_12B_2)
LSY4377_12B_4 <- add_DSB_kinetics_category(GRanges = LSY4377_12B_4)
LSY4377_15A_4 <- add_DSB_kinetics_category(GRanges = LSY4377_15A_4)


# calculate average (mean) as function of distance from DSB and grouped by DSB kinetics
agg_1 <- aggregate(score ~ distance_to_DSB + DSB_kinetics, data = LSY4377_12B_1, FUN = mean)
agg_2 <- aggregate(score ~ distance_to_DSB + DSB_kinetics, data = LSY4377_12B_2, FUN = mean)
agg_4 <- aggregate(score ~ distance_to_DSB + DSB_kinetics, data = LSY4377_12B_4, FUN = mean)


# plotting ================================================================
k <- 51
keep <- 2
my_colors <- rev(gray(level = c(0, 0.5, 0.75)))
# my_colors <- JFly_colors[c(1, 4, 5)]


# plot spreading from fast-cut DSBs ---------------------------------------
DSB_kinetics <- "fast"
y_shift <- 7

avg_dist <- data.frame(time = c(0, 1, 2, 4), dist = c(0, 0, 0, 0))

pdf(file = "tmp.pdf", width=3.5, height=3)
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(1.5, 0.5, 3.7, 1.7), tcl = -0.3, mgp = c(2.5, 0.6, 0), las = 1)

  x <- agg_1$distance_to_DSB[agg_1$DSB_kinetics == DSB_kinetics]
  y <- moving_average(x = agg_1$score[agg_1$DSB_kinetics == DSB_kinetics], k = k, keep = keep)
  max(y)
  max(y[-1:-10])
  th <- 55
  th_t <- 295
  f <- 1
  y <- trans(x = y, threshold = th, threshold_trans = th_t, factor = f)
  plot(x = x, y = y, ylim = c(0, trans(x = 310, threshold = th, threshold_trans = th_t, factor = f)), 
       xlim = c(0, 1800), ylab = "Average S1-seq [RPM]", xlab = "Distance from DSB [nt]", 
       type = "l", lwd = 1.5, col = my_colors[1], yaxt = "n")
  
  y <- y[x > keep]
  x <- x[x > keep]
  i <- find_x_where_half_AUC(x = x, y = y)
  yy <- par("usr")[4] - 0.2 * (par("usr")[4] - par("usr")[3]) + y_shift
  lines(x = c(x[i], x[i]), y = c(0, yy), col = my_colors[1], lty = "dashed")
  text(x = x[i], y = yy, labels = x[i], adj = c(0.5, -0.3), col = my_colors[1])
  avg_dist$dist[avg_dist$time == 1] <- x[i]
  
  x <- agg_2$distance_to_DSB[agg_1$DSB_kinetics == DSB_kinetics]
  y <- moving_average(x = agg_2$score[agg_2$DSB_kinetics == DSB_kinetics], k = k, keep = keep)
  max(y)
  max(y[-1])
  y <- trans(x = y, threshold = th, threshold_trans = th_t, factor = f)
  points(x = x, y = y, type = "l", lwd = 1.5, col = my_colors[2])
  
  y <- y[x > keep]
  x <- x[x > keep]
  i <- find_x_where_half_AUC(x = x, y = y)
  yy <- par("usr")[4] - 0.2 * (par("usr")[4] - par("usr")[3])
  lines(x = c(x[i], x[i]), y = c(0, yy), col = my_colors[2], lty = "dashed")
  text(x = x[i], y = yy, labels = x[i], adj = c(0.4, -0.3), col = my_colors[2])
  avg_dist$dist[avg_dist$time == 2] <- x[i]
  
  x <- agg_4$distance_to_DSB[agg_1$DSB_kinetics == DSB_kinetics]
  y <- moving_average(x = agg_4$score[agg_4$DSB_kinetics == DSB_kinetics], k = k, keep = keep)
  max(y)
  max(y[-1])
  y <- trans(x = y, threshold = th, threshold_trans = th_t, factor = f)
  points(x = x, y = y, type = "l", lwd = 1.5, col = my_colors[3])
  
  y <- y[x > keep]
  x <- x[x > keep]
  i <- find_x_where_half_AUC(x = x, y = y)
  yy <- par("usr")[4] - 0.2 * (par("usr")[4] - par("usr")[3]) - y_shift
  lines(x = c(x[i], x[i]), y = c(0, yy), col = my_colors[3], lty = "dashed")
  text(x = x[i], y = yy, labels = x[i], adj = c(0.5, -0.3), col = my_colors[3])
  avg_dist$dist[avg_dist$time == 4] <- x[i]
  
  below_break <- 0:5 * 10
  above_break <- c(300, 310)
  axis(side = 2, at = c(below_break, trans(x = above_break, threshold = th, threshold_trans = th_t, factor = f)), labels = c(below_break, above_break))
  
  axis.break(axis = 2, breakpos = th, style = "slash", bgcol = "white")
  rect(xleft = -10, ybottom = th - 0.5, xright = 10, ytop =  th + 0.5, col = "white", border = "white")
  
  legend(x = "topright", col = my_colors[1:3], lwd = 1.5, legend = paste0(c(1, 2, 4), " h"), inset = c(0.05, 0.05), bty = "n")

dev.off()
GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/Avg_S1-seq_spreading_", DSB_kinetics, ".pdf"))


# plot resection speed over time ------------------------------------------
avg_speed <- data.frame(time = c(1, 2, 4), speed = diff(avg_dist$dist) / diff(avg_dist$time))

pdf(file = "tmp.pdf", width=4, height=3)
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(2, 0.7, 3.7, -1.7), tcl = -0.25, mgp = c(2, 0.5, 0), las = 1)

plot(x = avg_dist$time, y = avg_dist$dist, xlim = c(0, 4), ylim = c(0, max(avg_dist$dist)),
     pch = 20, type = "o", col = "gray", xlab = "Time [h]", ylab = NA, yaxt = "n")
axis_col <- "gray"
axis(side = 4, at = pretty(c(0, max(avg_dist$dist))), col = axis_col, col.axis = axis_col)
text(x = par("usr")[2]*1.25, y = 0.9* mean(par("usr")[3:4]), labels = "Average Distance from DSB [nt]", srt = -90, pos = 3, col = axis_col, xpd = TRUE)

par(new = TRUE)
plot(x = avg_speed$time, y = avg_speed$speed, xlim = c(0, 4), ylim = c(0, max(avg_speed$speed)), 
     pch = 20, type = "o", axes = FALSE, ann = FALSE)
axis(side = 2, at = pretty(c(0, max(avg_speed$speed))))
title(ylab = "Average Resection Speed [nt/h]", line = 2.5)

dev.off()
GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/Avg_resection_speed.pdf"))


# plot spreading from middle-fast cut DSBs --------------------------------
DSB_kinetics <- "middle"
y_shift <- 5

pdf(file = "tmp.pdf", width=3.5, height=3)
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(1.5, 0.5, 3.7, 1.7), tcl = -0.3, mgp = c(2.5, 0.6, 0), las = 1)

  x <- agg_1$distance_to_DSB[agg_1$DSB_kinetics == DSB_kinetics]
  y <- moving_average(x = agg_1$score[agg_1$DSB_kinetics == DSB_kinetics], k = k, keep = keep)
  max(y)
  max(y[-1:-5])
  th <- 35
  th_t <- 50
  f <- 1/20
  y <- trans(x = y, threshold = th, threshold_trans = th_t, factor = f)
  plot(x = x, y = y, ylim = c(0, trans(x = 450, threshold = th, threshold_trans = th_t, factor = f)), 
       xlim = c(0, 1800), ylab = "Average S1-seq [RPM]", xlab = "Distance from DSB [nt]", 
       type = "l", lwd = 1.5, col = my_colors[1], yaxt = "n")
  
  y <- y[x > keep]
  x <- x[x > keep]
  i <- find_x_where_half_AUC(x = x, y = y)
  yy <- par("usr")[4] - 0.2 * (par("usr")[4] - par("usr")[3]) + y_shift
  lines(x = c(x[i], x[i]), y = c(0, yy), col = my_colors[1], lty = "dashed")
  text(x = x[i], y = yy, labels = x[i], adj = c(0.5, -0.3), col = my_colors[1])
  
  x <- agg_2$distance_to_DSB[agg_1$DSB_kinetics == DSB_kinetics]
  y <- moving_average(x = agg_2$score[agg_2$DSB_kinetics == DSB_kinetics], k = k, keep = keep)
  max(y)
  max(y[-1])
  y <- trans(x = y, threshold = th, threshold_trans = th_t, factor = f)
  points(x = x, y = y, type = "l", lwd = 1.5, col = my_colors[2])
  
  y <- y[x > keep]
  x <- x[x > keep]
  i <- find_x_where_half_AUC(x = x, y = y)
  yy <- par("usr")[4] - 0.2 * (par("usr")[4] - par("usr")[3])
  lines(x = c(x[i], x[i]), y = c(0, yy), col = my_colors[2], lty = "dashed")
  text(x = x[i], y = yy, labels = x[i], adj = c(0.5, -0.3), col = my_colors[2])
  
  x <- agg_4$distance_to_DSB[agg_1$DSB_kinetics == DSB_kinetics]
  y <- moving_average(x = agg_4$score[agg_4$DSB_kinetics == DSB_kinetics], k = k, keep = keep)
  max(y)
  max(y[-1])
  y <- trans(x = y, threshold = th, threshold_trans = th_t, factor = f)
  points(x = x, y = y, type = "l", lwd = 1.5, col = my_colors[3])
  
  y <- y[x > keep]
  x <- x[x > keep]
  i <- find_x_where_half_AUC(x = x, y = y)
  yy <- par("usr")[4] - 0.2 * (par("usr")[4] - par("usr")[3]) - y_shift
  lines(x = c(x[i], x[i]), y = c(0, yy), col = my_colors[3], lty = "dashed")
  text(x = x[i], y = yy, labels = x[i], adj = c(0.5, -0.3), col = my_colors[3])
  
  below_break <- 0:3 * 10
  above_break <- c(100, 200, 300, 400)
  axis(side = 2, at = c(below_break, trans(x = above_break, threshold = th, threshold_trans = th_t, factor = f)), labels = c(below_break, above_break))
  
  axis.break(axis = 2, breakpos = th, style = "slash", bgcol = "white")
  rect(xleft = -10, ybottom = th - 0.5, xright = 10, ytop =  th + 0.5, col = "white", border = "white")
  
  legend(x = "topright", col = my_colors[1:3], lwd = 1.5, legend = paste0(c(1, 2, 4), " h"), inset = c(0.05, 0.05), bty = "n")

dev.off()
GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/Avg_S1-seq_spreading_", DSB_kinetics, ".pdf"))


# plot spreading from slowly cut DSBs -------------------------------------
DSB_kinetics <- "slow"
y_shift <- 3.5

pdf(file = "tmp.pdf", width=3.5, height=3)
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(1.5, 0.5, 3.7, 1.7), tcl = -0.3, mgp = c(2.5, 0.6, 0), las = 1)

  x <- agg_1$distance_to_DSB[agg_1$DSB_kinetics == DSB_kinetics]
  xx <- x
  xx[1:5] <-  xx[1:5] - 10  # to prevent covering by other lines
  y <- moving_average(x = agg_1$score[agg_1$DSB_kinetics == DSB_kinetics], k = k, keep = keep)
  max(y)
  max(y[-1:-10])
  th <- 24.5
  th_t <- 25
  f <- 1/10
  y <- trans(x = y, threshold = th, threshold_trans = th_t, factor = f)
  plot(x = xx, y = y, ylim = c(0, trans(x = 160, threshold = th, threshold_trans = th_t, factor = f)), 
       xlim = c(0, 1800), ylab = "Average S1-seq [RPM]", xlab = "Distance from DSB [nt]", 
       type = "l", lwd = 1.5, col = my_colors[1], yaxt = "n")
  
  y <- y[x > keep]
  x <- x[x > keep]
  i <- find_x_where_half_AUC(x = x, y = y)
  yy <- par("usr")[4] - 0.2 * (par("usr")[4] - par("usr")[3]) + y_shift
  lines(x = c(x[i], x[i]), y = c(0, yy), col = my_colors[1], lty = "dashed")
  text(x = x[i], y = yy, labels = x[i], adj = c(0.5, -0.3), col = my_colors[1])
  
  x <- agg_2$distance_to_DSB[agg_1$DSB_kinetics == DSB_kinetics]
  y <- moving_average(x = agg_2$score[agg_2$DSB_kinetics == DSB_kinetics], k = k, keep = keep)
  max(y)
  max(y[-1])
  y <- trans(x = y, threshold = th, threshold_trans = th_t, factor = f)
  points(x = x, y = y, type = "l", lwd = 1.5, col = my_colors[2])
  
  y <- y[x > keep]
  x <- x[x > keep]
  i <- find_x_where_half_AUC(x = x, y = y)
  yy <- par("usr")[4] - 0.2 * (par("usr")[4] - par("usr")[3])
  lines(x = c(x[i], x[i]), y = c(0, yy), col = my_colors[2], lty = "dashed")
  text(x = x[i], y = yy, labels = x[i], adj = c(0.25, -0.3), col = my_colors[2])
  
  x <- agg_4$distance_to_DSB[agg_1$DSB_kinetics == DSB_kinetics]
  y <- moving_average(x = agg_4$score[agg_4$DSB_kinetics == DSB_kinetics], k = k, keep = keep)
  max(y)
  max(y[-1])
  y <- trans(x = y, threshold = th, threshold_trans = th_t, factor = f)
  points(x = x, y = y, type = "l", lwd = 1.5, col = my_colors[3])
  
  y <- y[x > keep]
  x <- x[x > keep]
  i <- find_x_where_half_AUC(x = x, y = y)
  yy <- par("usr")[4] - 0.2 * (par("usr")[4] - par("usr")[3]) - y_shift
  lines(x = c(x[i], x[i]), y = c(0, yy), col = my_colors[3], lty = "dashed")
  text(x = x[i], y = yy, labels = x[i], adj = c(0.5, -0.3), col = my_colors[3])
  
  below_break <- 0:4 * 5
  above_break <- c(50, 100, 150)
  axis(side = 2, at = c(below_break, trans(x = above_break, threshold = th, threshold_trans = th_t, factor = f)), labels = c(below_break, above_break))
  
  axis.break(axis = 2, breakpos = th, style = "slash", bgcol = "white")
  rect(xleft = -15, ybottom = th - 0.25, xright = 10, ytop =  th + 0.25, col = "white", border = "white")
  
  legend(x = "topright", col = my_colors[1:3], lwd = 1.5, legend = paste0(c(1, 2, 4), " h"), inset = c(0.05, 0.05), bty = "n")

dev.off()
GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/Avg_S1-seq_spreading_", DSB_kinetics, ".pdf"))


# plot spreading before DSB induction -------------------------------------
k <- 51
keep <- 2

agg <- aggregate(score ~ distance_to_DSB, data = LSY4377_12B_0, FUN = mean)

pdf(file = "tmp.pdf", width=3.5, height=3)
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(1.5, 0.5, 3.7, 1.7), tcl = -0.3, mgp = c(2.5, 0.6, 0), las = 1)
  x <- agg$distance_to_DSB
  y <- moving_average(x = agg$score, k = k, keep = keep)
  plot(x = x, y = y, xlim = c(0, 1800), ylab = "Average S1-seq [RPM]", xlab = "Distance from DSB [nt]", 
       type = "l", lwd = 1.5, col = gray(level = 0))
dev.off()
GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/Avg_S1-seq_spreading_t0.pdf"))


# plot spreading for mre11-nd ---------------------------------------------
k <- 51
keep <- 2

agg <- aggregate(score ~ distance_to_DSB, data = LSY4377_15A_4, FUN = mean)

pdf(file = "tmp.pdf", width=3.5, height=3)
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(1.5, 0, 3.7, 1.7), tcl = -0.3, mgp = c(2.5, 0.6, 0), las = 1)
  x <- agg$distance_to_DSB
  y <- moving_average(x = agg$score, k = k, keep = keep)
  plot(x = x, y = y, xlim = c(0, 1800), ylab = NA, xlab = "Distance from DSB [nt]", 
       type = "l", lwd = 1.5, col = my_colors[3])
  title(ylab = "Average S1-seq [RPM]", line = 3)
dev.off()
GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/Avg_S1-seq_spreading_mre11-nd.pdf"))

# zoom
pdf(file = "tmp.pdf", width=1.7, height=1.5)
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(3.6, 3, 4, 1.4), tcl = -0.3, mgp = c(2.5, 0.6, 0), las = 1)
  x <- agg$distance_to_DSB
  y <- moving_average(x = agg$score, k = k, keep = keep)
  plot(x = x, y = y, ylim = c(0, 5), xlim = c(0, 1000), xaxt = "n",
       xlab = NA, ylab = NA, type = "l", lwd = 1.5, col = my_colors[3])
  axis(side = 1, at = 0:2 * 500)
dev.off()
GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/Avg_S1-seq_spreading_mre11-nd_zoom.pdf"))