# info --------------------------------------------------------------------
# purpose: group S1-seq by DSB formation kinetics for subsequent meta analyses
# author: Robert Gnuegge (robert.gnuegge@gmail.com)
# date: 01/27/21
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

plot_dir <- "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Figures/05_DSB_formation/LSY4602-20C_1"
save_dir <- "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Data/S1-seq/LSY4602-20C_1"

# function definitions ====================================================

# add distance to nearest DSB
# argument: GRanges object
# result: GRanges object
add_distance_to_nearest_DSB <- function(subject, DSBs, add_DSB_id = TRUE){
  hits <- distanceToNearest(x = subject, subject = DSBs, ignore.strand = FALSE)
  out <- subject
  out$distance_to_DSB <- mcols(hits)$distance
  if(add_DSB_id){
    out$DSB_id <- as.character(DSBs)[subjectHits(hits)]
  }
  # we need to correct the distance for the + strand entries (they are shifted by -1, except for the ones at the DSBs)
  out[strand(out) == "+"]$distance_to_DSB <- out[strand(out) == "+"]$distance_to_DSB + 1
  strand(DSBs) <- "+"
  hits <- findOverlaps(query = DSBs, subject = out)
  out[subjectHits(hits)]$distance_to_DSB <- 0
  return(out)
}


# find resection extend
# argument: GRanges
# result: GRanges
find_resection_extend <- function(GRanges, DSBs, max_DSB_dist = 2500, threshold = 0.95){
  
  out <- DSBs
  out$DSB_id <- as.character(DSBs)
  
  for(n in 1:length(out)){
    # find tract end on - strand
    roi <- DSBs[n]
    strand(roi) <- "-"
    roi <- resize(x = roi, width = max_DSB_dist, fix = "start")
    tmp <- subsetByIntersect(subject =GRanges, query = roi)
    tmp <- rev(tmp)  # for increasing distance_to_DSB values
    idx <- which.max(cumsum(tmp$score) >= threshold * max(cumsum(tmp$score)))
    start(out[n]) <- start(tmp[idx])
    
    # find tract end on - strand
    roi <- DSBs[n]
    strand(roi) <- "+"
    roi <- resize(x = roi, width = max_DSB_dist, fix = "start")
    tmp <- subsetByIntersect(subject =GRanges, query = roi)
    idx <- which.max(cumsum(tmp$score) >= threshold * max(cumsum(tmp$score)))
    end(out[n]) <- end(tmp[idx])
  }
  
  return(out)
  
}


# get resection tract end distances from DSBs
# argument: GRanges
# result: data.frame
resection_extend_from_DSBs <- function(GRanges, DSBs){
  data.frame(DSB_id = as.character(DSBs),
             upstream = distance(x = DSBs, y = resize(x = GRanges, width = 1, fix = "start")),
             downstream = distance(x = DSBs, y = resize(x = GRanges, width = 1, fix = "end")))
}


# process data for plotting ===============================================

# read S1-seq coverage
load(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Data/S1-seq/LSY4602-20C_1/LSY4602-20C_1.RData")
load(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Data/S1-seq/LSY4377-12B_4377-15A_merged/LSY4377-12B_LSY4377-15A_merged.RData")

# define DSBs and surrounding regions
DSBs <- SrfIcs[- c(9, 17)]
# let's exclude chrIX:22006 and chrXV:27760 regions from analyses
# because they are in duplicates regions where some alignments could not be 
# unambiguously assigned and, thus, were filtered out
roi <- DSB_regions(DSBs = DSBs, region_width = 10000)

# retain only -/+ 5 kb regions around SrfIcs (for faster data processing)
LSY4602_20C_1_S1_seq <- subsetByIntersect(subject = LSY4602_20C_1_S1_seq, query = roi)
LSY4377_12B_1_S1_seq <- subsetByIntersect(subject = LSY4377_12B_1_merged_S1_seq, query = roi)

# make nt-resolved GRanges objects and sort
LSY4602_20C_1 <- sort(as_nt_resolved_GRanges(LSY4602_20C_1_S1_seq), ignore.strand = TRUE)
LSY4377_12B_1 <- sort(as_nt_resolved_GRanges(LSY4377_12B_1_S1_seq), ignore.strand = TRUE)


# add distance to nearest SrfIcs ------------------------------------------
# (for calculating averages over multiple DSBs)
LSY4602_20C_1 <- add_distance_to_nearest_DSB(subject = LSY4602_20C_1, DSBs = DSBs, add_DSB_id = TRUE)
LSY4377_12B_1 <- add_distance_to_nearest_DSB(subject = LSY4377_12B_1, DSBs = DSBs, add_DSB_id = TRUE)


# find furthest resection tract ends --------------------------------------
# (to avoid evaluating background signals e.g. for sequence-specific analyses [see later scripts])
resection_extend_LSY4602_20C_1 <- find_resection_extend(GRanges = LSY4602_20C_1, DSBs = DSBs, max_DSB_dist = 750)
resection_extend_LSY4377_12B_1 <- find_resection_extend(GRanges = LSY4377_12B_1, DSBs = DSBs, max_DSB_dist = 750)

# save
save(list = c("resection_extend_LSY4602_20C_1"), file = paste0(save_dir, "/LSY4602_20C_1_resection_extend.RData"))


# let's check resection extend over time
resection_extend_over_time <- resection_extend_from_DSBs(GRanges = c(resection_extend_LSY4602_20C_1, resection_extend_LSY4377_12B_1), DSBs = rep(DSBs, 2))
resection_extend_over_time <- cbind(resection_extend_over_time, sample = rep(c("LSY4602_20C_1", "LSY4377_12B_1"), rep(length(DSBs), 2)))
resection_extend_over_time$downstream <- - resection_extend_over_time$downstream

aggregate(upstream ~ sample, data = resection_extend_over_time, FUN = summary)
aggregate(downstream ~ sample, data = resection_extend_over_time, FUN = summary)


# plot furthest resection end points --------------------------------------

idx <- order(width(resection_extend_LSY4602_20C_1))

pdf(file = "tmp.pdf", width = 7, height = 4)
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(0.3, -0.5, 0.7, 2), tcl = -0.3, mgp = c(3.5, 0.6, 0), las = 1)
  plot(x = 1:length(DSBs), y = resection_extend_over_time$upstream[resection_extend_over_time$sample == "LSY4377_12B_1"][idx], ylim = range(resection_extend_over_time[, c(2:3)]),
       pch = 21, col = NA, bg = adjustcolor(col = JFly_colors[1], alpha.f = 0.5), 
       xlab = NA, xaxt = "n", ylab = "Resection Extend from DSBs [nt]")
  points(x = 1:length(DSBs), y = resection_extend_over_time$downstream[resection_extend_over_time$sample == "LSY4377_12B_1"][idx],
         pch = 21, col = NA, bg = adjustcolor(col = JFly_colors[1], alpha.f = 0.5))
  
  text(x = 1:length(DSBs), y = 1.05 * par("usr")[3], labels = resection_extend_over_time$DSB_id[resection_extend_over_time$sample == "LSY4377_12B_1"][idx], 
       srt = 45, xpd = TRUE, adj = c(1,1))
  
  points(x = 1:length(DSBs), y = resection_extend_over_time$upstream[resection_extend_over_time$sample == "LSY4602_20C_1"][idx],
         pch = 21, col = NA, bg = adjustcolor(col = JFly_colors[2], alpha.f = 0.5))
  points(x = 1:length(DSBs), y = resection_extend_over_time$downstream[resection_extend_over_time$sample == "LSY4602_20C_1"][idx],
         pch = 21, col = NA, bg = adjustcolor(col = JFly_colors[2], alpha.f = 0.5))
  
  abline(h = 0, lty = "dashed")
  
  legend(x = "top", inset = -0.27, ncol = 3, xpd = TRUE,
         pch = 21, col = NA, pt.bg = adjustcolor(col = JFly_colors[1:2], alpha.f = 0.5), 
         legend = c("WT", expression(italic("ku70"*Delta))))
dev.off()
GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/Resection_extend_from_DSBs.pdf"))


# S1-seq sum per DSB over time --------------------------------------------

tmp <- subsetByIntersect(subject = LSY4377_12B_1, query = resection_extend_LSY4377_12B_1)
tmp <- sapply(X = split(tmp, ~ DSB_id), FUN = function(x) {sum(x$score)})
S1_seq_frac <- data.frame(DSB_id = names(tmp), sample = "LSY4377_12B_1", score = tmp / sum(tmp))

tmp <- subsetByIntersect(subject = LSY4602_20C_1, query = resection_extend_LSY4602_20C_1)
tmp <- sapply(X = split(tmp, ~ DSB_id), FUN = function(x) {sum(x$score)})
S1_seq_frac <- rbind(S1_seq_frac, data.frame(DSB_id = names(tmp), sample = "LSY4602_20C_1", score = tmp / sum(tmp)))

# plot S1-seq sum grouped by DSB kinetics ---------------------------------

# Faster DSB formation will result in more ligatable ends and an over-representation of this region at early time points.
idx <- order(S1_seq_frac$score[S1_seq_frac$sample == "LSY4377_12B_1"])
# inspection of the plot below allows to group into slow, middle, and fast DSB formation kinetics.
idx_slow <- idx[1:3]
idx_middle <- idx[4:10]
idx_fast <- idx[11:19]

pdf(file = "tmp.pdf", width = 7, height = 4)
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(0.3, -0.5, 0.7, 2), tcl = -0.3, mgp = c(3.5, 0.6, 0), las = 1)
  shift <- 0.05
  plot(x = 1:length(DSBs) - 1 * shift, y = S1_seq_frac$score[S1_seq_frac$sample == "LSY4377_12B_1"][idx], ylim = range(S1_seq_frac$score),
       pch = 21, col = NA, bg = adjustcolor(col = JFly_colors[1], alpha.f = 0.5), 
       xlab = NA, xaxt = "n", ylab = "Fraction of Summed S1-seq")
  
  text(x = 1:length(DSBs), y = par("usr")[3] - 0.025 * (par("usr")[4] - par("usr")[3]), srt = 45, 
       labels = S1_seq_frac$DSB_id[S1_seq_frac$sample == "LSY4377_12B_1"][idx], xpd = TRUE, adj = c(1,1), font = c(2, rep(1, 6), 2, 1, 2, 2, rep(1,3), 2, rep(1, 4)))
  
  points(x = 1:length(DSBs) + 0 * shift, y = S1_seq_frac$score[S1_seq_frac$sample == "LSY4602_20C_1"][idx],
         pch = 21, col = NA, bg = adjustcolor(col = JFly_colors[2], alpha.f = 0.5))
  
  abline(h = 1 / 19, lty = "dashed", col = "gray")
  
  legend(x = "top", inset = -0.27, ncol = 3, xpd = TRUE,
         pch = 21, col = NA, pt.bg = adjustcolor(col = JFly_colors[1:2], alpha.f = 0.5), 
         legend = c("WT", expression(italic("ku70"*Delta))))
  
  abline(v = length(idx_slow) + 0.5, lty = "dashed")
  text(x = 0.5 * length(idx_slow) + 0.5, y = par("usr")[4] - 0.01 * (par("usr")[4] - par("usr")[3]), pos = 1, labels = "Slow")
  abline(v = length(idx_slow) + length(idx_middle) + 0.5, lty = "dashed")
  text(x = length(idx_slow) + 0.5 * length(idx_middle) + 0.5, y = par("usr")[4] - 0.01 * (par("usr")[4] - par("usr")[3]), pos = 1, labels = "Middle")
  text(x = length(idx_slow) + length(idx_middle) + 0.5 * length(idx_fast) + 0.5, 
       y = par("usr")[4] - 0.01 * (par("usr")[4] - par("usr")[3]), pos = 1, labels = "Fast")
  
dev.off()
GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/S1-seq_fraction_per_DSB.pdf"))


# add DSB kinetics categories to S1-seq data sets -------------------------

add_DSB_kinetics_rank <- function(GRanges, DSB_rank_map){
  GRanges$DSB_kinetics_rank <- NA
  for(n in 1:nrow(DSB_rank_map)){
    GRanges[GRanges$DSB_id == DSB_rank_map$DSB_id[n]]$DSB_kinetics_rank <- DSB_rank_map$rank[n]
  }
  return(GRanges)
}

DSB_rank_map <- S1_seq_frac[S1_seq_frac$sample == "LSY4602_20C_1", c(1, 3)]
DSB_rank_map <- DSB_rank_map[order(DSB_rank_map$score, decreasing = TRUE), ]
DSB_rank_map$rank <- 1:19
LSY4602_20C_1 <- add_DSB_kinetics_rank(GRanges = LSY4602_20C_1, DSB_rank_map = DSB_rank_map)

LSY4377_12B_1 <- add_DSB_kinetics_rank(GRanges = LSY4377_12B_1, DSB_rank_map = DSB_rank_map)

# save data
LSY4602_20C_1_S1_seq <- LSY4602_20C_1
save(list = c("LSY4602_20C_1_S1_seq"), file = paste0(save_dir, "/LSY4602_20C_1_around_DSBs.RData"))


# analyze fraction of unprocessed DSBs ------------------------------------

tmp <- subsetByIntersect(subject = LSY4377_12B_1, query = resection_extend_LSY4377_12B_1)
DSB_kinetics_rank <- sapply(X = split(tmp, ~ DSB_id), FUN = function(x) {unique(x$DSB_kinetics_rank)})
tmp <- sapply(X = split(tmp, ~ DSB_id), FUN = function(x) {sum(x$score[x$distance_to_DSB <= 1]) / sum(x$score)})
unprocessed_DSB_frac <- data.frame(DSB_id = names(tmp), sample = "LSY4377_12B_1", score = tmp, DSB_kinetics_rank = DSB_kinetics_rank)

tmp <- subsetByIntersect(subject = LSY4602_20C_1, query = resection_extend_LSY4602_20C_1)
DSB_kinetics_rank <- sapply(X = split(tmp, ~ DSB_id), FUN = function(x) {unique(x$DSB_kinetics_rank)})
tmp <- sapply(X = split(tmp, ~ DSB_id), FUN = function(x) {sum(x$score[x$distance_to_DSB <= 1]) / sum(x$score)})
unprocessed_DSB_frac <- rbind(unprocessed_DSB_frac, data.frame(DSB_id = names(tmp), sample = "LSY4602_20C_1", score = tmp, DSB_kinetics_rank = DSB_kinetics_rank))


# plot S1-seq frac vs. unprocessed DSBs -----------------------------------

all.equal(S1_seq_frac$DSB_id[S1_seq_frac$sample == "LSY4602_20C_1"], unprocessed_DSB_frac$DSB_id[unprocessed_DSB_frac$sample == "LSY4602_20C_1"])

DSB_metrics <- S1_seq_frac[S1_seq_frac$sample == "LSY4602_20C_1", ]
colnames(DSB_metrics)[3] <- "S1_seq_frac"
DSB_metrics$unprocessed <- unprocessed_DSB_frac$score[unprocessed_DSB_frac$sample == "LSY4602_20C_1"]

DSB_metrics <- DSB_metrics[order(DSB_metrics$S1_seq_frac), ]

pdf(file = "tmp.pdf", width = 4.5, height = 3)
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(1, -0.5, 3.8, 1.7), tcl = -0.3, mgp = c(3.5, 0.6, 0), las = 1)
  plot(x = DSB_metrics$S1_seq_frac, y = DSB_metrics$unprocessed, 
       pch = 21, col = NA, bg = adjustcolor(col = JFly_colors[1], alpha.f = 0.5),
       xlab = NA, ylab = "Unprocessed DSB Fraction")
  title(xlab = "Fraction of Summed S1-seq", line = 2.75)
  x_th <- c(mean(DSB_metrics$S1_seq_frac[3:4]), mean(DSB_metrics$S1_seq_frac[10:11]))
  abline(v = x_th, lty = "dashed", col = "gray")
  text(x = c(par("usr")[1] + 0.5 * (x_th[1] - par("usr")[1]), mean(x_th)), y = par("usr")[3] + 0.01 * (par("usr")[4] - par("usr")[3]), pos = 3, labels = c("Low", "Middle"), col = "gray")
  text(x = par("usr")[2] - 0.5 * (par("usr")[2] - x_th[2]), y = par("usr")[4] - 0.01 * (par("usr")[4] - par("usr")[3]), pos = 1, labels = "Fast", col = "gray")
dev.off()
GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/S1-seq_fraction_vs_unprocessed_DSB_fraction.pdf"))


# plot unprocessed DSB fraction grouped by DSB kinetics -------------------
idx_slow <- idx[1:3]
idx_middle <- idx[4:10]
idx_fast <- idx[11:19]
idx <- c(idx_slow, idx_middle, idx_fast)

all(unprocessed_DSB_frac$DSB_id[unprocessed_DSB_frac$sample == "LSY4377_12B_1"][idx] == unprocessed_DSB_frac$DSB_id[unprocessed_DSB_frac$sample == "LSY4602_20C_1"][idx])

pdf(file = "tmp.pdf", width = 7, height = 4)
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(0.3, -0.5, 0.7, 2), tcl = -0.3, mgp = c(3.5, 0.6, 0), las = 1)
shift <- 0.05
plot(x = 1:length(DSBs) - 1 * shift, y = unprocessed_DSB_frac$score[unprocessed_DSB_frac$sample == "LSY4377_12B_1"][idx], ylim = range(unprocessed_DSB_frac$score),
     pch = 21, col = NA, bg = adjustcolor(col = JFly_colors[1], alpha.f = 0.5), 
     xlab = NA, xaxt = "n", ylab = "Unprocessed DSB Fraction")

text(x = 1:length(DSBs), y = par("usr")[3] - 0.025 * (par("usr")[4] - par("usr")[3]), srt = 45, 
     labels = unprocessed_DSB_frac$DSB_id[unprocessed_DSB_frac$sample == "LSY4377_12B_1"][idx], xpd = TRUE, adj = c(1,1))

points(x = 1:length(DSBs) + 0 * shift, y = unprocessed_DSB_frac$score[unprocessed_DSB_frac$sample == "LSY4602_20C_1"][idx],
       pch = 21, col = NA, bg = adjustcolor(col = JFly_colors[2], alpha.f = 0.5))

legend(x = "top", inset = -0.27, ncol = 3, xpd = TRUE,
       pch = 21, col = NA, pt.bg = adjustcolor(col = JFly_colors[1:2], alpha.f = 0.5), 
       legend = c("WT", expression(italic("ku70"*Delta))))

abline(v = length(idx_slow) + 0.5, lty = "dashed")
text(x = 0.5 * length(idx_slow) + 0.5, y = par("usr")[4] - 0.01 * (par("usr")[4] - par("usr")[3]), pos = 1, labels = "Slow")
abline(v = length(idx_slow) + length(idx_middle) + 0.5, lty = "dashed")
text(x = length(idx_slow) + 0.5 * length(idx_middle) + 0.5, y = par("usr")[4] - 0.01 * (par("usr")[4] - par("usr")[3]), pos = 1, labels = "Middle")
text(x = length(idx_slow) + length(idx_middle) + 0.5 * length(idx_fast) + 0.5, 
     y = par("usr")[4] - 0.01 * (par("usr")[4] - par("usr")[3]), pos = 1, labels = "Fast")

dev.off()
GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/Unprocessed_DSB_fraction.pdf"))


# plot distributions grouped by DSB kinetics

# add slow/middle/fast category
unprocessed_DSB_frac$DSB_kinetics <- NA
unprocessed_DSB_frac$DSB_kinetics[unprocessed_DSB_frac$DSB_kinetics_rank %in% 19:17] <- "slow"
unprocessed_DSB_frac$DSB_kinetics[unprocessed_DSB_frac$DSB_kinetics_rank %in% 16:10] <- "middle"
unprocessed_DSB_frac$DSB_kinetics[unprocessed_DSB_frac$DSB_kinetics_rank %in% 9:1] <- "fast"

# convert group ids to factors to control plotting
unprocessed_DSB_frac$DSB_kinetics <- factor(x = unprocessed_DSB_frac$DSB_kinetics, levels = c("slow", "middle", "fast"))

pdf(file = "tmp.pdf", width = 4.5, height = 3.5)
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(0.5, -0.5, 3.9, 2), tcl = -0.3, mgp = c(3.5, 0.6, 0), las = 1)

at <- c(1:3, 5:7)
bp <- boxplot(at = at, score ~ DSB_kinetics + sample, data = unprocessed_DSB_frac, outline = FALSE, col = NA, border = "gray",
              xlab = "DSB Kinetics", ylab = "Unprocessed DSB Fraction", names = rep(NA, 6), xaxt = "n")
text(x = at, y = par("usr")[3] - 0.025 * (par("usr")[4] - par("usr")[3]), 
     srt = 45, adj = c(1,1), labels = rep(c("Slow", "Middle", "Fast"), 4), xpd = TRUE)

abline(v = 4)

stripchart(at = at, score ~ DSB_kinetics + sample, data = unprocessed_DSB_frac, 
           vertical = TRUE, pch = 20, method = "jitter", add = TRUE, col = rep(adjustcolor(col = JFly_colors[1:2], alpha.f = 0.5), rep(3, 2)))

text(x = c(2, 6), y = par("usr")[4] - 0.01 * (par("usr")[4] - par("usr")[3]), pos = 1, labels = c("WT", expression(italic("ku70"*Delta))))

dev.off()
GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/Unprocessed_DSB_fraction_box_plots.pdf"))