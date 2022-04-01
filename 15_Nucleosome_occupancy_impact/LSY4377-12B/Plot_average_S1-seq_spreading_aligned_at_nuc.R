# info --------------------------------------------------------------------
# purpose: plot average S1-seq spreading aligned at nucleosomes
# author: Robert Gnuegge (robert.gnuegge@gmail.com)
# date: 01/24/22
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

plot_dir <- "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Figures/13_Nucleosome_occupancy_impact/LSY4377-12B_4377-15A/S1-seq_aligned_at_nucleosomes"


# function definitions ----------------------------------------------------

# complete nucleosome center information (strand, nucleosome number relative to DSB)
# argument: GRanges objects
# result: GRanges object
add_nuc_numbers_and_strand <- function(nuc_centers, ideal_nuc_centers, min_dist_to_DSB = 83){
  tmp <- nuc_centers
  strand(tmp) <- ifelse(test = tmp$distance_to_DSB < 0, yes = "-", no = "+")
  tmp <- tmp[abs(tmp$distance_to_DSB) > min_dist_to_DSB]
  hits <- distanceToNearest(x = tmp, subject = ideal_nuc_centers)
  mcols(tmp)$nuc_number <- NA
  mcols(tmp[queryHits(hits)])$nuc_number <- mcols(ideal_nuc_centers[subjectHits(hits)])$nuc_number
  return(tmp)
}

# add distance to selected nucleosome (or nearest-to-selected nucleosome)
# argument: GRanges objects and numeric
# result: GRanges object
add_distance_to_nucleosome <- function(GRanges, nuc_centers, roi, nuc, nuc_dist = 165){
  out <- GRanges()
  for(n in 1:length(roi)){
    nuc_region <- subsetByOverlaps(x = nuc_centers, ranges = roi[n])
    nuc_pos <- nuc_region[which.min(abs(nuc_region$nuc_number) - nuc)]
    tmp <- subsetByOverlaps(x = GRanges, ranges = roi[n])
    if(unique(strand(tmp)) == "-"){
      tmp$dist_from_nuc <- start(nuc_pos) - start(tmp) + (abs(nuc_pos$nuc_number) - nuc) * nuc_dist
    }else{
      tmp$dist_from_nuc <- start(tmp) - start(nuc_pos) + (nuc_pos$nuc_number - nuc) * nuc_dist
    }
    tmp$nuc_number <- nuc_pos$nuc_number 
    out <- c(out, tmp)
  }
  return(out)
}

# process S1-seq and MNase-seq data ---------------------------------------

# read S1-seq coverage
load(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Data/S1-seq/LSY4377-12B_4377-15A_rep2/LSY4377-12B_LSY4377-15A_rep2_around_DSBs.RData")

# keep only regions with correct orientation w.r.t DSBs
DSBs <- SrfIcs[-c(9, 17)]
roi <- DSB_regions(DSBs = DSBs, region_width = 6000, up_rev_down_fw = TRUE)
S1_seq_0 <- subsetByIntersect(subject = LSY4377_12B_0_rep2_S1_seq, query = roi)
S1_seq_1 <- subsetByIntersect(subject = LSY4377_12B_1_rep2_S1_seq, query = roi)
S1_seq_2 <- subsetByIntersect(subject = LSY4377_12B_2_rep2_S1_seq, query = roi)
S1_seq_4 <- subsetByIntersect(subject = LSY4377_12B_4_rep2_S1_seq, query = roi)

# read MNase-seq coverage
load(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Data/MNase-seq/LSY4377-12B_4377-15A/LSY4377-12B_LSY4377-15A_trimmed.RData")

# retain only regions around SrfIcs (for faster data processing)
roi <- DSB_regions(DSBs = DSBs, region_width = 6000, up_rev_down_fw = TRUE)
MNase_seq_0 <- as_nt_resolved_GRanges(GRanges = subsetByIntersect(subject = LSY4377_12B_0_MNase_seq_trimmed, query = roi))
MNase_seq_1 <- as_nt_resolved_GRanges(GRanges = subsetByIntersect(subject = LSY4377_12B_1_MNase_seq_trimmed, query = roi))
MNase_seq_2 <- as_nt_resolved_GRanges(GRanges = subsetByIntersect(subject = LSY4377_12B_2_MNase_seq_trimmed, query = roi))
MNase_seq_4 <- as_nt_resolved_GRanges(GRanges = subsetByIntersect(subject = LSY4377_12B_4_MNase_seq_trimmed, query = roi))


# process nucleosome center data ------------------------------------------

# read nucleosome centers
load(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Data/MNase-seq/LSY4377-12B_4377-15A/Nuc_centers_around_SrfIcs.RData")

# nucleosome properties (according to Jansen et al., 2011; pmid: 21646431)
nuc_width <- 147
nuc_dist <- 165  # nuc_width + 18 (average linker length)

# define ideal nucleosome positions relative to DSBs
ideal_nuc_centers <- GRanges()
for(n in 0:19){
  tmp <- shift(x = DSBs, shift = -round((n + 0.5) * nuc_dist))
  tmp$nuc_number <- -(n + 1)
  ideal_nuc_centers <- c(ideal_nuc_centers, tmp)
  tmp <- shift(x = DSBs, shift = round((n + 0.5) * nuc_dist))
  tmp$nuc_number <- (n + 1)
  ideal_nuc_centers <- c(ideal_nuc_centers, tmp)
}
ideal_nuc_centers <- sort(ideal_nuc_centers)
strand(ideal_nuc_centers) <- ifelse(test = ideal_nuc_centers$nuc_number < 0, yes = "-", no = "+")

# add nucleosome numbers and strand information to Nuc_centers GRanges  
Nuc_centers_0 <- add_nuc_numbers_and_strand(nuc_centers = Nuc_centers_0, ideal_nuc_centers = ideal_nuc_centers)
Nuc_centers_1 <- add_nuc_numbers_and_strand(nuc_centers = Nuc_centers_1, ideal_nuc_centers = ideal_nuc_centers)
Nuc_centers_2 <- add_nuc_numbers_and_strand(nuc_centers = Nuc_centers_2, ideal_nuc_centers = ideal_nuc_centers)
Nuc_centers_4 <- add_nuc_numbers_and_strand(nuc_centers = Nuc_centers_4, ideal_nuc_centers = ideal_nuc_centers)


# plotting ================================================================

DSBs <- SrfIcs[-c(9, 17)]
roi <- DSB_regions(DSBs = DSBs, region_width = 6000, up_rev_down_fw = TRUE)

# function for plotting of nucleosome-aligned MNase-seq and S1-seq data in same plot
# arguments: numeric data.frames, doubles
# result: plot
plotting_function <- function(S1_seq, S1_seq_ylim, MNase_seq, MNase_seq_ylim, xlim, nuc, nuc_dist = 165){
  # start empty plot
  plot(x = NA, y = NA, xlim = xlim, ylim = MNase_seq_ylim, xaxt = "n",
       xlab = paste0("Distance from +", nuc, " dyade [nt]"), ylab = "Average MNase-seq [RPM]")
  axis(side = 1, at = 0:2 * 500)
  # plot MNase-seq
  y <- moving_average(x = MNase_seq$score, k = 51, keep = 0)
  polygon(x = c(MNase_seq$dist_from_nuc[1], MNase_seq$dist_from_nuc, MNase_seq$dist_from_nuc[length(MNase_seq$dist_from_nuc)]), y = c(0, y, 0), col = gray(level = 0.8), border = gray(level = 0.5))
  # add nucleosome labels 
  x <- 1:5
  segments(x0 = (x - nuc) * nuc_dist, y0 = 1.01 * par("usr")[4], x1 = (x - nuc) * nuc_dist, y1 = 0, col = gray(level = 0.4), lty = "dashed", xpd = TRUE)
  text(x = (x - nuc) * nuc_dist, y = par("usr")[4], labels = paste0("+", x), pos = 3, xpd = TRUE, col = gray(level = 0.4))
  # add S1-seq plot
  par(new = TRUE)
  y <- moving_average(x = S1_seq$score, k = 51, keep = 1)
  plot(x = S1_seq$dist_from_nuc, y = y, ylim = S1_seq_ylim, xlim = xlim, type = "l", col = JFly_colors[5], axes = FALSE, ann = FALSE)
  axis(side = 4, at = pretty(S1_seq_ylim), col = JFly_colors[5], col.axis = JFly_colors[5])
  mtext("Average S1-seq [RPM]", side = 4, line = 2.5, col = JFly_colors[5], las = 0)
}


# t = 1 -------------------------------------------------------------------
nuc <- 2
S1 <- add_distance_to_nucleosome(GRanges = S1_seq_1, nuc_centers = Nuc_centers_1, roi = roi, nuc = nuc)
MNase <- add_distance_to_nucleosome(GRanges = MNase_seq_1, nuc_centers = Nuc_centers_1, roi = roi, nuc = nuc)
idx <- S1$DSB_kinetics_rank < 20 & S1$dist_from_nuc > -(nuc - 0.5) * nuc_dist  # & S1$dist_from_nuc <= 1200
S1_agg_1 <- aggregate(score ~ dist_from_nuc, data = S1[idx], FUN = mean)
MNase_agg_1 <- aggregate(score ~ dist_from_nuc, data = MNase[idx], FUN = mean)

pdf(file = "tmp.pdf", width=5, height=3)
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(1.5, 0.6, 2.7, -1.5), las = 1, tcl = -0.3, mgp = c(2.5, 0.75, 0))
plotting_function(S1_seq = S1_agg_1, MNase_seq = MNase_agg_1, nuc = nuc, 
                  S1_seq_ylim = c(0,26), MNase_seq_ylim = c(0,32), xlim = c(-(nuc - 0.5) * nuc_dist, 1200))
dev.off()
GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/Aligned_S1_seq_and_MNase-seq_1h.pdf"))

# t = 2 -------------------------------------------------------------------
nuc <- 3
S1 <- add_distance_to_nucleosome(GRanges = S1_seq_2, nuc_centers = Nuc_centers_2, roi = roi, nuc = nuc)
MNase <- add_distance_to_nucleosome(GRanges = MNase_seq_2, nuc_centers = Nuc_centers_2, roi = roi, nuc = nuc)
idx <- S1$DSB_kinetics_rank < 20 & S1$dist_from_nuc > -(nuc - 0.5) * nuc_dist  # & S1$dist_from_nuc <= 1200
S1_agg_2 <- aggregate(score ~ dist_from_nuc, data = S1[idx], FUN = mean)
MNase_agg_2 <- aggregate(score ~ dist_from_nuc, data = MNase[idx], FUN = mean)

pdf(file = "tmp.pdf", width=5, height=3)
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(1.5, 0.6, 2.7, -1.5), las = 1, tcl = -0.3, mgp = c(2.5, 0.75, 0))
plotting_function(S1_seq = S1_agg_2, MNase_seq = MNase_agg_2, nuc = nuc, 
                  S1_seq_ylim = c(0,26), MNase_seq_ylim = c(0,32), xlim = c(-(nuc - 0.5) * nuc_dist, 1200))
dev.off()
GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/Aligned_S1_seq_and_MNase-seq_2h.pdf"))

# t = 4 -------------------------------------------------------------------
nuc <- 4
S1 <- add_distance_to_nucleosome(GRanges = S1_seq_4, nuc_centers = Nuc_centers_4, roi = roi, nuc = nuc)
MNase <- add_distance_to_nucleosome(GRanges = MNase_seq_4, nuc_centers = Nuc_centers_4, roi = roi, nuc = nuc)
idx <- S1$DSB_kinetics_rank < 20 & S1$dist_from_nuc > -(nuc - 0.5) * nuc_dist  # & S1$dist_from_nuc <= 1200
S1_agg_4 <- aggregate(score ~ dist_from_nuc, data = S1[idx], FUN = mean)
MNase_agg_4 <- aggregate(score ~ dist_from_nuc, data = MNase[idx], FUN = mean)

pdf(file = "tmp.pdf", width=5, height=3)
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(1.5, 0.6, 2.7, -1.5), las = 1, tcl = -0.3, mgp = c(2.5, 0.75, 0))
plotting_function(S1_seq = S1_agg_4, MNase_seq = MNase_agg_4, nuc = nuc, 
                  S1_seq_ylim = c(0,26), MNase_seq_ylim = c(0,32), xlim = c(-(nuc - 0.5) * nuc_dist, 1200))
dev.off()
GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/Aligned_S1_seq_and_MNase-seq_4h.pdf"))

# t = 0 -------------------------------------------------------------------
nuc <- 1
MNase <- add_distance_to_nucleosome(GRanges = MNase_seq_0, nuc_centers = Nuc_centers_0, roi = roi, nuc = nuc)
S1 <- add_distance_to_nucleosome(GRanges = S1_seq_0, nuc_centers = Nuc_centers_0, roi = roi, nuc = nuc)
idx <- S1$DSB_kinetics_rank < 20 & S1$dist_from_nuc > -(nuc - 0.5) * nuc_dist  # & S1$dist_from_nuc <= 1200
S1_agg_0 <- aggregate(score ~ dist_from_nuc, data = S1[idx], FUN = mean)
MNase_agg_0 <- aggregate(score ~ dist_from_nuc, data = MNase[idx], FUN = mean)

pdf(file = "tmp.pdf", width=5, height=3)
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(1.5, 0.6, 2.7, -1.5), las = 1, tcl = -0.3, mgp = c(2.5, 0.75, 0))
plotting_function(S1_seq = S1_agg_0, MNase_seq = MNase_agg_0, nuc = nuc, 
                  S1_seq_ylim = c(0,26), MNase_seq_ylim = c(0,32), xlim = c(-(nuc - 0.5) * nuc_dist, 1200))
dev.off()
GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/Aligned_S1_seq_and_MNase-seq_0h.pdf"))


# all t -------------------------------------------------------------------
k <- 51

pdf(file = "tmp.pdf", width=3.75, height=3)
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(2.3, 1, 2.9, -1), las = 1, tcl = -0.25, mgp = c(1.75, 0.4, 0))

# start empty plot
plot(x = NA, y = NA, xlim = c(-80, 1100), ylim = c(0,30.65), xaxt = "n",
     xlab = "Distance from +1 dyade [nt]", ylab = NA, yaxt = "n")
axis(side = 1, at = 0:4 * 250)
axis_col <- gray(level = 0.6)
axis(side = 4, at = pretty(c(0,30.65)), col = axis_col, col.axis = axis_col)
text(x = par("usr")[2]*1.2, y = 0.9* mean(par("usr")[3:4]), labels = "Average MNase-seq [RPM]", srt = -90, pos = 3, col = axis_col, xpd = TRUE)

# plot MNase-seq
MNase_seq <- MNase_agg_0
MNase_seq <- MNase_seq[MNase_seq$dist_from_nuc < 1100, ]
y <- moving_average(x = MNase_seq$score, k = k, keep = 0)
polygon(x = c(MNase_seq$dist_from_nuc[1], MNase_seq$dist_from_nuc, MNase_seq$dist_from_nuc[length(MNase_seq$dist_from_nuc)]), 
        y = c(0, y, 0), col = gray(level = 0.9), border = NA)

MNase_seq <- MNase_agg_1
MNase_seq <- MNase_seq[MNase_seq$dist_from_nuc < 1100 - nuc_dist, ]
y <- moving_average(x = MNase_seq$score, k = k, keep = 0)
polygon(x = c(MNase_seq$dist_from_nuc[1], MNase_seq$dist_from_nuc, MNase_seq$dist_from_nuc[length(MNase_seq$dist_from_nuc)]) + nuc_dist, 
        y = c(0, y, 0), col = gray(level = 0.8), border = NA)

MNase_seq <- MNase_agg_2
MNase_seq <- MNase_seq[MNase_seq$dist_from_nuc < 1100 - 2 * nuc_dist, ]
y <- moving_average(x = MNase_seq$score, k = k, keep = 0)
polygon(x = c(MNase_seq$dist_from_nuc[1], MNase_seq$dist_from_nuc, MNase_seq$dist_from_nuc[length(MNase_seq$dist_from_nuc)]) + 2 * nuc_dist, 
        y = c(0, y, 0), col = gray(level = 0.7), border = NA)

MNase_seq <- MNase_agg_4
MNase_seq <- MNase_seq[MNase_seq$dist_from_nuc < 1100 - 3 * nuc_dist, ]
y <- moving_average(x = MNase_seq$score, k = k, keep = 0)
polygon(x = c(MNase_seq$dist_from_nuc[1], MNase_seq$dist_from_nuc, MNase_seq$dist_from_nuc[length(MNase_seq$dist_from_nuc)]) + 3 * nuc_dist, 
        y = c(0, y, 0), col = gray(level = 0.6), border = NA)

# add nucleosome labels 
x <- 1:6
nuc <- 1
segments(x0 = (x - nuc) * nuc_dist, y0 = 1.01 * par("usr")[4], x1 = (x - nuc) * nuc_dist, y1 = 0, col = gray(level = 0.4), lty = "dashed", xpd = TRUE)
text(x = (x - nuc) * nuc_dist, y = par("usr")[4], labels = paste0("+", x), pos = 3, xpd = TRUE, col = gray(level = 0.4))

# add S1-seq plots
par(new = TRUE)

S1_seq <- S1_agg_1
S1_seq <- S1_seq[S1_seq$dist_from_nuc < 1100 - nuc_dist, ]
y <- moving_average(x = S1_seq$score, k = k, keep = 0)
plot(x = S1_seq$dist_from_nuc + nuc_dist, y = y, ylim = c(0, 26), xlim = c(-80, 1100), type = "l", col = JFly_colors[1], axes = FALSE, ann = FALSE)

S1_seq <- S1_agg_2
S1_seq <- S1_seq[S1_seq$dist_from_nuc < 1100 - 2 * nuc_dist, ]
y <- moving_average(x = S1_seq$score, k = k, keep = 0)
points(x = S1_seq$dist_from_nuc + 2 * nuc_dist, y = y, type = "l", col = JFly_colors[4])

S1_seq <- S1_agg_4
S1_seq <- S1_seq[S1_seq$dist_from_nuc < 1100 - 3 * nuc_dist, ]
y <- moving_average(x = S1_seq$score, k = k, keep = 0)
points(x = S1_seq$dist_from_nuc + 3 * nuc_dist, y = y, type = "l", col = JFly_colors[5])

axis(side = 2, at = pretty(c(0,26)))
mtext("Average S1-seq [RPM]", side = 2, line = 2, las = 0)

dev.off()
GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/Aligned_S1_seq_and_MNase-seq.pdf"))
