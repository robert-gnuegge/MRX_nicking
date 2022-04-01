# info --------------------------------------------------------------------
# purpose: generate S1-seq and MNase-seq coverage plots around genomic SrfIcs
# author: Robert Gnuegge (robert.gnuegge@gmail.com)
# date: 01/28/22
# version: 2.0


# preamble ----------------------------------------------------------------

# set wd to this file's location
wd_path <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(wd_path)

# libraries
library(GenomicRanges)
library(Gviz)
options(ucscChromosomeNames=FALSE)  # for using custom chromosome names (e.g. "micron")

# load helper files and functions
source(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Analysis_scripts/Misc/Genomic_helper_functions.R")
source(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Analysis_scripts/Misc/Misc_helper_functions.R")
source(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Analysis_scripts/Misc/JFly_colors.R")


# function definitions ====================================================

# make Gviz track with separately marked fw and rev scores
# argument: GRanges object, color definitions, ... (e.g. type, name)
# result: OverlayTrack
DataTrack_fw_rev <- function(GRanges, fw_col, rev_col, ...){
  fw <- DataTrack(range = GRanges[strand(GRanges) == "+"], col = fw_col, ...)
  rev <- DataTrack(range = GRanges[strand(GRanges) == "-"], col = rev_col, ...)
  OverlayTrack(trackList = list(fw, rev))
}


# process S1-seq data for plotting ========================================

DSB <- GRanges("chrIII:295660")
roi <- DSB_regions(DSBs = DSB, region_width = 10000)

# read S1-seq coverage
load(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Data/S1-seq/LSY4753-3B_4754-3A/LSY4753-3B_4754-3A.RData")

# retain only -/+ 5 kb regions around SrfIcs (for faster data processing)
S1_seq_HMR <- subsetByIntersect(subject = LSY4753_3B_2_S1_seq, query = roi)
S1_seq_hmr <- subsetByIntersect(subject = LSY4754_3A_2_S1_seq, query = roi)
S1_seq_HMR$score[S1_seq_HMR$score > 1500] <- 1500  # adjust for better visualization
S1_seq_hmr$score[S1_seq_hmr$score > 1500] <- 1500

# make nt-resolved GRanges objects and sort
S1_seq_HMR <- sort(as_nt_resolved_GRanges(S1_seq_HMR), ignore.strand = TRUE)
S1_seq_hmr <- sort(as_nt_resolved_GRanges(S1_seq_hmr), ignore.strand = TRUE)

# apply Hanning smoother
n <- 51  # Hanning window size
S1_seq_HMR_smooth <- as_smoothed_GRanges(GRanges = S1_seq_HMR, hanning_window_size = n)
S1_seq_hmr_smooth <- as_smoothed_GRanges(GRanges = S1_seq_hmr, hanning_window_size = n)


# process MNase-seq data for plotting =====================================

# read MNase-seq coverage
load(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Data/MNase-seq/LSY4753-3B_4754-3A/LSY4753-3B_4754-3A_trimmed.RData")

# retain only -/+ 5 kb regions around SrfIcs (for faster data processing)
MNase_seq_HMR_0 <- subsetByIntersect(subject = LSY4753_3B_0_MNase_seq_trimmed, query = roi)
MNase_seq_HMR_2 <- subsetByIntersect(subject = LSY4753_3B_2_MNase_seq_trimmed, query = roi)
MNase_seq_hmr_0 <- subsetByIntersect(subject = LSY4754_3A_0_MNase_seq_trimmed, query = roi)
MNase_seq_hmr_2 <- subsetByIntersect(subject = LSY4754_3A_2_MNase_seq_trimmed, query = roi)

# apply runmed smoother
MNase_seq_HMR_0$score <- runmed(x = MNase_seq_HMR_0$score, k = 31)
MNase_seq_HMR_2$score <- runmed(x = MNase_seq_HMR_2$score, k = 31)
MNase_seq_hmr_0$score <- runmed(x = MNase_seq_hmr_0$score, k = 31)
MNase_seq_hmr_2$score <- runmed(x = MNase_seq_hmr_2$score, k = 31)

# plotting ================================================================

plot_dir <- "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Figures/13_Nucleosome_occupancy_impact/LSY4753-3B_4754-3A"

# set global Gviz parameters for plotting
options(Gviz.scheme="default")
scheme <- getScheme()  # copy current scheme
scheme$GdObject$rotation.title <- 90
scheme$GdObject$fontcolor.title <- "black"
scheme$GdObject$background.title <- "white"
scheme$GdObject$cex.title <- 0.75
scheme$GdObject$col.axis <- "black"
addScheme(scheme, "MyScheme")  # define new scheme
options(Gviz.scheme = "MyScheme")  # set new scheme

# plotting colors
fw_col <- JFly_colors[2]
rev_col <- JFly_colors[3]
MNase_seq_0_col <- gray(level = 0.67)
MNase_seq_col <- gray(level = 0.9)
MNase_seq_border_col <- gray(level = 0.75)

Y_axis_ranges <- data.frame()  # initialize for collecting

# define name for data and figure saving
DSB <- GRanges("chrIII:295660")
DSB_locus <- gsub(pattern = ":", replacement = "_", x = as.character(DSB))
cat("\nPlotting", DSB_locus, "...")

# define plot area
roi <- DSB_regions(DSBs = DSB, region_width = 3000)

# AnnotationTrack ---------------------------------------------------------
# data downloaded from yeastmine.yeastgenome.org and curated with Adjust_chromosomal_features.R

load(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Data/S_cerevisiae_chromosomal_features/S_cerevisiae_genome_features.RData")

# adjust for use with AnnotationTrack
names(mcols(all_features))[1] <- "id"
all_features <- all_features[!(all_features$id == "") & all_features$type %in% c("ORF", "tRNA gene")]

# remove deleted MATA1 and YCR097W-A
all_features <- all_features[!all_features$id %in% c("HMRA1", "YCR097W-A")]

# add HMR-I and Citrine
tmp <- GRanges(seqnames = "chrIII", ranges = IRanges(start = c(293502, 294805), end = c(294221, 294893)), strand = c("+", "*"))
mcols(tmp) <- data.frame(id = c("Citrine", "HMR-I"), systematic_name = NA, type = c("ORF", "ARS"), qualifier = "Verified")
all_features <- sort(c(all_features, tmp), ignore.strand = TRUE)

subsetByOverlaps(x = all_features, ranges = roi)

AT <- AnnotationTrack(range = all_features, name = NULL, featureAnnotation = "id", cex = 0.5,
                      arrowHeadMaxWidth = 10, fill = "white", col = "gray", fontcolor.item = "black")
# add DSB to AnnotationTrack
AT_w_DSBs <- HighlightTrack(trackList = list(AT), range = DSB, col = JFly_colors[8], fill = NA, inBackground = FALSE)


# S1-seq DataTracks -------------------------------------------------------

# get GRanges in plot area, calculate ylim, and construct DataTracks
S1_seq_HMR_roi <- subsetByOverlaps(x = S1_seq_HMR, ranges = roi)
S1_seq_hmr_roi <- subsetByOverlaps(x = S1_seq_hmr, ranges = roi)

ylim <- range(c(S1_seq_HMR_roi, S1_seq_hmr_roi)$score)
ylim <- c(min(ylim) - 0.05 * max(pretty(ylim)), max(pretty(ylim)))  # adjust for prettier plotting

DT_S1_seq_HMR <- DataTrack_fw_rev(GRanges = GRanges_zero_to_NA(S1_seq_HMR), fw_col = adjustcolor(col = fw_col, alpha.f = 0.5), rev_col = adjustcolor(col = rev_col, alpha.f = 0.5), type = "h", name = "HMR", ylim = ylim)
DT_S1_seq_hmr <- DataTrack_fw_rev(GRanges = GRanges_zero_to_NA(S1_seq_hmr), fw_col = adjustcolor(col = fw_col, alpha.f = 0.5), rev_col = adjustcolor(col = rev_col, alpha.f = 0.5), type = "h", name = "hmr", ylim = ylim)

S1_seq_HMR_smooth_roi <- subsetByOverlaps(x = S1_seq_HMR_smooth, ranges = roi)
S1_seq_hmr_smooth_roi <- subsetByOverlaps(x = S1_seq_hmr_smooth, ranges = roi)

ylim_smooth <- range(c(S1_seq_HMR_smooth_roi, S1_seq_hmr_smooth_roi)$score)
ylim_smooth <- c(min(ylim_smooth) - 0.05 * max(pretty(ylim_smooth)), max(pretty(ylim_smooth)))  # adjust for prettier plotting

DT_S1_seq_HMR_smooth <- DataTrack_fw_rev(GRanges = GRanges_zero_to_NA(S1_seq_HMR_smooth), fw_col = fw_col, rev_col = rev_col, type = "l", name = "pho4-SA", ylim = ylim_smooth)
DT_S1_seq_hmr_smooth <- DataTrack_fw_rev(GRanges = GRanges_zero_to_NA(S1_seq_hmr_smooth), fw_col = fw_col, rev_col = rev_col, type = "l", name = "HMR", ylim = ylim_smooth)


# MNase-seq DataTracks ----------------------------------------------------

MNase_seq_hmr_0_roi <- subsetByOverlaps(x = MNase_seq_hmr_0, ranges = roi) 
MNase_seq_hmr_2_roi <- subsetByOverlaps(x = MNase_seq_hmr_2, ranges = roi)
MNase_seq_HMR_0_roi <- subsetByOverlaps(x = MNase_seq_HMR_0, ranges = roi)
MNase_seq_HMR_2_roi <- subsetByOverlaps(x = MNase_seq_HMR_2, ranges = roi)
 
ylim_MNase <- range(c(MNase_seq_hmr_0_roi, MNase_seq_HMR_0_roi, MNase_seq_hmr_2_roi, MNase_seq_HMR_2_roi)$score)
ylim_MNase <- c(min(ylim_MNase) - 0.05 * max(pretty(ylim_MNase)), max(pretty(ylim_MNase)))  # adjust for prettier plotting

DT_MNase_seq_HMR_0 <- DataTrack(range = MNase_seq_HMR_0, type = "polygon", col = MNase_seq_0_col, fill.mountain = rep(NA, 2), name = "HMR", ylim = ylim_MNase)
DT_MNase_seq_hmr_0 <- DataTrack(range = MNase_seq_hmr_0, type = "polygon", col = MNase_seq_0_col, fill.mountain = rep(NA, 2), name = "hmr", ylim = ylim_MNase)
DT_MNase_seq_HMR_2 <- DataTrack(range = MNase_seq_HMR_2, type = "polygon", col = MNase_seq_border_col, fill.mountain = rep(MNase_seq_col, 2), name = "HMR", ylim = ylim_MNase)
DT_MNase_seq_hmr_2 <- DataTrack(range = MNase_seq_hmr_2, type = "polygon", col = MNase_seq_border_col, fill.mountain = rep(MNase_seq_col, 2), name = "hmr", ylim = ylim_MNase)


# plot --------------------------------------------------------------------

pdf(file = "tmp.pdf", width = 3.5, height = 2.5)
  plotTracks(trackList = list(OverlayTrack(trackList = list(DT_MNase_seq_HMR_2, DT_MNase_seq_HMR_0, DT_S1_seq_HMR_smooth, DT_S1_seq_HMR)),
                              OverlayTrack(trackList = list(DT_MNase_seq_hmr_2, DT_MNase_seq_hmr_0, DT_S1_seq_hmr_smooth, DT_S1_seq_hmr)),
                              AT_w_DSBs), 
             from = start(roi), to = end(roi), chromosome = seqnames(roi), margin = 1, innerMargin = 0)
dev.off()
GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/", DSB_locus, ".pdf"))

# record ylims (for manual scale addition)
Y_axis_ranges <- rbind(Y_axis_ranges,
                       data.frame(DSB_locus = DSB_locus, 
                                  y_raw_min = ylim[1], y_raw_max = ylim[2],
                                  y_smooth_min = ylim_smooth[1], y_smooth_max = ylim_smooth[2],
                                  y_MNase_min = ylim_MNase[1], y_MNase_max = ylim_MNase[2]))

write.table(x = Y_axis_ranges, file = paste0(plot_dir, "/Y_axis_ranges.txt"), row.names = FALSE)
