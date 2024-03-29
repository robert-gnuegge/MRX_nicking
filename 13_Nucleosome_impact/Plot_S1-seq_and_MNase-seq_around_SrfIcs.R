# info --------------------------------------------------------------------
# purpose: generate S1-seq and MNase-seq coverage plots around genomic SrfIcs in compact form
# author: Robert Gnuegge (robert.gnuegge@gmail.com)
# date: 02/08/22
# version: 1.0


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
source(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Analysis_scripts/Misc/S_cerevisiae_SrfI_cut_sites.R")
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

# set 0 to NA in GRanges object
# argument: GRanges object
# result: GRanges object
GRanges_zero_to_NA <- function(GRanges){
  out <- GRanges
  out$score[out$score == 0] <- NA
  return(out)
}


# process S1-seq data for plotting ========================================

roi <- DSB_regions(DSBs = SrfIcs, region_width = 10000)

# read S1-seq coverage
load(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Data/S1-seq/LSY4377-12B_4377-15A/LSY4377-12B_LSY4377-15A.RData")

# retain only -/+ 5 kb regions around SrfIcs (for faster data processing)
S1_seq_0 <- subsetByIntersect(subject = LSY4377_12B_0_S1_seq, query = roi)
S1_seq_1 <- subsetByIntersect(subject = LSY4377_12B_1_S1_seq, query = roi)
S1_seq_2 <- subsetByIntersect(subject = LSY4377_12B_2_S1_seq, query = roi)
S1_seq_4 <- subsetByIntersect(subject = LSY4377_12B_4_S1_seq, query = roi)

# make nt-resolved GRanges objects and sort
S1_seq_0 <- sort(as_nt_resolved_GRanges(S1_seq_0), ignore.strand = TRUE)
S1_seq_1 <- sort(as_nt_resolved_GRanges(S1_seq_1), ignore.strand = TRUE)
S1_seq_2 <- sort(as_nt_resolved_GRanges(S1_seq_2), ignore.strand = TRUE)
S1_seq_4 <- sort(as_nt_resolved_GRanges(S1_seq_4), ignore.strand = TRUE)

# apply Hanning smoother
n <- 51  # Hanning window size
S1_seq_0_smooth <- as_smoothed_GRanges(GRanges = S1_seq_0, hanning_window_size = n)
S1_seq_1_smooth <- as_smoothed_GRanges(GRanges = S1_seq_1, hanning_window_size = n)
S1_seq_2_smooth <- as_smoothed_GRanges(GRanges = S1_seq_2, hanning_window_size = n)
S1_seq_4_smooth <- as_smoothed_GRanges(GRanges = S1_seq_4, hanning_window_size = n)


# process MNase-seq data for plotting =====================================

roi <- DSB_regions(DSBs = SrfIcs, region_width = 10000)

# read MNase-seq coverage
load(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Data/MNase-seq/LSY4377-12B_4377-15A/LSY4377-12B_LSY4377-15A_trimmed.RData")

# retain only -/+ 5 kb regions around SrfIcs (for faster data processing)
MNase_seq_0 <- subsetByIntersect(subject = LSY4377_12B_0_MNase_seq_trimmed, query = roi)
MNase_seq_1 <- subsetByIntersect(subject = LSY4377_12B_1_MNase_seq_trimmed, query = roi)
MNase_seq_2 <- subsetByIntersect(subject = LSY4377_12B_2_MNase_seq_trimmed, query = roi)
MNase_seq_4 <- subsetByIntersect(subject = LSY4377_12B_4_MNase_seq_trimmed, query = roi)

# apply runmed smoother
k <- 31
MNase_seq_0$score <- runmed(x = MNase_seq_0$score, k = k)
MNase_seq_1$score <- runmed(x = MNase_seq_1$score, k = k)
MNase_seq_2$score <- runmed(x = MNase_seq_2$score, k = k)
MNase_seq_4$score <- runmed(x = MNase_seq_4$score, k = k)


# plotting ================================================================

plot_dir <- "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Figures/13_Nucleosome_occupancy_impact/LSY4377-12B_4377-15A/S1-seq_and_MNase-seq_around_SrfIcs_compact"

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

# AnnotationTrack ---------------------------------------------------------
# data downloaded from yeastmine.yeastgenome.org and curated with Adjust_chromosomal_features.R

load(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Data/S_cerevisiae_chromosomal_features/S_cerevisiae_genome_features.RData")

# adjust for use with AnnotationTrack
names(mcols(all_features))[1] <- "id"
all_features <- all_features[!(all_features$id == "") & all_features$type == "ORF"]

AT <- AnnotationTrack(range = all_features, name = NULL, showFeatureId = FALSE, cex = 0.67,  # featureAnnotation = "id"
                      arrowHeadMaxWidth = 10, fill = "white", col = "gray", fontcolor.item = "black")


# plot regions around all SrfIcs ==========================================

Y_axis_ranges <- data.frame()  # initialize for collecting
n <- 20
for (n in 1:length(SrfIcs)){

  # define name for data and figure saving
  DSB_locus <- SrfIcs[n]
  DSB_locus <- gsub(pattern = ":", replacement = "_", x = as.character(DSB_locus))
  cat("\nPlotting", DSB_locus, "...")
  
  # define plot area
  roi <- DSB_regions(DSBs = SrfIcs, region_width = 3000)[n]
  
  # get GRanges in plot area, calculate ylim, and construct DataTracks
  S1_seq_4_smooth_roi <- subsetByOverlaps(x = S1_seq_4_smooth, ranges = roi)
  ylim_smooth <- range(S1_seq_4_smooth_roi$score)
  ylim_smooth <- c(min(ylim_smooth) - 0.05 * max(pretty(ylim_smooth)), max(pretty(ylim_smooth)))  # adjust for prettier plotting
  DT_S1_seq_4_smooth <- DataTrack_fw_rev(GRanges = S1_seq_4_smooth, fw_col = fw_col, rev_col = rev_col, type = "l", name = "S1-seq", ylim = ylim_smooth)
  
  MNase_seq_0_roi <- subsetByOverlaps(x = MNase_seq_0, ranges = roi)
  MNase_seq_4_roi <- subsetByOverlaps(x = MNase_seq_4, ranges = roi)
  ylim_MNase <- range(c(MNase_seq_0_roi, MNase_seq_4_roi)$score)
  ylim_MNase <- c(min(ylim_MNase) - 0.05 * max(pretty(ylim_MNase)), max(pretty(ylim_MNase)))  # adjust for prettier plotting
  DT_MNase_seq_0 <- DataTrack(range = MNase_seq_0, type = "polygon", col = NA, fill.mountain = rep(gray(level = 0.9), 2), name = "MNase-seq", ylim = ylim_MNase)
  DT_MNase_seq_4 <- DataTrack(range = MNase_seq_4, type = "polygon", col = NA, fill.mountain = rep(gray(level = 0.8), 2), name = "MNase-seq", ylim = ylim_MNase)
  
  # add DSB to AnnotationTrack
  AT_w_DSBs <- HighlightTrack(trackList = list(AT), range = SrfIcs[n], col = JFly_colors[8], fill = NA, inBackground = FALSE)
  
  # plot
  pdf(file = "tmp.pdf", width = 3.5, height = 2)
    plotTracks(trackList = list(OverlayTrack(trackList = list(DT_MNase_seq_0, DT_MNase_seq_4, DT_S1_seq_4_smooth)), AT), 
               from = start(roi), to = end(roi), chromosome = seqnames(roi), margin = 1, innerMargin = 0, sizes = c(0.875, 0.125))
  dev.off()
  GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/", DSB_locus, ".pdf"))
  
  # record ylims (for manual scale addition)
  Y_axis_ranges <- rbind(Y_axis_ranges,
                         data.frame(DSB_locus = DSB_locus, 
                                    y_smooth_min = ylim_smooth[1], y_smooth_max = ylim_smooth[2],
                                    y_MNase_min = ylim_MNase[1], y_MNase_max = ylim_MNase[2]))
  
}

write.table(x = Y_axis_ranges, file = paste0(plot_dir, "/Y_axis_ranges.txt"), row.names = FALSE)
