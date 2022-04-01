# info --------------------------------------------------------------------
# purpose: generate S1-seq coverage plots around genomic DSBs
# author: Robert Gnuegge (robert.gnuegge@gmail.com)
# date: 01/31/22
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

load(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Data/S1-seq/LSY4810-41D_4811-11B_4812-5A/LSY4810-41D_4811-11B_4812-5A.RData")

# function definitions ====================================================

# make Gviz track with separately marked fw and rev scores ----------------
# argument: GRanges object, color definitions, ... (e.g. type, name)
# result: OverlayTrack
DataTrack_fw_rev <- function(GRanges, fw_col, rev_col, ...){
  if(is.null(ylim)){
    ylim <- pretty(range(GRanges$score))
  }
  fw <- DataTrack(range = GRanges[strand(GRanges) == "+"], col = fw_col, ...)
  rev <- DataTrack(range = GRanges[strand(GRanges) == "-"], col = rev_col, ...)
  OverlayTrack(trackList = list(fw, rev))
}

# set 0 to NA in GRanges object -------------------------------------------
# argument: GRanges object
# result: GRanges object
GRanges_zero_to_NA <- function(GRanges){
  out <- GRanges
  out$score[out$score == 0] <- NA
  return(out)
}

# make Gviz track with raw and smooth fw and rev scores -----------------
# argument: GRanges objects with raw and smoothed data, ylims for raw and smoothed data
# result: OverlayTrack
DataTrack_raw_smooth_fw_rev <- function(GRanges_raw, GRanges_smooth, fw_col, rev_col, ylim_raw, ylim_smooth, alpha = 0.5, lwd = 1, name, zero_to_NA = TRUE){
  if(zero_to_NA){
    GRanges_raw <- GRanges_zero_to_NA(GRanges_raw)
    GRanges_smooth <- GRanges_zero_to_NA(GRanges_smooth)
  }
  fw_raw <- DataTrack(range = GRanges_raw[strand(GRanges_raw) == "+"], col = adjustcolor(col = fw_col, alpha.f = alpha), ylim = ylim_raw, type = "h", name = name)
  rev_raw <- DataTrack(range = GRanges_raw[strand(GRanges_raw) == "-"], col = adjustcolor(col = rev_col, alpha.f = alpha), ylim = ylim_raw, type = "h")
  fw_smooth <- DataTrack(range = GRanges_smooth[strand(GRanges_smooth) == "+"], col = fw_col, ylim = ylim_smooth, type = "l", lwd = lwd, name = name)
  rev_smooth <- DataTrack(range = GRanges_smooth[strand(GRanges_smooth) == "-"], col = rev_col, ylim = ylim_smooth, type = "l", lwd = lwd)
  OverlayTrack(trackList = list(fw_raw, rev_raw, fw_smooth, rev_smooth))
}

# make Gviz track with separately marked fw and rev scores
# argument: GRanges object, color definitions, ... (e.g. type, name)
# result: OverlayTrack
DataTrack_fw_rev <- function(GRanges, fw_col, rev_col, ...){
  if(is.null(ylim)){
    ylim <- pretty(range(GRanges$score))
  }
  fw <- DataTrack(range = GRanges[strand(GRanges) == "+"], col = fw_col, ...)
  rev <- DataTrack(range = GRanges[strand(GRanges) == "-"], col = rev_col, ...)
  OverlayTrack(trackList = list(fw, rev))
}


# plotting setup ==========================================================

plot_dir <- "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Figures/11_Transcription_impact/LSY4810-41D_4811-11B_4812-5A"

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

Y_axis_ranges <- data.frame()  # initialize for collecting


# AnnotationTrack ---------------------------------------------------------
# data downloaded from yeastmine.yeastgenome.org and curated with Adjust_chromosomal_features.R

load(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Data/S_cerevisiae_chromosomal_features/S_cerevisiae_genome_features.RData")

# adjust for use with AnnotationTrack
names(mcols(all_features))[1] <- "id"
all_features <- all_features[!(all_features$id == "") & all_features$type %in% c("ORF", "tRNA gene")]

# remove MET17, MRM1, HIS3, and YOR203W (shifted up by plasmid integration)
all_features <- all_features[!all_features$id %in% c("MET17", "MRM1", "HIS3", "YOR203W")]

# add transcription reporters
tmp <- GRanges(seqnames = c("chrXII", "chrXV"), ranges = IRanges(start = c(733579, 721161), end = c(734277, 722081)), strand = c("+", "+"))
mcols(tmp) <- data.frame(id = c("mKate2", "Citrine"), systematic_name = NA, type = c("ORF", "ORF"), qualifier = "Verified")
all_features <- sort(c(all_features, tmp), ignore.strand = TRUE)

# make AnnotationTrack
AT <- AnnotationTrack(range = all_features, name = NULL, featureAnnotation = "id", cex = 0.5,
                      arrowHeadMaxWidth = 10, fill = "white", col = "gray", fontcolor.item = "black")


# mKate2 reporter =========================================================
DSB <- GRanges("chrXII:734489")

# define name for data and figure saving
DSB_locus <- DSB
DSB_locus <- gsub(pattern = ":", replacement = "_", x = as.character(DSB_locus))
cat("\nPlotting", DSB_locus, "...")

# load and process data ---------------------------------------------------
# retain only -/+ 5 kb regions around DSBs (for faster data processing)
LSY4810_41D_2 <- subsetByIntersect(subject = LSY4810_41D_2_S1_seq, query = DSB_regions(DSBs = DSB, region_width = 10000))
LSY4811_11B_2  <- subsetByIntersect(subject = LSY4811_11B_2_S1_seq, query = DSB_regions(DSBs = DSB, region_width = 10000))
LSY4812_5A_2  <- subsetByIntersect(subject = LSY4812_5A_2_S1_seq, query = DSB_regions(DSBs = DSB, region_width = 10000))

# make nt-resolved GRanges objects and sort
LSY4810_41D_2 <- sort(as_nt_resolved_GRanges(LSY4810_41D_2), ignore.strand = TRUE)
LSY4811_11B_2 <- sort(as_nt_resolved_GRanges(LSY4811_11B_2), ignore.strand = TRUE)
LSY4812_5A_2 <- sort(as_nt_resolved_GRanges(LSY4812_5A_2), ignore.strand = TRUE)
LSY4810_41D_2$score[LSY4810_41D_2$score > 1500] <- 1500  # adjust for better visualization
LSY4811_11B_2$score[LSY4811_11B_2$score > 1500] <- 1500
LSY4812_5A_2$score[LSY4812_5A_2$score > 1500] <- 1500

# apply Hanning smoother
n <- 51  # Hanning window size
LSY4810_41D_2_smooth <- as_smoothed_GRanges(GRanges = LSY4810_41D_2, hanning_window_size = n)
LSY4811_11B_2_smooth <- as_smoothed_GRanges(GRanges = LSY4811_11B_2, hanning_window_size = n)
LSY4812_5A_2_smooth <- as_smoothed_GRanges(GRanges = LSY4812_5A_2, hanning_window_size = n)

# define plot area
roi <- DSB_regions(DSBs = DSB, region_width = 2000)

# get GRanges in plot area (for ylim calculation)
LSY4810_41D_2_roi <- subsetByOverlaps(x = LSY4810_41D_2, ranges = roi)
LSY4811_11B_2_roi <- subsetByOverlaps(x = LSY4811_11B_2, ranges = roi)
LSY4812_5A_2_roi <- subsetByOverlaps(x = LSY4812_5A_2, ranges = roi)

LSY4810_41D_2_smooth_roi <- subsetByOverlaps(x = LSY4810_41D_2_smooth, ranges = roi)
LSY4811_11B_2_smooth_roi <- subsetByOverlaps(x = LSY4811_11B_2_smooth, ranges = roi)
LSY4812_5A_2_smooth_roi <- subsetByOverlaps(x = LSY4812_5A_2_smooth, ranges = roi)

# find common ylim
ylim <- range(c(LSY4811_11B_2_roi, LSY4810_41D_2_roi, LSY4812_5A_2_roi)$score)
ylim <- c(min(ylim) - 0.05 * max(pretty(ylim)), max(pretty(ylim)))  # adjust for prettier plotting

ylim_smooth <- range(c(LSY4811_11B_2_smooth_roi, LSY4810_41D_2_smooth_roi, LSY4812_5A_2_smooth_roi)$score)
ylim_smooth <- c(min(ylim_smooth) - 0.05 * max(pretty(ylim_smooth)), max(pretty(ylim_smooth)))  # adjust for prettier plotting

Y_axis_ranges <- rbind(Y_axis_ranges,
                       data.frame(DSB_locus = DSB_locus, 
                                  y_raw_min = ylim[1], y_raw_max = ylim[2],
                                  y_smooth_min = ylim_smooth[1], y_smooth_max = ylim_smooth[2]))

# make OverlayTracks
DT_LSY4812_5A_2 <- DataTrack_raw_smooth_fw_rev(GRanges_raw = LSY4812_5A_2_roi, GRanges_smooth = LSY4812_5A_2_smooth_roi,
                                               fw_col = fw_col, rev_col = rev_col, ylim_raw = ylim, ylim_smooth = ylim_smooth, name = "no P")
DT_LSY4811_11B_2 <- DataTrack_raw_smooth_fw_rev(GRanges_raw = LSY4811_11B_2_roi, GRanges_smooth = LSY4811_11B_2_smooth_roi,
                                                fw_col = fw_col, rev_col = rev_col, ylim_raw = ylim, ylim_smooth = ylim_smooth, name = "P_ACT1")
DT_LSY4810_41D_2 <- DataTrack_raw_smooth_fw_rev(GRanges_raw = LSY4810_41D_2_roi, GRanges_smooth = LSY4810_41D_2_smooth_roi,
                                               fw_col = fw_col, rev_col = rev_col, ylim_raw = ylim, ylim_smooth = ylim_smooth, name = "P_TDH3")

# add DSB to AnnotationTrack
AT_w_DSBs <- HighlightTrack(trackList = list(AT), range = DSB, col = JFly_colors[8], fill = NA, inBackground = FALSE)

# plot
pdf(file = "tmp.pdf", width = 3.5, height = 2.5)
plotTracks(trackList = list(DT_LSY4812_5A_2, DT_LSY4811_11B_2, DT_LSY4810_41D_2, AT_w_DSBs),  # 
           chromosome = seqnames(roi), from = start(roi), to = end(roi), margin = 1, innerMargin = 0)
dev.off()
GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/", DSB_locus, "_mKate2.pdf"))


# mKate2 overlay ----------------------------------------------------------

# construct DataTracks for smoothed data
DT_LSY4812_5A_2 <- DataTrack_fw_rev(GRanges = GRanges_zero_to_NA(LSY4812_5A_2_smooth_roi), fw_col = adjustcolor(col = fw_col, alpha.f = 0.5), rev_col = adjustcolor(col = rev_col, alpha.f = 0.33), type = "l", name = "No P", ylim = ylim_smooth)
DT_LSY4810_41D_2 <- DataTrack_fw_rev(GRanges = GRanges_zero_to_NA(LSY4810_41D_2_smooth_roi), fw_col = adjustcolor(col = fw_col, alpha.f = 1), rev_col = adjustcolor(col = rev_col, alpha.f = 1), type = "l", name = "P_TDH3", ylim = ylim_smooth)

# plot
pdf(file = "tmp.pdf", width = 3.5, height = 1.5)
plotTracks(trackList = list(OverlayTrack(trackList = list(DT_LSY4812_5A_2, DT_LSY4810_41D_2), name = "S1-seq"), AT_w_DSBs),
           chromosome = seqnames(roi), from = start(roi), to = end(roi), margin = 1, innerMargin = 0, sizes = c(0.85, 0.15))
dev.off()
GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/", DSB_locus, "_mKate2_overlay.pdf"))


# Citrine reporter ========================================================
DSB <- GRanges("chrXV:722028")

# define name for data and figure saving
DSB_locus <- DSB
DSB_locus <- gsub(pattern = ":", replacement = "_", x = as.character(DSB_locus))
cat("\nPlotting", DSB_locus, "...")

# load and process data ---------------------------------------------------
# retain only -/+ 5 kb regions around DSBs (for faster data processing)
LSY4810_41D_2 <- subsetByIntersect(subject = LSY4810_41D_2_S1_seq, query = DSB_regions(DSBs = DSB, region_width = 10000))
LSY4811_11B_2  <- subsetByIntersect(subject = LSY4811_11B_2_S1_seq, query = DSB_regions(DSBs = DSB, region_width = 10000))
LSY4812_5A_2  <- subsetByIntersect(subject = LSY4812_5A_2_S1_seq, query = DSB_regions(DSBs = DSB, region_width = 10000))

# make nt-resolved GRanges objects and sort
LSY4810_41D_2 <- sort(as_nt_resolved_GRanges(LSY4810_41D_2), ignore.strand = TRUE)
LSY4811_11B_2 <- sort(as_nt_resolved_GRanges(LSY4811_11B_2), ignore.strand = TRUE)
LSY4812_5A_2 <- sort(as_nt_resolved_GRanges(LSY4812_5A_2), ignore.strand = TRUE)
LSY4810_41D_2$score[LSY4810_41D_2$score > 750] <- 750  # adjust for better visualization
LSY4811_11B_2$score[LSY4811_11B_2$score > 750] <- 750
LSY4812_5A_2$score[LSY4812_5A_2$score > 750] <- 750

# apply Hanning smoother
n <- 51  # Hanning window size
LSY4810_41D_2_smooth <- as_smoothed_GRanges(GRanges = LSY4810_41D_2, hanning_window_size = n)
LSY4811_11B_2_smooth <- as_smoothed_GRanges(GRanges = LSY4811_11B_2, hanning_window_size = n)
LSY4812_5A_2_smooth <- as_smoothed_GRanges(GRanges = LSY4812_5A_2, hanning_window_size = n)

# define plot area
roi <- DSB_regions(DSBs = DSB, region_width = 2000)

# get GRanges in plot area (for ylim calculation)
LSY4810_41D_2_roi <- subsetByOverlaps(x = LSY4810_41D_2, ranges = roi)
LSY4811_11B_2_roi <- subsetByOverlaps(x = LSY4811_11B_2, ranges = roi)
LSY4812_5A_2_roi <- subsetByOverlaps(x = LSY4812_5A_2, ranges = roi)

LSY4810_41D_2_smooth_roi <- subsetByOverlaps(x = LSY4810_41D_2_smooth, ranges = roi)
LSY4811_11B_2_smooth_roi <- subsetByOverlaps(x = LSY4811_11B_2_smooth, ranges = roi)
LSY4812_5A_2_smooth_roi <- subsetByOverlaps(x = LSY4812_5A_2_smooth, ranges = roi)

# find common ylim
ylim <- range(c(LSY4811_11B_2_roi, LSY4810_41D_2_roi, LSY4812_5A_2_roi)$score)
ylim <- c(min(ylim) - 0.05 * max(pretty(ylim)), max(pretty(ylim)))  # adjust for prettier plotting

ylim_smooth <- range(c(LSY4811_11B_2_smooth_roi, LSY4810_41D_2_smooth_roi, LSY4812_5A_2_smooth_roi)$score)
ylim_smooth <- c(min(ylim_smooth) - 0.05 * max(pretty(ylim_smooth)), max(pretty(ylim_smooth)))  # adjust for prettier plotting

Y_axis_ranges <- rbind(Y_axis_ranges,
                       data.frame(DSB_locus = DSB_locus, 
                                  y_raw_min = ylim[1], y_raw_max = ylim[2],
                                  y_smooth_min = ylim_smooth[1], y_smooth_max = ylim_smooth[2]))

# make OverlayTracks
DT_LSY4810_41D_2 <- DataTrack_raw_smooth_fw_rev(GRanges_raw = LSY4810_41D_2_roi, GRanges_smooth = LSY4810_41D_2_smooth_roi,
                                                fw_col = fw_col, rev_col = rev_col, ylim_raw = ylim, ylim_smooth = ylim_smooth, name = "no P")
DT_LSY4812_5A_2 <- DataTrack_raw_smooth_fw_rev(GRanges_raw = LSY4812_5A_2_roi, GRanges_smooth = LSY4812_5A_2_smooth_roi,
                                               fw_col = fw_col, rev_col = rev_col, ylim_raw = ylim, ylim_smooth = ylim_smooth, name = "P_ACT1")
DT_LSY4811_11B_2 <- DataTrack_raw_smooth_fw_rev(GRanges_raw = LSY4811_11B_2_roi, GRanges_smooth = LSY4811_11B_2_smooth_roi,
                                                fw_col = fw_col, rev_col = rev_col, ylim_raw = ylim, ylim_smooth = ylim_smooth, name = "P_TDH3")

# add DSB to AnnotationTrack
AT_w_DSBs <- HighlightTrack(trackList = list(AT), range = DSB, col = JFly_colors[8], fill = NA, inBackground = FALSE)

# plot
pdf(file = "tmp.pdf", width = 3.5, height = 2.5)
plotTracks(trackList = list(DT_LSY4810_41D_2, DT_LSY4812_5A_2, DT_LSY4811_11B_2, AT_w_DSBs),  # 
           chromosome = seqnames(roi), from = start(roi), to = end(roi), margin = 1, innerMargin = 0)
dev.off()
GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/", DSB_locus, "_Citrine.pdf"))


# Citrine overlay ----------------------------------------------------------

# construct DataTracks for smoothed data
DT_LSY4810_41D_2 <- DataTrack_fw_rev(GRanges = GRanges_zero_to_NA(LSY4810_41D_2_smooth_roi), fw_col = adjustcolor(col = fw_col, alpha.f = 0.5), rev_col = adjustcolor(col = rev_col, alpha.f = 1), type = "l", name = "No P", ylim = ylim_smooth)
DT_LSY4811_11B_2 <- DataTrack_fw_rev(GRanges = GRanges_zero_to_NA(LSY4811_11B_2_smooth_roi), fw_col = adjustcolor(col = fw_col, alpha.f = 1), rev_col = adjustcolor(col = rev_col, alpha.f = 0.33), type = "l", name = "P_TDH3", ylim = ylim_smooth)

# plot
pdf(file = "tmp.pdf", width = 3.5, height = 1.5)
plotTracks(trackList = list(OverlayTrack(trackList = list(DT_LSY4810_41D_2, DT_LSY4811_11B_2), name = "S1-seq"), AT_w_DSBs),
           chromosome = seqnames(roi), from = start(roi), to = end(roi), margin = 1, innerMargin = 0, sizes = c(0.85, 0.15))
dev.off()
GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/", DSB_locus, "_Citrine_overlay.pdf"))


write.table(x = Y_axis_ranges, file = paste0(plot_dir, "/Y_axis_ranges.txt"), row.names = FALSE)