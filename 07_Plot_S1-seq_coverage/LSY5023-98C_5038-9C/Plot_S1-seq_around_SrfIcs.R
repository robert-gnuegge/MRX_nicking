# info --------------------------------------------------------------------
# purpose: generate S1-seq coverage plots around genomic SrfIcs
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
source(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Analysis_scripts/Misc/S_cerevisiae_SrfI_cut_sites.R")
source(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Analysis_scripts/Misc/Misc_helper_functions.R")
source(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Analysis_scripts/Misc/JFly_colors.R")


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
DataTrack_raw_smooth_fw_rev <- function(GRanges_raw, GRanges_smooth, fw_col, rev_col, ylim_raw, ylim_smooth, alpha = 0.67, lwd = 1, name, zero_to_NA = TRUE){
  if(zero_to_NA){
    GRanges_raw <- GRanges_zero_to_NA(GRanges_raw)
    GRanges_smooth <- GRanges_zero_to_NA(GRanges_smooth)
  }
  fw_raw <- DataTrack(range = GRanges_raw[strand(GRanges_raw) == "+"], col = fw_col, ylim = ylim_raw, type = "h", name = name)
  rev_raw <- DataTrack(range = GRanges_raw[strand(GRanges_raw) == "-"], col = rev_col, ylim = ylim_raw, type = "h")
  fw_smooth <- DataTrack(range = GRanges_smooth[strand(GRanges_smooth) == "+"], col = adjustcolor(col = fw_col, alpha.f = alpha), ylim = ylim_smooth, type = "l", lwd = lwd, name = name)
  rev_smooth <- DataTrack(range = GRanges_smooth[strand(GRanges_smooth) == "-"], col = adjustcolor(col = rev_col, alpha.f = alpha), ylim = ylim_smooth, type = "l", lwd = lwd)
  OverlayTrack(trackList = list(fw_raw, rev_raw, fw_smooth, rev_smooth))
}


# process data for plotting ===============================================

# read S1-seq coverage
load(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Data/S1-seq/LSY5023-98C_5038-9C/LSY5023-98C_5038-9C.RData")
load(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Data/S1-seq/LSY4377-12B_4377-15A_merged/LSY4377-12B_LSY4377-15A_merged.RData")

# retain only -/+ 5 kb regions around SrfIcs (for faster data processing)
LSY5023_98C_2_S1_seq <- subsetByIntersect(subject = LSY5023_98C_2_S1_seq, query = DSB_regions(DSBs = SrfIcs, region_width = 10000))
LSY5038_9C_2_S1_seq  <- subsetByIntersect(subject = LSY5038_9C_2_S1_seq, query = DSB_regions(DSBs = SrfIcs, region_width = 10000))
LSY4377_12B_2_S1_seq <- subsetByIntersect(subject = LSY4377_12B_2_merged_S1_seq, query = DSB_regions(DSBs = SrfIcs, region_width = 10000))

# make nt-resolved GRanges objects and sort
LSY5023_98C_2 <- sort(as_nt_resolved_GRanges(LSY5023_98C_2_S1_seq), ignore.strand = TRUE)
LSY5038_9C_2 <- sort(as_nt_resolved_GRanges(LSY5038_9C_2_S1_seq), ignore.strand = TRUE)
LSY4377_12B_2 <- sort(as_nt_resolved_GRanges(LSY4377_12B_2_S1_seq), ignore.strand = TRUE)

# apply Hanning smoother
n <- 51  # Hanning window size
LSY5023_98C_2_smooth <- as_smoothed_GRanges(GRanges = LSY5023_98C_2, hanning_window_size = n)
LSY5038_9C_2_smooth <- as_smoothed_GRanges(GRanges = LSY5038_9C_2, hanning_window_size = n)
LSY4377_12B_2_smooth <- as_smoothed_GRanges(GRanges = LSY4377_12B_2, hanning_window_size = n)


# plotting ================================================================

plot_dir <- "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Figures/04_S1-seq_next_to_DSBs/LSY5023-98C_5038-9C"

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

for (n in 1:length(SrfIcs)){

  # define name for data and figure saving
  DSB_locus <- SrfIcs[n]
  DSB_locus <- gsub(pattern = ":", replacement = "_", x = as.character(DSB_locus))
  cat("\nPlotting", DSB_locus, "...")
  
  # define plot area
  roi <- DSB_regions(DSBs = SrfIcs, region_width = 3000)[n]
  
  # get GRanges in plot area (for ylim calculation)
  LSY5023_98C_2_roi <- subsetByOverlaps(x = LSY5023_98C_2, ranges = roi)
  LSY5038_9C_2_roi <- subsetByOverlaps(x = LSY5038_9C_2, ranges = roi)
  LSY4377_12B_2_roi <- subsetByOverlaps(x = LSY4377_12B_2, ranges = roi)
  
  LSY5023_98C_2_smooth_roi <- subsetByOverlaps(x = LSY5023_98C_2_smooth, ranges = roi)
  LSY5038_9C_2_smooth_roi <- subsetByOverlaps(x = LSY5038_9C_2_smooth, ranges = roi)
  LSY4377_12B_2_smooth_roi <- subsetByOverlaps(x = LSY4377_12B_2_smooth, ranges = roi)
  
  # find common ylim
  ylim <- range(c(LSY4377_12B_2_roi, LSY5038_9C_2_roi, LSY5023_98C_2_roi)$score)
  ylim <- c(min(ylim) - 0.05 * max(pretty(ylim)), max(pretty(ylim)))  # adjust for prettier plotting
  
  ylim_smooth <- range(c(LSY4377_12B_2_smooth_roi, LSY5038_9C_2_smooth_roi, LSY5023_98C_2_smooth_roi)$score)
  ylim_smooth <- c(min(ylim_smooth) - 0.05 * max(pretty(ylim_smooth)), max(pretty(ylim_smooth)))  # adjust for prettier plotting
  
  Y_axis_ranges <- rbind(Y_axis_ranges,
                         data.frame(DSB_locus = DSB_locus, 
                                    y_raw_min = ylim[1], y_raw_max = ylim[2],
                                    y_smooth_min = ylim_smooth[1], y_smooth_max = ylim_smooth[2]))
  
  
  # make OverlayTracks
  DT_LSY4377_12B_2 <- DataTrack_raw_smooth_fw_rev(GRanges_raw = LSY4377_12B_2_roi, GRanges_smooth = LSY4377_12B_2_smooth_roi,
                                                  fw_col = fw_col, rev_col = rev_col, ylim_raw = ylim, ylim_smooth = ylim_smooth, name = "WT\n2 h")
  DT_LSY5038_9C_2 <- DataTrack_raw_smooth_fw_rev(GRanges_raw = LSY5038_9C_2_roi, GRanges_smooth = LSY5038_9C_2_smooth_roi,
                                                  fw_col = fw_col, rev_col = rev_col, ylim_raw = ylim, ylim_smooth = ylim_smooth, name = "pho4Delta\n2 h")
  DT_LSY5023_98C_2 <- DataTrack_raw_smooth_fw_rev(GRanges_raw = LSY5023_98C_2_roi, GRanges_smooth = LSY5023_98C_2_smooth_roi,
                                                 fw_col = fw_col, rev_col = rev_col, ylim_raw = ylim, ylim_smooth = ylim_smooth, name = "pho4-SA\n2 h")
  
  # plot
  pdf(file = "tmp.pdf", width = 3.5, height = 3)
    plotTracks(trackList = list(DT_LSY4377_12B_2, DT_LSY5038_9C_2, DT_LSY5023_98C_2),
               chromosome = seqnames(roi), from = start(roi), to = end(roi), margin = 1, innerMargin = 0)
  dev.off()
  GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/", DSB_locus, ".pdf"))

}

write.table(x = Y_axis_ranges, file = paste0(plot_dir, "/Y_axis_ranges.txt"), row.names = FALSE)