# info --------------------------------------------------------------------
# purpose: identify nucleosome positions next to DSBs
# author: Robert Gnuegge (robert.gnuegge@gmail.com)
# date: 01/26/22
# version: 1.1


# preamble ----------------------------------------------------------------

# set wd to this file's location
wd_path <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(wd_path)

# libraries
library(GenomicRanges)
library(Gviz)
options(ucscChromosomeNames=FALSE)  # for using custom chromosome names (e.g. "micron")

# set wd to this file's location
wd.path <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(wd.path)

# load helper files and functions
source(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Analysis_scripts/Misc/Genomic_helper_functions.R")
source(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Analysis_scripts/Misc/S_cerevisiae_SrfI_cut_sites.R")
source(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Analysis_scripts/Misc/Misc_helper_functions.R")
source(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Analysis_scripts/Misc/JFly_colors.R")

save_dir <- "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Data/MNase-seq/LSY4377-12B_4377-15A_merged"
plot_dir <- "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Figures/13_Nucleosome_occupancy_impact/LSY4377-12B_4377-15A_merged/Nucleosome_centers"


# function definitions ====================================================

# find local maxima in noisy scatter plot
# argument: numeric vector (y coordinates)
# result: numeric vector (indices of local maximum y coordinates)
find_local_maxima <- function(y, k = NULL, smooth_spline = FALSE, peak_threshold = NULL, flank_threshold = NULL, peak_width = 147, plot = FALSE, merge_peak_threshold = NULL){
  if (plot){  # plot data
    plot(x = 1:length(y), y = y, type = "l", xlab = "x", ylab = "y")
  }
  if (!is.null(k)){  # smooth using moving average (rolling mean)
    y_2 <- moving_mean(x = y, k = k, keep = 0)
  }else{
    y_2 <- y
  }
  if (smooth_spline){  # smooth using spline
    y_3 <- smooth.spline(y_2)$y
  }else{
    y_3 <- y_2
  }
  local_max_idx <- which(diff(sign(diff(y_3))) == -2) + 1  # from https://stackoverflow.com/a/6836583
  if (!is.null(peak_threshold)){  # filter for local maxima with y values above threshold
    if (peak_threshold >= 0 & peak_threshold < 1){
      y_threshold <- max(y) * peak_threshold
    }else{
      y_threshold <- peak_threshold
    }
    stopifnot(peak_threshold > 0)
    local_max_idx <- local_max_idx[y[local_max_idx] >= y_threshold]
    if (plot){
      abline(h = y_threshold, col = "gray", lty = "dotted")  
    }
  }
  if (!is.null(flank_threshold)){  # filter for local maxima with decreasing flanking values
    left_flank_idx <- local_max_idx - round(0.5 * peak_width)
    left_flank_idx[left_flank_idx < 1] <- 1  # in case shifted out of range
    right_flank_idx <- local_max_idx + round(0.5 * peak_width)
    right_flank_idx[right_flank_idx > length(y)] <- length(y)  # in case shifted out of range
    idx_filter <- (y[left_flank_idx] <= flank_threshold * y[local_max_idx]) | (y[right_flank_idx] <= flank_threshold * y[local_max_idx])
    local_max_idx <- local_max_idx[idx_filter]
  }
  if(!is.null(merge_peak_threshold)){  # merge peaks
    idx_start <- which(diff(local_max_idx) < merge_peak_threshold)  # find index below threshold
    if(length(idx_start) > 0){
      idx_end <- idx_start + 1
      idx_all <- sort(unique(c(idx_start, idx_end)))
      idx_rm <- idx_start[which(idx_start %in% idx_end)]  # remove "middle" indices within consecutive stretches
      idx_first <- idx_start[!(idx_start %in% idx_rm)]
      idx_last <- idx_end[!(idx_end %in% idx_rm)]
      centers <- sapply(X = 1:length(idx_first), FUN = function(n) {round(mean(local_max_idx[idx_first[n]:idx_last[n]]))})
      local_max_idx <- sort(c(local_max_idx[-idx_all], centers))
    }
  }
  if (plot){
    abline(v = local_max_idx, col = "gray")
  }
  return(local_max_idx)
}



# find nucleosome centers =================================================

DSBs <- SrfIcs[-c(9, 17)]
roi <- DSB_regions(DSBs = DSBs, region_width = 6000)

# read MNase-seq coverage
load(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Data/MNase-seq/LSY4377-12B_4377-15A_merged/LSY4377-12B_LSY4377-15A_merged_trimmed.RData")

# retain only regions around SrfIcs (for faster data processing)
MNase_seq_0 <- as_nt_resolved_GRanges(GRanges = subsetByIntersect(subject = LSY4377_12B_0_merged_MNase_seq, query = roi))
MNase_seq_1 <- as_nt_resolved_GRanges(GRanges = subsetByIntersect(subject = LSY4377_12B_1_merged_MNase_seq, query = roi))
MNase_seq_2 <- as_nt_resolved_GRanges(GRanges = subsetByIntersect(subject = LSY4377_12B_2_merged_MNase_seq, query = roi))
MNase_seq_4 <- as_nt_resolved_GRanges(GRanges = subsetByIntersect(subject = LSY4377_12B_4_merged_MNase_seq, query = roi))


# 0 h ---------------------------------------------------------------------

Nuc_centers_0 <- GRanges()  # initialize
for (n in 1:length(DSBs)){
  tmp <- subsetByOverlaps(x = MNase_seq_0, ranges = roi[n])
  idx <- find_local_maxima(y = tmp$score, k = 25, smooth_spline = TRUE, peak_threshold = 0.025, flank_threshold = 0.9, merge_peak_threshold = 75)
  tmp <- tmp[idx]
  tmp$distance_to_DSB <- start(tmp) - start(DSBs[n])
  tmp$DSB_id <- as.character(DSBs[n])
  Nuc_centers_0 <- c(Nuc_centers_0, tmp)
}


# 1 h ---------------------------------------------------------------------

Nuc_centers_1 <- GRanges()  # initialize
for (n in 1:length(DSBs)){
  tmp <- subsetByOverlaps(x = MNase_seq_1, ranges = roi[n])
  idx <- find_local_maxima(y = tmp$score, k = 25, smooth_spline = TRUE, peak_threshold = 0.025, flank_threshold = 0.9, merge_peak_threshold = 75)
  tmp <- tmp[idx]
  tmp$distance_to_DSB <- start(tmp) - start(DSBs[n])
  tmp$DSB_id <- as.character(DSBs[n])
  Nuc_centers_1 <- c(Nuc_centers_1, tmp)
}


# 2 h ---------------------------------------------------------------------

Nuc_centers_2 <- GRanges()  # initialize
for (n in 1:length(DSBs)){
  tmp <- subsetByOverlaps(x = MNase_seq_2, ranges = roi[n])
  idx <- find_local_maxima(y = tmp$score, k = 25, smooth_spline = TRUE, peak_threshold = 0.025, flank_threshold = 0.9, merge_peak_threshold = 75)
  tmp <- tmp[idx]
  tmp$distance_to_DSB <- start(tmp) - start(DSBs[n])
  tmp$DSB_id <- as.character(DSBs[n])
  Nuc_centers_2 <- c(Nuc_centers_2, tmp)
}


# 4 h ---------------------------------------------------------------------

Nuc_centers_4 <- GRanges()  # initialize
for (n in 1:length(DSBs)){
  tmp <- subsetByOverlaps(x = MNase_seq_4, ranges = roi[n])
  idx <- find_local_maxima(y = tmp$score, k = 25, smooth_spline = TRUE, peak_threshold = 0.025, flank_threshold = 0.9, merge_peak_threshold = 75)
  tmp <- tmp[idx]
  tmp$distance_to_DSB <- start(tmp) - start(DSBs[n])
  tmp$DSB_id <- as.character(DSBs[n])
  Nuc_centers_4 <- c(Nuc_centers_4, tmp)
}


# save data
save(list = paste0("Nuc_centers_", c(0, 1, 2, 4)), file = paste0(save_dir, "/Nuc_centers_around_SrfIcs.RData"))


# plot identified nucleosome centers ======================================

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
MNase_seq_col <- gray(level = 0.9)
MNase_seq_border_col <- gray(level = 0.75)
Nuc_center_col = JFly_colors[8]


# AnnotationTrack ---------------------------------------------------------
# data downloaded from yeastmine.yeastgenome.org and curated with Adjust_chromosomal_features.R

load(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Data/S_cerevisiae_chromosomal_features/S_cerevisiae_genome_features.RData")

# adjust for use with AnnotationTrack
names(mcols(all_features))[1] <- "id"
all_features <- all_features[!(all_features$id == "") & all_features$type == "ORF"]

AT <- AnnotationTrack(range = all_features, name = NULL, featureAnnotation = "id", cex = 0.5,
                      arrowHeadMaxWidth = 10, fill = "white", col = "gray", fontcolor.item = "black")


# plot regions around DSBs ------------------------------------------------

for (n in 1:length(DSBs)){
  
  # define name for data and figure saving
  DSB_locus <- DSBs[n]
  DSB_locus <- gsub(pattern = ":", replacement = "_", x = as.character(DSB_locus))
  cat("\nPlotting", DSB_locus, "...")
  
  # define plot area
  roi <- DSB_regions(DSBs = DSBs, region_width = 4000)[n]
  
  # get GRanges in plot area, calculate ylim, and construct DataTracks
  MNase_seq_0_roi <- subsetByOverlaps(x = MNase_seq_0, ranges = roi)
  MNase_seq_1_roi <- subsetByOverlaps(x = MNase_seq_1, ranges = roi)
  MNase_seq_2_roi <- subsetByOverlaps(x = MNase_seq_2, ranges = roi)
  MNase_seq_4_roi <- subsetByOverlaps(x = MNase_seq_4, ranges = roi)
  
  ylim_MNase <- range(c(MNase_seq_0_roi, MNase_seq_1_roi, MNase_seq_2_roi, MNase_seq_4_roi)$score)
  ylim_MNase <- c(min(ylim_MNase) - 0.05 * max(pretty(ylim_MNase)), max(pretty(ylim_MNase)))  # adjust for prettier plotting
  
  DT_MNase_seq_0 <- DataTrack(range = MNase_seq_0, type = "polygon", col = MNase_seq_border_col, fill.mountain = rep(MNase_seq_col, 2), name = "0 h", ylim = ylim_MNase)
  DT_MNase_seq_1 <- DataTrack(range = MNase_seq_1, type = "polygon", col = MNase_seq_border_col, fill.mountain = rep(MNase_seq_col, 2), name = "1 h", ylim = ylim_MNase)
  DT_MNase_seq_2 <- DataTrack(range = MNase_seq_2, type = "polygon", col = MNase_seq_border_col, fill.mountain = rep(MNase_seq_col, 2), name = "2 h", ylim = ylim_MNase)
  DT_MNase_seq_4 <- DataTrack(range = MNase_seq_4, type = "polygon", col = MNase_seq_border_col, fill.mountain = rep(MNase_seq_col, 2), name = "4 h", ylim = ylim_MNase)
  
  # add DSB to AnnotationTrack
  AT_w_DSBs <- HighlightTrack(trackList = list(AT), range = DSBs[n], col = JFly_colors[8], fill = NA, inBackground = FALSE)
  
  # plot
  pdf(file = "tmp.pdf", width = 3.5, height = 4)
  plotTracks(trackList = list(HighlightTrack(trackList = list(DT_MNase_seq_0), 
                                             range = subsetByOverlaps(x = Nuc_centers_0, ranges = roi), 
                                             col = Nuc_center_col, fill = NA, inBackground = FALSE),
                              HighlightTrack(trackList = list(DT_MNase_seq_1), 
                                             range = subsetByOverlaps(x = Nuc_centers_1, ranges = roi), 
                                             col = Nuc_center_col, fill = NA, inBackground = FALSE),
                              HighlightTrack(trackList = list(DT_MNase_seq_2), 
                                             range = subsetByOverlaps(x = Nuc_centers_2, ranges = roi), 
                                             col = Nuc_center_col, fill = NA, inBackground = FALSE),
                              HighlightTrack(trackList = list(DT_MNase_seq_4), 
                                             range = subsetByOverlaps(x = Nuc_centers_4, ranges = roi), 
                                             col = Nuc_center_col, fill = NA, inBackground = FALSE),
                              AT_w_DSBs), 
             from = start(roi), to = end(roi), chromosome = seqnames(roi), margin = 1, innerMargin = 0)
  dev.off()
  GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/", DSB_locus, ".pdf"))
  
}