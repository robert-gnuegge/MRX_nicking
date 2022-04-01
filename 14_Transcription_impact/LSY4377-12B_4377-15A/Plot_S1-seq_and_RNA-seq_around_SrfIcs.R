# info --------------------------------------------------------------------
# purpose: generate S1-seq and RNA-seq coverage plots around genomic SrfIcs
# author: Robert Gnuegge (robert.gnuegge@gmail.com)
# date: 01/19/22
# version: 2.0


# preamble ----------------------------------------------------------------

# set wd to this file's location
wd_path <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(wd_path)

# libraries
library(GenomicRanges)
library(rtracklayer)
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
  if(is.null(ylim)){
    ylim <- pretty(range(GRanges$score))
  }
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

# make Gviz track with raw and smooth fw and rev scores
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


# process S1-seq data for plotting ========================================

roi <- DSB_regions(DSBs = SrfIcs, region_width = 10000)

# read S1-seq coverage
load(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Data/S1-seq/LSY4377-12B_4377-15A/LSY4377-12B_LSY4377-15A.RData")

# retain only -/+ 5 kb regions around SrfIcs (for faster data processing)
LSY4377_12B_2_S1_seq <- subsetByIntersect(subject = LSY4377_12B_2_S1_seq, query = roi)

# make nt-resolved GRanges objects and sort
S1_seq <- sort(as_nt_resolved_GRanges(LSY4377_12B_2_S1_seq), ignore.strand = TRUE)

# apply Hanning smoother
n <- 51  # Hanning window size
S1_seq_smooth <- as_smoothed_GRanges(GRanges = S1_seq, hanning_window_size = n)


# plotting ================================================================

plot_dir <- "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Figures/11_Transcription_impact/LSY4377-12B_4377-15A/S1-seq_and_RNA-seq_around_SrfIcs"

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
RNA_seq_col <- gray(level = 0.9)

# RNA-seq DataTrack -------------------------------------------------------
# data from Maya-Miles et al, 2019 (pmid: 31331360)

# there are two biological replicates
RNA_seq_1 <- import("/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Resources/Maya-Miles2019/GSM3567364_w303_rep1.bigwig",
                    which = roi)
RNA_seq_2 <- import("/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Resources/Maya-Miles2019/GSM3567365_w303_rep2.bigwig",
                    which = roi)

# let's check how they correlate
RNA_seq_1 <- subsetByIntersect(subject = RNA_seq_1, query = RNA_seq_2)
RNA_seq_2 <- subsetByIntersect(subject = RNA_seq_2, query = RNA_seq_1)
all(granges(RNA_seq_1) == granges(RNA_seq_2))
plot(x = RNA_seq_1$score, y = RNA_seq_2$score, pch = 20, col = gray(level = 0, alpha = 0.5))
cor(x = RNA_seq_1$score, y = RNA_seq_2$score, method = "pearson")
cor.test(x = RNA_seq_1$score, y = RNA_seq_2$score, method = "pearson", alternative = "greater")

# let's take the average RNA-seq score
RNA_seq <- RNA_seq_1
RNA_seq$score <- apply(X = cbind(RNA_seq_1$score, RNA_seq_2$score), MARGIN = 1, FUN = mean)


# AnnotationTrack ---------------------------------------------------------
# data downloaded from yeastmine.yeastgenome.org and curated with Adjust_chromosomal_features.R

load(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Resources/S_cerevisiae_chromosomal_features/S_cerevisiae_genome_features.RData")

# adjust for use with AnnotationTrack
names(mcols(all_features))[1] <- "id"
all_features <- all_features[!(all_features$id == "") & all_features$type == "ORF"]

AT <- AnnotationTrack(range = all_features, name = NULL, featureAnnotation = "id", cex = 0.5,
                      arrowHeadMaxWidth = 10, fill = "white", col = "gray", fontcolor.item = "black")


# plot regions around all SrfIcs ==========================================

Y_axis_ranges <- data.frame()  # initialize for collecting

for (n in 1:length(SrfIcs)){

  # define name for data and figure saving
  DSB_locus <- SrfIcs[n]
  DSB_locus <- gsub(pattern = ":", replacement = "_", x = as.character(DSB_locus))
  cat("\nPlotting", DSB_locus, "...")
  
  # define plot area
  roi <- DSB_regions(DSBs = SrfIcs, region_width = 3000)[n]
  
  # get GRanges in plot area, calculate ylim, and construct DataTracks
  S1_seq_roi <- subsetByOverlaps(x = S1_seq, ranges = roi)
  ylim <- range(pretty(range(S1_seq_roi$score)))
  DT_S1_seq <- DataTrack_fw_rev(GRanges = GRanges_zero_to_NA(S1_seq), fw_col = fw_col, rev_col = rev_col, type = "h", name = "S1-seq", ylim = ylim)
  
  S1_seq_smooth_roi <- subsetByOverlaps(x = S1_seq_smooth, ranges = roi)
  ylim_smooth <- range(pretty(range(S1_seq_smooth_roi$score)))
  DT_S1_seq_smooth <- DataTrack_fw_rev(GRanges = GRanges_zero_to_NA(S1_seq_smooth), type = "l", ylim = ylim_smooth, name = "smoothed S1-seq",
                                       fw_col = adjustcolor(col = fw_col, alpha.f = 0.67), rev_col = adjustcolor(col = rev_col, alpha.f = 0.67))
  
  RNA_seq_roi <- subsetByOverlaps(x = RNA_seq, ranges = roi)
  ylim_RNA <- range(pretty(range(RNA_seq_roi$score)))
  DT_RNA_seq <- DataTrack(range = RNA_seq, type = "polygon", ylim = ylim_RNA, name = "RNA-seq",
                          col = RNA_seq_col, fill.mountain = rep(RNA_seq_col, 2))
  
  # add DSB to AnnotationTrack
  AT_w_DSBs <- HighlightTrack(trackList = list(AT), range = SrfIcs[n], col = JFly_colors[8], fill = NA, inBackground = FALSE)
  
  # plot
  pdf(file = "tmp.pdf", width = 3.5, height = 2)
    plotTracks(trackList = list(OverlayTrack(trackList = list(DT_RNA_seq, DT_S1_seq_smooth, DT_S1_seq)), AT_w_DSBs), 
               from = start(roi), to = end(roi), chromosome = seqnames(roi), margin = 1, innerMargin = 0)
  dev.off()
  GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/", DSB_locus, ".pdf"))
  
  # record ylims (for manual scale addition)
  Y_axis_ranges <- rbind(Y_axis_ranges,
                         data.frame(DSB_locus = DSB_locus, 
                                    y_raw_min = ylim[1], y_raw_max = ylim[2],
                                    y_smooth_min = ylim_smooth[1], y_smooth_max = ylim_smooth[2],
                                    y_RNA_min = ylim_RNA[1], y_RNA_max = ylim_RNA[2]))
  
}

write.table(x = Y_axis_ranges, file = paste0(plot_dir, "/Y_axis_ranges.txt"), row.names = FALSE)
