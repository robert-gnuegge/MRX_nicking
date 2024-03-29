# info --------------------------------------------------------------------
# purpose: generate S1-seq coverage plots across all S. cerevisiae DNAs
# author: Robert Gnuegge (robert.gnuegge@gmail.com)
# date: 01/26/22
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


# sum scores in bins ------------------------------------------------------
# argument: GRanges object with mcols called "fw" and "rev", GRanges object defining bins
# result: GRanges object
sum_scores_in_bins <- function(GRanges, bins){
  out <- bins
  hits <- findOverlaps(query = out, subject = GRanges)
  fw <- aggregate(GRanges, hits, score = sum(fw))
  out$fw <- fw$score
  rev <- aggregate(GRanges, hits, score = sum(rev))
  out$rev <- rev$score
  return(out)
}


# process all samples =====================================================

# read S1-seq coverage
load(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Data/S1-seq/LSY4377-12B_4377-15A_merged/LSY4377-12B_LSY4377-15A_w_chrM_and_2micron_merged.RData")


# process data for plotting -----------------------------------------------

# derive GRanges objects with 1-kb binned separate fw and rev score mcols (required for grouped plotting with Gviz)
bins <- unlist(tile(x = SacCer_chromosomes, width = 1000))
LSY4377_12B_0 <- sum_scores_in_bins(GRanges = as_GRanges_with_fw_rev_scores(GRanges = LSY4377_12B_0_merged_S1_seq_w_chrM_and_2micron, rev_negative = TRUE), bins = bins)
LSY4377_12B_1 <- sum_scores_in_bins(GRanges = as_GRanges_with_fw_rev_scores(GRanges = LSY4377_12B_1_merged_S1_seq_w_chrM_and_2micron, rev_negative = TRUE), bins = bins)
LSY4377_12B_2 <- sum_scores_in_bins(GRanges = as_GRanges_with_fw_rev_scores(GRanges = LSY4377_12B_2_merged_S1_seq_w_chrM_and_2micron, rev_negative = TRUE), bins = bins)
LSY4377_12B_4 <- sum_scores_in_bins(GRanges = as_GRanges_with_fw_rev_scores(GRanges = LSY4377_12B_4_merged_S1_seq_w_chrM_and_2micron, rev_negative = TRUE), bins = bins)
LSY4377_15A_4 <- sum_scores_in_bins(GRanges = as_GRanges_with_fw_rev_scores(GRanges = LSY4377_15A_4_merged_S1_seq_w_chrM_and_2micron, rev_negative = TRUE), bins = bins)


# find y range that covers all sample score ranges ------------------------
y_range <- range(c(LSY4377_12B_0$fw, LSY4377_12B_0$rev, LSY4377_12B_1$fw, LSY4377_12B_1$rev,
                   LSY4377_12B_2$fw, LSY4377_12B_2$rev, LSY4377_12B_4$fw, LSY4377_12B_4$rev,
                   LSY4377_15A_4$fw, LSY4377_15A_4$rev))
y_range <- c(- max(abs(pretty(y_range))), max(abs(pretty(y_range))))  # makes y_range "nicer"


# make DataTracks ---------------------------------------------------------

# set global Gviz and DataTrack parameters for plotting
options(Gviz.scheme = "default")
scheme <- getScheme()  # copy current scheme
scheme$GdObject$showAxis <- FALSE
scheme$GdObject$showTitle <- FALSE
scheme$DataTrack$type <- "h"
scheme$DataTrack$groups <- c("fw", "rev")
scheme$DataTrack$col <- c(JFly_colors[2], JFly_colors[3])
scheme$DataTrack$legend <- FALSE
addScheme(scheme, "MyScheme")  # define new scheme
options(Gviz.scheme = "MyScheme")  # set new scheme

# construct DataTrack objects
LSY4377_12B_0_DT <- DataTrack(LSY4377_12B_0, ylim = y_range)
LSY4377_12B_1_DT <- DataTrack(LSY4377_12B_1, ylim = y_range)
LSY4377_12B_2_DT <- DataTrack(LSY4377_12B_2, ylim = y_range)
LSY4377_12B_4_DT <- DataTrack(LSY4377_12B_4, ylim = y_range)
LSY4377_15A_4_DT <- DataTrack(LSY4377_15A_4, ylim = y_range)


# plotting ----------------------------------------------------------------
# plot DataTrack for each sample and chromosome
# scale plot width to relative chromosome size
max_plot_width <- 4

# plotting dir
plot_dir <- "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Figures/03_S1-seq_genome-wide/Individual"

# iterate through all samples and chromosomes
for (sample in c(paste0("LSY4377_12B_", c(0, 1, 2, 4)), "LSY4377_15A_4")){
  dir.create(paste0(plot_dir, "/", sample, "_merged"), recursive = TRUE, showWarnings = FALSE)  # create sub directory for each sample
  for (chr in SacCer_chromosomes_df$chr){
    plot_width <- SacCer_chromosomes_df$length[SacCer_chromosomes_df$chr == chr] / max(SacCer_chromosomes_df$length) * max_plot_width
    pdf(file = "tmp.pdf", width = plot_width, height = 2)
      plotTracks(get(paste0(sample, "_DT")), chromosome = chr, margin = 0, innerMargin = 1, lwd = 1)
    dev.off()
    GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/", sample, "_merged/", sample, "_merged_", chr, ".pdf"))
  }
}