# info --------------------------------------------------------------------
# purpose: generate MNase-seq coverage plots around genomic SrfIcs for both replicates
# author: Robert Gnuegge (robert.gnuegge@gmail.com)
# date: 02/25/22
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


# process MNase-seq data for plotting =====================================

# read MNase-seq coverage
load(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Data/MNase-seq/LSY4377-12B_4377-15A/LSY4377-12B_LSY4377-15A_trimmed.RData")
load(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Data/MNase-seq/LSY4377-12B_4377-15A_rep2/LSY4377-12B_LSY4377-15A_rep2_trimmed.RData")

# retain only -/+ 5 kb regions around SrfIcs (for faster data processing)
roi <- DSB_regions(DSBs = SrfIcs, region_width = 10000)

MNase_seq_0_rep1 <- subsetByIntersect(subject = LSY4377_12B_0_MNase_seq_trimmed, query = roi)
MNase_seq_1_rep1 <- subsetByIntersect(subject = LSY4377_12B_1_MNase_seq_trimmed, query = roi)
MNase_seq_2_rep1 <- subsetByIntersect(subject = LSY4377_12B_2_MNase_seq_trimmed, query = roi)
MNase_seq_4_rep1 <- subsetByIntersect(subject = LSY4377_12B_4_MNase_seq_trimmed, query = roi)

MNase_seq_0_rep2 <- subsetByIntersect(subject = LSY4377_12B_0_rep2_MNase_seq_trimmed, query = roi)
MNase_seq_1_rep2 <- subsetByIntersect(subject = LSY4377_12B_1_rep2_MNase_seq_trimmed, query = roi)
MNase_seq_2_rep2 <- subsetByIntersect(subject = LSY4377_12B_2_rep2_MNase_seq_trimmed, query = roi)
MNase_seq_4_rep2 <- subsetByIntersect(subject = LSY4377_12B_4_rep2_MNase_seq_trimmed, query = roi)


# plotting ================================================================

plot_dir <- "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Figures/02_Replicate_correlation/MNase-seq/side_by_side"

# set global Gviz parameters for plotting
options(Gviz.scheme="default")
scheme <- getScheme()  # copy current scheme
scheme$GdObject$showTitle <- FALSE
scheme$GdObject$background.title <- "white"
scheme$GdObject$col.axis <- "black"
addScheme(scheme, "MyScheme")  # define new scheme
options(Gviz.scheme = "MyScheme")  # set new scheme

# plotting colors
MNase_seq_col <- gray(level = 0.7)
MNase_seq_border_col <- gray(level = 0.5)

# AnnotationTrack ---------------------------------------------------------
# data downloaded from yeastmine.yeastgenome.org and curated with Adjust_chromosomal_features.R
load(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Data/S_cerevisiae_chromosomal_features/S_cerevisiae_genome_features.RData")

# adjust for use with AnnotationTrack
names(mcols(all_features))[1] <- "id"
all_features <- all_features[!(all_features$id == "") & all_features$type == "ORF"]

AT <- AnnotationTrack(range = all_features, name = NULL, showFeatureId = FALSE, cex = 0.67,  # featureAnnotation = "id"
                      arrowHeadMaxWidth = 10, fill = "white", col = "gray", fontcolor.item = "black")

# plot regions around all SrfIcs ==========================================
n <- 10
for (n in 1:length(SrfIcs)){

  # define name for data and figure saving
  DSB_locus <- SrfIcs[n]
  DSB_locus <- gsub(pattern = ":", replacement = "_", x = as.character(DSB_locus))
  cat("\nPlotting", DSB_locus, "...")
  
  # define plot area
  roi <- DSB_regions(DSBs = SrfIcs, region_width = 3000)[n]
  
  # get GRanges in plot area, calculate ylim, and construct DataTracks
  MNase_seq_0_rep1_roi <- subsetByOverlaps(x = MNase_seq_0_rep1, ranges = roi)
  MNase_seq_1_rep1_roi <- subsetByOverlaps(x = MNase_seq_1_rep1, ranges = roi)
  MNase_seq_2_rep1_roi <- subsetByOverlaps(x = MNase_seq_2_rep1, ranges = roi)
  MNase_seq_4_rep1_roi <- subsetByOverlaps(x = MNase_seq_4_rep1, ranges = roi)
  
  MNase_seq_0_rep2_roi <- subsetByOverlaps(x = MNase_seq_0_rep2, ranges = roi)
  MNase_seq_1_rep2_roi <- subsetByOverlaps(x = MNase_seq_1_rep2, ranges = roi)
  MNase_seq_2_rep2_roi <- subsetByOverlaps(x = MNase_seq_2_rep2, ranges = roi)
  MNase_seq_4_rep2_roi <- subsetByOverlaps(x = MNase_seq_4_rep2, ranges = roi)
  
  ylim_MNase <- range(c(MNase_seq_0_rep1_roi, MNase_seq_1_rep1_roi, MNase_seq_2_rep1_roi, MNase_seq_4_rep1_roi,
                        MNase_seq_0_rep2_roi, MNase_seq_1_rep2_roi, MNase_seq_2_rep2_roi, MNase_seq_4_rep2_roi)$score)
  ylim_MNase <- c(min(ylim_MNase) - 0.05 * max(pretty(ylim_MNase)), max(pretty(ylim_MNase)))  # adjust for prettier plotting
  
  DT_MNase_seq_0_rep1 <- DataTrack(range = MNase_seq_0_rep1_roi, type = "polygon", col = MNase_seq_border_col, fill.mountain = rep(MNase_seq_col, 2), name = "0 h", ylim = ylim_MNase)
  DT_MNase_seq_1_rep1 <- DataTrack(range = MNase_seq_1_rep1_roi, type = "polygon", col = MNase_seq_border_col, fill.mountain = rep(MNase_seq_col, 2), name = "1 h", ylim = ylim_MNase)
  DT_MNase_seq_2_rep1 <- DataTrack(range = MNase_seq_2_rep1_roi, type = "polygon", col = MNase_seq_border_col, fill.mountain = rep(MNase_seq_col, 2), name = "2 h", ylim = ylim_MNase)
  DT_MNase_seq_4_rep1 <- DataTrack(range = MNase_seq_4_rep1_roi, type = "polygon", col = MNase_seq_border_col, fill.mountain = rep(MNase_seq_col, 2), name = "4 h", ylim = ylim_MNase)
  
  DT_MNase_seq_0_rep2 <- DataTrack(range = MNase_seq_0_rep2_roi, type = "polygon", col = MNase_seq_border_col, fill.mountain = rep(MNase_seq_col, 2), name = "0 h", ylim = ylim_MNase)
  DT_MNase_seq_1_rep2 <- DataTrack(range = MNase_seq_1_rep2_roi, type = "polygon", col = MNase_seq_border_col, fill.mountain = rep(MNase_seq_col, 2), name = "1 h", ylim = ylim_MNase)
  DT_MNase_seq_2_rep2 <- DataTrack(range = MNase_seq_2_rep2_roi, type = "polygon", col = MNase_seq_border_col, fill.mountain = rep(MNase_seq_col, 2), name = "2 h", ylim = ylim_MNase)
  DT_MNase_seq_4_rep2 <- DataTrack(range = MNase_seq_4_rep2_roi, type = "polygon", col = MNase_seq_border_col, fill.mountain = rep(MNase_seq_col, 2), name = "4 h", ylim = ylim_MNase)

  # plot
  pdf(file = "tmp.pdf", width = 3.5, height = 4)
    plotTracks(trackList = list(DT_MNase_seq_0_rep1, DT_MNase_seq_1_rep1, DT_MNase_seq_2_rep1, DT_MNase_seq_4_rep1, AT), 
               from = start(roi), to = end(roi), chromosome = seqnames(roi), margin = 0, sizes = c(rep(0.94/4, 4), 0.06))
  dev.off()
  GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/", DSB_locus, "_rep1.pdf"))
  
  pdf(file = "tmp.pdf", width = 3.5, height = 4)
  plotTracks(trackList = list(DT_MNase_seq_0_rep2, DT_MNase_seq_1_rep2, DT_MNase_seq_2_rep2, DT_MNase_seq_4_rep2, AT), 
             from = start(roi), to = end(roi), chromosome = seqnames(roi), margin = 0, sizes = c(rep(0.94/4, 4), 0.06))
  dev.off()
  GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/", DSB_locus, "_rep2.pdf"))

}