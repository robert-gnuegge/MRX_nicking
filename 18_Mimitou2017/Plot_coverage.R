# info --------------------------------------------------------------------
# purpose: plot S1-seq coverage for meiotic S. cerevisiae cells (Mimitou et al., 2017)
# author: Robert Gnuegge (robert.gnuegge@gmail.com)
# date: 02/15/22
# version: 1.0


# preamble ----------------------------------------------------------------

# set wd to this file's location
wd_path <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(wd_path)

# libraries
library(GenomicRanges)
library(rtracklayer)
library(Gviz)

# load helper files and functions
source(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Analysis_scripts/Misc/Genomic_helper_functions.R")
source(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Analysis_scripts/Misc/Misc_helper_functions.R")
source(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Analysis_scripts/Misc/JFly_colors.R")

plot_dir <- "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Figures/17_Mimitou2017/Coverage"


# function definitions ----------------------------------------------------

# make Gviz track with separately marked fw and rev scores
# argument: GRanges object, color definitions, ... (e.g. type, name)
# result: OverlayTrack
DataTrack_fw_rev <- function(GRanges, fw_col, rev_col, ...){
  fw <- DataTrack(range = GRanges[strand(GRanges) == "+"], col = fw_col, ...)
  rev <- DataTrack(range = GRanges[strand(GRanges) == "-"], col = rev_col, ...)
  OverlayTrack(trackList = list(fw, rev))
}


# read and process data ===================================================

# read S1-seq coverage data
load(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Data/Mimitou2017/Coverage/S1-seq_coverage.RData")

# read Spo11
Spo11_coverage <- GRanges(import.wig(con = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Resources/Mohibullah2017/GSM2247761_20160706_wt4_2_normalized_hitmap.wig"))

# read Spo11 hot spot locations
hot_spot_data <- read.csv(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Resources/Mohibullah2017/Supplemental_Table_S3.csv", header = TRUE)
# make GRanges
Spo11_hot_spots <- GRanges(seqnames = paste0("chr", as.roman(hot_spot_data$chr)), 
                           ranges = IRanges(start = hot_spot_data$Start, end = hot_spot_data$End),
                           midpoint = hot_spot_data$midpoint,
                           heat = apply(X = cbind(hot_spot_data$Spo11.oligo.hits.in.WT.4.h.sample.1, hot_spot_data$Spo11.oligo.hits.in.WT.4.h.sample.2), MARGIN = 1, FUN = mean))


# plot S1-sed and Spo11-seq data around Spo11 hot spots ===================

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

fw_col <- adjustcolor(col = JFly_colors[2], alpha.f = 0.5)
rev_col <- adjustcolor(col = JFly_colors[3], alpha.f = 0.5)

# smoothing function
hanning <- function(x, n = 51){
  hanning_window <- function(n){
    if(n == 1){
      c <- 1
    }else{
      n <- n - 1
      c <- 0.5 - 0.5 * cos(2 * pi * (0:n) / n)
    }
    return(c / sum(c))
  }
  return(weighted.mean(x = x, w = hanning_window(n)))
}

# make DataTracks
DT_Spo11 <- DataTrack(range = Spo11_coverage, col = "gray", name = "Spo11", type = "l", window = -1, windowSize = 51, aggregation = "hanning")
DT_wt_4 <- DataTrack_fw_rev(GRanges = wt_4_merged, fw_col = fw_col, rev_col = rev_col, name = "wt", type = "l", window = -1, windowSize = 51, aggregation = "hanning")
DT_exo1_4 <- DataTrack_fw_rev(GRanges = exo1_4_merged, fw_col = fw_col, rev_col = rev_col, name = "exo1", type = "l", window = -1, windowSize = 51, aggregation = "hanning")
GAT <- GenomeAxisTrack()

# define plotting roi
hot_spot <- subsetByOverlaps(x = Spo11_hot_spots, ranges = GRanges(seqnames = "chrIV", ranges = IRanges(start = 836000, width = 1000)))
roi <- resize(x = hot_spot, width = 5000, fix = "center")

# plot Mimitou et al., 2017 Fig. 1B region
pdf(file = "tmp.pdf", width = 3, height = 2)
plotTracks(trackList = HighlightTrack(trackList = list(GAT, DT_Spo11, DT_wt_4, DT_exo1_4),
                                      range = subsetByOverlaps(x = Spo11_hot_spots, ranges = roi),
                                      col = NA, fill = adjustcolor(col = JFly_colors[8], alpha.f = 0.25), inBackground = TRUE),
           chromosome = as.character(seqnames(roi)), from = start(roi), to = end(roi), margin = 1, innerMargin = 0)
dev.off()
GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/", gsub(pattern = ":", replacement = "_", x = as.character(hot_spot)), ".pdf"))


# find hot spots that are at least 3.0 kb apart ---------------------------
# max. resection after 4 h seems to be around 1.8 kb, but most resection tracts are < 1.5 kb in wt (Mimitou et al, 2017, Fig. 1E)

hits <- distanceToNearest(x = Spo11_hot_spots)
lonely_hot_spots <- Spo11_hot_spots[queryHits(hits)[mcols(hits)$distance >= 3000]]

# plot hottest lonely hot spots
idx <- order(lonely_hot_spots$heat, decreasing = TRUE)[1:3]

for(n in idx){
  hot_spot <- lonely_hot_spots[n]
  roi <- resize(x = hot_spot, width = 5000, fix = "center")
  
  pdf(file = "tmp.pdf", width = 3, height = 2)
  plotTracks(trackList = HighlightTrack(trackList = list(GAT, DT_Spo11, DT_wt_4, DT_exo1_4),
                                        range = subsetByOverlaps(x = Spo11_hot_spots, ranges = roi),
                                        col = NA, fill = adjustcolor(col = JFly_colors[8], alpha.f = 0.25), inBackground = TRUE),
             chromosome = as.character(seqnames(roi)), from = start(roi), to = end(roi), margin = 1, innerMargin = 0)
  dev.off()
  GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/", gsub(pattern = ":", replacement = "_", x = as.character(hot_spot)), ".pdf"))
}

