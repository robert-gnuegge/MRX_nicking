# info --------------------------------------------------------------------
# purpose:plot coverage
# author: Robert Gnuegge (robert.gnuegge@gmail.com)
# date: 02/13/22
# version: 1.0


# preamble ----------------------------------------------------------------

# set wd to this file's location
wd_path <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(wd_path)

# libraries
library(GenomicRanges)
library(Gviz)
options(ucscChromosomeNames=FALSE)  # for using custom chromosome names

# load helper files and functions
source(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Analysis_scripts/Misc/Genomic_helper_functions.R")
source(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Analysis_scripts/Misc/Misc_helper_functions.R")
source(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Analysis_scripts/Misc/JFly_colors.R")

# read coverage
load(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Data/Bazzano2021/Coverage/exo1_sgs1_coverage.RData")



# plotting ================================================================

plot_dir <- "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Figures/16_Bazzano2021/Coverage"

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


# plot without zoom -------------------------------------------------------

AT <- GenomeAxisTrack(range = GRanges(seqnames = "dsb", ranges = IRanges(start = c(530, 614), width = 1)), col.range = "red")
ylim <- c(0, 8e5)
DT_0 <- DataTrack(range = exo1_sgs1_T0, type = "h", name = "0 min", ylim = ylim)
DT_35 <- DataTrack(range = exo1_sgs1_T35, type = "h", name = "35 min", ylim = ylim)
DT_60 <- DataTrack(range = exo1_sgs1_T60, type = "h", name = "60 min", ylim = ylim)
DT_75 <- DataTrack(range = exo1_sgs1_T75, type = "h", name = "75 min", ylim = ylim)
DT_90 <- DataTrack(range = exo1_sgs1_T90, type = "h", name = "90 min", ylim = ylim)

pdf(file = "tmp.pdf", width = 3.5, height = 4)
plotTracks(trackList = list(AT, DT_0, DT_35, DT_60, DT_75, DT_90), chromosome = "dsb")
dev.off()
GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/Coverage.pdf"))


# plot with zoom ----------------------------------------------------------

ylim <- c(0, 5e4)
DT_0 <- DataTrack(range = exo1_sgs1_T0, type = "h", name = "0 min", ylim = ylim)
DT_35 <- DataTrack(range = exo1_sgs1_T35, type = "h", name = "35 min", ylim = ylim)
DT_60 <- DataTrack(range = exo1_sgs1_T60, type = "h", name = "60 min", ylim = ylim)
DT_75 <- DataTrack(range = exo1_sgs1_T75, type = "h", name = "75 min", ylim = ylim)
DT_90 <- DataTrack(range = exo1_sgs1_T90, type = "h", name = "90 min", ylim = ylim)

pdf(file = "tmp.pdf", width = 3.5, height = 4)
plotTracks(trackList = list(AT, DT_0, DT_35, DT_60, DT_75, DT_90), chromosome = "dsb")
dev.off()
GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/Coverage_zoomed.pdf"))