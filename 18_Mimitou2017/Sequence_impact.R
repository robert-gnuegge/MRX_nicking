# info --------------------------------------------------------------------
# purpose: analyze nt sequence impact on MRX nicking in meiotic S. cerevisiae cells
# author: Robert Gnuegge (robert.gnuegge@gmail.com)
# date: 02/15/22
# version: 2.0


# preamble ----------------------------------------------------------------

# set wd to this file's location
wd_path <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(wd_path)

# libraries
library(GenomicRanges)
library(BSgenome.Scerevisiae.UCSC.sacCer3)
genome <- BSgenome.Scerevisiae.UCSC.sacCer3  # alias
library(rtracklayer)


# load helper files and functions
source(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Analysis_scripts/Misc/Genomic_helper_functions.R")
source(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Analysis_scripts/Misc/Misc_helper_functions.R")
source(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Analysis_scripts/Misc/JFly_colors.R")

plot_dir <- "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Figures/17_Mimitou2017/Sequence_Impact"


# function definitions ----------------------------------------------------

# calculate position weight matrices
# argument: GRanges, even integer, genome (DNAStringSet)
# result: matrix
calc_PWM <- function(GRanges, width, genome, calc_beakground_PWM = FALSE){
  if(calc_beakground_PWM){
    nick_pos <- GRanges
  }else{
    nick_pos <- GRanges[rep(1:length(GRanges), GRanges$score)]
    # score counts how often a nick occurred at a specific genome position
    # let's repeat score times to retrieve the sequence context score times (see below)
  }
  if(width %% 2 == 1){
    width <- width + 1  
    warning("width must be an even integer. Changing width to ", width)
  }
  nick_context <- flank(x = nick_pos, width = 0.5 * width, both = TRUE)
  nt_sequences <- getSeq(x = genome, names = nick_context)
  PWM <- consensusMatrix(x = nt_sequences, as.prob = TRUE)[1:4, ]
  colnames(PWM) <- c(-(0.5 * width - 1):-0, 0:(0.5 * width - 1))
  return(PWM)
}

# derive sequence with highest frequency base at each position from PWM
# argument: matrix (PWM)
# result: character string
get_highest_freq_sequence <- function(PWM){
  paste0(rownames(PWM)[apply(X = PWM, MARGIN = 2, FUN = which.max)], collapse = "")
}

# plotting function
# argument: 70 nt PWMs
# result: plot
plot_PWM <- function(PWM, PWM_bg = NULL, ylim = NULL, my_colors = c(A = JFly_colors[1], T = JFly_colors[2], G = JFly_colors[3], C = JFly_colors[8]), xlab = "Distance from Nick Site [nt]"){
  if(is.null(ylim)){
    if(!is.null(PWM_bg)){
      ylim <- range(c(PWM, PWM_bg))
    }else{
      ylim <- range(PWM)
    }
  }
  
  # A
  plot(x = 1:dim(PWM)[2], y = PWM["A", ], type = "l", col = my_colors["A"], ylim = ylim,
       xlab = xlab, ylab = "Fraction", xaxt = "n")
  abline(v = 35.5, col = "gray", lty = "dashed")  # nick
  abline(h = c(0.192, 0.308), col = "gray")  # genome average
  
  # add customized x axis
  x <- c(-35:-1, 1:35)
  at <-  which(x %% 5 == 0)
  labels <- x[at]
  axis(side = 1, at = 1:70, labels = FALSE, tcl = -0.15)
  axis(side = 1, at = c(at, 35, 36), labels = NA, tcl = -0.3)
  axis(side = 1, at = at, labels = labels, tick = FALSE)
  
  # add other nts
  points(x = 1:dim(PWM)[2], y = PWM["T", ], type = "l", col = my_colors["T"])
  points(x = 1:dim(PWM)[2], y = PWM["G", ], type = "l", col = my_colors["G"])
  points(x = 1:dim(PWM)[2], y = PWM["C", ], type = "l", col = my_colors["C"])
  
  # add PWM_bg data
  if(!is.null(PWM_bg)){
    points(x = 1:dim(PWM_bg)[2], y = PWM_bg["A", ], type = "l", col = my_colors["A"], lty = "dashed")
    points(x = 1:dim(PWM_bg)[2], y = PWM_bg["T", ], type = "l", col = my_colors["T"], lty = "dashed")
    points(x = 1:dim(PWM_bg)[2], y = PWM_bg["G", ], type = "l", col = my_colors["G"], lty = "dashed")
    points(x = 1:dim(PWM_bg)[2], y = PWM_bg["C", ], type = "l", col = my_colors["C"], lty = "dashed")
    
    # add legend
    y_T <- mean(c(PWM["T", 70], PWM_bg["T", 70]))
    y_A <- mean(c(PWM["A", 70], PWM_bg["A", 70]))
    y_G <- mean(c(PWM["G", 70], PWM_bg["G", 70]))
    y_C <- mean(c(PWM["C", 70], PWM_bg["C", 70]))
    text(x = 70, y = c(y_T, y_A, y_G, y_C), pos = c(1, 3, 3, 1), labels = c("T", "A", "G", "C"), col = JFly_colors[c(2, 1, 3, 8)])
  }else{
    # add legend
    text(x = 70, y = PWM[c("T", "A", "G", "C"), 70], pos = 4, labels = c("T", "A", "G", "C"), col = JFly_colors[c(2, 1, 3, 8)])
  }
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


# plot Spo11 sequence preference ==========================================

# all Spo11 data ----------------------------------------------------------
Spo11_unnormalized <- Spo11_coverage
Spo11_unnormalized$score <- round(Spo11_unnormalized$score / min(Spo11_unnormalized$score))

Spo11_PWM <- calc_PWM(GRanges = Spo11_unnormalized, width = 70, genome = genome)

pdf(file = "tmp.pdf", width = 8, height = 3)
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(1.5, 0.5, 3.8, 2), tcl = -0.25, mgp = c(2.5, 0.5, 0), las = 1)
plot_PWM(PWM = Spo11_PWM, xlab = "Distance from Cleavage Site [nt]")
dev.off()
GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/Spo11-seq_all.pdf"))


# Spo11 data in hot spot regions ------------------------------------------
Spo11_in_hot_spots <- subsetByOverlaps(x = Spo11_unnormalized, ranges = Spo11_hot_spots)

Spo11_in_hot_spots_PWM <- calc_PWM(GRanges = Spo11_in_hot_spots, width = 70, genome = genome)

pdf(file = "tmp.pdf", width = 8, height = 3)
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(1.5, 0.5, 3.8, 2), tcl = -0.25, mgp = c(2.5, 0.5, 0), las = 1)
plot_PWM(PWM = Spo11_in_hot_spots_PWM, xlab = "Distance from Cleavage Site [nt]")
dev.off()
GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/Spo11-seq_in_hot_spots.pdf"))

Spo11_in_hot_spots_PWM


# plot MRX sequence preference ============================================

# all S1-seq data ---------------------------------------------------------
exo1_all_PWM <- calc_PWM(GRanges = exo1_4_1, width = 70, genome = genome)  # let's use 1st replicate data only, as they have more coverage
wt_all_PWM <- calc_PWM(GRanges = wt_4_1, width = 70, genome = genome)

pdf(file = "tmp.pdf", width = 8, height = 3)
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(1.5, 0.5, 3.8, 2), tcl = -0.25, mgp = c(2.5, 0.5, 0), las = 1)
plot_PWM(PWM = exo1_all_PWM, PWM_bg = wt_all_PWM)
dev.off()
GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/MRX_preference_all.pdf"))


# S1-seq data outside Spo11 hot spots -------------------------------------
hits <- findOverlaps(query = Spo11_hot_spots, subject = exo1_4_1)
exo1_cold <- exo1_4_1[-subjectHits(hits)]

hits <- findOverlaps(query = Spo11_hot_spots, subject = wt_4_1)
wt_cold <- wt_4_1[-subjectHits(hits)]

exo1_cold_PWM <- calc_PWM(GRanges = exo1_cold, width = 70, genome = genome)
wt_cold_PWM <- calc_PWM(GRanges = wt_cold, width = 70, genome = genome)

pdf(file = "tmp.pdf", width = 8, height = 3)
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(1.5, 0.5, 4, 2), tcl = -0.25, mgp = c(2.5, 0.5, 0), las = 1)
plot_PWM(PWM = exo1_cold_PWM, PWM_bg = wt_cold_PWM)
dev.off()
GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/MRX_preference_outside_Spo11_hot_spots.pdf"))


# S1-seq data around Spo11 hot spots ----------------------------------------
# max. resection after 4 h seems to be around 1.8 kb, but most resection tracts are < 1.5 kb in wt (Mimitou et al, 2017, Fig. 1E)
around_hot_spots <- Spo11_hot_spots
seqlengths(around_hot_spots) <- seqlengths(genome)[1:16]
start(around_hot_spots) <- around_hot_spots$midpoint - 1500
end(around_hot_spots) <- around_hot_spots$midpoint + 1500
around_hot_spots <- trim(around_hot_spots)

exo1_around <- subsetByOverlaps(x = exo1_4_1, ranges = around_hot_spots)
wt_around <- subsetByOverlaps(x = wt_4_1, ranges = around_hot_spots)

exo1_around_PWM <- calc_PWM(GRanges = exo1_around, width = 70, genome = genome)
wt_around_PWM <- calc_PWM(GRanges = wt_around, width = 70, genome = genome)

pdf(file = "tmp.pdf", width = 8, height = 3)
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(1.5, 0.5, 4, 2), tcl = -0.25, mgp = c(2.5, 0.5, 0), las = 1)
plot_PWM(PWM = exo1_around_PWM, PWM_bg = wt_around_PWM)
dev.off()
GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/MRX_preference_around_Spo11_hot_spots.pdf"))


# S1-seq data around lonely hot spots ---------------------------------------
hits <- distanceToNearest(x = Spo11_hot_spots)
lonely_hot_spots <- Spo11_hot_spots[queryHits(hits)[mcols(hits)$distance >= 3000]]
start(lonely_hot_spots) <- lonely_hot_spots$midpoint - 1500
end(lonely_hot_spots) <- lonely_hot_spots$midpoint + 1500

exo1_lonely <- subsetByOverlaps(x = exo1_4_1, ranges = lonely_hot_spots)
wt_lonely <- subsetByOverlaps(x = wt_4_1, ranges = lonely_hot_spots)

exo1_lonely_PWM <- calc_PWM(GRanges = exo1_lonely, width = 70, genome = genome)
wt_lonely_PWM <- calc_PWM(GRanges = wt_lonely, width = 70, genome = genome)

pdf(file = "tmp.pdf", width = 8, height = 3)
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(1.5, 0.5, 4, 2), tcl = -0.25, mgp = c(2.5, 0.5, 0), las = 1)
plot_PWM(PWM = exo1_lonely_PWM, PWM_bg = wt_lonely_PWM)
dev.off()
GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/MRX_preference_lonely_Spo11_hot_spots.pdf"))


# S1-seq data outside lonely hot spots ------------------------------------
lonely_hot_spots <- as_nt_resolved_GRanges(GRanges = lonely_hot_spots)
hits <- findOverlaps(query = Spo11_hot_spots, subject = lonely_hot_spots)
outside_lonely_hot_spots <- lonely_hot_spots[-subjectHits(hits)]

exo1_outside_lonely <- subsetByOverlaps(x = exo1_4_1, ranges = outside_lonely_hot_spots)
wt_outside_lonely <- subsetByOverlaps(x = wt_4_1, ranges = outside_lonely_hot_spots)

exo1_outside_lonely_PWM <- calc_PWM(GRanges = exo1_outside_lonely, width = 70, genome = genome)
wt_outside_lonely_PWM <- calc_PWM(GRanges = wt_outside_lonely, width = 70, genome = genome)

pdf(file = "tmp.pdf", width = 8, height = 3)
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(1.5, 0.5, 4, 2), tcl = -0.25, mgp = c(2.5, 0.5, 0), las = 1)
plot_PWM(PWM = exo1_outside_lonely_PWM, PWM_bg = wt_outside_lonely_PWM)
dev.off()
GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/MRX_preference_outside_lonely_Spo11_hot_spots.pdf"))


# S1-seq data around hottest spots ----------------------------------------
hottest_hot_spots <- Spo11_hot_spots[Spo11_hot_spots$heat >= quantile(x = Spo11_hot_spots$heat, probs = 0.75)]
start(hottest_hot_spots) <- hottest_hot_spots$midpoint - 1500
end(hottest_hot_spots) <- hottest_hot_spots$midpoint + 1500

exo1_hottest <- subsetByOverlaps(x = exo1_4_1, ranges = hottest_hot_spots)
wt_hottest <- subsetByOverlaps(x = wt_4_1, ranges = hottest_hot_spots)

exo1_hottest_PWM <- calc_PWM(GRanges = exo1_hottest, width = 70, genome = genome)
wt_hottest_PWM <- calc_PWM(GRanges = wt_hottest, width = 70, genome = genome)

pdf(file = "tmp.pdf", width = 8, height = 3)
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(1.5, 0.5, 4, 2), tcl = -0.25, mgp = c(2.5, 0.5, 0), las = 1)
plot_PWM(PWM = exo1_hottest_PWM, PWM_bg = wt_hottest_PWM)
dev.off()
GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/MRX_preference_hottest_Spo11_hot_spots.pdf"))


# S1-seq data outside hottest hot spots -----------------------------------
hottest_hot_spots <- as_nt_resolved_GRanges(GRanges = hottest_hot_spots)
hits <- findOverlaps(query = Spo11_hot_spots, subject = hottest_hot_spots)
outside_hottest_hot_spots <- hottest_hot_spots[-subjectHits(hits)]

exo1_outside_hottest <- subsetByOverlaps(x = exo1_4_1, ranges = outside_hottest_hot_spots)
wt_outside_hottest <- subsetByOverlaps(x = wt_4_1, ranges = outside_hottest_hot_spots)

exo1_outside_hottest_PWM <- calc_PWM(GRanges = exo1_outside_hottest, width = 70, genome = genome)
wt_outside_hottest_PWM <- calc_PWM(GRanges = wt_outside_hottest, width = 70, genome = genome)

pdf(file = "tmp.pdf", width = 8, height = 3)
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(1.5, 0.5, 4, 2), tcl = -0.25, mgp = c(2.5, 0.5, 0), las = 1)
plot_PWM(PWM = exo1_outside_hottest_PWM, PWM_bg = wt_outside_hottest_PWM)
dev.off()
GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/MRX_preference_outside_hottest_Spo11_hot_spots.pdf"))
