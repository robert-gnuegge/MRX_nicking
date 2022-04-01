# info --------------------------------------------------------------------
# purpose: analyze nt sequence impact on MRX nicking
# author: Robert Gnuegge (robert.gnuegge@gmail.com)
# date: 01/27/22
# version: 2.0


# preamble ----------------------------------------------------------------

# set wd to this file's location
wd_path <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(wd_path)

# libraries
library(GenomicRanges)
library(BSgenome)
library(ggplot2)
library(ggseqlogo)

# read modified S. cerevisiae genome
genome <- import(con = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Data/S_cerevisiae_reference_genomes/S288C_R64-2-1_W303_SNPs_MATa_hocs2SrfIcs_hml_hmr.fasta")

# load helper files and functions
source(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Analysis_scripts/Misc/Genomic_helper_functions.R")
source(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Analysis_scripts/Misc/S_cerevisiae_SrfI_cut_sites.R")
source(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Analysis_scripts/Misc/Misc_helper_functions.R")
source(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Analysis_scripts/Misc/JFly_colors.R")

plot_dir <- "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Figures/08_Sequence_impact/LSY4602-20C_1"


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


# process data for plotting ===============================================

# read S1-seq coverage
load(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Data/S1-seq/LSY4602-20C_1/LSY4602-20C_1_unnormalized.RData")
# use unnormalized coverage data (see explanation below)!
load(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Data/S1-seq/LSY4602-20C_1/LSY4602_20C_1_resection_extend.RData")

# keep only regions with correct orientation w.r.t DSBs
LSY4602_20C_1 <- subsetByIntersect(subject = LSY4602_20C_1_S1_seq_unnormalized, query = DSB_regions(DSBs = SrfIcs[-c(9, 17)], region_width = 4000, up_rev_down_fw = TRUE))

# retain GRanges in regions around DSBs where resection was detected
# other regions are excluded to prevent pattern "dilution" by background noise
LSY4602_20C_1 <- subsetByIntersect(subject = LSY4602_20C_1, query = resection_extend_LSY4602_20C_1)

# make nt-resolved GRanges
LSY4602_20C_1 <- as_nt_resolved_GRanges(GRanges = LSY4602_20C_1)

# remove SrfIcs sequence regions
# they give a lot of S1-seq signal, but derive from unprocessed DSBs and not from MRX nicking
DSBs <- DSB_regions(DSBs = SrfIcs, region_width = 8)
unique(getSeq(x = genome, names = DSBs))  # all SrfIcs sequence

LSY4602_20C_1 <- LSY4602_20C_1[-subjectHits(findOverlaps(query = DSBs, subject = LSY4602_20C_1))]

# calculate PWM and background PWM
PWM <- calc_PWM(GRanges = LSY4602_20C_1, width = 100, genome = genome)
PWM_bg <- calc_PWM(GRanges = LSY4602_20C_1, width = 100, genome = genome, calc_beakground_PWM = TRUE)


# plot nt fractions -------------------------------------------------------

# plotting function
plot_nt_fractions <- function(PWM, PWM_bg){
  plot(x = 1:dim(PWM)[2], y = PWM["A", ], type = "l", col = JFly_colors[1], ylim = range(c(PWM, PWM_bg)),
       xlab = "Distance from Nick Site [nt]", ylab = "Fraction", xaxt = "n")
  points(x = 1:dim(PWM_bg)[2], y = PWM_bg["A", ], type = "l", col = JFly_colors[1], lty = "dashed")
  
  axis(side = 1, at = 1:100, labels = FALSE, tcl = -0.2)
  axis(side = 1, at = c(0:5 * 10, 5:9 * 10 + 1), labels = c(NA, -4:-1 * 10, NA, NA, 1:4 * 10), tcl = -0.3)
  axis(side = 1, at = 50.5, labels = 0, tick = FALSE)
  
  points(x = 1:dim(PWM)[2], y = PWM["T", ], type = "l", col = JFly_colors[2])
  points(x = 1:dim(PWM_bg)[2], y = PWM_bg["T", ], type = "l", col = JFly_colors[2], lty = "dashed")
  
  points(x = 1:dim(PWM)[2], y = PWM["G", ], type = "l", col = JFly_colors[3])
  points(x = 1:dim(PWM_bg)[2], y = PWM_bg["G", ], type = "l", col = JFly_colors[3], lty = "dashed")
  
  points(x = 1:dim(PWM)[2], y = PWM["C", ], type = "l", col = JFly_colors[4])
  points(x = 1:dim(PWM_bg)[2], y = PWM_bg["C", ], type = "l", col = JFly_colors[4], lty = "dashed")
  
  abline(v = 50.5, col = "gray")
  
  legend(x = "top", legend = c("A", "T", "G", "C"), col = JFly_colors[1:4], inset = -0.275, lty = "solid", ncol = 4, xpd = TRUE)
}

pdf(file = "tmp.pdf", width = 6, height = 3.5)
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(1.5, 0.5, 1, 1.7), tcl = -0.3, mgp = c(2.5, 0.6, 0), las = 1)
 plot_nt_fractions(PWM = PWM, PWM_bg = PWM_bg)
dev.off()
GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/Sequence_preference.pdf"))


# plot sequence logo ------------------------------------------------------

col_scheme = make_col_scheme(chars = c("A", "T", "G", "C"), cols = JFly_colors[1:4])  # color assignment

# plotting function
plot_seq_logo <- function(PWM){
  ggplot() + geom_logo(PWM[, which(colnames(PWM) == "-8"):which(colnames(PWM) == "8")], method = 'prob', col_scheme = col_scheme) + theme_logo() +
    scale_x_continuous(breaks = c(1, 3, 5, 7, 12, 14, 16, 18), labels = c(-8, -6, -4, -2, 2, 4, 6, 8)) + 
    labs(x = "Distance from Nick Site [nt]", y = "Fraction") +
    theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
    theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
    geom_vline(xintercept = 9.5, linetype = "solid", col = "gray")
}

plot_seq_logo(PWM = PWM)
ggsave(filename = "tmp.pdf", device = "pdf", width = 3.25, height = 2, units = "in")
GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/Sequence_logo.pdf"))

plot_seq_logo(PWM = PWM_bg)
ggsave(filename = "tmp.pdf", device = "pdf", width = 3.25, height = 2, units = "in")
GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/Sequence_logo_background.pdf"))


# compare bg PWM to S. cerevisiae genome-wide nt frequencies --------------

# S. cerevisiae genome-wide N content
tmp <- alphabetFrequency(x = genome)
tmp <- tmp[c(1:16), c("A", "T", "G", "C")]
apply(X = tmp, MARGIN = 2, FUN = sum) / sum(tmp)
#         A         T         G         C 
# 0.3089952 0.3079593 0.1913379 0.1917076

# N content around SrfIcs
apply(X = PWM_bg, MARGIN = 1, FUN = mean)
#         A        T         G         C
# 0.2922961 0.3054438 0.1954275 0.2068327