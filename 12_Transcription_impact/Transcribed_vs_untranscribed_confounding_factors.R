# info --------------------------------------------------------------------
# purpose: exclude confounding factors (sequence, melting) in (un-)transcribed regions
# author: Robert Gnuegge (robert.gnuegge@gmail.com)
# date: 07/25/22
# version: 1.0


# preamble ----------------------------------------------------------------

# set wd to this file's location
wd_path <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(wd_path)

# libraries
library(GenomicRanges)
library(rtracklayer)
library(BSgenome)
library(beanplot)

# read modified S. cerevisiae genome
genome <- import(con = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Data/S_cerevisiae_reference_genomes/S288C_R64-2-1_W303_SNPs_MATa_hocs2SrfIcs_hml_hmr.fasta")

# load helper files and functions
source(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Analysis_scripts/Misc/Genomic_helper_functions.R")
source(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Analysis_scripts/Misc/S_cerevisiae_SrfI_cut_sites.R")
source(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Analysis_scripts/Misc/Misc_helper_functions.R")
source(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Analysis_scripts/Misc/JFly_colors.R")

plot_dir <- "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Figures/11_Transcription_impact/LSY4377-12B_4377-15A_merged/Exclude_confounding_factors"


# functions ---------------------------------------------------------------

add_significance <- function(x0, x1, y, text, text_offset = -2.5, x_shrink = 0.05){
  segments(x0 = (x0 + x_shrink), y0 = y, x1 = (x1 - x_shrink), y1 = y)
  text(x = x0 + diff(x = c(x0, x1))/2, y = y + text_offset, labels = text, pos = 3, xpd = TRUE)
}


# load S1-seq data --------------------------------------------------------

load(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Data/S1-seq/LSY4377-12B_4377-15A_merged/LSY4377-12B_LSY4377-15A_merged_around_DSBs.RData")

# only include DSBs that are within an ORF
DSBs <- SrfIcs[-c(9, 17)]  # exclude DSBs in duplicated genome region
# only keep regions with correct orientation w.r.t DSBs
roi <- DSB_regions(DSBs = DSBs, region_width = 5000, up_rev_down_fw = TRUE)

LSY4377_12B_4 <- subsetByIntersect(subject = LSY4377_12B_4_merged_S1_seq, query = roi)

# load RNA-seq data -------------------------------------------------------
# data are from Maya-Miles et al, 2019 (pmid: 31331360)

# there are two biological replicates
roi <- DSB_regions(DSBs = DSBs, region_width = 20000)
RNA_seq_1 <- import("/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Resources/Maya-Miles2019/GSM3567364_w303_rep1.bigwig",
                    which = roi)
RNA_seq_2 <- import("/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Resources/Maya-Miles2019/GSM3567365_w303_rep2.bigwig",
                    which = roi)
# let's take the average RNA-seq score
RNA_seq_1 <- subsetByIntersect(subject = RNA_seq_1, query = RNA_seq_2)  # to let both data sets have the same granges
RNA_seq_2 <- subsetByIntersect(subject = RNA_seq_2, query = RNA_seq_1)
all(granges(RNA_seq_1) == granges(RNA_seq_2))
RNA_seq <- RNA_seq_1
RNA_seq$score <- apply(X = cbind(RNA_seq_1$score, RNA_seq_2$score), MARGIN = 1, FUN = mean)


# load genome features ----------------------------------------------------
# data downloaded from yeastmine.yeastgenome.org and curated with Adjust_chromosomal_features.R

load(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Data/S_cerevisiae_chromosomal_features/S_cerevisiae_genome_features.RData")
all_features <- all_features[all_features$type == "ORF"]


# determine transcribed vs untranscribed ----------------------------------

tmp <- LSY4377_12B_4
hits <- findOverlaps(subject = RNA_seq, query = tmp, ignore.strand = TRUE)
tmp$RNA_seq <- mcols(RNA_seq[subjectHits(hits)])$score
tmp$transcribed <- tmp$RNA_seq > 0

transcribed <- reduce(tmp[tmp$transcribed])
mcols(transcribed) <- mcols(subsetByOverlaps(x = tmp, ranges = resize(x = transcribed, width = 1)))[, 3:4]

untranscribed <- reduce(tmp[!tmp$transcribed])
mcols(untranscribed) <- mcols(subsetByOverlaps(x = tmp, ranges = resize(x = untranscribed, width = 1)))[, 3:4]


# Difference in DSB formation kinetics? -----------------------------------

DSB_kinetics_ranks <- data.frame(Transcription = c(rep("Untranscribed", length(untranscribed)), rep("Transcribed", length(transcribed))),
                                 rank = c(untranscribed$DSB_kinetics_rank, transcribed$DSB_kinetics_rank))

# print to PDF ----------------------------------
pdf(file = "tmp.pdf", width=3, height=3)
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(3.7, 1.2, 2.4, 2), xpd = TRUE, las=1, tcl = -0.25, mgp = c(2, 0.5, 0))

at <- 1:2
bp <- boxplot(at = at, rank ~ Transcription, data = DSB_kinetics_ranks, outline = FALSE, col = NA, border = "gray",
              xlab = NA, ylab = "DSB Kinetic Rank")

stripchart(at = at, rank ~ Transcription, data = DSB_kinetics_ranks, 
           vertical = TRUE, pch = 20, method = "jitter", add = TRUE, col = gray(level = 0, alpha = 0.5))

wilcox.test(rank ~ Transcription, data = DSB_kinetics_ranks)

add_significance(x0 = 1, x1 = 2, y = 1.025 * par("usr")[4], text_offset = -0.25, x_shrink = 0, text = expression("p = 0.79"))

dev.off()
GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/DSB_formation.pdf"))


# Difference in sequence composition? -------------------------------------

# get genome sequences in (un)transcribed regions
seq_transcribed <- getSeq(x = genome, names = transcribed)
seq_untranscribed <- getSeq(x = genome, names = untranscribed)

# count appearances of motif (per kbp)
nick_motif <- data.frame(freq = c(vcountPattern(pattern = DNAString("S"), subject = seq_transcribed, fixed = FALSE) / width(seq_transcribed) * 1000, 
                                  vcountPattern(pattern = DNAString("S"), subject = seq_untranscribed, fixed = FALSE) / width(seq_untranscribed) * 1000,
                                  vcountPattern(pattern = DNAString("WSW"), subject = seq_transcribed, fixed = FALSE) / width(seq_transcribed) * 1000, 
                                  vcountPattern(pattern = DNAString("WSW"), subject = seq_untranscribed, fixed = FALSE) / width(seq_untranscribed) * 1000,
                                  vcountPattern(pattern = DNAString("WWSWW"), subject = seq_transcribed, fixed = FALSE) / width(seq_transcribed) * 1000, 
                                  vcountPattern(pattern = DNAString("WWSWW"), subject = seq_untranscribed, fixed = FALSE) / width(seq_untranscribed) * 1000,
                                  vcountPattern(pattern = DNAString("WWWSWWW"), subject = seq_transcribed, fixed = FALSE) / width(seq_transcribed) * 1000, 
                                  vcountPattern(pattern = DNAString("WWWSWWW"), subject = seq_untranscribed, fixed = FALSE) / width(seq_untranscribed) * 1000,
                                  vcountPattern(pattern = DNAString("WWWWSWWWW"), subject = seq_transcribed, fixed = FALSE) / width(seq_transcribed) * 1000, 
                                  vcountPattern(pattern = DNAString("WWWWSWWWW"), subject = seq_untranscribed, fixed = FALSE) / width(seq_untranscribed) * 1000), 
                         pattern = c(rep("S", length(c(seq_transcribed, seq_untranscribed))),
                                     rep("WSW", length(c(seq_transcribed, seq_untranscribed))),
                                     rep("WWSWW", length(c(seq_transcribed, seq_untranscribed))),
                                     rep("WWWSWWW", length(c(seq_transcribed, seq_untranscribed))),
                                     rep("WWWWSWWWW", length(c(seq_transcribed, seq_untranscribed)))),
                         Transcription = c(rep("Transcribed", length(seq_transcribed)), rep("Untranscribed", length(seq_untranscribed))))

# print to PDF ----------------------------------
pdf(file = "tmp.pdf", width=4, height=4)
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(-1.7, 0.7, 0.3, 1.6), xpd = TRUE, las=1, tcl = -0.25, mgp = c(2.5, 0.5, 0))

at <- c(1, 2, 3.5, 4.5, 6, 7, 8.5, 9.5, 11, 12)
bp <- boxplot(at = at, freq ~ Transcription + pattern, data = nick_motif, ylim = range(nick_motif$freq),
              outline = FALSE, col = NA, border = "gray", xlab = NA, ylab = "  Motif frequency per 1 kbp", names = NA)

stripchart(at = at, freq ~ Transcription + pattern, data = nick_motif, 
           vertical = TRUE, pch = 20, method = "jitter", add = TRUE, col = gray(level = 0, alpha = 0.5))

# labels
text(x = at, y = par("usr")[3] - 0 * (par("usr")[4] - par("usr")[3]), labels = c("+", "-"), pos = 1)
text(x = -1, y = par("usr")[3] - 0 * (par("usr")[4] - par("usr")[3]), labels = "Transcr.", pos = 1)

add_significance(x0 = at[1], x1 = at[2], y = par("usr")[3] - 0.14 * (par("usr")[4] - par("usr")[3]), x_shrink = -0.1, text = NA)
text(x = 0.5 * sum(at[1:2]), y = par("usr")[3] - 0.17 * (par("usr")[4] - par("usr")[3]),
     adj = c(1, 1), label = "S", srt = 45)

add_significance(x0 = at[3], x1 = at[4], y = par("usr")[3] - 0.14 * (par("usr")[4] - par("usr")[3]), x_shrink = -0.1, text = NA)
text(x = 0.5 * sum(at[3:4]), y = par("usr")[3] - 0.17 * (par("usr")[4] - par("usr")[3]),
     adj = c(1, 1), label = "WSW", srt = 45)

add_significance(x0 = at[5], x1 = at[6], y = par("usr")[3] - 0.14 * (par("usr")[4] - par("usr")[3]), x_shrink = -0.1, text = NA)
text(x = 0.5 * sum(at[5:6]), y = par("usr")[3] - 0.17 * (par("usr")[4] - par("usr")[3]),
     adj = c(1, 1), label = "WWSWW", srt = 45)

add_significance(x0 = at[7], x1 = at[8], y = par("usr")[3] - 0.14 * (par("usr")[4] - par("usr")[3]), x_shrink = -0.1, text = NA)
text(x = 0.5 * sum(at[7:8]), y = par("usr")[3] - 0.17 * (par("usr")[4] - par("usr")[3]),
     adj = c(1, 1), label = "WWWSWWW", srt = 45)

add_significance(x0 = at[9], x1 = at[10], y = par("usr")[3] - 0.14 * (par("usr")[4] - par("usr")[3]), x_shrink = -0.1, text = NA)
text(x = 0.5 * sum(at[9:10]), y = par("usr")[3] - 0.17 * (par("usr")[4] - par("usr")[3]),
     adj = c(1, 1), label = "WWWWSWWWW", srt = 45)

# significance
wilcox.test(freq ~ Transcription, data = nick_motif[nick_motif$pattern == "S", ])
add_significance(x0 = at[1], x1 = at[2], y = 1.025 * par("usr")[4], text_offset = -0.01, x_shrink = -0.15, text = NA)
text(x = 0.5 * sum(at[1:2]), y = 1.045 * par("usr")[4], adj = c(0, 0), label = expression("p = 0.0015"), srt = 45)

wilcox.test(freq ~ Transcription, data = nick_motif[nick_motif$pattern == "WSW", ])
add_significance(x0 = at[3], x1 = at[4], y = 1.025 * par("usr")[4], text_offset = -0.01, x_shrink = -0.15, text = NA)
text(x = 0.5 * sum(at[3:4]), y = 1.045 * par("usr")[4], adj = c(0, 0), label = expression("p = 0.99"), srt = 45)

wilcox.test(freq ~ Transcription, data = nick_motif[nick_motif$pattern == "WWSWW", ])
add_significance(x0 = at[5], x1 = at[6], y = 1.025 * par("usr")[4], text_offset = -0.01, x_shrink = -0.15, text = NA)
text(x = 0.5 * sum(at[5:6]), y = 1.045 * par("usr")[4], adj = c(0, 0), label = expression("p = 0.096"), srt = 45)

wilcox.test(freq ~ Transcription, data = nick_motif[nick_motif$pattern == "WWWSWWW", ])
add_significance(x0 = at[7], x1 = at[8], y = 1.025 * par("usr")[4], text_offset = -0.01, x_shrink = -0.15, text = NA)
text(x = 0.5 * sum(at[7:8]), y = 1.045 * par("usr")[4], adj = c(0, 0), label = expression("p = 0.48"), srt = 45)

wilcox.test(freq ~ Transcription, data = nick_motif[nick_motif$pattern == "WWWWSWWWW", ])
add_significance(x0 = at[9], x1 = at[10], y = 1.025 * par("usr")[4], text_offset = -0.01, x_shrink = -0.15, text = NA)
text(x = 0.5 * sum(at[9:10]), y = 1.045 * par("usr")[4], adj = c(0, 0), label = expression("p = 0.019"), srt = 45)

dev.off()
GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/Motif_freq.pdf"))


# Difference in meltability? ----------------------------------------------

# load calculated T_m values
load(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Data/Melting/Melting_10bp_v2.RData")

melt_transcribed <- subsetByIntersect(subject = melting_10bp, query = transcribed)
melt_untranscribed <- subsetByIntersect(subject = melting_10bp, query = untranscribed)

melting_grouped <- data.frame(Transcription = c(rep("Transcribed", length(melt_transcribed)), rep("Untranscribed", length(melt_untranscribed))),
                              T_m = c(melt_transcribed$T_m, melt_untranscribed$T_m))

# print to PDF ----------------------------------
pdf(file = "tmp.pdf", width=3, height=3)
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(3.7, 1, 2.4, 2), xpd = TRUE, las=1, tcl = -0.25, mgp = c(2, 0.5, 0))

beanplot(T_m ~ Transcription, data = melting_grouped, what = c(0, 1, 0, 0), col = gray(0.67), border = gray(level = 0.33), ylab = expression("Melting Temperature ["*degree*"C"*"]"))
boxplot(T_m ~ Transcription, data = melting_grouped, add = TRUE, outline = FALSE, col = NA, names = NA, axes = FALSE)

wilcox.test(T_m ~ Transcription, data = melting_grouped)

add_significance(x0 = 1, x1 = 2, y = 1.15 * max(melting_grouped$T_m), text_offset = -1, x_shrink = 0, text = expression("p < 10"^"-15"))

dev.off()
GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/Melting.pdf"))
