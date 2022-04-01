# info --------------------------------------------------------------------
# purpose: compare S1-seq up- and downstream of DSBs that interrupt a transcribed gene
# author: Robert Gnuegge (robert.gnuegge@gmail.com)
# date: 01/26/22
# version: 1.0
# note: DSB "shuts off" transcription on the promoter-distal site of the break,
# this creates a transcribed and an untranscribed region


# preamble ----------------------------------------------------------------

# set wd to this file's location
wd_path <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(wd_path)

# libraries
library(GenomicRanges)
library(rtracklayer)
library(plotrix)

# load helper files and functions
source(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Analysis_scripts/Misc/Genomic_helper_functions.R")
source(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Analysis_scripts/Misc/S_cerevisiae_SrfI_cut_sites.R")
source(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Analysis_scripts/Misc/Misc_helper_functions.R")

plot_dir <- "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Figures/11_Transcription_impact/LSY4377-12B_4377-15A_merged/Up_vs_downstream_of_DSB"


# load S1-seq data --------------------------------------------------------

load(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Data/S1-seq/LSY4377-12B_4377-15A_merged/LSY4377-12B_LSY4377-15A_merged_around_DSBs.RData")

# only include DSBs that are within an ORF
DSBs <- SrfIcs[c(3, 5, 8, 12, 13, 15, 18, 19, 21)]
# only keep regions with correct orientation w.r.t DSBs
roi <- DSB_regions(DSBs = DSBs, region_width = 5000, up_rev_down_fw = TRUE)

LSY4377_12B_1 <- subsetByIntersect(subject = LSY4377_12B_1_merged_S1_seq, query = roi)
LSY4377_12B_2 <- subsetByIntersect(subject = LSY4377_12B_2_merged_S1_seq, query = roi)
LSY4377_12B_4 <- subsetByIntersect(subject = LSY4377_12B_4_merged_S1_seq, query = roi)


# load RNA-seq data -------------------------------------------------------
# data are from Maya-Miles et al, 2019 (pmid: 31331360)

# there are two biological replicates
roi <- DSB_regions(DSBs = DSBs, region_width = 5000)
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


# make transcript GRanges -------------------------------------------------

transcripts <- GRanges()  # initialize

for(n in 1:length(DSBs)){
  # find cut gene
  gene <- subsetByOverlaps(x = all_features, ranges = DSBs[n])
  
  # find transcript UTR borders
  UTR <- GRanges(seqnames = seqnames(gene), 
                 ranges = IRanges(start = mean(c(start(gene), end(all_features[which(all_features$name == gene$name) - 1]))),
                                  end = start(gene)))
  UTR_RNA_seq <- subsetByOverlaps(x = RNA_seq, ranges = UTR)
  adjusted_start <- start(UTR_RNA_seq[which.min(UTR_RNA_seq$score)])
  
  UTR <- GRanges(seqnames = seqnames(gene), 
                 ranges = IRanges(start = end(gene),
                                  end = mean(c(end(gene), start(all_features[which(all_features$name == gene$name) + 1])))))
  UTR_RNA_seq <- subsetByOverlaps(x = RNA_seq, ranges = UTR)
  adjusted_end <- end(UTR_RNA_seq[which.min(UTR_RNA_seq$score)])
  
  # construct GRange containing the transcript (ORF + UTRs) coordinates
  roi <- GRanges(seqnames = seqnames(gene), ranges = IRanges(start = adjusted_start, end = adjusted_end))
  transcript <- subsetByOverlaps(x = RNA_seq, ranges = roi)
  strand(transcript) <- strand(gene)
  transcript$name <- gene$name

  transcripts <- c(transcripts, transcript)  
}


# plotting ================================================================

my_colors <- gray(level = c(0, 0.7))
k <- 51  # for moving average
keep <- 2

roi <- transcripts
strand(roi) <- "*"


# plot t = 1 h ------------------------------------------------------------

S1_seq <- subsetByIntersect(subject = LSY4377_12B_1, query = roi)
S1_seq$prom_distal <- strand(S1_seq) == strand(as_nt_resolved_GRanges(transcripts))  # co-directional = promoter-distal

agg <- aggregate(score ~ distance_to_DSB + prom_distal, data = S1_seq, FUN = mean)

pdf(file = "tmp.pdf", width=2.75, height=2.5)
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(1.8, 0.9, 4, 1.5), tcl = -0.25, mgp = c(2.25, 0.5, 0), las = 1)
  x <- agg$distance_to_DSB[!agg$prom_distal]
  y <- moving_average(x = agg$score[!agg$prom_distal], k = k, keep = keep)
  max(y)
  max(y[-1:-5])
  th <- 13
  th_t <- 310
  f <- 0.075
  y <- trans(x = y, threshold = th, threshold_trans = th_t, factor = f)
  plot(x = x, y = y, ylim = c(0, trans(x = 360, threshold = th, threshold_trans = th_t, factor = f)), 
       xlim = c(0, max(agg$distance_to_DSB[agg$prom_distal])), ylab = "Average S1-seq [RPM]", xlab = "Distance from DSB [nt]", 
       type = "l", lwd = 1.5, col = my_colors[1], yaxt = "n", xaxt = "n")
  axis(side = 1, at = c(0, 2, 4, 6) * 100)
  
  x <- agg$distance_to_DSB[agg$prom_distal]
  y <- moving_average(x = agg$score[agg$prom_distal], k = k, keep = keep)
  max(y)
  max(y[-1:-5])
  y <- trans(x = y, threshold = th, threshold_trans = th_t, factor = f)
  points(x = x, y = y, type = "l", lwd = 1.5, col = my_colors[2])
  
  below_break <- 0:6 * 2
  above_break <- c(320, 340, 360)
  axis(side = 2, at = c(below_break, trans(x = above_break, threshold = th, threshold_trans = th_t, factor = f)), labels = c(below_break, above_break))
  
  axis.break(axis = 2, breakpos = th, style = "slash", bgcol = "white")
  rect(xleft = -10, ybottom = th - 0.2, xright = 10, ytop =  th + 0.2, col = "white", border = "white")
  
  legend(x = "topright", legend = c("Promoter-proximal", "Promoter-distal"), col = my_colors, lwd = 1.5, seg.len = 1, bty = "n", inset = -0.02)
  
dev.off()
GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/Avg_S1-seq_spreading_1h.pdf"))


# plot t = 2 h ------------------------------------------------------------

S1_seq <- subsetByIntersect(subject = LSY4377_12B_2, query = roi)
S1_seq$prom_distal <- strand(S1_seq) == strand(as_nt_resolved_GRanges(transcripts))  # co-directional = promoter-distal

agg <- aggregate(score ~ distance_to_DSB + prom_distal, data = S1_seq, FUN = mean)

pdf(file = "tmp.pdf", width=2.75, height=2.5)
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(1.8, 0.9, 4, 1.5), tcl = -0.25, mgp = c(2.25, 0.5, 0), las = 1)
  x <- agg$distance_to_DSB[!agg$prom_distal]
  y <- moving_average(x = agg$score[!agg$prom_distal], k = k, keep = keep)
  max(y)
  max(y[-1:-5])
  th <- 23
  th_t <- 85
  f <- 0.15
  y <- trans(x = y, threshold = th, threshold_trans = th_t, factor = f)
  plot(x = x, y = y, ylim = c(0, trans(x = 125, threshold = th, threshold_trans = th_t, factor = f)), 
       xlim = c(0, max(agg$distance_to_DSB[agg$prom_distal])), ylab = "Average S1-seq [RPM]", xlab = "Distance from DSB [nt]", 
       type = "l", lwd = 1.5, col = my_colors[1], yaxt = "n", xaxt = "n")
  axis(side = 1, at = c(0, 2, 4, 6) * 100)
  
  x <- agg$distance_to_DSB[agg$prom_distal]
  y <- moving_average(x = agg$score[agg$prom_distal], k = k, keep = keep)
  max(y)
  max(y[-1:-5])
  y <- trans(x = y, threshold = th, threshold_trans = th_t, factor = f)
  points(x = x, y = y, type = "l", lwd = 1.5, col = my_colors[2])
  
  below_break <- 0:4 * 5
  above_break <- c(100, 120)
  axis(side = 2, at = c(below_break, trans(x = above_break, threshold = th, threshold_trans = th_t, factor = f)), labels = c(below_break, above_break))
  
  axis.break(axis = 2, breakpos = th, style = "slash", bgcol = "white")
  rect(xleft = -10, ybottom = th - 0.5, xright = 10, ytop =  th + 0.5, col = "white", border = "white")
dev.off()
GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/Avg_S1-seq_spreading_2h.pdf"))


# plot t = 4 h ------------------------------------------------------------

S1_seq <- subsetByIntersect(subject = LSY4377_12B_4, query = roi)
S1_seq$prom_distal <- strand(S1_seq) == strand(as_nt_resolved_GRanges(transcripts))  # co-directional = promoter-distal

agg <- aggregate(score ~ distance_to_DSB + prom_distal, data = S1_seq, FUN = mean)

pdf(file = "tmp.pdf", width=2.75, height=2.5)
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(1.8, 0.9, 4, 1.5), tcl = -0.25, mgp = c(2.25, 0.5, 0), las = 1)
  x <- agg$distance_to_DSB[!agg$prom_distal]
  y <- moving_average(x = agg$score[!agg$prom_distal], k = k, keep = keep)
  # max(y)
  # max(y[-1:-5])
  # th <- 22.5
  # th_t <- 37.5
  # f <- 1
  # y <- trans(x = y, threshold = th, threshold_trans = th_t, factor = f)
  plot(x = x, y = y, # ylim = c(0, trans(x = 45, threshold = th, threshold_trans = th_t, factor = f)), 
       xlim = c(0, max(agg$distance_to_DSB[agg$prom_distal])), ylab = "Average S1-seq [RPM]", xlab = "Distance from DSB [nt]", 
       type = "l", lwd = 1.5, col = my_colors[1], xaxt = "n")
  axis(side = 1, at = c(0, 2, 4, 6) * 100)
  
  x <- agg$distance_to_DSB[agg$prom_distal]
  y <- moving_average(x = agg$score[agg$prom_distal], k = k, keep = keep)
  max(y)
  max(y[-1:-5])
  # y <- trans(x = y, threshold = th, threshold_trans = th_t, factor = f)
  points(x = x, y = y, type = "l", lwd = 1.5, col = my_colors[2])
  
  # below_break <- c(0, 5, 10, 15, 20)
  # above_break <- c(40, 45)
  # axis(side = 2, at = c(below_break, trans(x = above_break, threshold = th, threshold_trans = th_t, factor = f)), labels = c(below_break, above_break))
  # axis.break(axis = 2, breakpos = th, style = "slash", bgcol = "white")
  # rect(xleft = -10, ybottom = th - 0.5, xright = 10, ytop =  th + 0.5, col = "white", border = "white")
dev.off()
GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/Avg_S1-seq_spreading_4h.pdf"))
