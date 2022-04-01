# info --------------------------------------------------------------------
# purpose: compare S1-seq of transcribed regions that are co-directional or converging with resection direction
# author: Robert Gnuegge (robert.gnuegge@gmail.com)
# date: 01/19/22
# version: 2.0
# note: co-directional, but not converging transcription will be inactivated by DSB formation or resection of promoter


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
source(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Analysis_scripts/Misc/JFly_colors.R")

plot_dir <- "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Figures/11_Transcription_impact/LSY4377-12B_4377-15A/Transcribed_vs_untranscribed"


# load S1-seq data --------------------------------------------------------

load(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Data/S1-seq/LSY4377-12B_4377-15A/LSY4377-12B_LSY4377-15A_around_DSBs.RData")

# only include DSBs that are within an ORF
DSBs <- SrfIcs[-c(9, 17)]  # exclude DSBs in duplicated genome region
# only keep regions with correct orientation w.r.t DSBs
roi <- DSB_regions(DSBs = DSBs, region_width = 5000, up_rev_down_fw = TRUE)

LSY4377_12B_1 <- subsetByIntersect(subject = LSY4377_12B_1_S1_seq, query = roi)
LSY4377_12B_2 <- subsetByIntersect(subject = LSY4377_12B_2_S1_seq, query = roi)
LSY4377_12B_4 <- subsetByIntersect(subject = LSY4377_12B_4_S1_seq, query = roi)


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


# plot t = 1 h ------------------------------------------------------------
tmp <- LSY4377_12B_1

hits <- findOverlaps(subject = RNA_seq, query = tmp, ignore.strand = TRUE)
tmp$RNA_seq <- mcols(RNA_seq[subjectHits(hits)])$score
tmp$transcribed <- tmp$RNA_seq > 0

mean(sort(unique(tmp[tmp$transcribed == TRUE]$DSB_kinetics_rank)))
# 10
mean(sort(unique(tmp[tmp$transcribed == FALSE]$DSB_kinetics_rank)))
# 10.625
# similar average DSB kinetics rank indicates that different nicking spreading is not due to different DSB kinetics

agg <- aggregate(score ~ distance_to_DSB + transcribed, data = tmp[tmp$DSB_kinetics_rank < 17], FUN = mean)
# exclude slowly formed DSBs

pdf(file = "tmp.pdf", width=3.5, height=3)
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(1.5, 0.6, 4, 1.6), tcl = -0.3, mgp = c(2.5, 0.6, 0), las = 1)
x <- agg$distance_to_DSB[agg$transcribed]
y <- moving_average(x = agg$score[agg$transcribed], k = k, keep = keep)
max(y)
max(y[-1:-5])
th <- 34
th_t <- 120
f <- 0.1
y <- trans(x = y, threshold = th, threshold_trans = th_t, factor = f)
plot(x = x, y = y, ylim = c(0, trans(x = 270, threshold = th, threshold_trans = th_t, factor = f)), 
     xlim = c(0, 1500), ylab = "Average S1-seq [RPM]", xlab = "Distance from DSB [nt]", 
     type = "l", lwd = 1.5, col = JFly_colors[1], yaxt = "n")

x <- agg$distance_to_DSB[!agg$transcribed]
y <- moving_average(x = agg$score[!agg$transcribed], k = k, keep = keep)
max(y)
max(y[-1:-5])
y <- trans(x = y, threshold = th, threshold_trans = th_t, factor = f)
points(x = x, y = y, type = "l", lwd = 1.5, col = JFly_colors[2])

below_break <- 0:6 * 5
above_break <- c(150, 200, 250)
axis(side = 2, at = c(below_break, trans(x = above_break, threshold = th, threshold_trans = th_t, factor = f)), labels = c(below_break, above_break))

axis.break(axis = 2, breakpos = th, style = "slash", bgcol = "white")
rect(xleft = -10, ybottom = th - 0.75, xright = 10, ytop =  th + 0.75, col = "white", border = "white")
dev.off()
GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/Transcribed_vs_untranscribed_1h.pdf"))


# plot t = 2 h ------------------------------------------------------------
tmp <- LSY4377_12B_2

hits <- findOverlaps(subject = RNA_seq, query = tmp, ignore.strand = TRUE)
tmp$RNA_seq <- mcols(RNA_seq[subjectHits(hits)])$score
tmp$transcribed <- tmp$RNA_seq > 0

agg <- aggregate(score ~ distance_to_DSB + transcribed, data = tmp[tmp$DSB_kinetics_rank < 17], FUN = mean)
# exclude slowly formed DSBs

pdf(file = "tmp.pdf", width=3.5, height=3)
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(1.5, 0.6, 4, 1.6), tcl = -0.3, mgp = c(2.5, 0.6, 0), las = 1)
x <- agg$distance_to_DSB[agg$transcribed]
y <- moving_average(x = agg$score[agg$transcribed], k = k, keep = keep)
plot(x = x, y = y,
     xlim = c(0, 1500), ylab = "Average S1-seq [RPM]", xlab = "Distance from DSB [nt]", 
     type = "l", lwd = 1.5, col = JFly_colors[1])

x <- agg$distance_to_DSB[!agg$transcribed]
y <- moving_average(x = agg$score[!agg$transcribed], k = k, keep = keep)
points(x = x, y = y, type = "l", lwd = 1.5, col = JFly_colors[2])
dev.off()
GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/Transcribed_vs_untranscribed_2h.pdf"))


# plot t = 4 h ------------------------------------------------------------
tmp <- LSY4377_12B_4

hits <- findOverlaps(subject = RNA_seq, query = tmp, ignore.strand = TRUE)
tmp$RNA_seq <- mcols(RNA_seq[subjectHits(hits)])$score
tmp$transcribed <- tmp$RNA_seq > 0

agg <- aggregate(score ~ distance_to_DSB + transcribed, data = tmp[tmp$DSB_kinetics_rank < 17], FUN = mean)
# exclude slowly formed DSBs

pdf(file = "tmp.pdf", width=3.5, height=3)
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(1.5, 0.6, 4, 1.6), tcl = -0.3, mgp = c(2.5, 0.6, 0), las = 1)
x <- agg$distance_to_DSB[agg$transcribed]
x[1:2] <- x[1:2] - 10
y <- moving_average(x = agg$score[agg$transcribed], k = k, keep = keep)
plot(x = x, y = y, ylim = c(0, 20), 
     xlim = c(0, 1500), ylab = "Average S1-seq [RPM]", xlab = "Distance from DSB [nt]", 
     type = "l", lwd = 1.5, col = JFly_colors[1])

x <- agg$distance_to_DSB[!agg$transcribed]
y <- moving_average(x = agg$score[!agg$transcribed], k = k, keep = keep)
points(x = x, y = y, type = "l", lwd = 1.5, col = JFly_colors[2])
dev.off()
GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/Transcribed_vs_untranscribed_4h.pdf"))


# plot legend -------------------------------------------------------------
pdf(file = "tmp.pdf", width=1.75, height=0.65)
par(cex = 1, mar = rep(0, 4))
plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
legend(x = 1, y = 1, col = JFly_colors[1:2], lwd = 1.5, legend = c("Transcribed", "Untranscribed"), xjust = 0.5, yjust = 0.5)

dev.off()
GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/Transcribed_vs_untranscribed_Legend.pdf"))


# determine codirectional vs. colliding ===================================

# make transcript GRanges -------------------------------------------------
# the RNA-seq data set does not contain strand information
# let's add it using strand of the corresponding gene

genes <- subsetByOverlaps(x = all_features, ranges = DSB_regions(DSBs = DSBs, region_width = 3000))

# remove unverified ORFs that overlap with verified ORFS
verified <- genes[genes$qualifier == "Verified"]
unverified <- genes[genes$qualifier != "Verified"]
hits <- findOverlaps(query = verified, subject = unverified, ignore.strand = TRUE)
genes <- genes[!genes$name %in% mcols(unverified[subjectHits(hits)])$name]
all_features <- all_features[!all_features$name %in% mcols(unverified[subjectHits(hits)])$name]

transcripts <- GRanges()  # initialize

for(n in 1:length(genes)){
  # find transcript UTR borders
  UTR <- GRanges(seqnames = seqnames(genes[n]), 
                 ranges = IRanges(start = mean(c(start(genes[n]), end(all_features[which(all_features$name == genes[n]$name) - 1]))),
                                  end = start(genes[n])))
  UTR_RNA_seq <- subsetByOverlaps(x = RNA_seq, ranges = UTR)
  adjusted_start <- start(UTR_RNA_seq[which.min(UTR_RNA_seq$score)])
  
  UTR <- GRanges(seqnames = seqnames(genes[n]), 
                 ranges = IRanges(start = end(genes[n]),
                                  end = mean(c(end(genes[n]), start(all_features[which(all_features$name == genes[n]$name) + 1])))))
  UTR_RNA_seq <- subsetByOverlaps(x = RNA_seq, ranges = UTR)
  adjusted_end <- end(UTR_RNA_seq[which.min(UTR_RNA_seq$score)])
  
  # construct GRange containing the transcript (ORF + UTRs) coordinates
  roi <- GRanges(seqnames = seqnames(genes[n]), ranges = IRanges(start = adjusted_start, end = adjusted_end))
  transcript <- subsetByOverlaps(x = RNA_seq, ranges = roi)
  strand(transcript) <- strand(genes[n])
  transcript$name <- genes[n]$name
  
  transcripts <- c(transcripts, transcript)  
}


# add codirectionality to GRanges -----------------------------------------

# add directionality (resection and transcription) to S1-seq GRanges object
# argument: GRanges, GRanges
# result: GRanges
add_directionality <- function(resection, transcription){
  out <- resection
  hits <- findOverlaps(query = transcription, subject = resection, ignore.strand = TRUE)
  resection$transcription_direction <- NA
  mcols(resection[subjectHits(hits)])$transcription_direction <- strand(transcription[queryHits(hits)])
  out$codirectional <- strand(resection) == mcols(resection)$transcription_direction
  return(out)
}

LSY4377_12B_1 <- add_directionality(resection = LSY4377_12B_1, transcription = transcripts)
LSY4377_12B_2 <- add_directionality(resection = LSY4377_12B_2, transcription = transcripts)
LSY4377_12B_4 <- add_directionality(resection = LSY4377_12B_4, transcription = transcripts)


mean(unique(LSY4377_12B_1[!is.na(LSY4377_12B_1$codirectional) & LSY4377_12B_1$codirectional == TRUE]$DSB_kinetics_rank))
# 10
mean(unique(LSY4377_12B_1[!is.na(LSY4377_12B_1$codirectional) & LSY4377_12B_1$codirectional == FALSE]$DSB_kinetics_rank))
# 10.11111
# similar average DSB kinetics rank indicates that different nicking spreading is not due to different DSB kinetics


# plot codirectional vs. colliding ========================================

k <- 51  # for moving average
keep <- 2


# plot t = 1 h ------------------------------------------------------------
tmp <- LSY4377_12B_1

agg <- aggregate(score ~ distance_to_DSB + codirectional, data = tmp, FUN = mean)
# exclude slowly formed DSBs

pdf(file = "tmp.pdf", width=3.5, height=3)
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(1.5, 0.6, 4, 1.6), tcl = -0.3, mgp = c(2.5, 0.6, 0), las = 1)
x <- agg$distance_to_DSB[!agg$codirectional]
y <- moving_average(x = agg$score[!agg$codirectional], k = k, keep = keep)
max(y)
max(y[-1:-5])
th <- 35
th_t <- 190
f <- 0.2
y <- trans(x = y, threshold = th, threshold_trans = th_t, factor = f)

plot(x = x, y = y, ylim = c(0, trans(x = 230, threshold = th, threshold_trans = th_t, factor = f)), 
     xlim = c(0, 1500), ylab = "Average S1-seq [RPM]", xlab = "Distance from DSB [nt]", 
     type = "l", lwd = 1.5, col = JFly_colors[1], yaxt = "n")

x <- agg$distance_to_DSB[agg$codirectional]
y <- moving_average(x = agg$score[agg$codirectional], k = k, keep = keep)
max(y)
max(y[-1:-5])
y <- trans(x = y, threshold = th, threshold_trans = th_t, factor = f)
points(x = x, y = y, type = "l", lwd = 1.5, col = JFly_colors[2])

below_break <- 0:3 * 10
above_break <- c(200, 220)
axis(side = 2, at = c(below_break, trans(x = above_break, threshold = th, threshold_trans = th_t, factor = f)), labels = c(below_break, above_break))

axis.break(axis = 2, breakpos = th, style = "slash", bgcol = "white")
rect(xleft = -10, ybottom = th - 0.5, xright = 10, ytop =  th + 0.5, col = "white", border = "white")
dev.off()
GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/Co-directional_vs_converging_1h.pdf"))


# plot t = 2 h ------------------------------------------------------------
tmp <- LSY4377_12B_2

agg <- aggregate(score ~ distance_to_DSB + codirectional, data = tmp, FUN = mean)
# exclude slowly formed DSBs

pdf(file = "tmp.pdf", width=3.5, height=3)
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(1.5, 0.6, 4, 1.6), tcl = -0.3, mgp = c(2.5, 0.6, 0), las = 1)
x <- agg$distance_to_DSB[!agg$codirectional]
y <- moving_average(x = agg$score[!agg$codirectional], k = k, keep = keep)
max(y)
max(y[-1:-5])
th <- 32
th_t <- 45
f <- 0.4
y <- trans(x = y, threshold = th, threshold_trans = th_t, factor = f)

plot(x = x, y = y, ylim = c(0, trans(x = 65, threshold = th, threshold_trans = th_t, factor = f)), 
     xlim = c(0, 1500), ylab = "Average S1-seq [RPM]", xlab = "Distance from DSB [nt]", 
     type = "l", lwd = 1.5, col = JFly_colors[1], yaxt = "n")

x <- agg$distance_to_DSB[agg$codirectional]
y <- moving_average(x = agg$score[agg$codirectional], k = k, keep = keep)
max(y)
max(y[-1:-5])
y <- trans(x = y, threshold = th, threshold_trans = th_t, factor = f)
points(x = x, y = y, type = "l", lwd = 1.5, col = JFly_colors[2])

below_break <- 0:6 * 5
above_break <- c(50, 60)
axis(side = 2, at = c(below_break, trans(x = above_break, threshold = th, threshold_trans = th_t, factor = f)), labels = c(below_break, above_break))

axis.break(axis = 2, breakpos = th, style = "slash", bgcol = "white")
rect(xleft = -10, ybottom = th - 0.5, xright = 10, ytop =  th + 0.5, col = "white", border = "white")
dev.off()
GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/Co-directional_vs_converging_2h.pdf"))


# plot t = 4 h ------------------------------------------------------------
tmp <- LSY4377_12B_4

agg <- aggregate(score ~ distance_to_DSB + codirectional, data = tmp, FUN = mean)
# exclude slowly formed DSBs

pdf(file = "tmp.pdf", width=3.5, height=3)
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(1.5, 0.6, 4, 1.6), tcl = -0.3, mgp = c(2.5, 0.6, 0), las = 1)
x <- agg$distance_to_DSB[!agg$codirectional]
y <- moving_average(x = agg$score[!agg$codirectional], k = k, keep = keep)
plot(x = x, y = y, 
     xlim = c(0, 1500), ylab = "Average S1-seq [RPM]", xlab = "Distance from DSB [nt]", 
     type = "l", lwd = 1.5, col = JFly_colors[1])

x <- agg$distance_to_DSB[agg$codirectional]
y <- moving_average(x = agg$score[agg$codirectional], k = k, keep = keep)
points(x = x, y = y, type = "l", lwd = 1.5, col = JFly_colors[2])
dev.off()
GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/Co-directional_vs_converging_4h.pdf"))


# plot codirectinal vs. colliding for (un)transcribed =====================

# plot t = 1 h ------------------------------------------------------------
tmp <- LSY4377_12B_1

hits <- findOverlaps(subject = RNA_seq, query = tmp, ignore.strand = TRUE)
tmp$RNA_seq <- mcols(RNA_seq[subjectHits(hits)])$score
tmp$transcribed <- tmp$RNA_seq > 0

agg <- aggregate(score ~ distance_to_DSB + codirectional + transcribed, data = tmp, FUN = mean)
# exclude slowly formed DSBs

# transcribed -------------------------------------------------------------
pdf(file = "tmp.pdf", width=3.5, height=3)
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(1.5, 0.6, 4, 1.6), tcl = -0.3, mgp = c(2.5, 0.6, 0), las = 1)
x <- agg$distance_to_DSB[!agg$codirectional & agg$transcribed]
y <- moving_average(x = agg$score[!agg$codirectional & agg$transcribed], k = k, keep = keep)
max(y)
max(y[-1:-5])
th <- 34
th_t <- 205
f <- 0.33
y <- trans(x = y, threshold = th, threshold_trans = th_t, factor = f)

plot(x = x, y = y, ylim = c(0, trans(x = 225, threshold = th, threshold_trans = th_t, factor = f)), 
     xlim = c(0, 1500), ylab = "Average S1-seq [RPM]", xlab = "Distance from DSB [nt]", 
     type = "l", lwd = 1.5, col = JFly_colors[1], yaxt = "n")

x <- agg$distance_to_DSB[agg$codirectional & agg$transcribed]
y <- moving_average(x = agg$score[agg$codirectional & agg$transcribed], k = k, keep = keep)
max(y)
max(y[-1:-5])
y <- trans(x = y, threshold = th, threshold_trans = th_t, factor = f)
points(x = x, y = y, type = "l", lwd = 1.5, col = JFly_colors[2])

below_break <- 0:6 * 5
above_break <- c(210, 220)
axis(side = 2, at = c(below_break, trans(x = above_break, threshold = th, threshold_trans = th_t, factor = f)), labels = c(below_break, above_break))

axis.break(axis = 2, breakpos = th, style = "slash", bgcol = "white")
rect(xleft = -10, ybottom = th - 0.5, xright = 10, ytop =  th + 0.5, col = "white", border = "white")
dev.off()
GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/Co-directional_vs_converging_1h_transcribed.pdf"))

# untranscribed -----------------------------------------------------------
pdf(file = "tmp.pdf", width=3.5, height=3)
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(1.5, 0.6, 4, 1.6), tcl = -0.3, mgp = c(2.5, 0.6, 0), las = 1)
x <- agg$distance_to_DSB[!agg$codirectional & !agg$transcribed]
y <- moving_average(x = agg$score[!agg$codirectional & !agg$transcribed], k = k, keep = keep)
max(y)
max(y[-1:-5])
th <- 27
th_t <- 25
f <- 0.03
y <- trans(x = y, threshold = th, threshold_trans = th_t, factor = f)

plot(x = x, y = y, ylim = c(0, trans(x = 300, threshold = th, threshold_trans = th_t, factor = f)), 
     xlim = c(0, 1500), ylab = "Average S1-seq [RPM]", xlab = "Distance from DSB [nt]", 
     type = "l", lwd = 1.5, col = JFly_colors[1], yaxt = "n")

x <- agg$distance_to_DSB[agg$codirectional & !agg$transcribed]
y <- moving_average(x = agg$score[agg$codirectional & !agg$transcribed], k = k, keep = keep)
max(y)
max(y[-1:-5])
y <- trans(x = y, threshold = th, threshold_trans = th_t, factor = f)
points(x = x, y = y, type = "l", lwd = 1.5, col = JFly_colors[2])

below_break <- 0:5 * 5
above_break <- c(100, 200, 300)
axis(side = 2, at = c(below_break, trans(x = above_break, threshold = th, threshold_trans = th_t, factor = f)), labels = c(below_break, above_break))

axis.break(axis = 2, breakpos = th, style = "slash", bgcol = "white")
rect(xleft = -10, ybottom = th - 0.75, xright = 10, ytop =  th + 0.75, col = "white", border = "white")
dev.off()
GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/Co-directional_vs_converging_1h_untranscribed.pdf"))


# plot t = 2 h ------------------------------------------------------------
tmp <- LSY4377_12B_2

hits <- findOverlaps(subject = RNA_seq, query = tmp, ignore.strand = TRUE)
tmp$RNA_seq <- mcols(RNA_seq[subjectHits(hits)])$score
tmp$transcribed <- tmp$RNA_seq > 0

agg <- aggregate(score ~ distance_to_DSB + codirectional + transcribed, data = tmp, FUN = mean)
# exclude slowly formed DSBs

# transcribed -------------------------------------------------------------
pdf(file = "tmp.pdf", width=3.5, height=3)
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(1.5, 0.6, 4, 1.6), tcl = -0.3, mgp = c(2.5, 0.6, 0), las = 1)
x <- agg$distance_to_DSB[!agg$codirectional & agg$transcribed]
y <- moving_average(x = agg$score[!agg$codirectional & agg$transcribed], k = k, keep = keep)
max(y)
max(y[-1:-5])
th <- 29
th_t <- 45
f <- 0.3
y <- trans(x = y, threshold = th, threshold_trans = th_t, factor = f)

plot(x = x, y = y, ylim = c(0, trans(x = 65, threshold = th, threshold_trans = th_t, factor = f)), 
     xlim = c(0, 1500), ylab = "Average S1-seq [RPM]", xlab = "Distance from DSB [nt]", 
     type = "l", lwd = 1.5, col = JFly_colors[1], yaxt = "n")

x <- agg$distance_to_DSB[agg$codirectional & agg$transcribed]
y <- moving_average(x = agg$score[agg$codirectional & agg$transcribed], k = k, keep = keep)
max(y)
max(y[-1:-5])
y <- trans(x = y, threshold = th, threshold_trans = th_t, factor = f)
points(x = x, y = y, type = "l", lwd = 1.5, col = JFly_colors[2])

below_break <- 0:5 * 5
above_break <- c(50, 60)
axis(side = 2, at = c(below_break, trans(x = above_break, threshold = th, threshold_trans = th_t, factor = f)), labels = c(below_break, above_break))

axis.break(axis = 2, breakpos = th, style = "slash", bgcol = "white")
rect(xleft = -10, ybottom = th - 0.5, xright = 10, ytop =  th + 0.5, col = "white", border = "white")
dev.off()
GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/Co-directional_vs_converging_2h_transcribed.pdf"))

# untranscribed -----------------------------------------------------------
pdf(file = "tmp.pdf", width=3.5, height=3)
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(1.5, 0.6, 4, 1.6), tcl = -0.3, mgp = c(2.5, 0.6, 0), las = 1)
x <- agg$distance_to_DSB[!agg$codirectional & !agg$transcribed]
y <- moving_average(x = agg$score[!agg$codirectional & !agg$transcribed], k = k, keep = keep)
max(y)
max(y[-1:-5])

plot(x = x, y = y, ylim = c(0, 57), 
     xlim = c(0, 1500), ylab = "Average S1-seq [RPM]", xlab = "Distance from DSB [nt]", 
     type = "l", lwd = 1.5, col = JFly_colors[1])

x <- agg$distance_to_DSB[agg$codirectional & !agg$transcribed]
y <- moving_average(x = agg$score[agg$codirectional & !agg$transcribed], k = k, keep = keep)
max(y)
max(y[-1:-5])
points(x = x, y = y, type = "l", lwd = 1.5, col = JFly_colors[2])
dev.off()
GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/Co-directional_vs_converging_2h_untranscribed.pdf"))


# plot t = 4 h ------------------------------------------------------------
tmp <- LSY4377_12B_4

hits <- findOverlaps(subject = RNA_seq, query = tmp, ignore.strand = TRUE)
tmp$RNA_seq <- mcols(RNA_seq[subjectHits(hits)])$score
tmp$transcribed <- tmp$RNA_seq > 0

agg <- aggregate(score ~ distance_to_DSB + codirectional + transcribed, data = tmp, FUN = mean)

# transcribed -------------------------------------------------------------
pdf(file = "tmp.pdf", width=3.5, height=3)
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(1.5, 0.6, 4, 1.6), tcl = -0.3, mgp = c(2.5, 0.6, 0), las = 1)
x <- agg$distance_to_DSB[!agg$codirectional & agg$transcribed]
y <- moving_average(x = agg$score[!agg$codirectional & agg$transcribed], k = k, keep = keep)
plot(x = x, y = y, 
     xlim = c(0, 1500), ylab = "Average S1-seq [RPM]", xlab = "Distance from DSB [nt]", 
     type = "l", lwd = 1.5, col = JFly_colors[1])

x <- agg$distance_to_DSB[agg$codirectional & agg$transcribed]
y <- moving_average(x = agg$score[agg$codirectional & agg$transcribed], k = k, keep = keep)
points(x = x, y = y, type = "l", lwd = 1.5, col = JFly_colors[2])
dev.off()
GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/Co-directional_vs_converging_4h_transcribed.pdf"))

# untranscribed -----------------------------------------------------------
pdf(file = "tmp.pdf", width=3.5, height=3)
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(1.5, 0.6, 4, 1.6), tcl = -0.3, mgp = c(2.5, 0.6, 0), las = 1)
x <- agg$distance_to_DSB[!agg$codirectional & !agg$transcribed]
y <- moving_average(x = agg$score[!agg$codirectional & !agg$transcribed], k = k, keep = keep)
plot(x = x, y = y,
     xlim = c(0, 1500), ylab = "Average S1-seq [RPM]", xlab = "Distance from DSB [nt]", 
     type = "l", lwd = 1.5, col = JFly_colors[1])

x <- agg$distance_to_DSB[agg$codirectional & !agg$transcribed]
y <- moving_average(x = agg$score[agg$codirectional & !agg$transcribed], k = k, keep = keep)
max(y)
points(x = x, y = y, type = "l", lwd = 1.5, col = JFly_colors[2])
dev.off()
GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/Co-directional_vs_converging_4h_untranscribed.pdf"))

# plot legend -------------------------------------------------------------
pdf(file = "tmp.pdf", width=1.8, height=0.65)
par(cex = 1, mar = rep(0, 4))
plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
legend(x = 1, y = 1, col = JFly_colors[1:2], lwd = 1.5, legend = c("Converging", "Co-directional"), xjust = 0.5, yjust = 0.5)

dev.off()
GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/Co-directional_vs_converging_Legend.pdf"))