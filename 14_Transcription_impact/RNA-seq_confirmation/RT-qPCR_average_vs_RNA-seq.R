# info --------------------------------------------------------------------
# purpose: calculate and plot average RT-qPCR levels and plot vs. RNA-seq
# author: Robert Gnuegge (robert.gnuegge@gmail.com)
# date: 08/16/22
# version: 1.0


# preamble ----------------------------------------------------------------

# set wd to this file's location
wd_path <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(wd_path)

# load libraries
library(propagate)
library(GenomicRanges)
library(rtracklayer)

# set up directory paths
plot_dir <- "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Figures/11_Transcription_impact/RNA-seq_confirmation/RT-qPCR_vs_RNA-seq/"
data_dir <- "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Data/RNA-seq_confirmation/Replicates"

# read in helper files and functions
source(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Analysis_scripts/Misc/Genomic_helper_functions.R")
source(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Analysis_scripts/Misc/Misc_helper_functions.R")
source(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Analysis_scripts/Misc/JFly_colors.R")


# collect and process RT-qPCR data ========================================

# read in data
files <- grep(pattern = ".txt", x = list.files(path = data_dir, full.names = TRUE), value = TRUE)

RT_qPCR_reps <- data.frame()
for(file in files){
  tmp <- read.table(file = file, header = TRUE)
  RT_qPCR_reps <- rbind(RT_qPCR_reps, tmp)
}

RT_qPCR_reps <- RT_qPCR_reps[RT_qPCR_reps$sample == "+RT", ]

# calc Cq mean, sd, and cv
RT_qPCR <- aggregate(mean ~ gene, data = RT_qPCR_reps, FUN = function(x) c(mean = mean(x), sd = sd(x)))

# mean and sd are in a single column (as a matrix)
# convert into separate columns
RT_qPCR <- data.frame(RT_qPCR[, 1:(ncol(RT_qPCR)-1)], as.data.frame(RT_qPCR$mean))
colnames(RT_qPCR) <- c("gene", "mean", "sd")


# plot RT-qPCR vs RNA-seq =================================================

# read RNA-seq data
RNA_seq_1 <- import("/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Resources/Maya-Miles2019/GSM3567364_w303_rep1.bigwig")
RNA_seq_2 <- import("/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Resources/Maya-Miles2019/GSM3567365_w303_rep2.bigwig")

# find ORFs of RT-qPCR genes in RNA-seq data sets
load(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Data/S_cerevisiae_chromosomal_features/S_cerevisiae_genome_features.RData")
ORFs <- all_features[all_features$name %in% RT_qPCR$gene]

# get RNA-seq score averaged over each ORF for both RNA-seq replicates
RNA_seq_score <- data.frame(rep1 = rep(0, length(ORFs)), rep2 = rep(0, length(ORFs)))

for(n in 1:nrow(RNA_seq_score)){
  RNA_seq_score$rep1[n] <- mean(mcols(subsetByIntersect(subject = RNA_seq_1, query = ORFs[n]))$score)
  RNA_seq_score$rep2[n] <- mean(mcols(subsetByIntersect(subject = RNA_seq_2, query = ORFs[n]))$score)
}

RNA_seq_score$mean <- apply(X = RNA_seq_score, MARGIN = 1, FUN = mean)
RNA_seq_score$sd <- apply(X = RNA_seq_score, MARGIN = 1, FUN = sd)
mcols(ORFs) <- cbind(mcols(ORFs), mean = RNA_seq_score$mean, sd = RNA_seq_score$sd)

ORFs$rel_mean <- ORFs$mean / ORFs$mean[ORFs$name == "ADH1"]
ORFs$rel_sd <- ORFs$rel_mean * ORFs$sd / ORFs$mean

# apply same gene order to both data sets
RNA_seq <- ORFs[order(ORFs$rel_mean)]  # order with increasing RNA-seq score
RT_qPCR <- RT_qPCR[match(x = RNA_seq$name, table = RT_qPCR$gene), ]  # apply same order to RT-qPCR data
all(RNA_seq$name == RT_qPCR$gene)


# print to PDF ------------------------------------------------------------
pdf(file = "tmp.pdf", width=3.75, height=3.5)
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(2.1, -0.1, 4, 2), tcl = -0.25, mgp = c(2, 0.5, 0), las = 1)

# plot means
plot(x = RNA_seq$rel_mean, y = RT_qPCR$mean, log = "xy", pch = 20, xlim = c(0.0005, 1.5), ylim = c(0.0005, 1.5),
     xlab = "RNA-seq", ylab = NA, xaxt = "n", yaxt = "n")
axis(side = 1, at = 10^(-3:0))
axis(side = 2, at = 10^(-3:0))
title(ylab = "RT-qPCR", line = 3.25)

# add RT-qPCR error bars
arrows(x0 = RNA_seq$rel_mean,
       y0 = RT_qPCR$mean - RT_qPCR$sd,
       x1 = RNA_seq$rel_mean,
       y1 = RT_qPCR$mean + RT_qPCR$sd,
       length = 0.03, # length of arrow head
       angle = 90, # angle of arrow head
       code = 3, # to draw arrow head on both ends
       col = "black")

# add linear regression line
lin_reg <- lm(RT_qPCR$mean ~ RNA_seq$rel_mean, weights = 1 / RT_qPCR$sd^2)
abline(lm(RT_qPCR$mean ~ RNA_seq$rel_mean), lty = "dashed", col = "gray")
# save regression results in file
sink(file = paste0(plot_dir, "Linear_regression_result.txt"))
summary(lin_reg)
sink()

# add gene names
adj_x <- c(0.7, 0.3, 1.2, 1.2, 1.2, 1.2, 1.2, -0.05, 0.3, 1.1, 0.9)
adj_y <- c(1.5, 1.7, 0.5, 0.5, 0.5, 0.2, 0.5, 0.8, -0.6, 0.3, -0.5)
for(n in 1:11){
  text(x =  RNA_seq$rel_mean[n], y =  RT_qPCR$mean[n], adj = c(adj_x[n],adj_y[n]), labels = RNA_seq$name[n], col = "gray", font = 3)
}

dev.off()
GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/Average_RT-qPCR_vs_RNA-seq.pdf"))
