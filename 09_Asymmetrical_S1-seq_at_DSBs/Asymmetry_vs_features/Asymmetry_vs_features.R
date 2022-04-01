# info --------------------------------------------------------------------
# purpose: check how much of the observed nicking can be predicted based on the identified impacting features
# author: Robert Gnuegge (robert.gnuegge@gmail.com)
# date: 01/31/22
# version: 2.0


# preamble ----------------------------------------------------------------

# set wd to this file's location
wd_path <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(wd_path)

# libraries
library(GenomicRanges)
library(BSgenome)
library(rtracklayer)

# load helper files and functions
source(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Analysis_scripts/Misc/Genomic_helper_functions.R")
source(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Analysis_scripts/Misc/S_cerevisiae_SrfI_cut_sites.R")
source(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Analysis_scripts/Misc/Misc_helper_functions.R")
source(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Analysis_scripts/Misc/JFly_colors.R")

plot_dir <- "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Figures/06_Asymmetry/Asymmetry_vs_features"


# function definition -----------------------------------------------------

# add RNA-seq to S1-seq data sets
# arguments: GRanges
# result: GRanges
add_RNA_seq <- function(S1_seq, RNA_seq){
  hits <- findOverlaps(query = RNA_seq, subject = S1_seq)
  S1_seq$RNA_seq[subjectHits(hits)] <- RNA_seq$score[queryHits(hits)]
  return(S1_seq)
}

# add directionality (resection and transcription) to S1-seq GRanges object
# arguments: GRanges
# result: GRanges
add_directionality <- function(S1_seq, transcripts){
  hits <- findOverlaps(query = transcripts, subject = S1_seq, ignore.strand = TRUE)
  S1_seq$co_direction[subjectHits(hits)] <- as.character(strand(S1_seq[subjectHits(hits)])) == as.character(strand(transcripts[queryHits(hits)]))
  return(S1_seq)
}

# add MNase-seq to S1-seq data sets
# arguments: GRanges
# result: GRanges
add_MNase_seq <- function(S1_seq, MNase_seq){
  hits <- findOverlaps(query = MNase_seq, subject = S1_seq)
  S1_seq$MNase_seq[subjectHits(hits)] <- MNase_seq$score[queryHits(hits)]
  return(S1_seq)
}

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

# calculate position weight matrices
# argument: GRanges, even integer, genome (DNAStringSet)
# result: matrix
calc_PWM_matching_score <- function(GRanges, PWM, genome){
  nick_pos <- granges(GRanges)
  width <- ncol(PWM)
  nick_context <- flank(x = nick_pos, width = 0.5 * width, both = TRUE)
  nt_sequences <- getSeq(x = genome, names = nick_context)
  nts <- strsplit(as.character(nt_sequences), "")
  sapply(X = nts, FUN = function(x){sum(diag(PWM[x, ]))})
}

# calculate average meltability profiles
# argument: list with GRanges and matrix (from get_meltability_profiles function)
# result: numeric vector
average_profile <- function(GRanges_and_profiles, GRanges, bg = FALSE){
  hits <- findOverlaps(query = GRanges, subject = GRanges_and_profiles$GRanges)
  if(bg){
    out <- apply(X = GRanges_and_profiles$meltability_profiles[subjectHits(hits), ], MARGIN = 2, FUN = mean)
  }else{
    out <- apply(X = GRanges_and_profiles$meltability_profiles[subjectHits(hits), ] * GRanges$score[queryHits(hits)], MARGIN = 2, FUN = sum) / sum(GRanges$score[queryHits(hits)])
  }
}


# load S1-seq data --------------------------------------------------------

load(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Data/S1-seq/LSY4377-12B_4377-15A_merged/LSY4377-12B_LSY4377-15A_merged_around_DSBs.RData")

# only include DSBs that are within an ORF
DSBs <- SrfIcs[-c(9, 17, 14, 15)]  # exclude DSBs in duplicated genome region and on chrXIII with low distance between DSBs
# only keep regions with correct orientation w.r.t DSBs
roi <- DSB_regions(DSBs = DSBs, region_width = 4000, up_rev_down_fw = TRUE)

LSY4377_12B_1 <- subsetByIntersect(subject = LSY4377_12B_1_merged_S1_seq, query = roi)
LSY4377_12B_2 <- subsetByIntersect(subject = LSY4377_12B_2_merged_S1_seq, query = roi)
LSY4377_12B_4 <- subsetByIntersect(subject = LSY4377_12B_4_merged_S1_seq, query = roi)


SrfIcs_df

# add RNA-seq data --------------------------------------------------------
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

LSY4377_12B_1 <- add_RNA_seq(S1_seq = LSY4377_12B_1, RNA_seq = RNA_seq)
LSY4377_12B_2 <- add_RNA_seq(S1_seq = LSY4377_12B_2, RNA_seq = RNA_seq)
LSY4377_12B_4 <- add_RNA_seq(S1_seq = LSY4377_12B_4, RNA_seq = RNA_seq)


# add codirectional vs. colliding -----------------------------------------
# the RNA-seq data set does not contain strand information
# let's add it using strand of the corresponding gene

# load genome feature
# data downloaded from yeastmine.yeastgenome.org and curated with Adjust_chromosomal_features.R
load(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Data/S_cerevisiae_chromosomal_features/S_cerevisiae_genome_features.RData")
all_features <- all_features[all_features$type == "ORF"]

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

LSY4377_12B_1 <- add_directionality(S1_seq = LSY4377_12B_1, transcripts = transcripts)
LSY4377_12B_2 <- add_directionality(S1_seq = LSY4377_12B_2, transcripts = transcripts)
LSY4377_12B_4 <- add_directionality(S1_seq = LSY4377_12B_4, transcripts = transcripts)


# add MNase-seq data ------------------------------------------------------
load(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Data/MNase-seq/LSY4377-12B_4377-15A_merged/LSY4377-12B_LSY4377-15A_merged_trimmed.RData")

LSY4377_12B_1 <- add_MNase_seq(S1_seq = LSY4377_12B_1, MNase_seq = LSY4377_12B_0_merged_MNase_seq)
LSY4377_12B_2 <- add_MNase_seq(S1_seq = LSY4377_12B_2, MNase_seq = LSY4377_12B_0_merged_MNase_seq)
LSY4377_12B_4 <- add_MNase_seq(S1_seq = LSY4377_12B_4, MNase_seq = LSY4377_12B_0_merged_MNase_seq)


# add similarity to preferred nt sequence ---------------------------------

# read modified S. cerevisiae genome
genome <- import(con = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Data/S_cerevisiae_reference_genomes/S288C_R64-2-1_W303_SNPs_MATa_hocs2SrfIcs_hml_hmr.fasta")

# calc PWM
PWM_all <- calc_PWM(GRanges = c(LSY4377_12B_1, LSY4377_12B_2, LSY4377_12B_4), width = 10, genome = genome)

LSY4377_12B_1$nt_seq_score <- calc_PWM_matching_score(GRanges = LSY4377_12B_1, PWM = PWM_all, genome = genome)
LSY4377_12B_2$nt_seq_score <- calc_PWM_matching_score(GRanges = LSY4377_12B_2, PWM = PWM_all, genome = genome)
LSY4377_12B_4$nt_seq_score <- calc_PWM_matching_score(GRanges = LSY4377_12B_4, PWM = PWM_all, genome = genome)


# add similarity to preferred melting profile -----------------------------
load(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Data/Melting/Extracted_melting_profiles.RData")

LSY4377_12B_1 <- subsetByIntersect(subject = LSY4377_12B_1, query = LSY4377_12B_4_10bp_melting_profiles$GRanges)
LSY4377_12B_2 <- subsetByIntersect(subject = LSY4377_12B_2, query = LSY4377_12B_4_10bp_melting_profiles$GRanges)
LSY4377_12B_4 <- subsetByIntersect(subject = LSY4377_12B_4, query = LSY4377_12B_4_10bp_melting_profiles$GRanges)

melt_profile <- average_profile(GRanges_and_profiles = LSY4377_12B_4_10bp_melting_profiles, GRanges = c(LSY4377_12B_1, LSY4377_12B_2, LSY4377_12B_4))

roi <- DSB_regions(DSBs = DSBs, region_width = 4000, up_rev_down_fw = TRUE)
hits <- findOverlaps(subject = LSY4377_12B_4_10bp_melting_profiles$GRanges, query = roi)
idx <- 46:56
tmp <- apply(X = (LSY4377_12B_4_10bp_melting_profiles$meltability_profiles[subjectHits(hits), idx] - melt_profile[idx])^2, MARGIN = 1, FUN = sum)
tmp <- max(tmp) - tmp

LSY4377_12B_1$melting <- tmp
LSY4377_12B_2$melting <- tmp
LSY4377_12B_4$melting <- tmp

LSY4377_12B_4_10bp_melting_profiles$GRanges


# calculate up/down ratios ================================================

calc_up_down_ratio <- function(x, mcol, f){ do.call(what = f, args = list(mcols(x[strand(x) == "-"])[names(mcols(x)) == mcol][, 1])) / do.call(what = f, args = list(mcols(x[strand(x) == "+"])[names(mcols(x)) == mcol][, 1])) }

LSY4377_12B_all_t <- c(LSY4377_12B_1, LSY4377_12B_2, LSY4377_12B_4)

up_down <- data.frame(DSB_id = unique(LSY4377_12B_all_t$DSB_id))
up_down <- cbind(up_down, S1_sum_ratio = sapply(X = split(LSY4377_12B_all_t, ~ DSB_id), FUN = calc_up_down_ratio, mcol = "score", f = sum))
up_down <- cbind(up_down, MNase_mean_ratio = sapply(X = split(LSY4377_12B_all_t, ~ DSB_id), FUN = calc_up_down_ratio, mcol = "MNase_seq", f = mean))
up_down <- cbind(up_down, MNase_var_ratio = sapply(X = split(LSY4377_12B_all_t, ~ DSB_id), FUN = calc_up_down_ratio, mcol = "MNase_seq", f = var))
up_down <- cbind(up_down, RNA_mean_ratio = sapply(X = split(LSY4377_12B_all_t, ~ DSB_id), FUN = calc_up_down_ratio, mcol = "RNA_seq", f = mean))
up_down <- cbind(up_down, RNA_var_ratio = sapply(X = split(LSY4377_12B_all_t, ~ DSB_id), FUN = calc_up_down_ratio, mcol = "RNA_seq", f = var))
up_down <- cbind(up_down, nt_mean_ratio = sapply(X = split(LSY4377_12B_all_t, ~ DSB_id), FUN = calc_up_down_ratio, mcol = "nt_seq_score", f = mean))
up_down <- cbind(up_down, nt_var_ratio = sapply(X = split(LSY4377_12B_all_t, ~ DSB_id), FUN = calc_up_down_ratio, mcol = "nt_seq_score", f = var))
up_down <- cbind(up_down, melt_mean_ratio = sapply(X = split(LSY4377_12B_all_t, ~ DSB_id), FUN = calc_up_down_ratio, mcol = "melting", f = mean))
up_down <- cbind(up_down, melt_var_ratio = sapply(X = split(LSY4377_12B_all_t, ~ DSB_id), FUN = calc_up_down_ratio, mcol = "melting", f = var))

# add prevailing directionality up/downstream of DSB
add_prevailing_directionality <- function(GRanges){
  direct <- rep(NA, length(GRanges))
  direct[GRanges$co_direction == TRUE] <- 1
  direct[GRanges$co_direction == FALSE] <- -1
  direct[is.na(direct)] <- 0
  up <- sum(direct[as.character(strand(GRanges)) == "-"] * GRanges$RNA_seq[as.character(strand(GRanges)) == "-"])
  down <- sum(direct[as.character(strand(GRanges)) == "+"] * GRanges$RNA_seq[as.character(strand(GRanges)) == "+"])
  return(data.frame(up = up, down = down))
}

tmp <- matrix(data = unlist(lapply(X = split(LSY4377_12B_all_t, ~ DSB_id), FUN = add_prevailing_directionality)), ncol = 2, byrow = TRUE)
up_down$up_direct <- tmp[, 1]
up_down$down_direct <- tmp[, 2]

# plotting function
plot_correlation <- function(x, y, ylab){
  plot(x = x, y = y, pch = 20, xlab = expression("Sum S1-seq Sum Ratio"), ylab = ylab)
  # abline(lm(y ~ x), col = "red")  # add linear regression line
  # add Pearson's correlation coefficient
  text(x = par("usr")[1] + 0.1 * (par("usr")[2] - par("usr")[1]), y = par("usr")[4] - 0.1 * (par("usr")[4]-par("usr")[3]), adj = c(0, 1),
       labels = paste0("r = ", round(cor(x = x, y = y, method = "spearman"), 2), " (p = ", round(cor.test(x = x, y = y, alternative = "two.sided", method = "spearman")$p.value, 2), ")"))
}


# plots
pdf(file = "tmp.pdf", width = 3, height = 3)
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(1.7, 0.5, 3.8, 1.8), tcl = -0.3, mgp = c(2.5, 0.6, 0), las = 1)
plot_correlation(x = up_down$S1_sum_ratio, y = up_down$MNase_mean_ratio, ylab = "MNase-seq Mean Ratio")
dev.off()
GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/MNase-seq_mean.pdf"))

pdf(file = "tmp.pdf", width = 3, height = 3)
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(1.7, 0.5, 3.8, 1.8), tcl = -0.3, mgp = c(2.5, 0.6, 0), las = 1)
plot_correlation(x = up_down$S1_sum_ratio, y = up_down$MNase_var_ratio, ylab = "MNase-seq Variance Ratio")
dev.off()
GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/MNase-seq_var.pdf"))

up_down[which.max(up_down$MNase_var_ratio), ]

pdf(file = "tmp.pdf", width = 3, height = 3)
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(1.7, 0.5, 3.8, 1.8), tcl = -0.3, mgp = c(2.5, 0.6, 0), las = 1)
plot_correlation(x = up_down$S1_sum_ratio, y = up_down$RNA_mean_ratio, ylab = "RNA-seq Mean Ratio")
dev.off()
GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/RNA-seq_mean.pdf"))

pdf(file = "tmp.pdf", width = 3, height = 3)
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(1.7, 0.5, 3.8, 1.8), tcl = -0.3, mgp = c(2.5, 0.6, 0), las = 1)
plot_correlation(x = up_down$S1_sum_ratio, y = up_down$RNA_var_ratio, ylab = "RNA-seq Variance Ratio")
dev.off()
GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/RNA-seq_var.pdf"))

pdf(file = "tmp.pdf", width = 3, height = 3)
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(1.7, 0.5, 3.8, 1.8), tcl = -0.3, mgp = c(2.5, 0.6, 0), las = 1)
plot_correlation(x = up_down$S1_sum_ratio, y = up_down$nt_mean_ratio, ylab = "Nt Similarity Mean Ratio")
dev.off()
GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/nt_mean.pdf"))

pdf(file = "tmp.pdf", width = 3, height = 3)
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(1.7, 0.5, 3.8, 1.8), tcl = -0.3, mgp = c(2.5, 0.6, 0), las = 1)
plot_correlation(x = up_down$S1_sum_ratio, y = up_down$nt_var_ratio, ylab = "Nt Similarity Variance Ratio")
dev.off()
GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/nt_var.pdf"))

pdf(file = "tmp.pdf", width = 3, height = 3)
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(1.7, 0.5, 3.8, 1.8), tcl = -0.3, mgp = c(2.5, 0.6, 0), las = 1)
plot_correlation(x = up_down$S1_sum_ratio, y = up_down$melt_mean_ratio, ylab = "Melting Similarity Mean Ratio")
dev.off()
GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/melt_mean.pdf"))

pdf(file = "tmp.pdf", width = 3, height = 3)
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(1.7, 0.5, 3.8, 1.8), tcl = -0.3, mgp = c(2.5, 0.6, 0), las = 1)
plot_correlation(x = up_down$S1_sum_ratio, y = up_down$melt_var_ratio, ylab = "Melting Similarity Variance Ratio      ")
dev.off()
GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/melt_var.pdf"))