# info --------------------------------------------------------------------
# purpose: derive S1-seq coverage for meiotic cells (from Mimitou et al., 2017)
# author: Robert Gnuegge (robert.gnuegge@gmail.com)
# date: 02/12/22
# version: 2.0


# preamble ----------------------------------------------------------------

# set wd to this file's location
wd_path <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(wd_path)

# libraries
library(GenomicRanges)

# load helper files and functions
source(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Analysis_scripts/Misc/Genomic_helper_functions.R")
source(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Analysis_scripts/Misc/Misc_helper_functions.R")
source(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Analysis_scripts/Misc/JFly_colors.R")

data_dir <- "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Resources/Mimitou2017/"
save_dir <- "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Data/Mimitou2017/Coverage"

# function definitions ----------------------------------------------------

# merge two replicates
# argument: GRanges objects
# result: GRanges object
merge_replicates <- function(rep1, rep2){
  merged <- as_nt_resolved_GRanges(GRanges = reduce(c(rep1, rep2)))
  
  score_1 <- rep(0, length(merged))
  hits <- findOverlaps(query = merged, subject = rep1)
  score_1[queryHits(hits)] <- rep1[subjectHits(hits)]$score
  
  score_2 <- rep(0, length(merged))
  hits <- findOverlaps(query = merged, subject = rep2)
  score_2[queryHits(hits)] <- rep2[subjectHits(hits)]$score
  merged$score <- apply(X = cbind(score_1, score_2), MARGIN = 1, FUN = mean)
  
  return(merged)
}

# round score valued n.5 half the time to n and half the time to n + 1 (not to change overall distribution)
# argument: GRanges
# result: GRanges
round_half_scores <- function(GRanges){
  idx <- which(GRanges$score %% 1 == 0.5)
  half_scores <- GRanges$score[idx]
  GRanges$score[idx] <- sapply(X = half_scores, FUN = function(i){ ifelse(test = sample(x = c(0,1), size = 1), yes = floor(i), no = ceiling(i)) })
  return(GRanges)
}


# process S1-seq data =====================================================

# wt ----------------------------------------------------------------------

# read S1-seq coverage from Mimitou et al., 2017 (PMID: 28059759)
data_1 <- read.table(file = paste0(data_dir, "GSE85253_wt1.txt"), header = TRUE)
data_2 <- read.table(file = paste0(data_dir, "GSE85253_wt2.txt"), header = TRUE)

colnames(data_1)
colnames(data_2)

# make GRanges
wt_4_1 <- GRanges(seqnames = rep(paste0("chr", as.roman(gsub(pattern = "chr", replacement = "", x = data_1$chr))), 2),
                    ranges = rep(IRanges(start = data_1$pos, width = 1), 2),
                    strand = c(rep("+", nrow(data_1)), rep("-", nrow(data_1))),
                    score = c(data_1$w4, data_1$c4))

wt_4_2 <- GRanges(seqnames = rep(paste0("chr", as.roman(gsub(pattern = "chr", replacement = "", x = data_2$chr))), 2),
                    ranges = rep(IRanges(start = data_2$pos, width = 1), 2),
                    strand = c(rep("+", nrow(data_2)), rep("-", nrow(data_2))),
                    score = c(data_2$w4, data_2$c4))


# reverse RPM scaling, such that sequence context can be recorded "score times"
wt_4_1$score <- round(wt_4_1$score / min(wt_4_1$score[wt_4_1$score > 0]))
wt_4_2$score <- round(wt_4_2$score / min(wt_4_2$score[wt_4_2$score > 0]))

# merge replicates
wt_4_merged <- merge_replicates(rep1 = wt_4_1, rep2 = wt_4_2)
wt_4_merged <- round_half_scores(wt_4_merged)


# exo1 --------------------------------------------------------------------

# read S1-seq coverage from Mimitou et al., 2017 (PMID: 28059759)
data_1 <- read.table(file = paste0(data_dir, "GSE85253_exo1.txt"), header = TRUE)
data_2 <- read.table(file = paste0(data_dir, "GSE85253_exo2.txt"), header = TRUE)

colnames(data_1)
colnames(data_2)

# make GRanges
exo1_4_1 <- GRanges(seqnames = rep(paste0("chr", as.roman(gsub(pattern = "chr", replacement = "", x = data_1$chr))), 2),
                    ranges = rep(IRanges(start = data_1$pos, width = 1), 2),
                    strand = c(rep("+", nrow(data_1)), rep("-", nrow(data_1))),
                    score = c(data_1$w4, data_1$c4))

exo1_4_2 <- GRanges(seqnames = rep(paste0("chr", as.roman(gsub(pattern = "chr", replacement = "", x = data_2$chr))), 2),
                    ranges = rep(IRanges(start = data_2$pos, width = 1), 2),
                    strand = c(rep("+", nrow(data_2)), rep("-", nrow(data_2))),
                    score = c(data_2$w4, data_2$c4))

# reverse RPM scaling, such that sequence context can be recorded "score times"
exo1_4_1$score <- round(exo1_4_1$score / min(exo1_4_1$score[exo1_4_1$score > 0]))
exo1_4_2$score <- round(exo1_4_2$score / min(exo1_4_2$score[exo1_4_2$score > 0]))

# merge replicates
exo1_4_merged <- merge_replicates(rep1 = exo1_4_1, rep2 = exo1_4_2)
exo1_4_merged <- round_half_scores(exo1_4_merged)


# mask high-background regions --------------------------------------------
mask_top <- read.csv(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Resources/Mimitou2017/mask_coordinates_top.csv", header = TRUE)
mask_bottom  <- read.csv(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Resources/Mimitou2017/mask_coordinates_bottom.csv", header = TRUE)

# make GRanges
masked <- GRanges(seqnames = paste0("chr", as.roman(c(mask_top$CHROM, mask_bottom$CHROM))), 
                  ranges = IRanges(start = c(mask_top$POS, mask_bottom$POS), width = 1), 
                  strand = c(rep("+", nrow(mask_top)), rep("-", nrow(mask_bottom))))

subsetByOverlaps(wt_4_merged, x = masked)  # masking was already applied in downloaded S1-seq data
subsetByOverlaps(exo1_4_merged, x = masked)


# save data ---------------------------------------------------------------
save(list = c(paste0("exo1_4_", c("1", "2", "merged")), paste0("wt_4_", c("1", "2", "merged"))), file = paste0(save_dir, "/S1-seq_coverage.RData"))