# info --------------------------------------------------------------------
# purpose: check if there is indication for SrfI off-target cutting
# author: Robert Gnuegge (robert.gnuegge@gmail.com)
# date: 01/14/22
# version: 2.0


# approach ----------------------------------------------------------------
# preferred off-targets (if there are any) would probably be 1-mismatch-SrfIcs
# (compare star acitivity of restriction enzymes and NdeI control in HydEN-seq)
# let's compare S1-seq levels at SrfIcs, SrfIcs* (containing 1 mismatch), and background
# let's use the mre11-nd data set, as S1-seq levels are most focussed here


# preamble ----------------------------------------------------------------

# set wd to this file's location
wd_path <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(wd_path)

# libraries
library(GenomicRanges)
library(BSgenome)

# read modified S. cerevisiae genome
genome <- import(con = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Data/Reference_genomes//home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Data/S_cerevisiae_reference_genomes/S288C_R64-2-1_W303_SNPs_MATa_hocs2SrfIcs_hml_hmr.fasta")

# load helper files and functions
source(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Analysis_scripts/Misc/Genomic_helper_functions.R")
source(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Analysis_scripts/Misc/S_cerevisiae_SrfI_cut_sites.R")
source(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Analysis_scripts/Misc/Misc_helper_functions.R")
source(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Analysis_scripts/Misc/JFly_colors.R")

plot_dir <- "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Figures/05_DSB_formation/LSY4377-12B_4377-15A"


# define SrfIcs -----------------------------------------------------------

SrfIcs <- GRanges()
for(chr in seqlevels(genome)[1:16]){  # iterate through genomic sequences, exclude chrM and 2micron
  tmp <- matchPattern(pattern = "GCCCGGGC", subject = genome[[chr]])  # find all patterns
  if(length(tmp) > 0){
    SrfIcs <- append(SrfIcs, GRanges(seqnames = chr, ranges = ranges(tmp), strand = "+", nt.seq = as.character(tmp)))  # save results
  }
}

SrfIcs <- SrfIcs[-c(9, 17)]  # remove SrfIcs in duplicated region


# define SrfIcs* ----------------------------------------------------------

mm_SrfIcs <- GRanges()
for(chr in seqlevels(genome)[1:16]){  # iterate through genomic sequences, exclude chrM and 2micron
  tmp <- matchPattern(pattern = "GCCCGGGC", subject = genome[[chr]], max.mismatch = 1, min.mismatch = 1)  # find all patterns
  if(length(tmp) > 0){
    mm_SrfIcs <- append(mm_SrfIcs, GRanges(seqnames = chr, ranges = ranges(tmp), strand = "+", nt.seq = as.character(tmp)))  # save results
  }
}

mm_SrfIcs
# found 600 sites

# Are there SrfIcs (with 1 mismatch) in 2micron?
matchPattern(pattern = "GCCCGGGC", subject = genome[["micron"]])
# none
matchPattern(pattern = "GCCCGGGC", subject = genome[["micron"]], max.mismatch = 1, min.mismatch = 1)
# none

# Are there SrfIcs (with 1 mismatch) in mtDNA?
matchPattern(pattern = "GCCCGGGC", subject = genome[["chrM"]])
# none
matchPattern(pattern = "GCCCGGGC", subject = genome[["chrM"]], max.mismatch = 1, min.mismatch = 1)
# 46
# however, SrfIcs should not localize to mitochondria



# define bg ---------------------------------------------------------------

# construct masking GRanges object
load(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Data/S_cerevisiae_chromosomal_features/S_cerevisiae_genome_features.RData")

masking <- GRanges()

# add TELs (regions up/downstream of first/last gene on each chr)
all_features_by_chr <- split(all_features[!is.na(all_features$qualifier) & all_features$qualifier == "Verified"], ~seqnames)
first_genes <- unlist(endoapply(X = all_features_by_chr, FUN = function(x) {x[1]}))
TELs_left <- GRanges(seqnames = seqnames(first_genes), ranges = IRanges(start = 1, end = (start(first_genes) - 1)))
last_genes <- unlist(endoapply(X = all_features_by_chr, FUN = function(x) {x[length(x)]}))
TELs_right <- GRanges(seqnames = seqnames(last_genes), ranges = IRanges(start = (end(last_genes) + 1), end = seqlengths(genome)[1:17]))
masking <- c(masking, TELs_left, TELs_right)

# add rDNA locus (region between flanking all_features)
rDNA <- GRanges(seqnames = "chrXII", ranges = IRanges(start = end(all_features[all_features$name == "RNH203"]) + 1, end = start(all_features[all_features$name == "MAS1"]) - 1))
masking <- c(masking, rDNA)

# add SrfIcs and mm_SrfIcs
masking <- c(masking, resize(x = SrfIcs, width = 100, fix = "center"), resize(x = mm_SrfIcs, width = 100, fix = "center"))
strand(masking) <- "*"
masking <- reduce(masking)

genome_for_sampling <- GRanges(seqnames = seqlevels(genome), ranges = IRanges(start = 1, end = seqlengths(genome)))[1:16]  # exclude chrM and 2micron
genome_for_sampling <- as_nt_resolved_GRanges(genome_for_sampling)
hits <- findOverlaps(query = masking, subject = genome_for_sampling)
genome_for_sampling <- genome_for_sampling[-subjectHits(hits)]

bg <- genome_for_sampling[sample(x = 1:length(genome_for_sampling), size = 1e4, replace = TRUE)]
strand(bg) <- "+"

# resize all GRanges to 2 nt on fw and rev strand
fw_rev_2nt_ranges <- function(GRanges){
  fw <- resize(x = shift(x = resize(x = GRanges, width = 1, fix = "center"), shift = 1), width = 1, fix = "start")
  rev <- fw
  strand(rev) <- "-"
  rev <- shift(x = rev, shift = -1)
  out <- sort(c(fw, rev), ignore.strand = TRUE)
  return(out)
}

SrfIcs <- fw_rev_2nt_ranges(GRanges = SrfIcs)
mm_SrfIcs <- fw_rev_2nt_ranges(GRanges = mm_SrfIcs)
bg <- fw_rev_2nt_ranges(bg)


# get S1-seq scores -------------------------------------------------------

load(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Data/S1-seq/LSY4377-12B_4377-15A/LSY4377-12B_LSY4377-15A.RData")

SrfIcs$score <- subsetByIntersect(query = SrfIcs, subject = LSY4377_15A_4_S1_seq)$score
mm_SrfIcs$score <- subsetByIntersect(query = mm_SrfIcs, subject = LSY4377_15A_4_S1_seq)$score
bg$score <- subsetByIntersect(query = bg, subject = LSY4377_15A_4_S1_seq)$score


# plot --------------------------------------------------------------------

pairwise.wilcox.test(x = c(bg$score, mm_SrfIcs$score, SrfIcs$score), 
                     g = c(rep("bg", length(bg)), rep("SrfIcs*", length(mm_SrfIcs)), rep("SrfIcs", length(SrfIcs))),
                     p.adjust.method = "bonferroni")

add_significance <- function(x0, x1, y, text, text_offset = -2.5, x_shrink = 0.05){
  segments(x0 = (x0 + x_shrink), y0 = y, x1 = (x1 - x_shrink), y1 = y)
  text(x = x0 + diff(x = c(x0, x1))/2, y = y + text_offset, labels = text, pos = 3, xpd = TRUE)
}

pdf(file = "tmp.pdf", width=3, height=3.5)
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(3.6, 0.6, 0.4, 2), xpd = TRUE, las=1, tcl = -0.25, mgp = c(2.25, 0.5, 0))
  boxplot(list(bg$score / 1000, mm_SrfIcs$score / 1000, SrfIcs$score / 1000),
          names = c(expression("No"~italic("SrfIcs")), expression(italic("SrfIcs*")), expression(italic("SrfIcs"))),
          ylab = expression("S1-seq [RPM] (x" * 10 ^ 3 *")"))
  
  y0 = 19
  y_step = 2.75
  add_significance(x0 = 1, x1 = 2, y = y0, text_offset = -0.5, x_shrink = 0.1, text = expression("p = 0.2"))
  add_significance(x0 = 2, x1 = 3, y = y0, text_offset = -0.5, x_shrink = 0.1, text = expression("p < 10"^-15))
  add_significance(x0 = 1, x1 = 3, y = y0 + y_step, text_offset = -0.5, x_shrink = 0.1, text = expression("p < 10"^-15))
dev.off()

GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/SrfI_off_target_cutting.pdf"))
