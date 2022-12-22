# info --------------------------------------------------------------------
# purpose: calculate DNA shape features along genome regions
# author: Robert Gnuegge (robert.gnuegge@gmail.com)
# date: 04/25/22
# version: 1.0


# preamble ----------------------------------------------------------------

# set wd to this file's location
wd_path <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(wd_path)

# libraries
library(GenomicRanges)
library(BSgenome)
library(DNAshapeR)
library(parallel)
library(tictoc)


# read modified S. cerevisiae genome
genome <- import(con = "~/Research/Resources/S_cerevisiae_genome/Modified_SacCer3/S288C_R64-2-1_W303_SNPs_MATa_hocs2SrfIcs_hml_hmr.fasta")

# load helper files and functions
source(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Analysis_scripts/Misc/Genomic_helper_functions.R")
source(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Analysis_scripts/Misc/S_cerevisiae_SrfI_cut_sites.R")
source(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Analysis_scripts/Misc/Misc_helper_functions.R")

save_dir <- "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Data/Shape"

# derive sequences where shape is to be calculated ------------------------

# define genome roi
load(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Data/S1-seq/LSY4377-12B_4377-15A_merged/LSY4377-12B_LSY4377-15A_merged_resection_extend.RData")
DSBs <- SrfIcs[- c(9, 17)]  # exclude SrfIcs in duplicated region
max(c(start(DSBs) - start(c(resection_extend_LSY4377_12B_1_merged, resection_extend_LSY4377_12B_2_merged, resection_extend_LSY4377_12B_4_merged)),
      end(c(resection_extend_LSY4377_12B_1_merged, resection_extend_LSY4377_12B_2_merged, resection_extend_LSY4377_12B_4_merged)) - end(DSBs)))
# 1317
roi <- DSB_regions(DSBs = DSBs, region_width = 2700, up_rev_down_fw = TRUE)

# get GRanges around all potential nick sites
shape_width <- 70
pos <- as_nt_resolved_GRanges(roi)
tiles <- resize(x = pos, width = shape_width, fix = "center")
# remove tiles that overlap with DSBs 
hits <- findOverlaps(query = DSBs, subject = tiles)
tiles <- tiles[-subjectHits(hits)]
pos <- pos[-subjectHits(hits)]

# get sequences
sequences <- as.character(getSeq(x = genome, names = tiles))

# save as fasta files (required for using getShape command) 
names(sequences) <- paste0(seqnames(pos), ":", start(pos))  # add names (will be fasta entry IDs)
sequences <- DNAStringSet(sequences)
file_path <- paste0(save_dir, "/Nt_sequences_70nt_windows_around_SrfIcs.fasta")
writeXStringSet(x = sequences, filepath = file_path)


# calculate DNA shape features for each sequence -------------------------
res <- getShape(filename = file_path)
lapply(X = res, FUN = dim)

# save results
DNA_shape <- pos
DNA_shape$MGW <- res$MGW
DNA_shape$Roll <- cbind(res$Roll, NA)  # add NA column for equal ncol
DNA_shape$HelT <- cbind(res$HelT, NA)  # add NA column for equal ncol
DNA_shape$ProT <- res$ProT
DNA_shape$EP <- res$EP

save(DNA_shape, file = paste0(save_dir, "/DNA_shape_around_SrfIcs.RData"))
