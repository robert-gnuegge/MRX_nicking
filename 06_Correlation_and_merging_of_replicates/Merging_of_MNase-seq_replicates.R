# info --------------------------------------------------------------------
# purpose: merge (average) data of biological replicates
# author: Robert Gnuegge (robert.gnuegge@gmail.com)
# date: 01/26/21
# version: 2.0

# preamble ----------------------------------------------------------------

# set wd to this file's location
wd_path <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(wd_path)

# libraries
library(GenomicAlignments)

# load helper files and functions
source(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Analysis_scripts/Misc/Genomic_helper_functions.R")
source(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Analysis_scripts/Misc/Misc_helper_functions.R")

rep1_dir <- "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Data/MNase-seq/LSY4377-12B_4377-15A"
rep2_dir <- "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Data/MNase-seq/LSY4377-12B_4377-15A_rep2"
save_dir <- "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Data/MNase-seq/LSY4377-12B_4377-15A_merged"

# function definitions ----------------------------------------------------

# remove chrM and micron GRanges
# argument: GRanges
# result: GRanges
remove_chrM_and_2micron <- function(GRanges){
  GRanges[!seqnames(GRanges) %in% c("chrM", "micron")]
}

# merge two replicates
# argument: GRanges
# result: GRanges
# note: only GRanges with numeric mcols are supported
merge_replicates <- function(rep1, rep2, verbose = FALSE){
  stopifnot(names(mcols(rep1)) %in% names(mcols(rep2)))
  if(verbose){
    cat("\nMaking sorted nt-resolved GRanges...")
  }
  rep1 <- sort(as_nt_resolved_GRanges(GRanges = rep1))
  rep2 <- sort(as_nt_resolved_GRanges(GRanges = rep2))
  if(verbose){
    cat("\nChecking identity of all granges of replicates...")
  }
  stopifnot(all(granges(rep1) == granges(rep2)))
  if(verbose){
    cat("\nSetting up merged GRanges object...")
  }
  out <- granges(rep1)
  mcol_names <- names(mcols(rep1))
  mcols(out) <- matrix(data = rep(NA, length(out) * length(mcol_names)), ncol = length(mcol_names))
  names(mcols(out)) <- mcol_names
  for(n in 1:length(mcol_names)){
    if(verbose){
      cat("\nMerging mcol '", mcol_names[n], "'...", sep = "")
    }
    mcol1 <- mcols(rep1)[names(mcols(rep1)) == mcol_names[n]]
    mcol2 <- mcols(rep2)[names(mcols(rep2)) == mcol_names[n]]
    mcols(out)[names(mcols(out)) == mcol_names[n]] <- apply(X = cbind(mcol1$score, mcol2$score), MARGIN = 1, FUN = mean)
  }
  return(out)
}


# merge MNase-seq files ---------------------------------------------------

load(file = paste0(rep1_dir, "/LSY4377-12B_LSY4377-15A_trimmed.RData"))
load(file = paste0(rep2_dir, "/LSY4377-12B_LSY4377-15A_rep2_trimmed.RData"))

LSY4377_12B_0_merged_MNase_seq <- merge_replicates(rep1 = remove_chrM_and_2micron(GRanges = LSY4377_12B_0_MNase_seq_trimmed),
                                                   rep2 = remove_chrM_and_2micron(GRanges = LSY4377_12B_0_rep2_MNase_seq_trimmed),
                                                   verbose = TRUE)

LSY4377_12B_1_merged_MNase_seq <- merge_replicates(rep1 = remove_chrM_and_2micron(GRanges = LSY4377_12B_1_MNase_seq_trimmed),
                                                   rep2 = remove_chrM_and_2micron(GRanges = LSY4377_12B_1_rep2_MNase_seq_trimmed),
                                                   verbose = TRUE)

LSY4377_12B_2_merged_MNase_seq <- merge_replicates(rep1 = remove_chrM_and_2micron(GRanges = LSY4377_12B_2_MNase_seq_trimmed),
                                                   rep2 = remove_chrM_and_2micron(GRanges = LSY4377_12B_2_rep2_MNase_seq_trimmed),
                                                   verbose = TRUE)

LSY4377_12B_4_merged_MNase_seq <- merge_replicates(rep1 = remove_chrM_and_2micron(GRanges = LSY4377_12B_4_MNase_seq_trimmed),
                                                   rep2 = remove_chrM_and_2micron(GRanges = LSY4377_12B_4_rep2_MNase_seq_trimmed),
                                                   verbose = TRUE)

save(list = paste0("LSY4377_12B_", c(0, 1, 2, 4), "_merged_MNase_seq"), file = paste0(save_dir, "/LSY4377-12B_LSY4377-15A_merged_trimmed.RData"))