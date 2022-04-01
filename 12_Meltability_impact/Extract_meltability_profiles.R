# info --------------------------------------------------------------------
# purpose: Extract meltability profiles around DSBs
# author: Robert Gnuegge (robert.gnuegge@gmail.com)
# date: 01/12/22
# version: 2.0


# preamble ----------------------------------------------------------------

# set wd to this file's location
wd_path <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(wd_path)

# libraries
library(GenomicRanges)


# load helper files and functions
source(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Analysis_scripts/Misc/Genomic_helper_functions.R")
source(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Analysis_scripts/Misc/S_cerevisiae_SrfI_cut_sites.R")
source(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Analysis_scripts/Misc/Misc_helper_functions.R")
source(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Analysis_scripts/Misc/JFly_colors.R")

save_dir <- "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Data/Melting"


# function definitions ----------------------------------------------------

# collect meltability profiles
# argument: GRanges
# result: list with GRanges and matrix
get_meltability_profiles <- function(GRanges, meltability, width){
  tiles <- resize(x = GRanges, width = width, fix = "center")
  hits <- findOverlaps(query = SrfIcs[-c(9, 17)], subject = tiles)
  tiles <- tiles[-subjectHits(hits)]  # remove tiles that overlap DSBs
  GRanges <- GRanges[-subjectHits(hits)]  # also remove corresponding GRanges
  profiles <- matrix(data = 0, nrow = length(tiles), ncol = unique(width(tiles)))
  pb <- txtProgressBar(min = 0, max = length(tiles), initial = 0, style = 3, width = 76)  # show progress bar
  for(n in 1:length(tiles)){
    hits <- findOverlaps(query = tiles[n], subject = meltability)
    profiles[n, ] <- mcols(meltability[subjectHits(hits)])$temperature.C
    setTxtProgressBar(pb, n)
  }
  profiles[which(strand(tiles) == "-"), ] <- t(apply(X = profiles[which(strand(tiles) == "-"), ], MARGIN = 1, FUN = rev))  # for profiles from opposite strand (rev direction relative to nick)
  out <- list(GRanges = GRanges, meltability_profiles = profiles)
  return(out)
}



# extract meltability profiles ============================================

load(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Data/S1-seq/LSY4377-12B_4377-15A/LSY4377-12B_LSY4377-15A.RData")
load(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Data/S1-seq/LSY4377-12B_4377-15A/LSY4377-12B_LSY4377-15A_resection_extend.RData")

load(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Data/Melting/Melting_temperatures/Melting_5bp.RData")
load(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Data/Melting/Melting_temperatures/Melting_10bp.RData")
load(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Data/Melting/Melting_temperatures/Melting_15bp.RData")

# check that resection extend at t = 4 contains earlier t
all(start(resection_extend_LSY4377_12B_2) < start(resection_extend_LSY4377_12B_1))
all(start(resection_extend_LSY4377_12B_4) < start(resection_extend_LSY4377_12B_2))
all(end(resection_extend_LSY4377_12B_1) < end(resection_extend_LSY4377_12B_2))
all(end(resection_extend_LSY4377_12B_2) < end(resection_extend_LSY4377_12B_4))
# TRUE
# this allow to get melting profiles for t = 4 and use them for calculations with all t

# keep only regions with correct orientation w.r.t DSBs
LSY4377_12B_4 <- subsetByIntersect(subject = LSY4377_12B_4_S1_seq, query = DSB_regions(DSBs = SrfIcs[-c(9, 17)], region_width = 4000, up_rev_down_fw = TRUE))

# retain GRanges in regions around DSBs where resection was detected
# other regions are excluded to prevent pattern "dilution" by background noise
LSY4377_12B_4 <- subsetByIntersect(subject = LSY4377_12B_4, query = resection_extend_LSY4377_12B_4)

# make nt-resolved GRanges
LSY4377_12B_4 <- as_nt_resolved_GRanges(GRanges = LSY4377_12B_4)


# calculate average melting profiles (scaled by S1-seq and bg) ------------

# get melting profiles for melting windows
LSY4377_12B_4_5bp_melting_profiles <- get_meltability_profiles(GRanges = LSY4377_12B_4, meltability = melting_5bp, width = 100)
LSY4377_12B_4_10bp_melting_profiles <- get_meltability_profiles(GRanges = LSY4377_12B_4, meltability = melting_10bp, width = 100)
LSY4377_12B_4_15bp_melting_profiles <- get_meltability_profiles(GRanges = LSY4377_12B_4, meltability = melting_15bp, width = 100)

save(list = paste0("LSY4377_12B_4_", c(5, 10, 15), "bp_melting_profiles"), file = paste0(save_dir, "/Extracted_melting_profiles.RData"))