# info --------------------------------------------------------------------
# purpose: get the cyclizability values around the SrfIcs on chrV
# author: Robert Gnuegge (robert.gnuegge@gmail.com)
# date: 01/18/21
# version: 2.0


# preamble ----------------------------------------------------------------

# set wd to this file's location
wd_path <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(wd_path)

# libraries
library(GenomicRanges)
library(BSgenome.Scerevisiae.UCSC.sacCer3)
genome <- BSgenome.Scerevisiae.UCSC.sacCer3

# load helper files and functions
source(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Analysis_scripts/Misc/Genomic_helper_functions.R")
source(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Analysis_scripts/Misc/S_cerevisiae_SrfI_cut_sites.R")
source(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Analysis_scripts/Misc/Misc_helper_functions.R")

save_dir <- "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Data/Cyclizability"


# read cyclizability data -------------------------------------------------
# from Basu et al., 2021 (pmid: 33328628)
cyclizability <- read.table(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Resources/Basu2021/41586_2020_3052_MOESM9_ESM.txt", header = TRUE, sep = "\t", comment.char="", stringsAsFactors = FALSE)

# remove flanking adapter sequences
cyclizability$Sequence <- substr(x = cyclizability$Sequence, start = 26, stop = 75)
# remove duplicates
cyclizability <- cyclizability[!duplicated(cyclizability$Sequence), ]


# get cyclizability values around SrfIcs at chrV:123366 -------------------

# define genome roi
roi <- DSB_regions(DSBs = SrfIcs[4], region_width = 10000)
pos_chrV_123366 <- as_nt_resolved_GRanges(roi)
tiles <- resize(x = pos_chrV_123366, width = nchar(cyclizability$Sequence[1]), fix = "center")
sequences <- as.character(getSeq(x = genome, names = tiles))

# cyclizability data set contains values for every 7th nt
# let's find first and last one matching roi GRanges

hit <- FALSE
first_idx_pos <- 0
while(!hit){
  first_idx_pos <- first_idx_pos + 1
  first_idx_cycl <- grep(pattern = sequences[first_idx_pos], x = cyclizability$Sequence, fixed = TRUE)
  if(length(first_idx_cycl) > 0){
    hit <- TRUE
  }
}
# sanity checks
sequences[first_idx_pos]
cyclizability$Sequence[first_idx_cycl]
sequences[first_idx_pos + 7]
cyclizability$Sequence[first_idx_cycl + 1]

hit <- FALSE
last_idx_pos <- length(sequences) + 1
while(!hit){
  last_idx_pos <- last_idx_pos - 1
  last_idx_cycl <- grep(pattern = sequences[last_idx_pos], x = cyclizability$Sequence, fixed = TRUE)
  if(length(last_idx_cycl) > 0){
    hit <- TRUE
  }
}
# sanity checks
sequences[last_idx_pos]
cyclizability$Sequence[last_idx_cycl]
sequences[last_idx_pos - 7]
cyclizability$Sequence[last_idx_cycl - 1]

pos_chrV_123366$value <- NA
pos_chrV_123366$value[seq(from = first_idx_pos, to = last_idx_pos, by = 7)] <- cyclizability$C0[first_idx_cycl:last_idx_cycl]

# fill in missing values by linear interpolation
pos_chrV_123366$interpolated <- pos_chrV_123366$value
idx <- which(!(is.na(pos_chrV_123366$value)))
n <- (idx[1]):(idx[length(idx)])
pos_chrV_123366$interpolated[n] <- approx(x = pos_chrV_123366$interpolated[n], n = length(n))$y  # interpolates between non-NA values

# get cyclizability values around SrfIcs at chrV:399646 -------------------

# define genome roi
roi <- DSB_regions(DSBs = SrfIcs[5], region_width = 10000)
pos_chrV_399646 <- as_nt_resolved_GRanges(roi)
tiles <- resize(x = pos_chrV_399646, width = nchar(cyclizability$Sequence[1]), fix = "center")
sequences <- as.character(getSeq(x = genome, names = tiles))

# cyclizability data set contains values for every 7th nt
# let's find first and last one matching roi GRanges

hit <- FALSE
first_idx_pos <- 0
while(!hit){
  first_idx_pos <- first_idx_pos + 1
  first_idx_cycl <- grep(pattern = sequences[first_idx_pos], x = cyclizability$Sequence, fixed = TRUE)
  if(length(first_idx_cycl) > 0){
    hit <- TRUE
  }
}
# sanity checks
sequences[first_idx_pos]
cyclizability$Sequence[first_idx_cycl]
sequences[first_idx_pos + 7]
cyclizability$Sequence[first_idx_cycl + 1]

hit <- FALSE
last_idx_pos <- length(sequences) + 1
while(!hit){
  last_idx_pos <- last_idx_pos - 1
  last_idx_cycl <- grep(pattern = sequences[last_idx_pos], x = cyclizability$Sequence, fixed = TRUE)
  if(length(last_idx_cycl) > 0){
    hit <- TRUE
  }
}
# sanity checks
sequences[last_idx_pos]
cyclizability$Sequence[last_idx_cycl]
sequences[last_idx_pos - 7]
cyclizability$Sequence[last_idx_cycl - 1]

pos_chrV_399646$value <- NA
pos_chrV_399646$value[seq(from = first_idx_pos, to = last_idx_pos, by = 7)] <- cyclizability$C0[first_idx_cycl:last_idx_cycl]

# fill in missing values by linear interpolation
pos_chrV_399646$interpolated <- pos_chrV_399646$value
idx <- which(!(is.na(pos_chrV_399646$value)))
n <- (idx[1]):(idx[length(idx)])
pos_chrV_399646$interpolated[n] <- approx(x = pos_chrV_399646$interpolated[n], n = length(n))$y  # interpolates between non-NA values

# save
Cyclizability <- c(pos_chrV_123366, pos_chrV_399646)
save(list = c("Cyclizability"), file = paste0(save_dir, "/Cyclizability.RData"))


# extract cyclizability profiles ==========================================

# collect cyclizability profiles function
# argument: GRanges
# result: list with GRanges and matrix
get_cyclizability_profiles <- function(GRanges, cyclizability, width, interpolated = FALSE){
  tiles <- resize(x = GRanges, width = width, fix = "center")
  hits <- findOverlaps(query = SrfIcs[4:5], subject = tiles)
  tiles <- tiles[-subjectHits(hits)]  # remove tiles that overlap DSBs
  GRanges <- GRanges[-subjectHits(hits)]  # also remove corresponding GRanges
  profiles <- matrix(data = NA, nrow = length(tiles), ncol = width)
  pb <- txtProgressBar(min = 0, max = length(tiles), initial = 0, style = 3, width = 76)  # show progress bar
  for(n in 1:length(tiles)){
    hits <- findOverlaps(query = tiles[n], subject = cyclizability)
    if(interpolated){
      profiles[n, ] <- mcols(cyclizability[subjectHits(hits)])$interpolated
    }else{
      profiles[n, ] <- mcols(cyclizability[subjectHits(hits)])$value
    }
    setTxtProgressBar(pb, n)
  }
  profiles[which(strand(tiles) == "-"), ] <- t(apply(X = profiles[which(strand(tiles) == "-"), ], MARGIN = 1, FUN = rev))  # for profiles from opposite strand (rev direction relative to nick)
  out <- list(GRanges = GRanges, cyclizability_profiles = profiles)
  return(out)
}


# load and process data ---------------------------------------------------

load(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Data/S1-seq/LSY4377-12B_4377-15A/LSY4377-12B_LSY4377-15A_unnormalized.RData")
# use unnormalized coverage data (for easy scaling with S1-seq score)
load(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Data/S1-seq/LSY4377-12B_4377-15A/LSY4377-12B_LSY4377-15A_resection_extend.RData")

# keep only regions with correct orientation w.r.t DSBs
LSY4377_12B_4 <- subsetByIntersect(subject = LSY4377_12B_4_S1_seq_unnormalized, query = DSB_regions(DSBs = SrfIcs[4:5], region_width = 4000, up_rev_down_fw = TRUE))

# retain GRanges in regions around DSBs where resection was detected
# other regions are excluded to prevent pattern "dilution" by background noise
LSY4377_12B_4 <- subsetByIntersect(subject = LSY4377_12B_4, query = resection_extend_LSY4377_12B_4)

# check that resection extend at t = 4 contains earlier t
all(start(resection_extend_LSY4377_12B_2) < start(resection_extend_LSY4377_12B_1))
all(start(resection_extend_LSY4377_12B_4) < start(resection_extend_LSY4377_12B_2))
all(end(resection_extend_LSY4377_12B_1) < end(resection_extend_LSY4377_12B_2))
all(end(resection_extend_LSY4377_12B_2) < end(resection_extend_LSY4377_12B_4))
# TRUE
# this allow to get cyclizability profiles for t = 4 and use them for calculations with all t

# make nt-resolved GRanges
LSY4377_12B_4 <- as_nt_resolved_GRanges(GRanges = LSY4377_12B_4)

# calculate average cyclizability profiles (scaled by S1-seq and bg) ------------

# get cyclizability profiles
LSY4377_12B_4_cyclizability_profiles_200bp <- get_cyclizability_profiles(GRanges = LSY4377_12B_4, cyclizability = Cyclizability, width = 200)
LSY4377_12B_4_cyclizability_profiles_500bp <- get_cyclizability_profiles(GRanges = LSY4377_12B_4, cyclizability = Cyclizability, width = 500)
LSY4377_12B_4_cyclizability_profiles_200bp_interpolated <- get_cyclizability_profiles(GRanges = LSY4377_12B_4, cyclizability = Cyclizability, width = 200, interpolated = TRUE)
LSY4377_12B_4_cyclizability_profiles_500bp_interpolated <- get_cyclizability_profiles(GRanges = LSY4377_12B_4, cyclizability = Cyclizability, width = 500, interpolated = TRUE)

save(list = paste0("LSY4377_12B_4_cyclizability_profiles_", c("200bp", "500bp", "200bp_interpolated", "500bp_interpolated")), file = paste0(save_dir, "/Cyclizability_chrV/Extracted_cyclizability_profiles.RData"))
