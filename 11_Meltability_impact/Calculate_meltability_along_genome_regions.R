# info --------------------------------------------------------------------
# purpose: calculate DNA meltability along genome regions
# author: Robert Gnuegge (robert.gnuegge@gmail.com)
# date: 04/30/22
# version: 2.0


# preamble ----------------------------------------------------------------

# set wd to this file's location
wd_path <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(wd_path)

# libraries
library(GenomicRanges)
library(BSgenome)
library(rmelting)
library(parallel)
library(tictoc)


# read modified S. cerevisiae genome
genome <- import(con = "~/Research/Resources/S_cerevisiae_genome/Modified_SacCer3/S288C_R64-2-1_W303_SNPs_MATa_hocs2SrfIcs_hml_hmr.fasta")

# load helper files and functions
source(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Analysis_scripts/Misc/Genomic_helper_functions.R")
source(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Analysis_scripts/Misc/S_cerevisiae_SrfI_cut_sites.R")
source(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Analysis_scripts/Misc/Misc_helper_functions.R")

save_dir <- "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Data/Melting"


# function definitions ----------------------------------------------------

# calculate meltability for DNA sequence
# argument: character string
# result: data.frame with rows enthalpy.J, entropy.J, and temperature.C
calc_melting <- function(sequence, nucleic.acid.conc = 5.54e-10, Na.conc = 0.02, K.conc = 0.3, Mg.conc = 0.002, hybridisation.type = "dnadna"){
  tmp <- melting(sequence = sequence, nucleic.acid.conc = nucleic.acid.conc, 
                 Na.conc = Na.conc, K.conc = K.conc, Mg.conc = Mg.conc, hybridisation.type = hybridisation.type)
  out <- unlist(tmp$Results[5])
  return(out)
}


# extract meltability profiles
# argument: GRanges
# result: list with GRanges and matrix
extract_T_m_profiles <- function(GRanges, GRanges_w_T_m, width, DSBs){
  tiles <- resize(x = GRanges, width = width, fix = "center")
  hits <- findOverlaps(query = DSBs, subject = tiles)
  if(length(hits) > 0){
    tiles <- tiles[-subjectHits(hits)]  # remove tiles that overlap DSBs
    GRanges <- GRanges[-subjectHits(hits)]  # also remove corresponding GRanges
  }
  profiles <- matrix(data = 0, nrow = length(tiles), ncol = width)
  pb <- txtProgressBar(min = 0, max = length(tiles), initial = 0, style = 3, width = 76)  # show progress bar
  for(n in 1:length(tiles)){
    hits <- findOverlaps(query = tiles[n], subject = GRanges_w_T_m)
    profiles[n, ] <- GRanges_w_T_m[subjectHits(hits)]$T_m
    setTxtProgressBar(pb, n)
  }
  GRanges$T_m_profile <- profiles
  return(GRanges)
}


# derive sequences where meltability is to be calculated ------------------

# define genome roi
load(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Data/S1-seq/LSY4377-12B_4377-15A_merged/LSY4377-12B_LSY4377-15A_merged_resection_extend.RData")
DSBs <- SrfIcs[- c(9, 17)]  # exclude SrfIcs in duplicated region
max(c(start(DSBs) - start(c(resection_extend_LSY4377_12B_1_merged, resection_extend_LSY4377_12B_2_merged, resection_extend_LSY4377_12B_4_merged)),
      end(c(resection_extend_LSY4377_12B_1_merged, resection_extend_LSY4377_12B_2_merged, resection_extend_LSY4377_12B_4_merged)) - end(DSBs)))
# 1317
roi <- DSB_regions(DSBs = DSBs, region_width = 2720, up_rev_down_fw = TRUE)


# 10 bp ===================================================================
# get GRanges around all potential nick sites
melting_window <- 10  # define how long the molten DNA stretch is
pos <- as_nt_resolved_GRanges(roi)
tiles <- resize(x = pos, width = melting_window, fix = "center")

# get sequences
sequences <- as.character(getSeq(x = genome, names = tiles))


# calculate melting temperature for each sequence -------------------------
# let's use parallel computation for speed-up
# WARNING: depending on your hardware and the chosen melting_window, this will run for several hours!
# E.g. it ran ca. 1 h on my laptop with melting_window = 10.
tic(msg = "Melting calculation")
n_cores <- detectCores()
cl <- makeCluster(spec = n_cores - 2)  # start cluster
clusterEvalQ(cl = cl, expr = {library(rmelting)})  # load required libraries
rmelting_out <- parSapply(cl = cl, X = sequences, FUN = calc_melting)
stopCluster(cl)  # stop cluster
toc()  # Melting calculation: 3611.069 sec elapsed


# construct GRanges object with calculated melting parameters for each nt
melting_results <- pos
melting_results$T_m <- rmelting_out


# save results
melting_10bp <- melting_results
save(melting_10bp, file = paste0(save_dir, "/Melting_10bp_v2.RData"))
load(file = paste0(save_dir, "/Melting_10bp_v2.RData"))

# derive T_m profiles around all potential nick sites ---------------------
pos <- as_nt_resolved_GRanges(DSB_regions(DSBs = DSBs, region_width = 2640, up_rev_down_fw = TRUE))
Melting_profiles_10bp <- extract_T_m_profiles(GRanges = pos, GRanges_w_T_m = melting_10bp, width = 70, DSBs = DSBs)



# 5 bp ====================================================================
# get GRanges around all potential nick sites
melting_window <- 5  # define how long the molten DNA stretch is
pos <- as_nt_resolved_GRanges(roi)
tiles <- resize(x = pos, width = melting_window, fix = "center")

# get sequences
sequences <- as.character(getSeq(x = genome, names = tiles))


# calculate melting temperature for each sequence -------------------------
# let's use parallel computation for speed-up
tic(msg = "Melting calculation")
n_cores <- detectCores()
cl <- makeCluster(spec = n_cores - 2)  # start cluster
clusterEvalQ(cl = cl, expr = {library(rmelting)})  # load required libraries
rmelting_out <- parSapply(cl = cl, X = sequences, FUN = calc_melting)
stopCluster(cl)  # stop cluster
toc()  # Melting calculation: 2597.922 sec elapsed


# construct GRanges object with calculated melting parameters for each nt
melting_results <- pos
melting_results$T_m <- rmelting_out


# save results
melting_5bp <- melting_results
save(melting_5bp, file = paste0(save_dir, "/Melting_5bp_v2.RData"))
load(file = paste0(save_dir, "/Melting_5bp_v2.RData"))

# derive T_m profiles around all potential nick sites ---------------------
pos <- as_nt_resolved_GRanges(DSB_regions(DSBs = DSBs, region_width = 2640, up_rev_down_fw = TRUE))
Melting_profiles_5bp <- extract_T_m_profiles(GRanges = pos, GRanges_w_T_m = melting_5bp, width = 70, DSBs = DSBs)

# 15 bp ===================================================================
# get GRanges around all potential nick sites
melting_window <- 15  # define how long the molten DNA stretch is
pos <- as_nt_resolved_GRanges(roi)
tiles <- resize(x = pos, width = melting_window, fix = "center")

# get sequences
sequences <- as.character(getSeq(x = genome, names = tiles))


# calculate melting temperature for each sequence -------------------------
# let's use parallel computation for speed-up
tic(msg = "Melting calculation")
n_cores <- detectCores()
cl <- makeCluster(spec = n_cores - 2)  # start cluster
clusterEvalQ(cl = cl, expr = {library(rmelting)})  # load required libraries
rmelting_out <- parSapply(cl = cl, X = sequences, FUN = calc_melting)
stopCluster(cl)  # stop cluster
toc()  # Melting calculation: 4559.781 sec elapsed


# construct GRanges object with calculated melting parameters for each nt
melting_results <- pos
melting_results$T_m <- rmelting_out


# save results
melting_15bp <- melting_results
save(melting_15bp, file = paste0(save_dir, "/Melting_15bp_v2.RData"))
load(file = paste0(save_dir, "/Melting_15bp_v2.RData"))

# derive T_m profiles around all potential nick sites ---------------------
pos <- as_nt_resolved_GRanges(DSB_regions(DSBs = DSBs, region_width = 2640, up_rev_down_fw = TRUE))
Melting_profiles_15bp <- extract_T_m_profiles(GRanges = pos, GRanges_w_T_m = melting_15bp, width = 70, DSBs = DSBs)


# save results ============================================================
save(list = c("Melting_profiles_5bp", "Melting_profiles_10bp", "Melting_profiles_15bp"), file = paste0(save_dir, "/Melting_profiles_v2.RData"))
