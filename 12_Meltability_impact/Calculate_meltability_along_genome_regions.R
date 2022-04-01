# info --------------------------------------------------------------------
# purpose: calculate DNA melting temperatures along genome
# author: Robert Gnuegge (robert.gnuegge@gmail.com)
# date: 01/13/21
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
  out <- unlist(tmp$Results[3:5])
  return(out)
}



# derive sequences where melting is to be calculated ----------------------

# define genome roi
load(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Data/S1-seq/LSY4377-12B_LSY4377-15A_resection_extend.RData")
DSBs <- SrfIcs[- c(9, 17)]  # exclude SrfIcs in duplicated region
max(c(start(DSBs) - start(c(resection_extend_LSY4377_12B_1, resection_extend_LSY4377_12B_2, resection_extend_LSY4377_12B_4)),
      end(c(resection_extend_LSY4377_12B_1, resection_extend_LSY4377_12B_2, resection_extend_LSY4377_12B_4)) - end(DSBs)))
# 1260
roi <- DSB_regions(DSBs = DSBs, region_width = 2650)

# get sequences
melting_window <- 10  # define how long the molten DNA stretch is
pos <- as_nt_resolved_GRanges(roi)
tiles <- resize(x = pos, width = melting_window, fix = "center")
sequences <- as.character(getSeq(x = genome, names = tiles))



# calculate meltablitity for each sequence --------------------------------

# let's use parallel computation for speed-up
# WARNING: depending on your hardware and the chosen melting_window, this will run for several hours!
# E.g. it ran ca. 4 h on my laptop with melting_window = 10.
n_cores <- detectCores()
cl <- makeCluster(spec = n_cores - 2)  # start cluster
clusterEvalQ(cl = cl, expr = {library(rmelting)})  # load required libraries
rmelting_out <- parSapply(cl = cl, X = sequences, FUN = calc_melting)
stopCluster(cl)  # stop cluster


# construct GRanges object with calculated melting parameters for each nt
melting_results <- pos
melting_results$enthalpy_J <- rmelting_out[1, ]
melting_results$entropy_J <- rmelting_out[2, ]
melting_results$temperature_C <- rmelting_out[3, ]


# save results
melting_10bp <- melting_results
save(melting_10bp, file = paste0(save_dir, "/Melting_10bp.RData"))