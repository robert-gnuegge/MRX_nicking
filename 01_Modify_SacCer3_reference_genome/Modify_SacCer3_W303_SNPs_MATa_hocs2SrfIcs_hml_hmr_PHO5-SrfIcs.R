# info --------------------------------------------------------------------
# purpose: modify S288C reference genome to improve mapping of W303-derived reads; add strain-specific mutations
# author: Robert Gnuegge (robert.gnuegge@gmail.com)
# date: 01/27/22
# version: 1.0

# preamble ----------------------------------------------------------------

# set wd to this file's location
wd_path <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(wd_path)

# libraries
library(GenomicRanges)
library(Biostrings)
library(rtracklayer)


# load S288C reference genome ---------------------------------------------
# sequence downloaded from http://sgd-archive.yeastgenome.org/sequence/S288C_reference/genome_releases/ on 05/11/20

genome <- import(con = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Data/S_cerevisiae_reference_genomes/S288C_R64-2-1_W303_SNPs_MATa_hocs2SrfIcs_hml_hmr.fasta")
names(genome)
seqinfo(genome)

matchPattern(pattern = DNAString(x = "GCCCGGGC"), subject = genome[["chrIII"]])


# add PHO5-SrfIcs modification --------------------------------------------

PHO5_mod <- DNAString(x = "TCTATGCGAAAACGTGCTAATTTAATATATATTTCTTTGTGCAGACAAAGAAAAAGCGCTATGAACCTTTTACCTTCGTTTGAAGTTTATAGACGACGTCCGCTTACATGAGAAAACATTAAAGCAGCGCACATAAGATGACTTCCAAATACGTTGACATATTTGCGCATTCTTGTTGAATAgccCgggCGGATCCTTcacATGCCGGTCTCCGGCGTCCTGTAAAGAGAGCGTGCGACACGCCGCTATTAGCGTGGGGACTGACATCAGGTCAACATTATTATTGTAGCGAGCTACTTTCGGTCAGTAAATAGAACATATGTAGAGATACAAGCGATTATTAGAAGGAGCGGACCTACAACAAAGTTGGGTCACGT")

# find corresponding coordinates in genome
start <- start(matchPattern(pattern = PHO5_mod[1:50], subject = genome[["chrII"]]))
end <- end(matchPattern(pattern = PHO5_mod[(length(PHO5_mod)-50):length(PHO5_mod)], subject = genome[["chrII"]]))

# replace
genome[["chrII"]][start:end] <- PHO5_mod

# are there two SrfIcs now on chrII?
matchPattern(pattern = "GCCCGGGC", subject = genome[["chrII"]])


# export modified reference genome to fast file ---------------------------

export(object = genome, con = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Data/S_cerevisiae_reference_genomes/S288C_R64-2-1_W303_SNPs_MATa_hocs2SrfIcs_hml_hmr_PHO5-SrfIcs.fasta")
