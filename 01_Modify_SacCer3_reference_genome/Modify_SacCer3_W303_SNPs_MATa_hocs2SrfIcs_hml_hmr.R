# info --------------------------------------------------------------------
# purpose: modify S288C reference genome to improve mapping of W303-derived reads; add strain-specific mutations
# author: Robert Gnuegge (robert.gnuegge@gmail.com)
# date: 12/16/21
# version: 1.0

# preamble ----------------------------------------------------------------

# set wd to this file's location
wd_path <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(wd_path)

# libraries
library(GenomicRanges)
library(Biostrings)
library(rtracklayer)

source(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Analysis_scripts/Misc/S_cerevisiae_SrfI_cut_sites.R")


# load S288C reference genome ---------------------------------------------
# sequence downloaded from http://sgd-archive.yeastgenome.org/sequence/S288C_reference/genome_releases/ on 05/11/20

genome <- import(con = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Resources/S_cerevisiae_genome/S288C_reference_sequence_R64-2-1_20150113.fasta")

# change names
names(genome)
names(genome) <- c(paste0("chr", as.roman(1:16)), "chrM")
names(genome)


# add 2 micron plasmid ----------------------------------------------------
# sequence downloaded from https://www.ncbi.nlm.nih.gov/nuccore/CM001822.1 on 05/11/20

micron <- import(con = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Resources/S_cerevisiae_genome/2micron_W303.fasta")

# the 2 micron plasmid is circular
# to support read mapping near its 3' end, add a sequence stretch from its 5' end
read_length <- 75
micron$micron <- c(micron$micron, micron$micron[1:(read_length - 1)])

# add to genome
genome <- DNAStringSet(x = c(as.character(genome), as.character(micron)))


# inject W303 SNPs --------------------------------------------------------
# SNPs downloaded from https://www.g3journal.org/highwire/filestream/485032/field_highwire_adjunct_files/0/FileS1.xlsx on 05/11/20
# (Matheson et al., 2017)

SNPs <- read.csv(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Resources/S_cerevisiae_genome/W303_SNPs.csv", header = TRUE, stringsAsFactors = FALSE)

# change chrmt to chrM
SNPs$Chromosome[SNPs$Chromosome == "chrmt"] <- "chrM"

# check that the S288C nts in the SNPs file are consistent with the genome sequence
all(sapply(X = names(genome)[1:17], FUN = function(c) {as.character(genome[[c]][SNPs$Position.in.S288c[SNPs$Chromosome == c]]) == paste0(SNPs$Base.in.S288c[SNPs$Chromosome == c], collapse = "")}))

# substitute with W303 SNPs
for(c in c(paste0("chr", as.roman(1:16)), "chrM")){
  genome[[c]][SNPs$Position.in.S288c[SNPs$Chromosome == c]] <- DNAString(paste0(SNPs$Base.in.W303[SNPs$Chromosome == c], collapse = ""))
}

# check that SNPs have been injected
all(sapply(X = names(genome)[1:17], FUN = function(c) {as.character(genome[[c]][SNPs$Position.in.S288c[SNPs$Chromosome == c]]) == paste0(SNPs$Base.in.W303[SNPs$Chromosome == c], collapse = "")}))


# change MATalpha to MATa -------------------------------------------------
# the strains used for S1-seq library preparation are all MATa

# correct HMRa1 stk mutation (see https://www.yeastgenome.org/locus/S000029699#paragraph for details)
matchPattern(pattern = DNAString("CGCAACAGTAAAATTTTATAA"), subject = genome[["chrIII"]])
matchPattern(pattern = DNAString("CGCAACAGTATAATTTTATAA"), subject = genome[["chrIII"]])
# already corrected during the W303 SNP injection above

# define MATalpha Z region (coordinates from https://www.yeastgenome.org/locus/S000029699)
Z1 <- genome[["chrIII"]][200851:201089]

# find Z regions in HML, MAT, and HMR
all.Z1 <- matchPattern(pattern = Z1, subject = genome[["chrIII"]], max.mismatch = 1)

# turns out that there is a 1 bp difference between HML/MAT and HMR (C->T)
mismatch(pattern = Z1, x = all.Z1)

# define MATalpha X region
X <- genome[["chrIII"]][199401:200103] 

# find X regions in HML, MAT, and HMR
all.X <- matchPattern(pattern = X, subject = genome[["chrIII"]], max.mismatch = 2, with.indels = TRUE)

# MATalpha has a 2 bp deletion at position 199512:199513 (in a poly-TA stretch)
mismatch(pattern = X, x = all.X)  # gives mismatches for HML and HMR at relative position > 111

# let's only use the first 111 bp of the MATalpha X
X <- genome[["chrIII"]][199401:199511] 

# find X regions in HML, MAT, and HMR
all.X <- matchPattern(pattern = X, subject = genome[["chrIII"]], max.mismatch = 0)

# MATalpha has a 2 bp deletion at position 199512:199513 (in a poly-TA stretch)
mismatch(pattern = X, x = all.X)  # gives mismatches for HML and HMR at relative position > 111

# check size difference of regions
end(all.Z1[2]) - start(all.X[2]) - ( end(all.Z1[3]) - start(all.X[3]) )
# MATalpha X-Z region is 103 bp longer than MATa X-Z region

# replace MATalpha with HMR and add 103 Ns downstream of TAF1 gene (ca. 5 kb downstream of HO cut site)
TAF1.end <- 206260  # upstream of YCR043C, which is downstream of TAF1 (https://www.yeastgenome.org/locus/YCR043C)
MATalpha <- genome[["chrIII"]][start(all.X[2]):TAF1.end]

HMR <- genome[["chrIII"]][start(all.X[3]):end(all.Z1[3])]
TAF1 <- genome[["chrIII"]][(end(all.Z1[2]) + 1):TAF1.end]
Ns <- DNAString(x = paste0(rep("N", 103), collapse = ""))
replacement <- c(HMR, TAF1, Ns)

# sanity check
length(MATalpha) == length(replacement)

# replace
genome[["chrIII"]][start(all.X[2]):TAF1.end] <- replacement

# check if HMR sequence was copied to MAT
matchPattern(pattern = HMR, subject = genome[["chrIII"]])


# change MATa HO cut site to SrfI cut site --------------------------------

# find HO cut sites
HOcs <- matchPattern(pattern = DNAString(x = "GCAACAGTATAATTTTATAAACCCT"), subject = genome[["chrIII"]])
# the middle hit is MAT

# replace with SrfI cut site
genome[["chrIII"]][start(HOcs[2]):(start(HOcs[2]) + 7)] <- DNAString(x = "GCCCGGGC")

# sanity check
genome[["chrIII"]][(start(HOcs[2]) - 10):(start(HOcs[2]) + 17)]



# HML deletion ------------------------------------------------------------
# the HML deletion was generated by C. Zierhut (https://www.ncbi.nlm.nih.gov/pubmed/18511906 and personal communication)
# HML was replaced by URA3, which was then in turn replaced by a PCR fragment derived from pRS306 (ori) using the following primers
OCZ330 <- DNAString(x = "gacaagtagcgcagttattttcttattttcatCATCAGAGCAGATTGTACTGAGAGTGCAtaattgcgttgcgctcactg")
OCZ331 <- reverseComplement(x = DNAString(x = "CAAATTTACTAGCTcgagagttagttaatccttactaagtgaagaaaagcaactatagtcTCAAGAACTCTGTAGCACCG"))

# find replaced genome region
HML.start <- start(matchPattern(pattern = OCZ330[1:32], subject = genome[["chrIII"]]))
HML.end <- end(matchPattern(pattern = OCZ331[21:80], subject = genome[["chrIII"]]))
HML <- genome[["chrIII"]][HML.start:HML.end]

# find amplified pRS306 region
# pRS306 sequence downloaded from https://www.ncbi.nlm.nih.gov/nuccore/U03438.1/ on 05/11/20
pRS306 <- import(con = "~/Research/Resources/pRS306.fasta")[[1]]
ori.start <- start(matchPattern(pattern = OCZ330[61:80], subject = pRS306))
ori.end <- end(matchPattern(pattern = OCZ331[1:20], subject = pRS306))
ori <- pRS306[ori.start:ori.end]
# expand with primer overhangs
replacement <- c(OCZ330[1:60], ori, OCZ331[21:80])

# length of replacement is 3495 bp shorter than replaced region
length(HML) - length(replacement)
# let's fill up with Ns in the center of the replacement
replacement <- c(replacement[1:(round(0.5 * length(replacement)) - 1)],
                 DNAString(x = paste0(rep("N", 3495), collapse = "")),
                 replacement[round(0.5 * length(replacement)):length(replacement)])

# sanity check
length(HML) == length(replacement)

# replace
genome[["chrIII"]][HML.start:HML.end] <- replacement

# check if ori parts were inserted into genome
matchPattern(pattern = ori[1:floor(0.5 * length(ori) - 1)], subject = genome[["chrIII"]])
matchPattern(pattern = ori[ceiling(0.5 * length(ori) + 1):length(ori)], subject = genome[["chrIII"]])
matchPattern(pattern = HML, subject = genome[["chrIII"]])


# HMR deletion ------------------------------------------------------------
# the HMR deletion was generated by C. Zierhut (https://www.ncbi.nlm.nih.gov/pubmed/18511906 and personal communication).
# HMR was replaced by URA3, which was then in turn replaced by a PCR fragment derived from pRS306 (bla) using the following primers
OCZ332 <- DNAString(x = "gatttaagcgtgcgtGAAGATAACACTACAATCCATTTTAAAGCAACATCCACATTGAGTcgttcatccatagttgcctg")
OCZ333 <- reverseComplement(x = DNAString(x = "ggtatatagaatataCTAGAAGTTCTCCTCGAGGATATAGGAATCCTCAAAAGGGAATCTttcaacatttccgtgtcgcc"))

# find replaced genome region
HMR.start <- start(matchPattern(pattern = OCZ332[1:60], subject = genome[["chrIII"]]))
HMR.end <- end(matchPattern(pattern = OCZ333[21:80], subject = genome[["chrIII"]]))
HMR <- genome[["chrIII"]][HMR.start:HMR.end]

# find amplified pRS306 region
pRS306 <- import(con = "~/Research/Resources/pRS306.fasta")[[1]]
bla.start <- start(matchPattern(pattern = OCZ332[61:80], subject = pRS306))
bla.end <- end(matchPattern(pattern = OCZ333[1:20], subject = pRS306))
bla <- pRS306[bla.start:bla.end]

# expand with primer overhangs
replacement <- c(OCZ332[1:60], bla, OCZ333[21:80])

# length of replacement is 3341 bp shorter than replaced region
length(HMR) - length(replacement)
# let's fill up with Ns in the center of the replacement
replacement <- c(replacement[1:(round(0.5 * length(replacement)) - 1)],
                 DNAString(x = paste0(rep("N", 3341), collapse = "")),
                 replacement[round(0.5 * length(replacement)):length(replacement)])

# sanity check
length(HMR) == length(replacement)

# replace
genome[["chrIII"]][HMR.start:HMR.end] <- replacement

# check if bla parts were inserted into the genome
matchPattern(pattern = bla[1:floor(0.5 * length(bla) - 1)], subject = genome[["chrIII"]])
matchPattern(pattern = bla[ceiling(0.5 * length(bla) + 1):length(bla)], subject = genome[["chrIII"]])
matchPattern(pattern = HMR, subject = genome[["chrIII"]])


# check if sequences around SrfI cut sites correspond to W303 reference sequence --------

# load W303 reference genome 
# downloaded from https://www.ncbi.nlm.nih.gov/assembly/GCA_002163515.1 on 05/11/20
# (Matheson et al., 2017) 
W303 <- import(con = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Resources/S_cerevisiae_genome/W303_Matheson2017.fasta")

# change names
names(W303) <- c(paste0("chr", as.roman(1:16)), "chrM", "micron", paste0("chrXII-rep_seq-", 1:3))

# function to find W303 sequence corresponding to S288C sequence, align, and return mismatches and indels
FindAndAlign <- function(chr, coord, nt.window = 4000, genome1 = genome, genome1name = "S288C", genome2 = W303, genome2name = "W303", return.alignment = FALSE){
  
  cat("\n\n------------------------------------------")
  cat("\nSrfIcs at ", chr, ":", coord, sep = "")
  
  coords <- coord + (-0.5 * nt.window):(0.5 * nt.window)
  
  cat("\nAligning", genome1name, "and", genome2name, "at", chr, ":", min(coords), "-", max(coords), "...")
  
  # extract sequence 1
  seq1 <- genome1[[chr]][coords]
  
  # find and extract sequence 2
  hit <- matchPattern(pattern = seq1, subject = genome2[[chr]], max.mismatch = 10, with.indels = TRUE)
  
  if(length(hit) == 0){
    cat("\nSelected ", genome1name, " sequence could not be found in ", genome2name, ".")
  }else{
    
    seq2 <- genome2[[chr]][start(hit):end(hit)] 
    
    # align both sequences
    algnmnt <- pairwiseAlignment(pattern = seq1, subject = seq2)
    
    # analyze alignment
    cat("\n\nNumber of mismatches:", nmismatch(algnmnt))
    
    # Insertions
    ins <- insertion(algnmnt)[[1]]
    
    if(length(ins) == 0){
      cat("\n\nNumber of insertions: 0")
    }else{
      cat("\n\nNumber of insertions:", length(ins))
      
      for(n in 1:length(ins)){
        cat("\n\nInsertion", n)
        cat("\n", genome1name, ": ", as.character(seq1[(start(ins)[n]-5):(end(ins)[n]+5)]))
        cat("\n", genome2name, ":  ", as.character(seq2[(start(ins)[n]-5):(end(ins)[n]+5)]))
        
      }
      
    }
    
    # Deletions
    dels <- deletion(algnmnt)[[1]]
    
    if(length(dels) == 0){
      cat("\n\nNumber of deletions: 0")
    }else{
      cat("\n\nNumber of deletions:", length(dels))
      
      for(n in 1:length(dels)){
        cat("\n\nDeletion", n)
        cat("\n", genome1name, ": ", as.character(seq1[(start(dels)[n]-5):(end(dels)[n]+5)]))
        cat("\n", genome2name, ":  ", as.character(seq2[(start(dels)[n]-5):(end(dels)[n]+5)]))
        
      }
      
    }
    
    if(return.alignment){
      return(algnmnt)
    }
    
  }
  
}

# channel output to file and console
sink(file = "Mod_S288C_W303_differences_around_SrfIcs.txt", split = TRUE)

# check all SrfIcs
for(n in 1:nrow(SrfIcs_df)){
  FindAndAlign(chr = SrfIcs_df$seqnames[n], coord = SrfIcs_df$pos[n])
}

# close output file connection
sink()


# evaluate region around SrfIcs at chrIV:370477 ---------------------------
# sequence around SrfIcs  at chrIV:370477 could not be mapped to W303
# let's try nearer to SrfIcs
n <- 2
FindAndAlign(chr = SrfIcs_df$seqnames[n], coord = SrfIcs_df$pos[n], nt.window = 1000)
# successful mapping!
# let's expand the window around the SrfIcs and find out, when mapping becomes unsuccessful
for(nt.window in 1000 + 1:30 * 100){
  tmp <- FindAndAlign(chr = SrfIcs_df$seqnames[n], coord = SrfIcs_df$pos[n], nt.window = nt.window, return.alignment = TRUE)
  if(is.null(tmp)){
    cat("nt.windwo =", nt.window)
    break
  }
}

FindAndAlign(chr = SrfIcs_df$seqnames[n], coord = SrfIcs_df$pos[n], nt.window = 3260)
# probably there is a GCR or a stretch of > 10 mismatches at ca. 1.6 kb from SrfIcs
# seems to be distant enough to allow lots of MRX nick site mapping


# evaluate region around SrfIcs at MAT (chrIII:200908) --------------------
# there were 6 mismatches
# let's inspect them
n <- 21
tmp <- FindAndAlign(chr = SrfIcs_df$seqnames[n], coord = SrfIcs_df$pos[n], return.alignment = TRUE)

# 5/6 mismatches are due to the inserted hocs::SrfIcs
mismatchTable(tmp)
genome[[SrfIcs_df$seqnames[n]]][SrfIcs_df$pos[n] - 2000 + 1835:1855]
# can't figure out where the 6th mutation (SNP) is
SNPs[SNPs$Chromosome == "chrIII" & SNPs$Position.in.S288c %in% (SrfIcs_df$pos[n] + (-2000:2000)), ]
# might be shifted in its location due to the accompanying insertion in this region; 
# it's a single mutation, so let's neglect it


# export modified reference genome to fast file ---------------------------

export(object = genome, con = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Data/S_cerevisiae_reference_genomes/S288C_R64-2-1_W303_SNPs_MATa_hocs2SrfIcs_hml_hmr.fasta")
