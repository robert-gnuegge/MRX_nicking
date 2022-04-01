# info --------------------------------------------------------------------
# purpose: define SrfI cut sites (and windows around these) in the S. cerevisiae genome
# author: Robert Gnuegge (robert.gnuegge@gmail.com)
# date: 12/23/21
# version: 0.1

# -------------------------------------------------------------------------
# SrfIcs as data.frame (SrfIcs on chrIII:200753 [MAT] is not endogenous, but engineered)
SrfIcs_df <- data.frame(seqnames = paste0("chr", as.roman(c(2, 3, 4, 5, 5, 7, 7, 7, 9, 10, 12, 13, 13, 13, 13, 14, 15, 15, 15, 15, 16))),
                        pos = c(256173, 200753, 370477, 123366, 399646, 398052, 517001, 642694, 22006, 360907, 498747, 129799, 566584, 664938, 676597, 301414, 27760, 370687, 756594, 1039563, 431983),
                        stringsAsFactors = FALSE)


# -------------------------------------------------------------------------
# SrfIcs as GRanges object
SrfIcs <- GRanges(seqnames = SrfIcs_df$seqnames, ranges = IRanges(start = SrfIcs_df$pos, width = 1))
S_cerevisiae_seqinfo <- Seqinfo(seqnames = c(paste0("chr", as.roman(1:16)), "chrM", "micron"), 
                                seqlengths = c(230218, 813184, 316620, 1531933, 576874, 270161, 1090940, 562643, 439888, 
                                               745751, 666816, 1078177, 924431, 784333, 1091291, 948066, 85779, 6392))
seqlevels(SrfIcs) <- seqlevels(S_cerevisiae_seqinfo)
seqinfo(SrfIcs) <- S_cerevisiae_seqinfo
