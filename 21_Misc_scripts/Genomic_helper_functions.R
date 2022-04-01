# info --------------------------------------------------------------------
# purpose: helper functions to process and plot genomic data
# author: Robert Gnuegge (robert.gnuegge@gmail.com)
# date: 12/16/21
# version: 0.1


# find by intersection ----------------------------------------------------
# argument: GRanges object
# result: GRanges object
# note: based on https://support.bioconductor.org/p/9141415/#9141492
subsetByIntersect <- function(subject, query){
  hits <- findOverlaps(subject, query)
  out <- pintersect(subject[queryHits(hits)], query[subjectHits(hits)])
  mcols(out) <- mcols(out)[names(mcols(out)) != "hit"]
  return(out)
}


# make GRanges nt-resolved ------------------------------------------------
# argument: GRanges object 
# result: GRanges object
as_nt_resolved_GRanges <- function(GRanges){
  out <- GPos(GRanges)
  mcols(out) <- rep(mcols(GRanges), width(GRanges))  # add back metadata columns, which are dropped when converting to GPos (bug?)
  out <- GRanges(out)  # convert back to GRanges, as in anyway no storage difference for my data, but missing functionality of GPos compared to GRanges objects
  return(out) 
}


# make GRanges with fw/rev score mcols ------------------------------------
# argument: GRanges object
# result: GRanges object
as_GRanges_with_fw_rev_scores <- function(GRanges, rev_negative = FALSE, rm_0_coverage = FALSE){
  GRanges_nt_resolved <- as_nt_resolved_GRanges(GRanges)
  fw_rev_intersect <- intersect(GRanges_nt_resolved[strand(GRanges_nt_resolved) == "+"], GRanges_nt_resolved[strand(GRanges_nt_resolved) == "-"], ignore.strand = TRUE)
  GRanges_nt_resolved <- subsetByOverlaps(x = GRanges_nt_resolved, ranges = fw_rev_intersect)  # remove ranges unique to one strand
  tmp <- GRanges_nt_resolved
  strand(tmp) <- "*"
  tmp <- reduce(tmp)  # removes duplicated ranges
  tmp <- as_nt_resolved_GRanges(tmp)
  tmp$fw <- GRanges_nt_resolved[strand(GRanges_nt_resolved) == "+"]$score
  tmp$rev <- GRanges_nt_resolved[strand(GRanges_nt_resolved) == "-"]$score * ifelse(test = rev_negative, yes = -1, no = 1)
  if(rm_0_coverage){
    tmp <- tmp[!(tmp$fw == 0 & tmp$rev == 0)]
  }
  return(tmp)
}


# apply Hanning smoother to GRanges score data ----------------------------
# argument: GRanges object
# result: GRanges object
# note: artifacts will be produced when neighboring entries are from different loci
as_smoothed_GRanges <- function(GRanges, hanning_window_size){
  require(zoo)
  hanning_window <- function(n = 51){
    if(n == 1){
      c <- 1
    }else{
      n <- n - 1
      c <- 0.5 - 0.5 * cos(2 * pi * (0:n) / n)
    }
    return(c / sum(c))
  }
  GRanges[strand(GRanges) == "+"]$score <- rollapply(data = GRanges[strand(GRanges) == "+"]$score, width = hanning_window_size, FUN = weighted.mean, w = hanning_window(n = hanning_window_size), fill = NA)
  GRanges[strand(GRanges) == "-"]$score <- rollapply(data = GRanges[strand(GRanges) == "-"]$score, width = hanning_window_size, FUN = weighted.mean, w = hanning_window(n = hanning_window_size), fill = NA)
  GRanges[strand(GRanges) == "*"]$score <- rollapply(data = GRanges[strand(GRanges) == "*"]$score, width = hanning_window_size, FUN = weighted.mean, w = hanning_window(n = hanning_window_size), fill = NA)
  return(GRanges)
}


# S. cerevisiae chromosome information ------------------------------------
SacCer_chromosomes_df <- data.frame(chr = c(paste0("chr", as.roman(1:16)), "chrM", "micron"),
                                    length = c(230218, 813184, 316620, 1531933, 576874, 270161, 1090940, 562643, 439888, 745751, 666816, 1078177, 924431, 784333, 1091291, 948066, 85779, 6392),
                                    CEN = c(151524, 238265, 114443, 449766, 152046, 148569, 496979, 105645, 355687, 436366, 440187, 150888, 268090, 628817, 326643, 556015, NA, NA))

SacCer_chromosomes <- GRanges(seqnames = SacCer_chromosomes_df$chr, ranges = IRanges(start = 1, end = SacCer_chromosomes_df$length))
seqlengths(SacCer_chromosomes) <- SacCer_chromosomes_df$length


# -------------------------------------------------------------------------
# get genomic regions around DSBs
# argument: GRanges object, integer
# result: GRanges object
DSB_regions <- function(DSBs, region_width = 4000, up_rev_down_fw = FALSE){
  if (up_rev_down_fw){
    fw <- DSBs
    strand(fw) <- "+"
    fw <- resize(x = fw, width = 0.5 * region_width, fix = "start")
    rev <- DSBs
    strand(rev) <- "-"
    rev <- shift(x = rev, shift = -1)
    rev <- resize(x = rev, width = 0.5 * region_width, fix = "start")
    out <- sort(x = c(fw, rev), ignore.strand = TRUE)
    return(out)
  }else{
    out <- resize(x = DSBs, width = region_width, fix = "center")
    return(out)
  }
}