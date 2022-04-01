# info --------------------------------------------------------------------
# purpose: fcalculate coverage
# author: Robert Gnuegge (robert.gnuegge@gmail.com)
# date: 02/13/22
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


# function definitions ====================================================

# calculate coverage
# argument: GAlignments object
# result: GRanges object with coverage (mcol "score")
calc_coverage <- function(GAlignments, rpm_normalization = TRUE){
  out <- GRanges(coverage(x = GAlignments[strand(GAlignments)=="-"]))  # based on library prep, only alignments to - strand are senseful
  strand(out) <- "-"
  if(rpm_normalization){
    out$score <- out$score / length(GAlignments) * 1e6
  }
  return(out)
}


# process all samples =====================================================

# file base paths
BAM_dir <- "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Data/Bazzano2021/BAM"
save_dir <- "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Data/Bazzano2021/Coverage"
plot_dir <- "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Figures/16_Bazzano2021/Mapping_edit_distances"

# initialize data.frames for result recording
mapping_stats <- data.frame()

# iterate through samples
samples <- paste0("exo1_sgs1_T", c(0, 35, 60, 75, 90))
for(sample in samples){
  
  cat("\n\nProcessing ", sample, "...", sep = "")
  
  # read BAM file -----------------------------------------------------------
  cat("\nReading BAM file...")
  tmp <- readGAlignments(file = paste0(BAM_dir, "/", sample, ".bam"), param = ScanBamParam(tag = c("AS", "XS", "NM")))
  # AS: alignment score (max. read length * match bonus [--ma, default: 2])
  # XS: alignment score for 2nd best alignment
  # NM: edit distance
  mapped <- length(tmp)  # for mapping statistics
  cat(mapped, "alignments read.")
  
  # how many alignments to alleles? -----------------------------------------
  dsb <- sum(seqnames(tmp) == "dsb")
  ctl <- sum(seqnames(tmp) == "ctl")
  cat("\n", dsb, " alignments to dsb allele (", round(100 * dsb / length(tmp), 2), "%).", sep = "")
  cat("\n", ctl, " alignments to ctl allele (", round(100 * ctl / length(tmp), 2), "%).", sep = "")
  
  # plot edit distance distribution -----------------------------------------
  cat("\nPlotting edit distance distribution...")
  pdf(file = "tmp.pdf", width=3, height=2.5)
  par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(1.6, 0.6, 3.6, 2.1), las = 1, tcl = -0.3, mgp = c(2.5, 0.6, 0))
    # plot fractions of alignments with edit distance = 0, 1, 2, ..., 5, >5
    h <- hist(x = mcols(tmp)$NM, breaks = 0:max(mcols(tmp)$NM), plot = FALSE)  # use hist function to calculate densities (fractions)
    bp <- barplot(height = c(h$density[1:6], sum(h$density[7:length(h$density)])), ylim = c(0, 1),
                  names.arg = c(as.character(0:5), ">5"), ylab = "Fraction of Alignments", xlab = "Edit Distance")  # plot density values as bar plot
    # add edit distance <= 3 threshold
    th <- mean(bp[4:5, ])
    abline(v = th, col = "red")
    text(x = 0.95 * th, y = 1, adj = c(1,1), xpd = TRUE, col = "red", labels = paste0(round(x = 100 * ecdf(mcols(tmp)$NM)(3), digits = 2), "%"))
  dev.off()
  GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/", sample, ".pdf"))
  
  # calculate coverage ----------------------------------------------------
  cat("\nCalculating coverage...")
  tmp <- resize(x = GRanges(tmp), width = 1)  # only keep 5' ends
  # save RPM normalized coverage
  tmp_coverage <- calc_coverage(GAlignments = tmp, rpm_normalization = TRUE)
  # transfer coverage from ctrl to dsb after background subtraction
  ctl_coverage <- as_nt_resolved_GRanges(tmp_coverage[seqnames(tmp_coverage) == "ctl"])
  dsb_coverage <- as_nt_resolved_GRanges(tmp_coverage[seqnames(tmp_coverage) == "dsb"])
  homology <- GRanges(seqnames = c("dsb", "ctl"), ranges = IRanges(start = 1, end = 187))
  heterology <- GRanges(seqnames = c("dsb", "ctl"), ranges = IRanges(start = 188, end = 525))
  # subtract background
  ctl_unique <- subsetByOverlaps(x = ctl_coverage, ranges = heterology)
  lm <- lm(score ~ start, data = ctl_unique)
  bg <- predict(object = lm, newdata = data.frame(start = 1:716))
  idx <- which(ctl_coverage$score != 0)
  ctl_coverage$score[idx] <- ctl_coverage$score[idx] - bg[idx]
  ctl_coverage$score[ctl_coverage$score < 0] <- 0
  idx <- which(dsb_coverage$score != 0)
  dsb_coverage$score[idx] <- dsb_coverage$score[idx] - bg[idx]
  dsb_coverage$score[dsb_coverage$score < 0] <- 0
  ctl_shared <- subsetByOverlaps(x = ctl_coverage, ranges = homology)
  hits <- findOverlaps(query = homology, subject = dsb_coverage)
  dsb_coverage[subjectHits(hits)]$score <- dsb_coverage[subjectHits(hits)]$score + ctl_shared$score
  hits <- findOverlaps(query = homology, subject = dsb_coverage)
  ctl_coverage[subjectHits(hits)]$score <- ctl_coverage[subjectHits(hits)]$score - ctl_shared$score
  corrected_coverage <- c(dsb_coverage, ctl_coverage)
  assign(x = sample, value = corrected_coverage)
  
  # save unnormalized coverage
  tmp_coverage <- calc_coverage(GAlignments = tmp, rpm_normalization = FALSE)
  # transfer coverage from ctrl to dsb after background subtraction
  ctl_coverage <- as_nt_resolved_GRanges(tmp_coverage[seqnames(tmp_coverage) == "ctl"])
  dsb_coverage <- as_nt_resolved_GRanges(tmp_coverage[seqnames(tmp_coverage) == "dsb"])
  homology <- GRanges(seqnames = c("dsb", "ctl"), ranges = IRanges(start = 1, end = 187))
  heterology <- GRanges(seqnames = c("dsb", "ctl"), ranges = IRanges(start = 188, end = 525))
  # subtract background
  ctl_unique <- subsetByOverlaps(x = ctl_coverage, ranges = heterology)
  lm <- lm(score ~ start, data = ctl_unique)
  bg <- predict(object = lm, newdata = data.frame(start = 1:716))
  idx <- which(ctl_coverage$score != 0)
  ctl_coverage$score[idx] <- round(ctl_coverage$score[idx] - bg[idx])
  ctl_coverage$score[ctl_coverage$score < 0] <- 0
  idx <- which(dsb_coverage$score != 0)
  dsb_coverage$score[idx] <- round(dsb_coverage$score[idx] - bg[idx])
  dsb_coverage$score[dsb_coverage$score < 0] <- 0
  ctl_shared <- subsetByOverlaps(x = ctl_coverage, ranges = homology)
  hits <- findOverlaps(query = homology, subject = dsb_coverage)
  dsb_coverage[subjectHits(hits)]$score <- dsb_coverage[subjectHits(hits)]$score + ctl_shared$score
  hits <- findOverlaps(query = homology, subject = dsb_coverage)
  ctl_coverage[subjectHits(hits)]$score <- ctl_coverage[subjectHits(hits)]$score - ctl_shared$score
  corrected_coverage <- c(dsb_coverage, ctl_coverage)
  assign(x = paste0(sample, "_unnormalized"), value = corrected_coverage)
  
  # record mapping statistics data ------------------------------------------
  mapping_stats <- rbind(mapping_stats,
                         data.frame(sample = sample,
                                    reads_mapped = mapped,
                                    dsb_alignments = dsb,
                                    ctl_alignments = ctl))

}


# save data ---------------------------------------------------------------
save(list = samples, file = paste0(save_dir, "/exo1_sgs1_coverage.RData"))
save(list = paste0(samples, "_unnormalized"), file = paste0(save_dir, "/exo1_sgs1_coverage_unnormalized.RData"))

write.table(x = mapping_stats, file = paste0(save_dir, "/exo1_sgs1_mapping_stats.txt"), row.names = FALSE)
