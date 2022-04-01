# info --------------------------------------------------------------------
# purpose: filter alignments, calculate S1-seq coverage, and extract mapping stats
# author: Robert Gnuegge (robert.gnuegge@gmail.com)
# date: 01/28/22
# version: 1.1


# preamble ----------------------------------------------------------------

# set wd to this file's location
wd_path <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(wd_path)

# libraries
library(GenomicAlignments)

# load helper files and functions
source(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Analysis_scripts/Misc/Genomic_helper_functions.R")
source(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Analysis_scripts/Misc/Misc_helper_functions.R")


# process all samples =====================================================

# file base paths
BAM_dir <- "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Data/BAM/MNase-seq"
save_dir <- "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Data/MNase-seq/LSY4753-3B_4754-3A"
plot_dir_edit <- "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Figures/12_MNase-seq_edit_distance_distributions/LSY4753-3B_4754-3A"
plot_dir_insert <- "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Figures/12_MNase-seq_insert_size_distributions/LSY4753-3B_4754-3A"

# initialize data.frames for result recording
mapping_stats <- data.frame()

# iterate through samples
samples <- c("LSY4753-3B_0", "LSY4753-3B_2", "LSY4754-3A_0", "LSY4754-3A_2")
for(sample in samples){
  
  cat("\n\nProcessing ", sample, "...", sep = "")
  
  # read BAM file -----------------------------------------------------------
  cat("\nReading BAM file...")
  tmp <- readGAlignmentPairs(file = paste0(BAM_dir, "/", sample, "/", sample, ".bam"), param = ScanBamParam(tag = c("AS", "YS", "NM"), what = "isize"))
  # AS: alignment score for first mate (max. read length * match bonus [--ma, default: 2])
  # YS:alignment score for opposite mate
  # NM: edit distance
  # isize: insert size
  seqlevels(tmp)[seqlevels(tmp) == "mikron"] <- "micron"
  all_mapped <- length(tmp)  # for mapping statistics
  cat(all_mapped, "alignments read.")
  
  # Remove alignments with too large insert size ----------------------------
  cat("\nRemoving alignments with insert size >250 bp...")
  tmp <- tmp[mcols(first(tmp))$isize <= 250]
  mapped <- length(tmp)
  cat(" kept ", mapped, " alignments (",  round(100 * mapped / all_mapped, digits = 2), "%).", sep = "")

  # how do alignments distribute on yeast DNA types? ------------------------
  chrM <- sum(seqnames(tmp) == "chrM")
  micron <- sum(seqnames(tmp) == "micron")
  nuc_genome <- sum(!seqnames(tmp) %in% c("chrM", "micron"))
  cat("\n", nuc_genome, " alignments to nuclear genome (", round(100 * nuc_genome / length(tmp), 2), "%).", sep = "")
  cat("\n", chrM, " alignments to mitochondrial genome (", round(100 * chrM / length(tmp), 2), "%).", sep = "")
  cat("\n", micron, " alignments to 2-micron plasmid (", round(100 * micron / length(tmp), 2), "%).", sep = "")
  
  # plot edit distance distribution -----------------------------------------
  cat("\nPlotting edit distance distribution...")
  pdf(file = "tmp.pdf", width=2.5, height=2.5)
  par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(1.6, 0.6, 3.6, 2.1), las = 1, tcl = -0.3, mgp = c(2.5, 0.6, 0))
  # plot fractions of alignments with edit distance = 0, 1, 2, ..., 5, >5
  edit_dist <- c(mcols(first(tmp))$NM, mcols(last(tmp))$NM)
  h <- hist(x = edit_dist, breaks = 0:max(edit_dist), plot = FALSE)  # use hist function to calculate densities (fractions)
  bp <- barplot(height = c(h$density[1:6], sum(h$density[7:length(h$density)])), ylim = c(0, 1),
                names.arg = c(as.character(0:5), ">5"), ylab = "Fraction of Alignments", xlab = "Edit Distance")  # plot density values as bar plot
  # add edit distance <= 3 threshold
  th <- mean(bp[4:5, ])
  abline(v = th, col = "red")
  text(x = 0.95 * th, y = 1, adj = c(1,1), xpd = TRUE, col = "red", labels = paste0(round(x = 100 * ecdf(edit_dist)(3), digits = 2), "%"))
  dev.off()
  GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir_edit, "/", sample, ".pdf"))
  
  # plot insert size distribution -----------------------------------------
  cat("\nPlotting insert size distribution...")
  pdf(file = "tmp.pdf", width=2.75, height=2.5)
  par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(2, 0.2, 4.1, 1.6), las = 1, tcl = -0.3, mgp = c(2, 0.5, 0))
  insert_size <- abs(mcols(first(tmp))$isize)
  h <- hist(x = insert_size, breaks = 100, plot = FALSE)  # to calculate ylim of histogram
  hist(x = insert_size, breaks = 100, xlim = c(50, 250), probability = TRUE, 
       ylim = c(0, 1.2 * max(h$density)), xlab = "Insert size [bp]", ylab = NA, main = NA)
  title(ylab = "Probability", line = 3)
  med <- median(insert_size)
  segments(x0 = med, y0 = 0, y1 = 1.1 * max(h$density), col = "red")
  text(x = med, y = 1.1 * max(h$density), labels = med, pos = 3, offset = 0.4, col = "red")
  dev.off()
  GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir_insert, "/", sample, ".pdf"))
  
  # Convert to GRanges ------------------------------------------------------
  # conversion to GRanges is necessary for trimming (setting start and end is not possible for GAlignmentPairs class)
  # and for coverage calculation to consider the complete insert sequence
  cat("\nConverting to GRanges...")
  tmp <- GRanges(tmp)
  
  # Adjust to median insert size = 147 bp -----------------------------------
  cat("\nMedian insert size is ", med , " bp.", sep = "")
  tmp_trimmed <- tmp
  # trim
  if(med != 147){
      cat(" Adjusting to median insert size = 147 bp...")
      trim_length <- (med - 147)
      if(trim_length %% 2 != 0){  # in case trim_length is uneven
        # add ceiling of half trim_length randomly to left or right end of insert
        trim_left <- ifelse(test = sample(x = c(0, 1), size = length(tmp), replace = TRUE), yes = ceiling(0.5 * trim_length), no = floor(0.5 * trim_length))
      }else{
        trim_left <- 0.5 * trim_length
      }
      trim_right <- trim_length - trim_left
      start(tmp_trimmed) <- start(tmp_trimmed) + trim_left
      end(tmp_trimmed) <- end(tmp_trimmed) - trim_right
      tmp_trimmed <- trim(tmp_trimmed)  # in case of out-of-bond ranges
  }

  # calculate MNase-seq coverage --------------------------------------------
  cat("\nCalculating coverage...")
  tmp_coverage <- GRanges(coverage(tmp))
  tmp_coverage$score <- tmp_coverage$score / length(tmp_coverage) * 1e6  # convert to RPM 
  assign(x = paste0(dash_to_underscore(sample), "_MNase_seq"), value = tmp_coverage)
  tmp_trimmed_coverage <- GRanges(coverage(tmp_trimmed))
  tmp_trimmed_coverage$score <- tmp_trimmed_coverage$score / length(tmp_trimmed_coverage) * 1e6  # convert to RPM 
  assign(x = paste0(dash_to_underscore(sample), "_MNase_seq_trimmed"), value = tmp_trimmed_coverage)
  
  # record mapping statistics data ------------------------------------------
  mapping_stats <- rbind(mapping_stats,
                         data.frame(sample = sample,
                                    reads_mapped = mapped,
                                    nuc_genome = nuc_genome,
                                    mt_genome = chrM,
                                    micron = micron,
                                    insert_size_median = med))
  
}


# save data ---------------------------------------------------------------
save(list = paste0(dash_to_underscore(samples), "_MNase_seq"), file = paste0(save_dir, "/LSY4753-3B_4754-3A.RData"))
save(list = paste0(dash_to_underscore(samples), "_MNase_seq_trimmed"), file = paste0(save_dir, "/LSY4753-3B_4754-3A_trimmed.RData"))
write.table(x = mapping_stats, file = paste0(save_dir, "/LSY4753-3B_4754-3A_mapping_stats.txt"), row.names = FALSE)