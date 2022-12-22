# info --------------------------------------------------------------------
# purpose: plot density plots for flow cytometry data
# author: Robert Gnuegge (robert.gnuegge@gmail.com)
# date: 02/25/22
# version: 1.0


# preamble ----------------------------------------------------------------

# set wd to this file's location
wd_path <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(wd_path)

# load helper files and functions
source(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Analysis_scripts/Misc/Misc_helper_functions.R")

# process data ============================================================

# read extracted fluorescence data
df <- read.table(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Data/Flow_cytometry/LSY4810-41D_4811-11B_4812-5A/Extracted_fluorescence_data.txt", header = TRUE, stringsAsFactors = FALSE)
str(df)

# aggregate data per sample replicate
fluor.modes <- aggregate(formula = cbind(Citrine, mKate2) ~ sample, data = df,
                         FUN = function(x) c(mean = mean(x, na.rm = TRUE), sd = sd(x, na.rm = TRUE), n = length(x)))

# convert into separate columns
fluor.modes <- data.frame(fluor.modes$sample, as.data.frame(fluor.modes$Citrine), as.data.frame(fluor.modes$mKate2), stringsAsFactors = FALSE)
colnames(fluor.modes) <- c("strain", "Citrine.mean", "Citrine.sd", "Citrine.n", "mKate2.mean", "mKate2.sd", "n")
# keep only one "n" column
fluor.modes <- fluor.modes[, !(colnames(fluor.modes) == "Citrine.n")]

# extract replicate number and add as column
fluor.modes$rep <- as.integer(substr(x = fluor.modes$strain, start = 6, stop = 6))
fluor.modes$strain <- substr(x = fluor.modes$strain, start = 1, stop = 4)

# aggregate per sample
avgs <- aggregate(formula = cbind(Citrine.mean, mKate2.mean) ~ strain, data = fluor.modes,
                  FUN = function(x) c(mean = mean(x, na.rm = TRUE), sd = sd(x, na.rm = TRUE)))
# convert into separate columns
avgs <- data.frame(avgs$strain, as.data.frame(avgs$Citrine.mean), as.data.frame(avgs$mKate2.mean), stringsAsFactors = FALSE)

colnames(avgs) <- c("strain", "Citrine.mean", "Citrine.sd", "mKate2.mean", "mKate2.sd")


# plotting ====================================================================

plot_dir <- "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Figures/19_Flow_cytometry/LSY4810-41D_4811-11B_4812-5A"
tmp <- avgs


# plot mKate2 -------------------------------------------------------------
idx <- c(1, 4, 2)

pdf(file = "tmp.pdf", width=1.75, height=2.5)
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(0.3, 0.2, 4, 2.1), tcl = -0.25, mgp = c(3, 0.5, 0), las = 1)

# plot means
bp <- barplot(height = tmp$mKate2.mean[idx], ylim = c(0, max(tmp$mKate2.mean + tmp$mKate2.sd)), ylab = "Fluorescence [AU]")

# add error bars
arrows(x0 = bp[, 1],
       y0 = tmp$mKate2.mean[idx] - tmp$mKate2.sd[idx],
       x1 = bp[, 1],
       y1 = tmp$mKate2.mean[idx] + tmp$mKate2.sd[idx],
       length = 0.03, # length of arrow head
       angle = 90, # angle of arrow head
       code = 3, # to draw arrow head on both ends
       xpd = TRUE)

# x labels
LabelTxt <- c("Background", "No promoter", expression(italic("TDH3")~"promoter"))
text(x = bp[, 1], y = -0.05 * par("usr")[4], labels = LabelTxt, srt = 45, xpd = TRUE, adj = c(1,1))

dev.off()

GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/mKate2.pdf"))

# -----------------------------------------------------------------------------
# Citrine
idx <- c(1, 2, 3)

pdf(file = "tmp.pdf", width=1.75, height=2.5)
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(0.3, 0.2, 4, 2.1), tcl = -0.25, mgp = c(3, 0.5, 0), las = 1)

# plot means
bp <- barplot(height = tmp$Citrine.mean[idx], ylim = c(0, max(tmp$Citrine.mean + tmp$Citrine.sd)), ylab = "Fluorescence [AU]")

# add error bars
arrows(x0 = bp[, 1],
       y0 = tmp$Citrine.mean[idx] - tmp$Citrine.sd[idx],
       x1 = bp[, 1],
       y1 = tmp$Citrine.mean[idx] + tmp$Citrine.sd[idx],
       length = 0.03, # length of arrow head
       angle = 90, # angle of arrow head
       code = 3, # to draw arrow head on both ends
       xpd = TRUE)

# x labels
LabelTxt <- c("Background", "No promoter", expression(italic("TDH3")~"promoter"))
text(x = bp[, 1], y = -0.05 * par("usr")[4], labels = LabelTxt, srt = 45, xpd = TRUE, adj = c(1,1))

dev.off()

GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/Citrine.pdf"))

