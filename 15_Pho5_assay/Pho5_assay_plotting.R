# info --------------------------------------------------------------------
# purpose: plot Pho5 assay averages
# author: Robert Gnuegge (robert.gnuegge@gmail.com)
# date: 02/25/22
# version: 1.0


# preamble ----------------------------------------------------------------

# set wd to this file's location
wd_path <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(wd_path)

# load helper files and functions
source(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Analysis_scripts/Misc/Misc_helper_functions.R")


# process data ------------------------------------------------------------

# read data
ODs <- read.table(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Data/Pho5_assay/OD_405nm.txt", header = TRUE)

replicate_means <- aggregate(OD.405nm ~ Strain + replicate, data = ODs, FUN = mean)

# calculate PHO5 units
coeff <- 4519  # OD to conc (M) conversion
time <- 16  # min
OD <- 1.5  # OD_600nm (cell titer)
V <- 0.25e-3  # reaction volume
conversion <- 1e9  # to nM
replicate_means <- cbind(replicate_means, U = replicate_means$OD.405nm / coeff / time / OD * conversion * V)

# calc mean and sd
modes <- aggregate(U ~ Strain, data = replicate_means, FUN = function(x) c(mean = mean(x), sd = sd(x)))
# mean and sd are in a single column (as a matrix); convert into separate columns
modes <- data.frame(modes[, 1:(ncol(modes)-1)], as.data.frame(modes$U))
colnames(modes)[1] <- "Strain"

# plotting ----------------------------------------------------------------

pdf(file = "tmp.pdf", width=2, height=3.5)
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(-0.7, 0.6, 3.7, 2.1), xpd = TRUE, las=1, tcl = -0.3, mgp = c(2.5, 0.6, 0))

idx <- c(1, 3, 2)
bp <- barplot(height = modes$mean[idx], ylab = "Pho5 Activity [U]", ylim = c(0, 6))

# add data points
raw_data <- split(x = replicate_means, f = ~ Strain)
raw_data <- raw_data[idx]
stripchart(at = bp[, 1], x = lapply(X = raw_data, FUN = function(x){x$U}), 
           vertical = TRUE, pch = 21, method = "jitter", jitter = 0.33, add = TRUE, col = gray(level = 0.4), cex = 0.5)

# add error bars
arrows(x0 = bp[, 1],
       y0 = modes$mean[idx] - modes$sd[idx],
       x1 = bp[, 1],
       y1 = modes$mean[idx] + modes$sd[idx],
       length = 0.03, # length of arrow head    
       angle = 90, # angle of arrow head
       code = 3 # to draw arrow head on both ends
)

LabelTxt <- c("Wild type",
              expression(italic("pho4")*Delta),
              expression(italic("pho4-SA1234PA6")))

text(x = bp[, 1], y = par("usr")[3] - 0.3, labels = LabelTxt, srt = 45, xpd = TRUE, adj = c(1, 0.5))

dev.off()

# embed fonts
GS_embed_fonts(input = "tmp.pdf", output = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Figures/18_Pho5_assay/PHO5_activity.pdf")