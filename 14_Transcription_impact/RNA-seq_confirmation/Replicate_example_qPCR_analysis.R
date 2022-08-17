# info --------------------------------------------------------------------
# purpose: calculate and plot relative transcript levels
# author: Robert Gnuegge (robert.gnuegge@gmail.com)
# date: 08/16/22
# version: 1.0


# preamble ----------------------------------------------------------------

# set wd to this file's location
wd_path <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(wd_path)

# load libraries
library(propagate)
library(GenomicRanges)
library(rtracklayer)

# set up directory paths
plot_dir <- "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Figures/11_Transcription_impact/RNA-seq_confirmation/Replicate_example"
data_dir <- "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Data/RNA-seq_confirmation/Replicates"

# read in helper files and functions
source(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Analysis_scripts/Misc/Genomic_helper_functions.R")
source(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Analysis_scripts/Misc/Misc_helper_functions.R")
source(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Analysis_scripts/Misc/JFly_colors.R")

# read primer efficiencies
efficiencies <- read.table(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Resources/qPCR_primer_efficiencies/primer_efficiencies.txt", header = TRUE)


# read, check, and summarize Cq values ====================================

# define plate layout -----------------------------------------------------
qPCR <- data.frame(well = c(paste0(rep(LETTERS[c(2, 3, 6)], rep(18, 3)), formatC(x = 2:19, width = 2, format = "d", flag = "0")),
                            paste0(rep(LETTERS[c(7, 8, 11)], rep(15, 3)), formatC(x = 2:16, width = 2, format = "d", flag = "0"))),
                   sample = c(rep("+RT", 18), rep("-RT", 18), rep("NTC", 18), rep("+RT", 15), rep("-RT", 15), rep("NTC", 15)),
                   primers = c(rep(rep(c("oRG876_oRG877", "oRG878_oRG879", "oRG882_oRG883", "oRG888_oRG889", "oRG892_oRG893", "oRG896_oRG897"), rep(3, 6)), 3),
                               rep(rep(c("oRG898_oRG899", "oRG903_oRG904", "oRG907_oRG908", "oRG909_oRG910", "oRG54_oRG55"), rep(3, 5)), 3)),
                   Cq = 0)

nrow(qPCR) == 18 * 3 + 15 * 3
head(qPCR)
tail(qPCR)

# read Cq Results file --------------------------------------------------------
file.path <- "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Data/RNA-seq_confirmation/qPCR_data/admin_2022-08-12 11-55-31_CT017724 -  Quantification Cq Results_0.txt"
tmp <- read.table(file = file.path, sep = "\t", header = TRUE)

# save Cq values in qPCR data frame
for (well in qPCR$well){
  qPCR$Cq[qPCR$well == well] <- tmp$Cq[tmp$Well == well & tmp$Fluor == "SYBR"]
}

head(qPCR)

# check NTCs ------------------------------------------------------------------

# min NTC Cq
min(qPCR$Cq[qPCR$sample == "NTC"], na.rm = TRUE)
qPCR[qPCR$sample == "NTC", ]
# NaN

# max non-NTC Cq
max(qPCR$Cq[qPCR$sample != "NTC"], na.rm = TRUE)
# 37.99326

# Confirm that the DNase treatment was sufficient by checking the Cq difference between +/-RT samples!
summary(qPCR$Cq[qPCR$sample == "-RT"] - qPCR$Cq[qPCR$sample == "+RT"])
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#   10.99   13.94   15.19   15.97   19.01   20.35      15

# calc Cq mean, sd, and cv ----------------------------------------------------
Cq.modes <- aggregate(Cq ~ sample + primers, data = qPCR,
                      FUN = function(x) c(mean = mean(x, na.rm = TRUE), sd = sd(x, na.rm = TRUE), cv = sd(x, na.rm = TRUE)/mean(x, na.rm = TRUE)),
                      na.action = na.pass)

# sort
Cq.modes <- Cq.modes[order(Cq.modes$sample, Cq.modes$primers, decreasing = c(FALSE, TRUE), method = "radix"),]

# mean, sd, and cv are in a single column (as a matrix)
# convert into separate columns
Cq.modes <- data.frame(Cq.modes[, 1:(ncol(Cq.modes)-1)], as.data.frame(Cq.modes$Cq))

# what is the largest cv?
max(Cq.modes$cv, na.rm = TRUE)
# 0.04123197

# which are the samples with large cv?
Cq.modes[Cq.modes$cv > 0.025, ]
# only "-RT" samples


# calc relative transcript levels =========================================

# function to calculate relative transcript levels
Pfaffl <- function(Cq_target_triplicate, Cq_ref_triplicate, E_target = 2, E_ref = 2){
  
  data <- cbind(Cq_target = c(mean(Cq_target_triplicate, na.rm = TRUE), sd(Cq_target_triplicate, na.rm = TRUE)),
                Cq_ref = c(mean(Cq_ref_triplicate, na.rm = TRUE), sd(Cq_ref_triplicate, na.rm = TRUE)))
  
  expr <- as.expression(substitute(expr = E_target^-Cq_target / E_ref^-Cq_ref, env = list(E_target = E_target, E_ref = E_ref)))
  
  res <- propagate(expr = expr, data = data, do.sim = FALSE)
  
  return(res$prop[c("Mean.2", "sd.2", "2.5%", "97.5%")])
  
}

# set up data.frame
RNA <- Cq.modes[Cq.modes$sample != "NTC", c("sample", "primers")]
RNA$gene <- efficiencies$amplicon[match(x = RNA$primers, table = efficiencies$primers)]
RNA$mean <- 0
RNA$sd <- 0
RNA$conf97.5 <- 0
RNA$conf2.5 <- 0


# fill in transcript levels relative to ADH1 (for same concentration)
for(n in 1:nrow(RNA)){
  
  sample <- RNA$sample[n]
  conc <- RNA$conc[n]
  primers <- RNA$primers[n]
  
  tmp <- Pfaffl(Cq_target_triplicate = qPCR$Cq[qPCR$sample == sample & qPCR$primers == primers],
                Cq_ref_triplicate = qPCR$Cq[qPCR$sample == "+RT" & qPCR$primers == "oRG54_oRG55"],
                E_target = efficiencies$efficiency[efficiencies$primers == primers],
                E_ref = efficiencies$efficiency[efficiencies$primers == "oRG54_oRG55"])
  
  RNA$mean[n] <- tmp["Mean.2"]
  RNA$sd[n] <- tmp["sd.2"]
  RNA$conf97.5[n] <- tmp["97.5%"]
  RNA$conf2.5[n] <- tmp["2.5%"]
  
}

# save data
write.table(x = RNA, file = paste0(data_dir, "/RT-qPCR_levels_rep2.txt"), row.names = FALSE)


# plot relative transcript levels -----------------------------------------

# assemble data for grouped bar plot
tmp <- RNA[RNA$sample == "+RT", ]
gene_order <- RNA$gene[order(tmp$mean, decreasing = TRUE)]

means <- tmp$mean[match(x = gene_order, table = tmp$gene)]
sd <- tmp$sd[match(x = gene_order, table = tmp$gene)]

tmp <- RNA[RNA$sample == "-RT", ]
means <- rbind(means, tmp$mean[match(x = gene_order, table = tmp$gene)])
sd <- cbind(sd, tmp$sd[match(x = gene_order, table = tmp$gene)])

sd <- t(sd)  # transpose to make it parallel to means

# print to PDF
pdf(file = "tmp.pdf", width=4.5, height=3)
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(2.4, 0.4, 3.4, 2.1), tcl = -0.25, mgp = c(2.75, 0.5, 0), las = 1)

range(c(means + sd, means - sd), na.rm = TRUE)
bp <- barplot(height = means, beside = TRUE, log = "y", ylim = c(0.67e-7, 1.5), yaxt = "n", ylab = "Relative mRNA level")
axis(side = 2, at = 10^(0:-7), 
     labels = c(expression(10^0), expression(10^-1), expression(10^-2), expression(10^-3),
                expression(10^-4), expression(10^-5), expression(10^-6), expression(10^-7)))

arrows(x0 = c(bp),  # c(matrix) transforms into vector
       y0 = c(means) + c(sd),
       x1 = c(bp),
       y1 = c(means) - c(sd),
       length = 0.03, # length of arrow head
       angle = 90, # angle of arrow head
       code = 3, # to draw arrow head on both ends
       col = "black")

text(x = apply(X = bp, MARGIN = 2, FUN = mean), y = 0.4e-7, srt = 45, 
     labels = gene_order, xpd = TRUE, adj = c(1,1), font = 3)

legend(x = "topright", inset = c(0.03, -0.05), fill = gray.colors(2), # col = gray.colors(n) is used by barplot by default
       ncol = 2, legend = c("+RT", "-RT"), xpd = TRUE)

dev.off()

GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/RT-qPCR_levels_rep2.pdf"))