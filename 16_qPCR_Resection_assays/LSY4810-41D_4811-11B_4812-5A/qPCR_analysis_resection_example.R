# info --------------------------------------------------------------------
# purpose: analyze and plot qPCR-based resection assay data (see Gnugge et al., 2018 [pmid 29458754] for details)
# author: Robert Gnuegge (robert.gnuegge@gmail.com)
# date: 02/09/22
# version: 1.0

# preamble ----------------------------------------------------------------


# set working directory to this file's location
wd.path <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(wd.path)

# directories for plots and distribution data
plot_dir <- "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Figures/15_qPCR_Resection_assays/LSY4810-41D_4811-11B_4812-5A/Replicate_example"
moment_dir <- "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Data/qPCR_Resection_assays/LSY4810-41D_4811-11B_4812-5A/Replicate_example/Moments"
raw_data_dir <- "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Data/qPCR_Resection_assays/LSY4810-41D_4811-11B_4812-5A/Replicate_example/qPCR_data"

# read in helper files and functions
source(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Analysis_scripts/Misc/JFly_colors.R")
source(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Analysis_scripts/Misc/Misc_helper_functions.R")
source(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Analysis_scripts/Misc/CalcCutFraction2.R")
source(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Analysis_scripts/Misc/CalcResectedFraction2.R")

# read primer efficiencies
efficiencies <- read.table(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Resources/qPCR_primer_efficiencies/primer_efficiencies.txt", header = TRUE)

# read, check, and summarize Cq values ========================================

# define plate layout ---------------------------------------------------------
qPCR <- data.frame(well = c(paste0(rep(LETTERS[1:15], rep(21, 15)), formatC(x = 2:22, width = 2, format = "d", flag = "0")),
                            paste0("P", formatC(x = 2:16, width = 2, format = "d", flag = "0")),
                            paste0(LETTERS[2:7], formatC(x = 23, width = 2, format = "d", flag = "0"))),
                   sample = c(rep(c("4810-41D", "4811-11B", "4812-5A"), rep(21 * 5, 3)), rep("NTC", 21)),
                   time = c(rep(rep(c(0, 1, 2, 4, 6), rep(21, 5)), 3), rep(0, 21)),
                   digest = c(rep(c(rep("none", 15), rep("AluI", 3), rep("StyI", 3)), 5), rep(c(rep("none", 15), rep("AluI", 6)), 10), rep("none", 21)),
                   primers = c(rep(rep(c("oRG54_oRG55", "oRG641_oRG642", "oRG645_oRG646", "oRG711_oRG712", "oRG721_oRG722", "oRG711_oRG712", "oRG721_oRG722"), rep(3, 7)), 5),
                               rep(rep(c("oRG54_oRG55", "oRG641_oRG642", "oRG645_oRG646", "oRG711_oRG712", "oRG723_oRG724", "oRG711_oRG712", "oRG723_oRG724"), rep(3, 7)), 5),
                               rep(rep(c("oRG54_oRG55", "oRG641_oRG642", "oRG645_oRG646", "oRG709_oRG710", "oRG723_oRG724", "oRG709_oRG710", "oRG723_oRG724"), rep(3, 7)), 5),
                               rep(c("oRG54_oRG55", "oRG641_oRG642", "oRG645_oRG646", "oRG709_oRG710", "oRG723_oRG724", "oRG711_oRG712", "oRG721_oRG722"), rep(3, 7))),
                   Cq = 0)

head(qPCR)
tail(qPCR)

# read Cq Results file --------------------------------------------------------
file.path <- grep(pattern = "Cq", x = list.files(path = raw_data_dir, full.names = TRUE, recursive = TRUE), value = TRUE)

# save Cq values in qPCR data frame
tmp <- read.table(file = file.path, sep = "\t", header = TRUE)

for (well in qPCR$well){
  qPCR$Cq[qPCR$well == well] <- tmp$Cq[tmp$Well == well & tmp$Fluor == "SYBR"]
}

head(qPCR)

# check NTCs ------------------------------------------------------------------

# min NTC Cq
min(qPCR$Cq[qPCR$sample == "NTC"], na.rm = TRUE)
# 37.67751
qPCR[qPCR$sample == "NTC", ]

# max non-NTC Cq
max(qPCR$Cq[qPCR$sample != "NTC"], na.rm = TRUE)
# 28.5616

# check non-NTC samples with high Cq values
qPCR[qPCR$Cq > 25 & qPCR$sample != "NTC",]

# calc Cq mean, sd, and cv ----------------------------------------------------
Cq.modes <- aggregate(Cq ~ sample + digest + time + primers, data = qPCR,
                      FUN = function(x) c(mean = mean(x, na.rm = TRUE),
                                          sd = sd(x, na.rm = TRUE),
                                          cv = sd(x, na.rm = TRUE)/mean(x, na.rm = TRUE)
                                          )
                      )

# sort
Cq.modes <- Cq.modes[order(Cq.modes$sample, Cq.modes$primers, Cq.modes$digest, Cq.modes$time),]

# mean, sd, and cv are in a single column (as a matrix)
# convert into separate columns
Cq.modes <- data.frame(Cq.modes[, 1:(ncol(Cq.modes)-1)], as.data.frame(Cq.modes$Cq))

# what is the largest cv?
max(Cq.modes$cv, na.rm = TRUE)
# 0.02013872

# which are the samples with large cv?
Cq.modes[Cq.modes$cv > 0.025, ]

# what is the ADH1 fold spread?
tmp <- Cq.modes$mean[Cq.modes$primers == "oRG54_oRG55" & Cq.modes$sample != "NTC"]
efficiencies$efficiency[efficiencies$primers == "oRG54_oRG55"]^(diff(range(tmp)))
# 2.358656

# what is the ADH1 fold spread within each strain
aggregate(mean ~ sample, data = Cq.modes[Cq.modes$primers == "oRG54_oRG55" & Cq.modes$sample != "NTC", ],
          FUN = function(x) efficiencies$efficiency[efficiencies$primers == "oRG54_oRG55"]^(diff(range(x))))
#     sample     mean
# 1 4810-41D 1.730975
# 2 4811-11B 1.385996
# 3  4812-5A 2.170406

# calc cut fraction over time =================================================

# initialize data.frame to collect results
Cut.fraction <- data.frame()

ref.primers <- "oRG54_oRG55"

# iterate through strains
strains <- c("4810-41D", "4811-11B", "4812-5A")
for (strain in strains){
  
  for(cut.primers in c("oRG641_oRG642", "oRG645_oRG646")){
    # calc cut fraction
    tmp <- CalcCutFraction2(cut = Cq.modes[Cq.modes$sample == strain & Cq.modes$primers == cut.primers, c("time", "mean", "sd")], 
                            ref = Cq.modes[Cq.modes$sample == strain & Cq.modes$digest == "none" & Cq.modes$primers == ref.primers, c("time", "mean", "sd")], 
                            cut.primer.eff = efficiencies$efficiency[efficiencies$primers == cut.primers],
                            ref.primer.eff = efficiencies$efficiency[efficiencies$primers == ref.primers])
    
    # save results
    Cut.fraction <- rbind(Cut.fraction, cbind(strain = strain, primers = cut.primers, tmp, stringsAsFactors = FALSE))
  }
  
}

# write data to file
write.table(x = Cut.fraction, file = paste0(moment_dir, "/Cutting.txt"), row.names = FALSE, col.names = TRUE)

# plot cutting ================================================================

# mKate2 ----------------------------------------------------------------------
strains <- c("4812-5A", "4811-11B", "4810-41D")
MyColors <- JFly_colors[1:length(strains)]

tmp <- Cut.fraction[Cut.fraction$primers == "oRG641_oRG642", ]

# print to PDF
pdf(file = "tmp.pdf", width=3.5, height=3)
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(1, 0, 4, 2), las=1)

# start empty plot
y.range <- range(c(tmp$mean + tmp$sd, tmp$mean - tmp$sd, 0, 1))
# y.range <- c(-0.5, 1.5)
plot(x = NA, y = NA, ylim = y.range, xlim = range(tmp$time), xlab = "Time [h]", ylab = "Cut Fraction")

# iterate through all samples
for (n in 1:length(strains)){
  
  strain <- strains[n]
  
  # time values with jitter
  t <- jitter(x = tmp$time[tmp$strain == strain], factor = 0)
  
  # line plot of means
  points(x = t, 
         y = tmp$mean[tmp$strain == strain],
         pch = 20, col = MyColors[n], type = "o"
  )
  
  # add error bars
  arrows(x0 = t,
         y0 = tmp$mean[tmp$strain == strain] - tmp$sd[tmp$strain == strain],
         x1 = t,
         y1 = tmp$mean[tmp$strain == strain] + tmp$sd[tmp$strain == strain],
         length = 0.03, # length of arrow head
         angle = 90, # angle of arrow head
         code = 3, # to draw arrow head on both ends
         col = MyColors[n])
}

dev.off()

GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/mKate2_cut_fraction.pdf"))

# plot legend -----------------------------------------------------------------
LegendTxt <- c(expression(italic("mKate2-HOcs")),
               expression(italic("P")[italic("ACT1")]*italic("-mKate2-HOcs")),
               expression(italic("P")[italic("TDH3")]*italic("-mKate2-HOcs")))

MyColors <- c(JFly_colors[1:3])

pdf(file = "tmp.pdf", width=2.25, height=0.85)
par(cex = 1, mar = rep(0, 4))

plot(1, type="n", axes=FALSE, xlab="", ylab="")
legend(1, 1, pch = 20, col = MyColors, lty = c(rep("solid", 3), "dashed"), legend = LegendTxt, xjust=0.5, yjust=0.5, ncol = 1)

dev.off()

GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/mKate2_Legend.pdf"))

# Citrine ---------------------------------------------------------------------
strains <- c("4810-41D", "4812-5A", "4811-11B")
MyColors <- JFly_colors[1:3]

tmp <- Cut.fraction[Cut.fraction$primers == "oRG645_oRG646", ]

# print to PDF
pdf(file = "tmp.pdf", width=3.5, height=3)
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(1, 0, 4, 2), las=1)

# start empty plot
y.range <- range(c(tmp$mean + tmp$sd, tmp$mean - tmp$sd, 0, 1))
# y.range <- c(-0.5, 1.5)
plot(x = NA, y = NA, ylim = y.range, xlim = range(tmp$time), xlab = "Time [h]", ylab = "Cut Fraction")

# iterate through all samples
for (n in 1:length(strains)){
  
  strain <- strains[n]
  
  # time values with jitter
  t <- jitter(x = tmp$time[tmp$strain == strain], factor = 0)
  
  # line plot of means
  points(x = t, 
         y = tmp$mean[tmp$strain == strain],
         pch = 20, col = MyColors[n], type = "o"
  )
  
  # add error bars
  arrows(x0 = t,
         y0 = tmp$mean[tmp$strain == strain] - tmp$sd[tmp$strain == strain],
         x1 = t,
         y1 = tmp$mean[tmp$strain == strain] + tmp$sd[tmp$strain == strain],
         length = 0.03, # length of arrow head
         angle = 90, # angle of arrow head
         code = 3, # to draw arrow head on both ends
         col = MyColors[n])
}

dev.off()

GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/Citrine_cut_fraction.pdf"))

# plot legend -----------------------------------------------------------------
LegendTxt <- c(expression(italic("Citrine-HOcs")),
               expression(italic("P")[italic("ACT1")]*italic("-Citrine-hisG-HOcs")),
               expression(italic("P")[italic("TDH3")]*italic("-Citrine-hisG-HOcs")))

MyColors <- c(JFly_colors[1:3])

pdf(file = "tmp.pdf", width=2.6, height=0.85)
par(cex = 1, mar = rep(0, 4))

plot(1, type="n", axes=FALSE, xlab="", ylab="")
legend(1, 1, pch = 20, col = MyColors, lty = c(rep("solid", 3), "dashed"), legend = LegendTxt, xjust=0.5, yjust=0.5, ncol = 1)

dev.off()

GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/Citrine_Legend.pdf"))

# calc resected fraction over time ============================================

# initialize data.frame to collect results
Resected.fraction <- data.frame()

# mKate2 ----------------------------------------------------------------------
cut.primers <- "oRG641_oRG642"

strains <- c("4810-41D", "4811-11B")
resect.primers <- "oRG711_oRG712"

# iterate through strains
for (strain in strains){
  
  # calc resected fraction
  tmp <- CalcResectedFraction2(resect.mock = Cq.modes[Cq.modes$sample == strain & Cq.modes$primers == resect.primers & Cq.modes$digest == "none", c("time", "mean", "sd")], 
                               resect.digest = Cq.modes[Cq.modes$sample == strain & Cq.modes$primers == resect.primers & Cq.modes$digest != "none", c("time", "mean", "sd")],
                               resect.primer.eff = efficiencies$efficiency[efficiencies$primers == resect.primers],
                               ref.mock = Cq.modes[Cq.modes$sample == strain & Cq.modes$primers == "oRG54_oRG55" & Cq.modes$digest == "none", c("time", "mean", "sd")], 
                               ref.digest = Cq.modes[Cq.modes$sample == strain & Cq.modes$primers == "oRG54_oRG55" & Cq.modes$digest == "none", c("time", "mean", "sd")],
                               ref.primer.eff = efficiencies$efficiency[efficiencies$primers == "oRG54_oRG55"],
                               cut.fraction = Cut.fraction[Cut.fraction$strain == strain & Cut.fraction$primers == cut.primers, c("time", "mean", "sd")]
  )
    
  # save results
  Resected.fraction <- rbind(Resected.fraction, cbind(strain = strain, primers = resect.primers, tmp, stringsAsFactors = FALSE))
}

strains <- c("4812-5A")
resect.primers <- "oRG709_oRG710"

# iterate through strains
for (strain in strains){
  
  # calc resected fraction
  tmp <- CalcResectedFraction2(resect.mock = Cq.modes[Cq.modes$sample == strain & Cq.modes$primers == resect.primers & Cq.modes$digest == "none", c("time", "mean", "sd")], 
                               resect.digest = Cq.modes[Cq.modes$sample == strain & Cq.modes$primers == resect.primers & Cq.modes$digest != "none", c("time", "mean", "sd")],
                               resect.primer.eff = efficiencies$efficiency[efficiencies$primers == resect.primers],
                               ref.mock = Cq.modes[Cq.modes$sample == strain & Cq.modes$primers == "oRG54_oRG55" & Cq.modes$digest == "none", c("time", "mean", "sd")], 
                               ref.digest = Cq.modes[Cq.modes$sample == strain & Cq.modes$primers == "oRG54_oRG55" & Cq.modes$digest == "none", c("time", "mean", "sd")],
                               ref.primer.eff = efficiencies$efficiency[efficiencies$primers == "oRG54_oRG55"],
                               cut.fraction = Cut.fraction[Cut.fraction$strain == strain & Cut.fraction$primers == cut.primers, c("time", "mean", "sd")]
  )
  
  # save results
  Resected.fraction <- rbind(Resected.fraction, cbind(strain = strain, primers = resect.primers, tmp, stringsAsFactors = FALSE))
}

# Citrine ----------------------------------------------------------------------
cut.primers <- "oRG645_oRG646"

strains <- c("4810-41D")
resect.primers <- "oRG721_oRG722"

# iterate through strains
for (strain in strains){
  
  # calc resected fraction
  tmp <- CalcResectedFraction2(resect.mock = Cq.modes[Cq.modes$sample == strain & Cq.modes$primers == resect.primers & Cq.modes$digest == "none", c("time", "mean", "sd")], 
                               resect.digest = Cq.modes[Cq.modes$sample == strain & Cq.modes$primers == resect.primers & Cq.modes$digest != "none", c("time", "mean", "sd")],
                               resect.primer.eff = efficiencies$efficiency[efficiencies$primers == resect.primers],
                               ref.mock = Cq.modes[Cq.modes$sample == strain & Cq.modes$primers == "oRG54_oRG55" & Cq.modes$digest == "none", c("time", "mean", "sd")], 
                               ref.digest = Cq.modes[Cq.modes$sample == strain & Cq.modes$primers == "oRG54_oRG55" & Cq.modes$digest == "none", c("time", "mean", "sd")],
                               ref.primer.eff = efficiencies$efficiency[efficiencies$primers == "oRG54_oRG55"],
                               cut.fraction = Cut.fraction[Cut.fraction$strain == strain & Cut.fraction$primers == cut.primers, c("time", "mean", "sd")]
  )
  
  # save results
  Resected.fraction <- rbind(Resected.fraction, cbind(strain = strain, primers = resect.primers, tmp, stringsAsFactors = FALSE))
}

strains <- c("4811-11B", "4812-5A")
resect.primers <- "oRG723_oRG724"

# iterate through strains
for (strain in strains){
  
  # calc resected fraction
  tmp <- CalcResectedFraction2(resect.mock = Cq.modes[Cq.modes$sample == strain & Cq.modes$primers == resect.primers & Cq.modes$digest == "none", c("time", "mean", "sd")], 
                               resect.digest = Cq.modes[Cq.modes$sample == strain & Cq.modes$primers == resect.primers & Cq.modes$digest != "none", c("time", "mean", "sd")],
                               resect.primer.eff = efficiencies$efficiency[efficiencies$primers == resect.primers],
                               ref.mock = Cq.modes[Cq.modes$sample == strain & Cq.modes$primers == "oRG54_oRG55" & Cq.modes$digest == "none", c("time", "mean", "sd")], 
                               ref.digest = Cq.modes[Cq.modes$sample == strain & Cq.modes$primers == "oRG54_oRG55" & Cq.modes$digest == "none", c("time", "mean", "sd")],
                               ref.primer.eff = efficiencies$efficiency[efficiencies$primers == "oRG54_oRG55"],
                               cut.fraction = Cut.fraction[Cut.fraction$strain == strain & Cut.fraction$primers == cut.primers, c("time", "mean", "sd")]
  )
  
  # save results
  Resected.fraction <- rbind(Resected.fraction, cbind(strain = strain, primers = resect.primers, tmp, stringsAsFactors = FALSE))
}

# write data to file
write.table(x = Resected.fraction, file = paste0(moment_dir, "/Resection.txt"), row.names = FALSE, col.names = TRUE)


# plot resected fraction ======================================================

# mKate2 -------------------------------------------------------------------------
strains <- c("4812-5A", "4811-11B", "4810-41D")
MyColors <- JFly_colors[1:length(strains)]

tmp <- Resected.fraction[Resected.fraction$strain %in% strains & Resected.fraction$primers %in% c("oRG711_oRG712", "oRG709_oRG710"), ]

# print to PDF
pdf(file = "tmp.pdf", width=3.5, height=3)
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(1, 0, 4, 2), las=1)

# start empty plot
y.range <- range(c(tmp$mean + tmp$sd, tmp$mean - tmp$sd, 0, 1), na.rm = TRUE)
# y.range <- c(-0.5, 2)
plot(x = NA, y = NA, ylim = y.range, xlim = range(tmp$time), xlab = "Time [h]", ylab = "ssDNA Fraction")

# iterate through all samples
for (n in 1:length(strains)){
  
  strain <- strains[n]
  
  # time values with jitter
  t <- jitter(x = tmp$time[tmp$strain == strain], factor = 0)
  
  # line plot of means
  points(x = t, 
         y = tmp$mean[tmp$strain == strain],
         pch = 20, col = MyColors[n], type = "o")
  
  # add error bars
  arrows(x0 = t,
         y0 = tmp$mean[tmp$strain == strain] - tmp$sd[tmp$strain == strain],
         x1 = t,
         y1 = tmp$mean[tmp$strain == strain] + tmp$sd[tmp$strain == strain],
         length = 0.03, # length of arrow head
         angle = 90, # angle of arrow head
         code = 3, # to draw arrow head on both ends
         col = MyColors[n])
}

dev.off()

GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/mKate2_resection_-430_bp.pdf"))


# Citrine -------------------------------------------------------------------------
strains <- c("4810-41D", "4812-5A", "4811-11B")
MyColors <- JFly_colors[1:length(strains)]

tmp <- Resected.fraction[Resected.fraction$strain %in% strains & Resected.fraction$primers %in% c("oRG721_oRG722", "oRG723_oRG724"), ]

# print to PDF
pdf(file = "tmp.pdf", width=3.5, height=3)
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(1, 0, 4, 2), las=1)

# start empty plot
y.range <- range(c(tmp$mean + tmp$sd, tmp$mean - tmp$sd, 0, 1), na.rm = TRUE)
# y.range <- c(-0.5, 2)
plot(x = NA, y = NA, ylim = y.range, xlim = range(tmp$time), xlab = "Time [h]", ylab = "ssDNA Fraction")

# iterate through all samples
for (n in 1:length(strains)){
  
  strain <- strains[n]
  
  # time values with jitter
  t <- jitter(x = tmp$time[tmp$strain == strain], factor = 0)
  
  # line plot of means
  points(x = t, 
         y = tmp$mean[tmp$strain == strain],
         pch = 20, col = MyColors[n], type = "o")
  
  # add error bars
  arrows(x0 = t,
         y0 = tmp$mean[tmp$strain == strain] - tmp$sd[tmp$strain == strain],
         x1 = t,
         y1 = tmp$mean[tmp$strain == strain] + tmp$sd[tmp$strain == strain],
         length = 0.03, # length of arrow head
         angle = 90, # angle of arrow head
         code = 3, # to draw arrow head on both ends
         col = MyColors[n])
}

dev.off()

GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/Citrine_resection_-565_bp.pdf"))