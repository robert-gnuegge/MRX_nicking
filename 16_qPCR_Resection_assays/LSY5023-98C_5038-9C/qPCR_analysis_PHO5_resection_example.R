# info --------------------------------------------------------------------
# purpose: analyze and plot qPCR-based resection assay data (see Gnugge et al., 2018 [pmid 29458754] for details)
# author: Robert Gnuegge (robert.gnuegge@gmail.com)
# date: 02/08/22
# version: 1.0


# preamble ----------------------------------------------------------------

# set working directory to this file's location
wd.path <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(wd.path)

# directories for plots and distribution data
plot_dir <- "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Figures/15_qPCR_Resection_assays/LSY5023-98C_5038-9C/Replicate_example"
moment_dir <- "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Data/qPCR_Resection_assays/LSY5023-98C_5038-9C/Replicate_example/Moments"
raw_data_dir <- "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Data/qPCR_Resection_assays/LSY5023-98C_5038-9C/Replicate_example/qPCR_data"

# read in helper files and functions
source(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Analysis_scripts/Misc/JFly_colors.R")
source(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Analysis_scripts/Misc/Misc_helper_functions.R")
source(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Analysis_scripts/Misc/CalcCutFraction2.R")
source(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Analysis_scripts/Misc/CalcResectedFraction2.R")

# read primer efficiencies
efficiencies <- read.table(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Resources/qPCR_primer_efficiencies/primer_efficiencies.txt", header = TRUE)

# read, check, and summarize Cq values ========================================

# define plate layout ---------------------------------------------------------
qPCR <- data.frame(well = c(paste0(rep(LETTERS[2:11], rep(24, 10)), formatC(x = 1:24, width = 2, format = "d", flag = "0")),
                            paste0(rep(LETTERS[12], 15), formatC(x = 1:15, width = 2, format = "d", flag = "0"))),
                   sample = c(rep(c("5023-98C", "5038-9C"), rep(24 * 5, 2)), rep("NTC", 15)),
                   time = c(rep(rep(c(0, 1, 2, 4, 6), rep(24, 5)), 2), rep(0, 15)),
                   digest = c(rep(c(rep("none", 15), rep("MseI", 9)), 10), rep("none", 15)),
                   primers = c(rep(rep(c("oRG54_oRG55", "oRG790_oRG791", "oRG782_oRG783", "oRG802_oRG803", "oRG806_oRG807", "oRG782_oRG783", "oRG802_oRG803", "oRG806_oRG807"), rep(3, 8)), 10),
                               rep(c("oRG54_oRG55", "oRG790_oRG791", "oRG782_oRG783", "oRG802_oRG803", "oRG806_oRG807"), rep(3, 5))),
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
# 33.53631
qPCR[qPCR$sample == "NTC", ]

# max non-NTC Cq
max(qPCR$Cq[qPCR$sample != "NTC"], na.rm = TRUE)
# 27.33716

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
# 0.04780513

# which are the samples with large cv?
Cq.modes[Cq.modes$cv > 0.025, ]
# only NTCs

# what is the ADH1 fold spread?
tmp <- Cq.modes$mean[Cq.modes$primers == "oRG54_oRG55" & Cq.modes$sample != "NTC"]
efficiencies$efficiency[efficiencies$primers == "oRG54_oRG55"]^(diff(range(tmp)))
# 1.886178

# what is the ADH1 fold spread within each strain
aggregate(mean ~ sample, data = Cq.modes[Cq.modes$primers == "oRG54_oRG55" & Cq.modes$sample != "NTC", ],
          FUN = function(x) efficiencies$efficiency[efficiencies$primers == "oRG54_oRG55"]^(diff(range(x))))
#     sample     mean
# 1 5023-98C 1.241087
# 2  5038-9C 1.886178

# calc cut fraction over time =================================================
cut.primers <- "oRG790_oRG791"
ref.primers <- "oRG54_oRG55"

# initialize data.frame to collect results
Cut.fraction <- data.frame()

# iterate through strains
strains <- unique(Cq.modes$sample[Cq.modes$sample != "NTC"])
for (strain in strains){
  
  # calc cut fraction
  tmp <- CalcCutFraction2(cut = Cq.modes[Cq.modes$sample == strain & Cq.modes$primers == cut.primers, c("time", "mean", "sd")], 
                          ref = Cq.modes[Cq.modes$sample == strain & Cq.modes$digest == "none" & Cq.modes$primers == ref.primers, c("time", "mean", "sd")], 
                          cut.primer.eff = efficiencies$efficiency[efficiencies$primers == cut.primers],
                          ref.primer.eff = efficiencies$efficiency[efficiencies$primers == ref.primers])
  
  # save results
  Cut.fraction <- rbind(Cut.fraction, cbind(strain = strain, primers = "oRG790_oRG791", tmp, stringsAsFactors = FALSE))  
}

# write data to file
write.table(x = Cut.fraction, file = paste0(moment_dir, "/Cutting.txt"), row.names = FALSE, col.names = TRUE)

# plot cutting ================================================================
strains <- c("5038-9C", "5023-98C")
MyColors <- JFly_colors[1:length(strains)]

# print to PDF
pdf(file = "tmp.pdf", width=3.5, height=3)
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(1, 0, 4, 2), las=1)

# start empty plot
y.range <- range(c(Cut.fraction$mean + Cut.fraction$sd, Cut.fraction$mean - Cut.fraction$sd, 0, 1))
# y.range <- c(-0.5, 1.5)
plot(x = NA, y = NA, ylim = y.range, xlim = c(0, 6), xlab = "Time [h]", ylab = "Cut Fraction")

# iterate through all samples
for (n in 1:length(strains)){
  
  strain <- strains[n]
  
  # time values with jitter
  t <- jitter(x = Cut.fraction$time[Cut.fraction$strain == strain], factor = 0)
  
  # line plot of means
  points(x = t, 
         y = Cut.fraction$mean[Cut.fraction$strain == strain],
         pch = 20, col = MyColors[n], type = "o"
  )
  
  # add error bars
  arrows(x0 = t,
         y0 = Cut.fraction$mean[Cut.fraction$strain == strain] - Cut.fraction$sd[Cut.fraction$strain == strain],
         x1 = t,
         y1 = Cut.fraction$mean[Cut.fraction$strain == strain] + Cut.fraction$sd[Cut.fraction$strain == strain],
         length = 0.03, # length of arrow head
         angle = 90, # angle of arrow head
         code = 3, # to draw arrow head on both ends
         col = MyColors[n])
}

legend(x = "bottomright", legend = c(expression(italic("pho4"*Delta)), expression(italic("pho4-SA1234PA6"))), 
       col = MyColors, pch = 20, inset = 0.05)

dev.off()

GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/Cut_fraction.pdf"))

# plot legend =================================================================
pdf(file = "tmp.pdf", width=1.75, height=0.65)
par(cex = 1, mar = rep(0, 4))
plot(1, type="n", axes=FALSE, xlab="", ylab="")
legend(1, 1, xjust=0.5, yjust=0.5, legend = c(expression(italic("pho4"*Delta)), expression(italic("pho4-SA1234PA6"))), 
       col = MyColors, pch = 20)
dev.off()
GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/Legend.pdf"))

# calc resected fraction over time ============================================
strains <- unique(Cq.modes$sample[Cq.modes$sample != "NTC"])

# initialize data.frame to collect results
Resected.fraction <- data.frame()

# iterate through strains
for (strain in strains){
  
  # iterate through resection evaluation amplicons
  for (resect.primers in c("oRG782_oRG783", "oRG802_oRG803", "oRG806_oRG807")){
    
    # calc resected fraction
    tmp <- CalcResectedFraction2(resect.mock = Cq.modes[Cq.modes$sample == strain & Cq.modes$primers == resect.primers & Cq.modes$digest == "none", c("time", "mean", "sd")], 
                                 resect.digest = Cq.modes[Cq.modes$sample == strain & Cq.modes$primers == resect.primers & Cq.modes$digest != "none", c("time", "mean", "sd")],
                                 resect.primer.eff = efficiencies$efficiency[efficiencies$primers == resect.primers],
                                 ref.mock = Cq.modes[Cq.modes$sample == strain & Cq.modes$primers == "oRG54_oRG55" & Cq.modes$digest == "none", c("time", "mean", "sd")], 
                                 ref.digest = Cq.modes[Cq.modes$sample == strain & Cq.modes$primers == "oRG54_oRG55" & Cq.modes$digest == "none", c("time", "mean", "sd")],
                                 ref.primer.eff = efficiencies$efficiency[efficiencies$primers == "oRG54_oRG55"],
                                 cut.fraction = Cut.fraction[Cut.fraction$strain == strain, c("time", "mean", "sd")]
    )
    
    # save results
    Resected.fraction <- rbind(Resected.fraction, cbind(strain = strain, primers = resect.primers, tmp, stringsAsFactors = FALSE))
  }

}

# write data to file
write.table(x = Resected.fraction, file = paste0(moment_dir, "/Resection.txt"), row.names = FALSE, col.names = TRUE)


# plot resected fraction ======================================================
strains <- c("5038-9C", "5023-98C")
MyColors <- JFlyColors[1:length(strains)]

# -66 bp ---------------------------------------------------------------------
tmp <- Resected.fraction[Resected.fraction$primers == "oRG782_oRG783", ]

# print to PDF
pdf(file = "tmp.pdf", width=3.5, height=3)
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(1, 0, 4, 2), las=1)

# start empty plot
y.range <- range(c(tmp$mean + tmp$sd, tmp$mean - tmp$sd, 0, 1), na.rm = TRUE)
# y.range <- c(-0.5, 2)
plot(x = NA, y = NA, ylim = y.range, xlim = c(0, 6), xlab = "Time [h]", ylab = "ssDNA Fraction")

# iterate through all samples
for (n in 1:length(strains)){
  
  strain <- strains[n]
  
  # time values with jitter
  t <- jitter(x = Cut.fraction$time[Cut.fraction$strain == strain], factor = 0)
  
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

GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/Resection_-66_bp.pdf"))


# -223 bp ---------------------------------------------------------------------
tmp <- Resected.fraction[Resected.fraction$primers == "oRG802_oRG803", ]

# print to PDF
pdf(file = "tmp.pdf", width=3.5, height=3)
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(1, 0, 4, 2), las=1)

# start empty plot
y.range <- range(c(tmp$mean + tmp$sd, tmp$mean - tmp$sd, 0, 1), na.rm = TRUE)
# y.range <- c(-0.5, 2)
plot(x = NA, y = NA, ylim = y.range, xlim = c(0, 6), xlab = "Time [h]", ylab = "ssDNA Fraction")

# iterate through all samples
for (n in 1:length(strains)){
  
  strain <- strains[n]
  
  # time values with jitter
  t <- jitter(x = Cut.fraction$time[Cut.fraction$strain == strain], factor = 0)
  
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

GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/Resection_-223_bp.pdf"))

# -538 bp ---------------------------------------------------------------------
tmp <- Resected.fraction[Resected.fraction$primers == "oRG806_oRG807", ]

# print to PDF
pdf(file = "tmp.pdf", width=3.5, height=3)
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(1, 0, 4, 2), las=1)

# start empty plot
y.range <- range(c(tmp$mean + tmp$sd, tmp$mean - tmp$sd, 0, 1), na.rm = TRUE)
# y.range <- c(-0.5, 2)
plot(x = NA, y = NA, ylim = y.range, xlim = c(0, 6), xlab = "Time [h]", ylab = "ssDNA Fraction")

# iterate through all samples
for (n in 1:length(strains)){
  
  strain <- strains[n]
  
  # time values with jitter
  t <- jitter(x = Cut.fraction$time[Cut.fraction$strain == strain], factor = 0)
  
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

GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/Resection_-538_bp.pdf"))