# info --------------------------------------------------------------------
# purpose: analyze and plot qPCR-based resection assay data (see Gnugge et al., 2018 [pmid 29458754] for details)
# author: Robert Gnuegge (robert.gnuegge@gmail.com)
# date: 02/22/22
# version: 1.0


# preamble ----------------------------------------------------------------

# set working directory to this file's location
wd.path <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(wd.path)

# directories for plots and distribution data
plot_dir <- "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Figures/15_qPCR_Resection_assays/LSY4377-12B_4377-15A/Replicate_example"
moment_dir <- "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Data/qPCR_Resection_assays/LSY4377-12B_4377-15A/Replicate_example/Moments"
raw_data_dir <- "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Data/qPCR_Resection_assays/LSY4377-12B_4377-15A/Replicate_example/qPCR_data"

# read in helper files and functions
source(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Analysis_scripts/Misc/JFly_colors.R")
source(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Analysis_scripts/Misc/Misc_helper_functions.R")
source(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Analysis_scripts/Misc/CalcCutFraction2.R")
source(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Analysis_scripts/Misc/CalcResectedFraction2.R")

# read primer efficiencies
efficiencies <- read.table(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Resources/qPCR_primer_efficiencies/primer_efficiencies.txt", header = TRUE)

# read, check, and summarize Cq values ========================================

# set up data frame
qPCR <- data.frame(well = c(paste0(rep(LETTERS[1:16], rep(21, 16)), formatC(x = 2:22, width = 2, format = "d", flag = "0")), # samples
                            paste0(LETTERS[3:14], 23)), # NTCs
                   sample = c(rep(c("4376-3C", "4376-4A", "4377-12B", "4377-15A"), rep(21*4, 4)), # samples
                              rep("NTC", 12)), # NTCs
                   digest = c(rep(c(rep("none", 12), rep("RsaI", 9)), 16), # samples
                              rep("none", 12)), # NTCs
                   time = c(rep(rep(c(0, 1, 2, 4), rep(21, 4)), 4), # samples
                            rep(0, 12)), # NTCs
                   primers = c(rep(rep(c("oRG46_oRG47", "oRG54_oRG55", "oRG50_oRG51", "oRG52_oRG53", "oRG54_oRG55", "oRG50_oRG51", "oRG52_oRG53"), rep(3, 7)), 16), # samples
                              rep(c("oRG46_oRG47", "oRG54_oRG55", "oRG50_oRG51", "oRG52_oRG53"), rep(3, 4))), # NTCs
                   Cq = 0)

# sanity checks
nrow(qPCR) == 21*16+12
head(qPCR)

# read Cq Results file ----------------------------------------------------
file.path <- grep(pattern = "Cq", x = list.files(path = raw_data_dir, full.names = TRUE, recursive = TRUE), value = TRUE)
tmp <- read.table(file = file.path, sep = "\t", header = TRUE)
  
# save Cq values in qPCR data frame
for (well in qPCR$well){
  qPCR$Cq[qPCR$well == well] <- tmp$Cq[tmp$Well == well & tmp$Fluor == "SYBR"]
}

head(qPCR)

# check NTCs --------------------------------------------------------------

# min NTC Cq
min(qPCR$Cq[qPCR$sample == "NTC"], na.rm = TRUE)
# 31.96858

# max non-NTC Cq
max(qPCR$Cq[qPCR$sample != "NTC"], na.rm = TRUE)
# 26.93768

# check non-NTC samples with high Cq values
qPCR[qPCR$Cq > 25 & qPCR$sample != "NTC",]
# RsaI digests and HO cutting

# calc Cq mean, sd, and cv ------------------------------------------------
# options(scipen = 999) # switch off scientific notation 
# options(scipen = 0) # switch on scientific notation 

Cq.modes <- aggregate(Cq ~ sample + digest + time + primers, data = qPCR,
                      FUN = function(x) c(mean = mean(x, na.rm = TRUE),
                                          sd = sd(x, na.rm = TRUE),
                                          cv = sd(x, na.rm = TRUE)/mean(x, na.rm = TRUE)
                                          )
                      )

# sort
Cq.modes <- Cq.modes[order(Cq.modes$sample, Cq.modes$primers, Cq.modes$time),]

# mean, sd, and cv are in a single column (as a matrix)
# convert into separate columns
Cq.modes <- data.frame(Cq.modes[, 1:(ncol(Cq.modes)-1)], as.data.frame(Cq.modes$Cq))

# what is the largest cv?
max(Cq.modes$cv, na.rm = TRUE)
# 0.02936018

# which are the samples with large cv?
Cq.modes[Cq.modes$cv > 0.025, ]
# only NTC

# what is the ADH1 fold spread?
# ADH1 Cq mean values
tmp <- Cq.modes$mean[Cq.modes$primers == "oRG54_oRG55" & Cq.modes$sample != "NTC"]
# fold spread
efficiencies$efficiency[efficiencies$primers == "oRG54_oRG55"]^(diff(range(tmp)))
# 2.357728

# what is the ADH1 fold spread within each strain
aggregate(mean ~ sample, data = Cq.modes[Cq.modes$primers == "oRG54_oRG55" & Cq.modes$sample != "NTC", ],
          FUN = function(x) efficiencies$efficiency[efficiencies$primers == "oRG54_oRG55"]^(diff(range(x))))
#     sample     mean
# 1  4376-3C 1.605572
# 2  4376-4A 1.244854
# 3 4377-12B 1.755001
# 4 4377-15A 1.944129


# calc HOcs cut fraction over time ========================================
cut.primers <- "oRG46_oRG47"
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
  Cut.fraction <- rbind(Cut.fraction, cbind(strain = strain, primers = cut.primers, tmp, stringsAsFactors = FALSE))  
}

# write data to file
write.table(x = Cut.fraction, file = paste0(moment_dir, "/Cutting.txt"), row.names = FALSE, col.names = TRUE)


# plot cutting ================================================================
strains <- c("4376-3C", "4376-4A", "4377-12B", "4377-15A")
MyColors <- JFly_colors[1:length(strains)]

# print to PDF
pdf(file = "tmp.pdf", width=3.5, height=3)
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(1, 0, 4, 2), las=1)

# start empty plot
y.range <- range(c(Cut.fraction$mean + Cut.fraction$sd, Cut.fraction$mean - Cut.fraction$sd, 0, 1))
# y.range <- c(-0.5, 1.5)
plot(x = NA, y = NA, ylim = y.range, xlim = range(Cut.fraction$time), xlab = "Time [h]", ylab = "Cut Fraction")

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

dev.off()

GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/Cut_fraction.pdf"))


# plot legend =================================================================
LegTxt <- c(expression("wt"),
            expression(italic("exo1")*Delta~italic("sgs1")*Delta),
            expression(italic("exo1")*Delta~italic("sgs1-aa dna2-aa")),
            expression(italic("exo1")*Delta~italic("sgs1-aa dna2-aa mre11-H125N"))
)

pdf(file = "tmp.pdf", width=3.4, height=1.05)
par(cex = 1, mar = rep(0, 4))
plot(1, type="n", axes=FALSE, xlab="", ylab="")
legend(1, 1, xjust=0.5, yjust=0.5, legend = LegTxt, 
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
  for (resect.primers in c("oRG50_oRG51", "oRG52_oRG53")){
    
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
strains <- c("4376-3C", "4376-4A", "4377-12B", "4377-15A")
MyColors <- JFlyColors[1:length(strains)]

# +98 bp ---------------------------------------------------------------------
tmp <- Resected.fraction[Resected.fraction$primers == "oRG50_oRG51", ]

# print to PDF
pdf(file = "tmp.pdf", width=3.5, height=3)
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(1, 0, 4, 2), las=1)

# start empty plot
y.range <- range(c(tmp$mean + tmp$sd, tmp$mean - tmp$sd, 0, 1), na.rm = TRUE)
# y.range <- c(-0.5, 2)
plot(x = NA, y = NA, ylim = y.range, xlim = range(Resected.fraction$time), xlab = "Time [h]", ylab = "ssDNA Fraction")

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

GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/Resection_98_bp.pdf"))


# +640 bp ---------------------------------------------------------------------
tmp <- Resected.fraction[Resected.fraction$primers == "oRG52_oRG53", ]

# print to PDF
pdf(file = "tmp.pdf", width=3.5, height=3)
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(1, 0, 4, 2), las=1)

# start empty plot
y.range <- range(c(tmp$mean + tmp$sd, tmp$mean - tmp$sd, 0, 1), na.rm = TRUE)
# y.range <- c(-0.5, 2)
plot(x = NA, y = NA, ylim = y.range, xlim = range(Resected.fraction$time), xlab = "Time [h]", ylab = "ssDNA Fraction")

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

GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/Resection_640_bp.pdf"))
