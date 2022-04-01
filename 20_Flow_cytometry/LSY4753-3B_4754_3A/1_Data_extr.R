# info --------------------------------------------------------------------
# purpose: Extract and gate flow cytometry data from FCS files
# author: Robert Gnuegge (robert.gnuegge@gmail.com)
# date: 02/25/22
# version: 1.0

# preamble ----------------------------------------------------------------

# set wd to this file's location
wd_path <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(wd_path)

# load helper files and functions
source(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Analysis_scripts/Misc/Misc_helper_functions.R")

# load libraries
library(flowCore)
library(flowViz)

# read the fcs files ----------------------------------------------------------
folder.name <- "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Data/Flow_cytometry/LSY4753-3B_4754-3A/Flow_cytometry_data/"
gating_dir <- "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Figures/19_Flow_cytometry/LSY4753-3B_4754-3A/Gating"

# extract fcs.file.names
fcs.file.names <- list.files(path = folder.name)

# read files into flowSet
fcs.data <- read.flowSet(path = folder.name)
sampleNames(fcs.data) <- sub(pattern = ".fcs", replacement = "", x = sampleNames(fcs.data))
colnames(fcs.data)

# plot SSC-H vs. FSC-H overlayed for all samples ------------------------------

# extract all FSC-H and SSC-H data
for(gal in c("Gal_0_01", "Gal_2")){

  idx <- grep(pattern = gal, x = sampleNames(fcs.data))
  all.fsc.ssc <- data.frame(fsApply(x = fcs.data[idx, c("FSC-H", "SSC-H")], FUN = exprs, simplify = TRUE))

  # plot to png -----------------------------------------------------------------
  plot.name <- paste0(gating_dir, "/all_SSC-H_vs_FSC-H_", gal, ".png")
  png(plot.name, width = 500, height = 500)
  par(mar = rep(0,4))  # suppress outer margins

  idx <- seq(from = 1, to = nrow(all.fsc.ssc), by = 2)  # plot only every "nth" element

  plot(x = all.fsc.ssc$FSC.H[idx], y = all.fsc.ssc$SSC.H[idx], xlim = c(0, 200000), ylim=c(0, 200000),
       pch = ".", col = adjustcolor("black", alpha.f = 0.1), axes = FALSE, xlab = NA, ylab = NA)  # suppress axis for easier measurement in ImageJ

  # add axis lines to measure scale
  points(x = c(0,200000), y = c(0,0), type = "l")
  points(x = c(0,0), y = c(0,200000), type = "l")

  dev.off()

}

# get parameters of ellipsoid gate using ImageJ -------------------------------

# 0. (if necessary) Analyze -> Set Measurments -> select "Center of mass" and "Fit ellipse"
# 1. open PNG image
# 2. place elliptical selection
# 3. measure
# 4. save results as "ElliGate.txt"

# construct gating filter object ----------------------------------------------
elli.data <- read.table(file = paste0(gating_dir, "/ElliGate_Gal_0_01.txt"), header = TRUE)

x.center <- elli.data$XM
y.center <- elli.data$YM
major.axis <- elli.data$Major
minor.axis <- elli.data$Minor
h.min.ax <- minor.axis / 2
angle <- elli.data$Angle

# define measures for scaling
real.dis <- 463
ass.dis <- 200000
real.x.ori <- 18
real.y.ori <- 481
ass.x.ori <- 0
ass.y.ori <- 0 	

# calculate center
f.center.x <- ass.x.ori + (x.center - real.x.ori) * ass.dis / real.dis
f.center.y <- ass.y.ori + (real.y.ori - y.center) * ass.dis / real.dis
center <- c("FSC-H" = f.center.x, "SSC-H" = f.center.y)

# calc inverse of covariance matrix
h.maj.ax <- major.axis / 2 * ass.dis / real.dis
h.min.ax <- minor.axis / 2 * ass.dis / real.dis
angle.rad <- angle * pi / 180

# calc matrix elements
m.1.1 <- cos(angle.rad)^2 / h.maj.ax^2 + sin(angle.rad)^2 / h.min.ax^2
m.1.2 <- sin(angle.rad) * cos(angle.rad) *  (1 / h.maj.ax^2 - 1/h.min.ax^2)
m.2.1 <- m.1.2
m.2.2 <- sin(angle.rad)^2 / h.maj.ax^2 + cos(angle.rad)^2 / h.min.ax^2

# calc inverse and non-inverse matrix
m.i <- matrix(data = c(m.1.1, m.1.2, m.2.1, m.2.2), nrow=2,
              dimnames = list(c("FSC-H", "SSC-H"), c("FSC-H", "SSC-H")))
m <- solve(m.i)

# define the gating filter
MyGate <- ellipsoidGate(.gate = m, mean = center)

# plot gating filter on all SSC-H vs FSC-H plots
for (n in 1:2){
  
  idx <- (1 + (n - 1) * 10):(n * 10)
  pdf(file = paste0(gating_dir, "/Gate_on_all_SSC-H_vs_FSC-H_Gal_0_01_", n, ".pdf"), width = 8.5, height = 11)
  print(xyplot(x = `SSC-H` ~ `FSC-H`, data = fcs.data[idx], filter = MyGate, strip = strip.custom(factor.levels = sampleNames(fcs.data)[idx])))
  dev.off()
  
}

# gate and save data ----------------------------------------------------------
gated <- Subset(x = fcs.data[grep(pattern = "Gal_0_01", x = sampleNames(fcs.data))], subset = MyGate)
sampleNames(gated)

# construct gating filter object ----------------------------------------------
elli.data <- read.table(file = paste0(gating_dir, "/ElliGate_Gal_2.txt"), header = TRUE)

x.center <- elli.data$XM
y.center <- elli.data$YM
major.axis <- elli.data$Major
minor.axis <- elli.data$Minor
h.min.ax <- minor.axis / 2
angle <- elli.data$Angle

# define measures for scaling
real.dis <- 463
ass.dis <- 200000
real.x.ori <- 18
real.y.ori <- 481
ass.x.ori <- 0
ass.y.ori <- 0 	

# calculate center
f.center.x <- ass.x.ori + (x.center - real.x.ori) * ass.dis / real.dis
f.center.y <- ass.y.ori + (real.y.ori - y.center) * ass.dis / real.dis
center <- c("FSC-H" = f.center.x, "SSC-H" = f.center.y)

# calc inverse of covariance matrix
h.maj.ax <- major.axis / 2 * ass.dis / real.dis
h.min.ax <- minor.axis / 2 * ass.dis / real.dis
angle.rad <- angle * pi / 180

# calc matrix elements
m.1.1 <- cos(angle.rad)^2 / h.maj.ax^2 + sin(angle.rad)^2 / h.min.ax^2
m.1.2 <- sin(angle.rad) * cos(angle.rad) *  (1 / h.maj.ax^2 - 1/h.min.ax^2)
m.2.1 <- m.1.2
m.2.2 <- sin(angle.rad)^2 / h.maj.ax^2 + cos(angle.rad)^2 / h.min.ax^2

# calc inverse and non-inverse matrix
m.i <- matrix(data = c(m.1.1, m.1.2, m.2.1, m.2.2), nrow=2,
              dimnames = list(c("FSC-H", "SSC-H"), c("FSC-H", "SSC-H")))
m <- solve(m.i)

# define the gating filter
MyGate <- ellipsoidGate(.gate = m, mean = center)

# plot gating filter on all SSC-H vs FSC-H plots
idx <- grep(pattern = "Gal_2", x = sampleNames(fcs.data))

for (n in 1:2){
  
  idx <- (1 + (n - 1) * 10):(n * 10)
  pdf(file = paste0(gating_dir, "/Gate_on_all_SSC-H_vs_FSC-H_Gal_2_", n, ".pdf"), width = 8.5, height = 11)
  print(xyplot(x = `SSC-H` ~ `FSC-H`, data = fcs.data[idx], filter = MyGate, strip = strip.custom(factor.levels = sampleNames(fcs.data)[idx])))
  dev.off()
  
}

# gate and save data ----------------------------------------------------------
tmp <- Subset(x = fcs.data[grep(pattern = "Gal_2", x = sampleNames(fcs.data))], subset = MyGate)
sampleNames(tmp)

gated <- rbind2(x = gated, y = tmp)
sampleNames(gated)

# extract gated Citrine values ------------------------------------------------
gated.Cit.df <- data.frame()

for(gal in c(0.01, 2)){
  
  for (strain in c("4436", "4518", "4753", "4754")){
    
    for (t in c(0, 2, 4, 6, 8)){
      
      name <- paste0("Gal_", sub(pattern = "\\.", replacement = "_", x = gal), "_", strain, "_", t)  # sample name
      tmp <- exprs(gated[[name, "FITC-H"]])  # get Citrine data
      tmp.df <- data.frame(gal = gal, strain = strain, time = t, Citrine = tmp[, 1], stringsAsFactors = FALSE)  # construct df for sample
      gated.Cit.df <- rbind(gated.Cit.df, tmp.df)  # append to collection df
      
    }
    
  }
  
}

str(gated.Cit.df)

# save data
write.table(x = gated.Cit.df, file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Data/Flow_cytometry/LSY4753-3B_4754-3A/Extracted_Citrine_data.txt")
