# info --------------------------------------------------------------------
# purpose: Extract and gate flow cytometry data from FCS files
# author: Robert Gnuegge (robert.gnuegge@gmail.com)
# date: 02/26/22
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
folder.name <- "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Data/Flow_cytometry/LSY4810-41D_4811-11B_4812-5A/Flow_cytometry_data"
gating_dir <- "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Figures/19_Flow_cytometry/LSY4810-41D_4811-11B_4812-5A/Gating"

# extract fcs.file.names
fcs.file.names <- list.files(path = folder.name)
fcs.file.names

# sample.names <- fcs.file.names and some cosmetics
sample.names <- sub(pattern = "Measurements_", replacement = "", x = fcs.file.names)
sample.names <- sub(pattern = ".fcs", replacement = "", x = sample.names)
sample.names

# read files into flowSet
fcs.data <- read.flowSet(path = folder.name)
sampleNames(fcs.data) <- sample.names
colnames(fcs.data)

# =============================================================================
# gating

tmp <- fcs.data

xyplot(`SSC-H` ~ `FSC-H`, data = tmp[3])

# plot SSC-H vs. FSC-H overlayed for all samples ------------------------------

# extract all FSC-H and SSC-H data
all.fsc.ssc <- data.frame(fsApply(x = tmp[, c("FSC-H", "SSC-H")], FUN = exprs, simplify = TRUE))

# plot to png -----------------------------------------------------------------

plot.name <- paste0(gating_dir, "/all_SSC-H_vs_FSC-H.png")
png(plot.name, width = 500, height = 500)
par(mar = rep(0,4))  # suppress outer margins

idx <- seq(from = 1, to = nrow(all.fsc.ssc), by = 1)  # plot only every "nth" element

plot(x = all.fsc.ssc$FSC.H[idx], y = all.fsc.ssc$SSC.H[idx], xlim = c(0, 200000), ylim=c(0, 200000),
     pch = ".", col = adjustcolor("black", alpha.f = 0.1), axes = FALSE, xlab = NA, ylab = NA)  # suppress axis for easier measurement in ImageJ

# add axis lines to measure scale
points(x = c(0,200000), y = c(0,0), type = "l")
points(x = c(0,0), y = c(0,200000), type = "l")

dev.off()

# get parameters of ellipsoid gate using ImageJ -------------------------------

# 0. (if necessary) Analyze -> Set Measurments -> select "Center of mass" and "Fit ellipse"
# 1. open PNG image
# 2. place elliptical selection
# 3. measure
# 4. save results as "ElliGate.txt"

# construct gating filter object ----------------------------------------------
elli.data <- read.table(file = paste0(gating_dir, "/ElliGate.txt"), header = TRUE)

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
pdf(file = paste0(gating_dir, "/Gate_on_all_SSC-H_vs_FSC-H.pdf"), width = 8.5, height = 11)
print(xyplot(x = `SSC-H` ~ `FSC-H`, data = tmp, filter = MyGate, strip = strip.custom(factor.levels = sampleNames(tmp))))
dev.off()
  
# gate and save data ----------------------------------------------------------
gated <- Subset(x = tmp, subset = MyGate)

gated.df <- data.frame()  # initialize collection df

# iterate through all samples
for (strain in sampleNames(tmp)){
  
  Citrine <- exprs(gated[[strain, "FITC-H"]])  # get Citrine data
  mKate2 <- exprs(gated[[strain, "PE-Texas Red-H"]])  # get mKate2 data
  tmp.df <- data.frame(sample = strain, Citrine = Citrine[, 1], mKate2 = mKate2[, 1])  # construct df for sample
  gated.df <- rbind(gated.df, tmp.df)  # append to collection df
  
}

str(gated.df)

# save data
write.table(x = gated.df, file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Data/Flow_cytometry/LSY4810-41D_4811-11B_4812-5A/Extracted_fluorescence_data.txt")
