# info --------------------------------------------------------------------
# purpose: miscellaneous helper functions
# author: Robert Gnuegge (robert.gnuegge@gmail.com)
# date: 12/16/21
# version: 0.1


# convert dash to underscore ----------------------------------------------
# argument: string
# result: string
dash_to_underscore <- function(string){
  gsub(pattern = "-", replacement = "_", x = string)
}


# embed fonts in PDF (using ghostscript) ----------------------------------
# argument: PDF input and output file names (strings)
# result: none
GS_embed_fonts <- function(input, output, remove_input_file = TRUE){
  success <- system2(command = "gs", args = paste0("-q -dNOPAUSE -dBATCH -dPDFSETTINGS=/prepress -sDEVICE=pdfwrite -sOutputFile=", output, " ", input))
  if(remove_input_file & success == 0){  # success == 0 only if system2 command ran successfully
    file.remove(input)
  }
}


# moving average (median) ---------------------------------------------------
# argument: numeric vector, integer
# result: numeric vector
# note: runmed, except for first n = keep elements, where original element values are used
moving_average <- function(x, k, keep = 2){
  x_length <- length(x)
  out <- rep(0, x_length)
  if(keep > 0){
    out[1:keep] <- x[1:keep]
  }
  out[(keep + 1):x_length] <- runmed(x = x[(keep + 1):x_length], k = k)
  return(out)
}


# moving mean -------------------------------------------------------------
# argument: numeric vector, integer
# result: numeric vector
# note: rolling mean, except for first and last elements (defined by keep), where original element values are used
# note: for elements near beginning and end of x, value of k is adjusted
moving_mean <- function(x, k, keep = 3){
  require(zoo)
  out <- rollmean(x = x, k = k, fill = NA)
  if(keep > 0){
    out[1:keep] <- x[1:keep]
    out[(length(x) - keep + 1):length(x)] <- x[(length(x) - keep + 1):length(x)]
  }
  for(i in (keep + 1):(floor(k / 2) + keep)){
    out[i] <- mean(x[(keep + 1):(2 * i - keep - 1)])
    out[(length(out) - i + 1)] <- mean(x[length(out) - (keep):(2 * i - keep - 2)])
  }
  return(out)
}


# transform data for axis break with scale change -------------------------
# argument: numeric vector
# result: numeric vector
trans <- function(x, threshold, threshold_trans, factor){
  pmin(x, threshold) + factor * pmax(x - threshold_trans, 0)
}
