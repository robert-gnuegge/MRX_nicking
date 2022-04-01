###############################################################################
# Functions to analyze resection assay qPCR data

# necessary libraries
require("propagate")

# CalcCutFraction2 function ---------------------------------------------------
CalcResectedFraction2 <- function(resect.mock, resect.digest, resect.primer.eff = 2, ref.mock, ref.digest, ref.primer.eff = 2, cut.fraction, t){
  
  # check that all necessary data are supplied
  if (is.null(resect.mock) | is.null(resect.digest) | is.null(ref.mock) | is.null(ref.digest) | is.null(cut.fraction)){
    stop("resect.mock, resect.digest, ref.mock, ref.digest, and cut.fraction data need to be supplied.")
  }
  
  # check that time is consistent between all supplied data frames
  t <- resect.mock$time
  t2 <- resect.digest$time
  t3 <- ref.mock$time
  t4 <- ref.digest$time
  t5 <- cut.fraction$time
  if(any(sort(t) != sort(t2), sort(t) != sort(t3), sort(t) != sort(t4), sort(t) != sort(t5))){
    stop("Time points are not consistent among all supplied data frames. Supply time points to be evaluated as function argument!")
  }
  # check if zero time point is supplied
  if (!(0 %in% t)){
    stop("No zero time point in supplied data.")
  }
  
  # create resected.fraction expression
  resected.fraction.expr <- as.expression(substitute(expr = 2 / ((resect.primer.eff^(resect.digest - resect.mock) / ref.primer.eff^(ref.digest - ref.mock) + 1) * cut.fraction),
                                                     env = list(resect.primer.eff = resect.primer.eff, ref.primer.eff = ref.primer.eff)
                                                     )
                                          )
  
  # initialize result df
  out <- data.frame()
  
  # calc resected fraction for t
  for (t in t){
    
    # calc mean and sd of cut fraction with error propagation
    out.propagate <- propagate(expr = resected.fraction.expr,
                               data = cbind(ref.digest = as.numeric(ref.digest[ref.digest$time == t, c("mean", "sd")]),
                                            ref.mock = as.numeric(ref.mock[ref.mock$time == t, c("mean", "sd")]),
                                            resect.digest = as.numeric(resect.digest[resect.digest$time == t, c("mean", "sd")]),
                                            resect.mock = as.numeric(resect.mock[resect.mock$time == t, c("mean", "sd")]),
                                            cut.fraction = as.numeric(cut.fraction[cut.fraction$time == t, c("mean", "sd")])
                                            ),
                               do.sim = FALSE)
    
    # save results to result df
    out <- rbind(out, data.frame(time = t, 
                                 mean = out.propagate$prop["Mean.1"], sd = out.propagate$prop["sd.1"], 
                                 conf.2.5 = out.propagate$prop["2.5%"], conf.97.5 = out.propagate$prop["97.5%"])
                 )
    
  }
  
  # return result df
  out$mean[out$time == 0 & out$mean == Inf] <- 0  # replace Inf with 0 for t = 0
  out <- data.frame(out, row.names = NULL, stringsAsFactors = FALSE)  # remove row.names
  return(out)
  
}