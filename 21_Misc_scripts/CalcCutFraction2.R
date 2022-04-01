###############################################################################
# Functions do analyze gDNA cutting qPCR data

# necessary libraries
require("propagate")

# CalcCutFraction2 function ----------------------------------------------------
CalcCutFraction2 <- function(cut, ref, cut.primer.eff = 2, ref.primer.eff = 2, t){
  
  # check that all necessary data are supplied
  if (is.null(cut) | is.null(ref)){
    stop("cut and ref data need to be supplied.")
  }
  
  # check that time is consistent between all supplied data frames
  if(any(sort(cut$time) != sort(ref$time))){
    stop("Time points are not consistent among supplied data frames. Supply time points to be evaluated as function argument!")
  }else{
    t <- cut$time
  }
  # check if zero time point is supplied
  if (!(0 %in% t)){
    stop("No zero time point in supplied data.")
  }
  
  # create cut.fraction expression
  cut.fraction.expr <- as.expression(substitute(expr = 1 - cut.primer.eff^(cut.0 - cut.t) / ref.primer.eff^(ref.0 - ref.t),
                                                env = list(cut.primer.eff = cut.primer.eff, ref.primer.eff = ref.primer.eff)
                                                )
                                     )

  # initialize result df
  out <- data.frame()
  
  # calc cut fraction for all rows in df
  for (t in t){
    
    # calc mean and sd of cut fraction with error propagation
    out.propagate <- propagate(expr = cut.fraction.expr,
                               data = cbind(cut.0 = as.numeric(cut[cut$time == 0, c("mean", "sd")]),
                                            cut.t = as.numeric(cut[cut$time == t, c("mean", "sd")]),
                                            ref.0 = as.numeric(ref[ref$time == 0, c("mean", "sd")]),
                                            ref.t = as.numeric(ref[ref$time == t, c("mean", "sd")])
                                            ),
                               do.sim = FALSE
                               )
    
    # save results to result df
    out <- rbind(out, data.frame(time = t, 
                                 mean = out.propagate$prop["Mean.1"], sd = out.propagate$prop["sd.1"], 
                                 conf.2.5 = out.propagate$prop["2.5%"], conf.97.5 = out.propagate$prop["97.5%"])
                 )
    
  }
  
  # return result df
  out <- data.frame(out, row.names = NULL)
  return(out)
  
}