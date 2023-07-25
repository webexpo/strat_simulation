#' function that estimates the compliance probability based on the expostats bayesian model 
#' 
#' compliance is defined as : the X% upper confidence limit on the 95th percentile is below the OEL
#' 
#' inputs include true GSD, true 95th percentile relative to OEL, level of confidence and proportion of non detects
#' 
#' 

##### Libraries #####  

library(rjags)

##### External scripts #####

source("data generation/load.webexpo.data_gen.functions.R")
source("parameter estimation/bayesian/load.webexpo.SEG.functions.R")

## data generation

sample1 <-  webexpo.seg.gener.LN <-function(
    # Sample size
  n = 10,
  # indicator for the presence of censoring
  no.censoring = FALSE,
  # Percentage Left censored 
  perc.lowerthan = 20,
  # Percentage right censored
  perc.greaterthan = 10,
  # Percentage interval censored
  perc.between = 20,
  # Distributional parameters
  gm = 0.3,
  gsd = 2.5,
  # Censoring parameters
  left_factor = 1.5,
  right_factor = 1/1.5,
  int_factor = 1.5,
  # Measurement error structure
  error = "none",   #(or "cv" or "sd")
  me.cv = 0.20 ,
  me.sd = 0.05)