# Script for performing a simulation study with real GSD : 5000 iterations per scenario

#  errors encountered with JAGS on certain scenarios (see debug section)

#  solved with STAN

# LIBRARIES -------------------------------------

library(rjags)
library(rstan)

# DATA ---------------------------

real_gsds <- readRDS( "created data/real_gsd_values.RDS")


# SCRIPTS ---------------------------

source("full simulations/TD EXIL 2024 measurement error impact/scripts/sample_analysis_functions.R")

source("full simulations/TD EXIL 2024 measurement error impact/scripts/simulation_interpretation_functions.R")

source("other/support_functions.R")

source("other/performance_metrics.R")

source("parameter estimation/bayesian/load.webexpo.SEG.functions.R")

source("data generation/SEG data simulation.R")


# SIMULATiON SPACE --------------------

scenarios <- expand.grid(
  true_exceedance_perc = c(0.1, 1, 3, 7, 10, 25),
  sample_size      = c(3L, 6L , 9L, 12L),
  proportion_censored        = c(0, 0.3, 0.6),
  stringsAsFactors = FALSE)

  ## scenarios for n=6 or 9, F=10% 
  my_scenarios <- scenarios[c(11,35,59,17,41,65),]

  # fixed parameters
  
  n_sim <- 500
  
  true_p95 <- 100
  
  index <- 65
  
  sample_size <- scenarios$sample_size[index]

# ANALYSES --------------------------------------------------------

## parameter vectors due to real GSDs
  
  #true_gsd <- rep(2.5, n_sim)

  true_gsd <- sample( real_gsds[ real_gsds>=quantile(real_gsds,0.025) & real_gsds<=quantile(real_gsds,0.975)] , n_sim , replace = TRUE )
  
  true_gm <- exp( log(true_p95) - qnorm(0.95)*log(true_gsd) )
  
  oel <- exp( qnorm(1 - scenarios$true_exceedance_perc[index]/100, mean = log(true_gm) , sd = log(true_gsd) ))
  
  loq <- signif(exp(qnorm(scenarios$proportion_censored[index], mean = log(true_gm) , sd = log(true_gsd) )) ,3)
  
  
  ## testing whether issue in vercorized generation of OEL, LOQ : NO
  
  oel_man <- numeric(n_sim)
  loq_man <- numeric(n_sim)
  
  for (i in 1:n_sim) { 
    
    oel_man[i] <- exp( qnorm(1 - scenarios$true_exceedance_perc[index]/100, mean = log(true_gm[i]) , sd = log(true_gsd[i]) ) )
    
    loq_man[i] <- exp(qnorm(scenarios$proportion_censored[index], mean = log(true_gm[i]) , sd = log(true_gsd[i]) ))
    
    }
  
  
  ## simulated exposure values - approach A
  
  simulated_data_object <- data.simulation.seg.gen(  true_gm = true_gm, 
                                                     true_gsd = true_gsd,
                                                     sample_size = sample_size ,
                                                           me_cv_inf = 0,
                                                           me_cv_sup = 0,
                                                           censor_level = loq )
  
  
  
  ## simulated exposure values - approach B
  
  simulated_data_object_prime <- matrix(nrow=sample_size, ncol=n_sim)
  
  for (i in 1:n_sim) {
    
    simulated_data_object_prime[,i] <-data.simulation.seg(true_gsd = true_gsd[i], 
                                                          true_gm = true_gm[i], 
                                                          n_sim = 1, 
                                                          sample_size = sample_size ,
                                                          me_cv_inf = 0.15,
                                                          me_cv_sup = 0.15,
                                                          censor_level = loq[i] )$true[,1] 
    
  }
  
  ## simulated exposure values - approach C
  
  simulated_data_object_ter <- matrix(nrow=sample_size, ncol=n_sim)
  
  censor_level = loq
  
  data_matrix <- matrix(nrow=9, ncol=n_sim)
  
  for (i in 1:n_sim) {
    
    obs_vec <- exp( rnorm( n = sample_size , mean = log(true_gm[i]) , sd = log(true_gsd[i]) ) )

    obs_vec[ obs_vec < exp(-15) ] <- exp(-15)
    obs_vec[ obs_vec > exp(15) ] <- exp(15)
    
      data_matrix[,i] 
    
      concentrations <- obs_vec
      
      concentrations[ concentrations < censor_level[i] ] <- paste("<", signif(censor_level[i],3) , sep="")
      
      data_matrix[,i] <- concentrations
      
      
  
  simulated_data_object_ter[,i]<-data_matrix[,i]
  
  }
  
  ## mean n censored : ## OK for B and C but NOT A
  
  mean(apply(simulated_data_object$true, 2, function(x){mean(is.na(as.numeric(x)))}))
  mean(apply(simulated_data_object_prime, 2, function(x){mean(is.na(as.numeric(x)))}))
  mean(apply(simulated_data_object_ter, 2, function(x){mean(is.na(as.numeric(x)))}))
  
  ## running analysis function in a loop
  
  mistake_a <- logical(n_sim)
  mistake_b <- logical(n_sim)
  mistake_c <- logical(n_sim)
  
  for (iter_i in 1:n_sim) {
  
  # data preparation
  
  true_data_a <- simulated_data_object$true[,iter_i]
  true_data_b <- simulated_data_object_prime[,iter_i]
  true_data_c <- simulated_data_object_ter[,iter_i]
  
  ## analysis
  
  res_ideal_a = expostats.naive( true_data_a , oel[iter_i])["F_70ucl"]
  res_ideal_b = expostats.naive( true_data_b , oel[iter_i])["F_70ucl"]
  res_ideal_c = expostats.naive( true_data_c , oel[iter_i])["F_70ucl"]
  
  ## decision
  
  mistake_a[iter_i] <- res_ideal_a < 5 
  mistake_b[iter_i] <- res_ideal_b < 5 
  mistake_c[iter_i] <- res_ideal_c < 5 
  
  }

  mean(mistake_a)
  mean(mistake_b)
  mean(mistake_c)


  