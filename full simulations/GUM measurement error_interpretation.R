

##### DATA ####

## true parameters

true_p95 <- 100

true_gsd <- 2.5

true_gm <- exp( log(true_p95) - qnorm(0.95)*log(true_gsd) )


## simulation parameters

expanded_uncertainty <- 0.50

coverage_factor <- qnorm(0.975)

sample_size <- c(6,3,6,3)

gsd <- c(2.5,2.5,1.5,1.5)

simulation_summary1 <- readRDS("C:/jerome/Dropbox/temp/aioh2023-S2_1.RDS")
simulation_summary2 <- readRDS("C:/jerome/Dropbox/temp/aioh2023-S2_2.RDS")
simulation_summary3 <- readRDS("C:/jerome/Dropbox/temp/aioh2023-S2_3.RDS")
simulation_summary4 <- readRDS("C:/jerome/Dropbox/temp/aioh2023-S2_4.RDS")

    
#### FUNCTIONS ####

## relative rmse results function

rmse.result <- function( simulation_summary , true_p95 , true_gsd , true_gm ) { 
  
  
  rmse_table <- data.frame( method = c("ideal_F" , "ideal_B" ,
                                       "naive_F" , "naive_B",
                                       "me_F" , "me_B"),
                            gm = numeric(6),
                            gsd = numeric(6),
                            p95 = numeric(6))
  
  rmse_table$gm[1] <- compute_relrmse( simulation_summary$ideal_gm_freq , true_gm ) 
  rmse_table$gm[2] <- compute_relrmse( simulation_summary$ideal_gm , true_gm )
  rmse_table$gm[3] <- compute_relrmse( simulation_summary$naive_gm_freq , true_gm )
  rmse_table$gm[4] <- compute_relrmse( simulation_summary$naive_gm , true_gm )
  rmse_table$gm[5] <- compute_relrmse( simulation_summary$gum_gm , true_gm )
  rmse_table$gm[6] <- compute_relrmse( simulation_summary$me_gm , true_gm )
  
  rmse_table$gsd[1] <- compute_relrmse( simulation_summary$ideal_gsd_freq , true_gsd )
  rmse_table$gsd[2] <- compute_relrmse( simulation_summary$ideal_gsd , true_gsd )
  rmse_table$gsd[3] <- compute_relrmse( simulation_summary$naive_gsd_freq , true_gsd )
  rmse_table$gsd[4] <- compute_relrmse( simulation_summary$naive_gsd , true_gsd )
  rmse_table$gsd[5] <- compute_relrmse( simulation_summary$gum_gsd , true_gsd )
  rmse_table$gsd[6] <- compute_relrmse( simulation_summary$me_gsd , true_gsd )
  
  rmse_table$p95[1] <- compute_relrmse( simulation_summary$ideal_p95_freq , true_p95 )
  rmse_table$p95[2] <- compute_relrmse( simulation_summary$ideal_p95 , true_p95 )
  rmse_table$p95[3] <- compute_relrmse( simulation_summary$naive_p95_freq , true_p95 )
  rmse_table$p95[4] <- compute_relrmse( simulation_summary$naive_p95 , true_p95 )
  rmse_table$p95[5] <- compute_relrmse( simulation_summary$gum_p95 , true_p95 )
  rmse_table$p95[6] <- compute_relrmse( simulation_summary$me_p95 , true_p95 )
  
  return(rmse_table)
  
  }
  
  
## relative precision results function

precision.result <- function( simulation_summary , true_p95 , true_gsd , true_gm ) { 
  
  
  precision_table <- data.frame( method = c("ideal_F" , "ideal_B" ,
                                       "naive_F" , "naive_B",
                                       "me_F" , "me_B"),
                            gm = numeric(6),
                            gsd = numeric(6),
                            p95 = numeric(6))
  
  precision_table$gm[1] <- compute_relprecision( simulation_summary$ideal_gm_freq , true_gm ) 
  precision_table$gm[2] <- compute_relprecision( simulation_summary$ideal_gm , true_gm )
  precision_table$gm[3] <- compute_relprecision( simulation_summary$naive_gm_freq , true_gm )
  precision_table$gm[4] <- compute_relprecision( simulation_summary$naive_gm , true_gm )
  precision_table$gm[5] <- compute_relprecision( simulation_summary$gum_gm , true_gm )
  precision_table$gm[6] <- compute_relprecision( simulation_summary$me_gm , true_gm )
  
  precision_table$gsd[1] <- compute_relprecision( simulation_summary$ideal_gsd_freq , true_gsd )
  precision_table$gsd[2] <- compute_relprecision( simulation_summary$ideal_gsd , true_gsd )
  precision_table$gsd[3] <- compute_relprecision( simulation_summary$naive_gsd_freq , true_gsd )
  precision_table$gsd[4] <- compute_relprecision( simulation_summary$naive_gsd , true_gsd )
  precision_table$gsd[5] <- compute_relprecision( simulation_summary$gum_gsd , true_gsd )
  precision_table$gsd[6] <- compute_relprecision( simulation_summary$me_gsd , true_gsd )
  
  precision_table$p95[1] <- compute_relprecision( simulation_summary$ideal_p95_freq , true_p95 )
  precision_table$p95[2] <- compute_relprecision( simulation_summary$ideal_p95 , true_p95 )
  precision_table$p95[3] <- compute_relprecision( simulation_summary$naive_p95_freq , true_p95 )
  precision_table$p95[4] <- compute_relprecision( simulation_summary$naive_p95 , true_p95 )
  precision_table$p95[5] <- compute_relprecision( simulation_summary$gum_p95 , true_p95 )
  precision_table$p95[6] <- compute_relprecision( simulation_summary$me_p95 , true_p95 )

  
  return(precision_table)
  
}


## relative bias results function

bias.result <- function( simulation_summary , true_p95 , true_gsd , true_gm ) { 
  
  
  bias_table <- data.frame( method = c("ideal_F" , "ideal_B" ,
                                       "naive_F" , "naive_B",
                                       "me_F" , "me_B"),
                            gm = numeric(6),
                            gsd = numeric(6),
                            p95 = numeric(6))
  
  bias_table$gm[1] <- compute_relbias( simulation_summary$ideal_gm_freq , true_gm ) 
  bias_table$gm[2] <- compute_relbias( simulation_summary$ideal_gm , true_gm )
  bias_table$gm[3] <- compute_relbias( simulation_summary$naive_gm_freq , true_gm )
  bias_table$gm[4] <- compute_relbias( simulation_summary$naive_gm , true_gm )
  bias_table$gm[5] <- compute_relbias( simulation_summary$gum_gm , true_gm )
  bias_table$gm[6] <- compute_relbias( simulation_summary$me_gm , true_gm )
  
  bias_table$gsd[1] <- compute_relbias( simulation_summary$ideal_gsd_freq , true_gsd )
  bias_table$gsd[2] <- compute_relbias( simulation_summary$ideal_gsd , true_gsd )
  bias_table$gsd[3] <- compute_relbias( simulation_summary$naive_gsd_freq , true_gsd )
  bias_table$gsd[4] <- compute_relbias( simulation_summary$naive_gsd , true_gsd )
  bias_table$gsd[5] <- compute_relbias( simulation_summary$gum_gsd , true_gsd )
  bias_table$gsd[6] <- compute_relbias( simulation_summary$me_gsd , true_gsd )
  
  bias_table$p95[1] <- compute_relbias( simulation_summary$ideal_p95_freq , true_p95 )
  bias_table$p95[2] <- compute_relbias( simulation_summary$ideal_p95 , true_p95 )
  bias_table$p95[3] <- compute_relbias( simulation_summary$naive_p95_freq , true_p95 )
  bias_table$p95[4] <- compute_relbias( simulation_summary$naive_p95 , true_p95 )
  bias_table$p95[5] <- compute_relbias( simulation_summary$gum_p95 , true_p95 )
  bias_table$p95[6] <- compute_relbias( simulation_summary$me_p95 , true_p95 )

  return(bias_table)
  
}



## coverage results function

coverage.result <- function( simulation_summary , true_p95 ) { 
  
  
  coverage_table <- data.frame( method = c("ideal_F" , "ideal_B" ,
                                       "naive_F" , "naive_B",
                                       "me_F" , "me_B"),
                            p95 = numeric(6))
  
  coverage_table$p95[1] <-  100*sum(simulation_summary$ideal_p95_ucl_freq>=true_p95)/length(simulation_summary$ideal_p95_ucl_freq)
  coverage_table$p95[2] <-  100*sum(simulation_summary$ideal_p95_ucl>=true_p95)/length(simulation_summary$ideal_p95_ucl)
  coverage_table$p95[3] <-  100*sum(simulation_summary$naive_p95_ucl_freq>=true_p95)/length(simulation_summary$naive_p95_ucl_freq)
  coverage_table$p95[4] <-  100*sum(simulation_summary$naive_p95_ucl>=true_p95)/length(simulation_summary$naive_p95_ucl)
  coverage_table$p95[5] <-  100*sum(simulation_summary$gum_p95_ucl>=true_p95)/length(simulation_summary$gum_p95_ucl)
  coverage_table$p95[6] <-  100*sum(simulation_summary$me_p95_ucl>=true_p95)/length(simulation_summary$me_p95_ucl)
  
  return(coverage_table)
  
}

#### RESULTS ####

## summary of coverage results (ideal = 70%)

coverage_table <- data.frame( approach = c("ideal_F" , "ideal_B" ,
                                           "naive_F" , "naive_B",
                                           "me_F" , "me_B"),
                              n6_2.5 = coverage.result( simulation_summary1 , true_p95 )$p95,
                              n3_2.5 = coverage.result( simulation_summary2 , true_p95 )$p95,
                              n6_1.5 = coverage.result( simulation_summary3 , true_p95 )$p95,
                              n3_1.5 = coverage.result( simulation_summary4 , true_p95 )$p95)

  
#### SUPPORT #####

compute_rmse = function(x, theta){
  
  mean.x = mean(x)
  
  N = length(x)
  
  if(is.null(N))
    stop("x must be a vector of minimum length 1")
  
  rmse = sqrt((mean.x - theta) ^ 2 + sum((x - mean.x) ^ 2) / (N - 1))
  return(rmse)
}

compute_relrmse = function(x, theta){
  
  mean.x = mean(x)
  
  N = length(x)
  
  if(is.null(N))
    stop("x must be a vector of minimum length 1")
  
  rmse = sqrt((mean.x - theta) ^ 2 + sum((x - mean.x) ^ 2) / (N - 1))/theta
  return(rmse)
}


compute_bias = function(x, theta){
  
  bias = mean(x) - theta 
  
  return(bias)
}


compute_relbias = function(x, theta){
  
  bias = 100 * (mean(x) - theta) / theta
  
  return(bias)
}


compute_relprecision = function(x, theta){
  
  sd.x = sd(x)
  
  rp = 100 * sd.x / theta
  
  return(rp)
}
