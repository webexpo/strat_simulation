#' exploring the impact of measurement error industrial hygiene measurement data interpretation

#' question : for one example of true distribution, what is the impact of a 50% expanded measurement error (CV=50/1.96) on P95 and its UCL

#' This is an extension of the AIOH2023.r analysis, triggered by exchanges with Theo Scheffers who described a GUM approach to measurement error handling

#' This script presents the interpretation of the simulation, which was run 5 times


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

source("other/support_functions.R")

setwd("C:/jerome/Dropbox/GITHUB/WEBEXPO/sampling_strats/GUM measurement error 2024")

simulation_summary1 <- list ( sa = readRDS("aioh2023-S2_1.RDS"),
                              sb = readRDS("aioh2023-S2_1b.RDS"),
                              sc = readRDS("aioh2023-S2_1c.RDS"),
                              sd = readRDS("aioh2023-S2_1d.RDS"),
                              se = readRDS("aioh2023-S2_1e.RDS")
                                           )

simulation_summary2 <- list ( sa = readRDS("aioh2023-S2_2.RDS"),
                              sb = readRDS("aioh2023-S2_2b.RDS"),
                              sc = readRDS("aioh2023-S2_2c.RDS"),
                              sd = readRDS("aioh2023-S2_2d.RDS"),
                              se = readRDS("aioh2023-S2_2e.RDS")
                           )

simulation_summary3 <- list ( sa = readRDS("aioh2023-S2_3.RDS"),
                              sb = readRDS("aioh2023-S2_3b.RDS"),
                              sc = readRDS("aioh2023-S2_3c.RDS"),
                              sd = readRDS("aioh2023-S2_3d.RDS"),
                              se = readRDS("aioh2023-S2_3e.RDS")
                           )

simulation_summary4 <- list ( sa = readRDS("aioh2023-S2_4.RDS"),
                              sb = readRDS("aioh2023-S2_4b.RDS"),
                              sc = readRDS("aioh2023-S2_4c.RDS"),
                              sd = readRDS("aioh2023-S2_4d.RDS"),
                              se = readRDS("aioh2023-S2_4e.RDS")
                           )



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

###### coverage ######

coverage_matrix_n6_2.5 <- matrix( nrow = 6 , ncol = 5)

for (i in 1:5) { coverage_matrix_n6_2.5[,i] <- coverage.result( simulation_summary1[[i]] , true_p95 )$p95 }
  
coverage_matrix_n3_2.5 <- matrix( nrow = 6 , ncol = 5) 

for (i in 1:5) { coverage_matrix_n3_2.5[,i] <- coverage.result( simulation_summary2[[i]] , true_p95 )$p95 }
  
coverage_matrix_n6_1.5 <- matrix( nrow = 6 , ncol = 5)

for (i in 1:5) { coverage_matrix_n6_1.5[,i] <- coverage.result( simulation_summary3[[i]] , true_p95 )$p95 }

coverage_matrix_n3_1.5 <- matrix( nrow = 6 , ncol = 5)

for (i in 1:5) { coverage_matrix_n3_1.5[,i] <- coverage.result( simulation_summary4[[i]] , true_p95 )$p95 }
  
  
  
  coverage_table <- data.frame( approach = c("ideal_F" , "ideal_B" ,
                                           "naive_F" , "naive_B",
                                           "me_F" , "me_B"),
                              n6_2.5_mean =  apply( coverage_matrix_n6_2.5 , 1 , mean ),
                              n6_2.5_sd =  apply( coverage_matrix_n6_2.5 , 1 , sd ),
                              n3_2.5_mean =  apply( coverage_matrix_n3_2.5 , 1 , mean ),
                              n3_2.5_sd =  apply( coverage_matrix_n3_2.5 , 1 , sd ),
                              n6_1.5_mean =  apply( coverage_matrix_n6_1.5 , 1 , mean ),
                              n6_1.5_sd =  apply( coverage_matrix_n6_1.5 , 1 , sd ),
                              n3_1.5_mean =  apply( coverage_matrix_n3_1.5 , 1 , mean ),
                              n3_1.5_sd =  apply( coverage_matrix_n3_1.5 , 1 , sd ))

###### bias ######                                
  
       bias_array_n6_2.5 <- array( dim=c(6,3,5) )
  
      for (i in 1:5) { bias_array_n6_2.5[,,i] <- as.matrix(bias.result( simulation_summary1[[i]] , true_p95 , true_gsd , true_gm )[,2:4]) }
  
       
      bias_table_n6_2.5 <- data.frame( approach = c("ideal_F" , "ideal_B" ,
                                                 "naive_F" , "naive_B",
                                                 "me_F" , "me_B"),
                                    gm_mean =  apply( bias_array_n6_2.5[,1,] , 1 , mean ),
                                    gm_sd =  apply( bias_array_n6_2.5[,1,] , 1 , sd ),
                                    gsd_mean =  apply( bias_array_n6_2.5[,2,] , 1 , mean ),
                                    gsd_sd =  apply( bias_array_n6_2.5[,2,] , 1 , sd ),
                                    p95_mean =  apply( bias_array_n6_2.5[,3,] , 1 , mean ),
                                    p95_sd =  apply( bias_array_n6_2.5[,3,] , 1 , sd ))
      
      bias_array_n3_2.5 <- array( dim=c(6,3,5) )
      
      for (i in 1:5) { bias_array_n3_2.5[,,i] <- as.matrix(bias.result( simulation_summary2[[i]] , true_p95 , true_gsd , true_gm )[,2:4]) }
      
      bias_table_n3_2.5 <- data.frame( approach = c("ideal_F" , "ideal_B" ,
                                                 "naive_F" , "naive_B",
                                                 "me_F" , "me_B"),
                                    gm_mean =  apply( bias_array_n3_2.5[,1,] , 1 , mean ),
                                    gm_sd =  apply( bias_array_n3_2.5[,1,] , 1 , sd ),
                                    gsd_mean =  apply( bias_array_n3_2.5[,2,] , 1 , mean ),
                                    gsd_sd =  apply( bias_array_n3_2.5[,2,] , 1 , sd ),
                                    p95_mean =  apply( bias_array_n3_2.5[,3,] , 1 , mean ),
                                    p95_sd =  apply( bias_array_n3_2.5[,3,] , 1 , sd ))
  
      bias_array_n6_1.5 <- array( dim=c(6,3,5) )
      
      for (i in 1:5) { bias_array_n6_1.5[,,i] <- as.matrix(bias.result( simulation_summary3[[i]] , true_p95 , true_gsd , true_gm )[,2:4]) }
      
      bias_table_n6_1.5 <- data.frame( approach = c("ideal_F" , "ideal_B" ,
                                                 "naive_F" , "naive_B",
                                                 "me_F" , "me_B"),
                                    gm_mean =  apply( bias_array_n6_1.5[,1,] , 1 , mean ),
                                    gm_sd =  apply( bias_array_n6_1.5[,1,] , 1 , sd ),
                                    gsd_mean =  apply( bias_array_n6_1.5[,2,] , 1 , mean ),
                                    gsd_sd =  apply( bias_array_n6_1.5[,2,] , 1 , sd ),
                                    p95_mean =  apply( bias_array_n6_1.5[,3,] , 1 , mean ),
                                    p95_sd =  apply( bias_array_n6_1.5[,3,] , 1 , sd ))
      
      bias_array_n3_1.5 <- array( dim=c(6,3,5) )
      
      for (i in 1:5) { bias_array_n3_1.5[,,i] <- as.matrix(bias.result( simulation_summary4[[i]] , true_p95 , true_gsd , true_gm )[,2:4]) }
      
      bias_table_n3_1.5 <- data.frame( approach = c("ideal_F" , "ideal_B" ,
                                                 "naive_F" , "naive_B",
                                                 "me_F" , "me_B"),
                                    gm_mean =  apply( bias_array_n3_1.5[,1,] , 1 , mean ),
                                    gm_sd =  apply( bias_array_n3_1.5[,1,] , 1 , sd ),
                                    gsd_mean =  apply( bias_array_n3_1.5[,2,] , 1 , mean ),
                                    gsd_sd =  apply( bias_array_n3_1.5[,2,] , 1 , sd ),
                                    p95_mean =  apply( bias_array_n3_1.5[,3,] , 1 , mean ),
                                    p95_sd =  apply( bias_array_n3_1.5[,3,] , 1 , sd ))
      
      
###### precision ######
      
      precision_array_n6_2.5 <- array( dim=c(6,3,5) )
      
      for (i in 1:5) { precision_array_n6_2.5[,,i] <- as.matrix(precision.result( simulation_summary1[[i]] , true_p95 , true_gsd , true_gm )[,2:4]) }
      
      precision_table_n6_2.5 <- data.frame( approach = c("ideal_F" , "ideal_B" ,
                                                     "naive_F" , "naive_B",
                                                     "me_F" , "me_B"),
                                        gm_mean =  apply( precision_array_n6_2.5[,1,] , 1 , mean ),
                                        gm_sd =  apply( precision_array_n6_2.5[,1,] , 1 , sd ),
                                        gsd_mean =  apply( precision_array_n6_2.5[,2,] , 1 , mean ),
                                        gsd_sd =  apply( precision_array_n6_2.5[,2,] , 1 , sd ),
                                        p95_mean =  apply( precision_array_n6_2.5[,3,] , 1 , mean ),
                                        p95_sd =  apply( precision_array_n6_2.5[,3,] , 1 , sd ))
      
      precision_array_n3_2.5 <- array( dim=c(6,3,5) )
      
      for (i in 1:5) { precision_array_n3_2.5[,,i] <- as.matrix(precision.result( simulation_summary2[[i]] , true_p95 , true_gsd , true_gm )[,2:4]) }
      
      precision_table_n3_2.5 <- data.frame( approach = c("ideal_F" , "ideal_B" ,
                                                     "naive_F" , "naive_B",
                                                     "me_F" , "me_B"),
                                        gm_mean =  apply( precision_array_n3_2.5[,1,] , 1 , mean ),
                                        gm_sd =  apply( precision_array_n3_2.5[,1,] , 1 , sd ),
                                        gsd_mean =  apply( precision_array_n3_2.5[,2,] , 1 , mean ),
                                        gsd_sd =  apply( precision_array_n3_2.5[,2,] , 1 , sd ),
                                        p95_mean =  apply( precision_array_n3_2.5[,3,] , 1 , mean ),
                                        p95_sd =  apply( precision_array_n3_2.5[,3,] , 1 , sd ))
      
      precision_array_n6_1.5 <- array( dim=c(6,3,5) )
      
      for (i in 1:5) { precision_array_n6_1.5[,,i] <- as.matrix(precision.result( simulation_summary3[[i]] , true_p95 , true_gsd , true_gm )[,2:4]) }
      
      precision_table_n6_1.5 <- data.frame( approach = c("ideal_F" , "ideal_B" ,
                                                     "naive_F" , "naive_B",
                                                     "me_F" , "me_B"),
                                        gm_mean =  apply( precision_array_n6_1.5[,1,] , 1 , mean ),
                                        gm_sd =  apply( precision_array_n6_1.5[,1,] , 1 , sd ),
                                        gsd_mean =  apply( precision_array_n6_1.5[,2,] , 1 , mean ),
                                        gsd_sd =  apply( precision_array_n6_1.5[,2,] , 1 , sd ),
                                        p95_mean =  apply( precision_array_n6_1.5[,3,] , 1 , mean ),
                                        p95_sd =  apply( precision_array_n6_1.5[,3,] , 1 , sd ))
      
      precision_array_n3_1.5 <- array( dim=c(6,3,5) )
      
      for (i in 1:5) { precision_array_n3_1.5[,,i] <- as.matrix(precision.result( simulation_summary4[[i]] , true_p95 , true_gsd , true_gm )[,2:4]) }
      
      precision_table_n3_1.5 <- data.frame( approach = c("ideal_F" , "ideal_B" ,
                                                     "naive_F" , "naive_B",
                                                     "me_F" , "me_B"),
                                        gm_mean =  apply( precision_array_n3_1.5[,1,] , 1 , mean ),
                                        gm_sd =  apply( precision_array_n3_1.5[,1,] , 1 , sd ),
                                        gsd_mean =  apply( precision_array_n3_1.5[,2,] , 1 , mean ),
                                        gsd_sd =  apply( precision_array_n3_1.5[,2,] , 1 , sd ),
                                        p95_mean =  apply( precision_array_n3_1.5[,3,] , 1 , mean ),
                                        p95_sd =  apply( precision_array_n3_1.5[,3,] , 1 , sd ))
  
      
###### rmse ######
      
      rmse_array_n6_2.5 <- array( dim=c(6,3,5) )
      
      for (i in 1:5) { rmse_array_n6_2.5[,,i] <- as.matrix(rmse.result( simulation_summary1[[i]] , true_p95 , true_gsd , true_gm )[,2:4]) }
      
      rmse_table_n6_2.5 <- data.frame( approach = c("ideal_F" , "ideal_B" ,
                                                     "naive_F" , "naive_B",
                                                     "me_F" , "me_B"),
                                        gm_mean =  apply( rmse_array_n6_2.5[,1,] , 1 , mean ),
                                        gm_sd =  apply( rmse_array_n6_2.5[,1,] , 1 , sd ),
                                        gsd_mean =  apply( rmse_array_n6_2.5[,2,] , 1 , mean ),
                                        gsd_sd =  apply( rmse_array_n6_2.5[,2,] , 1 , sd ),
                                        p95_mean =  apply( rmse_array_n6_2.5[,3,] , 1 , mean ),
                                        p95_sd =  apply( rmse_array_n6_2.5[,3,] , 1 , sd ))
      
      rmse_array_n3_2.5 <- array( dim=c(6,3,5) )
      
      for (i in 1:5) { rmse_array_n3_2.5[,,i] <- as.matrix(rmse.result( simulation_summary2[[i]] , true_p95 , true_gsd , true_gm )[,2:4]) }
      
      rmse_table_n3_2.5 <- data.frame( approach = c("ideal_F" , "ideal_B" ,
                                                     "naive_F" , "naive_B",
                                                     "me_F" , "me_B"),
                                        gm_mean =  apply( rmse_array_n3_2.5[,1,] , 1 , mean ),
                                        gm_sd =  apply( rmse_array_n3_2.5[,1,] , 1 , sd ),
                                        gsd_mean =  apply( rmse_array_n3_2.5[,2,] , 1 , mean ),
                                        gsd_sd =  apply( rmse_array_n3_2.5[,2,] , 1 , sd ),
                                        p95_mean =  apply( rmse_array_n3_2.5[,3,] , 1 , mean ),
                                        p95_sd =  apply( rmse_array_n3_2.5[,3,] , 1 , sd ))
      
      rmse_array_n6_1.5 <- array( dim=c(6,3,5) )
      
      for (i in 1:5) { rmse_array_n6_1.5[,,i] <- as.matrix(rmse.result( simulation_summary3[[i]] , true_p95 , true_gsd , true_gm )[,2:4]) }
      
      rmse_table_n6_1.5 <- data.frame( approach = c("ideal_F" , "ideal_B" ,
                                                     "naive_F" , "naive_B",
                                                     "me_F" , "me_B"),
                                        gm_mean =  apply( rmse_array_n6_1.5[,1,] , 1 , mean ),
                                        gm_sd =  apply( rmse_array_n6_1.5[,1,] , 1 , sd ),
                                        gsd_mean =  apply( rmse_array_n6_1.5[,2,] , 1 , mean ),
                                        gsd_sd =  apply( rmse_array_n6_1.5[,2,] , 1 , sd ),
                                        p95_mean =  apply( rmse_array_n6_1.5[,3,] , 1 , mean ),
                                        p95_sd =  apply( rmse_array_n6_1.5[,3,] , 1 , sd ))
      
      rmse_array_n3_1.5 <- array( dim=c(6,3,5) )
      
      for (i in 1:5) { rmse_array_n3_1.5[,,i] <- as.matrix(rmse.result( simulation_summary4[[i]] , true_p95 , true_gsd , true_gm )[,2:4]) }
      
      rmse_table_n3_1.5 <- data.frame( approach = c("ideal_F" , "ideal_B" ,
                                                     "naive_F" , "naive_B",
                                                     "me_F" , "me_B"),
                                        gm_mean =  apply( rmse_array_n3_1.5[,1,] , 1 , mean ),
                                        gm_sd =  apply( rmse_array_n3_1.5[,1,] , 1 , sd ),
                                        gsd_mean =  apply( rmse_array_n3_1.5[,2,] , 1 , mean ),
                                        gsd_sd =  apply( rmse_array_n3_1.5[,2,] , 1 , sd ),
                                        p95_mean =  apply( rmse_array_n3_1.5[,3,] , 1 , mean ),
                                        p95_sd =  apply( rmse_array_n3_1.5[,3,] , 1 , sd ))
      
###### export ######
      
  result <- list( coverage = coverage_table,
                  bias = list( n6_2.5 = bias_table_n6_2.5 , n3_2.5 = bias_table_n3_2.5 , n6_1.5 = bias_table_n6_1.5 , n3_1.5 = bias_table_n3_1.5 ),
                  precision = list( n6_2.5 = precision_table_n6_2.5 , n3_2.5 = precision_table_n3_2.5 , n6_1.5 = precision_table_n6_1.5 , n3_1.5 = precision_table_n3_1.5 ),
                  rmse = list( n6_2.5 = rmse_table_n6_2.5 , n3_2.5 = rmse_table_n3_2.5 , n6_1.5 = rmse_table_n6_1.5 , n3_1.5 = rmse_table_n3_1.5 ) )

      # Make sure the rstudioapi package is installed
      if (!require(rstudioapi)) {
        install.packages("rstudioapi")
      }
      
      # Get the path to the active project
      project_directory <- rstudioapi::getActiveProject()
      
      # Set the working directory to the project directory
      setwd(project_directory)      
      
            
      saveRDS(result , "created data/GUM measurement error.RDS")
      
      