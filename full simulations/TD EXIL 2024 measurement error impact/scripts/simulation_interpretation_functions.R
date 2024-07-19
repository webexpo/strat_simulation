#'  script for interpreting the simulation results for one scenario


#### FUNCTIONS ####

#' function computes rmse for GM, GSD, P95 and Exceedance from the results of the simulation for one scenario  
#'
#' @param results_one_scenario index of the simulation 
#' @param true_gm list containing the simulated data for a single scenario 
#' @param true_gsd value for the measurement error as standard error CV, not in % and not expanded uncertainty CV
#' @param true_p95 occupational exposure limit
#' @param true_exceedance_perc number of iteration for the GUM approach
#'
#' @return data.frame of results with one line per approach, and columns for GM, GSD, P95 and Exceedance
#' @return column "perc_estimable" contains the percentage of iterations with estimable values
#' @return result is zero if percentage of estimable values is zero (e.g. all NDS and ROS approach)


rmse.result <- function( results_one_scenario , true_gm , true_gsd , true_p95 , true_exceedance_perc ) { 
  
  
  rmse_table <- data.frame( method = c("ideal_b" , "ideal_f" ,
                                       "naive_b" , "naive_f",
                                       "me_b" , "me_f_mean" , "me_f_median" , "me_f_quantile"),
                            gm = numeric(8),
                            gsd = numeric(8),
                            p95 = numeric(8))
  
  rmse_table$gm <- apply(results_one_scenario$array[1,,], 1, compute.rmse , theta = true_gm ) 
  
  rmse_table$gsd <- apply(results_one_scenario$array[2,,], 1, compute.rmse , theta = true_gsd )
  
  rmse_table$p95 <- apply(results_one_scenario$array[3,,], 1, compute.rmse , theta = true_p95 )
  
  rmse_table$exceedance <- apply(results_one_scenario$array[6,,], 1, compute.rmse , theta = true_exceedance_perc )
  
  rmse_table$perc_estimable <- apply(results_one_scenario$array[1,,], 1, function(x) { 100*sum(!is.na(x))/length(x) } )
  
  return(rmse_table)
  
}


#' function computes precision for GM, GSD, P95 and Exceedance from the results of the simulation for one scenario  
#'
#' @param results_one_scenario index of the simulation 
#' @param true_gm list containing the simulated data for a single scenario 
#' @param true_gsd value for the measurement error as standard error CV, not in % and not expanded uncertainty CV
#' @param true_p95 occupational exposure limit
#' @param true_exceedance_perc number of iteration for the GUM approach
#'
#' @return data.frame of results with one line per approach, and columns for GM, GSD, P95 and Exceedance
#' @return column "perc_estimable" contains the percentage of iterations with estimable values
#' @return result is zero if percentage of estimable values is zero (e.g. all NDS and ROS approach)

precision.result <- function( results_one_scenario , true_p95 , true_gsd , true_gm , true_exceedance_perc ) { 
  
  
  precision_table <- data.frame( method = c("ideal_b" , "ideal_f" ,
                                            "naive_b" , "naive_f",
                                            "me_b" , "me_f_mean" , "me_f_median" , "me_f_quantile"),
                                 gm = numeric(8),
                                 gsd = numeric(8),
                                 p95 = numeric(8))
  
  precision_table$gm <- apply(results_one_scenario$array[1,,], 1, compute.precision , theta = true_gm ) 
  
  precision_table$gsd <- apply(results_one_scenario$array[2,,], 1, compute.precision , theta = true_gsd )
  
  precision_table$p95 <- apply(results_one_scenario$array[3,,], 1, compute.precision , theta = true_p95 )
  
  precision_table$exceedance <- apply(results_one_scenario$array[6,,], 1, compute.precision , theta = true_exceedance_perc )
  
  precision_table$perc_estimable <- apply(results_one_scenario$array[1,,], 1, function(x) { 100*sum(!is.na(x))/length(x) } )
  
  
  return(precision_table)
  
}


#' function computes bias for GM, GSD, P95 and Exceedance from the results of the simulation for one scenario  
#'
#' @param results_one_scenario index of the simulation 
#' @param true_gm list containing the simulated data for a single scenario 
#' @param true_gsd value for the measurement error as standard error CV, not in % and not expanded uncertainty CV
#' @param true_p95 occupational exposure limit
#' @param true_exceedance_perc number of iteration for the GUM approach
#'
#' @return data.frame of results with one line per approach, and columns for GM, GSD, P95 and Exceedance
#' @return column "perc_estimable" contains the percentage of iterations with estimable values
#' @return result is zero if percentage of estimable values is zero (e.g. all NDS and ROS approach)


bias.result <- function( results_one_scenario , true_p95 , true_gsd , true_gm , true_exceedance_perc ) { 
  
  
  bias_table <- data.frame( method = c("ideal_b" , "ideal_f" ,
                                       "naive_b" , "naive_f",
                                       "me_b" , "me_f_mean" , "me_f_median" , "me_f_quantile"),
                            gm = numeric(8),
                            gsd = numeric(8),
                            p95 = numeric(8))
  
  bias_table$gm <- apply(results_one_scenario$array[1,,], 1, compute.bias , theta = true_gm ) 
  
  bias_table$gsd <- apply(results_one_scenario$array[2,,], 1, compute.bias , theta = true_gsd )
  
  bias_table$p95 <- apply(results_one_scenario$array[3,,], 1, compute.bias , theta = true_p95 )
  
  bias_table$exceedance <- apply(results_one_scenario$array[6,,], 1, compute.bias , theta = true_exceedance_perc )
  
  bias_table$perc_estimable <- apply(results_one_scenario$array[1,,], 1, function(x) { 100*sum(!is.na(x))/length(x) } )
  
  
  return(bias_table)
  
}

#' function computes median error for GM, GSD, P95 and Exceedance from the results of the simulation for one scenario  
#'
#' @param results_one_scenario index of the simulation 
#' @param true_gm list containing the simulated data for a single scenario 
#' @param true_gsd value for the measurement error as standard error CV, not in % and not expanded uncertainty CV
#' @param true_p95 occupational exposure limit
#' @param true_exceedance_perc number of iteration for the GUM approach
#'
#' @return data.frame of results with one line per approach, and columns for GM, GSD, P95 and Exceedance
#' @return column "perc_estimable" contains the percentage of iterations with estimable values
#' @return result is zero if percentage of estimable values is zero (e.g. all NDS and ROS approach)


median.error.result <- function( results_one_scenario , true_p95 , true_gsd , true_gm , true_exceedance_perc ) { 
  
  
  median_error_table <- data.frame( method = c("ideal_b" , "ideal_f" ,
                                               "naive_b" , "naive_f",
                                               "me_b" , "me_f_mean" , "me_f_median" , "me_f_quantile"),
                                    gm = numeric(8),
                                    gsd = numeric(8),
                                    p95 = numeric(8))
  
  median_error_table$gm <- apply(results_one_scenario$array[1,,], 1, compute.median.error , theta = true_gm ) 
  
  median_error_table$gsd <- apply(results_one_scenario$array[2,,], 1, compute.median.error , theta = true_gsd )
  
  median_error_table$p95 <- apply(results_one_scenario$array[3,,], 1, compute.median.error , theta = true_p95 )
  
  median_error_table$exceedance <- apply(results_one_scenario$array[6,,], 1, compute.median.error , theta = true_exceedance_perc )
  
  median_error_table$perc_estimable <- apply(results_one_scenario$array[1,,], 1, function(x) { 100*sum(!is.na(x))/length(x) } )
  
  
  return(median_error_table)
  
}

#' function computes rmsle for GM, GSD, P95 and Exceedance from the results of the simulation for one scenario  
#'
#' @param results_one_scenario index of the simulation 
#' @param true_gm list containing the simulated data for a single scenario 
#' @param true_gsd value for the measurement error as standard error CV, not in % and not expanded uncertainty CV
#' @param true_p95 occupational exposure limit
#' @param true_exceedance_perc number of iteration for the GUM approach
#'
#' @return data.frame of results with one line per approach, and columns for GM, GSD, P95 and Exceedance
#' @return column "perc_estimable" contains the percentage of iterations with estimable values
#' @return result is zero if percentage of estimable values is zero (e.g. all NDS and ROS approach)


rmsle.result <- function( results_one_scenario , true_gm , true_gsd , true_p95 , true_exceedance_perc ) { 
  
  
  rmsle_table <- data.frame( method = c("ideal_b" , "ideal_f" ,
                                        "naive_b" , "naive_f",
                                        "me_b" , "me_f_mean" , "me_f_median" , "me_f_quantile"),
                             gm = numeric(8),
                             gsd = numeric(8),
                             p95 = numeric(8))
  
  rmsle_table$gm <- apply(results_one_scenario$array[1,,], 1, compute.rmsle , theta = true_gm ) 
  
  rmsle_table$gsd <- apply(results_one_scenario$array[2,,], 1, compute.rmsle , theta = true_gsd )
  
  rmsle_table$p95 <- apply(results_one_scenario$array[3,,], 1, compute.rmsle , theta = true_p95 )
  
  rmsle_table$exceedance <- apply(results_one_scenario$array[6,,], 1, compute.rmsle , theta = true_exceedance_perc )
  
  rmsle_table$perc_estimable <- apply(results_one_scenario$array[1,,], 1, function(x) { 100*sum(!is.na(x))/length(x) } )
  
  
  return(rmsle_table)
  
}



#' function computes median absolute deviation  for GM, GSD, P95 and Exceedance from the results of the simulation for one scenario  
#'
#' @param results_one_scenario index of the simulation 
#' @param true_gm list containing the simulated data for a single scenario 
#' @param true_gsd value for the measurement error as standard error CV, not in % and not expanded uncertainty CV
#' @param true_p95 occupational exposure limit
#' @param true_exceedance_perc number of iteration for the GUM approach
#'
#' @return data.frame of results with one line per approach, and columns for GM, GSD, P95 and Exceedance
#' @return column "perc_estimable" contains the percentage of iterations with estimable values
#' @return result is zero if percentage of estimable values is zero (e.g. all NDS and ROS approach)


mad.result <- function( results_one_scenario , true_gm , true_gsd , true_p95 , true_exceedance_perc ) { 
  
  
  mad_table <- data.frame( method = c("ideal_b" , "ideal_f" ,
                                      "naive_b" , "naive_f",
                                      "me_b" , "me_f_mean" , "me_f_median" , "me_f_quantile"),
                           gm = numeric(8),
                           gsd = numeric(8),
                           p95 = numeric(8))
  
  mad_table$gm <- apply(results_one_scenario$array[1,,], 1, compute.mad , theta = true_gm ) 
  
  mad_table$gsd <- apply(results_one_scenario$array[2,,], 1, compute.mad , theta = true_gsd )
  
  mad_table$p95 <- apply(results_one_scenario$array[3,,], 1, compute.mad , theta = true_p95 )
  
  mad_table$exceedance <- apply(results_one_scenario$array[6,,], 1, compute.mad , theta = true_exceedance_perc )
  
  mad_table$perc_estimable <- apply(results_one_scenario$array[1,,], 1, function(x) { 100*sum(!is.na(x))/length(x) } )
  
  return(mad_table)
  
}


#' function computes coverage for the 70 and 95% UCLs for P95 and exceedance from the results of the simulation for one scenario  
#'
#' @param results_one_scenario index of the simulation 
#' @param true_p95 occupational exposure limit
#' @param true_exceedance_perc number of iteration for the GUM approach
#'
#' @return data.frame of results with one line per approach, and columns for 70% and 95% UCLs for P95 and Exceedance
#' @return column "perc_estimable" contains the percentage of iterations with estimable values
#' @return result is zero if percentage of estimable values is zero (e.g. all NDS and ROS approach)


coverage.result <- function( results_one_scenario , true_p95 , true_exceedance_perc ) { 
  
  #results_one_scenario <- test_parallel
  
  #x <- results_one_scenario$array[4,6,]
  
  coverage_table <- data.frame( method = c("ideal_b" , "ideal_f" ,
                                           "naive_b" , "naive_f",
                                           "me_b" , "me_f_mean" , "me_f_median" , "me_f_quantile"),
                                p95_ucl70 = numeric(8),
                                p95_ucl95 = numeric(8),
                                exceedance_ucl70 = numeric(8),
                                exceedance_ucl95 = numeric(8))
  
  coverage_table$p95_ucl70 <- apply(results_one_scenario$array[4,,], 1, function(x) {
    
    if ( sum(is.na(x)) == length(x)) return(0) else {   x <- x[!is.na(x)]
                                                        return(100*sum(x>=true_p95)/length(x))} } )
  
  coverage_table$p95_ucl95 <- apply(results_one_scenario$array[5,,], 1, function(x) {
    
    if ( sum(is.na(x)) == length(x)) return(0) else {   x <- x[!is.na(x)]
                                                        return(100*sum(x>=true_p95)/length(x))}} )
  
  coverage_table$exceedance_ucl70 <- apply(results_one_scenario$array[7,,], 1, function(x) {
    
    if ( sum(is.na(x)) == length(x)) return(0) else {   x <- x[!is.na(x)]
                                                        return(100*sum(x>=true_exceedance_perc)/length(x))}} ) 
  
  coverage_table$exceedance_ucl95 <- apply(results_one_scenario$array[8,,], 1, function(x) {
    
    if ( sum(is.na(x)) == length(x)) return(0) else {   x <- x[!is.na(x)]
                                                        return(100*sum(x>=true_exceedance_perc)/length(x))}} ) 
  
  coverage_table$perc_estimable <- apply(results_one_scenario$array[1,,], 1, function(x) { 100*sum(!is.na(x))/length(x) } )
  
  return(coverage_table)
  
}

#' function computes percentage of  errors for the four strategies UTL95,70 /  UTL95,95 /  F_UCL70 and F_UCL95 (as F and P95 are estimated differently)  
#' for truly OK situations, error = 100-specificity, for truly not OK situations, error = 100-sensitivity
#'
#' @param results_one_scenario index of the simulation 
#' @param true_p95 occupational exposure limit
#' @param true_exceedance_perc number of iteration for the GUM approach
#' @param oel occupational exposure limit
#'
#' @return data.frame of results with one line per approach, and columns for UTL95,70 /  UTL95,95 /  F_UCL70 and F_UCL95
#' @return column "perc_estimable" contains the percentage of iterations with estimable values
#' @return result is zero if percentage of estimable values is zero (e.g. all NDS and ROS approach)


perc.mistake.result <- function( results_one_scenario , true_p95 , true_exceedance_perc , oel ) { 
  
  perc_mistake_table <- data.frame( method = c("ideal_b" , "ideal_f" ,
                                               "naive_b" , "naive_f",
                                               "me_b" , "me_f_mean" , "me_f_median" , "me_f_quantile"),
                                    p95_ucl70 = numeric(8),
                                    p95_ucl95 = numeric(8),
                                    exceedance_ucl70 = numeric(8),
                                    exceedance_ucl95 = numeric(8))
  
  
  perc_mistake_table$p95_ucl70 <- apply(results_one_scenario$array[4,,], 1, function(x) { 
    if ( sum(is.na(x)) == length(x)) return(0) else {
      x <- x[!is.na(x)]
      if (true_p95>=oel) results <- 100*sum(x<oel)/length(x)
      else results <- 100*sum(x>=oel)/length(x) 
      return(results)}
  }) 
  
  perc_mistake_table$p95_ucl95 <- apply(results_one_scenario$array[5,,], 1, function(x) { 
    if ( sum(is.na(x)) == length(x)) return(0) else {
      x <- x[!is.na(x)]
      if (true_p95>=oel) results <- 100*sum(x<oel)/length(x)
      else results <- 100*sum(x>=oel)/length(x) 
      return(results) } 
  })
  
  perc_mistake_table$exceedance_ucl70 <- apply(results_one_scenario$array[7,,], 1, function(x) { 
    if ( sum(is.na(x)) == length(x)) return(0) else {
      x <- x[!is.na(x)]
      if (true_p95>=oel) results <- 100*sum(x<5)/length(x)
      else results <- 100*sum(x>=5)/length(x) 
      return(results) }
  })
  
  perc_mistake_table$exceedance_ucl95 <- apply(results_one_scenario$array[8,,], 1, function(x) {
    if ( sum(is.na(x)) == length(x)) return(0) else {
      x <- x[!is.na(x)]
      if (true_p95>=oel) results <- 100*sum(x<5)/length(x)
      else results <- 100*sum(x>=5)/length(x) 
      return(results) }
  }) 
  
  perc_mistake_table$perc_estimable <- apply(results_one_scenario$array[1,,], 1, function(x) { 100*sum(!is.na(x))/length(x) } )
  
  return(perc_mistake_table)
  
}



