#'  script for interpreting the simulation results for one scenario


#### FUNCTIONS ####

#' function computes rmse for GM, GSD, P95 and Exceedance from the results of the simulation for one scenario  
#'
#' @param results_one_scenario list containing the simulated data for a single scenario 
#' @param true_gm true GM, NA in the "real gsd" approach 
#' @param true_gsd true GSD, NA in the "real gsd" approach 
#' @param true_p95 true 95th percentile
#' @param true_exceedance_perc true exceedance fraction in percentage
#'
#' @return data.frame of results with one line per approach, and columns for GM, GSD, P95 and Exceedance
#' @return column "perc_estimable" contains the percentage of iterations with estimable values
#' @return result is zero if percentage of estimable values is zero (e.g. all NDS and ROS approach)


rmse.result <- function( results_one_scenario , true_gm , true_gsd , true_p95 , true_exceedance_perc ) { 
  
  
  rmse_table <- data.frame( method = c("ideal_b","ideal_f","naive_b","naive_f","me_b","me_f_mean","me_f_median","me_f_q2.5","me_f_q5","me_f_q95","me_f_q97.5"),
                            gm = numeric(11),
                            gsd = numeric(11),
                            p95 = numeric(11))
  
  if (is.na(true_gm[1])) rmse_table$gm <- NA else rmse_table$gm <- apply(results_one_scenario$array[1,,], 1, compute.rmse , theta = true_gm ) 
  
  if (is.na(true_gsd[1])) rmse_table$gsd <- NA else rmse_table$gsd <- apply(results_one_scenario$array[2,,], 1, compute.rmse , theta = true_gsd )
  
  rmse_table$p95 <- apply(results_one_scenario$array[3,,], 1, compute.rmse , theta = true_p95 )
  
  rmse_table$exceedance <- apply(results_one_scenario$array[6,,], 1, compute.rmse , theta = true_exceedance_perc )
  
  rmse_table$perc_estimable <- apply(results_one_scenario$array[1,,], 1, function(x) { 100*sum(!is.na(x))/length(x) } )
  
  return(rmse_table)
  
}


#' function computes precision for GM, GSD, P95 and Exceedance from the results of the simulation for one scenario  
#'
#' @param results_one_scenario list containing the simulated data for a single scenario 
#' @param true_gm true GM, NA in the "real gsd" approach 
#' @param true_gsd true GSD, NA in the "real gsd" approach 
#' @param true_p95 true 95th percentile
#' @param true_exceedance_perc true exceedance fraction in percentage
#'
#' @return data.frame of results with one line per approach, and columns for GM, GSD, P95 and Exceedance
#' @return column "perc_estimable" contains the percentage of iterations with estimable values
#' @return result is zero if percentage of estimable values is zero (e.g. all NDS and ROS approach)

precision.result <- function( results_one_scenario , true_p95 , true_gsd , true_gm , true_exceedance_perc ) { 
  
  
  precision_table <- data.frame( method = c("ideal_b","ideal_f","naive_b","naive_f","me_b","me_f_mean","me_f_median","me_f_q2.5","me_f_q5","me_f_q95","me_f_q97.5"),
                                 gm = numeric(11),
                                 gsd = numeric(11),
                                 p95 = numeric(11))
  
  if (is.na(true_gm[1])) precision_table$gm <- NA else precision_table$gm <- apply(results_one_scenario$array[1,,], 1, compute.precision , theta = true_gm ) 
  
  if (is.na(true_gsd[1])) precision_table$gsd <- NA else precision_table$gsd <- apply(results_one_scenario$array[2,,], 1, compute.precision , theta = true_gsd )
  
  precision_table$p95 <- apply(results_one_scenario$array[3,,], 1, compute.precision , theta = true_p95 )
  
  precision_table$exceedance <- apply(results_one_scenario$array[6,,], 1, compute.precision , theta = true_exceedance_perc )
  
  precision_table$perc_estimable <- apply(results_one_scenario$array[1,,], 1, function(x) { 100*sum(!is.na(x))/length(x) } )
  
  
  return(precision_table)
  
}


#' function computes bias for GM, GSD, P95 and Exceedance from the results of the simulation for one scenario  
#'
#' @param results_one_scenario list containing the simulated data for a single scenario 
#' @param true_gm true GM, NA in the "real gsd" approach 
#' @param true_gsd true GSD, NA in the "real gsd" approach 
#' @param true_p95 true 95th percentile
#' @param true_exceedance_perc true exceedance fraction in percentage
#'
#' @return data.frame of results with one line per approach, and columns for GM, GSD, P95 and Exceedance
#' @return column "perc_estimable" contains the percentage of iterations with estimable values
#' @return result is zero if percentage of estimable values is zero (e.g. all NDS and ROS approach)


bias.result <- function( results_one_scenario , true_p95 , true_gsd , true_gm , true_exceedance_perc ) { 
  
  
  bias_table <- data.frame( method = c("ideal_b","ideal_f","naive_b","naive_f","me_b","me_f_mean","me_f_median","me_f_q2.5","me_f_q5","me_f_q95","me_f_q97.5"),
                            gm = numeric(11),
                            gsd = numeric(11),
                            p95 = numeric(11))
  
  if (is.na(true_gm[1])) bias_table$gm <- NA else bias_table$gm <- apply(results_one_scenario$array[1,,], 1, compute.bias , theta = true_gm ) 
  
  if (is.na(true_gsd[1])) bias_table$gsd <- NA else bias_table$gsd <- apply(results_one_scenario$array[2,,], 1, compute.bias , theta = true_gsd )
  
  bias_table$p95 <- apply(results_one_scenario$array[3,,], 1, compute.bias , theta = true_p95 )
  
  bias_table$exceedance <- apply(results_one_scenario$array[6,,], 1, compute.bias , theta = true_exceedance_perc )
  
  bias_table$perc_estimable <- apply(results_one_scenario$array[1,,], 1, function(x) { 100*sum(!is.na(x))/length(x) } )
  
  
  return(bias_table)
  
}

#' function computes median error for GM, GSD, P95 and Exceedance from the results of the simulation for one scenario  
#'
#' @param results_one_scenario list containing the simulated data for a single scenario 
#' @param true_gm true GM, NA in the "real gsd" approach 
#' @param true_gsd true GSD, NA in the "real gsd" approach 
#' @param true_p95 true 95th percentile
#' @param true_exceedance_perc true exceedance fraction in percentage
#'
#' @return data.frame of results with one line per approach, and columns for GM, GSD, P95 and Exceedance
#' @return column "perc_estimable" contains the percentage of iterations with estimable values
#' @return result is zero if percentage of estimable values is zero (e.g. all NDS and ROS approach)


median.error.result <- function( results_one_scenario , true_p95 , true_gsd , true_gm , true_exceedance_perc ) { 
  
  
  median_error_table <- data.frame( method = c("ideal_b","ideal_f","naive_b","naive_f","me_b","me_f_mean","me_f_median","me_f_q2.5","me_f_q5","me_f_q95","me_f_q97.5"),
                                    gm = numeric(11),
                                    gsd = numeric(11),
                                    p95 = numeric(11))
  
  if (is.na(true_gm[1])) median_error_table$gm <- NA else median_error_table$gm <- apply(results_one_scenario$array[1,,], 1, compute.median.error , theta = true_gm ) 
  
  if (is.na(true_gsd[1])) median_error_table$gsd <- NA else median_error_table$gsd <- apply(results_one_scenario$array[2,,], 1, compute.median.error , theta = true_gsd )
  
  median_error_table$p95 <- apply(results_one_scenario$array[3,,], 1, compute.median.error , theta = true_p95 )
  
  median_error_table$exceedance <- apply(results_one_scenario$array[6,,], 1, compute.median.error , theta = true_exceedance_perc )
  
  median_error_table$perc_estimable <- apply(results_one_scenario$array[1,,], 1, function(x) { 100*sum(!is.na(x))/length(x) } )
  
  
  return(median_error_table)
  
}

#' function computes rmsle for GM, GSD, P95 and Exceedance from the results of the simulation for one scenario  
#'
#' @param results_one_scenario list containing the simulated data for a single scenario 
#' @param true_gm true GM, NA in the "real gsd" approach 
#' @param true_gsd true GSD, NA in the "real gsd" approach 
#' @param true_p95 true 95th percentile
#' @param true_exceedance_perc true exceedance fraction in percentage
#'
#' @return data.frame of results with one line per approach, and columns for GM, GSD, P95 and Exceedance
#' @return column "perc_estimable" contains the percentage of iterations with estimable values
#' @return result is zero if percentage of estimable values is zero (e.g. all NDS and ROS approach)


rmsle.result <- function( results_one_scenario , true_gm , true_gsd , true_p95 , true_exceedance_perc ) { 
  
  
  rmsle_table <- data.frame( method = c("ideal_b","ideal_f","naive_b","naive_f","me_b","me_f_mean","me_f_median","me_f_q2.5","me_f_q5","me_f_q95","me_f_q97.5"),
                             gm = numeric(11),
                             gsd = numeric(11),
                             p95 = numeric(11))
  
  if (is.na(true_gm[1])) rmsle_table$gm <- NA else rmsle_table$gm <- apply(results_one_scenario$array[1,,], 1, compute.rmsle , theta = true_gm ) 
  
  if (is.na(true_gsd[1])) rmsle_table$gsd <- NA else rmsle_table$gsd <- apply(results_one_scenario$array[2,,], 1, compute.rmsle , theta = true_gsd )
  
  rmsle_table$p95 <- apply(results_one_scenario$array[3,,], 1, compute.rmsle , theta = true_p95 )
  
  rmsle_table$exceedance <- apply(results_one_scenario$array[6,,], 1, compute.rmsle , theta = true_exceedance_perc )
  
  rmsle_table$perc_estimable <- apply(results_one_scenario$array[1,,], 1, function(x) { 100*sum(!is.na(x))/length(x) } )
  
  
  return(rmsle_table)
  
}



#' function computes median absolute deviation  for GM, GSD, P95 and Exceedance from the results of the simulation for one scenario  
#'
#' @param results_one_scenario list containing the simulated data for a single scenario 
#' @param true_gm true GM, NA in the "real gsd" approach 
#' @param true_gsd true GSD, NA in the "real gsd" approach 
#' @param true_p95 true 95th percentile
#' @param true_exceedance_perc true exceedance fraction in percentage
#'
#' @return data.frame of results with one line per approach, and columns for GM, GSD, P95 and Exceedance
#' @return column "perc_estimable" contains the percentage of iterations with estimable values
#' @return result is zero if percentage of estimable values is zero (e.g. all NDS and ROS approach)


mad.result <- function( results_one_scenario , true_gm , true_gsd , true_p95 , true_exceedance_perc ) { 
  
  
  mad_table <- data.frame( method = c("ideal_b","ideal_f","naive_b","naive_f","me_b","me_f_mean","me_f_median","me_f_q2.5","me_f_q5","me_f_q95","me_f_q97.5"),
                           gm = numeric(11),
                           gsd = numeric(11),
                           p95 = numeric(11))
  
  if (is.na(true_gm[1])) mad_table$gm <- NA else mad_table$gm <- apply(results_one_scenario$array[1,,], 1, compute.mad , theta = true_gm ) 
  
  if (is.na(true_gsd[1])) mad_table$gsd <- NA else mad_table$gsd <- apply(results_one_scenario$array[2,,], 1, compute.mad , theta = true_gsd )
  
  mad_table$p95 <- apply(results_one_scenario$array[3,,], 1, compute.mad , theta = true_p95 )
  
  mad_table$exceedance <- apply(results_one_scenario$array[6,,], 1, compute.mad , theta = true_exceedance_perc )
  
  mad_table$perc_estimable <- apply(results_one_scenario$array[1,,], 1, function(x) { 100*sum(!is.na(x))/length(x) } )
  
  return(mad_table)
  
}


#' function computes coverage for the 70 and 95% UCLs for P95 and exceedance from the results of the simulation for one scenario  
#'
#' @param results_one_scenario index of the simulation 
#' @param true_p95 true 95th percentile
#' @param true_exceedance_perc true exceedance fraction in percentage
#'
#' @return data.frame of results with one line per approach, and columns for 70% and 95% UCLs for P95 and Exceedance
#' @return column "perc_estimable" contains the percentage of iterations with estimable values
#' @return result is zero if percentage of estimable values is zero (e.g. all NDS and ROS approach)


coverage.result <- function( results_one_scenario , true_p95 , true_exceedance_perc ) { 
  
  coverage_table <- data.frame( method = c("ideal_b","ideal_f","naive_b","naive_f","me_b","me_f_mean","me_f_median","me_f_q2.5","me_f_q5","me_f_q95","me_f_q97.5"),
                                p95_ucl70 = numeric(11),
                                p95_ucl95 = numeric(11),
                                exceedance_ucl70 = numeric(11),
                                exceedance_ucl95 = numeric(11))
  
  

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
#' @param true_p95 true 95th percentile
#' @param oel vector of occupational exposure limit across iterations
#'
#' @return data.frame of results with one line per approach, and columns for UTL95,70 /  UTL95,95 /  F_UCL70 and F_UCL95
#' @return column "perc_estimable" contains the percentage of iterations with estimable values
#' @return result is zero if percentage of estimable values is zero (e.g. all NDS and ROS approach)


perc.mistake.result <- function( results_one_scenario , true_p95 , true_exceedance_perc , oel ) { 
  
  #results_one_scenario <- test_parallel
  
  
  perc_mistake_table <- data.frame( method = c("ideal_b","ideal_f","naive_b","naive_f","me_b","me_f_mean","me_f_median","me_f_q2.5","me_f_q5","me_f_q95","me_f_q97.5"),
                                    p95_ucl70 = numeric(11),
                                    p95_ucl95 = numeric(11),
                                    exceedance_ucl70 = numeric(11),
                                    exceedance_ucl95 = numeric(11))
  
  # matrix of OELs, with the same values for each methods, but different across iteration (real gsd approach), for comparison with the UCLs
  
  oel_matrix <- matrix(nrow = dim(results_one_scenario$array[4,,])[1], ncol = dim(results_one_scenario$array[4,,])[2])
  
  for (i in 1:dim(results_one_scenario$array[4,,])[1]) oel_matrix[i,] <- oel
  
  # dichotomous matrix of the UCLs, with TRUE if UCL is above OEL, FALSE otherwise
  
  p95_ucl70_array_over_oel <- results_one_scenario$array[4,,] >= oel_matrix
  p95_ucl95_array_over_oel <- results_one_scenario$array[5,,] >= oel_matrix
  
  # calculation of the mistake rate

  perc_mistake_table$p95_ucl70 <- apply(p95_ucl70_array_over_oel, 1, function(x) { 
    if ( sum(is.na(x)) == length(x)) return(0) else {
      x <- x[!is.na(x)]
      if (true_exceedance_perc>=5) results <- 100*sum(!x)/length(x)
      else results <- 100*sum(x)/length(x) 
      return(results)}
  }) 
  
  perc_mistake_table$p95_ucl95 <- apply(p95_ucl95_array_over_oel, 1, function(x) {
    if ( sum(is.na(x)) == length(x)) return(0) else {
      x <- x[!is.na(x)]
      if (true_exceedance_perc>=5) results <- 100*sum(!x)/length(x)
      else results <- 100*sum(x)/length(x) 
      return(results)}
  })
    

  
  perc_mistake_table$exceedance_ucl70 <- apply(results_one_scenario$array[7,,], 1, function(x) { 
    if ( sum(is.na(x)) == length(x)) return(0) else {
      x <- x[!is.na(x)]
      if (true_exceedance_perc>=5) results <- 100*sum(x<5)/length(x)
      else results <- 100*sum(x>=5)/length(x) 
      return(results) }
  })
  
  perc_mistake_table$exceedance_ucl95 <- apply(results_one_scenario$array[8,,], 1, function(x) {
    if ( sum(is.na(x)) == length(x)) return(0) else {
      x <- x[!is.na(x)]
      if (true_exceedance_perc>=5) results <- 100*sum(x<5)/length(x)
      else results <- 100*sum(x>=5)/length(x) 
      return(results) }
  }) 
  
  perc_mistake_table$perc_estimable <- apply(results_one_scenario$array[1,,], 1, function(x) { 100*sum(!is.na(x))/length(x) } )
  
  return(perc_mistake_table)
  
}


#' function which creates the interpretation of the simulations for a single scenario and one run 
#'
#' @param index index of the scenario of interest
#' @param simulation_result object of the results of a full run of the simulation 
#'
#' @return list of results with all performance metrics calculated across the iterations
#'


one_run_summary <- function( index , simulation_result ) {
  
  
  run_inter <- list( rmse = rmse.result( results_one_scenario = simulation_result[[index]] , 
                                         true_gm = scenarios$true_gm[index] ,
                                         true_gsd = scenarios$true_gsd[index] ,
                                         true_p95 = scenarios$true_p95[index] ,
                                         true_exceedance_perc = scenarios$true_exceedance_perc[index] ),
                     
                     precision = precision.result( results_one_scenario = simulation_result[[index]] , 
                                                   true_gm = scenarios$true_gm[index] , 
                                                   true_gsd = scenarios$true_gsd[index] ,
                                                   true_p95 = scenarios$true_p95[index] ,
                                                   true_exceedance_perc = scenarios$true_exceedance_perc[index] ),
                     
                     bias = bias.result( results_one_scenario = simulation_result[[index]] , 
                                         true_gm = scenarios$true_gm[index] , 
                                         true_gsd = scenarios$true_gsd[index] ,
                                         true_p95 = scenarios$true_p95[index] ,
                                         true_exceedance_perc = scenarios$true_exceedance_perc[index] ),
                     
                     median_error = median.error.result( results_one_scenario = simulation_result[[index]] , 
                                                         true_gm = scenarios$true_gm[index] , 
                                                         true_gsd = scenarios$true_gsd[index] ,
                                                         true_p95 = scenarios$true_p95[index] , 
                                                         true_exceedance_perc = scenarios$true_exceedance_perc[index] ),
                     
                     rmsle = rmsle.result( results_one_scenario = simulation_result[[index]] , 
                                           true_gm = scenarios$true_gm[index] , 
                                           true_gsd = scenarios$true_gsd[index] , 
                                           true_p95 = scenarios$true_p95[index] , 
                                           true_exceedance_perc = scenarios$true_exceedance_perc[index] ),
                     
                     mad = mad.result( results_one_scenario = simulation_result[[index]] , 
                                       true_gm = scenarios$true_gm[index] , 
                                       true_gsd = scenarios$true_gsd[index] , 
                                       true_p95 = scenarios$true_p95[index] , 
                                       true_exceedance_perc = scenarios$true_exceedance_perc[index] ),
                     
                     coverage = coverage.result( results_one_scenario = simulation_result[[index]] , 
                                                 true_p95 = scenarios$true_p95[index] , 
                                                 true_exceedance_perc = scenarios$true_exceedance_perc[index] ),
                     
                     perc_mistake = perc.mistake.result( results_one_scenario = simulation_result[[index]] , 
                                                         true_p95 = scenarios$true_p95[index] , 
                                                         true_exceedance_perc = scenarios$true_exceedance_perc[index] ,
                                                         oel = scenarios$oel[index] )
  )
  
  return(run_inter)
  
  
  
}



#' function which creates the interpretation of the simulations for a single scenario and three runs 
#'
#' @param run1 first run, object output of the on_run_summary
#' @param run2 second run, object output of the on_run_summary
#' @param run3 three run, object output of the on_run_summary
#'  
#'
#' @return list of results stratfied by performance metric then exposure metric, then : 3 runs, mean and CV
#'


three_run_summary <- function( run1 , run2 , run3 ) {
  
  results <- list()
  
  # rmse ------------------
  
  results$rmse <-  list()
  
  
  #gm
  
  results$rmse$gm <- data.frame( method = run1$rmse$method ,
                                 rep1 = run1$rmse$gm , 
                                 rep2 = run2$rmse$gm , 
                                 rep3 = run3$rmse$gm  )
  
  
  # ading the mean
  
  results$rmse$gm$mean <- rowMeans( results$rmse$gm[ , 2:4 ] )
  
  # adding the sd
  
  results$rmse$gm$cv_perc <- apply( results$rmse$gm[ , 2:4 ] , 1 , sd )*100 / results$rmse$gm$mean
  
  #gsd
  
  results$rmse$gsd <- data.frame( method = run1$rmse$method ,
                                  rep1 = run1$rmse$gsd , 
                                  rep2 = run2$rmse$gsd , 
                                  rep3 = run3$rmse$gsd  )
  
  
  # ading the mean
  results$rmse$gsd$mean <- rowMeans( results$rmse$gsd[ , 2:4 ] )
  
  # adding the sd
  results$rmse$gsd$cv <- apply( results$rmse$gsd[ , 2:4 ] , 1 , sd )*100 / results$rmse$gsd$mean   
  
  #p95
  
  results$rmse$p95 <- data.frame( method = run1$rmse$method ,
                                  rep1 = run1$rmse$p95 , 
                                  rep2 = run2$rmse$p95 , 
                                  rep3 = run3$rmse$p95  )
  
  
  # ading the mean
  results$rmse$p95$mean <- rowMeans( results$rmse$p95[ , 2:4 ] )
  
  # adding the sd
  results$rmse$p95$cv <- apply( results$rmse$p95[ , 2:4 ] , 1 , sd )*100 / results$rmse$p95$mean        
  
  # exceedance
  
  results$rmse$exceedance <- data.frame( method = run1$rmse$method ,
                                         rep1 = run1$rmse$exceedance , 
                                         rep2 = run2$rmse$exceedance , 
                                         rep3 = run3$rmse$exceedance  )
  
  
  # ading the mean
  results$rmse$exceedance$mean <- rowMeans( results$rmse$exceedance[ , 2:4 ] )
  
  # adding the sd
  results$rmse$exceedance$cv <- apply( results$rmse$exceedance[ , 2:4 ] , 1 , sd )*100 / results$rmse$exceedance$mean    
  
  # perc_estimable
  
  results$rmse$perc_estimable <- data.frame( method = run1$rmse$method ,
                                             rep1 = run1$rmse$perc_estimable , 
                                             rep2 = run2$rmse$perc_estimable , 
                                             rep3 = run3$rmse$perc_estimable  )
  
  
  # ading the mean
  results$rmse$perc_estimable$mean <- rowMeans( results$rmse$perc_estimable[ , 2:4 ] )
  
  # adding the sd
  results$rmse$perc_estimable$cv <- apply( results$rmse$perc_estimable[ , 2:4 ] , 1 , sd )*100 / results$rmse$perc_estimable$mean        
  
  # precision  ------------------
  
  results$precision <- list()
  
  
  #gm
  
  results$precision$gm <- data.frame( method = run1$precision$method ,
                                      rep1 = run1$precision$gm , 
                                      rep2 = run2$precision$gm , 
                                      rep3 = run3$precision$gm  )
  
  
  # ading the mean
  
  results$precision$gm$mean <- rowMeans( results$precision$gm[ , 2:4 ] )
  
  # adding the sd
  
  results$precision$gm$cv_perc <- apply( results$precision$gm[ , 2:4 ] , 1 , sd )*100 / results$precision$gm$mean
  
  #gsd
  
  results$precision$gsd <- data.frame( method = run1$precision$method ,
                                       rep1 = run1$precision$gsd , 
                                       rep2 = run2$precision$gsd , 
                                       rep3 = run3$precision$gsd  )
  
  
  # ading the mean
  results$precision$gsd$mean <- rowMeans( results$precision$gsd[ , 2:4 ] )
  
  # adding the sd
  results$precision$gsd$cv <- apply( results$precision$gsd[ , 2:4 ] , 1 , sd )*100 / results$precision$gsd$mean   
  
  #p95
  
  results$precision$p95 <- data.frame( method = run1$precision$method ,
                                       rep1 = run1$precision$p95 , 
                                       rep2 = run2$precision$p95 , 
                                       rep3 = run3$precision$p95  )
  
  
  # ading the mean
  results$precision$p95$mean <- rowMeans( results$precision$p95[ , 2:4 ] )
  
  # adding the sd
  results$precision$p95$cv <- apply( results$precision$p95[ , 2:4 ] , 1 , sd )*100 / results$precision$p95$mean        
  
  # exceedance
  
  results$precision$exceedance <- data.frame( method = run1$precision$method ,
                                              rep1 = run1$precision$exceedance , 
                                              rep2 = run2$precision$exceedance , 
                                              rep3 = run3$precision$exceedance  )
  
  
  # ading the mean
  results$precision$exceedance$mean <- rowMeans( results$precision$exceedance[ , 2:4 ] )
  
  # adding the sd
  results$precision$exceedance$cv <- apply( results$precision$exceedance[ , 2:4 ] , 1 , sd )*100 / results$precision$exceedance$mean    
  
  # perc_estimable
  
  results$precision$perc_estimable <- data.frame( method = run1$precision$method ,
                                                  rep1 = run1$precision$perc_estimable , 
                                                  rep2 = run2$precision$perc_estimable , 
                                                  rep3 = run3$precision$perc_estimable  )
  
  
  # ading the mean
  results$precision$perc_estimable$mean <- rowMeans( results$precision$perc_estimable[ , 2:4 ] )
  
  # adding the sd
  results$precision$perc_estimable$cv <- apply( results$precision$perc_estimable[ , 2:4 ] , 1 , sd )*100 / results$precision$perc_estimable$mean                
  
  
  
  # bias  ------------------
  
  results$bias <- list()
  
  
  #gm
  
  results$bias$gm <- data.frame( method = run1$bias$method ,
                                 rep1 = run1$bias$gm , 
                                 rep2 = run2$bias$gm , 
                                 rep3 = run3$bias$gm  )
  
  
  # ading the mean
  
  results$bias$gm$mean <- rowMeans( results$bias$gm[ , 2:4 ] )
  
  # adding the sd
  
  results$bias$gm$cv_perc <- apply( results$bias$gm[ , 2:4 ] , 1 , sd )*100 / results$bias$gm$mean
  
  #gsd
  
  results$bias$gsd <- data.frame( method = run1$bias$method ,
                                  rep1 = run1$bias$gsd , 
                                  rep2 = run2$bias$gsd , 
                                  rep3 = run3$bias$gsd  )
  
  
  # ading the mean
  results$bias$gsd$mean <- rowMeans( results$bias$gsd[ , 2:4 ] )
  
  # adding the sd
  results$bias$gsd$cv <- apply( results$bias$gsd[ , 2:4 ] , 1 , sd )*100 / results$bias$gsd$mean   
  
  #p95
  
  results$bias$p95 <- data.frame( method = run1$bias$method ,
                                  rep1 = run1$bias$p95 , 
                                  rep2 = run2$bias$p95 , 
                                  rep3 = run3$bias$p95  )
  
  
  # ading the mean
  results$bias$p95$mean <- rowMeans( results$bias$p95[ , 2:4 ] )
  
  # adding the sd
  results$bias$p95$cv <- apply( results$bias$p95[ , 2:4 ] , 1 , sd )*100 / results$bias$p95$mean        
  
  # exceedance
  
  results$bias$exceedance <- data.frame( method = run1$bias$method ,
                                         rep1 = run1$bias$exceedance , 
                                         rep2 = run2$bias$exceedance , 
                                         rep3 = run3$bias$exceedance  )
  
  
  # ading the mean
  results$bias$exceedance$mean <- rowMeans( results$bias$exceedance[ , 2:4 ] )
  
  # adding the sd
  results$bias$exceedance$cv <- apply( results$bias$exceedance[ , 2:4 ] , 1 , sd )*100 / results$bias$exceedance$mean    
  
  # perc_estimable
  
  results$bias$perc_estimable <- data.frame( method = run1$bias$method ,
                                             rep1 = run1$bias$perc_estimable , 
                                             rep2 = run2$bias$perc_estimable , 
                                             rep3 = run3$bias$perc_estimable  )
  
  
  # ading the mean
  results$bias$perc_estimable$mean <- rowMeans( results$bias$perc_estimable[ , 2:4 ] )
  
  # adding the sd
  results$bias$perc_estimable$cv <- apply( results$bias$perc_estimable[ , 2:4 ] , 1 , sd )*100 / results$bias$perc_estimable$mean                
  
  
  # median_error  ------------------
  
  results$median_error <- list()
  
  
  #gm
  
  results$median_error$gm <- data.frame( method = run1$median_error$method ,
                                         rep1 = run1$median_error$gm , 
                                         rep2 = run2$median_error$gm , 
                                         rep3 = run3$median_error$gm  )
  
  
  # ading the mean
  
  results$median_error$gm$mean <- rowMeans( results$median_error$gm[ , 2:4 ] )
  
  # adding the sd
  
  results$median_error$gm$cv_perc <- apply( results$median_error$gm[ , 2:4 ] , 1 , sd )*100 / results$median_error$gm$mean
  
  #gsd
  
  results$median_error$gsd <- data.frame( method = run1$median_error$method ,
                                          rep1 = run1$median_error$gsd , 
                                          rep2 = run2$median_error$gsd , 
                                          rep3 = run3$median_error$gsd  )
  
  
  # ading the mean
  results$median_error$gsd$mean <- rowMeans( results$median_error$gsd[ , 2:4 ] )
  
  # adding the sd
  results$median_error$gsd$cv <- apply( results$median_error$gsd[ , 2:4 ] , 1 , sd )*100 / results$median_error$gsd$mean   
  
  #p95
  
  results$median_error$p95 <- data.frame( method = run1$median_error$method ,
                                          rep1 = run1$median_error$p95 , 
                                          rep2 = run2$median_error$p95 , 
                                          rep3 = run3$median_error$p95  )
  
  
  # ading the mean
  results$median_error$p95$mean <- rowMeans( results$median_error$p95[ , 2:4 ] )
  
  # adding the sd
  results$median_error$p95$cv <- apply( results$median_error$p95[ , 2:4 ] , 1 , sd )*100 / results$median_error$p95$mean        
  
  # exceedance
  
  results$median_error$exceedance <- data.frame( method = run1$median_error$method ,
                                                 rep1 = run1$median_error$exceedance , 
                                                 rep2 = run2$median_error$exceedance , 
                                                 rep3 = run3$median_error$exceedance  )
  
  
  # ading the mean
  results$median_error$exceedance$mean <- rowMeans( results$median_error$exceedance[ , 2:4 ] )
  
  # adding the sd
  results$median_error$exceedance$cv <- apply( results$median_error$exceedance[ , 2:4 ] , 1 , sd )*100 / results$median_error$exceedance$mean    
  
  # perc_estimable
  
  results$median_error$perc_estimable <- data.frame( method = run1$median_error$method ,
                                                     rep1 = run1$median_error$perc_estimable , 
                                                     rep2 = run2$median_error$perc_estimable , 
                                                     rep3 = run3$median_error$perc_estimable  )
  
  
  # ading the mean
  results$median_error$perc_estimable$mean <- rowMeans( results$median_error$perc_estimable[ , 2:4 ] )
  
  # adding the sd
  results$median_error$perc_estimable$cv <- apply( results$median_error$perc_estimable[ , 2:4 ] , 1 , sd )*100 / results$median_error$perc_estimable$mean                
  
  
  # rmsle  ------------------
  
  results$rmsle <- list()
  
  
  #gm
  
  results$rmsle$gm <- data.frame( method = run1$rmsle$method ,
                                  rep1 = run1$rmsle$gm , 
                                  rep2 = run2$rmsle$gm , 
                                  rep3 = run3$rmsle$gm  )
  
  
  # ading the mean
  
  results$rmsle$gm$mean <- rowMeans( results$rmsle$gm[ , 2:4 ] )
  
  # adding the sd
  
  results$rmsle$gm$cv_perc <- apply( results$rmsle$gm[ , 2:4 ] , 1 , sd )*100 / results$rmsle$gm$mean
  
  #gsd
  
  results$rmsle$gsd <- data.frame( method = run1$rmsle$method ,
                                   rep1 = run1$rmsle$gsd , 
                                   rep2 = run2$rmsle$gsd , 
                                   rep3 = run3$rmsle$gsd  )
  
  
  # ading the mean
  results$rmsle$gsd$mean <- rowMeans( results$rmsle$gsd[ , 2:4 ] )
  
  # adding the sd
  results$rmsle$gsd$cv <- apply( results$rmsle$gsd[ , 2:4 ] , 1 , sd )*100 / results$rmsle$gsd$mean   
  
  #p95
  
  results$rmsle$p95 <- data.frame( method = run1$rmsle$method ,
                                   rep1 = run1$rmsle$p95 , 
                                   rep2 = run2$rmsle$p95 , 
                                   rep3 = run3$rmsle$p95  )
  
  
  # ading the mean
  results$rmsle$p95$mean <- rowMeans( results$rmsle$p95[ , 2:4 ] )
  
  # adding the sd
  results$rmsle$p95$cv <- apply( results$rmsle$p95[ , 2:4 ] , 1 , sd )*100 / results$rmsle$p95$mean        
  
  # exceedance
  
  results$rmsle$exceedance <- data.frame( method = run1$rmsle$method ,
                                          rep1 = run1$rmsle$exceedance , 
                                          rep2 = run2$rmsle$exceedance , 
                                          rep3 = run3$rmsle$exceedance  )
  
  
  # ading the mean
  results$rmsle$exceedance$mean <- rowMeans( results$rmsle$exceedance[ , 2:4 ] )
  
  # adding the sd
  results$rmsle$exceedance$cv <- apply( results$rmsle$exceedance[ , 2:4 ] , 1 , sd )*100 / results$rmsle$exceedance$mean    
  
  # perc_estimable
  
  results$rmsle$perc_estimable <- data.frame( method = run1$rmsle$method ,
                                              rep1 = run1$rmsle$perc_estimable , 
                                              rep2 = run2$rmsle$perc_estimable , 
                                              rep3 = run3$rmsle$perc_estimable  )
  
  
  # ading the mean
  results$rmsle$perc_estimable$mean <- rowMeans( results$rmsle$perc_estimable[ , 2:4 ] )
  
  # adding the sd
  results$rmsle$perc_estimable$cv <- apply( results$rmsle$perc_estimable[ , 2:4 ] , 1 , sd )*100 / results$rmsle$perc_estimable$mean                
  
  
  
  # mad  ------------------
  
  results$mad <- list()
  
  
  #gm
  
  results$mad$gm <- data.frame( method = run1$mad$method ,
                                rep1 = run1$mad$gm , 
                                rep2 = run2$mad$gm , 
                                rep3 = run3$mad$gm  )
  
  
  # ading the mean
  
  results$mad$gm$mean <- rowMeans( results$mad$gm[ , 2:4 ] )
  
  # adding the sd
  
  results$mad$gm$cv_perc <- apply( results$mad$gm[ , 2:4 ] , 1 , sd )*100 / results$mad$gm$mean
  
  #gsd
  
  results$mad$gsd <- data.frame( method = run1$mad$method ,
                                 rep1 = run1$mad$gsd , 
                                 rep2 = run2$mad$gsd , 
                                 rep3 = run3$mad$gsd  )
  
  
  # ading the mean
  results$mad$gsd$mean <- rowMeans( results$mad$gsd[ , 2:4 ] )
  
  # adding the sd
  results$mad$gsd$cv <- apply( results$mad$gsd[ , 2:4 ] , 1 , sd )*100 / results$mad$gsd$mean   
  
  #p95
  
  results$mad$p95 <- data.frame( method = run1$mad$method ,
                                 rep1 = run1$mad$p95 , 
                                 rep2 = run2$mad$p95 , 
                                 rep3 = run3$mad$p95  )
  
  
  # ading the mean
  results$mad$p95$mean <- rowMeans( results$mad$p95[ , 2:4 ] )
  
  # adding the sd
  results$mad$p95$cv <- apply( results$mad$p95[ , 2:4 ] , 1 , sd )*100 / results$mad$p95$mean        
  
  # exceedance
  
  results$mad$exceedance <- data.frame( method = run1$mad$method ,
                                        rep1 = run1$mad$exceedance , 
                                        rep2 = run2$mad$exceedance , 
                                        rep3 = run3$mad$exceedance  )
  
  
  # ading the mean
  results$mad$exceedance$mean <- rowMeans( results$mad$exceedance[ , 2:4 ] )
  
  # adding the sd
  results$mad$exceedance$cv <- apply( results$mad$exceedance[ , 2:4 ] , 1 , sd )*100 / results$mad$exceedance$mean    
  
  # perc_estimable
  
  results$mad$perc_estimable <- data.frame( method = run1$mad$method ,
                                            rep1 = run1$mad$perc_estimable , 
                                            rep2 = run2$mad$perc_estimable , 
                                            rep3 = run3$mad$perc_estimable  )
  
  
  # ading the mean
  results$mad$perc_estimable$mean <- rowMeans( results$mad$perc_estimable[ , 2:4 ] )
  
  # adding the sd
  results$mad$perc_estimable$cv <- apply( results$mad$perc_estimable[ , 2:4 ] , 1 , sd )*100 / results$mad$perc_estimable$mean                
  
  # coverage  ------------------
  
  results$coverage <-list()
  
  
  #p95_ucl70
  
  results$coverage$p95_ucl70 <- data.frame( method = run1$coverage$method ,
                                            rep1 = run1$coverage$p95_ucl70 , 
                                            rep2 = run2$coverage$p95_ucl70 , 
                                            rep3 = run3$coverage$p95_ucl70  )
  
  
  # ading the mean
  
  results$coverage$p95_ucl70$mean <- rowMeans( results$coverage$p95_ucl70[ , 2:4 ] )
  
  # adding the sd
  
  results$coverage$p95_ucl70$cv_perc <- apply( results$coverage$p95_ucl70[ , 2:4 ] , 1 , sd )*100 / results$coverage$p95_ucl70$mean
  
  #p95_ucl95
  
  results$coverage$p95_ucl95 <- data.frame( method = run1$coverage$method ,
                                            rep1 = run1$coverage$p95_ucl95 , 
                                            rep2 = run2$coverage$p95_ucl95 , 
                                            rep3 = run3$coverage$p95_ucl95  )
  
  
  # ading the mean
  results$coverage$p95_ucl95$mean <- rowMeans( results$coverage$p95_ucl95[ , 2:4 ] )
  
  # adding the sd
  results$coverage$p95_ucl95$cv <- apply( results$coverage$p95_ucl95[ , 2:4 ] , 1 , sd )*100 / results$coverage$p95_ucl95$mean   
  
  #exceedance_ucl70
  
  results$coverage$exceedance_ucl70 <- data.frame( method = run1$coverage$method ,
                                                   rep1 = run1$coverage$exceedance_ucl70 , 
                                                   rep2 = run2$coverage$exceedance_ucl70 , 
                                                   rep3 = run3$coverage$exceedance_ucl70  )
  
  
  # ading the mean
  results$coverage$exceedance_ucl70$mean <- rowMeans( results$coverage$exceedance_ucl70[ , 2:4 ] )
  
  # adding the sd
  results$coverage$exceedance_ucl70$cv <- apply( results$coverage$exceedance_ucl70[ , 2:4 ] , 1 , sd )*100 / results$coverage$exceedance_ucl70$mean        
  
  # exceedance_ucl95
  
  results$coverage$exceedance_ucl95 <- data.frame( method = run1$coverage$method ,
                                                   rep1 = run1$coverage$exceedance_ucl95 , 
                                                   rep2 = run2$coverage$exceedance_ucl95 , 
                                                   rep3 = run3$coverage$exceedance_ucl95  )
  
  
  # ading the mean
  results$coverage$exceedance_ucl95$mean <- rowMeans( results$coverage$exceedance_ucl95[ , 2:4 ] )
  
  # adding the sd
  results$coverage$exceedance_ucl95$cv <- apply( results$coverage$exceedance_ucl95[ , 2:4 ] , 1 , sd )*100 / results$coverage$exceedance_ucl95$mean    
  
  # perc_estimable
  
  results$coverage$perc_estimable <- data.frame( method = run1$coverage$method ,
                                                 rep1 = run1$coverage$perc_estimable , 
                                                 rep2 = run2$coverage$perc_estimable , 
                                                 rep3 = run3$coverage$perc_estimable  )
  
  
  # ading the mean
  results$coverage$perc_estimable$mean <- rowMeans( results$coverage$perc_estimable[ , 2:4 ] )
  
  # adding the sd
  results$coverage$perc_estimable$cv <- apply( results$coverage$perc_estimable[ , 2:4 ] , 1 , sd )*100 / results$coverage$perc_estimable$mean                
  
  
  
  
  
  
  
  
  # perc_mistake  ------------------
  
  results$perc_mistake <- list()
  
  
  #p95_ucl70
  
  results$perc_mistake$p95_ucl70 <- data.frame( method = run1$perc_mistake$method ,
                                                rep1 = run1$perc_mistake$p95_ucl70 , 
                                                rep2 = run2$perc_mistake$p95_ucl70 , 
                                                rep3 = run3$perc_mistake$p95_ucl70  )
  
  
  # ading the mean
  
  results$perc_mistake$p95_ucl70$mean <- rowMeans( results$perc_mistake$p95_ucl70[ , 2:4 ] )
  
  # adding the sd
  
  results$perc_mistake$p95_ucl70$cv_perc <- apply( results$perc_mistake$p95_ucl70[ , 2:4 ] , 1 , sd )*100 / results$perc_mistake$p95_ucl70$mean
  
  #p95_ucl95
  
  results$perc_mistake$p95_ucl95 <- data.frame( method = run1$perc_mistake$method ,
                                                rep1 = run1$perc_mistake$p95_ucl95 , 
                                                rep2 = run2$perc_mistake$p95_ucl95 , 
                                                rep3 = run3$perc_mistake$p95_ucl95  )
  
  
  # ading the mean
  results$perc_mistake$p95_ucl95$mean <- rowMeans( results$perc_mistake$p95_ucl95[ , 2:4 ] )
  
  # adding the sd
  results$perc_mistake$p95_ucl95$cv <- apply( results$perc_mistake$p95_ucl95[ , 2:4 ] , 1 , sd )*100 / results$perc_mistake$p95_ucl95$mean   
  
  #exceedance_ucl70
  
  results$perc_mistake$exceedance_ucl70 <- data.frame( method = run1$perc_mistake$method ,
                                                       rep1 = run1$perc_mistake$exceedance_ucl70 , 
                                                       rep2 = run2$perc_mistake$exceedance_ucl70 , 
                                                       rep3 = run3$perc_mistake$exceedance_ucl70  )
  
  
  # ading the mean
  results$perc_mistake$exceedance_ucl70$mean <- rowMeans( results$perc_mistake$exceedance_ucl70[ , 2:4 ] )
  
  # adding the sd
  results$perc_mistake$exceedance_ucl70$cv <- apply( results$perc_mistake$exceedance_ucl70[ , 2:4 ] , 1 , sd )*100 / results$perc_mistake$exceedance_ucl70$mean        
  
  # exceedance_ucl95
  
  results$perc_mistake$exceedance_ucl95 <- data.frame( method = run1$perc_mistake$method ,
                                                       rep1 = run1$perc_mistake$exceedance_ucl95 , 
                                                       rep2 = run2$perc_mistake$exceedance_ucl95 , 
                                                       rep3 = run3$perc_mistake$exceedance_ucl95  )
  
  
  # ading the mean
  results$perc_mistake$exceedance_ucl95$mean <- rowMeans( results$perc_mistake$exceedance_ucl95[ , 2:4 ] )
  
  # adding the sd
  results$perc_mistake$exceedance_ucl95$cv <- apply( results$perc_mistake$exceedance_ucl95[ , 2:4 ] , 1 , sd )*100 / results$perc_mistake$exceedance_ucl95$mean    
  
  # perc_estimable
  
  results$perc_mistake$perc_estimable <- data.frame( method = run1$perc_mistake$method ,
                                                     rep1 = run1$perc_mistake$perc_estimable , 
                                                     rep2 = run2$perc_mistake$perc_estimable , 
                                                     rep3 = run3$perc_mistake$perc_estimable  )
  
  
  # ading the mean
  results$perc_mistake$perc_estimable$mean <- rowMeans( results$perc_mistake$perc_estimable[ , 2:4 ] )
  
  # adding the sd
  results$perc_mistake$perc_estimable$cv <- apply( results$perc_mistake$perc_estimable[ , 2:4 ] , 1 , sd )*100 / results$perc_mistake$perc_estimable$mean                
  
  
  
  
  
  
  
  # finalization --------------
  
  
  
  return(results)
  
  
  
}

#### FUNCTIONS - REAL GSD ####


#' function which creates the interpretation of the simulations for a single scenario and one run - REAL GSD
#'
#' @param index index of the scenario of interest
#' @param simulation_result object of the results of a full run of the simulation 
#'
#' @return list of results with all performance metrics calculated across the iterations
#'


one_run_summary_gsd <- function( index , simulation_result , simulation_data ) {
  
  #index=i , simulation_result=run3 , simulation_data = run3_data
  
  
  run_inter <- list( rmse = rmse.result( results_one_scenario = simulation_result[[index]] , 
                                         true_gm = NA ,
                                         true_gsd = NA ,
                                         true_p95 = scenarios$true_p95[index] ,
                                         true_exceedance_perc = scenarios$true_exceedance_perc[index] ),
                     
                     precision = precision.result( results_one_scenario = simulation_result[[index]] , 
                                                   true_gm = NA , 
                                                   true_gsd = NA ,
                                                   true_p95 = scenarios$true_p95[index] ,
                                                   true_exceedance_perc = scenarios$true_exceedance_perc[index] ),
                     
                     bias = bias.result( results_one_scenario = simulation_result[[index]] , 
                                         true_gm = NA , 
                                         true_gsd = NA ,
                                         true_p95 = scenarios$true_p95[index] ,
                                         true_exceedance_perc = scenarios$true_exceedance_perc[index] ),
                     
                     median_error = median.error.result( results_one_scenario = simulation_result[[index]] , 
                                                         true_gm = NA , 
                                                         true_gsd = NA ,
                                                         true_p95 = scenarios$true_p95[index] , 
                                                         true_exceedance_perc = scenarios$true_exceedance_perc[index] ),
                     
                     rmsle = rmsle.result( results_one_scenario = simulation_result[[index]] , 
                                           true_gm = NA , 
                                           true_gsd = NA , 
                                           true_p95 = scenarios$true_p95[index] , 
                                           true_exceedance_perc = scenarios$true_exceedance_perc[index] ),
                     
                     mad = mad.result( results_one_scenario = simulation_result[[index]] , 
                                       true_gm = NA , 
                                       true_gsd = NA , 
                                       true_p95 = scenarios$true_p95[index] , 
                                       true_exceedance_perc = scenarios$true_exceedance_perc[index] ),
                     
                     coverage = coverage.result( results_one_scenario = simulation_result[[index]] , 
                                                 true_p95 = scenarios$true_p95[index] , 
                                                 true_exceedance_perc = scenarios$true_exceedance_perc[index] ),
                     
                     perc_mistake = perc.mistake.result( results_one_scenario = simulation_result[[index]] , 
                                                         true_p95 = scenarios$true_p95[index] , 
                                                         true_exceedance_perc = scenarios$true_exceedance_perc[index] ,
                                                         oel = simulation_data[[index]]$oel )
  )
  
  return(run_inter)
  
  
  
}




#' function which creates the interpretation of the simulations for a single scenario and three runs - real GSD
#'
#' @param run1 first run, object output of the on_run_summary
#' @param run2 second run, object output of the on_run_summary
#' @param run3 three run, object output of the on_run_summary
#'  
#'
#' @return list of results stratfied by performance metric then exposure metric, then : 3 runs, mean and CV
#'


three_run_summary_gsd <- function( run1 , run2 , run3 ) {
  
  results <- list()
  
  # rmse ------------------
  
  results$rmse <-  list()
  
  #p95
  
  results$rmse$p95 <- data.frame( method = run1$rmse$method ,
                                  rep1 = run1$rmse$p95 , 
                                  rep2 = run2$rmse$p95 , 
                                  rep3 = run3$rmse$p95  )
  
  
  # ading the mean
  results$rmse$p95$mean <- rowMeans( results$rmse$p95[ , 2:4 ] )
  
  # adding the sd
  results$rmse$p95$cv <- apply( results$rmse$p95[ , 2:4 ] , 1 , sd )*100 / results$rmse$p95$mean        
  
  # exceedance
  
  results$rmse$exceedance <- data.frame( method = run1$rmse$method ,
                                         rep1 = run1$rmse$exceedance , 
                                         rep2 = run2$rmse$exceedance , 
                                         rep3 = run3$rmse$exceedance  )
  
  
  # ading the mean
  results$rmse$exceedance$mean <- rowMeans( results$rmse$exceedance[ , 2:4 ] )
  
  # adding the sd
  results$rmse$exceedance$cv <- apply( results$rmse$exceedance[ , 2:4 ] , 1 , sd )*100 / results$rmse$exceedance$mean    
  
  # perc_estimable
  
  results$rmse$perc_estimable <- data.frame( method = run1$rmse$method ,
                                             rep1 = run1$rmse$perc_estimable , 
                                             rep2 = run2$rmse$perc_estimable , 
                                             rep3 = run3$rmse$perc_estimable  )
  
  
  # ading the mean
  results$rmse$perc_estimable$mean <- rowMeans( results$rmse$perc_estimable[ , 2:4 ] )
  
  # adding the sd
  results$rmse$perc_estimable$cv <- apply( results$rmse$perc_estimable[ , 2:4 ] , 1 , sd )*100 / results$rmse$perc_estimable$mean        
  
  # precision  ------------------
  
  results$precision <- list()
  #p95
  
  results$precision$p95 <- data.frame( method = run1$precision$method ,
                                       rep1 = run1$precision$p95 , 
                                       rep2 = run2$precision$p95 , 
                                       rep3 = run3$precision$p95  )
  
  
  # ading the mean
  results$precision$p95$mean <- rowMeans( results$precision$p95[ , 2:4 ] )
  
  # adding the sd
  results$precision$p95$cv <- apply( results$precision$p95[ , 2:4 ] , 1 , sd )*100 / results$precision$p95$mean        
  
  # exceedance
  
  results$precision$exceedance <- data.frame( method = run1$precision$method ,
                                              rep1 = run1$precision$exceedance , 
                                              rep2 = run2$precision$exceedance , 
                                              rep3 = run3$precision$exceedance  )
  
  
  # ading the mean
  results$precision$exceedance$mean <- rowMeans( results$precision$exceedance[ , 2:4 ] )
  
  # adding the sd
  results$precision$exceedance$cv <- apply( results$precision$exceedance[ , 2:4 ] , 1 , sd )*100 / results$precision$exceedance$mean    
  
  # perc_estimable
  
  results$precision$perc_estimable <- data.frame( method = run1$precision$method ,
                                                  rep1 = run1$precision$perc_estimable , 
                                                  rep2 = run2$precision$perc_estimable , 
                                                  rep3 = run3$precision$perc_estimable  )
  
  
  # ading the mean
  results$precision$perc_estimable$mean <- rowMeans( results$precision$perc_estimable[ , 2:4 ] )
  
  # adding the sd
  results$precision$perc_estimable$cv <- apply( results$precision$perc_estimable[ , 2:4 ] , 1 , sd )*100 / results$precision$perc_estimable$mean                
  
  
  
  # bias  ------------------
  
  results$bias <- list()
  
  
  #p95
  
  results$bias$p95 <- data.frame( method = run1$bias$method ,
                                  rep1 = run1$bias$p95 , 
                                  rep2 = run2$bias$p95 , 
                                  rep3 = run3$bias$p95  )
  
  
  # ading the mean
  results$bias$p95$mean <- rowMeans( results$bias$p95[ , 2:4 ] )
  
  # adding the sd
  results$bias$p95$cv <- apply( results$bias$p95[ , 2:4 ] , 1 , sd )*100 / results$bias$p95$mean        
  
  # exceedance
  
  results$bias$exceedance <- data.frame( method = run1$bias$method ,
                                         rep1 = run1$bias$exceedance , 
                                         rep2 = run2$bias$exceedance , 
                                         rep3 = run3$bias$exceedance  )
  
  
  # ading the mean
  results$bias$exceedance$mean <- rowMeans( results$bias$exceedance[ , 2:4 ] )
  
  # adding the sd
  results$bias$exceedance$cv <- apply( results$bias$exceedance[ , 2:4 ] , 1 , sd )*100 / results$bias$exceedance$mean    
  
  # perc_estimable
  
  results$bias$perc_estimable <- data.frame( method = run1$bias$method ,
                                             rep1 = run1$bias$perc_estimable , 
                                             rep2 = run2$bias$perc_estimable , 
                                             rep3 = run3$bias$perc_estimable  )
  
  
  # ading the mean
  results$bias$perc_estimable$mean <- rowMeans( results$bias$perc_estimable[ , 2:4 ] )
  
  # adding the sd
  results$bias$perc_estimable$cv <- apply( results$bias$perc_estimable[ , 2:4 ] , 1 , sd )*100 / results$bias$perc_estimable$mean                
  
  
  # median_error  ------------------
  
  results$median_error <- list()
  
  
  #p95
  
  results$median_error$p95 <- data.frame( method = run1$median_error$method ,
                                          rep1 = run1$median_error$p95 , 
                                          rep2 = run2$median_error$p95 , 
                                          rep3 = run3$median_error$p95  )
  
  
  # ading the mean
  results$median_error$p95$mean <- rowMeans( results$median_error$p95[ , 2:4 ] )
  
  # adding the sd
  results$median_error$p95$cv <- apply( results$median_error$p95[ , 2:4 ] , 1 , sd )*100 / results$median_error$p95$mean        
  
  # exceedance
  
  results$median_error$exceedance <- data.frame( method = run1$median_error$method ,
                                                 rep1 = run1$median_error$exceedance , 
                                                 rep2 = run2$median_error$exceedance , 
                                                 rep3 = run3$median_error$exceedance  )
  
  
  # ading the mean
  results$median_error$exceedance$mean <- rowMeans( results$median_error$exceedance[ , 2:4 ] )
  
  # adding the sd
  results$median_error$exceedance$cv <- apply( results$median_error$exceedance[ , 2:4 ] , 1 , sd )*100 / results$median_error$exceedance$mean    
  
  # perc_estimable
  
  results$median_error$perc_estimable <- data.frame( method = run1$median_error$method ,
                                                     rep1 = run1$median_error$perc_estimable , 
                                                     rep2 = run2$median_error$perc_estimable , 
                                                     rep3 = run3$median_error$perc_estimable  )
  
  
  # ading the mean
  results$median_error$perc_estimable$mean <- rowMeans( results$median_error$perc_estimable[ , 2:4 ] )
  
  # adding the sd
  results$median_error$perc_estimable$cv <- apply( results$median_error$perc_estimable[ , 2:4 ] , 1 , sd )*100 / results$median_error$perc_estimable$mean                
  
  
  # rmsle  ------------------
  
  results$rmsle <- list()
  
  
  #p95
  
  results$rmsle$p95 <- data.frame( method = run1$rmsle$method ,
                                   rep1 = run1$rmsle$p95 , 
                                   rep2 = run2$rmsle$p95 , 
                                   rep3 = run3$rmsle$p95  )
  
  
  # ading the mean
  results$rmsle$p95$mean <- rowMeans( results$rmsle$p95[ , 2:4 ] )
  
  # adding the sd
  results$rmsle$p95$cv <- apply( results$rmsle$p95[ , 2:4 ] , 1 , sd )*100 / results$rmsle$p95$mean        
  
  # exceedance
  
  results$rmsle$exceedance <- data.frame( method = run1$rmsle$method ,
                                          rep1 = run1$rmsle$exceedance , 
                                          rep2 = run2$rmsle$exceedance , 
                                          rep3 = run3$rmsle$exceedance  )
  
  
  # ading the mean
  results$rmsle$exceedance$mean <- rowMeans( results$rmsle$exceedance[ , 2:4 ] )
  
  # adding the sd
  results$rmsle$exceedance$cv <- apply( results$rmsle$exceedance[ , 2:4 ] , 1 , sd )*100 / results$rmsle$exceedance$mean    
  
  # perc_estimable
  
  results$rmsle$perc_estimable <- data.frame( method = run1$rmsle$method ,
                                              rep1 = run1$rmsle$perc_estimable , 
                                              rep2 = run2$rmsle$perc_estimable , 
                                              rep3 = run3$rmsle$perc_estimable  )
  
  
  # ading the mean
  results$rmsle$perc_estimable$mean <- rowMeans( results$rmsle$perc_estimable[ , 2:4 ] )
  
  # adding the sd
  results$rmsle$perc_estimable$cv <- apply( results$rmsle$perc_estimable[ , 2:4 ] , 1 , sd )*100 / results$rmsle$perc_estimable$mean                
  
  
  
  # mad  ------------------
  
  results$mad <- list()
  
  
  #p95
  
  results$mad$p95 <- data.frame( method = run1$mad$method ,
                                 rep1 = run1$mad$p95 , 
                                 rep2 = run2$mad$p95 , 
                                 rep3 = run3$mad$p95  )
  
  
  # ading the mean
  results$mad$p95$mean <- rowMeans( results$mad$p95[ , 2:4 ] )
  
  # adding the sd
  results$mad$p95$cv <- apply( results$mad$p95[ , 2:4 ] , 1 , sd )*100 / results$mad$p95$mean        
  
  # exceedance
  
  results$mad$exceedance <- data.frame( method = run1$mad$method ,
                                        rep1 = run1$mad$exceedance , 
                                        rep2 = run2$mad$exceedance , 
                                        rep3 = run3$mad$exceedance  )
  
  
  # ading the mean
  results$mad$exceedance$mean <- rowMeans( results$mad$exceedance[ , 2:4 ] )
  
  # adding the sd
  results$mad$exceedance$cv <- apply( results$mad$exceedance[ , 2:4 ] , 1 , sd )*100 / results$mad$exceedance$mean    
  
  # perc_estimable
  
  results$mad$perc_estimable <- data.frame( method = run1$mad$method ,
                                            rep1 = run1$mad$perc_estimable , 
                                            rep2 = run2$mad$perc_estimable , 
                                            rep3 = run3$mad$perc_estimable  )
  
  
  # ading the mean
  results$mad$perc_estimable$mean <- rowMeans( results$mad$perc_estimable[ , 2:4 ] )
  
  # adding the sd
  results$mad$perc_estimable$cv <- apply( results$mad$perc_estimable[ , 2:4 ] , 1 , sd )*100 / results$mad$perc_estimable$mean                
  
  # coverage  ------------------
  
  results$coverage <-list()
  
  
  #p95_ucl70
  
  results$coverage$p95_ucl70 <- data.frame( method = run1$coverage$method ,
                                            rep1 = run1$coverage$p95_ucl70 , 
                                            rep2 = run2$coverage$p95_ucl70 , 
                                            rep3 = run3$coverage$p95_ucl70  )
  
  
  # ading the mean
  
  results$coverage$p95_ucl70$mean <- rowMeans( results$coverage$p95_ucl70[ , 2:4 ] )
  
  # adding the sd
  
  results$coverage$p95_ucl70$cv_perc <- apply( results$coverage$p95_ucl70[ , 2:4 ] , 1 , sd )*100 / results$coverage$p95_ucl70$mean
  
  #p95_ucl95
  
  results$coverage$p95_ucl95 <- data.frame( method = run1$coverage$method ,
                                            rep1 = run1$coverage$p95_ucl95 , 
                                            rep2 = run2$coverage$p95_ucl95 , 
                                            rep3 = run3$coverage$p95_ucl95  )
  
  
  # ading the mean
  results$coverage$p95_ucl95$mean <- rowMeans( results$coverage$p95_ucl95[ , 2:4 ] )
  
  # adding the sd
  results$coverage$p95_ucl95$cv <- apply( results$coverage$p95_ucl95[ , 2:4 ] , 1 , sd )*100 / results$coverage$p95_ucl95$mean   
  
  #exceedance_ucl70
  
  results$coverage$exceedance_ucl70 <- data.frame( method = run1$coverage$method ,
                                                   rep1 = run1$coverage$exceedance_ucl70 , 
                                                   rep2 = run2$coverage$exceedance_ucl70 , 
                                                   rep3 = run3$coverage$exceedance_ucl70  )
  
  
  # ading the mean
  results$coverage$exceedance_ucl70$mean <- rowMeans( results$coverage$exceedance_ucl70[ , 2:4 ] )
  
  # adding the sd
  results$coverage$exceedance_ucl70$cv <- apply( results$coverage$exceedance_ucl70[ , 2:4 ] , 1 , sd )*100 / results$coverage$exceedance_ucl70$mean        
  
  # exceedance_ucl95
  
  results$coverage$exceedance_ucl95 <- data.frame( method = run1$coverage$method ,
                                                   rep1 = run1$coverage$exceedance_ucl95 , 
                                                   rep2 = run2$coverage$exceedance_ucl95 , 
                                                   rep3 = run3$coverage$exceedance_ucl95  )
  
  
  # ading the mean
  results$coverage$exceedance_ucl95$mean <- rowMeans( results$coverage$exceedance_ucl95[ , 2:4 ] )
  
  # adding the sd
  results$coverage$exceedance_ucl95$cv <- apply( results$coverage$exceedance_ucl95[ , 2:4 ] , 1 , sd )*100 / results$coverage$exceedance_ucl95$mean    
  
  # perc_estimable
  
  results$coverage$perc_estimable <- data.frame( method = run1$coverage$method ,
                                                 rep1 = run1$coverage$perc_estimable , 
                                                 rep2 = run2$coverage$perc_estimable , 
                                                 rep3 = run3$coverage$perc_estimable  )
  
  
  # ading the mean
  results$coverage$perc_estimable$mean <- rowMeans( results$coverage$perc_estimable[ , 2:4 ] )
  
  # adding the sd
  results$coverage$perc_estimable$cv <- apply( results$coverage$perc_estimable[ , 2:4 ] , 1 , sd )*100 / results$coverage$perc_estimable$mean                
  
  
  
  
  
  
  
  
  # perc_mistake  ------------------
  
  results$perc_mistake <- list()
  
  
  #p95_ucl70
  
  results$perc_mistake$p95_ucl70 <- data.frame( method = run1$perc_mistake$method ,
                                                rep1 = run1$perc_mistake$p95_ucl70 , 
                                                rep2 = run2$perc_mistake$p95_ucl70 , 
                                                rep3 = run3$perc_mistake$p95_ucl70  )
  
  
  # ading the mean
  
  results$perc_mistake$p95_ucl70$mean <- rowMeans( results$perc_mistake$p95_ucl70[ , 2:4 ] )
  
  # adding the sd
  
  results$perc_mistake$p95_ucl70$cv_perc <- apply( results$perc_mistake$p95_ucl70[ , 2:4 ] , 1 , sd )*100 / results$perc_mistake$p95_ucl70$mean
  
  #p95_ucl95
  
  results$perc_mistake$p95_ucl95 <- data.frame( method = run1$perc_mistake$method ,
                                                rep1 = run1$perc_mistake$p95_ucl95 , 
                                                rep2 = run2$perc_mistake$p95_ucl95 , 
                                                rep3 = run3$perc_mistake$p95_ucl95  )
  
  
  # ading the mean
  results$perc_mistake$p95_ucl95$mean <- rowMeans( results$perc_mistake$p95_ucl95[ , 2:4 ] )
  
  # adding the sd
  results$perc_mistake$p95_ucl95$cv <- apply( results$perc_mistake$p95_ucl95[ , 2:4 ] , 1 , sd )*100 / results$perc_mistake$p95_ucl95$mean   
  
  #exceedance_ucl70
  
  results$perc_mistake$exceedance_ucl70 <- data.frame( method = run1$perc_mistake$method ,
                                                       rep1 = run1$perc_mistake$exceedance_ucl70 , 
                                                       rep2 = run2$perc_mistake$exceedance_ucl70 , 
                                                       rep3 = run3$perc_mistake$exceedance_ucl70  )
  
  
  # ading the mean
  results$perc_mistake$exceedance_ucl70$mean <- rowMeans( results$perc_mistake$exceedance_ucl70[ , 2:4 ] )
  
  # adding the sd
  results$perc_mistake$exceedance_ucl70$cv <- apply( results$perc_mistake$exceedance_ucl70[ , 2:4 ] , 1 , sd )*100 / results$perc_mistake$exceedance_ucl70$mean        
  
  # exceedance_ucl95
  
  results$perc_mistake$exceedance_ucl95 <- data.frame( method = run1$perc_mistake$method ,
                                                       rep1 = run1$perc_mistake$exceedance_ucl95 , 
                                                       rep2 = run2$perc_mistake$exceedance_ucl95 , 
                                                       rep3 = run3$perc_mistake$exceedance_ucl95  )
  
  
  # ading the mean
  results$perc_mistake$exceedance_ucl95$mean <- rowMeans( results$perc_mistake$exceedance_ucl95[ , 2:4 ] )
  
  # adding the sd
  results$perc_mistake$exceedance_ucl95$cv <- apply( results$perc_mistake$exceedance_ucl95[ , 2:4 ] , 1 , sd )*100 / results$perc_mistake$exceedance_ucl95$mean    
  
  # perc_estimable
  
  results$perc_mistake$perc_estimable <- data.frame( method = run1$perc_mistake$method ,
                                                     rep1 = run1$perc_mistake$perc_estimable , 
                                                     rep2 = run2$perc_mistake$perc_estimable , 
                                                     rep3 = run3$perc_mistake$perc_estimable  )
  
  
  # ading the mean
  results$perc_mistake$perc_estimable$mean <- rowMeans( results$perc_mistake$perc_estimable[ , 2:4 ] )
  
  # adding the sd
  results$perc_mistake$perc_estimable$cv <- apply( results$perc_mistake$perc_estimable[ , 2:4 ] , 1 , sd )*100 / results$perc_mistake$perc_estimable$mean                
  
  
  
  
  
  
  
  # finalization --------------
  
  
  
  return(results)
  
  
  
}
