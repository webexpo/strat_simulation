# Script preparing performance results for the entire main analysis, with only the outcome : performance
# Main results for Sarah Exil's study

# 300 hours for me.sd=0.125
# 290 hours for me.sd=0.25


##### LIBRARIES ####

library(readxl)
library(writexl)


##### DATA ####
    
    #init_path <- "F:/Dropbox/"
    init_path <- "C:/jerome/Dropbox/"
    
    
    # results for me.cv=0.125
    me0.125 <- readRDS( file = paste( init_path ,"GITHUB/WEBEXPO/sampling_strats/EXIL TD 2024/aggregated results/cv0.125.RDS", sep="") )    
    
    # results for me.cv=0.25
    me0.25 <- readRDS( file = paste( init_path ,"GITHUB/WEBEXPO/sampling_strats/EXIL TD 2024/aggregated results/cv0.25.RDS", sep="") )    

##### Preparation ####
  
  # initialisation of one data.frame for each of the too me.sd values
  
  scenarios.0.25 <- me0.25$scenarios
  scenarios.0.125 <- me0.125$scenarios
  
  # for 0.25, adding all results as new columns
  
  # initializing additional columns
  scenarios.0.25$perf_freq_ideal_mean <- numeric(nrow(scenarios.0.25))
  scenarios.0.25$perf_freq_ideal_cv_perc <- numeric(nrow(scenarios.0.25))
  scenarios.0.25$perf_freq_naive_mean <- numeric(nrow(scenarios.0.25))
  scenarios.0.25$perf_freq_naive_cv_perc <- numeric(nrow(scenarios.0.25))
  scenarios.0.25$perf_freq_me_mean <- numeric(nrow(scenarios.0.25))
  scenarios.0.25$perf_freq_me_cv_perc <- numeric(nrow(scenarios.0.25))
  scenarios.0.25$perf_bayes_ideal_mean <- numeric(nrow(scenarios.0.25))
  scenarios.0.25$perf_bayes_ideal_cv_perc <- numeric(nrow(scenarios.0.25))
  scenarios.0.25$perf_bayes_naive_mean <- numeric(nrow(scenarios.0.25))
  scenarios.0.25$perf_bayes_naive_cv_perc <- numeric(nrow(scenarios.0.25))
  scenarios.0.25$perf_bayes_me_mean <- numeric(nrow(scenarios.0.25))
  scenarios.0.25$perf_bayes_me_cv_perc <- numeric(nrow(scenarios.0.25))
  
  scenarios.0.25$freq_estimable_perc_ideal <- numeric(nrow(scenarios.0.25))
  scenarios.0.25$freq_estimable_perc_naive <- numeric(nrow(scenarios.0.25))
  scenarios.0.25$freq_estimable_perc_me <- numeric(nrow(scenarios.0.25))
  
  # filling in the values
  
  for ( i in 1:nrow(scenarios.0.25) ) {
    
    # ideal
    scenarios.0.25$perf_freq_ideal_mean[i] <- me0.25$results[[i]]$perc_mistake$p95_ucl70$mean[2]
    scenarios.0.25$perf_freq_ideal_cv_perc[i] <- me0.25$results[[i]]$perc_mistake$p95_ucl70$cv_perc[2]
    
    scenarios.0.25$perf_bayes_ideal_mean[i] <- me0.25$results[[i]]$perc_mistake$p95_ucl70$mean[1]
    scenarios.0.25$perf_bayes_ideal_cv_perc[i] <- me0.25$results[[i]]$perc_mistake$p95_ucl70$cv_perc[1]
    
    # naive
    scenarios.0.25$perf_freq_naive_mean[i] <- me0.25$results[[i]]$perc_mistake$p95_ucl70$mean[4]
    scenarios.0.25$perf_freq_naive_cv_perc[i] <- me0.25$results[[i]]$perc_mistake$p95_ucl70$cv_perc[4]
    
    scenarios.0.25$perf_bayes_naive_mean[i] <- me0.25$results[[i]]$perc_mistake$p95_ucl70$mean[3]
    scenarios.0.25$perf_bayes_naive_cv_perc[i] <- me0.25$results[[i]]$perc_mistake$p95_ucl70$cv_perc[3]
    
    # me
    
    scenarios.0.25$perf_freq_me_mean[i] <- me0.25$results[[i]]$perc_mistake$p95_ucl70$mean[10]
    scenarios.0.25$perf_freq_me_cv_perc[i] <- me0.25$results[[i]]$perc_mistake$p95_ucl70$cv_perc[10]
    
    scenarios.0.25$perf_bayes_me_mean[i] <- me0.25$results[[i]]$perc_mistake$p95_ucl70$mean[5]
    scenarios.0.25$perf_bayes_me_cv_perc[i] <- me0.25$results[[i]]$perc_mistake$p95_ucl70$cv_perc[5]
   
    # frequentist estimability
    
    scenarios.0.25$freq_estimable_perc_ideal[i] <- me0.25$results[[i]]$perc_mistake$perc_estimable$mean[2]
    scenarios.0.25$freq_estimable_perc_naive[i] <- me0.25$results[[i]]$perc_mistake$perc_estimable$mean[4]
    scenarios.0.25$freq_estimable_perc_me[i] <- me0.25$results[[i]]$perc_mistake$perc_estimable$mean[10]
     
  }
  
  # for 0.125, adding all results as new columns
  
  # initializing additional columns
  
  scenarios.0.125$perf_freq_ideal_mean <- numeric(nrow(scenarios.0.125))
  scenarios.0.125$perf_freq_ideal_cv_perc <- numeric(nrow(scenarios.0.125))
  scenarios.0.125$perf_freq_naive_mean <- numeric(nrow(scenarios.0.125))
  scenarios.0.125$perf_freq_naive_cv_perc <- numeric(nrow(scenarios.0.125))
  scenarios.0.125$perf_freq_me_mean <- numeric(nrow(scenarios.0.125))
  scenarios.0.125$perf_freq_me_cv_perc <- numeric(nrow(scenarios.0.125))
  scenarios.0.125$perf_bayes_ideal_mean <- numeric(nrow(scenarios.0.125))
  scenarios.0.125$perf_bayes_ideal_cv_perc <- numeric(nrow(scenarios.0.125))
  scenarios.0.125$perf_bayes_naive_mean <- numeric(nrow(scenarios.0.125))
  scenarios.0.125$perf_bayes_naive_cv_perc <- numeric(nrow(scenarios.0.125))
  scenarios.0.125$perf_bayes_me_mean <- numeric(nrow(scenarios.0.125))
  scenarios.0.125$perf_bayes_me_cv_perc <- numeric(nrow(scenarios.0.125))
  
  scenarios.0.125$freq_estimable_perc_ideal <- numeric(nrow(scenarios.0.125))
  scenarios.0.125$freq_estimable_perc_naive <- numeric(nrow(scenarios.0.125))
  scenarios.0.125$freq_estimable_perc_me <- numeric(nrow(scenarios.0.125))
  
  # filling in the values
  
  for ( i in 1:nrow(scenarios.0.125) ) {
    
    # ideal
    scenarios.0.125$perf_freq_ideal_mean[i] <- me0.125$results[[i]]$perc_mistake$p95_ucl70$mean[2]
    scenarios.0.125$perf_freq_ideal_cv_perc[i] <- me0.125$results[[i]]$perc_mistake$p95_ucl70$cv_perc[2]
    
    scenarios.0.125$perf_bayes_ideal_mean[i] <- me0.125$results[[i]]$perc_mistake$p95_ucl70$mean[1]
    scenarios.0.125$perf_bayes_ideal_cv_perc[i] <- me0.125$results[[i]]$perc_mistake$p95_ucl70$cv_perc[1]
    
    # naive
    scenarios.0.125$perf_freq_naive_mean[i] <- me0.125$results[[i]]$perc_mistake$p95_ucl70$mean[4]
    scenarios.0.125$perf_freq_naive_cv_perc[i] <- me0.125$results[[i]]$perc_mistake$p95_ucl70$cv_perc[4]
    
    scenarios.0.125$perf_bayes_naive_mean[i] <- me0.125$results[[i]]$perc_mistake$p95_ucl70$mean[3]
    scenarios.0.125$perf_bayes_naive_cv_perc[i] <- me0.125$results[[i]]$perc_mistake$p95_ucl70$cv_perc[3]
    
    # me
    
    scenarios.0.125$perf_freq_me_mean[i] <- me0.125$results[[i]]$perc_mistake$p95_ucl70$mean[10]
    scenarios.0.125$perf_freq_me_cv_perc[i] <- me0.125$results[[i]]$perc_mistake$p95_ucl70$cv_perc[10]
    
    scenarios.0.125$perf_bayes_me_mean[i] <- me0.125$results[[i]]$perc_mistake$p95_ucl70$mean[5]
    scenarios.0.125$perf_bayes_me_cv_perc[i] <- me0.125$results[[i]]$perc_mist$p95_ucl70$cv_perc[5]
    
    # frequentist estimability
    
    scenarios.0.125$freq_estimable_perc_ideal[i] <- me0.125$results[[i]]$perc_mistake$perc_estimable$mean[2]
    scenarios.0.125$freq_estimable_perc_naive[i] <- me0.125$results[[i]]$perc_mistake$perc_estimable$mean[4]
    scenarios.0.125$freq_estimable_perc_me[i] <- me0.125$results[[i]]$perc_mistake$perc_estimable$mean[10]
    
  }
  
  # merging the two datasets
  
  main_dataset <- rbind( cbind( scenarios.0.25 , "me.sd" = rep("0.25",nrow(scenarios.0.25))),
                         cbind( scenarios.0.125 , "me.sd" = rep("0.125",nrow(scenarios.0.125)) ) )
  
  # additional parameters
  
  main_dataset$n_sim <- 5000
  main_dataset$n_iterations_gum <- 5000
  main_dataset$n_replications <- 3
  
  
  # reordering variables
  
  main_dataset <- main_dataset[,c(1:6,22:25,7:21) ]
  
  # creating a variable description table
  
  dictionnary <- data.frame( variable = names(main_dataset) , description = c("True GSD",
                                                                              "True exceedance fraction in %",
                                                                              "Sample size",
                                                                              "True 95th percentile",
                                                                              "True GM",
                                                                              "Occupational exposure limit",
                                                                              "Measurement error standard deviation",
                                                                              "Number of simulated datasets for each scenario",
                                                                              "Number of iterations for the GUM method",
                                                                              "Number of replications of the entire simulation",
                                                                             
                                                                             "Mistake rate of the frequentist method with ideal case (mean)",
                                                                             "Mistake rate of the frequentist method with ideal case (CV%)",
                                                                             "Mistake rate of the frequentist method with naive case (mean)",
                                                                             "Mistake rate of the frequentist method with naive case (CV%)",
                                                                             "Mistake rate of the frequentist method with measurement error (mean)",
                                                                             "Mistake rate of the frequentist method with measurement error (CV%)",
                                                                             "Mistake rate of the Bayesian method with ideal case (mean)",
                                                                             "Mistake rate of the Bayesian method with ideal case (CV%)",
                                                                             "Mistake rate of the Bayesian method with naive case (mean)",
                                                                             "Mistake rate of the Bayesian method with naive case (CV%)",
                                                                             "Mistake rate of the Bayesian method with measurement error (mean)",
                                                                             "Mistake rate of the Bayesian method with measurement error (CV%)",
                                                                             "Percentage of estimable values for the frequentist method with ideal case",
                                                                             "Percentage of estimable values for the frequentist method with naive case",
                                                                             "Percentage of estimable values for the frequentist method with measurement error"
                                                                              ) ) 


####### Data EXPORT as EXCEL files

write_xlsx( main_dataset , path = paste( init_path ,"GITHUB/WEBEXPO/sampling_strats/EXIL TD 2024/aggregated results/main_analysis.xlsx", sep="") )

write_xlsx( dictionnary , path = paste( init_path ,"GITHUB/WEBEXPO/sampling_strats/EXIL TD 2024/aggregated results/main_analysis_dictionnary.xlsx", sep="")  )
