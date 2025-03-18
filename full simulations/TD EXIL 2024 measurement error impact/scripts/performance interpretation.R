######################
#
#  Measurement error simulation - generating report format results
#  
#
#######################



# Load libraries ----------------------------

  library(readxl)
  library(writexl)

# Load data --------------------------------

    #init_path <- "F:/Dropbox/"
    init_path <- "C:/jerome/Dropbox/"
    
    # main analysis
    main_dataset <- read_xlsx( path = paste( init_path ,"GITHUB/WEBEXPO/sampling_strats/EXIL TD 2024/aggregated results/main_analysis.xlsx", sep="") )
    
    main_dataset_dict <- read_xlsx( path = paste( init_path ,"GITHUB/WEBEXPO/sampling_strats/EXIL TD 2024/aggregated results/main_analysis_dictionnary.xlsx", sep="")  )
    
    # real gsd analysis with non detects
    gsd_dataset <- read_xlsx( path = paste( init_path ,"GITHUB/WEBEXPO/sampling_strats/EXIL TD 2024/aggregated results/realGSD_analysis.xlsx", sep="") )
    
    gsd_dataset_dict <- read_xlsx(path = paste( init_path ,"GITHUB/WEBEXPO/sampling_strats/EXIL TD 2024/aggregated results/realGSD_analysis_dictionnary.xlsx", sep="")  )

    
# Reproductibility - main analyis ----------------------------------------------
    
    # creating a long format dataset 
    
    main_dataset_long <- main_dataset[ rep( 1:144 , 6 ) , c(1:10)]
    
    main_dataset_long$mean <- c( main_dataset$perf_bayes_ideal_mean , 
                                 main_dataset$perf_bayes_naive_mean,
                                 main_dataset$perf_bayes_me_mean,
                                 main_dataset$perf_freq_ideal_mean , 
                                 main_dataset$perf_freq_naive_mean,
                                 main_dataset$perf_freq_me_mean )
    
    main_dataset_long$cv_perc <- c( main_dataset$perf_bayes_ideal_cv_perc , 
                                 main_dataset$perf_bayes_naive_cv_perc,
                                 main_dataset$perf_bayes_me_cv_perc,
                                 main_dataset$perf_freq_ideal_cv_perc , 
                                 main_dataset$perf_freq_naive_cv_perc,
                                 main_dataset$perf_freq_me_cv_perc )
    
    main_dataset_long$method <- rep( c("bayes_ideal", "bayes_naive", "bayes_me", "freq_ideal", "freq_naive", "freq_me") , 144 )
                                 
    # some CVs equal to 0 are NA
    
    main_dataset_long$cv_perc[ is.na(main_dataset_long$cv_perc) ] <- 0
    
    
    # creating a flag for cv_perc < 5%
    
    main_dataset_long$cv_flag <- ifelse( main_dataset_long$cv_perc < 5 , 1 , 0 )
    
    # creating a flag for mean*cv_perc/100 < 1%
    
    main_dataset_long$cv_mean_flag <- ifelse( main_dataset_long$mean * main_dataset_long$cv_perc / 100 < 1 , 1 , 0 )

    
    # some results
    
    quantile( main_dataset_long$cv_perc , probs = c(0,0.05, 0.25, 0.5, 0.75, 0.95,1) )
    
    quantile( main_dataset_long$cv_perc*main_dataset_long$mean/100 , probs = c(0,0.05, 0.25, 0.5, 0.75, 0.95,1) )
    
    # proportion of CV_perc below 5%
    
    sum( main_dataset_long$cv_flag ) / nrow(main_dataset_long)
        
    # proprotion of absolute SD below 1%
    
    sum( main_dataset_long$cv_mean_flag ) / nrow(main_dataset_long)

    # proportion of CV_perc below 5% or absolute SD below 1%
    
    sum( main_dataset_long$cv_flag | main_dataset_long$cv_mean_flag ) / nrow(main_dataset_long)

    
    
# Reproductibility - GSD analysis ------------------------------------
    
    # creating a long format dataset 
    
    gsd_dataset_long <- gsd_dataset[ rep( 1:144 , 6 ) , c(1:10)]
    
    gsd_dataset_long$mean <- c( gsd_dataset$perf_bayes_ideal_mean , 
                                 gsd_dataset$perf_bayes_naive_mean,
                                 gsd_dataset$perf_bayes_me_mean,
                                 gsd_dataset$perf_freq_ideal_mean , 
                                 gsd_dataset$perf_freq_naive_mean,
                                 gsd_dataset$perf_freq_me_mean )
    
    gsd_dataset_long$cv_perc <- c( gsd_dataset$perf_bayes_ideal_cv_perc , 
                                 gsd_dataset$perf_bayes_naive_cv_perc,
                                 gsd_dataset$perf_bayes_me_cv_perc,
                                 gsd_dataset$perf_freq_ideal_cv_perc , 
                                 gsd_dataset$perf_freq_naive_cv_perc,
                                 gsd_dataset$perf_freq_me_cv_perc )
    
    gsd_dataset_long$method <- rep( c("bayes_ideal", "bayes_naive", "bayes_me", "freq_ideal", "freq_naive", "freq_me") , 144 )
                                 
    # some CVs equal to 0 are NA
    
    gsd_dataset_long$cv_perc[ is.na(gsd_dataset_long$cv_perc) ] <- 0
    
    
    # creating a flag for cv_perc < 5%
    
    gsd_dataset_long$cv_flag <- ifelse( gsd_dataset_long$cv_perc < 5 , 1 , 0 )
    
    # creating a flag for mean*cv_perc/100 < 1%
    
    gsd_dataset_long$cv_mean_flag <- ifelse( gsd_dataset_long$mean * gsd_dataset_long$cv_perc / 100 < 1 , 1 , 0 )

    
    # some results
    
    quantile( gsd_dataset_long$cv_perc , probs = c(0,0.05, 0.25, 0.5, 0.75, 0.95,1) )
    
    quantile( gsd_dataset_long$cv_perc*gsd_dataset_long$mean/100 , probs = c(0,0.05, 0.25, 0.5, 0.75, 0.95,1) )
    
    # proportion of CV_perc below 5%
    
    sum( gsd_dataset_long$cv_flag ) / nrow(gsd_dataset_long)
    
    # proprotion of absolute SD below 1%
    
    sum( gsd_dataset_long$cv_mean_flag ) / nrow(gsd_dataset_long)
    
    # proportion of CV_perc below 5% or absolute SD below 1%
    
    sum( gsd_dataset_long$cv_flag | gsd_dataset_long$cv_mean_flag ) / nrow(gsd_dataset_long)
    
    
# Estimability ---------------------------------------------------------------------------
    
    quantile(gsd_dataset$freq_estimable_perc_ideal[ gsd_dataset$proportion_censored==0.3 ], probs = c(0,0.05, 0.25, 0.5, 0.75, 0.95,1) )
    
    quantile(gsd_dataset$freq_estimable_perc_ideal[ gsd_dataset$proportion_censored==0.6 ], probs = c(0,0.05, 0.25, 0.5, 0.75, 0.95,1) )
    
    
    quantile(gsd_dataset$freq_estimable_perc_naive[ gsd_dataset$proportion_censored==0.3 ], probs = c(0,0.05, 0.25, 0.5, 0.75, 0.95,1) )
    
    quantile(gsd_dataset$freq_estimable_perc_naive[ gsd_dataset$proportion_censored==0.6 ], probs = c(0,0.05, 0.25, 0.5, 0.75, 0.95,1) )
    
    
    quantile(gsd_dataset$freq_estimable_perc_me[ gsd_dataset$proportion_censored==0.3 ], probs = c(0,0.05, 0.25, 0.5, 0.75, 0.95,1) )
    
    quantile(gsd_dataset$freq_estimable_perc_me[ gsd_dataset$proportion_censored==0.6 ], probs = c(0,0.05, 0.25, 0.5, 0.75, 0.95,1) )
    
    