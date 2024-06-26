#' exploring the impact of measurement error industrial hygiene measurement data interpretation

#' question : for one example of true distribution, what is the impact of a 50% expanded measurement error (CV=50/1.96) on P95 and its UCL

#' This is an extension of the AIOH2023.r analysis, triggered by exchanges with Theo Scheffers who described a GUM approach to measurement error handling

#' This script presents the running of the simulation, which was run 5 times

##### Data #####

    ## libraries  
    
    library(ggplot2)
    library(ggthemes)
    library(writexl)
    library(scales)
    library(readxl)
    library(rjags)
    library(parallel)
    library(tolerance)
    devtools::source_url("https://github.com/lhimp/scripts/raw/master/chemin.R")


    ##  Webexpo scripts      
    
    # data preparation / generation
    
    chemin(
      fileName = "webexpo.seg.randomgeneration.R",
      relPath = c("RANDOM SAMPLE GENERATION"),
      githubCredentials = list(userName = "webexpo", repoName = "webexpo_r_lib")
    )
    
    chemin(
      fileName = "webexpo.seg.dataprep.R",
      relPath = c("DATA PREPARATION", "SEG ANALYSIS"),
      githubCredentials = list(userName = "webexpo", repoName = "webexpo_r_lib")
    )
    
    # JAGS models
    
    chemin(
      fileName = "webexpo.seg.mainbayesian.R",
      relPath = c("jags models", "SEG ANALYSIS"),
      githubCredentials = list(userName = "webexpo", repoName = "webexpo_r_lib")
    )
    
    
    chemin(
      fileName = "webexpo.seg.informedvarbayesian.R",
      relPath = c("jags models", "SEG ANALYSIS"),
      githubCredentials = list(userName = "webexpo", repoName = "webexpo_r_lib")
    )
    
    chemin(
      fileName = "webexpo.seg.informedvarbayesian.models.R",
      relPath = c("jags models", "SEG ANALYSIS"),
      githubCredentials = list(userName = "webexpo", repoName = "webexpo_r_lib")
    )   
    
    # Data interpretation
    
    
    chemin(
      fileName = "webexpo.seg.interpretation.R",
      relPath = c("RESULT INTERPRETATION", "SEG ANALYSIS"),
      githubCredentials = list(userName = "webexpo", repoName = "webexpo_r_lib")
    ) 
    
    
    # ROS
    
    chemin(
      fileName = "function NDexpo_EN.R",
      relPath = c("NDEXPO"),
      githubCredentials = list(userName = "webexpo", repoName = "webexpo_r_lib")
    ) 

    
    # frequentit estimation
    
    source("parameter estimation/frequentist/percentile.R")
    
##### Functions #####

    #' smal function for generating MCMC results for an analysis assuming no error
    #'
    #' @param true_gsd 
    #' @param true_gm 
    #' @param sample_size 
    #' @param n_sim number of simulations 
    #' @return matrix of generated data, one column per iteration
    #'
    
    
    data.simulation <- function(   true_gsd = 2.5, 
                                   true_gm = 30, 
                                   n_sim = 50, 
                                   sample_size = 9 ,
                                   me.cv = 0.15) {
      
      ###### data generation
      
      
      data_vector <- exp(rnorm( n = sample_size*n_sim , mean = log(true_gm) , sd = log(true_gsd) ) )           
      
      data_vector_me <- pmax( rep(0.1,sample_size*n_sim)  , rnorm( sample_size*n_sim , data_vector , data_vector*me.cv ) )
      
      #formatting pre treatment
      
      data_matrix <- matrix( data = data_vector , nrow = sample_size , ncol = n_sim)
      
      data_matrix_me <- matrix( data = data_vector_me , nrow = sample_size , ncol = n_sim)
      
      return(list(true = data_matrix, observed = data_matrix_me))
      
    }
    
    
    
     

    #' smal function for generating MCMC results for an analysis assuming no error
    #'
    #' @param mysample sample to be analysed
    #' @param oel occupational exposure limit
    #'
    #' @return MCMC chains for GM, GSD and P95
    #'
    
    myfunction.naive <- function( mysample , oel) {
      
      mcmc <- Webexpo.seg.globalbayesian.jags( data.sample = mysample ,
                                               is.lognormal = TRUE , 
                                               error.type = "none" ,
                                               me.range = c(0.3,0.3) , 
                                               oel = oel ,
                                               prior.model = "informedvar")
      
      result <- list( gm_chain = exp(mcmc$mu.chain),
                      gsd_chain = exp(mcmc$sigma.chain),
                      p95_chain = exp(mcmc$mu.chain+qnorm(0.95)*mcmc$sigma.chain))
      
      return(result)
      
    }
    
    #' small function for generating MCMC results for an analysis assuming CV error within a range
    #'
    #' @param mysample sample to be analysed
    #' @param me.range range of error standard measurement error as CV
    #' @param oel occupational exposure limit
    #'
    #' @return MCMC chains for GM, GSD and P95
    #'
    
    myfunction.me <- function( mysample , me.range = c(0.2,0.4) , oel ) {
      
      mcmc <- Webexpo.seg.globalbayesian.jags( data.sample = mysample ,
                                               is.lognormal = TRUE , 
                                               error.type = "CV" ,
                                               me.range = me.range , 
                                               oel = oel ,
                                               prior.model = "informedvar")
      
      result <- list( gm_chain = exp(mcmc$mu.chain),
                      gsd_chain = exp(mcmc$sigma.chain),
                      p95_chain = exp(mcmc$mu.chain+qnorm(0.95)*mcmc$sigma.chain))

      
      return(result)
      
    }
    
    #' small function for using monte carlo simulation as described by Scheffers et al.
    #'
    #' @param mysample sample to be analysed
    #' @param me.cv standard measurement error as CV
    #' @param n.simul simulation size
    #'
    #' @return MCMC chains for GM, GSD and P95
    #'
    
    myfunction.gum <- function( mysample , n.simul , me.cv ) {
      
      # update june 25, correction of P95 UCL equation to mimick EN689, using the tolerance package
      # as of this update the simulations wer enot run with it. So frequentist coverage needs correction 
      
      
      # mysample <- c("28.9","19.40","<8.22","149.9","26.42","56.1")
      # n.simul <- 10000
      
      # ROS
      
      mysample.ros <- fun.NdExpo.lognorm( mysample )$data$xfin
      
      
      little.ros.gum <- function(i) {
      
      observed_data <- pmax( 0.1 , rnorm( length(mysample.ros) , mysample.ros , mysample.ros*me.cv ) ) #truncation at 0.1
        
      gm.ros <- exp(mean(log(observed_data)))
      
      gsd.ros <- exp(sd(log(observed_data)))
      
      p95.ros <- fun.perc(observed_data,alpha=0.05,perc=0.95)$est
      
      p95ucl70.ros <- fun.perc.en689(observed_data,alpha=0.30,perc=0.95)$uc
      
      return( c(gm.ros,gsd.ros,p95.ros,p95ucl70.ros) )
      
      }
      
      repetition <- sapply(1:n.simul, little.ros.gum)
      
      raw_result <- apply( repetition , 1 , mean)
      
      result <- list( gm = raw_result[1],
                      gsd = raw_result[2],
                      p95 = raw_result[3],
                      p95ucl70 = raw_result[4] )
      
      return(result)
      
    }
    
    
    #' small function for frequentist estimation of GM, GSD and P95 
    #'
    #' @param mysample sample to be analysed
    #'
    #' @return GM, GSD and P95
    #'
    
    myfunction.freq <- function( mysample  ) {
      
      # update june 25, correction of P95 UCL equation to mimick EN689, using the tolerance package
      # as of this update the simulations wer enot run with it. So frequentist coverage needs correction 
      
      # mysample <- c("28.9","19.40","<8.22","149.9","26.42","56.1")
      # n.simul <- 10000
      
      # ROS
      
      mysample.ros <- fun.NdExpo.lognorm( mysample )$data$xfin
      
      observed_data <- pmax( 0.1 , rnorm( length(mysample.ros) , mysample.ros , mysample.ros*me.cv ) ) #truncation at 0.1
        
        gm.ros <- exp(mean(log(mysample.ros)))
        
        gsd.ros <- exp(sd(log(mysample.ros)))
        
        p95.ros <- fun.perc(mysample.ros,alpha=0.05,perc=0.95)$est
        
        p95ucl70.ros <- fun.perc.en689(mysample.ros,alpha=0.30,perc=0.95)$uc
        
      result <- list( gm = gm.ros,
                      gsd = gsd.ros,
                      p95 = p95.ros,
                      p95ucl70 = p95ucl70.ros )
      
      return(result)
      
    }
    
    
    
    
    #' function repeatedly estimating P95 and its 70% UTL for ideal, naive and measurement error analysis
    #'
    #' @param sample_size 
    #' @param true_p95 
    #' @param true_gsd
    #' @param me.cv value for the measurement error as standard error CV, not in % and not expanded uncertainty CV
    #' @param oel occupational exposure limit
    #'
    #' @return point estimate and 70% UCL for p95 using the three approaches
    #'
    
    
    
    my.parallel.function <- function( index ) {
      
      # data generation
      
      true_gm <- exp( log(true_p95) - qnorm(0.95)*log(true_gsd) )
      
      true_data <- exp( rnorm( sample_size , log(true_gm) , log(true_gsd) ) )   
      
      observed_data <- pmax( rep(0.1,sample_size) , rnorm( sample_size , true_data , true_data*me.cv ) ) #truncation at 0.1
      
      ## analysis
      
      ideal_analysis <- myfunction.naive( true_data , oel)
      
      naive_analysis <- myfunction.naive( observed_data , oel)
      
      me_analysis <- myfunction.me( observed_data , me.range = c(me.cv.range.factor*me.cv,me.cv/me.cv.range.factor) ,oel)
      
      gum_analysis <- myfunction.gum( observed_data , 10000 , me.cv)
      
      ideal_analysis_freq <- myfunction.freq( true_data )
      
      naive_analysis_freq <- myfunction.freq( observed_data )
      
      
      
      ## data interpretation
      
      results <- c( median(ideal_analysis$p95_chain),quantile(ideal_analysis$p95_chain,0.7),
                    median(naive_analysis$p95_chain),quantile(naive_analysis$p95_chain,0.7),
                    median(me_analysis$p95_chain),quantile(me_analysis$p95_chain,0.7),
                    gum_analysis$p95,gum_analysis$p95ucl70,
                    ideal_analysis_freq$p95,ideal_analysis_freq$p95ucl70,
                    naive_analysis_freq$p95,naive_analysis_freq$p95ucl70,
                    
                    median(ideal_analysis$gm_chain),
                    median(naive_analysis$gm_chain),
                    median(me_analysis$gm_chain),
                    gum_analysis$gm,
                    ideal_analysis_freq$gm,
                    naive_analysis_freq$gm,
                    
                    median(ideal_analysis$gsd_chain),
                    median(naive_analysis$gsd_chain),
                    median(me_analysis$gsd_chain),
                    gum_analysis$gsd,
                    ideal_analysis_freq$gsd,
                    naive_analysis_freq$gsd  )
      
      names(results) <- c("ideal_p95","ideal_p95UCL70","naive_p95","naive_p95UCL70","me_p95","me_p95UCL70","gum_p95","gum_p95UCL70",
                          "ideal_p95_freq","ideal_p95UCL70_freq","naive_p95_freq","naive_p95UCL70_freq",
                          "ideal_gm","naive_gm","me_gm","gum_gm","ideal_gm_freq","naive_gm_freq",
                          "ideal_gsd","naive_gsd","me_gsd","gum_gsd","ideal_gsd_freq","naive_gsd_freq")
      
      return( results) }
    
#### Analysis #####
    
  # note : each of the simulation below was run 5 times
    
  # the raw simulation files are located on dropbox : Dropbox\GITHUB\WEBEXPO\sampling_strats\GUM measurement error 2024   
  
  ###### Example 1 ######  
    
    
    ## parameters
    
    expanded_uncertainty <- 0.50
    
    coverage_factor <- qnorm(0.975)
    
    me.cv.range.factor <- 1 #( <1, if one wants to input uncertain uncertainty in the ME analysis)
    
    me.cv <- me.cv.range.factor*expanded_uncertainty/coverage_factor
    
    sample_size <- 6
    
    true_p95 <- 100
    
    true_gsd <- 2.5
    
    oel <- 300
    
    true_gm <- exp( log(true_p95) - qnorm(0.95)*log(true_gsd) )
    
    ## data generated
    
    true_data <- c(15.446956, 121.058372,  39.942332,   9.226556, 197.908315,  22.717817 )   
    
    observed_data <- c(12.793548, 155.902654,  40.918517,   3.873618, 219.697256,  29.277506 )
    
    ## Bayesian analysis
    
    ideal_analysis_b <- myfunction.naive( true_data , oel)
    
    ideal_analysis_f <- myfunction.freq( true_data )
    
    naive_analysis_b <- myfunction.naive( observed_data , oel)
    
    naive_analysis_f <- myfunction.freq( observed_data )
    
    me_analysis_b <- myfunction.me( observed_data , me.range = c(me.cv.range.factor*expanded_uncertainty/coverage_factor,
                                                               expanded_uncertainty/(coverage_factor*me.cv.range.factor)) ,oel)
    
    me_analysis_f <- myfunction.gum( observed_data , 10000 , expanded_uncertainty/coverage_factor)
    
    ##numerical results - P95 estimates
    
    exemple1_table <- data.frame( type = c("true","ideal_b","ideal_f","naive_b","naive_f","me_b","me_f"),
                                  
                                  gm = c( true_gm,
                                          median(ideal_analysis_b$gm_chain),
                                          ideal_analysis_f$gm,
                                          median(naive_analysis_b$gm_chain),
                                          naive_analysis_f$gm,
                                          median(me_analysis_b$gm_chain),
                                          me_analysis_f$gm),
                                  
                                  gsd = c( true_gsd,
                                           median(ideal_analysis_b$gsd_chain),
                                           ideal_analysis_f$gsd,
                                           median(naive_analysis_b$gsd_chain),
                                           naive_analysis_f$gsd,
                                           median(me_analysis_b$gsd_chain),
                                           me_analysis_f$gsd),
                                  
                                  
                                  
                                  p95 = c( true_p95,
                                           quantile(ideal_analysis_b$p95_chain,0.5),
                                           ideal_analysis_f$p95,
                                           quantile(naive_analysis_b$p95_chain,0.5),
                                           naive_analysis_f$p95,
                                           quantile(me_analysis_b$p95_chain,0.5),
                                           me_analysis_f$p95),
                                  
                                  p95_ucl = c( NA,
                                               quantile(ideal_analysis_b$p95_chain,0.7),
                                               ideal_analysis_f$p95ucl70,
                                               quantile(naive_analysis_b$p95_chain,0.7),
                                               naive_analysis_f$p95ucl70,
                                               quantile(me_analysis_b$p95_chain,0.7),
                                               me_analysis_f$p95ucl70)
                                  
    )
    
    
    exemple1.object <- list( true_data = true_data,
                             observed_data = observed_data,
                             results = exemple1_table)
    
    saveRDS( exemple1.object, "created data/GUM measurement error_ex1.RDS")
    

###### Example 2 ######  
    
    ## example from Theo et al. : default example from expostats, with <5.5 replaced with 8.22 (unlike NDexpo)
    
    ## parameters
    
    me.cv <- 0.15
    
    sample_size <- 6
    
    true_data <- c(28.9, 19.4, 8.22, 149.9, 26.42, 56.1 )
  
    observed_data <- c(35.54, 20.73, 9.36, 136.83, 28.84, 47.27 )
    
    oel <-100
    
    ## Bayesian analysis
    
    ideal_analysis_b <- myfunction.naive( true_data , oel)
    
    ideal_analysis_f <- myfunction.freq( true_data )
    
    naive_analysis_b <- myfunction.naive( observed_data , oel)
    
    naive_analysis_f <- myfunction.freq( observed_data )
    
    me_analysis_b <- myfunction.me( observed_data , me.range = c(me.cv,me.cv) ,oel)
    
    me_analysis_f <- myfunction.gum( observed_data , 10000 , me.cv)
    
    ##numerical results - P95 estimates
    
    exemple2_table <- data.frame( type = c("true","ideal_b","ideal_f","naive_b","naive_f","me_b","me_f"),
                                  
                                  gm = c( true_gm,
                                          median(ideal_analysis_b$gm_chain),
                                          ideal_analysis_f$gm,
                                          median(naive_analysis_b$gm_chain),
                                          naive_analysis_f$gm,
                                          median(me_analysis_b$gm_chain),
                                          me_analysis_f$gm),
                                  
                                  gsd = c( true_gsd,
                                           median(ideal_analysis_b$gsd_chain),
                                           ideal_analysis_f$gsd,
                                           median(naive_analysis_b$gsd_chain),
                                           naive_analysis_f$gsd,
                                           median(me_analysis_b$gsd_chain),
                                           me_analysis_f$gsd),
                                  
                                  
                                  
                                  p95 = c( true_p95,
                                           quantile(ideal_analysis_b$p95_chain,0.5),
                                           ideal_analysis_f$p95,
                                           quantile(naive_analysis_b$p95_chain,0.5),
                                           naive_analysis_f$p95,
                                           quantile(me_analysis_b$p95_chain,0.5),
                                           me_analysis_f$p95),
                                  
                                  p95_ucl = c( NA,
                                               quantile(ideal_analysis_b$p95_chain,0.7),
                                               ideal_analysis_f$p95ucl70,
                                               quantile(naive_analysis_b$p95_chain,0.7),
                                               naive_analysis_f$p95ucl70,
                                               quantile(me_analysis_b$p95_chain,0.7),
                                               me_analysis_f$p95ucl70)
                                  
    )
    
    
    exemple2.object <- list( true_data = true_data,
                             observed_data = observed_data,
                             results = exemple2_table)
    
    saveRDS( exemple2.object, "created data/GUM measurement error_ex2.RDS")
    
    
    ## replication of the results provided by Theo :  GM, GSD and UTL95/70 replicated
    
    ## GUM UTL : apply the GUM procedure to the TRUE data, Theo's results replicated
    
          
  ###### simulation 1 ##### 
    
    ## impact on estimation of gm, gsd, p95, and the p95 70% UCL
    
    expanded_uncertainty <- 0.50
    
    coverage_factor <- qnorm(0.975)
    
    me.cv.range.factor <- 1 #( <1, if one wants to input uncertain uncertainty in the ME analysis)
    
    sample_size <- 6
    
    true_p95 <- 100
    
    true_gsd <- 2.5
    
    oel <- 300
    
    me.cv <- expanded_uncertainty/coverage_factor
   
    n_simul <- 10000 
    
    ##  traditional simulation
    
        # compteur de temps initialisé
        #start_time <- Sys.time()
    
        #simulation_result <- lapply( X = as.list(1:n_simul) , FUN = my.parallel.function)
    
        # estimation of computing time ( 2 min on my computer for n=1000)
        #end_time <- Sys.time()
        #mytime <- end_time - start_time 
    
    
    ##  parallel simulation (already performed and saved)
    
        # compteur de temps initialisé
        start_time <- Sys.time()
        
        
        #procédure parallele
        
        # creating the clusters ( use detectCores() to count the cores available, choose this number minus 2 or 4 below) 
        
        cl <- makeCluster(9)
        
        # libraries and scripts to be used in each cluster
        
        clusterEvalQ(cl, library(rjags))
        
        
        # sending objects to clusters
        clusterExport( cl , "jags.model.informedvar" , envir=environment())
        clusterExport( cl , "webexpo.seg.datapreparation" , envir=environment())
        clusterExport( cl , "Webexpo.seg.globalbayesian.jags" , envir=environment())
        clusterExport( cl , "fun.jags.informedvar" , envir=environment())
        clusterExport( cl , "myfunction.naive" , envir=environment())
        clusterExport( cl , "myfunction.me" , envir=environment())
        clusterExport( cl , "myfunction.me" , envir=environment())
        clusterExport( cl , "myfunction.gum" , envir=environment())
        clusterExport( cl , "myfunction.freq" , envir=environment())
        clusterExport( cl , "fun.NdExpo.lognorm" , envir=environment())
        clusterExport( cl , "fun.perc" , envir=environment())
        
        clusterExport( cl , "sample_size" , envir=environment())
        clusterExport( cl , "me.cv.range.factor" , envir=environment())
        clusterExport( cl , "true_p95" , envir=environment())
        clusterExport( cl , "true_gsd" , envir=environment())
        clusterExport( cl , "oel" , envir=environment())
        clusterExport( cl , "me.cv" , envir=environment())

        # calculations
        
        simulation_result_par <- parLapply(cl, X = as.list(1:n_simul) , fun = my.parallel.function) 
        
        # recommendation from the net: close the clusters
        stopCluster(cl)
        
        
        # estimation of computing time ( 9 min on my computer for 5000 iterations)
        end_time <- Sys.time()
        mytime <- end_time - start_time 
    
        
      
    ## interpretation of results 
        
        
    simulation_summary <- data.frame( ideal_p95 = numeric(n_simul),
                                               ideal_p95_ucl = numeric(n_simul),
                                               naive_p95 = numeric(n_simul),
                                               naive_p95_ucl = numeric(n_simul),
                                               me_p95 = numeric(n_simul),
                                               me_p95_ucl = numeric(n_simul),
                                               gum_p95 = numeric(n_simul),
                                               gum_p95_ucl = numeric(n_simul),
                                              ideal_p95_freq = numeric(n_simul),
                                              ideal_p95_ucl_freq = numeric(n_simul),
                                              naive_p95_freq = numeric(n_simul),
                                              naive_p95_ucl_freq = numeric(n_simul),
                                                        
                                                ideal_gm = numeric(n_simul),
                                                naive_gm = numeric(n_simul),
                                                me_gm = numeric(n_simul),
                                                gum_gm = numeric(n_simul),
                                                ideal_gm_freq = numeric(n_simul),
                                                naive_gm_freq = numeric(n_simul),
                                                  
                                                ideal_gsd = numeric(n_simul),
                                                naive_gsd = numeric(n_simul),
                                                me_gsd = numeric(n_simul),
                                                gum_gsd = numeric(n_simul),
                                                ideal_gsd_freq = numeric(n_simul),
                                                naive_gsd_freq = numeric(n_simul))
        for ( i in 1:n_simul) { 
          
          simulation_summary$ideal_p95[i] <- simulation_result_par[[i]][1]
          simulation_summary$ideal_p95_ucl[i] <- simulation_result_par[[i]][2]
          simulation_summary$naive_p95[i] <- simulation_result_par[[i]][3]
          simulation_summary$naive_p95_ucl[i] <- simulation_result_par[[i]][4]
          simulation_summary$me_p95[i] <- simulation_result_par[[i]][5]
          simulation_summary$me_p95_ucl[i] <- simulation_result_par[[i]][6]
          simulation_summary$gum_p95[i] <- simulation_result_par[[i]][7]
          simulation_summary$gum_p95_ucl[i] <- simulation_result_par[[i]][8]
          simulation_summary$ideal_p95_freq[i] <- simulation_result_par[[i]][9]
          simulation_summary$ideal_p95_ucl_freq[i] <- simulation_result_par[[i]][10]
          simulation_summary$naive_p95_freq[i] <- simulation_result_par[[i]][11]
          simulation_summary$naive_p95_ucl_freq[i] <- simulation_result_par[[i]][12]
          
          
          
          simulation_summary$ideal_gm[i] = simulation_result_par[[i]][13]
          simulation_summary$naive_gm[i] = simulation_result_par[[i]][14]
          simulation_summary$me_gm[i] = simulation_result_par[[i]][15]
          simulation_summary$gum_gm[i] = simulation_result_par[[i]][16]
          simulation_summary$ideal_gm_freq[i] = simulation_result_par[[i]][17]
          simulation_summary$naive_gm_freq[i] = simulation_result_par[[i]][18]
          
          
          simulation_summary$ideal_gsd[i] = simulation_result_par[[i]][19]
          simulation_summary$naive_gsd[i] = simulation_result_par[[i]][20]
          simulation_summary$me_gsd[i] = simulation_result_par[[i]][21]
          simulation_summary$gum_gsd[i] = simulation_result_par[[i]][22]
          simulation_summary$ideal_gsd_freq[i] = simulation_result_par[[i]][23]
          simulation_summary$naive_gsd_freq[i] = simulation_result_par[[i]][24]
          
        }

    
    ## saving simulation results
    
    #saveRDS( simulation_summary, "F:/Dropbox/temp/aioh2023-S2_1e.RDS")
    
    


###### simulation 2 ###### 
    
    ## impact on estimation of gm, gsd, p95, and the p95 70% UCL
    
    expanded_uncertainty <- 0.50
    
    coverage_factor <- qnorm(0.975)
    
    me.cv.range.factor <- 1 #( <1, if one wants to input uncertain uncertainty in the ME analysis)
    
    sample_size <- 3
    
    true_p95 <- 100
    
    true_gsd <- 2.5
    
    oel <- 300
    
    me.cv <- expanded_uncertainty/coverage_factor
    
    n_simul <- 10000 
    
    
    
    ##  parallel simulation (already performed and saved)
    
    # compteur de temps initialisé
    start_time <- Sys.time()
    
    
    #procédure parallele
    
    # creating the clusters ( use detectCores() to count the cores available, choose this number minus 2 or 4 below) 
    
    cl <- makeCluster(9)
    
    # libraries and scripts to be used in each cluster
    
    clusterEvalQ(cl, library(rjags))
    
    
    # sending objects to clusters
    clusterExport( cl , "jags.model.informedvar" , envir=environment())
    clusterExport( cl , "webexpo.seg.datapreparation" , envir=environment())
    clusterExport( cl , "Webexpo.seg.globalbayesian.jags" , envir=environment())
    clusterExport( cl , "fun.jags.informedvar" , envir=environment())
    clusterExport( cl , "myfunction.naive" , envir=environment())
    clusterExport( cl , "myfunction.me" , envir=environment())
    clusterExport( cl , "myfunction.me" , envir=environment())
    clusterExport( cl , "myfunction.gum" , envir=environment())
    clusterExport( cl , "myfunction.freq" , envir=environment())
    
    clusterExport( cl , "fun.NdExpo.lognorm" , envir=environment())
    clusterExport( cl , "fun.perc" , envir=environment())
    
    clusterExport( cl , "sample_size" , envir=environment())
    clusterExport( cl , "me.cv.range.factor" , envir=environment())
    clusterExport( cl , "true_p95" , envir=environment())
    clusterExport( cl , "true_gsd" , envir=environment())
    clusterExport( cl , "oel" , envir=environment())
    clusterExport( cl , "me.cv" , envir=environment())
    
    # calculations
    
    simulation_result_par <- parLapply(cl, X = as.list(1:n_simul) , fun = my.parallel.function) 
    
    # recommendation from the net: close the clusters
    stopCluster(cl)
    
    
    # estimation of computing time ( 9 min on my computer for 5000 iterations)
    end_time <- Sys.time()
    mytime <- end_time - start_time 
    
    
    
    ## interpretation of results 
    
    simulation_summary <- data.frame( ideal_p95 = numeric(n_simul),
                                      ideal_p95_ucl = numeric(n_simul),
                                      naive_p95 = numeric(n_simul),
                                      naive_p95_ucl = numeric(n_simul),
                                      me_p95 = numeric(n_simul),
                                      me_p95_ucl = numeric(n_simul),
                                      gum_p95 = numeric(n_simul),
                                      gum_p95_ucl = numeric(n_simul),
                                      ideal_p95_freq = numeric(n_simul),
                                      ideal_p95_ucl_freq = numeric(n_simul),
                                      naive_p95_freq = numeric(n_simul),
                                      naive_p95_ucl_freq = numeric(n_simul),
                                      
                                      ideal_gm = numeric(n_simul),
                                      naive_gm = numeric(n_simul),
                                      me_gm = numeric(n_simul),
                                      gum_gm = numeric(n_simul),
                                      ideal_gm_freq = numeric(n_simul),
                                      naive_gm_freq = numeric(n_simul),
                                      
                                      ideal_gsd = numeric(n_simul),
                                      naive_gsd = numeric(n_simul),
                                      me_gsd = numeric(n_simul),
                                      gum_gsd = numeric(n_simul),
                                      ideal_gsd_freq = numeric(n_simul),
                                      naive_gsd_freq = numeric(n_simul))
    for ( i in 1:n_simul) { 
      
      simulation_summary$ideal_p95[i] <- simulation_result_par[[i]][1]
      simulation_summary$ideal_p95_ucl[i] <- simulation_result_par[[i]][2]
      simulation_summary$naive_p95[i] <- simulation_result_par[[i]][3]
      simulation_summary$naive_p95_ucl[i] <- simulation_result_par[[i]][4]
      simulation_summary$me_p95[i] <- simulation_result_par[[i]][5]
      simulation_summary$me_p95_ucl[i] <- simulation_result_par[[i]][6]
      simulation_summary$gum_p95[i] <- simulation_result_par[[i]][7]
      simulation_summary$gum_p95_ucl[i] <- simulation_result_par[[i]][8]
      simulation_summary$ideal_p95_freq[i] <- simulation_result_par[[i]][9]
      simulation_summary$ideal_p95_ucl_freq[i] <- simulation_result_par[[i]][10]
      simulation_summary$naive_p95_freq[i] <- simulation_result_par[[i]][11]
      simulation_summary$naive_p95_ucl_freq[i] <- simulation_result_par[[i]][12]
      
      
      
      simulation_summary$ideal_gm[i] = simulation_result_par[[i]][13]
      simulation_summary$naive_gm[i] = simulation_result_par[[i]][14]
      simulation_summary$me_gm[i] = simulation_result_par[[i]][15]
      simulation_summary$gum_gm[i] = simulation_result_par[[i]][16]
      simulation_summary$ideal_gm_freq[i] = simulation_result_par[[i]][17]
      simulation_summary$naive_gm_freq[i] = simulation_result_par[[i]][18]
      
      
      simulation_summary$ideal_gsd[i] = simulation_result_par[[i]][19]
      simulation_summary$naive_gsd[i] = simulation_result_par[[i]][20]
      simulation_summary$me_gsd[i] = simulation_result_par[[i]][21]
      simulation_summary$gum_gsd[i] = simulation_result_par[[i]][22]
      simulation_summary$ideal_gsd_freq[i] = simulation_result_par[[i]][23]
      simulation_summary$naive_gsd_freq[i] = simulation_result_par[[i]][24]
      
    }
    
    
    ## saving simulation results
    
    #saveRDS( simulation_summary, "F:/Dropbox/temp/aioh2023-S2_2e.RDS")    
    
    
    
    ###### simulation 3 ##### 
    
    ## impact on estimation of gm, gsd, p95, and the p95 70% UCL
    
    expanded_uncertainty <- 0.50
    
    coverage_factor <- qnorm(0.975)
    
    me.cv.range.factor <- 1 #( <1, if one wants to input uncertain uncertainty in the ME analysis)
    
    sample_size <- 6
    
    true_p95 <- 100
    
    true_gsd <- 1.5
    
    oel <- 300
    
    me.cv <- expanded_uncertainty/coverage_factor
    
    n_simul <- 10000 
    
    ##  traditional simulation
    
    # compteur de temps initialisé
    #start_time <- Sys.time()
    
    #simulation_result <- lapply( X = as.list(1:n_simul) , FUN = my.parallel.function)
    
    # estimation of computing time ( 2 min on my computer for n=1000)
    #end_time <- Sys.time()
    #mytime <- end_time - start_time 
    
    
    ##  parallel simulation (already performed and saved)
    
    # compteur de temps initialisé
    start_time <- Sys.time()
    
    
    #procédure parallele
    
    # creating the clusters ( use detectCores() to count the cores available, choose this number minus 2 or 4 below) 
    
    cl <- makeCluster(9)
    
    # libraries and scripts to be used in each cluster
    
    clusterEvalQ(cl, library(rjags))
    
    
    # sending objects to clusters
    clusterExport( cl , "jags.model.informedvar" , envir=environment())
    clusterExport( cl , "webexpo.seg.datapreparation" , envir=environment())
    clusterExport( cl , "Webexpo.seg.globalbayesian.jags" , envir=environment())
    clusterExport( cl , "fun.jags.informedvar" , envir=environment())
    clusterExport( cl , "myfunction.naive" , envir=environment())
    clusterExport( cl , "myfunction.me" , envir=environment())
    clusterExport( cl , "myfunction.me" , envir=environment())
    clusterExport( cl , "myfunction.gum" , envir=environment())
    clusterExport( cl , "myfunction.freq" , envir=environment())
    clusterExport( cl , "fun.NdExpo.lognorm" , envir=environment())
    clusterExport( cl , "fun.perc" , envir=environment())
    
    clusterExport( cl , "sample_size" , envir=environment())
    clusterExport( cl , "me.cv.range.factor" , envir=environment())
    clusterExport( cl , "true_p95" , envir=environment())
    clusterExport( cl , "true_gsd" , envir=environment())
    clusterExport( cl , "oel" , envir=environment())
    clusterExport( cl , "me.cv" , envir=environment())
    
    # calculations
    
    simulation_result_par <- parLapply(cl, X = as.list(1:n_simul) , fun = my.parallel.function) 
    
    # recommendation from the net: close the clusters
    stopCluster(cl)
    
    
    # estimation of computing time ( 9 min on my computer for 5000 iterations)
    end_time <- Sys.time()
    mytime <- end_time - start_time 
    
    
    
    ## interpretation of results 
    
    
    simulation_summary <- data.frame( ideal_p95 = numeric(n_simul),
                                      ideal_p95_ucl = numeric(n_simul),
                                      naive_p95 = numeric(n_simul),
                                      naive_p95_ucl = numeric(n_simul),
                                      me_p95 = numeric(n_simul),
                                      me_p95_ucl = numeric(n_simul),
                                      gum_p95 = numeric(n_simul),
                                      gum_p95_ucl = numeric(n_simul),
                                      ideal_p95_freq = numeric(n_simul),
                                      ideal_p95_ucl_freq = numeric(n_simul),
                                      naive_p95_freq = numeric(n_simul),
                                      naive_p95_ucl_freq = numeric(n_simul),
                                      
                                      ideal_gm = numeric(n_simul),
                                      naive_gm = numeric(n_simul),
                                      me_gm = numeric(n_simul),
                                      gum_gm = numeric(n_simul),
                                      ideal_gm_freq = numeric(n_simul),
                                      naive_gm_freq = numeric(n_simul),
                                      
                                      ideal_gsd = numeric(n_simul),
                                      naive_gsd = numeric(n_simul),
                                      me_gsd = numeric(n_simul),
                                      gum_gsd = numeric(n_simul),
                                      ideal_gsd_freq = numeric(n_simul),
                                      naive_gsd_freq = numeric(n_simul))
    for ( i in 1:n_simul) { 
      
      simulation_summary$ideal_p95[i] <- simulation_result_par[[i]][1]
      simulation_summary$ideal_p95_ucl[i] <- simulation_result_par[[i]][2]
      simulation_summary$naive_p95[i] <- simulation_result_par[[i]][3]
      simulation_summary$naive_p95_ucl[i] <- simulation_result_par[[i]][4]
      simulation_summary$me_p95[i] <- simulation_result_par[[i]][5]
      simulation_summary$me_p95_ucl[i] <- simulation_result_par[[i]][6]
      simulation_summary$gum_p95[i] <- simulation_result_par[[i]][7]
      simulation_summary$gum_p95_ucl[i] <- simulation_result_par[[i]][8]
      simulation_summary$ideal_p95_freq[i] <- simulation_result_par[[i]][9]
      simulation_summary$ideal_p95_ucl_freq[i] <- simulation_result_par[[i]][10]
      simulation_summary$naive_p95_freq[i] <- simulation_result_par[[i]][11]
      simulation_summary$naive_p95_ucl_freq[i] <- simulation_result_par[[i]][12]
      
      
      
      simulation_summary$ideal_gm[i] = simulation_result_par[[i]][13]
      simulation_summary$naive_gm[i] = simulation_result_par[[i]][14]
      simulation_summary$me_gm[i] = simulation_result_par[[i]][15]
      simulation_summary$gum_gm[i] = simulation_result_par[[i]][16]
      simulation_summary$ideal_gm_freq[i] = simulation_result_par[[i]][17]
      simulation_summary$naive_gm_freq[i] = simulation_result_par[[i]][18]
      
      
      simulation_summary$ideal_gsd[i] = simulation_result_par[[i]][19]
      simulation_summary$naive_gsd[i] = simulation_result_par[[i]][20]
      simulation_summary$me_gsd[i] = simulation_result_par[[i]][21]
      simulation_summary$gum_gsd[i] = simulation_result_par[[i]][22]
      simulation_summary$ideal_gsd_freq[i] = simulation_result_par[[i]][23]
      simulation_summary$naive_gsd_freq[i] = simulation_result_par[[i]][24]
      
    }
    
    
    ## saving simulation results
    
    #saveRDS( simulation_summary, "F:/Dropbox/temp/aioh2023-S2_3e.RDS")
    
    
    
    
    ###### simulation 4 ###### 
    
    ## impact on estimation of gm, gsd, p95, and the p95 70% UCL
    
    expanded_uncertainty <- 0.50
    
    coverage_factor <- qnorm(0.975)
    
    me.cv.range.factor <- 1 #( <1, if one wants to input uncertain uncertainty in the ME analysis)
    
    sample_size <- 3
    
    true_p95 <- 100
    
    true_gsd <- 1.5
    
    oel <- 300
    
    me.cv <- expanded_uncertainty/coverage_factor
    
    n_simul <- 10000 
    
    
    
    ##  parallel simulation (already performed and saved)
    
    # compteur de temps initialisé
    start_time <- Sys.time()
    
    
    #procédure parallele
    
    # creating the clusters ( use detectCores() to count the cores available, choose this number minus 2 or 4 below) 
    
    cl <- makeCluster(9)
    
    # libraries and scripts to be used in each cluster
    
    clusterEvalQ(cl, library(rjags))
    
    
    # sending objects to clusters
    clusterExport( cl , "jags.model.informedvar" , envir=environment())
    clusterExport( cl , "webexpo.seg.datapreparation" , envir=environment())
    clusterExport( cl , "Webexpo.seg.globalbayesian.jags" , envir=environment())
    clusterExport( cl , "fun.jags.informedvar" , envir=environment())
    clusterExport( cl , "myfunction.naive" , envir=environment())
    clusterExport( cl , "myfunction.me" , envir=environment())
    clusterExport( cl , "myfunction.me" , envir=environment())
    clusterExport( cl , "myfunction.gum" , envir=environment())
    clusterExport( cl , "myfunction.freq" , envir=environment())
    
    clusterExport( cl , "fun.NdExpo.lognorm" , envir=environment())
    clusterExport( cl , "fun.perc" , envir=environment())
    
    clusterExport( cl , "sample_size" , envir=environment())
    clusterExport( cl , "me.cv.range.factor" , envir=environment())
    clusterExport( cl , "true_p95" , envir=environment())
    clusterExport( cl , "true_gsd" , envir=environment())
    clusterExport( cl , "oel" , envir=environment())
    clusterExport( cl , "me.cv" , envir=environment())
    
    # calculations
    
    simulation_result_par <- parLapply(cl, X = as.list(1:n_simul) , fun = my.parallel.function) 
    
    # recommendation from the net: close the clusters
    stopCluster(cl)
    
    
    # estimation of computing time ( 9 min on my computer for 5000 iterations)
    end_time <- Sys.time()
    mytime <- end_time - start_time 
    
    
    
    ## interpretation of results 
    
    simulation_summary <- data.frame( ideal_p95 = numeric(n_simul),
                                      ideal_p95_ucl = numeric(n_simul),
                                      naive_p95 = numeric(n_simul),
                                      naive_p95_ucl = numeric(n_simul),
                                      me_p95 = numeric(n_simul),
                                      me_p95_ucl = numeric(n_simul),
                                      gum_p95 = numeric(n_simul),
                                      gum_p95_ucl = numeric(n_simul),
                                      ideal_p95_freq = numeric(n_simul),
                                      ideal_p95_ucl_freq = numeric(n_simul),
                                      naive_p95_freq = numeric(n_simul),
                                      naive_p95_ucl_freq = numeric(n_simul),
                                      
                                      ideal_gm = numeric(n_simul),
                                      naive_gm = numeric(n_simul),
                                      me_gm = numeric(n_simul),
                                      gum_gm = numeric(n_simul),
                                      ideal_gm_freq = numeric(n_simul),
                                      naive_gm_freq = numeric(n_simul),
                                      
                                      ideal_gsd = numeric(n_simul),
                                      naive_gsd = numeric(n_simul),
                                      me_gsd = numeric(n_simul),
                                      gum_gsd = numeric(n_simul),
                                      ideal_gsd_freq = numeric(n_simul),
                                      naive_gsd_freq = numeric(n_simul))
    for ( i in 1:n_simul) { 
      
      simulation_summary$ideal_p95[i] <- simulation_result_par[[i]][1]
      simulation_summary$ideal_p95_ucl[i] <- simulation_result_par[[i]][2]
      simulation_summary$naive_p95[i] <- simulation_result_par[[i]][3]
      simulation_summary$naive_p95_ucl[i] <- simulation_result_par[[i]][4]
      simulation_summary$me_p95[i] <- simulation_result_par[[i]][5]
      simulation_summary$me_p95_ucl[i] <- simulation_result_par[[i]][6]
      simulation_summary$gum_p95[i] <- simulation_result_par[[i]][7]
      simulation_summary$gum_p95_ucl[i] <- simulation_result_par[[i]][8]
      simulation_summary$ideal_p95_freq[i] <- simulation_result_par[[i]][9]
      simulation_summary$ideal_p95_ucl_freq[i] <- simulation_result_par[[i]][10]
      simulation_summary$naive_p95_freq[i] <- simulation_result_par[[i]][11]
      simulation_summary$naive_p95_ucl_freq[i] <- simulation_result_par[[i]][12]
      
      
      
      simulation_summary$ideal_gm[i] = simulation_result_par[[i]][13]
      simulation_summary$naive_gm[i] = simulation_result_par[[i]][14]
      simulation_summary$me_gm[i] = simulation_result_par[[i]][15]
      simulation_summary$gum_gm[i] = simulation_result_par[[i]][16]
      simulation_summary$ideal_gm_freq[i] = simulation_result_par[[i]][17]
      simulation_summary$naive_gm_freq[i] = simulation_result_par[[i]][18]
      
      
      simulation_summary$ideal_gsd[i] = simulation_result_par[[i]][19]
      simulation_summary$naive_gsd[i] = simulation_result_par[[i]][20]
      simulation_summary$me_gsd[i] = simulation_result_par[[i]][21]
      simulation_summary$gum_gsd[i] = simulation_result_par[[i]][22]
      simulation_summary$ideal_gsd_freq[i] = simulation_result_par[[i]][23]
      simulation_summary$naive_gsd_freq[i] = simulation_result_par[[i]][24]
      
    }
    
    
    ## saving simulation results
    
    #saveRDS( simulation_summary, "F:/Dropbox/temp/aioh2023-S2_4e.RDS")    
    
###### simulation 5 ###### 
    
    ## impact on estimation of gm, gsd, p95, and the p95 70% UCL
    
    expanded_uncertainty <- 0.50
    
    coverage_factor <- qnorm(0.975)
    
    me.cv.range.factor <- 1 #( <1, if one wants to input uncertain uncertainty in the ME analysis)
    
    sample_size <- 9
    
    true_p95 <- 100
    
    true_gsd <- 2.5
    
    oel <- 300
    
    me.cv <- expanded_uncertainty/coverage_factor
    
    n_simul <- 10000 
    
    
    
    ##  parallel simulation (already performed and saved)
    
    # compteur de temps initialisé
    start_time <- Sys.time()
    
    
    #procédure parallele
    
    # creating the clusters ( use detectCores() to count the cores available, choose this number minus 2 or 4 below) 
    
    cl <- makeCluster(9)
    
    # libraries and scripts to be used in each cluster
    
    clusterEvalQ(cl, library(rjags))
    
    
    # sending objects to clusters
    clusterExport( cl , "jags.model.informedvar" , envir=environment())
    clusterExport( cl , "webexpo.seg.datapreparation" , envir=environment())
    clusterExport( cl , "Webexpo.seg.globalbayesian.jags" , envir=environment())
    clusterExport( cl , "fun.jags.informedvar" , envir=environment())
    clusterExport( cl , "myfunction.naive" , envir=environment())
    clusterExport( cl , "myfunction.me" , envir=environment())
    clusterExport( cl , "myfunction.me" , envir=environment())
    clusterExport( cl , "myfunction.gum" , envir=environment())
    clusterExport( cl , "myfunction.freq" , envir=environment())
    
    clusterExport( cl , "fun.NdExpo.lognorm" , envir=environment())
    clusterExport( cl , "fun.perc" , envir=environment())
    
    clusterExport( cl , "sample_size" , envir=environment())
    clusterExport( cl , "me.cv.range.factor" , envir=environment())
    clusterExport( cl , "true_p95" , envir=environment())
    clusterExport( cl , "true_gsd" , envir=environment())
    clusterExport( cl , "oel" , envir=environment())
    clusterExport( cl , "me.cv" , envir=environment())
    
    # calculations
    
    simulation_result_par <- parLapply(cl, X = as.list(1:n_simul) , fun = my.parallel.function) 
    
    # recommendation from the net: close the clusters
    stopCluster(cl)
    
    
    # estimation of computing time ( 9 min on my computer for 5000 iterations)
    end_time <- Sys.time()
    mytime <- end_time - start_time 
    
    
    
    ## interpretation of results 
    
    simulation_summary <- data.frame( ideal_p95 = numeric(n_simul),
                                      ideal_p95_ucl = numeric(n_simul),
                                      naive_p95 = numeric(n_simul),
                                      naive_p95_ucl = numeric(n_simul),
                                      me_p95 = numeric(n_simul),
                                      me_p95_ucl = numeric(n_simul),
                                      gum_p95 = numeric(n_simul),
                                      gum_p95_ucl = numeric(n_simul),
                                      ideal_p95_freq = numeric(n_simul),
                                      ideal_p95_ucl_freq = numeric(n_simul),
                                      naive_p95_freq = numeric(n_simul),
                                      naive_p95_ucl_freq = numeric(n_simul),
                                      
                                      ideal_gm = numeric(n_simul),
                                      naive_gm = numeric(n_simul),
                                      me_gm = numeric(n_simul),
                                      gum_gm = numeric(n_simul),
                                      ideal_gm_freq = numeric(n_simul),
                                      naive_gm_freq = numeric(n_simul),
                                      
                                      ideal_gsd = numeric(n_simul),
                                      naive_gsd = numeric(n_simul),
                                      me_gsd = numeric(n_simul),
                                      gum_gsd = numeric(n_simul),
                                      ideal_gsd_freq = numeric(n_simul),
                                      naive_gsd_freq = numeric(n_simul))
    for ( i in 1:n_simul) { 
      
      simulation_summary$ideal_p95[i] <- simulation_result_par[[i]][1]
      simulation_summary$ideal_p95_ucl[i] <- simulation_result_par[[i]][2]
      simulation_summary$naive_p95[i] <- simulation_result_par[[i]][3]
      simulation_summary$naive_p95_ucl[i] <- simulation_result_par[[i]][4]
      simulation_summary$me_p95[i] <- simulation_result_par[[i]][5]
      simulation_summary$me_p95_ucl[i] <- simulation_result_par[[i]][6]
      simulation_summary$gum_p95[i] <- simulation_result_par[[i]][7]
      simulation_summary$gum_p95_ucl[i] <- simulation_result_par[[i]][8]
      simulation_summary$ideal_p95_freq[i] <- simulation_result_par[[i]][9]
      simulation_summary$ideal_p95_ucl_freq[i] <- simulation_result_par[[i]][10]
      simulation_summary$naive_p95_freq[i] <- simulation_result_par[[i]][11]
      simulation_summary$naive_p95_ucl_freq[i] <- simulation_result_par[[i]][12]
      
      
      
      simulation_summary$ideal_gm[i] = simulation_result_par[[i]][13]
      simulation_summary$naive_gm[i] = simulation_result_par[[i]][14]
      simulation_summary$me_gm[i] = simulation_result_par[[i]][15]
      simulation_summary$gum_gm[i] = simulation_result_par[[i]][16]
      simulation_summary$ideal_gm_freq[i] = simulation_result_par[[i]][17]
      simulation_summary$naive_gm_freq[i] = simulation_result_par[[i]][18]
      
      
      simulation_summary$ideal_gsd[i] = simulation_result_par[[i]][19]
      simulation_summary$naive_gsd[i] = simulation_result_par[[i]][20]
      simulation_summary$me_gsd[i] = simulation_result_par[[i]][21]
      simulation_summary$gum_gsd[i] = simulation_result_par[[i]][22]
      simulation_summary$ideal_gsd_freq[i] = simulation_result_par[[i]][23]
      simulation_summary$naive_gsd_freq[i] = simulation_result_par[[i]][24]
      
    }
    
    
    ## saving simulation results
    
    #saveRDS( simulation_summary, "F:/Dropbox/temp/aioh2023-S2_5e.RDS")    
    
###### simulation 6 ###### 
    
    ## impact on estimation of gm, gsd, p95, and the p95 70% UCL
    
    expanded_uncertainty <- 0.50
    
    coverage_factor <- qnorm(0.975)
    
    me.cv.range.factor <- 1 #( <1, if one wants to input uncertain uncertainty in the ME analysis)
    
    sample_size <- 9
    
    true_p95 <- 100
    
    true_gsd <- 1.5
    
    oel <- 300
    
    me.cv <- expanded_uncertainty/coverage_factor
    
    n_simul <- 10000 
    
    
    
    ##  parallel simulation (already performed and saved)
    
    # compteur de temps initialisé
    start_time <- Sys.time()
    
    
    #procédure parallele
    
    # creating the clusters ( use detectCores() to count the cores available, choose this number minus 2 or 4 below) 
    
    cl <- makeCluster(9)
    
    # libraries and scripts to be used in each cluster
    
    clusterEvalQ(cl, library(rjags))
    
    
    # sending objects to clusters
    clusterExport( cl , "jags.model.informedvar" , envir=environment())
    clusterExport( cl , "webexpo.seg.datapreparation" , envir=environment())
    clusterExport( cl , "Webexpo.seg.globalbayesian.jags" , envir=environment())
    clusterExport( cl , "fun.jags.informedvar" , envir=environment())
    clusterExport( cl , "myfunction.naive" , envir=environment())
    clusterExport( cl , "myfunction.me" , envir=environment())
    clusterExport( cl , "myfunction.me" , envir=environment())
    clusterExport( cl , "myfunction.gum" , envir=environment())
    clusterExport( cl , "myfunction.freq" , envir=environment())
    
    clusterExport( cl , "fun.NdExpo.lognorm" , envir=environment())
    clusterExport( cl , "fun.perc" , envir=environment())
    
    clusterExport( cl , "sample_size" , envir=environment())
    clusterExport( cl , "me.cv.range.factor" , envir=environment())
    clusterExport( cl , "true_p95" , envir=environment())
    clusterExport( cl , "true_gsd" , envir=environment())
    clusterExport( cl , "oel" , envir=environment())
    clusterExport( cl , "me.cv" , envir=environment())
    
    # calculations
    
    simulation_result_par <- parLapply(cl, X = as.list(1:n_simul) , fun = my.parallel.function) 
    
    # recommendation from the net: close the clusters
    stopCluster(cl)
    
    
    # estimation of computing time ( 9 min on my computer for 5000 iterations)
    end_time <- Sys.time()
    mytime <- end_time - start_time 
    
    
    
    ## interpretation of results 
    
    simulation_summary <- data.frame( ideal_p95 = numeric(n_simul),
                                      ideal_p95_ucl = numeric(n_simul),
                                      naive_p95 = numeric(n_simul),
                                      naive_p95_ucl = numeric(n_simul),
                                      me_p95 = numeric(n_simul),
                                      me_p95_ucl = numeric(n_simul),
                                      gum_p95 = numeric(n_simul),
                                      gum_p95_ucl = numeric(n_simul),
                                      ideal_p95_freq = numeric(n_simul),
                                      ideal_p95_ucl_freq = numeric(n_simul),
                                      naive_p95_freq = numeric(n_simul),
                                      naive_p95_ucl_freq = numeric(n_simul),
                                      
                                      ideal_gm = numeric(n_simul),
                                      naive_gm = numeric(n_simul),
                                      me_gm = numeric(n_simul),
                                      gum_gm = numeric(n_simul),
                                      ideal_gm_freq = numeric(n_simul),
                                      naive_gm_freq = numeric(n_simul),
                                      
                                      ideal_gsd = numeric(n_simul),
                                      naive_gsd = numeric(n_simul),
                                      me_gsd = numeric(n_simul),
                                      gum_gsd = numeric(n_simul),
                                      ideal_gsd_freq = numeric(n_simul),
                                      naive_gsd_freq = numeric(n_simul))
    for ( i in 1:n_simul) { 
      
      simulation_summary$ideal_p95[i] <- simulation_result_par[[i]][1]
      simulation_summary$ideal_p95_ucl[i] <- simulation_result_par[[i]][2]
      simulation_summary$naive_p95[i] <- simulation_result_par[[i]][3]
      simulation_summary$naive_p95_ucl[i] <- simulation_result_par[[i]][4]
      simulation_summary$me_p95[i] <- simulation_result_par[[i]][5]
      simulation_summary$me_p95_ucl[i] <- simulation_result_par[[i]][6]
      simulation_summary$gum_p95[i] <- simulation_result_par[[i]][7]
      simulation_summary$gum_p95_ucl[i] <- simulation_result_par[[i]][8]
      simulation_summary$ideal_p95_freq[i] <- simulation_result_par[[i]][9]
      simulation_summary$ideal_p95_ucl_freq[i] <- simulation_result_par[[i]][10]
      simulation_summary$naive_p95_freq[i] <- simulation_result_par[[i]][11]
      simulation_summary$naive_p95_ucl_freq[i] <- simulation_result_par[[i]][12]
      
      
      
      simulation_summary$ideal_gm[i] = simulation_result_par[[i]][13]
      simulation_summary$naive_gm[i] = simulation_result_par[[i]][14]
      simulation_summary$me_gm[i] = simulation_result_par[[i]][15]
      simulation_summary$gum_gm[i] = simulation_result_par[[i]][16]
      simulation_summary$ideal_gm_freq[i] = simulation_result_par[[i]][17]
      simulation_summary$naive_gm_freq[i] = simulation_result_par[[i]][18]
      
      
      simulation_summary$ideal_gsd[i] = simulation_result_par[[i]][19]
      simulation_summary$naive_gsd[i] = simulation_result_par[[i]][20]
      simulation_summary$me_gsd[i] = simulation_result_par[[i]][21]
      simulation_summary$gum_gsd[i] = simulation_result_par[[i]][22]
      simulation_summary$ideal_gsd_freq[i] = simulation_result_par[[i]][23]
      simulation_summary$naive_gsd_freq[i] = simulation_result_par[[i]][24]
      
    }
    
    
    ## saving simulation results
    
    #saveRDS( simulation_summary, "F:/Dropbox/temp/aioh2023-S2_6e.RDS")    
    
###### simulation 7 : correction of the UTL95/70 calculation for coverage results ######  
    
    
    # constants
    
    oel     <- 100L   # Default Occupational Exposure Limit (OEL).
    n_sim   <-  10000L   # Number of simulations to use.
    expanded_uncertainty <- 0.50
    coverage_factor <- qnorm(0.975)
    me.cv.range.factor <- 1 #( <1, if one wants to input uncertain uncertainty in the ME analysis)
    me.cv <- expanded_uncertainty/coverage_factor
    true_p95 <- 100
    
    # Scenarios.
   
         # A scenario is a unique combination of constants. We consider 3025 scenarios.
         
         scenarios <- data.frame( sample_size = c(6,3,6,3,9,9),
                                  true_gsd = c(2.5,2.5,1.5,1.5,2.5,1.5))
        
        
         scenarios$true_gm <- exp( log(true_p95) - qnorm(0.95)*log(scenarios$true_gsd) )
     
    # testing the various steps    
           
        # calls to functions in data_simulation.R   
        
        
        simulated_data_list <- data.simulation( true_gm = scenarios$true_gm[1]  , 
                                                     true_gsd = scenarios$true_gsd[1]  , 
                                                     n_sim = n_sim , 
                                                     sample_size = scenarios$sample_size[1] ) 
        
        
    
    
        # Analysis of the generated samples 
    
        
    
         test <-  myfunction.gum( simulated_data_list$true[,1] , n.simul=10000 , me.cv=me.cv )$p95ucl70 > true_p95
    
    
        # Looping across the n_sim iterations for one scenario and one method out of 3 : ideal naive me
          
          one_scenario_ideal <- apply( X = simulated_data_list$true ,
                                      MARGIN = 2,
                                      FUN = function(x){ myfunction.freq( x )$p95ucl70 > true_p95} ,
                                      simplify = TRUE) 
          
          one_scenario_naive <- apply( X = simulated_data_list$observed ,
                                       MARGIN = 2,
                                       FUN = function(x){ myfunction.freq( x )$p95ucl70 > true_p95} ,
                                       simplify = TRUE) 
          
          one_scenario_me <- apply( X = simulated_data_list$true ,
                                       MARGIN = 2,
                                       FUN = function(x){ myfunction.gum( x , n.simul=10000 , me.cv=me.cv )$p95ucl70 > true_p95} ,
                                       simplify = TRUE) 
                                        
 
          mean(one_scenario_ideal)
          mean(one_scenario_naive)
          mean(one_scenario_me)
          
          

     # function  for the lapply loop
          
          myfunction <- function( x ) {
            
            simulated_data_list <- data.simulation( true_gm = x$true_gm[1]  , 
                                                         true_gsd = x$true_gsd[1]  , 
                                                         n_sim = n_sim , 
                                                         sample_size = x$sample_size[1] )
            
            one_scenario_ideal <- apply( X = simulated_data_list$true ,
                                         MARGIN = 2,
                                         FUN = function(x){ myfunction.freq( x )$p95ucl70 > true_p95} ,
                                         simplify = TRUE) 
            
            one_scenario_naive <- apply( X = simulated_data_list$observed ,
                                         MARGIN = 2,
                                         FUN = function(x){ myfunction.freq( x )$p95ucl70 > true_p95} ,
                                         simplify = TRUE) 
            
            one_scenario_me <- apply( X = simulated_data_list$true ,
                                      MARGIN = 2,
                                      FUN = function(x){ myfunction.gum( x , n.simul=10000 , me.cv=me.cv )$p95ucl70 > true_p95} ,
                                      simplify = TRUE) 
            
            return( list( ideal = mean(one_scenario_ideal),
                          naive = mean(one_scenario_naive),
                          me = mean(one_scenario_me)) ) }      
          
          
          
    # calculation loop #3.2 hours
    
    full_scenario_list <- scenarios
    
    start_time <- Sys.time()
    
      list_of_scenarios <- vector(  mode = "list" , length = length(scenarios[,1]))
      
      for (j in 1:length(scenarios[,1])) list_of_scenarios[[j]] <-scenarios[j,]  
      
      
      # creating the clusters  
      cl <- makeCluster(6)
      
      # libraries and scripts
      
      # sending objects to clusters
 
      clusterExport( cl , "myfunction.gum" , envir=environment())
      clusterExport( cl , "myfunction.freq" , envir=environment())
      clusterExport( cl , "fun.NdExpo.lognorm" , envir=environment())
      clusterExport( cl , "fun.perc.en689" , envir=environment())
      clusterExport( cl , "fun.perc" , envir=environment())
      clusterExport( cl , "data.simulation" , envir=environment())
      clusterExport( cl , "K.factor" , envir=environment())
      
      
      clusterExport( cl , "true_p95" , envir=environment())
      clusterExport( cl , "me.cv" , envir=environment())
      clusterExport( cl , "n_sim" , envir=environment())
      
      
      simulation_result <- parLapply(cl, X = list_of_scenarios , fun = myfunction) 
      
      # recommendation from the net
      stopCluster(cl)
    
      
    end_time <- Sys.time()
    mytime <- end_time - start_time 
    
    
    saveRDS( list(sim=simulation_result,
                  param=scenarios), "created data/GUM measurement error_sim7.RDS")    
    
    