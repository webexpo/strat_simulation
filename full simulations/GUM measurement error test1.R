#' exploring the impact of measurement error industrial hygiene measurement data interpretation

#" question : for one example of true distribution, what is the impact of a 50% expanded measurement error (CV=50/1.96) on P95 and its UCL


##### Data #####

    ## libraries  
    
    library(ggplot2)
    library(ggthemes)
    library(writexl)
    library(scales)
    library(readxl)
    library(rjags)
    library(parallel)
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
 
    ### smal function for summary results for a simple uninformative model
    
    
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
      
      # mysample <- c("28.9","19.40","<8.22","149.9","26.42","56.1")
      # n.simul <- 10000
      
      # ROS
      
      mysample.ros <- fun.NdExpo.lognorm( mysample )$data$xfin
      
      
      little.ros.gum <- function(i) {
      
      observed_data <- pmax( 0.1 , rnorm( length(mysample.ros) , mysample.ros , mysample.ros*me.cv ) ) #truncation at 0.1
        
      gm.ros <- exp(mean(log(observed_data)))
      
      gsd.ros <- exp(sd(log(observed_data)))
      
      p95.ros <- fun.perc(observed_data,alpha=0.05,perc=0.95)$est
      
      p95ucl70.ros <- fun.perc(observed_data,alpha=0.30,perc=0.95)$uc
      
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
      
      # mysample <- c("28.9","19.40","<8.22","149.9","26.42","56.1")
      # n.simul <- 10000
      
      # ROS
      
      mysample.ros <- fun.NdExpo.lognorm( mysample )$data$xfin
      
      observed_data <- pmax( 0.1 , rnorm( length(mysample.ros) , mysample.ros , mysample.ros*me.cv ) ) #truncation at 0.1
        
        gm.ros <- exp(mean(log(mysample.ros)))
        
        gsd.ros <- exp(sd(log(mysample.ros)))
        
        p95.ros <- fun.perc(mysample.ros,alpha=0.05,perc=0.95)$est
        
        p95ucl70.ros <- fun.perc(mysample.ros,alpha=0.30,perc=0.95)$uc
        
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
      
      observed_data <- pmax( 0.1 , rnorm( sample_size , true_data , true_data*me.cv ) ) #truncation at 0.1
      
      ## Bayesian analysis
      
      ideal_analysis <- myfunction.naive( true_data , oel)
      
      naive_analysis <- myfunction.naive( observed_data , oel)
      
      me_analysis <- myfunction.me( observed_data , me.range = c(me.cv.range.factor*me.cv,me.cv/me.cv.range.factor) ,oel)
      
      gum_analysis <- myfunction.gum( observed_data , 10000 , me.cv)
      
      ## data interpretation
      
      results <- c( median(ideal_analysis$p95_chain),quantile(ideal_analysis$p95_chain,0.7),
                    median(naive_analysis$p95_chain),quantile(naive_analysis$p95_chain,0.7),
                    median(me_analysis$p95_chain),quantile(me_analysis$p95_chain,0.7),
                    gum_analysis$p95,gum_analysis$p95ucl70,
                    
                    median(ideal_analysis$gm_chain),
                    median(naive_analysis$gm_chain),
                    median(me_analysis$gm_chain),
                    gum_analysis$gm,
                    
                    median(ideal_analysis$gsd_chain),
                    median(naive_analysis$gsd_chain),
                    median(me_analysis$gsd_chain),
                    gum_analysis$gsd)
      
      names(results) <- c("ideal_p95","ideal_p95UCL70","naive_p95","naive_p95UCL70","me_p95","me_p95UCL70","gum_p95","gum_p95UCL70",
                          "ideal_gm","naive_gm","me_gm","gum_gm",
                          "ideal_gsd","naive_gsd","me_gsd","gum_gsd")
      
      return( results) }
    
#### Analysis #####
    
  ###### Example 1 #####
    
          ## parameters
          
          expanded_uncertainty <- 0.50
          
          coverage_factor <- qnorm(0.975)
          
          me.cv.range.factor <- 1 #( <1, if one wants to input uncertain uncertainty in the ME analysis)
          
          sample_size <- 6
          
          true_p95 <- 100
          
          true_gsd <- 2.5
          
          oel <- 300
          
          true_gm <- exp( log(true_p95) - qnorm(0.95)*log(true_gsd) )
          
          ## data generation
          
          true_data <- exp( rnorm( sample_size , log(true_gm) , log(true_gsd) ) )   
          
          observed_data <- pmax( 0.1 , rnorm( sample_size , true_data , true_data*expanded_uncertainty/coverage_factor ) ) #truncation at 0.1
          
          ## Bayesian analysis
          
          ideal_analysis <- myfunction.naive( true_data , oel)
          
          naive_analysis <- myfunction.naive( observed_data , oel)
          
          me_analysis <- myfunction.me( observed_data , me.range = c(me.cv.range.factor*expanded_uncertainty/coverage_factor,
                                                                     expanded_uncertainty/(coverage_factor*me.cv.range.factor)) ,oel)
          
          gum_analysis <- myfunction.gum( observed_data , 10000 , expanded_uncertainty/coverage_factor)
          
          ##numerical results - P95 estimates
          
          exemple1_table <- data.frame( type = c("true","ideal","naive","me","gum"),
            
                                        gm = c( true_gm,
                                                median(ideal_analysis$gm_chain),
                                                median(naive_analysis$gm_chain),
                                                median(me_analysis$gm_chain),
                                                gum_analysis$gm),
                                        
                                        gsd = c( true_gsd,
                                                 median(ideal_analysis$gsd_chain),
                                                 median(naive_analysis$gsd_chain),
                                                 median(me_analysis$gsd_chain),
                                                 gum_analysis$gsd),
                                        
                                        p95 = c( true_p95,
                                                 quantile(ideal_analysis$p95_chain,0.5),
                                                 quantile(naive_analysis$p95_chain,0.5),
                                                 quantile(me_analysis$p95_chain,0.5),
                                                 gum_analysis$p95),
                                        
                                        p95_ucl = c( NA,
                                                    quantile(ideal_analysis$p95_chain,0.7),
                                                    quantile(naive_analysis$p95_chain,0.7),
                                                    quantile(me_analysis$p95_chain,0.7),
                                                    gum_analysis$p95ucl70)
                                       
                                       )
          
          
          exemple1.object <- list( true_data = true_data,
                                   observed_data = observed_data,
                                   results = exemple1_table)
          
          saveRDS( exemple1.object, "created data/aioh2023-t1.RDS")
          
          ## plot
      
          mylength <- 1000
          
          myindex <- sample( x=1:25000 , size=mylength , replace = FALSE)
          
          
          mcmc.data <- data.frame( p95 = c( ideal_analysis$p95_chain[myindex],
                                            naive_analysis$p95_chain[myindex],
                                            me_analysis$p95_chain[myindex]),
          
                                   type = c( rep("p95_ideal", mylength), rep("p95_naive", mylength) , rep("p95_me", mylength)))
                                   
          mcmc.data$index <- sample( x = 1:length(mcmc.data[,1]) , size = mylength*3 , replace = FALSE)
          
           
          p <- ggplot(data = mcmc.data, aes(x=index,y=p95,color = type) ) 
          
          p <- p + geom_point(alpha = 0.3, size = 3  )
          
          p <- p + theme_calc()
          
          
          p <- p +  theme(axis.text.x=element_text(size=14),axis.title.x=element_text(size=16,vjust=+0.55))+
            theme(axis.text.y=element_text(size=14),axis.title.y=element_text(size=16,vjust=+0.55))+
            labs(x=expression ( x=""), y=expression ( x="Plausible values for P95"))
          
          p <- p + scale_y_continuous(expand = c(0, 0) , 
                                      limits = c(1, quantile(mcmc.data$p95,0.95)),
                                      labels = scales::number_format(accuracy = 0.1,
                                                                     decimal.mark = ',') )
          p <- p + scale_x_continuous(expand = c(0, 20) , 
                                      limits = c(0, mylength*3) )
          
          p <- p + theme( plot.background = element_rect(color = NA) )
          
          
          p <- p + geom_hline(yintercept=quantile(ideal_analysis$p95_chain,0.7) , color = "red" , linetype="dashed",size = 1.5)
          p <- p + geom_hline(yintercept=quantile(naive_analysis$p95_chain,0.7) , color = "blue" , linetype="dashed",size = 1.5)
          p <- p + geom_hline(yintercept=quantile(me_analysis$p95_chain,0.7) , color = "green" , linetype="dashed",size = 1.5)
          
          
          p

    
  ###### simulation 1 R1 ##### 
    
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
        
        cl <- makeCluster(17)
        
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
                                                
                                                ideal_gm = numeric(n_simul),
                                                naive_gm = numeric(n_simul),
                                                me_gm = numeric(n_simul),
                                                gum_gm = numeric(n_simul),
                                                  
                                                ideal_gsd = numeric(n_simul),
                                                naive_gsd = numeric(n_simul),
                                                me_gsd = numeric(n_simul),
                                                gum_gsd = numeric(n_simul))
        for ( i in 1:n_simul) { 
          
          simulation_summary$ideal_p95[i] <- simulation_result_par[[i]][1]
          simulation_summary$ideal_p95_ucl[i] <- simulation_result_par[[i]][2]
          simulation_summary$naive_p95[i] <- simulation_result_par[[i]][3]
          simulation_summary$naive_p95_ucl[i] <- simulation_result_par[[i]][4]
          simulation_summary$me_p95[i] <- simulation_result_par[[i]][5]
          simulation_summary$me_p95_ucl[i] <- simulation_result_par[[i]][6]
          simulation_summary$gum_p95[i] <- simulation_result_par[[i]][7]
          simulation_summary$gum_p95_ucl[i] <- simulation_result_par[[i]][8]
         
           simulation_summary$ideal_gm[i] = simulation_result_par[[i]][9]
          simulation_summary$naive_gm[i] = simulation_result_par[[i]][10]
          simulation_summary$me_gm[i] = simulation_result_par[[i]][11]
          simulation_summary$gum_gm[i] = simulation_result_par[[i]][12]
          
          simulation_summary$ideal_gsd[i] = simulation_result_par[[i]][13]
          simulation_summary$naive_gsd[i] = simulation_result_par[[i]][14]
          simulation_summary$me_gsd[i] = simulation_result_par[[i]][15]
          simulation_summary$gum_gsd[i] = simulation_result_par[[i]][16]
          
        }

    
    ## saving simulation results
    
    saveRDS( simulation_summary, "created data/aioh2023-r1.RDS")
    
    
###### simulation 1 R2 ##### 
    
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
    

    ##  parallel simulation (already performed and saved)
    
    # compteur de temps initialisé
    start_time <- Sys.time()
    
    
    #procédure parallele
    
    # creating the clusters ( use detectCores() to count the cores available, choose this number minus 2 or 4 below) 
    
    cl <- makeCluster(10)
    
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
                                      
                                      ideal_gm = numeric(n_simul),
                                      naive_gm = numeric(n_simul),
                                      me_gm = numeric(n_simul),
                                      gum_gm = numeric(n_simul),
                                      
                                      ideal_gsd = numeric(n_simul),
                                      naive_gsd = numeric(n_simul),
                                      me_gsd = numeric(n_simul),
                                      gum_gsd = numeric(n_simul))
    for ( i in 1:n_simul) { 
      
      simulation_summary$ideal_p95[i] <- simulation_result_par[[i]][1]
      simulation_summary$ideal_p95_ucl[i] <- simulation_result_par[[i]][2]
      simulation_summary$naive_p95[i] <- simulation_result_par[[i]][3]
      simulation_summary$naive_p95_ucl[i] <- simulation_result_par[[i]][4]
      simulation_summary$me_p95[i] <- simulation_result_par[[i]][5]
      simulation_summary$me_p95_ucl[i] <- simulation_result_par[[i]][6]
      simulation_summary$gum_p95[i] <- simulation_result_par[[i]][7]
      simulation_summary$gum_p95_ucl[i] <- simulation_result_par[[i]][8]
      
      simulation_summary$ideal_gm[i] = simulation_result_par[[i]][9]
      simulation_summary$naive_gm[i] = simulation_result_par[[i]][10]
      simulation_summary$me_gm[i] = simulation_result_par[[i]][11]
      simulation_summary$gum_gm[i] = simulation_result_par[[i]][12]
      
      simulation_summary$ideal_gsd[i] = simulation_result_par[[i]][13]
      simulation_summary$naive_gsd[i] = simulation_result_par[[i]][14]
      simulation_summary$me_gsd[i] = simulation_result_par[[i]][15]
      simulation_summary$gum_gsd[i] = simulation_result_par[[i]][16]
      
    }
    
    
    ## saving simulation results
    
    saveRDS( simulation_summary, "created data/aioh2023-r1_2.RDS")    
    
    
    
    
###### simulation 1 R3 ##### 
    
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
    
    
    ##  parallel simulation (already performed and saved)
    
    # compteur de temps initialisé
    start_time <- Sys.time()
    
    
    #procédure parallele
    
    # creating the clusters ( use detectCores() to count the cores available, choose this number minus 2 or 4 below) 
    
    cl <- makeCluster(10)
    
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
                                      
                                      ideal_gm = numeric(n_simul),
                                      naive_gm = numeric(n_simul),
                                      me_gm = numeric(n_simul),
                                      gum_gm = numeric(n_simul),
                                      
                                      ideal_gsd = numeric(n_simul),
                                      naive_gsd = numeric(n_simul),
                                      me_gsd = numeric(n_simul),
                                      gum_gsd = numeric(n_simul))
    for ( i in 1:n_simul) { 
      
      simulation_summary$ideal_p95[i] <- simulation_result_par[[i]][1]
      simulation_summary$ideal_p95_ucl[i] <- simulation_result_par[[i]][2]
      simulation_summary$naive_p95[i] <- simulation_result_par[[i]][3]
      simulation_summary$naive_p95_ucl[i] <- simulation_result_par[[i]][4]
      simulation_summary$me_p95[i] <- simulation_result_par[[i]][5]
      simulation_summary$me_p95_ucl[i] <- simulation_result_par[[i]][6]
      simulation_summary$gum_p95[i] <- simulation_result_par[[i]][7]
      simulation_summary$gum_p95_ucl[i] <- simulation_result_par[[i]][8]
      
      simulation_summary$ideal_gm[i] = simulation_result_par[[i]][9]
      simulation_summary$naive_gm[i] = simulation_result_par[[i]][10]
      simulation_summary$me_gm[i] = simulation_result_par[[i]][11]
      simulation_summary$gum_gm[i] = simulation_result_par[[i]][12]
      
      simulation_summary$ideal_gsd[i] = simulation_result_par[[i]][13]
      simulation_summary$naive_gsd[i] = simulation_result_par[[i]][14]
      simulation_summary$me_gsd[i] = simulation_result_par[[i]][15]
      simulation_summary$gum_gsd[i] = simulation_result_par[[i]][16]
      
    }
    
    
    ## saving simulation results
    
    saveRDS( simulation_summary, "created data/aioh2023-r1_3.RDS")    
    
    
    

###### simulation 2 R1 ###### 
    
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
    
    cl <- makeCluster(10)
    
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
                                      
                                      ideal_gm = numeric(n_simul),
                                      naive_gm = numeric(n_simul),
                                      me_gm = numeric(n_simul),
                                      gum_gm = numeric(n_simul),
                                      
                                      ideal_gsd = numeric(n_simul),
                                      naive_gsd = numeric(n_simul),
                                      me_gsd = numeric(n_simul),
                                      gum_gsd = numeric(n_simul))
    for ( i in 1:n_simul) { 
      
      simulation_summary$ideal_p95[i] <- simulation_result_par[[i]][1]
      simulation_summary$ideal_p95_ucl[i] <- simulation_result_par[[i]][2]
      simulation_summary$naive_p95[i] <- simulation_result_par[[i]][3]
      simulation_summary$naive_p95_ucl[i] <- simulation_result_par[[i]][4]
      simulation_summary$me_p95[i] <- simulation_result_par[[i]][5]
      simulation_summary$me_p95_ucl[i] <- simulation_result_par[[i]][6]
      simulation_summary$gum_p95[i] <- simulation_result_par[[i]][7]
      simulation_summary$gum_p95_ucl[i] <- simulation_result_par[[i]][8]
      
      simulation_summary$ideal_gm[i] = simulation_result_par[[i]][9]
      simulation_summary$naive_gm[i] = simulation_result_par[[i]][10]
      simulation_summary$me_gm[i] = simulation_result_par[[i]][11]
      simulation_summary$gum_gm[i] = simulation_result_par[[i]][12]
      
      simulation_summary$ideal_gsd[i] = simulation_result_par[[i]][13]
      simulation_summary$naive_gsd[i] = simulation_result_par[[i]][14]
      simulation_summary$me_gsd[i] = simulation_result_par[[i]][15]
      simulation_summary$gum_gsd[i] = simulation_result_par[[i]][16]
      
    }
    
    
    ## saving simulation results
    
    saveRDS( simulation_summary, "created data/aioh2023-r2_1.RDS")    
    
    
###### simulation 2 R2 ###### 
    
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
    
    cl <- makeCluster(10)
    
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
                                      
                                      ideal_gm = numeric(n_simul),
                                      naive_gm = numeric(n_simul),
                                      me_gm = numeric(n_simul),
                                      gum_gm = numeric(n_simul),
                                      
                                      ideal_gsd = numeric(n_simul),
                                      naive_gsd = numeric(n_simul),
                                      me_gsd = numeric(n_simul),
                                      gum_gsd = numeric(n_simul))
    for ( i in 1:n_simul) { 
      
      simulation_summary$ideal_p95[i] <- simulation_result_par[[i]][1]
      simulation_summary$ideal_p95_ucl[i] <- simulation_result_par[[i]][2]
      simulation_summary$naive_p95[i] <- simulation_result_par[[i]][3]
      simulation_summary$naive_p95_ucl[i] <- simulation_result_par[[i]][4]
      simulation_summary$me_p95[i] <- simulation_result_par[[i]][5]
      simulation_summary$me_p95_ucl[i] <- simulation_result_par[[i]][6]
      simulation_summary$gum_p95[i] <- simulation_result_par[[i]][7]
      simulation_summary$gum_p95_ucl[i] <- simulation_result_par[[i]][8]
      
      simulation_summary$ideal_gm[i] = simulation_result_par[[i]][9]
      simulation_summary$naive_gm[i] = simulation_result_par[[i]][10]
      simulation_summary$me_gm[i] = simulation_result_par[[i]][11]
      simulation_summary$gum_gm[i] = simulation_result_par[[i]][12]
      
      simulation_summary$ideal_gsd[i] = simulation_result_par[[i]][13]
      simulation_summary$naive_gsd[i] = simulation_result_par[[i]][14]
      simulation_summary$me_gsd[i] = simulation_result_par[[i]][15]
      simulation_summary$gum_gsd[i] = simulation_result_par[[i]][16]
      
    }
    
    
    ## saving simulation results
    
    saveRDS( simulation_summary, "created data/aioh2023-r2_2.RDS")    
  
    
###### simulation 2 R3 ###### 
    
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
    
    cl <- makeCluster(10)
    
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
                                      
                                      ideal_gm = numeric(n_simul),
                                      naive_gm = numeric(n_simul),
                                      me_gm = numeric(n_simul),
                                      gum_gm = numeric(n_simul),
                                      
                                      ideal_gsd = numeric(n_simul),
                                      naive_gsd = numeric(n_simul),
                                      me_gsd = numeric(n_simul),
                                      gum_gsd = numeric(n_simul))
    for ( i in 1:n_simul) { 
      
      simulation_summary$ideal_p95[i] <- simulation_result_par[[i]][1]
      simulation_summary$ideal_p95_ucl[i] <- simulation_result_par[[i]][2]
      simulation_summary$naive_p95[i] <- simulation_result_par[[i]][3]
      simulation_summary$naive_p95_ucl[i] <- simulation_result_par[[i]][4]
      simulation_summary$me_p95[i] <- simulation_result_par[[i]][5]
      simulation_summary$me_p95_ucl[i] <- simulation_result_par[[i]][6]
      simulation_summary$gum_p95[i] <- simulation_result_par[[i]][7]
      simulation_summary$gum_p95_ucl[i] <- simulation_result_par[[i]][8]
      
      simulation_summary$ideal_gm[i] = simulation_result_par[[i]][9]
      simulation_summary$naive_gm[i] = simulation_result_par[[i]][10]
      simulation_summary$me_gm[i] = simulation_result_par[[i]][11]
      simulation_summary$gum_gm[i] = simulation_result_par[[i]][12]
      
      simulation_summary$ideal_gsd[i] = simulation_result_par[[i]][13]
      simulation_summary$naive_gsd[i] = simulation_result_par[[i]][14]
      simulation_summary$me_gsd[i] = simulation_result_par[[i]][15]
      simulation_summary$gum_gsd[i] = simulation_result_par[[i]][16]
      
    }
    
    
    ## saving simulation results
    
    saveRDS( simulation_summary, "created data/aioh2023-r2_3.RDS")    
    
    
