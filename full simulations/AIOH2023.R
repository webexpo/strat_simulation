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
    
    #' smal function for generating MCMC results for an analysis assuming CV error within a range
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
      
      me_analysis <- myfunction.me( observed_data , me.range = c(me.cv,me.cv) ,oel)
      
      ## data interpretation
      
      results <- c( median(ideal_analysis$p95_chain),quantile(ideal_analysis$p95_chain,0.7),
                    median(naive_analysis$p95_chain),quantile(naive_analysis$p95_chain,0.7),
                    median(me_analysis$p95_chain),quantile(me_analysis$p95_chain,0.7),
                    
                    median(ideal_analysis$gm_chain),
                    median(naive_analysis$gm_chain),
                    median(me_analysis$gm_chain),
                    
                    median(ideal_analysis$gsd_chain),
                    median(naive_analysis$gsd_chain),
                    median(me_analysis$gsd_chain))
      
      names(results) <- c("ideal_p95","ideal_p95UCL70","naive_p95","naive_p95UCL70","me_p95","me_p95UCL70",
                          "ideal_gm","naive_gm","me_gm",
                          "ideal_gsd","naive_gsd","me_gsd")
      
      return( results) }
    
#### Simulation ####
    
  ##### Example 1 ##### 
    
    ## parameters
    
    expanded_uncertainty <- 0.50
    
    coverage_factor <- qnorm(0.975)
    
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
    
    me_analysis <- myfunction.me( observed_data , me.range = c(0.5*expanded_uncertainty/coverage_factor,expanded_uncertainty/coverage_factor) ,oel)
    
    

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

    
##### Example 2 ##### 
    
    ## impact on estimation of gm, gsd, p95, and the p95 70% UCL
    
    expanded_uncertainty <- 0.50
    
    coverage_factor <- qnorm(0.975)
    
    sample_size <- 6
    
    true_p95 <- 100
    
    true_gsd <- 2.5
    
    oel <- 300
    
    me.cv <- 0.4
   
    n_simul <- 5000 
    
    ##  traditional simulation
    
        # compteur de temps initialisé
    start_time <- Sys.time()
    
    simulation_result <- lapply( X = as.list(1:n_simul) , FUN = my.parallel.function)
    
        # estimation of computing time ( 2 min on my computer for n=1000)
    end_time <- Sys.time()
    mytime <- end_time - start_time 
    
    
    
    ##  parallel simulation
    
        # compteur de temps initialisé
        start_time <- Sys.time()
        
        
        #procédure parallele
        
        # creating the clusters ( use detectCores() to count the cores available, choose this number minus 2 or 4 below) 
        
        cl <- makeCluster(8)
        
        # libraries and scripts to be used in each cluster
        
        clusterEvalQ(cl, library(rjags))
        
        
        # sending objects to clusters
        clusterExport( cl , "jags.model.informedvar" , envir=environment())
        clusterExport( cl , "webexpo.seg.datapreparation" , envir=environment())
        clusterExport( cl , "Webexpo.seg.globalbayesian.jags" , envir=environment())
        clusterExport( cl , "fun.jags.informedvar" , envir=environment())
        clusterExport( cl , "myfunction.naive" , envir=environment())
        clusterExport( cl , "myfunction.me" , envir=environment())
        
        
        
        clusterExport( cl , "sample_size" , envir=environment())
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
                                                
                                                ideal_gm = numeric(n_simul),
                                                naive_gm = numeric(n_simul),
                                                me_gm = numeric(n_simul),
                                                  
                                                ideal_gsd = numeric(n_simul),
                                                naive_gsd = numeric(n_simul),
                                                me_gsd = numeric(n_simul))
        for ( i in 1:n_simul) { 
          
          simulation_summary$ideal_p95[i] <- simulation_result_par[[i]][1]
          simulation_summary$ideal_p95_ucl[i] <- simulation_result_par[[i]][2]
          simulation_summary$naive_p95[i] <- simulation_result_par[[i]][3]
          simulation_summary$naive_p95_ucl[i] <- simulation_result_par[[i]][4]
          simulation_summary$me_p95[i] <- simulation_result_par[[i]][5]
          simulation_summary$me_p95_ucl[i] <- simulation_result_par[[i]][6]
         
           simulation_summary$ideal_gm[i] = simulation_result_par[[i]][7]
          simulation_summary$naive_gm[i] = simulation_result_par[[i]][8]
          simulation_summary$me_gm[i] = simulation_result_par[[i]][9]
          
          simulation_summary$ideal_gsd[i] = simulation_result_par[[i]][10]
          simulation_summary$naive_gsd[i] = simulation_result_par[[i]][11]
          simulation_summary$me_gsd[i] = simulation_result_par[[i]][12]
        }

    
  ###### visuals 1 : estimates ####
        
        mcmc.data <- data.frame( p95 = c( simulation_summary$ideal_p95_ucl,
                                          simulation_summary$naive_p95_ucl,
                                          simulation_summary$me_p95_ucl),
                                 
                                 type = c( rep("p95_ideal", mylength), rep("p95_naive", mylength) , rep("p95_me", mylength)))
        
        mcmc.data$index <- sample( x = 1:length(mcmc.data[,1]) , size = n_simul*3 , replace = FALSE)
        
        
        p <- ggplot(data = mcmc.data, aes(x=index,y=p95,color = type) ) 
        
        p <- p + geom_point(alpha = 0.3, size = 3  )
        
        p <- p + theme_calc()
        
        
        p <- p +  theme(axis.text.x=element_text(size=14),axis.title.x=element_text(size=16,vjust=+0.55))+
          theme(axis.text.y=element_text(size=14),axis.title.y=element_text(size=16,vjust=+0.55))+
          labs(x=expression ( x=""), y=expression ( x="Observed values for P95"))
        
        p <- p + scale_y_continuous(expand = c(0, 0) , 
                                    limits = c(1, quantile(mcmc.data$p95,0.8)),
                                    labels = scales::number_format(accuracy = 0.1,
                                                                   decimal.mark = ',') )
        p <- p + scale_x_continuous(expand = c(0, 20) , 
                                    limits = c(0, n_simul*3) )
        
        p <- p + theme( plot.background = element_rect(color = NA) )
        
        
        p <- p + geom_hline(yintercept=quantile(simulation_summary$ideal_p95,0.5) , color = "red" , linetype="dashed",size = 1.5)
        p <- p + geom_hline(yintercept=quantile(simulation_summary$naive_p95,0.5) , color = "blue" , linetype="dashed",size = 1.5)
        p <- p + geom_hline(yintercept=quantile(simulation_summary$me_p95,0.5) , color = "green" , linetype="dashed",size = 1.5)
        
        
        p  
        
        
        ###### visuals 2 : ucl ####
        
        mcmc.data <- data.frame( p95 = c( simulation_summary$ideal_p95_ucl,
                                          simulation_summary$naive_p95_ucl,
                                          simulation_summary$me_p95_ucl),
                                 
                                 type = c( rep("p95_ideal", mylength), rep("p95_naive", mylength) , rep("p95_me", mylength)))
        
        mcmc.data$index <- sample( x = 1:length(mcmc.data[,1]) , size = n_simul*3 , replace = FALSE)
        
        
        p <- ggplot(data = mcmc.data, aes(x=index,y=p95,color = type) ) 
        
        p <- p + geom_point(alpha = 0.3, size = 3  )
        
        p <- p + theme_calc()
        
        
        p <- p +  theme(axis.text.x=element_text(size=14),axis.title.x=element_text(size=16,vjust=+0.55))+
          theme(axis.text.y=element_text(size=14),axis.title.y=element_text(size=16,vjust=+0.55))+
          labs(x=expression ( x=""), y=expression ( x="Observed values for P95"))
        
        p <- p + scale_y_continuous(expand = c(0, 0) , 
                                    limits = c(1, quantile(mcmc.data$p95,0.8)),
                                    labels = scales::number_format(accuracy = 0.1,
                                                                   decimal.mark = ',') )
        p <- p + scale_x_continuous(expand = c(0, 20) , 
                                    limits = c(0, n_simul*3) )
        
        p <- p + theme( plot.background = element_rect(color = NA) )
        
        
        p <- p + geom_hline(yintercept=quantile(simulation_summary$ideal_p95_ucl,0.5) , color = "red" , linetype="dashed",size = 1.5)
        p <- p + geom_hline(yintercept=quantile(simulation_summary$naive_p95_ucl,0.5) , color = "blue" , linetype="dashed",size = 1.5)
        p <- p + geom_hline(yintercept=quantile(simulation_summary$me_p95_ucl,0.5) , color = "green" , linetype="dashed",size = 1.5)
        

        
fivenum( simulation_summary$ideal_p95)                
fivenum( simulation_summary$naive_p95)
fivenum( simulation_summary$me_p95)


fivenum( simulation_summary$ideal_p95_ucl)                
fivenum( simulation_summary$naive_p95_ucl)
fivenum( simulation_summary$me_p95_ucl)


    ###### visuals 3 : EST + IC ####


mydata1 <- data.frame( model = ordered(c("ideal","naive","ME"),
                                       levels = rev(c("ideal","naive","ME"))),
                       lcl = c( quantile(simulation_summary$ideal_p95,0.1),
                                quantile(simulation_summary$naive_p95,0.1),
                                quantile(simulation_summary$me_p95,0.1)),
                       est = c( quantile(simulation_summary$ideal_p95,0.5),
                                quantile(simulation_summary$naive_p95,0.5),
                                quantile(simulation_summary$me_p95,0.5)),
                       ucl = c( quantile(simulation_summary$ideal_p95,0.9),
                                quantile(simulation_summary$naive_p95,0.9),
                                quantile(simulation_summary$me_p95,0.9)))

p <- ggplot(mydata1) 

p <- p +  geom_segment(aes( x = model , xend = model , y = lcl , yend = ucl )) + coord_flip() + ylab("Concentration (ppm)")

p <- p + geom_point( aes( x = model , y = est) , size = 2)

p <- p + geom_hline( yintercept=100, color="red", linetype="dashed", size=1)

p <- p + theme_solarized()

p


###### visuals 3 : 70% UCL + IC ####


mydata1 <- data.frame( model = ordered(c("ideal","naive","ME"),
                                       levels = rev(c("ideal","naive","ME"))),
                       lcl = c( quantile(simulation_summary$ideal_p95_ucl,0.1),
                                quantile(simulation_summary$naive_p95_ucl,0.1),
                                quantile(simulation_summary$me_p95_ucl,0.1)),
                       est = c( quantile(simulation_summary$ideal_p95_ucl,0.5),
                                quantile(simulation_summary$naive_p95_ucl,0.5),
                                quantile(simulation_summary$me_p95_ucl,0.5)),
                       ucl = c( quantile(simulation_summary$ideal_p95_ucl,0.9),
                                quantile(simulation_summary$naive_p95_ucl,0.9),
                                quantile(simulation_summary$me_p95_ucl,0.9)))

p <- ggplot(mydata1) 

p <- p +  geom_segment(aes( x = model , xend = model , y = lcl , yend = ucl )) + coord_flip() + ylab("Concentration (ppm)")

p <- p + geom_point( aes( x = model , y = est) , size = 2)

p <- p + theme_solarized()

p


###### visuals 4 : gm + boxplot ####


mcmc.data <- data.frame( p95_ucl = c( simulation_summary$ideal_p95_ucl,
                                  simulation_summary$naive_p95_ucl,
                                  simulation_summary$me_p95_ucl),
                         p95 = c( simulation_summary$ideal_p95,
                                      simulation_summary$naive_p95,
                                      simulation_summary$me_p95),
                         gm = c( simulation_summary$ideal_gm,
                                  simulation_summary$naive_gm,
                                  simulation_summary$me_gm),
                         gsd = c( simulation_summary$ideal_gsd,
                                 simulation_summary$naive_gsd,
                                 simulation_summary$me_gsd),
                         type = c( rep("p95_ideal", mylength), rep("p95_naive", mylength) , rep("p95_me", mylength)))
f <- function(x) {
  r <- quantile(x, probs = c(0.1, 0.25, 0.5, 0.75, 0.9))
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}



p <- ggplot( mcmc.data, aes(x=type, y=p95_ucl)) 

p <- p +  stat_summary(fun.data = f, geom="boxplot")

p <- p +  scale_y_log10()

p



########################simulated data

data.g <- bayesian.output.D$group.id

data.g$mu <-numeric(length(data.g[,1]))

for (i in 1:length(data.g[,1])) data.g$mu[i]<-median(bayesian.output.D$mu.chain[i,])

data.g$sw <- median(bayesian.output.D$sigma.chain[i,])

dat <-data.frame(value=numeric(0), group=character(0),stringsAsFactors=F)

for (i in 1:length(data.g[,1])) dat <- rbind(dat, data.frame(value=exp(rnorm(1000,data.g$mu[i],data.g$sw[i])),      group=rep(data.g$name[i],1000),stringsAsFactors=F))


######observed data


data.fin <-data.frame(x.orig=data.formatted$data)

data.fin$x <-data.simply.imputed$imputed$data$xfin

data.fin$censored=!data.formatted$notcensored

data.fin$id <-data.formatted$var


#####min max for graph

min.val <-exp(min(data.g$mu-1.5*data.g$sw))
max.val <-exp(max(data.g$mu-1.5*data.g$sw))


####graph





p1 <- ggplot(dat, aes(x=group, y=value, fill = group, color=group, group=group),ylim=c(min.val,max.val))

if(pal_col == TRUE){
  
  p1 <- p1 + geom_boxplot(lwd = 0.8, alpha = 0.1) +
    geom_point(position=position_jitter(width=0.3), alpha=0.2) +
    geom_hline(yintercept = c.oel,colour = "red", size = 2)
  
  
}else{
  
  p1 <- p1 + geom_boxplot(color= "gray20", lwd = 0.8, fill= alpha(paste0("gray", round(seq(5, 25, length.out = length(unique(dat$group))))), 0.3)) +
    geom_point(position=position_jitter(width=0.3), alpha=0.2, color = "black") +
    geom_hline(yintercept = c.oel,colour = "black", size = 2)
}


p1 <- p1 + geom_point(data=data.fin, fill="white", colour = "black", shape=20, size = 4, aes(x=id, y=x, pch=censored, group=id)) +
  scale_y_log10(breaks=c(.01,.1,1,10,100),labels=c(.01,.1,1,10,100)) +
  labs(x = boxplot.cat.1, y = boxplot.cat.2) +
  theme(axis.title.x=element_text(size=16,vjust=-1)) +
  theme(axis.text.x=element_text(size=14)) +
  theme(axis.title.y=element_text(size=16,angle=90)) +
  theme(axis.text.y=element_text(size=14,angle=0)) +
  theme(aspect.ratio=0.6) +
  theme(legend.position = "none") +
  annotate("text", 2.5, c.oel*1.5, label = boxplot.cat.3, size=5 , color="red") +
  theme(axis.line = element_line(size = 3, colour = "grey80")) +
  theme(axis.ticks = element_line(size = 2))   +
  coord_flip(xlim = NULL, ylim = NULL, expand = TRUE)

return(p1)
