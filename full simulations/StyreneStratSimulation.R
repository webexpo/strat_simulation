#################################################################################
#
#
#   simulating the outcome of risk analysis (CEN) using Rappaport Styrene data
#
##################################################################################


# 1. libraries ---------------------------------------------

    #library(usethis)
    #library(here)
    #library(httr)
    #library(magrittr)
    
    
    library(ggplot2)
    library(ggthemes)
    library(writexl)
    library(scales)
    devtools::source_url("https://github.com/lhimp/scripts/raw/master/chemin.R")
    
    ##  Webexpo functions
    source("parameter estimation/bayesian/load.webexpo.SEG.functions.R")

    ##  arithmetic mean (MVUE + Land UCL95) function
    source("parameter estimation/frequentist/arithmetic_mean.R")
    

        
# 2. data -----------------------------------------------------------------------------------


        mydata <- read.csv( "raw data/rappaport/styr.rappap.csv" , header = TRUE , sep = "," , dec = "." )
          

# 3. functions -----------------------------------------------

        # function calculating P95 estimate + 70%UC + 95%UCL
        
       fun.expostats <- function( x , oel) {
         
         # x is a numerical vector
         
         x <- as.character(x)
         
         mcmc.seg.1 <- Webexpo.seg.globalbayesian.jags( data.sample = x ,
                                                        is.lognormal = TRUE , 
                                                        error.type = "none" ,
                                                        oel = oel )
         
         #Numerical interpretation of the MCMC chains _ 70%
         
         num.seg.1 <- all.numeric(mu.chain = mcmc.seg.1$mu.chain,
                                  sigma.chain = mcmc.seg.1$sigma.chain,
                                  probacred = 40 ,
                                  oel = oel,
                                  frac_threshold =5 ,
                                  target_perc = 95)
         
         #Numerical interpretation of the MCMC chains _ 95%
         
         num.seg.2 <- all.numeric(mu.chain = mcmc.seg.1$mu.chain,
                                  sigma.chain = mcmc.seg.1$sigma.chain,
                                  probacred = 90 ,
                                  oel = oel,
                                  frac_threshold =5 ,
                                  target_perc = 95)
         
         result <- c(num.seg.1$perc$est,
                     ucl70 = num.seg.1$perc$ucl,
                     ucl95 = num.seg.2$perc$ucl)
         
         names(result) <- c("est","UCL70","UCL95")
         
         return(result)
         
         
       }
      
       
       # function to apply 5 strategies
       
       fun.five.strats <- function( x , oel) {
         
            result <- vector( mode = "numeric" , length = 5)
       
            bayes.res <- fun.expostats( x = x , oel = oel )
            
            result[1] <- max( x )
            
            result[2] <- x[1]
            
            result[3] <- bayes.res[1]
            
            result[4] <- bayes.res[2]
            
            result[5] <- bayes.res[3]
            
            names(result) <- c( "max" , "AL" , "P95.est" , "P95.UCL70" , "P95.UCL95")
            
            return(result) 
            
       }
            
            
            
         

       # function to simulate a sample and apply  5 strategies
       
       fun.simulate <- function( n.sim , oel , n.sample , truepop  ) {
         
         #n.sim <- 100
         
         #oel <-350
         
         #n.sample <- 6
         
         #truepop <- mydata$x
         
         #matrix of samples
         
         samp.mat <- matrix( nrow = n.sample , ncol = n.sim)
         
         for (i in 1:n.sim ) samp.mat[ , i ] <- sample( x = truepop , size = n.sample , replace = TRUE )
         
         res.mat <- apply( samp.mat , 2 , fun.five.strats , oel = oel) 
         
         return(res.mat)
         
       }
         
         
      
       
# 4. MAIN analysis-----------------------------------------
        
        
      ##### 4.1 OEL
        
        oel <- 350
        
      #### 4.2 whole data set analysis
        
       main <- fun.expostats( x =  mydata$x , oel = oel) 
        
       main.fivestrats <- fun.five.strats( x =  mydata$x  , oel = oel ) 
        
      ##### 4.3 simulation : n=5000 6 samples OEL = 350
       
       results.1 <- fun.simulate( n.sim = 5000  , oel = 350 , n.sample = 6 , truepop = mydata$x  ) 
       
       
       
       ######## 4.4 interpretation for OEL = 350
       
       oel <- 350
       
       strat1 <- logical( length( results.1[ 1 , ] ) )
       
       for (i in 1:length( results.1[ 1 , ] ) ) strat1[i] <- ( results.1[ 1 , i ] < oel )
       
       strat2 <- logical( length( results.1[ 1 , ] ) )
       
       for (i in 1:length( results.1[ 1 , ] ) ) strat2[i] <- ( results.1[ 2 , i ] < oel/2 )
       
       strat3 <- logical( length( results.1[ 1 , ] ) )
       
       for (i in 1:length( results.1[ 1 , ] ) ) strat3[i] <- ( results.1[ 3 , i ] < oel )
       
       strat4 <- logical( length( results.1[ 1 , ] ) )
       
       for (i in 1:length( results.1[ 1 , ] ) ) strat4[i] <- ( results.1[ 4 , i ] < oel )
       
       strat5 <- logical( length( results.1[ 1 , ] ) )
       
       for (i in 1:length( results.1[ 1 , ] ) ) strat5[i] <- ( results.1[ 5 , i ] < oel )
       
       100*length(strat1[!strat1])/length(strat1) #16
       100*length(strat1[!strat2])/length(strat2) #34
       100*length(strat1[!strat3])/length(strat3) #63
       100*length(strat1[!strat4])/length(strat4) #86
       100*length(strat1[!strat5])/length(strat5) #99
       
       
       
       ######## 4.5 interpretation for OEL = 700
       
       oel <- 700
       
       strat1 <- logical( length( results.1[ 1 , ] ) )
       
       for (i in 1:length( results.1[ 1 , ] ) ) strat1[i] <- ( results.1[ 1 , i ] < oel )
       
       strat2 <- logical( length( results.1[ 1 , ] ) )
       
       for (i in 1:length( results.1[ 1 , ] ) ) strat2[i] <- ( results.1[ 2 , i ] < oel/2 )
       
       strat3 <- logical( length( results.1[ 1 , ] ) )
       
       for (i in 1:length( results.1[ 1 , ] ) ) strat3[i] <- ( results.1[ 3 , i ] < oel )
       
       strat4 <- logical( length( results.1[ 1 , ] ) )
       
       for (i in 1:length( results.1[ 1 , ] ) ) strat4[i] <- ( results.1[ 4 , i ] < oel )
       
       strat5 <- logical( length( results.1[ 1 , ] ) )
       
       for (i in 1:length( results.1[ 1 , ] ) ) strat5[i] <- ( results.1[ 5 , i ] < oel )
       
       100*length(strat1[!strat1])/length(strat1) #0
       100*length(strat1[!strat2])/length(strat2) #3
       100*length(strat1[!strat3])/length(strat3) #7
       100*length(strat1[!strat4])/length(strat4) #25
       100*length(strat1[!strat5])/length(strat5) #87
       
       
##################################   5. AIHA risk band
       
       
       ##################  5.1 FR
         
         c.oel <- 700
         
         perc.chain <- results.1[ 3 ,]
         
         C1 <-100*length(perc.chain[perc.chain<0.01*c.oel])/length(perc.chain)
         C2 <-100*length(perc.chain[perc.chain>=0.01*c.oel & perc.chain<0.1*c.oel])/length(perc.chain)
         C3 <-100*length(perc.chain[perc.chain>=0.1*c.oel & perc.chain<0.5*c.oel])/length(perc.chain)
         C4 <-100*length(perc.chain[perc.chain>=0.5*c.oel & perc.chain<1*c.oel])/length(perc.chain)
         C5 <-100*length(perc.chain[perc.chain>=1*c.oel ])/length(perc.chain)

         riskplot.2="Proportion des conclusions"
         riskplot.3="<1%\nVLE"
         riskplot.4="1-10%\nVLE"
         riskplot.5="10-50%\nVLE"
         riskplot.6="50-100%\nVLE"
         riskplot.7=">VLE"         
         
         cats <- factor(c('C1','C2','C3','C4','C5'),labels=c(riskplot.3,riskplot.4,riskplot.5,riskplot.6,riskplot.7))
         
         data <-data.frame(perc=c(C1,C2,C3,C4,C5),cat=cats)

         graph9<- ggplot(data,aes(x=cats,y=perc))
         graph9 <-graph9+
           geom_bar(stat="identity",fill=c('green3','cyan3','yellow','orange','red'))+
           theme(aspect.ratio=0.6)+
           xlab(expression ( "Catégorie du 95"^"e" * "centile"))+
           ylab(riskplot.2) +
           theme(axis.title.x=element_text(size=16,vjust=-1))+
           theme(axis.text.x=element_text(size=13))+
           theme(axis.title.y=element_text(size=16,angle=90))+
           theme(axis.text.y=element_text(size=13,angle=90,hjust=0.5))+
           theme(legend.position = "none")+
           geom_text(x=1:5,y=c(C1,C2,C3,C4,C5)+5,label=paste(signif(c(C1,C2,C3,C4,C5),3),'%',sep=''),size=5,colour='grey28') +
           scale_y_continuous(breaks=c(0,20,40,60,80,100),limits=c(0,110))
         
        
         
         ##################  5.2 EN
         
         c.oel <- 700
         
         perc.chain <- results.1[ 3 ,]
         
         C1 <-100*length(perc.chain[perc.chain<0.01*c.oel])/length(perc.chain)
         C2 <-100*length(perc.chain[perc.chain>=0.01*c.oel & perc.chain<0.1*c.oel])/length(perc.chain)
         C3 <-100*length(perc.chain[perc.chain>=0.1*c.oel & perc.chain<0.5*c.oel])/length(perc.chain)
         C4 <-100*length(perc.chain[perc.chain>=0.5*c.oel & perc.chain<1*c.oel])/length(perc.chain)
         C5 <-100*length(perc.chain[perc.chain>=1*c.oel ])/length(perc.chain)
         
         riskplot.2="Proportion of decisions"
         riskplot.3="<1%\nOEL"
         riskplot.4="1-10%\nOEL"
         riskplot.5="10-50%\nOEL"
         riskplot.6="50-100%\nOEL"
         riskplot.7=">OEL"         
         
         cats <- factor(c('C1','C2','C3','C4','C5'),labels=c(riskplot.3,riskplot.4,riskplot.5,riskplot.6,riskplot.7))
         
         data <-data.frame(perc=c(C1,C2,C3,C4,C5),cat=cats)
         
         graph9<- ggplot(data,aes(x=cats,y=perc))
         graph9 <-graph9+
           geom_bar(stat="identity",fill=c('green3','cyan3','yellow','orange','red'))+
           theme(aspect.ratio=0.6)+
           xlab(expression ( "95"^"th" * "percentile category"))+
           ylab(riskplot.2) +
           theme(axis.title.x=element_text(size=16,vjust=-1))+
           theme(axis.text.x=element_text(size=13))+
           theme(axis.title.y=element_text(size=16,angle=90))+
           theme(axis.text.y=element_text(size=13,angle=90,hjust=0.5))+
           theme(legend.position = "none")+
           geom_text(x=1:5,y=c(C1,C2,C3,C4,C5)+5,label=paste(signif(c(C1,C2,C3,C4,C5),3),'%',sep=''),size=5,colour='grey28') +
           scale_y_continuous(breaks=c(0,20,40,60,80,100),limits=c(0,110))
         
         
         
# 6. graphiques -----------------------------------------
       
       
       data.to.plot <- data.frame( max = res.mat[ 1 , ],
                                   al = res.mat[ 2 , ],
                                   P95est = res.mat[ 3 , ],
                                   P95ucl70 = res.mat[ 4 , ],
                                   P95ucl95 = res.mat[ 5 , ] )
       
       #P95 est
       
       p <-ggplot(data=data.to.plot, aes(x=P95ucl70)) + geom_histogram(aes(y =..density.. ) , color="white", fill="darkgray") #+ xlim(c(0,120)) 
       
       p <- p + theme_calc()
       
       p <-p + labs(x=expression ( "Concentration (mg/m"^"3" * ")") , y="Probability density")
       
       p <- p + geom_segment(aes(x = oel , y = 0, 
                                 xend = oel , 
                                 yend = 0.004 ),size=1.3,color="black")
       
       p <-p +theme(axis.title.x=element_text(vjust=-0.5,size=16))+
         theme(axis.title.y=element_text(size=12,angle=90))+
         theme(axis.text.x=element_text(size=14))+
         theme(axis.text.y=element_text(size=10,angle=90, hjust=0.3))
       
       #p <- p + scale_y_continuous(expand = c(0, 0) , 
        #                           limits = c(0,0.0025),
        #                           labels = scales::number_format(accuracy = 0.001,
        #                                                          decimal.mark = ',') ) 
       
       p <- p + scale_x_continuous( expand = c(0, 0) , limits = c(0,1500), breaks = seq( from = 0 , to = 1500, by = 500) )  
       
       
       p <- p + theme( plot.background = element_rect(color = NA) )
       
fivenum(data.to.plot$P95est)       
fivenum(data.to.plot$P95ucl70)       
fivenum(data.to.plot$P95ucl95)       

       
       
       
        

# 7. LEAP-TECK strategy -----------------------------------------
#
#   Decision rule : compliant ("green") if the 95% UCL on the arithmetic mean
#   (Land's exact method, fun.am) is below the OEL ; otherwise "not green"
#

       ##### 7.0 settings

       oel.lt <- 350          # OEL used for the LEAP-TECK assessment (TO CONFIRM)

       n.sample.lt <- 6       # number of measurements per participant in the simulation (TO CONFIRM)

       n.sim.lt <- 5000       # number of simulated participants

       # helper : LEAP-TECK decision for one vector of measurements

       fun.leapteck <- function( x , oel ) {

         am <- fun.am( x = x , L = oel )

         c( mvue = am$mvue , ucl95 = am$ucl , green = as.numeric( am$ucl < oel ) )

       }


       ##### 7.1 participants data

       # one row per participant, 6 measurements each (mg/m3)

       participants <- read.csv( text = "
id,x1,x2,x3,x4,x5,x6
1,179,154,256,72.6,221,70.5
2,188,6.64,361,384,131,222
3,241,66.7,96.8,72.6,82.3,19.4
4,88,74.2,74.2,304,208,221
5,2.78,285,180,185,18.7,19.4
6,6.64,58.5,35.1,88,176,287
7,190,158,208,241,184,56
8,285,82.7,71,222,39.5,66.7
9,108,82,210,158,42.4,195
10,154,189,118,179,74.4,18.7
12,211,214,72.6,210,67.9,66.9
13,287,199,10.2,139,6.32,48
14,86.2,2.78,49,157,361,31.6
15,82.7,66.7,153,69.2,139,118
16,43.6,195,47,56,264,30.6
17,6.64,118,56,195,234,176
18,276,66.6,48,131,131,69.6
22,42.4,50.7,180,48,31.6,109
23,158,210,154,133,109,109
24,181,188,391,82.3,189,35.1
25,195,118,115,67.9,47,256
26,153,72.6,139,74.2,109,245
27,2.78,10.2,176,35.1,32.1,181
28,71,155,140,32.1,204,188
29,21.6,361,118,66.6,6.32,139
30,140,43.6,82.3,19.4,89.5,234
31,154,56,18.7,157,49,304
32,176,84.7,81.1,199,69.2,96.8
33,10.2,21.6,31.6,32.1,140,10.2
34,74.4,304,31.6,208,304,287
35,214,6.64,88,72.6,35.1,285
36,96.8,70.5,66.9,70.5,181,176
" )

       res.participants <- data.frame( id = participants$id ,
                                       t( apply( participants[ , paste0( "x" , 1:6 ) ] , 1 ,
                                                 fun.leapteck , oel = oel.lt ) ) )

       res.participants

       # number and % of participants concluding green

       sum( res.participants$green )

       100 * mean( res.participants$green )


       ##### 7.2 simulation : n.sim.lt participants drawing n.sample.lt values from the styrene dataset

       set.seed( 12345 )

       samp.mat.lt <- matrix( sample( x = mydata$x , size = n.sample.lt * n.sim.lt , replace = TRUE ) ,
                              nrow = n.sample.lt , ncol = n.sim.lt )

       res.sim.lt <- t( apply( samp.mat.lt , 2 , fun.leapteck , oel = oel.lt ) )

       # number and % of simulated participants concluding green

       sum( res.sim.lt[ , "green" ] )

       100 * mean( res.sim.lt[ , "green" ] )

       # reference : the full styrene dataset assessed with the same rule

       fun.leapteck( x = mydata$x , oel = oel.lt )

       # distribution of the UCL95 across simulated participants

       fivenum( res.sim.lt[ , "ucl95" ] )
