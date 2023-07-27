###################################################################################################
#
#
#      Strat sim project : Using real values of GSD picked at random
#
#      preparation of the gsd data from the Kroumhout et al's paper 
#
#
##########################################################################################################


# libraries --------------------------

library(readxl)

paper_data <- read_xlsx("raw data/kromhout databasefromPaperforR.xlsx")

# gsd calculation

paper_data$gsd <- exp( sqrt( paper_data$sw^2 + paper_data$sb^2 ) )


saveRDS( paper_data$gsd , "created data/real_gsd_values.RDS")
