#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#

library(shiny)
library(utils)
library(shinythemes)

# Define UI for application that draws a histogram
shinyUI(fluidPage(
  
  theme = shinytheme("flatly"),
  
  titlePanel("Performance of screening tests in industrial hygiene sampling strategies"),
  
  mainPanel(
                  br(),
                  
                  p("The question answered here through simulation is : given a true exposure distribution, what are the chances of the screening test yielding a PASS"),
                  
                  p("principle of a screening test studied here : PASS if the maximum of a small number of measurements is lower than a fraction of the OEL"),
                  
                  br(),
                  
                  p("Choose the number of samples in the screening procedure"),
                 
                  sliderInput("sample_size", label = ("Number of measurements"),min = 1, max = 10, value = 3),
                  
                  p("Select the fraction of the OEL used for the test"),
                  
                  numericInput("screening_threshold", label = ("Screening threshold"),min = 0.001, max = 0.999, value = 0.1),
                  
                  p("Select true exceedance of the OEL (in %)"),
                  
                  numericInput("true_exceedance", label = ("True exceedance(%)"),min = 0.01, max = 99.9, value = 10),
                  
                  p("Select true GSD"),
                  
                  numericInput("true_gsd", label = ("True GSD"),min = 1.01, max = 20, value = 2.5),
                  
                  br(),
                  
                  h3("Proportion of time the screening test will yield a PASS (in %, across 100 000 trials)"),
                  
                  br(),
                  
                  h4("For the selected GSD : ",strong(textOutput("known_gsd_result" , inline = TRUE)),"%"),
                  
                  br(),
                  
                  h4("GSD unknown* : ",strong(textOutput("unknown_gsd_result" , inline = TRUE)),"%"),
                  
                  br(),
                  
                  p("*: GSDs are sampled for each iteration of the simulation among the dataset reported by Krouomout et al")
            
                
         
                
               
                
                
       )          

    

))
                
                
  

