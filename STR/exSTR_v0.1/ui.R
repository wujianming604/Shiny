library(ggplot2)
library(plotly)
library(shiny)
require(shinyBS)
require(shinydashboard)
require(shinyjs)
library(tidyr)
options(encoding = 'UTF-8')


locus_type = c("DM1","DRPLA","HD","HDL2","SBMA","SCA1","SCA2","SCA3","SCA6","SCA7","SCA12","SCA17","GDPAG")
sampleList <- c()
for(i in readLines(file("/share_data/wujm/project/run_exSTR/run_WES_T345V1/sampleList.txt"))){
  sampleList <- c(sampleList,i)
}



shinyUI(fluidPage(
  # Application title
  titlePanel("repeat STR"),
  sidebarLayout(
    sidebarPanel(width = 3,
      selectInput("product_type", label="类型:", choices = c("WES_T345V1","WES_T192V1","DD","诺禾","WGS"), selected="WES_T345V1"),
      selectInput("sampleName", label="样本名称:" , choices = sampleList),
      #textInput("sampleName", label="样本名称:" , value = ""),
      selectInput("locus", label="Locus", choices = locus_type),
      actionButton("goButton",label = "Go!"),
      downloadButton('downloadData', 'Download', style="float:right;")
    ), 
    
    #Show plot 
    mainPanel(width = 9,
      navbarPage(
        title = ">  ",
        tabPanel("single locus",fluidPage(
          imageOutput("single_locus")
        )),
        tabPanel("seven loci",fluidPage(
          imageOutput("seven_loci")
        )),
        tabPanel("six loci",fluidPage(
          imageOutput("six_loci")
        )),
        tabPanel("significant table",fluidPage(
          DT::dataTableOutput("table")
        ))
      )
    )
  )
))
