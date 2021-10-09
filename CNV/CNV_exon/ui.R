library(ggplot2)
library(plotly)
library(shiny)
require(shinyBS)
require(shinydashboard)
require(shinyjs)
library(tidyr)
library(data.table)
options(encoding = 'UTF-8')

geneData <- fread("/share_data/wujm/project/CNV/gene_exon/exonList/Gene/allGeneList.txt",header=FALSE)
geneList <- as.list(geneData)


shinyUI(fluidPage(
  # Application title
  titlePanel("CNV exon"),
  # includeCSS("style.css"),
  sidebarLayout(
    sidebarPanel(width = 3,
      selectInput('selectGene',label = 'Select Gene :', choices = geneList),
      textInput("sampleName", label="样本名称:" , value = ""),
      actionButton("goButton",label = "Go!")
    ),
    #Show plot 
    mainPanel(width = 9,
        plotlyOutput("exoncnv",width ='auto')

    )
  )
))
