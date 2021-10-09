library(ggplot2)
library(plotly)
library(shiny)
require(shinyBS)
require(shinydashboard)
require(shinyjs)
library(tidyr)
options(encoding = 'UTF-8')


chrList = c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY")

shinyUI(fluidPage(
  # Application title
  titlePanel("CNV_DD"),
  includeCSS("style.css"),
  sidebarLayout(
    sidebarPanel(width = 3,
      selectInput("product_type", label="产品类型:", choices = c("WES")),
      selectInput("sex", label="性别:", choices = c("男","女")),
      textInput("sampleName", label="样本名称:" , value = ""),
      actionButton("goButton",label = "Go!"),
      downloadButton('downloadData', 'Download', style="float:right;")
    ),
    
    #Show plot 
    mainPanel(width = 9,
      navbarPage(
        title = ">  ",
        tabPanel("1",fluidPage(
          plotOutput("chr1_chrom",height="200px"),
          plotlyOutput("chr1")
        )),
        tabPanel("2",fluidPage(
          plotOutput("chr2_chrom",height="200px"),
          plotlyOutput("chr2")
        )),
        tabPanel("3",fluidPage(
          plotOutput("chr3_chrom",height="200px"),
          plotlyOutput("chr3")
        )),
        tabPanel("4",fluidPage(
          plotOutput("chr4_chrom",height="200px"),
          plotlyOutput("chr4")
        )),
        tabPanel("5",fluidPage(
          plotOutput("chr5_chrom",height="200px"),
          plotlyOutput("chr5")
        )),
        tabPanel("6",fluidPage(
          plotOutput("chr6_chrom",height="200px"),
          plotlyOutput("chr6")
        )),
        tabPanel("7",fluidPage(
          plotOutput("chr7_chrom",height="200px"),
          plotlyOutput("chr7")
        )),
        tabPanel("8",fluidPage(
          plotOutput("chr8_chrom",height="200px"),
          plotlyOutput("chr8")
        )),
        tabPanel("9",fluidPage(
          plotOutput("chr9_chrom",height="200px"),
          plotlyOutput("chr9")
        )),
        tabPanel("10",fluidPage(
          plotOutput("chr10_chrom",height="200px"),
          plotlyOutput("chr10")
        )),
        tabPanel("11",fluidPage(
          plotOutput("chr11_chrom",height="200px"),
          plotlyOutput("chr11")
        )),
        tabPanel("12",fluidPage(
          plotOutput("chr12_chrom",height="200px"),
          plotlyOutput("chr12")
        )),
        tabPanel("13",fluidPage(
          plotOutput("chr13_chrom",height="200px"),
          plotlyOutput("chr13")
        )),
        tabPanel("14",fluidPage(
          plotOutput("chr14_chrom",height="200px"),
          plotlyOutput("chr14")
        )),
        tabPanel("15",fluidPage(
          plotOutput("chr15_chrom",height="200px"),
          plotlyOutput("chr15")
        )),
        tabPanel("16",fluidPage(
          plotOutput("chr16_chrom",height="200px"),
          plotlyOutput("chr16")
        )),
        tabPanel("17",fluidPage(
          plotOutput("chr17_chrom",height="200px"),
          plotlyOutput("chr17")
        )),
        tabPanel("18",fluidPage(
          plotOutput("chr18_chrom",height="200px"),
          plotlyOutput("chr18")
        )),
        tabPanel("19",fluidPage(
          plotOutput("chr19_chrom",height="200px"),
          plotlyOutput("chr19")
        )),
        tabPanel("20",fluidPage(
          plotOutput("chr20_chrom",height="200px"),
          plotlyOutput("chr20")
        )),
        tabPanel("21",fluidPage(
          plotOutput("chr21_chrom",height="200px"),
          plotlyOutput("chr21")
        )),
        tabPanel("22",fluidPage(
          plotOutput("chr22_chrom",height="200px"),
          plotlyOutput("chr22")
        )),
        tabPanel("X",fluidPage(
          plotOutput("chrX_chrom",height="200px"),
          plotlyOutput("chrX")
        )),
        tabPanel("Y",fluidPage(
          plotOutput("chrY_chrom",height="200px"),
          plotlyOutput("chrY")
        )),
        tabPanel("TABLE",fluidPage(
          DT::dataTableOutput("table")
        ))
      )
    )
  )
))
