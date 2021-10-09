library(ggplot2)
library(RMySQL)
library(yaml)
library(Hmisc)
library(dplyr)
library(data.table)


server <- function(input, output, session){
  for(con in dbListConnections(MySQL())) dbDisconnect(con)
  observeEvent(input$goButton,{
    #处理input sample
    sampleList = strsplit(input$sampleName,",")[[1]]
    inputGeneName <- input$selectGene
    rawData <- as.data.frame(matrix(nrow=0,ncol=9))
    reRawData <- select(rawData,chrom=1,start=2,end=3,symbol=4,transcript=5,exon_number=6,depth=7,corr_depth=8,sample_name=9)

    for(testSample in sampleList){
      conn <- dbConnect(MySQL(), host="192.168.99.7", dbname="nyuen_clinepilepsy", user='wujm', password='wjM123456++')
      cmd <- paste0("select depth_file,exon_file from shiny_exon_cnv where sample_name = '",testSample,"';")
      data <- dbGetQuery(conn, cmd)
      #get each position depth
      #data$dept_file
      depthFile <- data$depth_file
      #get echo exon region depth
      #data$exon_file
      exonDepthFile <- data$exon_file
      #read data
      depthData <- fread(depthFile, data.table = FALSE,sep="\t",header = FALSE)
      exonDepthData <- read.table(exonDepthFile,sep="\t", header = FALSE, col.names = c("chrom","start","end","symbol","transcript","exon_number","depth"))
      exonData <- filter(exonDepthData,symbol == inputGeneName)
      
      # get chrom of input gene
      thisChrom <- exonData[1,"chrom"]
      #染色体的总深度
      chromDepth <- colSums(depthData["V3"])
      #get Correct depth
      exonData["corr_depth"] = exonData['depth'] / (exonData['end'] - exonData['start']) / chromDepth
      exonData["sample_name"] = testSample

      reRawData <- rbind(reRawData, exonData)
    }

    #plotly output
    output$exoncnv <- renderPlotly({
      ggplotly(
        ggplot(reRawData, aes(x=exon_number,y=corr_depth,fill=sample_name))+
        geom_bar(position = "dodge",stat="identity",width=0.5) + 
        scale_x_continuous(breaks = reRawData$exon_number) + 
        theme(plot.title = element_text(hjust = 0.5),axis.text.y = element_blank(),axis.ticks.y=element_blank(),legend.title=element_blank())+
        labs(title=inputGeneName, y="")
      )
    })

  })
}