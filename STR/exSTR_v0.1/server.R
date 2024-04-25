library(ggplot2)
library(plotly)
library(data.table)
library(exSTRa)

wes_locus_type = c("DM1","DRPLA","HD","HDL2","SBMA","SCA1","SCA2","SCA3","SCA6","SCA7","SCA8","SCA12","SCA17","GDPAG")
wgs_locus_type = c('DM1','DM2','DRPLA','EPM1A','FRAXA','FRAXE','FTDALS1','HD','HDL2','SBMA','SCA1','SCA2','SCA3','SCA7','SCA8','SCA10','SCA12','SCA17','SCA36','FECD3','GDPAG','NIID','OPDM1','OPML1','DBQD2','SCA6','FAME1')
nuohe_locus_type = c("DRPLA","HD","HDL2","SBMA","SCA1","SCA2","SCA3","SCA6","SCA7","SCA12","SCA17")

sample_T345V1_List <- c()
for(i in readLines(file("/share_data/wujm/project/run_exSTR/run_WES_T345V1/sampleList.txt"))){
  sample_T345V1_List <- c(sample_T345V1_List,i)
}
sample_T192V1_List <- c()
for(i in readLines(file("/share_data/wujm/project/run_exSTR/run_WES_T192V1/sampleList.txt"))){
  sample_T192V1_List <- c(sample_T192V1_List,i)
}
sample_T192V1Plus_List <- c()
for(i in readLines(file("/share_data/wujm/project/run_exSTR/run_WES_T192V1Plus/sampleList.txt"))){
  sample_T192V1Plus_List <- c(sample_T192V1Plus_List,i)
}
sample_NUOHE_List <- c()
for(i in readLines(file("/share_data/wujm/project/run_exSTR/run_NUOHE/sampleList.txt"))){
  sample_NUOHE_List <- c(sample_NUOHE_List,i)
}
sample_WGS_List <- c()
for(i in readLines(file("/share_data/wujm/project/run_exSTR/run_WGS/sampleList.txt"))){
  sample_WGS_List <- c(sample_WGS_List,i)
}


server <- function(input, output, session){
    #联动
    observe({
      x <- input$product_type
      if(x == "WES_T345V1"){
        y <- sample_T345V1_List
        locus_type <- wes_locus_type
      } else if(x == "DD"){
        y <- sample_T192V1Plus_List
        locus_type <- wes_locus_type
      } else if(x == "WES_T192V1"){
        y <- sample_T192V1_List
        locus_type <- wes_locus_type
      } else if(x == "诺禾"){
        y <- sample_NUOHE_List
        locus_type <- nuohe_locus_type
      } else if(x == "WGS"){
        y <- sample_WGS_List
        locus_type <- wgs_locus_type
      }
      if(x != "EE"){
        updateSelectInput(session, "sampleName",
                          label = "样本名称",
                          choices = y)
        updateSelectInput(session, "locus",
                          label = "Locus",
                          choices = locus_type)
      }# 根据不同的样本类型，更新不同的样本名称

    })
    observeEvent(input$goButton,{
        if(input$product_type == "WES_T345V1"){
            sampleName <- input$sampleName
            locusName <- input$locus
            workDir <- "/share_data/wujm/project/run_exSTR/run_WES_T345V1/"
            # set var
            colo <- c(sample = "red")
            names(colo) <- c(sampleName)
            ##
            str_score <- read_score(
              file <- paste0(workDir,"samples_T345V1.exSTR2.txt"), 
              database = system.file("extdata", "repeat_expansion_disorders_hg19.txt", package = "exSTRa"),
              groups.regex = c(control = "^WGSrpt_0[24]$", case = "")
            )

            str_score_7 <- str_score[c("DM1","DRPLA","HD","HDL2","SBMA","SCA1","SCA2")]
            str_score_6 <- str_score[ c("SCA3","SCA6","SCA7","SCA8","SCA12","SCA17","GDPAG")]


            output$single_locus <- renderImage({
              # A temp file to save the output.
              # This file will be removed later by renderImage
              outfile <- tempfile(fileext = '.png')
              png(outfile, width = 1232,height = 665)
              plot(str_score, locusName, sample_col = colo)
              dev.off()
              # Return a list containing the filename
              list(src = outfile,
                  contentType = 'image/png',
                  width = 1232,
                  height = 665)
            }, deleteFile = TRUE)

            output$seven_loci <- renderImage({
              # A temp file to save the output.
              # This file will be removed later by renderImage
              outfile <- tempfile(fileext = '.png')
              png(outfile, width = 1232,height = 665)
              par(mfrow = c(3, 3))
              for(i in c("DM1","DRPLA","HD","HDL2","SBMA","SCA1","SCA2")){
                plot(str_score_7,i,sample_col = colo)}
              dev.off()
              # Return a list containing the filename
              list(src = outfile,
                  contentType = 'image/png',
                  width = 1232,
                  height = 665)
            }, deleteFile = TRUE)

            output$six_loci <- renderImage({
              # A temp file to save the output.
              # This file will be removed later by renderImage
              outfile <- tempfile(fileext = '.png')
              png(outfile, width = 1232,height = 665)
              par(mfrow = c(3, 3))
              for(i in c("SCA3","SCA6","SCA7","SCA8","SCA12","SCA17","GDPAG")){
                plot(str_score_6,i,sample_col = colo)}
              dev.off()
              # Return a list containing the filename
              list(src = outfile,
                  contentType = 'image/png',
                  width = 1232,
                  height = 665)
            }, deleteFile = TRUE)



            ##output table
            #用于输出significant table，会少一个SCA2
            str_score_12 <- str_score[c("DM1","DRPLA","HD","HDL2","SBMA","SCA1","SCA3","SCA6","SCA7","SCA8","SCA12","SCA17","GDPAG")]
            tsum <- tsum_test(str_score_12)
            ps <- p_values(tsum, only.signif = TRUE, correction = "samples")
            output$table <- DT::renderDataTable({ps})

            outFileName = paste0(sampleName,".T345V1.significant.csv")
            output$downloadData = downloadHandler(
              filename = function() {
                outFileName
                },
              content = function(file) {
                write.table(ps, file,sep=",", row.names= FALSE)
                }
            )
        }
        else if(input$product_type == "WES_T192V1"){
            sampleName <- input$sampleName
            locusName <- input$locus
            workDir <- "/share_data/wujm/project/run_exSTR/run_WES_T192V1/"
            # set var
            colo <- c(sample = "red")
            names(colo) <- c(sampleName)
            ##
            str_score <- read_score(
              file <- paste0(workDir,"samples_T192V1.exSTR2.txt"), 
              database = system.file("extdata", "repeat_expansion_disorders_hg19.txt", package = "exSTRa"),
              groups.regex = c(control = "^WGSrpt_0[24]$", case = "")
            )

            str_score_7 <- str_score[c("DM1","DRPLA","HD","HDL2","SBMA","SCA1","SCA2")]
            str_score_6 <- str_score[ c("SCA3","SCA6","SCA7","SCA8","SCA12","SCA17","GDPAG")]


            output$single_locus <- renderImage({
              # A temp file to save the output.
              # This file will be removed later by renderImage
              outfile <- tempfile(fileext = '.png')
              png(outfile, width = 1232,height = 665)
              plot(str_score, locusName, sample_col = colo)
              dev.off()
              # Return a list containing the filename
              list(src = outfile,
                  contentType = 'image/png',
                  width = 1232,
                  height = 665)
            }, deleteFile = TRUE)

            output$seven_loci <- renderImage({
              # A temp file to save the output.
              # This file will be removed later by renderImage
              outfile <- tempfile(fileext = '.png')
              png(outfile, width = 1232,height = 665)
              par(mfrow = c(3, 3))
              for(i in c("DM1","DRPLA","HD","HDL2","SBMA","SCA1","SCA2")){
                plot(str_score_7,i,sample_col = colo)}
              dev.off()
              # Return a list containing the filename
              list(src = outfile,
                  contentType = 'image/png',
                  width = 1232,
                  height = 665)
            }, deleteFile = TRUE)

            output$six_loci <- renderImage({
              # A temp file to save the output.
              # This file will be removed later by renderImage
              outfile <- tempfile(fileext = '.png')
              png(outfile, width = 1232,height = 665)
              par(mfrow = c(3, 3))
              for(i in c("SCA3","SCA6","SCA7","SCA8","SCA12","SCA17","GDPAG")){
                plot(str_score_6,i,sample_col = colo)}
              dev.off()
              # Return a list containing the filename
              list(src = outfile,
                  contentType = 'image/png',
                  width = 1232,
                  height = 665)
            }, deleteFile = TRUE)



            ##output table
            #用于输出significant table，会少一个SCA2
            str_score_12 <- str_score[c("DM1","DRPLA","HD","HDL2","SBMA","SCA1","SCA3","SCA6","SCA7","SCA8","SCA12","SCA17","GDPAG")]
            tsum <- tsum_test(str_score_12)
            ps <- p_values(tsum, only.signif = TRUE, correction = "samples")
            output$table <- DT::renderDataTable({ps})

            outFileName = paste0(sampleName,".T192V1.significant.csv")
            output$downloadData = downloadHandler(
              filename = function() {
                outFileName
                },
              content = function(file) {
                write.table(ps, file,sep=",", row.names= FALSE)
                }
            )
        }
        else if(input$product_type == "诺禾"){
            sampleName <- input$sampleName
            locusName <- input$locus
            workDir <- "/share_data/wujm/project/run_exSTR/run_NUOHE/"
            # set var
            colo <- c(sample = "red")
            names(colo) <- c(sampleName)
            str_score <- read_score(
              file <- paste0(workDir,"samples_NUOHE.exSTR2.txt"), 
              database = system.file("extdata", "repeat_expansion_disorders_hg19.txt", package = "exSTRa"),
              #database <- "/share_data/wujm/software/exSTRa/inst/extdata/nuohe_repeat_expansion_disorders_hg19.txt",
              groups.regex = c(control = "^WGSrpt_0[24]$", case = "")
            )

            str_score_6 <- str_score[c("DRPLA","HD","HDL2","SBMA","SCA1","SCA2")]
            str_score_5 <- str_score[ c("SCA3","SCA6","SCA7","SCA12","SCA17")]


            output$single_locus <- renderImage({
              # A temp file to save the output.
              # This file will be removed later by renderImage
              outfile <- tempfile(fileext = '.png')
              png(outfile, width = 1232,height = 665)
              plot(str_score, locusName, sample_col = colo)
              dev.off()
              # Return a list containing the filename
              list(src = outfile,
                  contentType = 'image/png',
                  width = 1232,
                  height = 665)
            }, deleteFile = TRUE)

            output$seven_loci <- renderImage({
              # A temp file to save the output.
              # This file will be removed later by renderImage
              outfile <- tempfile(fileext = '.png')
              png(outfile, width = 1232,height = 665)
              par(mfrow = c(3, 3))
              for(i in c("DRPLA","HD","HDL2","SBMA","SCA1","SCA2")){
                plot(str_score_6,i,sample_col = colo)}
              dev.off()
              # Return a list containing the filename
              list(src = outfile,
                  contentType = 'image/png',
                  width = 1232,
                  height = 665)
            }, deleteFile = TRUE)

            output$six_loci <- renderImage({
              # A temp file to save the output.
              # This file will be removed later by renderImage
              outfile <- tempfile(fileext = '.png')
              png(outfile, width = 1232,height = 665)
              par(mfrow = c(3, 3))
              for(i in c("SCA3","SCA6","SCA7","SCA12","SCA17")){
                plot(str_score_5,i,sample_col = colo)}
              dev.off()
              # Return a list containing the filename
              list(src = outfile,
                  contentType = 'image/png',
                  width = 1232,
                  height = 665)
            }, deleteFile = TRUE)



            ##output table
            #用于输出significant table
            str_score_11 <- str_score[c("DRPLA","HD","HDL2","SBMA","SCA1","SCA3","SCA6","SCA7","SCA12")]
            tsum <- tsum_test(str_score_11)
            ps <- p_values(tsum, only.signif = TRUE, correction = "samples")
            output$table <- DT::renderDataTable({ps})

            outFileName = paste0(sampleName,".NUOHE.significant.csv")
            output$downloadData = downloadHandler(
              filename = function() {
                outFileName
                },
              content = function(file) {
                write.table(ps, file,sep=",", row.names= FALSE)
                }
            )
        }
        else if(input$product_type == "DD"){
            sampleName <- input$sampleName
            locusName <- input$locus
            workDir <- "/share_data/wujm/project/run_exSTR/run_WES_T192V1Plus/"
            # set var
            colo <- c(sample = "red")
            names(colo) <- c(sampleName)
            ##
            str_score <- read_score(
              file <- paste0(workDir,"samples_T192V1Plus.exSTR2.txt"), 
              database = system.file("extdata", "repeat_expansion_disorders_hg19.txt", package = "exSTRa"),
              groups.regex = c(control = "^WGSrpt_0[24]$", case = "")
            )

            str_score_7 <- str_score[c("DM1","DRPLA","HD","HDL2","SBMA","SCA1","SCA2")]
            str_score_6 <- str_score[ c("SCA3","SCA6","SCA7","SCA8","SCA12","SCA17","GDPAG")]


            output$single_locus <- renderImage({
              # A temp file to save the output.
              # This file will be removed later by renderImage
              outfile <- tempfile(fileext = '.png')
              png(outfile, width = 1232,height = 665)
              plot(str_score, locusName, sample_col = colo)
              dev.off()
              # Return a list containing the filename
              list(src = outfile,
                  contentType = 'image/png',
                  width = 1232,
                  height = 665)
            }, deleteFile = TRUE)

            output$seven_loci <- renderImage({
              # A temp file to save the output.
              # This file will be removed later by renderImage
              outfile <- tempfile(fileext = '.png')
              png(outfile, width = 1232,height = 665)
              par(mfrow = c(3, 3))
              for(i in c("DM1","DRPLA","HD","HDL2","SBMA","SCA1","SCA2")){
                plot(str_score_7,i,sample_col = colo)}
              dev.off()
              # Return a list containing the filename
              list(src = outfile,
                  contentType = 'image/png',
                  width = 1232,
                  height = 665)
            }, deleteFile = TRUE)

            output$six_loci <- renderImage({
              # A temp file to save the output.
              # This file will be removed later by renderImage
              outfile <- tempfile(fileext = '.png')
              png(outfile, width = 1232,height = 665)
              par(mfrow = c(3, 3))
              for(i in c("SCA3","SCA6","SCA7","SCA8","SCA12","SCA17","GDPAG")){
                plot(str_score_6,i,sample_col = colo)}
              dev.off()
              # Return a list containing the filename
              list(src = outfile,
                  contentType = 'image/png',
                  width = 1232,
                  height = 665)
            }, deleteFile = TRUE)



            ##output table
            #用于输出significant table，会少一个SCA2
            str_score_12 <- str_score[c("DM1","DRPLA","HD","HDL2","SBMA","SCA1","SCA3","SCA6","SCA7","SCA8","SCA12","SCA17","GDPAG")]
            tsum <- tsum_test(str_score_12)
            ps <- p_values(tsum, only.signif = TRUE, correction = "samples")
            output$table <- DT::renderDataTable({ps})

            outFileName = paste0(sampleName,".T192V1Plus.significant.csv")
            output$downloadData = downloadHandler(
              filename = function() {
                outFileName
                },
              content = function(file) {
                write.table(ps, file,sep=",", row.names= FALSE)
                }
            )
        }
        else if(input$product_type == "WGS"){
            sampleName <- input$sampleName
            locusName <- input$locus
            workDir <- "/share_data/wujm/project/run_exSTR/run_WGS/"
            # set var
            colo <- c(sample = "red")
            names(colo) <- c(sampleName)
            ##
            str_score <- read_score(
              file <- paste0(workDir,"samples_WGS.exSTR2.txt"), 
              database = system.file("extdata", "repeat_expansion_disorders_WGS_hg19.txt", package = "exSTRa"),
              groups.regex = c(control = "^WGSrpt_0[24]$", case = "")
            )

            str_score_14 <- str_score[c('DM1','DM2','DRPLA','EPM1A','FRAXA','FRAXE','FTDALS1','HD','HDL2','SBMA','SCA1','SCA2','SCA3','SCA7')]
            str_score_13 <- str_score[c('SCA8','SCA6','SCA10','SCA12','SCA17','SCA36','FECD3','GDPAG','NIID','OPDM1','OPML1','DBQD2','FAME1')]


            output$single_locus <- renderImage({
              # A temp file to save the output.
              # This file will be removed later by renderImage
              outfile <- tempfile(fileext = '.png')
              png(outfile, width = 1232,height = 665)
              plot(str_score, locusName, sample_col = colo)
              dev.off()
              # Return a list containing the filename
              list(src = outfile,
                  contentType = 'image/png',
                  width = 1232,
                  height = 665)
            }, deleteFile = TRUE)

            output$seven_loci <- renderImage({
              # A temp file to save the output.
              # This file will be removed later by renderImage
              outfile <- tempfile(fileext = '.png')
              png(outfile, width = 1232,height = 665)
              par(mfrow = c(4, 4))
              for(i in c('DM1','DM2','DRPLA','EPM1A','FRAXA','FRAXE','FTDALS1','HD','HDL2','SBMA','SCA1','SCA2','SCA3','SCA7')){
                plot(str_score_14,i,sample_col = colo)}
              dev.off()
              # Return a list containing the filename
              list(src = outfile,
                  contentType = 'image/png',
                  width = 1232,
                  height = 665)
            }, deleteFile = TRUE)

            output$six_loci <- renderImage({
              # A temp file to save the output.
              # This file will be removed later by renderImage
              outfile <- tempfile(fileext = '.png')
              png(outfile, width = 1232,height = 665)
              par(mfrow = c(4, 4))
              for(i in c('SCA8','SCA6','SCA10','SCA12','SCA17','SCA36','FECD3','GDPAG','NIID','OPDM1','OPML1','DBQD2','FAME1')){
                plot(str_score_13,i,sample_col = colo)}
              dev.off()
              # Return a list containing the filename
              list(src = outfile,
                  contentType = 'image/png',
                  width = 1232,
                  height = 665)
            }, deleteFile = TRUE)



            ##output table
            #用于输出significant table，会少一个SCA2
            str_score_29 <- str_score[c('DM1','DM2','DRPLA','EPM1A','FRAXA','FRAXE','FTDALS1','HD','HDL2','SBMA','SCA1','SCA2','SCA3','SCA7','SCA8','SCA10','SCA29','SCA17','SCA36','FECD3','GDPAG','NIID','OPDM1','OPML1','DBQD2')]
            tsum <- tsum_test(str_score_29)
            ps <- p_values(tsum, only.signif = TRUE, correction = "samples")
            output$table <- DT::renderDataTable({ps})

            outFileName = paste0(sampleName,".WGS.significant.csv")
            output$downloadData = downloadHandler(
              filename = function() {
                outFileName
                },
              content = function(file) {
                write.table(ps, file,sep=",", row.names= FALSE)
                }
            )
        }
    })
}
