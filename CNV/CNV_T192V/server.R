library(ggplot2)
library(RMySQL)
library(yaml)
library(Hmisc)
library(data.table)


cytoFile <- "/share_data/wujm/project/chromPlot/cytoBandIdeo_hg19.txt"
chromNumber <- function(chrom){
  ideoFile<-read.table(cytoFile,header = FALSE,col.names=c("chrom","chromStart","chromEnd","name","gieStain"),sep = "\t")
  chr <- paste("chr",chrom,sep="")
  ideos <- ideoFile[ideoFile$chrom == chr,]
  centro_idx <- which(ideos$gieStain=="acen")
  centro_coord <- range(c(ideos[centro_idx, "chromStart"], ideos[centro_idx, "chromEnd"]))
  bplim <- range(c(ideos[, "chromStart"], ideos[, "chromEnd"]))
  ideo <- ideos[-centro_idx,]
  plot(type='n', yaxt='n',xpd=FALSE,xaxt='n',xlab=' ', ylab=' ', bty='n',x=0,y=0,cex.axis=1,font.axis=2,font=2,xlim = bplim,xaxs='i',ylim=c(-0.1,1.1))
  polygon(c(centro_coord[1], centro_coord[2], centro_coord[2], centro_coord[1]), c(0, 1, 0, 1), col="darkred", border="black", lwd=2) 
  col <- character(nrow(ideo))
  col[1:length(col)] <- "#000000"
  col[grep("gneg", ideo$gieStain)] <- "gray90"
  col[grep("gpos25", ideo$gieStain)] <- "gray75"
  col[grep("gpos33", ideo$gieStain)] <- "gray66"
  col[grep("gpos50", ideo$gieStain)] <- "gray50"
  col[grep("stalk", ideo$gieStain)] <- "darkred" 
  col[grep("gpos66", ideo$gieStain)] <- "gray33"
  col[grep("gpos75", ideo$gieStain)] <- "gray25"
  for (k in 1:nrow(ideo)){
    rect(xleft=ideo[k,2],ybottom=0,xright=ideo[k,3],ytop=1,col=col[k], border=col[k])
  }
  mtext(capitalize(chr), side=3,line=2,outer=FALSE, col='black', cex=1.5,font=1,las=1,adj=0)
  p_arm <- grep("p", ideo$name)
  if(any(p_arm)){
  p_rect_coord <- range(c(ideo[p_arm, "chromStart"], 
                          ideo[p_arm, "chromEnd"]))
     
    rect(p_rect_coord[1], 0, p_rect_coord[2], 1, 
         col=NA, border="black", lwd=2)
    rm(p_rect_coord)
  }
  q_arm <- grep("q", ideo$name)
  if(any(q_arm)){
  q_rect_coord <- range(c(ideo[q_arm, "chromStart"], 
                          ideo[q_arm, "chromEnd"]))
    rect(q_rect_coord[1], 0, q_rect_coord[2], 1, col=NA , 
         border="black", lwd=2)
    rm(q_rect_coord)
  }
  xy <- par("usr")
  n <- 0
  for(i in 1:nrow(ideos)){
    pos <- (ideos[i,3]-ideos[i,2])/2+ideos[i,2]+500000
    if(n == 0){
      legend(x=pos,y=xy[3]-yinch(0.5),legend="|",xjust=1,yjust=-0.8,xpd=T,ncol=1,bty="n")
      legend(x=pos,y=xy[3]-yinch(0.5),legend=ideos[i,4],xjust=0.8,yjust=-0.5,xpd=T,cex=0.8,ncol=1,bty="n")
      n <- n + 1
    }else{
      legend(x=pos,y=xy[3]+yinch(0.5),legend="|",xjust=1,yjust=-0.54,xpd=T,ncol=1,bty="n")
      legend(x=pos,y=xy[3]+yinch(0.5),legend=ideos[i,4],xjust=0.8,yjust=-1.5,xpd=T,cex=0.8,ncol=1,bty="n")
      n <- 0
    }
  }
  
}

dbconfig <- yaml.load_file("/home/wujianming/.ssh/mysql.yaml")
server <- function(input, output, session){
    for(con in dbListConnections(MySQL())) dbDisconnect(con)
    observeEvent(input$goButton,{
            if(input$product_type == "WES Plus 8225"){
              if(input$sex == "男"){
                sampleDir <- "/share_data/wujm/project/Shiny/CNV/reference/WES/man/"
                ##
                dat <- data.table::fread(paste0(sampleDir,"baseLineCorPercent.txt"),sep='\t',header=TRUE)
                gc_file <- "/share_data/wujm/project/getFastaGC/WES.cnn.gc2.txt"
                ##连接数据库，根据检测样本，获取对应的cnn文件
                testSample <- input$sampleName
                conn <- dbConnect(MySQL(), host=dbconfig$host, dbname="nyuen_clinepilepsy", user=dbconfig$username, password=dbconfig$password)
                cmd <- paste0("select cnn from CNVBaseLine_WES where sample_name = '",testSample,"';")
                cnnFile = dbGetQuery(conn, cmd)
                newData <- data.table::fread(cnnFile$cnn,sep='\t',header=TRUE)
                names(newData) <- c("chrom","start","end","gene","depth","log2")
                # newData$percent <- newData$depth / sum(newData$depth)
                GC <- data.table::fread(gc_file,sep="\t",header=TRUE)
                mergeData <- merge(newData, GC, by=c("chrom","start","end"))
                filter0 <-  mergeData[which(depth != 0),]
                ##计算所有bin的平均覆盖度（Median）
                allMedianDepth <- median(filter0$depth)

                #计算相同GC含量的bin的覆盖度，计算这些覆盖度的median，每个GC含量梯度对应一个median值
                GC_group = group_by(filter0, gc)
                delay <- summarise(GC_group, median_sales = median(depth))
                corData = merge(filter0,delay,by='gc')
                corData$corDepth = corData$depth * (allMedianDepth / corData$median_sales)
                corData$corPercent = corData$corDepth / sum(corData$corDepth)
                corData <- corData[,c("chrom","start","end","gene","corDepth","corPercent")]


                #再将基线样本和检测样本合并
                percent_dat <- merge(corData, dat, by=c("chrom","start","end","gene"))
                percent_dat$t <- (percent_dat$corPercent - percent_dat$mean) / percent_dat$sd
                #chromStr <- input$chrom
                dataChr1 <- percent_dat[percent_dat$chrom == 'chr1',]

                output$chr1 <- renderPlotly({ggplotly(ggplot(data=dataChr1, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr1$gene,"\n","area:",dataChr1$start,"-",dataChr1$end))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,249250621) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr1_chrom <- renderPlot(chromNumber("1"))
                dataChr2 <- percent_dat[percent_dat$chrom == 'chr2',]
                output$chr2 <- renderPlotly({ggplotly(ggplot(data=dataChr2, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr2$gene,"\n","area:",dataChr2$start,"-",dataChr2$end))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,243199373) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr2_chrom <- renderPlot(chromNumber("2"))
                dataChr3 <- percent_dat[percent_dat$chrom == 'chr3',]
                output$chr3 <- renderPlotly({ggplotly(ggplot(data=dataChr3, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr3$gene,"\n","area:",dataChr3$start,"-",dataChr3$end))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,198022430) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr3_chrom <- renderPlot(chromNumber("3"))
                dataChr4 <- percent_dat[percent_dat$chrom == 'chr4',]
                output$chr4 <- renderPlotly({ggplotly(ggplot(data=dataChr4, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr4$gene,"\n","area:",dataChr4$start,"-",dataChr4$end))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,191154276) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr4_chrom <- renderPlot(chromNumber("4"))
                dataChr5 <- percent_dat[percent_dat$chrom == 'chr5',]
                output$chr5 <- renderPlotly({ggplotly(ggplot(data=dataChr5, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr5$gene,"\n","area:",dataChr5$start,"-",dataChr5$end))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,180915260) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr5_chrom <- renderPlot(chromNumber("5"))
                dataChr6 <- percent_dat[percent_dat$chrom == 'chr6',]
                output$chr6 <- renderPlotly({ggplotly(ggplot(data=dataChr6, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr6$gene,"\n","area:",dataChr6$start,"-",dataChr6$end))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,171115067) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr6_chrom <- renderPlot(chromNumber("6"))
                dataChr7 <- percent_dat[percent_dat$chrom == 'chr7',]
                output$chr7 <- renderPlotly({ggplotly(ggplot(data=dataChr7, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr7$gene,"\n","area:",dataChr7$start,"-",dataChr7$end))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,159138663) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr7_chrom <- renderPlot(chromNumber("7"))
                dataChr8 <- percent_dat[percent_dat$chrom == 'chr8',]
                output$chr8 <- renderPlotly({ggplotly(ggplot(data=dataChr8, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr8$gene,"\n","area:",dataChr8$start,"-",dataChr8$end))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,146364022) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr8_chrom <- renderPlot(chromNumber("8"))
                dataChr9 <- percent_dat[percent_dat$chrom == 'chr9',]
                output$chr9 <- renderPlotly({ggplotly(ggplot(data=dataChr9, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr9$gene,"\n","area:",dataChr9$start,"-",dataChr9$end))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,141213431) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr9_chrom <- renderPlot(chromNumber("9"))
                dataChr10 <- percent_dat[percent_dat$chrom == 'chr10',]
                output$chr10 <- renderPlotly({ggplotly(ggplot(data=dataChr10, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr10$gene,"\n","area:",dataChr10$start,"-",dataChr10$end))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,135534747) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr10_chrom <- renderPlot(chromNumber("10"))
                dataChr11 <- percent_dat[percent_dat$chrom == 'chr11',]
                output$chr11 <- renderPlotly({ggplotly(ggplot(data=dataChr11, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr11$gene,"\n","area:",dataChr11$start,"-",dataChr11$end))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,135006516) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr11_chrom <- renderPlot(chromNumber("11"))
                dataChr12 <- percent_dat[percent_dat$chrom == 'chr12',]
                output$chr12 <- renderPlotly({ggplotly(ggplot(data=dataChr12, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr12$gene,"\n","area:",dataChr12$start,"-",dataChr12$end))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,133851895) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr12_chrom <- renderPlot(chromNumber("12"))
                dataChr13 <- percent_dat[percent_dat$chrom == 'chr13',]
                output$chr13 <- renderPlotly({ggplotly(ggplot(data=dataChr13, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr13$gene,"\n","area:",dataChr13$start,"-",dataChr13$end))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,115169878) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr13_chrom <- renderPlot(chromNumber("13"))
                dataChr14 <- percent_dat[percent_dat$chrom == 'chr14',]
                output$chr14 <- renderPlotly({ggplotly(ggplot(data=dataChr14, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr14$gene,"\n","area:",dataChr14$start,"-",dataChr14$end))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,107349540) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr14_chrom <- renderPlot(chromNumber("14"))
                dataChr15 <- percent_dat[percent_dat$chrom == 'chr15',]
                output$chr15 <- renderPlotly({ggplotly(ggplot(data=dataChr15, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr15$gene,"\n","area:",dataChr15$start,"-",dataChr15$end))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,102531392) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr15_chrom <- renderPlot(chromNumber("15"))
                dataChr16 <- percent_dat[percent_dat$chrom == 'chr16',]
                output$chr16 <- renderPlotly({ggplotly(ggplot(data=dataChr16, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr16$gene,"\n","area:",dataChr16$start,"-",dataChr16$end))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,90354753) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr16_chrom <- renderPlot(chromNumber("16"))
                dataChr17 <- percent_dat[percent_dat$chrom == 'chr17',]
                output$chr17 <- renderPlotly({ggplotly(ggplot(data=dataChr17, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr17$gene,"\n","area:",dataChr17$start,"-",dataChr17$end))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,81195210) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr17_chrom <- renderPlot(chromNumber("17"))
                dataChr18 <- percent_dat[percent_dat$chrom == 'chr18',]
                output$chr18 <- renderPlotly({ggplotly(ggplot(data=dataChr18, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr18$gene,"\n","area:",dataChr18$start,"-",dataChr18$end))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,78077248) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr18_chrom <- renderPlot(chromNumber("18"))
                dataChr19 <- percent_dat[percent_dat$chrom == 'chr19',]
                output$chr19 <- renderPlotly({ggplotly(ggplot(data=dataChr19, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr19$gene,"\n","area:",dataChr19$start,"-",dataChr19$end))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,59128983) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr19_chrom <- renderPlot(chromNumber("19"))
                dataChr20 <- percent_dat[percent_dat$chrom == 'chr20',]
                output$chr20 <- renderPlotly({ggplotly(ggplot(data=dataChr20, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr20$gene,"\n","area:",dataChr20$start,"-",dataChr20$end))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,63025520) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr20_chrom <- renderPlot(chromNumber("20"))
                dataChr21 <- percent_dat[percent_dat$chrom == 'chr21',]
                output$chr21 <- renderPlotly({ggplotly(ggplot(data=dataChr21, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr21$gene,"\n","area:",dataChr21$start,"-",dataChr21$end))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,48129895) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr21_chrom <- renderPlot(chromNumber("21"))
                dataChr22 <- percent_dat[percent_dat$chrom == 'chr22',]
                output$chr22 <- renderPlotly({ggplotly(ggplot(data=dataChr22, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr22$gene,"\n","area:",dataChr22$start,"-",dataChr22$end))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,51304566) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr22_chrom <- renderPlot(chromNumber("22"))
                dataChrX <- percent_dat[percent_dat$chrom == 'chrX',]
                output$chrX <- renderPlotly({ggplotly(ggplot(data=dataChrX, mapping=aes(x=start, y=t, text = paste0("gene:",dataChrX$gene,"\n","area:",dataChrX$start,"-",dataChrX$end))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,155270560) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chrX_chrom <- renderPlot(chromNumber("X"))
                dataChrY <- percent_dat[percent_dat$chrom == 'chrY',]
                output$chrY <- renderPlotly({ggplotly(ggplot(data=dataChrY, mapping=aes(x=start, y=t, text = paste0("gene:",dataChrY$gene,"\n","area:",dataChrY$start,"-",dataChrY$end))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,59373566) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chrY_chrom <- renderPlot(chromNumber("Y"))

                ##output table
                DPPercentData <- data.table::fread(paste0(sampleDir,"baseLineCorDepthPercent.txt"),sep='\t',header=TRUE)
                MergeData <- merge(DPPercentData, corData, by=c("chrom","start","end","gene"))
                output$table <- DT::renderDataTable({MergeData})
                dbDisconnect(conn)
                #download dataset
                outFileName = paste0(testSample,"depth.csv")
                output$downloadData = downloadHandler(
                  filename = function() {
                    outFileName
                   },
                  content = function(file) {
                    write.table(MergeData, file,sep=",", row.names= FALSE)
                   }
                )
              }
              else if(input$sex == "女"){
                sampleDir <- "/share_data/wujm/project/Shiny/CNV/reference/WES/woman/"
                ##
                dat <- data.table::fread(paste0(sampleDir,"baseLineCorPercent.txt"),sep='\t',header=TRUE)
                gc_file <- "/share_data/wujm/project/getFastaGC/WES.cnn.gc2.txt"
                ##连接数据库，根据检测样本，获取对应的cnn文件
                testSample <- input$sampleName
                conn <- dbConnect(MySQL(), host=dbconfig$host, dbname="nyuen_clinepilepsy", user=dbconfig$username, password=dbconfig$password)
                cmd <- paste0("select cnn from CNVBaseLine_WES where sample_name = '",testSample,"';")
                cnnFile = dbGetQuery(conn, cmd)
                newData <- data.table::fread(cnnFile$cnn,sep='\t',header=TRUE)
                names(newData) <- c("chrom","start","end","gene","depth","log2")
                # newData$percent <- newData$depth / sum(newData$depth)
                GC <- data.table::fread(gc_file,sep="\t",header=TRUE)
                mergeData <- merge(newData, GC, by=c("chrom","start","end"))
                filter0 <-  mergeData[which(depth != 0),]
                ##计算所有bin的平均覆盖度（Median）
                allMedianDepth <- median(filter0$depth)

                #计算相同GC含量的bin的覆盖度，计算这些覆盖度的median，每个GC含量梯度对应一个median值
                GC_group = group_by(filter0, gc)
                delay <- summarise(GC_group, median_sales = median(depth))
                corData = merge(filter0,delay,by='gc')
                corData$corDepth = corData$depth * (allMedianDepth / corData$median_sales)
                corData$corPercent = corData$corDepth / sum(corData$corDepth)
                corData <- corData[,c("chrom","start","end","gene","corDepth","corPercent")]
                

                #再将基线样本和检测样本合并
                percent_dat <- merge(corData, dat, by=c("chrom","start","end","gene"))
                percent_dat$t <- (percent_dat$corPercent - percent_dat$mean) / percent_dat$sd
                #chromStr <- input$chrom
                dataChr1 <- percent_dat[percent_dat$chrom == 'chr1',]
                output$chr1 <- renderPlotly({ggplotly(ggplot(data=dataChr1, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr1$gene,"\n","area:",dataChr1$start,"-",dataChr1$end))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,249250621) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr1_chrom <- renderPlot(chromNumber("1"))
                dataChr2 <- percent_dat[percent_dat$chrom == 'chr2',]
                output$chr2 <- renderPlotly({ggplotly(ggplot(data=dataChr2, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr2$gene,"\n","area:",dataChr2$start,"-",dataChr2$end))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,243199373) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr2_chrom <- renderPlot(chromNumber("2"))
                dataChr3 <- percent_dat[percent_dat$chrom == 'chr3',]
                output$chr3 <- renderPlotly({ggplotly(ggplot(data=dataChr3, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr3$gene,"\n","area:",dataChr3$start,"-",dataChr3$end))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,198022430) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr3_chrom <- renderPlot(chromNumber("3"))
                dataChr4 <- percent_dat[percent_dat$chrom == 'chr4',]
                output$chr4 <- renderPlotly({ggplotly(ggplot(data=dataChr4, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr4$gene,"\n","area:",dataChr4$start,"-",dataChr4$end))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,191154276) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr4_chrom <- renderPlot(chromNumber("4"))
                dataChr5 <- percent_dat[percent_dat$chrom == 'chr5',]
                output$chr5 <- renderPlotly({ggplotly(ggplot(data=dataChr5, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr5$gene,"\n","area:",dataChr5$start,"-",dataChr5$end))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,180915260) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr5_chrom <- renderPlot(chromNumber("5"))
                dataChr6 <- percent_dat[percent_dat$chrom == 'chr6',]
                output$chr6 <- renderPlotly({ggplotly(ggplot(data=dataChr6, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr6$gene,"\n","area:",dataChr6$start,"-",dataChr6$end))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,171115067) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr6_chrom <- renderPlot(chromNumber("6"))
                dataChr7 <- percent_dat[percent_dat$chrom == 'chr7',]
                output$chr7 <- renderPlotly({ggplotly(ggplot(data=dataChr7, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr7$gene,"\n","area:",dataChr7$start,"-",dataChr7$end))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,159138663) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr7_chrom <- renderPlot(chromNumber("7"))
                dataChr8 <- percent_dat[percent_dat$chrom == 'chr8',]
                output$chr8 <- renderPlotly({ggplotly(ggplot(data=dataChr8, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr8$gene,"\n","area:",dataChr8$start,"-",dataChr8$end))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,146364022) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr8_chrom <- renderPlot(chromNumber("8"))
                dataChr9 <- percent_dat[percent_dat$chrom == 'chr9',]
                output$chr9 <- renderPlotly({ggplotly(ggplot(data=dataChr9, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr9$gene,"\n","area:",dataChr9$start,"-",dataChr9$end))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,141213431) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr9_chrom <- renderPlot(chromNumber("9"))
                dataChr10 <- percent_dat[percent_dat$chrom == 'chr10',]
                output$chr10 <- renderPlotly({ggplotly(ggplot(data=dataChr10, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr10$gene,"\n","area:",dataChr10$start,"-",dataChr10$end))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,135534747) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr10_chrom <- renderPlot(chromNumber("10"))
                dataChr11 <- percent_dat[percent_dat$chrom == 'chr11',]
                output$chr11 <- renderPlotly({ggplotly(ggplot(data=dataChr11, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr11$gene,"\n","area:",dataChr11$start,"-",dataChr11$end))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,135006516) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr11_chrom <- renderPlot(chromNumber("11"))
                dataChr12 <- percent_dat[percent_dat$chrom == 'chr12',]
                output$chr12 <- renderPlotly({ggplotly(ggplot(data=dataChr12, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr12$gene,"\n","area:",dataChr12$start,"-",dataChr12$end))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,133851895) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr12_chrom <- renderPlot(chromNumber("12"))
                dataChr13 <- percent_dat[percent_dat$chrom == 'chr13',]
                output$chr13 <- renderPlotly({ggplotly(ggplot(data=dataChr13, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr13$gene,"\n","area:",dataChr13$start,"-",dataChr13$end))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,115169878) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr13_chrom <- renderPlot(chromNumber("13"))
                dataChr14 <- percent_dat[percent_dat$chrom == 'chr14',]
                output$chr14 <- renderPlotly({ggplotly(ggplot(data=dataChr14, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr14$gene,"\n","area:",dataChr14$start,"-",dataChr14$end))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,107349540) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr14_chrom <- renderPlot(chromNumber("14"))
                dataChr15 <- percent_dat[percent_dat$chrom == 'chr15',]
                output$chr15 <- renderPlotly({ggplotly(ggplot(data=dataChr15, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr15$gene,"\n","area:",dataChr15$start,"-",dataChr15$end))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,102531392) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr15_chrom <- renderPlot(chromNumber("15"))
                dataChr16 <- percent_dat[percent_dat$chrom == 'chr16',]
                output$chr16 <- renderPlotly({ggplotly(ggplot(data=dataChr16, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr16$gene,"\n","area:",dataChr16$start,"-",dataChr16$end))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,90354753) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr16_chrom <- renderPlot(chromNumber("16"))
                dataChr17 <- percent_dat[percent_dat$chrom == 'chr17',]
                output$chr17 <- renderPlotly({ggplotly(ggplot(data=dataChr17, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr17$gene,"\n","area:",dataChr17$start,"-",dataChr17$end))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,81195210) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr17_chrom <- renderPlot(chromNumber("17"))
                dataChr18 <- percent_dat[percent_dat$chrom == 'chr18',]
                output$chr18 <- renderPlotly({ggplotly(ggplot(data=dataChr18, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr18$gene,"\n","area:",dataChr18$start,"-",dataChr18$end))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,78077248) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr18_chrom <- renderPlot(chromNumber("18"))
                dataChr19 <- percent_dat[percent_dat$chrom == 'chr19',]
                output$chr19 <- renderPlotly({ggplotly(ggplot(data=dataChr19, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr19$gene,"\n","area:",dataChr19$start,"-",dataChr19$end))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,59128983) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr19_chrom <- renderPlot(chromNumber("19"))
                dataChr20 <- percent_dat[percent_dat$chrom == 'chr20',]
                output$chr20 <- renderPlotly({ggplotly(ggplot(data=dataChr20, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr20$gene,"\n","area:",dataChr20$start,"-",dataChr20$end))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,63025520) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr20_chrom <- renderPlot(chromNumber("20"))
                dataChr21 <- percent_dat[percent_dat$chrom == 'chr21',]
                output$chr21 <- renderPlotly({ggplotly(ggplot(data=dataChr21, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr21$gene,"\n","area:",dataChr21$start,"-",dataChr21$end))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,48129895) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr21_chrom <- renderPlot(chromNumber("21"))
                dataChr22 <- percent_dat[percent_dat$chrom == 'chr22',]
                output$chr22 <- renderPlotly({ggplotly(ggplot(data=dataChr22, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr22$gene,"\n","area:",dataChr22$start,"-",dataChr22$end))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,51304566) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr22_chrom <- renderPlot(chromNumber("22"))
                dataChrX <- percent_dat[percent_dat$chrom == 'chrX',]
                output$chrX <- renderPlotly({ggplotly(ggplot(data=dataChrX, mapping=aes(x=start, y=t, text = paste0("gene:",dataChrX$gene,"\n","area:",dataChrX$start,"-",dataChrX$end))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,155270560) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chrX_chrom <- renderPlot(chromNumber("X"))
                dataChrY <- percent_dat[percent_dat$chrom == 'chrY',]
                output$chrY <- renderPlotly({ggplotly(ggplot(data=dataChrY, mapping=aes(x=start, y=t, text = paste0("gene:",dataChrY$gene,"\n","area:",dataChrY$start,"-",dataChrY$end))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,59373566) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chrY_chrom <- renderPlot(chromNumber("Y"))

                ##output table
                DPPercentData <- data.table::fread(paste0(sampleDir,"baseLineCorDepthPercent.txt"),sep='\t',header=TRUE)
                MergeData <- merge(DPPercentData, corData, by=c("chrom","start","end","gene"))
                output$table <- DT::renderDataTable({MergeData})
                dbDisconnect(conn)

                outFileName = paste0(testSample,"depth.csv")
                output$downloadData = downloadHandler(
                  filename = function() {
                    outFileName
                   },
                  content = function(file) {
                    write.table(MergeData, file,sep=",", row.names= FALSE)
                   }
                )
              }
            }
            else if(input$product_type == "WES Plus 8226"){
              if(input$sex == "男"){
                sampleDir <- "/share_data/wujm/project/Shiny/CNV/reference_V2/WES/man/"
                ##
                dat <- data.table::fread(paste0(sampleDir,"baseLineCorPercent.txt"),sep='\t',header=TRUE)
                gc_file <- "/share_data/wujm/project/getFastaGC/WES.cnn.gc2.txt"
                ##连接数据库，根据检测样本，获取对应的cnn文件
                testSample <- input$sampleName
                conn <- dbConnect(MySQL(), host=dbconfig$host, dbname="nyuen_clinepilepsy", user=dbconfig$username, password=dbconfig$password)
                cmd <- paste0("select cnn from CNVBaseLine_WES where sample_name = '",testSample,"';")
                cnnFile = dbGetQuery(conn, cmd)
                newData <- data.table::fread(cnnFile$cnn,sep='\t',header=TRUE)
                names(newData) <- c("chrom","start","end","gene","depth","log2")
                # newData$percent <- newData$depth / sum(newData$depth)
                GC <- data.table::fread(gc_file,sep="\t",header=TRUE)
                mergeData <- merge(newData, GC, by=c("chrom","start","end"))
                filter0 <-  mergeData[which(depth != 0),]
                ##计算所有bin的平均覆盖度（Median）
                allMedianDepth <- median(filter0$depth)

                #计算相同GC含量的bin的覆盖度，计算这些覆盖度的median，每个GC含量梯度对应一个median值
                GC_group = group_by(filter0, gc)
                delay <- summarise(GC_group, median_sales = median(depth))
                corData = merge(filter0,delay,by='gc')
                corData$corDepth = corData$depth * (allMedianDepth / corData$median_sales)
                corData$corPercent = corData$corDepth / sum(corData$corDepth)
                corData <- corData[,c("chrom","start","end","gene","corDepth","corPercent")]


                #再将基线样本和检测样本合并
                percent_dat <- merge(corData, dat, by=c("chrom","start","end","gene"))
                percent_dat$t <- (percent_dat$corPercent - percent_dat$mean) / percent_dat$sd
                #chromStr <- input$chrom
                dataChr1 <- percent_dat[percent_dat$chrom == 'chr1',]

                output$chr1 <- renderPlotly({ggplotly(ggplot(data=dataChr1, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr1$gene))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,249250621) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr1_chrom <- renderPlot(chromNumber("1"))
                dataChr2 <- percent_dat[percent_dat$chrom == 'chr2',]
                output$chr2 <- renderPlotly({ggplotly(ggplot(data=dataChr2, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr2$gene))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,243199373) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr2_chrom <- renderPlot(chromNumber("2"))
                dataChr3 <- percent_dat[percent_dat$chrom == 'chr3',]
                output$chr3 <- renderPlotly({ggplotly(ggplot(data=dataChr3, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr3$gene))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,198022430) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr3_chrom <- renderPlot(chromNumber("3"))
                dataChr4 <- percent_dat[percent_dat$chrom == 'chr4',]
                output$chr4 <- renderPlotly({ggplotly(ggplot(data=dataChr4, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr4$gene))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,191154276) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr4_chrom <- renderPlot(chromNumber("4"))
                dataChr5 <- percent_dat[percent_dat$chrom == 'chr5',]
                output$chr5 <- renderPlotly({ggplotly(ggplot(data=dataChr5, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr5$gene))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,180915260) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr5_chrom <- renderPlot(chromNumber("5"))
                dataChr6 <- percent_dat[percent_dat$chrom == 'chr6',]
                output$chr6 <- renderPlotly({ggplotly(ggplot(data=dataChr6, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr6$gene))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,171115067) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr6_chrom <- renderPlot(chromNumber("6"))
                dataChr7 <- percent_dat[percent_dat$chrom == 'chr7',]
                output$chr7 <- renderPlotly({ggplotly(ggplot(data=dataChr7, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr7$gene))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,159138663) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr7_chrom <- renderPlot(chromNumber("7"))
                dataChr8 <- percent_dat[percent_dat$chrom == 'chr8',]
                output$chr8 <- renderPlotly({ggplotly(ggplot(data=dataChr8, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr8$gene))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,146364022) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr8_chrom <- renderPlot(chromNumber("8"))
                dataChr9 <- percent_dat[percent_dat$chrom == 'chr9',]
                output$chr9 <- renderPlotly({ggplotly(ggplot(data=dataChr9, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr9$gene))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,141213431) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr9_chrom <- renderPlot(chromNumber("9"))
                dataChr10 <- percent_dat[percent_dat$chrom == 'chr10',]
                output$chr10 <- renderPlotly({ggplotly(ggplot(data=dataChr10, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr10$gene))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,135534747) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr10_chrom <- renderPlot(chromNumber("10"))
                dataChr11 <- percent_dat[percent_dat$chrom == 'chr11',]
                output$chr11 <- renderPlotly({ggplotly(ggplot(data=dataChr11, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr11$gene))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,135006516) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr11_chrom <- renderPlot(chromNumber("11"))
                dataChr12 <- percent_dat[percent_dat$chrom == 'chr12',]
                output$chr12 <- renderPlotly({ggplotly(ggplot(data=dataChr12, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr12$gene))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,133851895) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr12_chrom <- renderPlot(chromNumber("12"))
                dataChr13 <- percent_dat[percent_dat$chrom == 'chr13',]
                output$chr13 <- renderPlotly({ggplotly(ggplot(data=dataChr13, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr13$gene))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,115169878) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr13_chrom <- renderPlot(chromNumber("13"))
                dataChr14 <- percent_dat[percent_dat$chrom == 'chr14',]
                output$chr14 <- renderPlotly({ggplotly(ggplot(data=dataChr14, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr14$gene))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,107349540) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr14_chrom <- renderPlot(chromNumber("14"))
                dataChr15 <- percent_dat[percent_dat$chrom == 'chr15',]
                output$chr15 <- renderPlotly({ggplotly(ggplot(data=dataChr15, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr15$gene))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,102531392) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr15_chrom <- renderPlot(chromNumber("15"))
                dataChr16 <- percent_dat[percent_dat$chrom == 'chr16',]
                output$chr16 <- renderPlotly({ggplotly(ggplot(data=dataChr16, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr16$gene))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,90354753) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr16_chrom <- renderPlot(chromNumber("16"))
                dataChr17 <- percent_dat[percent_dat$chrom == 'chr17',]
                output$chr17 <- renderPlotly({ggplotly(ggplot(data=dataChr17, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr17$gene))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,81195210) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr17_chrom <- renderPlot(chromNumber("17"))
                dataChr18 <- percent_dat[percent_dat$chrom == 'chr18',]
                output$chr18 <- renderPlotly({ggplotly(ggplot(data=dataChr18, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr18$gene))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,78077248) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr18_chrom <- renderPlot(chromNumber("18"))
                dataChr19 <- percent_dat[percent_dat$chrom == 'chr19',]
                output$chr19 <- renderPlotly({ggplotly(ggplot(data=dataChr19, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr19$gene))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,59128983) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr19_chrom <- renderPlot(chromNumber("19"))
                dataChr20 <- percent_dat[percent_dat$chrom == 'chr20',]
                output$chr20 <- renderPlotly({ggplotly(ggplot(data=dataChr20, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr20$gene))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,63025520) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr20_chrom <- renderPlot(chromNumber("20"))
                dataChr21 <- percent_dat[percent_dat$chrom == 'chr21',]
                output$chr21 <- renderPlotly({ggplotly(ggplot(data=dataChr21, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr21$gene))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,48129895) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr21_chrom <- renderPlot(chromNumber("21"))
                dataChr22 <- percent_dat[percent_dat$chrom == 'chr22',]
                output$chr22 <- renderPlotly({ggplotly(ggplot(data=dataChr22, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr22$gene))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,51304566) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr22_chrom <- renderPlot(chromNumber("22"))
                dataChrX <- percent_dat[percent_dat$chrom == 'chrX',]
                output$chrX <- renderPlotly({ggplotly(ggplot(data=dataChrX, mapping=aes(x=start, y=t, text = paste0("gene:",dataChrX$gene))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,155270560) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chrX_chrom <- renderPlot(chromNumber("X"))
                dataChrY <- percent_dat[percent_dat$chrom == 'chrY',]
                output$chrY <- renderPlotly({ggplotly(ggplot(data=dataChrY, mapping=aes(x=start, y=t, text = paste0("gene:",dataChrY$gene))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,59373566) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chrY_chrom <- renderPlot(chromNumber("Y"))

                ##output table
                DPPercentData <- data.table::fread(paste0(sampleDir,"baseLineCorDepthPercent.txt"),sep='\t',header=TRUE)
                MergeData <- merge(DPPercentData, corData, by=c("chrom","start","end","gene"))
                output$table <- DT::renderDataTable({MergeData})
                dbDisconnect(conn)

                outFileName = paste0(testSample,".depth.csv")
                output$downloadData = downloadHandler(
                  filename = function() {
                    outFileName
                   },
                  content = function(file) {
                    write.table(MergeData, file,sep=",", row.names= FALSE)
                   }
                )
              }else if(input$sex == "女"){
                sampleDir <- "/share_data/wujm/project/Shiny/CNV/reference_V2/WES/woman/"
                ##
                dat <- data.table::fread(paste0(sampleDir,"baseLineCorPercent.txt"),sep='\t',header=TRUE)
                gc_file <- "/share_data/wujm/project/getFastaGC/WES.cnn.gc2.txt"
                ##连接数据库，根据检测样本，获取对应的cnn文件
                testSample <- input$sampleName
                conn <- dbConnect(MySQL(), host=dbconfig$host, dbname="nyuen_clinepilepsy", user=dbconfig$username, password=dbconfig$password)
                cmd <- paste0("select cnn from CNVBaseLine_WES where sample_name = '",testSample,"';")
                cnnFile = dbGetQuery(conn, cmd)
                newData <- data.table::fread(cnnFile$cnn,sep='\t',header=TRUE)
                names(newData) <- c("chrom","start","end","gene","depth","log2")
                # newData$percent <- newData$depth / sum(newData$depth)
                GC <- data.table::fread(gc_file,sep="\t",header=TRUE)
                mergeData <- merge(newData, GC, by=c("chrom","start","end"))
                filter0 <-  mergeData[which(depth != 0),]
                ##计算所有bin的平均覆盖度（Median）
                allMedianDepth <- median(filter0$depth)

                #计算相同GC含量的bin的覆盖度，计算这些覆盖度的median，每个GC含量梯度对应一个median值
                GC_group = group_by(filter0, gc)
                delay <- summarise(GC_group, median_sales = median(depth))
                corData = merge(filter0,delay,by='gc')
                corData$corDepth = corData$depth * (allMedianDepth / corData$median_sales)
                corData$corPercent = corData$corDepth / sum(corData$corDepth)
                corData <- corData[,c("chrom","start","end","gene","corDepth","corPercent")]
                

                #再将基线样本和检测样本合并
                percent_dat <- merge(corData, dat, by=c("chrom","start","end","gene"))
                percent_dat$t <- (percent_dat$corPercent - percent_dat$mean) / percent_dat$sd
                #chromStr <- input$chrom
                dataChr1 <- percent_dat[percent_dat$chrom == 'chr1',]
                output$chr1 <- renderPlotly({ggplotly(ggplot(data=dataChr1, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr1$gene))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,249250621) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr1_chrom <- renderPlot(chromNumber("1"))
                dataChr2 <- percent_dat[percent_dat$chrom == 'chr2',]
                output$chr2 <- renderPlotly({ggplotly(ggplot(data=dataChr2, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr2$gene))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,243199373) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr2_chrom <- renderPlot(chromNumber("2"))
                dataChr3 <- percent_dat[percent_dat$chrom == 'chr3',]
                output$chr3 <- renderPlotly({ggplotly(ggplot(data=dataChr3, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr3$gene))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,198022430) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr3_chrom <- renderPlot(chromNumber("3"))
                dataChr4 <- percent_dat[percent_dat$chrom == 'chr4',]
                output$chr4 <- renderPlotly({ggplotly(ggplot(data=dataChr4, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr4$gene))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,191154276) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr4_chrom <- renderPlot(chromNumber("4"))
                dataChr5 <- percent_dat[percent_dat$chrom == 'chr5',]
                output$chr5 <- renderPlotly({ggplotly(ggplot(data=dataChr5, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr5$gene))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,180915260) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr5_chrom <- renderPlot(chromNumber("5"))
                dataChr6 <- percent_dat[percent_dat$chrom == 'chr6',]
                output$chr6 <- renderPlotly({ggplotly(ggplot(data=dataChr6, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr6$gene))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,171115067) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr6_chrom <- renderPlot(chromNumber("6"))
                dataChr7 <- percent_dat[percent_dat$chrom == 'chr7',]
                output$chr7 <- renderPlotly({ggplotly(ggplot(data=dataChr7, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr7$gene))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,159138663) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr7_chrom <- renderPlot(chromNumber("7"))
                dataChr8 <- percent_dat[percent_dat$chrom == 'chr8',]
                output$chr8 <- renderPlotly({ggplotly(ggplot(data=dataChr8, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr8$gene))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,146364022) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr8_chrom <- renderPlot(chromNumber("8"))
                dataChr9 <- percent_dat[percent_dat$chrom == 'chr9',]
                output$chr9 <- renderPlotly({ggplotly(ggplot(data=dataChr9, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr9$gene))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,141213431) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr9_chrom <- renderPlot(chromNumber("9"))
                dataChr10 <- percent_dat[percent_dat$chrom == 'chr10',]
                output$chr10 <- renderPlotly({ggplotly(ggplot(data=dataChr10, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr10$gene))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,135534747) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr10_chrom <- renderPlot(chromNumber("10"))
                dataChr11 <- percent_dat[percent_dat$chrom == 'chr11',]
                output$chr11 <- renderPlotly({ggplotly(ggplot(data=dataChr11, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr11$gene))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,135006516) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr11_chrom <- renderPlot(chromNumber("11"))
                dataChr12 <- percent_dat[percent_dat$chrom == 'chr12',]
                output$chr12 <- renderPlotly({ggplotly(ggplot(data=dataChr12, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr12$gene))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,133851895) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr12_chrom <- renderPlot(chromNumber("12"))
                dataChr13 <- percent_dat[percent_dat$chrom == 'chr13',]
                output$chr13 <- renderPlotly({ggplotly(ggplot(data=dataChr13, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr13$gene))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,115169878) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr13_chrom <- renderPlot(chromNumber("13"))
                dataChr14 <- percent_dat[percent_dat$chrom == 'chr14',]
                output$chr14 <- renderPlotly({ggplotly(ggplot(data=dataChr14, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr14$gene))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,107349540) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr14_chrom <- renderPlot(chromNumber("14"))
                dataChr15 <- percent_dat[percent_dat$chrom == 'chr15',]
                output$chr15 <- renderPlotly({ggplotly(ggplot(data=dataChr15, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr15$gene))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,102531392) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr15_chrom <- renderPlot(chromNumber("15"))
                dataChr16 <- percent_dat[percent_dat$chrom == 'chr16',]
                output$chr16 <- renderPlotly({ggplotly(ggplot(data=dataChr16, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr16$gene))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,90354753) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr16_chrom <- renderPlot(chromNumber("16"))
                dataChr17 <- percent_dat[percent_dat$chrom == 'chr17',]
                output$chr17 <- renderPlotly({ggplotly(ggplot(data=dataChr17, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr17$gene))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,90354753) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr17_chrom <- renderPlot(chromNumber("17"))
                dataChr18 <- percent_dat[percent_dat$chrom == 'chr18',]
                output$chr18 <- renderPlotly({ggplotly(ggplot(data=dataChr18, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr18$gene))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,90354753) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr18_chrom <- renderPlot(chromNumber("18"))
                dataChr19 <- percent_dat[percent_dat$chrom == 'chr19',]
                output$chr19 <- renderPlotly({ggplotly(ggplot(data=dataChr19, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr19$gene))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,59128983) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr19_chrom <- renderPlot(chromNumber("19"))
                dataChr20 <- percent_dat[percent_dat$chrom == 'chr20',]
                output$chr20 <- renderPlotly({ggplotly(ggplot(data=dataChr20, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr20$gene))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,63025520) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr20_chrom <- renderPlot(chromNumber("20"))
                dataChr21 <- percent_dat[percent_dat$chrom == 'chr21',]
                output$chr21 <- renderPlotly({ggplotly(ggplot(data=dataChr21, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr21$gene))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,48129895) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr21_chrom <- renderPlot(chromNumber("21"))
                dataChr22 <- percent_dat[percent_dat$chrom == 'chr22',]
                output$chr22 <- renderPlotly({ggplotly(ggplot(data=dataChr22, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr22$gene))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,51304566) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr22_chrom <- renderPlot(chromNumber("22"))
                dataChrX <- percent_dat[percent_dat$chrom == 'chrX',]
                output$chrX <- renderPlotly({ggplotly(ggplot(data=dataChrX, mapping=aes(x=start, y=t, text = paste0("gene:",dataChrX$gene))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,155270560) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chrX_chrom <- renderPlot(chromNumber("X"))
                dataChrY <- percent_dat[percent_dat$chrom == 'chrY',]
                output$chrY <- renderPlotly({ggplotly(ggplot(data=dataChrY, mapping=aes(x=start, y=t, text = paste0("gene:",dataChrY$gene))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,59373566) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chrY_chrom <- renderPlot(chromNumber("Y"))

                ##output table
                DPPercentData <- data.table::fread(paste0(sampleDir,"baseLineCorDepthPercent.txt"),sep='\t',header=TRUE)
                MergeData <- merge(DPPercentData, corData, by=c("chrom","start","end","gene"))
                output$table <- DT::renderDataTable({MergeData})
                dbDisconnect(conn)

                outFileName = paste0(testSample,".depth.csv")
                output$downloadData = downloadHandler(
                  filename = function() {
                    outFileName
                   },
                  content = function(file) {
                    write.table(MergeData, file,sep=",", row.names= FALSE)
                   }
                )
              }
            }
            else if(input$product_type == "JM_WES_Plus"){
              if(input$sex == "男"){
                sampleDir <- "/share_data/wujm/project/Shiny/CNV/reference_JM/man/"
                ##
                dat <- data.table::fread(paste0(sampleDir,"baseLineCorPercent.txt"),sep='\t',header=TRUE)
                gc_file <- "/share_data/wujm/project/getFastaGC/WES.cnn.gc2.txt"
                ##连接数据库，根据检测样本，获取对应的cnn文件
                testSample <- input$sampleName
                conn <- dbConnect(MySQL(), host=dbconfig$host, dbname="nyuen_clinepilepsy", user=dbconfig$username, password=dbconfig$password)
                cmd <- paste0("select cnn from CNVBaseLine_JM_WP where sample_name = '",testSample,"';")
                cnnFile = dbGetQuery(conn, cmd)
                newData <- data.table::fread(cnnFile$cnn,sep='\t',header=TRUE)
                names(newData) <- c("chrom","start","end","gene","depth","log2")
                # newData$percent <- newData$depth / sum(newData$depth)
                GC <- data.table::fread(gc_file,sep="\t",header=TRUE)
                mergeData <- merge(newData, GC, by=c("chrom","start","end"))
                filter0 <-  mergeData[which(depth != 0),]
                ##计算所有bin的平均覆盖度（Median）
                allMedianDepth <- median(filter0$depth)

                #计算相同GC含量的bin的覆盖度，计算这些覆盖度的median，每个GC含量梯度对应一个median值
                GC_group = group_by(filter0, gc)
                delay <- summarise(GC_group, median_sales = median(depth))
                corData = merge(filter0,delay,by='gc')
                corData$corDepth = corData$depth * (allMedianDepth / corData$median_sales)
                corData$corPercent = corData$corDepth / sum(corData$corDepth)
                corData <- corData[,c("chrom","start","end","gene","corDepth","corPercent")]


                #再将基线样本和检测样本合并
                percent_dat <- merge(corData, dat, by=c("chrom","start","end","gene"))
                percent_dat$t <- (percent_dat$corPercent - percent_dat$mean) / percent_dat$sd
                #chromStr <- input$chrom
                dataChr1 <- percent_dat[percent_dat$chrom == 'chr1',]

                output$chr1 <- renderPlotly({ggplotly(ggplot(data=dataChr1, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr1$gene))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,249250621) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr1_chrom <- renderPlot(chromNumber("1"))
                dataChr2 <- percent_dat[percent_dat$chrom == 'chr2',]
                output$chr2 <- renderPlotly({ggplotly(ggplot(data=dataChr2, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr2$gene))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,243199373) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr2_chrom <- renderPlot(chromNumber("2"))
                dataChr3 <- percent_dat[percent_dat$chrom == 'chr3',]
                output$chr3 <- renderPlotly({ggplotly(ggplot(data=dataChr3, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr3$gene))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,198022430) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr3_chrom <- renderPlot(chromNumber("3"))
                dataChr4 <- percent_dat[percent_dat$chrom == 'chr4',]
                output$chr4 <- renderPlotly({ggplotly(ggplot(data=dataChr4, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr4$gene))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,191154276) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr4_chrom <- renderPlot(chromNumber("4"))
                dataChr5 <- percent_dat[percent_dat$chrom == 'chr5',]
                output$chr5 <- renderPlotly({ggplotly(ggplot(data=dataChr5, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr5$gene))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,180915260) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr5_chrom <- renderPlot(chromNumber("5"))
                dataChr6 <- percent_dat[percent_dat$chrom == 'chr6',]
                output$chr6 <- renderPlotly({ggplotly(ggplot(data=dataChr6, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr6$gene))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,171115067) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr6_chrom <- renderPlot(chromNumber("6"))
                dataChr7 <- percent_dat[percent_dat$chrom == 'chr7',]
                output$chr7 <- renderPlotly({ggplotly(ggplot(data=dataChr7, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr7$gene))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,159138663) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr7_chrom <- renderPlot(chromNumber("7"))
                dataChr8 <- percent_dat[percent_dat$chrom == 'chr8',]
                output$chr8 <- renderPlotly({ggplotly(ggplot(data=dataChr8, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr8$gene))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,146364022) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr8_chrom <- renderPlot(chromNumber("8"))
                dataChr9 <- percent_dat[percent_dat$chrom == 'chr9',]
                output$chr9 <- renderPlotly({ggplotly(ggplot(data=dataChr9, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr9$gene))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,141213431) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr9_chrom <- renderPlot(chromNumber("9"))
                dataChr10 <- percent_dat[percent_dat$chrom == 'chr10',]
                output$chr10 <- renderPlotly({ggplotly(ggplot(data=dataChr10, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr10$gene))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,135534747) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr10_chrom <- renderPlot(chromNumber("10"))
                dataChr11 <- percent_dat[percent_dat$chrom == 'chr11',]
                output$chr11 <- renderPlotly({ggplotly(ggplot(data=dataChr11, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr11$gene))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,135006516) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr11_chrom <- renderPlot(chromNumber("11"))
                dataChr12 <- percent_dat[percent_dat$chrom == 'chr12',]
                output$chr12 <- renderPlotly({ggplotly(ggplot(data=dataChr12, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr12$gene))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,133851895) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr12_chrom <- renderPlot(chromNumber("12"))
                dataChr13 <- percent_dat[percent_dat$chrom == 'chr13',]
                output$chr13 <- renderPlotly({ggplotly(ggplot(data=dataChr13, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr13$gene))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,115169878) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr13_chrom <- renderPlot(chromNumber("13"))
                dataChr14 <- percent_dat[percent_dat$chrom == 'chr14',]
                output$chr14 <- renderPlotly({ggplotly(ggplot(data=dataChr14, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr14$gene))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,107349540) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr14_chrom <- renderPlot(chromNumber("14"))
                dataChr15 <- percent_dat[percent_dat$chrom == 'chr15',]
                output$chr15 <- renderPlotly({ggplotly(ggplot(data=dataChr15, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr15$gene))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,102531392) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr15_chrom <- renderPlot(chromNumber("15"))
                dataChr16 <- percent_dat[percent_dat$chrom == 'chr16',]
                output$chr16 <- renderPlotly({ggplotly(ggplot(data=dataChr16, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr16$gene))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,90354753) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr16_chrom <- renderPlot(chromNumber("16"))
                dataChr17 <- percent_dat[percent_dat$chrom == 'chr17',]
                output$chr17 <- renderPlotly({ggplotly(ggplot(data=dataChr17, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr17$gene))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,81195210) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr17_chrom <- renderPlot(chromNumber("17"))
                dataChr18 <- percent_dat[percent_dat$chrom == 'chr18',]
                output$chr18 <- renderPlotly({ggplotly(ggplot(data=dataChr18, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr18$gene))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,78077248) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr18_chrom <- renderPlot(chromNumber("18"))
                dataChr19 <- percent_dat[percent_dat$chrom == 'chr19',]
                output$chr19 <- renderPlotly({ggplotly(ggplot(data=dataChr19, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr19$gene))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,59128983) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr19_chrom <- renderPlot(chromNumber("19"))
                dataChr20 <- percent_dat[percent_dat$chrom == 'chr20',]
                output$chr20 <- renderPlotly({ggplotly(ggplot(data=dataChr20, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr20$gene))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,63025520) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr20_chrom <- renderPlot(chromNumber("20"))
                dataChr21 <- percent_dat[percent_dat$chrom == 'chr21',]
                output$chr21 <- renderPlotly({ggplotly(ggplot(data=dataChr21, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr21$gene))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,48129895) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr21_chrom <- renderPlot(chromNumber("21"))
                dataChr22 <- percent_dat[percent_dat$chrom == 'chr22',]
                output$chr22 <- renderPlotly({ggplotly(ggplot(data=dataChr22, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr22$gene))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,51304566) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr22_chrom <- renderPlot(chromNumber("22"))
                dataChrX <- percent_dat[percent_dat$chrom == 'chrX',]
                output$chrX <- renderPlotly({ggplotly(ggplot(data=dataChrX, mapping=aes(x=start, y=t, text = paste0("gene:",dataChrX$gene))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,155270560) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chrX_chrom <- renderPlot(chromNumber("X"))
                dataChrY <- percent_dat[percent_dat$chrom == 'chrY',]
                output$chrY <- renderPlotly({ggplotly(ggplot(data=dataChrY, mapping=aes(x=start, y=t, text = paste0("gene:",dataChrY$gene))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,59373566) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chrY_chrom <- renderPlot(chromNumber("Y"))

                ##output table
                DPPercentData <- data.table::fread(paste0(sampleDir,"baseLineCorDepthPercent.txt"),sep='\t',header=TRUE)
                MergeData <- merge(DPPercentData, corData, by=c("chrom","start","end","gene"))
                output$table <- DT::renderDataTable({MergeData})
                dbDisconnect(conn)

                outFileName = paste0(testSample,".depth.csv")
                output$downloadData = downloadHandler(
                  filename = function() {
                    outFileName
                   },
                  content = function(file) {
                    write.table(MergeData, file,sep=",", row.names= FALSE)
                   }
                )
              }else if(input$sex == "女"){
                sampleDir <- "/share_data/wujm/project/Shiny/CNV/reference_JM/woman/"
                ##
                dat <- data.table::fread(paste0(sampleDir,"baseLineCorPercent.txt"),sep='\t',header=TRUE)
                gc_file <- "/share_data/wujm/project/getFastaGC/WES.cnn.gc2.txt"
                ##连接数据库，根据检测样本，获取对应的cnn文件
                testSample <- input$sampleName
                conn <- dbConnect(MySQL(), host=dbconfig$host, dbname="nyuen_clinepilepsy", user=dbconfig$username, password=dbconfig$password)
                cmd <- paste0("select cnn from CNVBaseLine_JM_WP where sample_name = '",testSample,"';")
                cnnFile = dbGetQuery(conn, cmd)
                newData <- data.table::fread(cnnFile$cnn,sep='\t',header=TRUE)
                names(newData) <- c("chrom","start","end","gene","depth","log2")
                # newData$percent <- newData$depth / sum(newData$depth)
                GC <- data.table::fread(gc_file,sep="\t",header=TRUE)
                mergeData <- merge(newData, GC, by=c("chrom","start","end"))
                filter0 <-  mergeData[which(depth != 0),]
                ##计算所有bin的平均覆盖度（Median）
                allMedianDepth <- median(filter0$depth)

                #计算相同GC含量的bin的覆盖度，计算这些覆盖度的median，每个GC含量梯度对应一个median值
                GC_group = group_by(filter0, gc)
                delay <- summarise(GC_group, median_sales = median(depth))
                corData = merge(filter0,delay,by='gc')
                corData$corDepth = corData$depth * (allMedianDepth / corData$median_sales)
                corData$corPercent = corData$corDepth / sum(corData$corDepth)
                corData <- corData[,c("chrom","start","end","gene","corDepth","corPercent")]
                

                #再将基线样本和检测样本合并
                percent_dat <- merge(corData, dat, by=c("chrom","start","end","gene"))
                percent_dat$t <- (percent_dat$corPercent - percent_dat$mean) / percent_dat$sd
                #chromStr <- input$chrom
                dataChr1 <- percent_dat[percent_dat$chrom == 'chr1',]
                output$chr1 <- renderPlotly({ggplotly(ggplot(data=dataChr1, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr1$gene))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,249250621) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr1_chrom <- renderPlot(chromNumber("1"))
                dataChr2 <- percent_dat[percent_dat$chrom == 'chr2',]
                output$chr2 <- renderPlotly({ggplotly(ggplot(data=dataChr2, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr2$gene))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,243199373) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr2_chrom <- renderPlot(chromNumber("2"))
                dataChr3 <- percent_dat[percent_dat$chrom == 'chr3',]
                output$chr3 <- renderPlotly({ggplotly(ggplot(data=dataChr3, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr3$gene))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,198022430) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr3_chrom <- renderPlot(chromNumber("3"))
                dataChr4 <- percent_dat[percent_dat$chrom == 'chr4',]
                output$chr4 <- renderPlotly({ggplotly(ggplot(data=dataChr4, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr4$gene))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,191154276) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr4_chrom <- renderPlot(chromNumber("4"))
                dataChr5 <- percent_dat[percent_dat$chrom == 'chr5',]
                output$chr5 <- renderPlotly({ggplotly(ggplot(data=dataChr5, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr5$gene))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,180915260) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr5_chrom <- renderPlot(chromNumber("5"))
                dataChr6 <- percent_dat[percent_dat$chrom == 'chr6',]
                output$chr6 <- renderPlotly({ggplotly(ggplot(data=dataChr6, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr6$gene))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,171115067) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr6_chrom <- renderPlot(chromNumber("6"))
                dataChr7 <- percent_dat[percent_dat$chrom == 'chr7',]
                output$chr7 <- renderPlotly({ggplotly(ggplot(data=dataChr7, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr7$gene))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,159138663) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr7_chrom <- renderPlot(chromNumber("7"))
                dataChr8 <- percent_dat[percent_dat$chrom == 'chr8',]
                output$chr8 <- renderPlotly({ggplotly(ggplot(data=dataChr8, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr8$gene))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,146364022) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr8_chrom <- renderPlot(chromNumber("8"))
                dataChr9 <- percent_dat[percent_dat$chrom == 'chr9',]
                output$chr9 <- renderPlotly({ggplotly(ggplot(data=dataChr9, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr9$gene))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,141213431) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr9_chrom <- renderPlot(chromNumber("9"))
                dataChr10 <- percent_dat[percent_dat$chrom == 'chr10',]
                output$chr10 <- renderPlotly({ggplotly(ggplot(data=dataChr10, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr10$gene))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,135534747) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr10_chrom <- renderPlot(chromNumber("10"))
                dataChr11 <- percent_dat[percent_dat$chrom == 'chr11',]
                output$chr11 <- renderPlotly({ggplotly(ggplot(data=dataChr11, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr11$gene))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,135006516) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr11_chrom <- renderPlot(chromNumber("11"))
                dataChr12 <- percent_dat[percent_dat$chrom == 'chr12',]
                output$chr12 <- renderPlotly({ggplotly(ggplot(data=dataChr12, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr12$gene))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,133851895) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr12_chrom <- renderPlot(chromNumber("12"))
                dataChr13 <- percent_dat[percent_dat$chrom == 'chr13',]
                output$chr13 <- renderPlotly({ggplotly(ggplot(data=dataChr13, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr13$gene))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,115169878) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr13_chrom <- renderPlot(chromNumber("13"))
                dataChr14 <- percent_dat[percent_dat$chrom == 'chr14',]
                output$chr14 <- renderPlotly({ggplotly(ggplot(data=dataChr14, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr14$gene))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,107349540) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr14_chrom <- renderPlot(chromNumber("14"))
                dataChr15 <- percent_dat[percent_dat$chrom == 'chr15',]
                output$chr15 <- renderPlotly({ggplotly(ggplot(data=dataChr15, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr15$gene))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,102531392) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr15_chrom <- renderPlot(chromNumber("15"))
                dataChr16 <- percent_dat[percent_dat$chrom == 'chr16',]
                output$chr16 <- renderPlotly({ggplotly(ggplot(data=dataChr16, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr16$gene))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,90354753) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr16_chrom <- renderPlot(chromNumber("16"))
                dataChr17 <- percent_dat[percent_dat$chrom == 'chr17',]
                output$chr17 <- renderPlotly({ggplotly(ggplot(data=dataChr17, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr17$gene))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,90354753) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr17_chrom <- renderPlot(chromNumber("17"))
                dataChr18 <- percent_dat[percent_dat$chrom == 'chr18',]
                output$chr18 <- renderPlotly({ggplotly(ggplot(data=dataChr18, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr18$gene))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,90354753) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr18_chrom <- renderPlot(chromNumber("18"))
                dataChr19 <- percent_dat[percent_dat$chrom == 'chr19',]
                output$chr19 <- renderPlotly({ggplotly(ggplot(data=dataChr19, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr19$gene))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,59128983) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr19_chrom <- renderPlot(chromNumber("19"))
                dataChr20 <- percent_dat[percent_dat$chrom == 'chr20',]
                output$chr20 <- renderPlotly({ggplotly(ggplot(data=dataChr20, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr20$gene))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,63025520) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr20_chrom <- renderPlot(chromNumber("20"))
                dataChr21 <- percent_dat[percent_dat$chrom == 'chr21',]
                output$chr21 <- renderPlotly({ggplotly(ggplot(data=dataChr21, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr21$gene))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,48129895) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr21_chrom <- renderPlot(chromNumber("21"))
                dataChr22 <- percent_dat[percent_dat$chrom == 'chr22',]
                output$chr22 <- renderPlotly({ggplotly(ggplot(data=dataChr22, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr22$gene))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,51304566) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr22_chrom <- renderPlot(chromNumber("22"))
                dataChrX <- percent_dat[percent_dat$chrom == 'chrX',]
                output$chrX <- renderPlotly({ggplotly(ggplot(data=dataChrX, mapping=aes(x=start, y=t, text = paste0("gene:",dataChrX$gene))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,155270560) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chrX_chrom <- renderPlot(chromNumber("X"))
                dataChrY <- percent_dat[percent_dat$chrom == 'chrY',]
                output$chrY <- renderPlotly({ggplotly(ggplot(data=dataChrY, mapping=aes(x=start, y=t, text = paste0("gene:",dataChrY$gene))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,59373566) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chrY_chrom <- renderPlot(chromNumber("Y"))

                ##output table
                DPPercentData <- data.table::fread(paste0(sampleDir,"baseLineCorDepthPercent.txt"),sep='\t',header=TRUE)
                MergeData <- merge(DPPercentData, corData, by=c("chrom","start","end","gene"))
                output$table <- DT::renderDataTable({MergeData})
                dbDisconnect(conn)

                outFileName = paste0(testSample,".depth.csv")
                output$downloadData = downloadHandler(
                  filename = function() {
                    outFileName
                   },
                  content = function(file) {
                    write.table(MergeData, file,sep=",", row.names= FALSE)
                   }
                )
              }
            }
            else if(input$product_type == "单人WES/Trio WES"){
              if(input$sex == "男"){
                sampleDir <- "/share_data/wujm/project/Shiny/CNV/reference_T192V/man/"
                ##
                dat <- data.table::fread(paste0(sampleDir,"baseLineCorPercent.txt"),sep='\t',header=TRUE)
                gc_file <- "/share_data/wujm/project/getFastaGC/WES.T192V1.gc2.txt"
                ##连接数据库，根据检测样本，获取对应的cnn文件
                testSample <- input$sampleName
                conn <- dbConnect(MySQL(), host=dbconfig$host, dbname="nyuen_clinepilepsy", user=dbconfig$username, password=dbconfig$password)
                cmd <- paste0("select cnn from CNVBaseLine_T192V_WES where sample_name = '",testSample,"';")
                cnnFile = dbGetQuery(conn, cmd)
                newData <- data.table::fread(cnnFile$cnn,sep='\t',header=TRUE)
                names(newData) <- c("chrom","start","end","gene","depth","log2")
                # newData$percent <- newData$depth / sum(newData$depth)
                GC <- data.table::fread(gc_file,sep="\t",header=TRUE)
                mergeData <- merge(newData, GC, by=c("chrom","start","end"))
                filter0 <-  mergeData[which(depth != 0),]
                ##计算所有bin的平均覆盖度（Median）
                allMedianDepth <- median(filter0$depth)

                #计算相同GC含量的bin的覆盖度，计算这些覆盖度的median，每个GC含量梯度对应一个median值
                GC_group = group_by(filter0, gc)
                delay <- summarise(GC_group, median_sales = median(depth))
                corData = merge(filter0,delay,by='gc')
                corData$corDepth = corData$depth * (allMedianDepth / corData$median_sales)
                corData$corPercent = corData$corDepth / sum(corData$corDepth)
                corData <- corData[,c("chrom","start","end","gene","corDepth","corPercent")]


                #再将基线样本和检测样本合并
                percent_dat <- merge(corData, dat, by=c("chrom","start","end","gene"))
                percent_dat$t <- (percent_dat$corPercent - percent_dat$mean) / percent_dat$sd
                #chromStr <- input$chrom
                dataChr1 <- percent_dat[percent_dat$chrom == 'chr1',]

                output$chr1 <- renderPlotly({ggplotly(ggplot(data=dataChr1, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr1$gene))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,249250621) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr1_chrom <- renderPlot(chromNumber("1"))
                dataChr2 <- percent_dat[percent_dat$chrom == 'chr2',]
                output$chr2 <- renderPlotly({ggplotly(ggplot(data=dataChr2, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr2$gene))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,243199373) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr2_chrom <- renderPlot(chromNumber("2"))
                dataChr3 <- percent_dat[percent_dat$chrom == 'chr3',]
                output$chr3 <- renderPlotly({ggplotly(ggplot(data=dataChr3, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr3$gene))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,198022430) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr3_chrom <- renderPlot(chromNumber("3"))
                dataChr4 <- percent_dat[percent_dat$chrom == 'chr4',]
                output$chr4 <- renderPlotly({ggplotly(ggplot(data=dataChr4, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr4$gene))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,191154276) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr4_chrom <- renderPlot(chromNumber("4"))
                dataChr5 <- percent_dat[percent_dat$chrom == 'chr5',]
                output$chr5 <- renderPlotly({ggplotly(ggplot(data=dataChr5, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr5$gene))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,180915260) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr5_chrom <- renderPlot(chromNumber("5"))
                dataChr6 <- percent_dat[percent_dat$chrom == 'chr6',]
                output$chr6 <- renderPlotly({ggplotly(ggplot(data=dataChr6, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr6$gene))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,171115067) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr6_chrom <- renderPlot(chromNumber("6"))
                dataChr7 <- percent_dat[percent_dat$chrom == 'chr7',]
                output$chr7 <- renderPlotly({ggplotly(ggplot(data=dataChr7, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr7$gene))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,159138663) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr7_chrom <- renderPlot(chromNumber("7"))
                dataChr8 <- percent_dat[percent_dat$chrom == 'chr8',]
                output$chr8 <- renderPlotly({ggplotly(ggplot(data=dataChr8, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr8$gene))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,146364022) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr8_chrom <- renderPlot(chromNumber("8"))
                dataChr9 <- percent_dat[percent_dat$chrom == 'chr9',]
                output$chr9 <- renderPlotly({ggplotly(ggplot(data=dataChr9, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr9$gene))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,141213431) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr9_chrom <- renderPlot(chromNumber("9"))
                dataChr10 <- percent_dat[percent_dat$chrom == 'chr10',]
                output$chr10 <- renderPlotly({ggplotly(ggplot(data=dataChr10, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr10$gene))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,135534747) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr10_chrom <- renderPlot(chromNumber("10"))
                dataChr11 <- percent_dat[percent_dat$chrom == 'chr11',]
                output$chr11 <- renderPlotly({ggplotly(ggplot(data=dataChr11, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr11$gene))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,135006516) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr11_chrom <- renderPlot(chromNumber("11"))
                dataChr12 <- percent_dat[percent_dat$chrom == 'chr12',]
                output$chr12 <- renderPlotly({ggplotly(ggplot(data=dataChr12, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr12$gene))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,133851895) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr12_chrom <- renderPlot(chromNumber("12"))
                dataChr13 <- percent_dat[percent_dat$chrom == 'chr13',]
                output$chr13 <- renderPlotly({ggplotly(ggplot(data=dataChr13, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr13$gene))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,115169878) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr13_chrom <- renderPlot(chromNumber("13"))
                dataChr14 <- percent_dat[percent_dat$chrom == 'chr14',]
                output$chr14 <- renderPlotly({ggplotly(ggplot(data=dataChr14, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr14$gene))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,107349540) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr14_chrom <- renderPlot(chromNumber("14"))
                dataChr15 <- percent_dat[percent_dat$chrom == 'chr15',]
                output$chr15 <- renderPlotly({ggplotly(ggplot(data=dataChr15, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr15$gene))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,102531392) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr15_chrom <- renderPlot(chromNumber("15"))
                dataChr16 <- percent_dat[percent_dat$chrom == 'chr16',]
                output$chr16 <- renderPlotly({ggplotly(ggplot(data=dataChr16, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr16$gene))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,90354753) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr16_chrom <- renderPlot(chromNumber("16"))
                dataChr17 <- percent_dat[percent_dat$chrom == 'chr17',]
                output$chr17 <- renderPlotly({ggplotly(ggplot(data=dataChr17, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr17$gene))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,81195210) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr17_chrom <- renderPlot(chromNumber("17"))
                dataChr18 <- percent_dat[percent_dat$chrom == 'chr18',]
                output$chr18 <- renderPlotly({ggplotly(ggplot(data=dataChr18, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr18$gene))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,78077248) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr18_chrom <- renderPlot(chromNumber("18"))
                dataChr19 <- percent_dat[percent_dat$chrom == 'chr19',]
                output$chr19 <- renderPlotly({ggplotly(ggplot(data=dataChr19, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr19$gene))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,59128983) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr19_chrom <- renderPlot(chromNumber("19"))
                dataChr20 <- percent_dat[percent_dat$chrom == 'chr20',]
                output$chr20 <- renderPlotly({ggplotly(ggplot(data=dataChr20, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr20$gene))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,63025520) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr20_chrom <- renderPlot(chromNumber("20"))
                dataChr21 <- percent_dat[percent_dat$chrom == 'chr21',]
                output$chr21 <- renderPlotly({ggplotly(ggplot(data=dataChr21, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr21$gene))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,48129895) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr21_chrom <- renderPlot(chromNumber("21"))
                dataChr22 <- percent_dat[percent_dat$chrom == 'chr22',]
                output$chr22 <- renderPlotly({ggplotly(ggplot(data=dataChr22, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr22$gene))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,51304566) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr22_chrom <- renderPlot(chromNumber("22"))
                dataChrX <- percent_dat[percent_dat$chrom == 'chrX',]
                output$chrX <- renderPlotly({ggplotly(ggplot(data=dataChrX, mapping=aes(x=start, y=t, text = paste0("gene:",dataChrX$gene))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,155270560) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chrX_chrom <- renderPlot(chromNumber("X"))
                dataChrY <- percent_dat[percent_dat$chrom == 'chrY',]
                output$chrY <- renderPlotly({ggplotly(ggplot(data=dataChrY, mapping=aes(x=start, y=t, text = paste0("gene:",dataChrY$gene))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,59373566) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chrY_chrom <- renderPlot(chromNumber("Y"))

                ##output table
                DPPercentData <- data.table::fread(paste0(sampleDir,"baseLineCorDepthPercent.txt"),sep='\t',header=TRUE)
                MergeData <- merge(DPPercentData, corData, by=c("chrom","start","end","gene"))
                output$table <- DT::renderDataTable({MergeData})
                dbDisconnect(conn)

                outFileName = paste0(testSample,".depth.csv")
                output$downloadData = downloadHandler(
                  filename = function() {
                    outFileName
                   },
                  content = function(file) {
                    write.table(MergeData, file,sep=",", row.names= FALSE)
                   }
                )
              }
              else if(input$sex == "女"){
                sampleDir <- "/share_data/wujm/project/Shiny/CNV/reference_T192V/woman/"
                ##
                dat <- data.table::fread(paste0(sampleDir,"baseLineCorPercent.txt"),sep='\t',header=TRUE)
                gc_file <- "/share_data/wujm/project/getFastaGC/WES.T192V1.gc2.txt"
                ##连接数据库，根据检测样本，获取对应的cnn文件
                testSample <- input$sampleName
                conn <- dbConnect(MySQL(), host=dbconfig$host, dbname="nyuen_clinepilepsy", user=dbconfig$username, password=dbconfig$password)
                cmd <- paste0("select cnn from CNVBaseLine_T192V_WES where sample_name = '",testSample,"';")
                cnnFile = dbGetQuery(conn, cmd)
                newData <- data.table::fread(cnnFile$cnn,sep='\t',header=TRUE)
                names(newData) <- c("chrom","start","end","gene","depth","log2")
                # newData$percent <- newData$depth / sum(newData$depth)
                GC <- data.table::fread(gc_file,sep="\t",header=TRUE)
                mergeData <- merge(newData, GC, by=c("chrom","start","end"))
                filter0 <-  mergeData[which(depth != 0),]
                ##计算所有bin的平均覆盖度（Median）
                allMedianDepth <- median(filter0$depth)

                #计算相同GC含量的bin的覆盖度，计算这些覆盖度的median，每个GC含量梯度对应一个median值
                GC_group = group_by(filter0, gc)
                delay <- summarise(GC_group, median_sales = median(depth))
                corData = merge(filter0,delay,by='gc')
                corData$corDepth = corData$depth * (allMedianDepth / corData$median_sales)
                corData$corPercent = corData$corDepth / sum(corData$corDepth)
                corData <- corData[,c("chrom","start","end","gene","corDepth","corPercent")]


                #再将基线样本和检测样本合并
                percent_dat <- merge(corData, dat, by=c("chrom","start","end","gene"))
                percent_dat$t <- (percent_dat$corPercent - percent_dat$mean) / percent_dat$sd
                #chromStr <- input$chrom
                dataChr1 <- percent_dat[percent_dat$chrom == 'chr1',]

                output$chr1 <- renderPlotly({ggplotly(ggplot(data=dataChr1, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr1$gene))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,249250621) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr1_chrom <- renderPlot(chromNumber("1"))
                dataChr2 <- percent_dat[percent_dat$chrom == 'chr2',]
                output$chr2 <- renderPlotly({ggplotly(ggplot(data=dataChr2, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr2$gene))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,243199373) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr2_chrom <- renderPlot(chromNumber("2"))
                dataChr3 <- percent_dat[percent_dat$chrom == 'chr3',]
                output$chr3 <- renderPlotly({ggplotly(ggplot(data=dataChr3, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr3$gene))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,198022430) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr3_chrom <- renderPlot(chromNumber("3"))
                dataChr4 <- percent_dat[percent_dat$chrom == 'chr4',]
                output$chr4 <- renderPlotly({ggplotly(ggplot(data=dataChr4, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr4$gene))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,191154276) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr4_chrom <- renderPlot(chromNumber("4"))
                dataChr5 <- percent_dat[percent_dat$chrom == 'chr5',]
                output$chr5 <- renderPlotly({ggplotly(ggplot(data=dataChr5, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr5$gene))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,180915260) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr5_chrom <- renderPlot(chromNumber("5"))
                dataChr6 <- percent_dat[percent_dat$chrom == 'chr6',]
                output$chr6 <- renderPlotly({ggplotly(ggplot(data=dataChr6, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr6$gene))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,171115067) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr6_chrom <- renderPlot(chromNumber("6"))
                dataChr7 <- percent_dat[percent_dat$chrom == 'chr7',]
                output$chr7 <- renderPlotly({ggplotly(ggplot(data=dataChr7, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr7$gene))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,159138663) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr7_chrom <- renderPlot(chromNumber("7"))
                dataChr8 <- percent_dat[percent_dat$chrom == 'chr8',]
                output$chr8 <- renderPlotly({ggplotly(ggplot(data=dataChr8, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr8$gene))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,146364022) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr8_chrom <- renderPlot(chromNumber("8"))
                dataChr9 <- percent_dat[percent_dat$chrom == 'chr9',]
                output$chr9 <- renderPlotly({ggplotly(ggplot(data=dataChr9, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr9$gene))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,141213431) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr9_chrom <- renderPlot(chromNumber("9"))
                dataChr10 <- percent_dat[percent_dat$chrom == 'chr10',]
                output$chr10 <- renderPlotly({ggplotly(ggplot(data=dataChr10, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr10$gene))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,135534747) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr10_chrom <- renderPlot(chromNumber("10"))
                dataChr11 <- percent_dat[percent_dat$chrom == 'chr11',]
                output$chr11 <- renderPlotly({ggplotly(ggplot(data=dataChr11, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr11$gene))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,135006516) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr11_chrom <- renderPlot(chromNumber("11"))
                dataChr12 <- percent_dat[percent_dat$chrom == 'chr12',]
                output$chr12 <- renderPlotly({ggplotly(ggplot(data=dataChr12, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr12$gene))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,133851895) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr12_chrom <- renderPlot(chromNumber("12"))
                dataChr13 <- percent_dat[percent_dat$chrom == 'chr13',]
                output$chr13 <- renderPlotly({ggplotly(ggplot(data=dataChr13, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr13$gene))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,115169878) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr13_chrom <- renderPlot(chromNumber("13"))
                dataChr14 <- percent_dat[percent_dat$chrom == 'chr14',]
                output$chr14 <- renderPlotly({ggplotly(ggplot(data=dataChr14, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr14$gene))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,107349540) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr14_chrom <- renderPlot(chromNumber("14"))
                dataChr15 <- percent_dat[percent_dat$chrom == 'chr15',]
                output$chr15 <- renderPlotly({ggplotly(ggplot(data=dataChr15, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr15$gene))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,102531392) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr15_chrom <- renderPlot(chromNumber("15"))
                dataChr16 <- percent_dat[percent_dat$chrom == 'chr16',]
                output$chr16 <- renderPlotly({ggplotly(ggplot(data=dataChr16, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr16$gene))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,90354753) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr16_chrom <- renderPlot(chromNumber("16"))
                dataChr17 <- percent_dat[percent_dat$chrom == 'chr17',]
                output$chr17 <- renderPlotly({ggplotly(ggplot(data=dataChr17, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr17$gene))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,81195210) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr17_chrom <- renderPlot(chromNumber("17"))
                dataChr18 <- percent_dat[percent_dat$chrom == 'chr18',]
                output$chr18 <- renderPlotly({ggplotly(ggplot(data=dataChr18, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr18$gene))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,78077248) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr18_chrom <- renderPlot(chromNumber("18"))
                dataChr19 <- percent_dat[percent_dat$chrom == 'chr19',]
                output$chr19 <- renderPlotly({ggplotly(ggplot(data=dataChr19, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr19$gene))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,59128983) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr19_chrom <- renderPlot(chromNumber("19"))
                dataChr20 <- percent_dat[percent_dat$chrom == 'chr20',]
                output$chr20 <- renderPlotly({ggplotly(ggplot(data=dataChr20, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr20$gene))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,63025520) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr20_chrom <- renderPlot(chromNumber("20"))
                dataChr21 <- percent_dat[percent_dat$chrom == 'chr21',]
                output$chr21 <- renderPlotly({ggplotly(ggplot(data=dataChr21, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr21$gene))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,48129895) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr21_chrom <- renderPlot(chromNumber("21"))
                dataChr22 <- percent_dat[percent_dat$chrom == 'chr22',]
                output$chr22 <- renderPlotly({ggplotly(ggplot(data=dataChr22, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr22$gene))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,51304566) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr22_chrom <- renderPlot(chromNumber("22"))
                dataChrX <- percent_dat[percent_dat$chrom == 'chrX',]
                output$chrX <- renderPlotly({ggplotly(ggplot(data=dataChrX, mapping=aes(x=start, y=t, text = paste0("gene:",dataChrX$gene))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,155270560) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chrX_chrom <- renderPlot(chromNumber("X"))
                dataChrY <- percent_dat[percent_dat$chrom == 'chrY',]
                output$chrY <- renderPlotly({ggplotly(ggplot(data=dataChrY, mapping=aes(x=start, y=t, text = paste0("gene:",dataChrY$gene))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,59373566) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chrY_chrom <- renderPlot(chromNumber("Y"))

                ##output table
                DPPercentData <- data.table::fread(paste0(sampleDir,"baseLineCorDepthPercent.txt"),sep='\t',header=TRUE)
                MergeData <- merge(DPPercentData, corData, by=c("chrom","start","end","gene"))
                output$table <- DT::renderDataTable({MergeData})
                dbDisconnect(conn)

                outFileName = paste0(testSample,".depth.csv")
                output$downloadData = downloadHandler(
                  filename = function() {
                    outFileName
                   },
                  content = function(file) {
                    write.table(MergeData, file,sep=",", row.names= FALSE)
                   }
                )
              }
            }
            else if(input$product_type == "Panel"){
              if(input$sex == "男"){
                sampleDir <- "/share_data/wujm/project/Shiny/CNV/reference/Panel/man/"
                ##
                dat <- data.table::fread(paste0(sampleDir,"baseLineCorPercent.txt"),sep='\t',header=TRUE)
                gc_file <- "/share_data/wujm/project/getFastaGC/Panel.extra150.gc.txt"
                ##连接数据库，根据检测样本，获取对应的cnn文件
                testSample <- input$sampleName
                conn <- dbConnect(MySQL(), host=dbconfig$host, dbname="nyuen_clinepilepsy", user=dbconfig$username, password=dbconfig$password)
                cmd <- paste0("select cnn from CNVBaseLine_Panel where sample_name = '",testSample,"';")
                cnnFile = dbGetQuery(conn, cmd)
                newData <- data.table::fread(cnnFile$cnn,sep='\t',header=TRUE)
                names(newData) <- c("chrom","start","end","gene","depth","log2")
                # newData$percent <- newData$depth / sum(newData$depth)
                GC <- data.table::fread(gc_file,sep="\t",header=TRUE)
                mergeData <- merge(newData, GC, by=c("chrom","start","end"))
                filter0 <-  mergeData[which(depth != 0),]
                ##计算所有bin的平均覆盖度（Median）
                allMedianDepth <- median(filter0$depth)

                #计算相同GC含量的bin的覆盖度，计算这些覆盖度的median，每个GC含量梯度对应一个median值
                GC_group = group_by(filter0, gc)
                delay <- summarise(GC_group, median_sales = median(depth))
                corData = merge(filter0,delay,by='gc')
                corData$corDepth = corData$depth * (allMedianDepth / corData$median_sales)
                corData$corPercent = corData$corDepth / sum(corData$corDepth)
                corData <- corData[,c("chrom","start","end","gene","corDepth","corPercent")]
                

                #再将基线样本和检测样本合并
                percent_dat <- merge(corData, dat, by=c("chrom","start","end","gene"))
                percent_dat$t <- (percent_dat$corPercent - percent_dat$mean) / percent_dat$sd
                #chromStr <- input$chrom
                dataChr1 <- percent_dat[percent_dat$chrom == 'chr1',]
                output$chr1 <- renderPlotly({ggplotly(ggplot(data=dataChr1, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr1$gene,"\n","area:",dataChr1$start,"-",dataChr1$end))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,249250621) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr1_chrom <- renderPlot(chromNumber("1"))
                dataChr2 <- percent_dat[percent_dat$chrom == 'chr2',]
                output$chr2 <- renderPlotly({ggplotly(ggplot(data=dataChr2, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr2$gene,"\n","area:",dataChr2$start,"-",dataChr2$end))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,243199373) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr2_chrom <- renderPlot(chromNumber("2"))
                dataChr3 <- percent_dat[percent_dat$chrom == 'chr3',]
                output$chr3 <- renderPlotly({ggplotly(ggplot(data=dataChr3, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr3$gene,"\n","area:",dataChr3$start,"-",dataChr3$end))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,198022430) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr3_chrom <- renderPlot(chromNumber("3"))
                dataChr4 <- percent_dat[percent_dat$chrom == 'chr4',]
                output$chr4 <- renderPlotly({ggplotly(ggplot(data=dataChr4, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr4$gene,"\n","area:",dataChr4$start,"-",dataChr4$end))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,191154276) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr4_chrom <- renderPlot(chromNumber("4"))
                dataChr5 <- percent_dat[percent_dat$chrom == 'chr5',]
                output$chr5 <- renderPlotly({ggplotly(ggplot(data=dataChr5, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr5$gene,"\n","area:",dataChr5$start,"-",dataChr5$end))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,180915260) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr5_chrom <- renderPlot(chromNumber("5"))
                dataChr6 <- percent_dat[percent_dat$chrom == 'chr6',]
                output$chr6 <- renderPlotly({ggplotly(ggplot(data=dataChr6, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr6$gene,"\n","area:",dataChr6$start,"-",dataChr6$end))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,171115067) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr6_chrom <- renderPlot(chromNumber("6"))
                dataChr7 <- percent_dat[percent_dat$chrom == 'chr7',]
                output$chr7 <- renderPlotly({ggplotly(ggplot(data=dataChr7, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr7$gene,"\n","area:",dataChr7$start,"-",dataChr7$end))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,159138663) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr7_chrom <- renderPlot(chromNumber("7"))
                dataChr8 <- percent_dat[percent_dat$chrom == 'chr8',]
                output$chr8 <- renderPlotly({ggplotly(ggplot(data=dataChr8, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr8$gene,"\n","area:",dataChr8$start,"-",dataChr8$end))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,146364022) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr8_chrom <- renderPlot(chromNumber("8"))
                dataChr9 <- percent_dat[percent_dat$chrom == 'chr9',]
                output$chr9 <- renderPlotly({ggplotly(ggplot(data=dataChr9, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr9$gene,"\n","area:",dataChr9$start,"-",dataChr9$end))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,141213431) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr9_chrom <- renderPlot(chromNumber("9"))
                dataChr10 <- percent_dat[percent_dat$chrom == 'chr10',]
                output$chr10 <- renderPlotly({ggplotly(ggplot(data=dataChr10, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr10$gene,"\n","area:",dataChr10$start,"-",dataChr10$end))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,135534747) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr10_chrom <- renderPlot(chromNumber("10"))
                dataChr11 <- percent_dat[percent_dat$chrom == 'chr11',]
                output$chr11 <- renderPlotly({ggplotly(ggplot(data=dataChr11, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr11$gene,"\n","area:",dataChr11$start,"-",dataChr11$end))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,135006516) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr11_chrom <- renderPlot(chromNumber("11"))
                dataChr12 <- percent_dat[percent_dat$chrom == 'chr12',]
                output$chr12 <- renderPlotly({ggplotly(ggplot(data=dataChr12, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr12$gene,"\n","area:",dataChr12$start,"-",dataChr12$end))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,133851895) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr12_chrom <- renderPlot(chromNumber("12"))
                dataChr13 <- percent_dat[percent_dat$chrom == 'chr13',]
                output$chr13 <- renderPlotly({ggplotly(ggplot(data=dataChr13, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr13$gene,"\n","area:",dataChr13$start,"-",dataChr13$end))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,115169878) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr13_chrom <- renderPlot(chromNumber("13"))
                dataChr14 <- percent_dat[percent_dat$chrom == 'chr14',]
                output$chr14 <- renderPlotly({ggplotly(ggplot(data=dataChr14, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr14$gene,"\n","area:",dataChr14$start,"-",dataChr14$end))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,107349540) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr14_chrom <- renderPlot(chromNumber("14"))
                dataChr15 <- percent_dat[percent_dat$chrom == 'chr15',]
                output$chr15 <- renderPlotly({ggplotly(ggplot(data=dataChr15, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr15$gene,"\n","area:",dataChr15$start,"-",dataChr15$end))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,102531392) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr15_chrom <- renderPlot(chromNumber("15"))
                dataChr16 <- percent_dat[percent_dat$chrom == 'chr16',]
                output$chr16 <- renderPlotly({ggplotly(ggplot(data=dataChr16, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr16$gene,"\n","area:",dataChr16$start,"-",dataChr16$end))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,90354753) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr16_chrom <- renderPlot(chromNumber("16"))
                dataChr17 <- percent_dat[percent_dat$chrom == 'chr17',]
                output$chr17 <- renderPlotly({ggplotly(ggplot(data=dataChr17, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr17$gene,"\n","area:",dataChr17$start,"-",dataChr17$end))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,81195210) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr17_chrom <- renderPlot(chromNumber("17"))
                dataChr18 <- percent_dat[percent_dat$chrom == 'chr18',]
                output$chr18 <- renderPlotly({ggplotly(ggplot(data=dataChr18, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr18$gene,"\n","area:",dataChr18$start,"-",dataChr18$end))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,78077248) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr18_chrom <- renderPlot(chromNumber("18"))
                dataChr19 <- percent_dat[percent_dat$chrom == 'chr19',]
                output$chr19 <- renderPlotly({ggplotly(ggplot(data=dataChr19, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr19$gene,"\n","area:",dataChr19$start,"-",dataChr19$end))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,59128983) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr19_chrom <- renderPlot(chromNumber("19"))
                dataChr20 <- percent_dat[percent_dat$chrom == 'chr20',]
                output$chr20 <- renderPlotly({ggplotly(ggplot(data=dataChr20, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr20$gene,"\n","area:",dataChr20$start,"-",dataChr20$end))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,63025520) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr20_chrom <- renderPlot(chromNumber("20"))
                dataChr21 <- percent_dat[percent_dat$chrom == 'chr21',]
                output$chr21 <- renderPlotly({ggplotly(ggplot(data=dataChr21, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr21$gene,"\n","area:",dataChr21$start,"-",dataChr21$end))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,48129895) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr21_chrom <- renderPlot(chromNumber("21"))
                dataChr22 <- percent_dat[percent_dat$chrom == 'chr22',]
                output$chr22 <- renderPlotly({ggplotly(ggplot(data=dataChr22, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr22$gene,"\n","area:",dataChr22$start,"-",dataChr22$end))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,51304566) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr22_chrom <- renderPlot(chromNumber("22"))
                dataChrX <- percent_dat[percent_dat$chrom == 'chrX',]
                output$chrX <- renderPlotly({ggplotly(ggplot(data=dataChrX, mapping=aes(x=start, y=t, text = paste0("gene:",dataChrX$gene,"\n","area:",dataChrX$start,"-",dataChrX$end))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,155270560) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chrX_chrom <- renderPlot(chromNumber("X"))
                dataChrY <- percent_dat[percent_dat$chrom == 'chrY',]
                output$chrY <- renderPlotly({ggplotly(ggplot(data=dataChrY, mapping=aes(x=start, y=t, text = paste0("gene:",dataChrY$gene,"\n","area:",dataChrY$start,"-",dataChrY$end))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,59373566) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chrY_chrom <- renderPlot(chromNumber("Y"))

                ##output table
                DPPercentData <- data.table::fread(paste0(sampleDir,"baseLineCorDepthPercent.txt"),sep='\t',header=TRUE)
                MergeData <- merge(DPPercentData, corData, by=c("chrom","start","end","gene"))
                output$table <- DT::renderDataTable({MergeData})
                dbDisconnect(conn)

                outFileName = paste0(testSample,"depth.csv")
                output$downloadData = downloadHandler(
                  filename = function() {
                    outFileName
                   },
                  content = function(file) {
                    write.table(MergeData, file,sep=",", row.names= FALSE)
                   }
                )
              }else if(input$sex == "女"){
                sampleDir <- "/share_data/wujm/project/Shiny/CNV/reference/Panel/woman/"
                ##
                dat <- data.table::fread(paste0(sampleDir,"baseLineCorPercent.txt"),sep='\t',header=TRUE)
                gc_file <- "/share_data/wujm/project/getFastaGC/Panel.extra150.gc.txt"
                ##连接数据库，根据检测样本，获取对应的cnn文件
                testSample <- input$sampleName
                conn <- dbConnect(MySQL(), host=dbconfig$host, dbname="nyuen_clinepilepsy", user=dbconfig$username, password=dbconfig$password)
                cmd <- paste0("select cnn from CNVBaseLine_Panel where sample_name = '",testSample,"';")
                cnnFile = dbGetQuery(conn, cmd)
                newData <- data.table::fread(cnnFile$cnn,sep='\t',header=TRUE)
                names(newData) <- c("chrom","start","end","gene","depth","log2")
                # newData$percent <- newData$depth / sum(newData$depth)
                GC <- data.table::fread(gc_file,sep="\t",header=TRUE)
                mergeData <- merge(newData, GC, by=c("chrom","start","end"))
                filter0 <-  mergeData[which(depth != 0),]
                ##计算所有bin的平均覆盖度（Median）
                allMedianDepth <- median(filter0$depth)

                #计算相同GC含量的bin的覆盖度，计算这些覆盖度的median，每个GC含量梯度对应一个median值
                GC_group = group_by(filter0, gc)
                delay <- summarise(GC_group, median_sales = median(depth))
                corData = merge(filter0,delay,by='gc')
                corData$corDepth = corData$depth * (allMedianDepth / corData$median_sales)
                corData$corPercent = corData$corDepth / sum(corData$corDepth)
                corData <- corData[,c("chrom","start","end","gene","corDepth","corPercent")]
                

                #再将基线样本和检测样本合并
                percent_dat <- merge(corData, dat, by=c("chrom","start","end","gene"))
                percent_dat$t <- (percent_dat$corPercent - percent_dat$mean) / percent_dat$sd
                dataChr1 <- percent_dat[percent_dat$chrom == 'chr1',]
                output$chr1 <- renderPlotly({ggplotly(ggplot(data=dataChr1, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr1$gene,"\n","area:",dataChr1$start,"-",dataChr1$end))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,249250621) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr1_chrom <- renderPlot(chromNumber("1"))
                dataChr2 <- percent_dat[percent_dat$chrom == 'chr2',]
                output$chr2 <- renderPlotly({ggplotly(ggplot(data=dataChr2, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr2$gene,"\n","area:",dataChr2$start,"-",dataChr2$end))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,243199373) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr2_chrom <- renderPlot(chromNumber("2"))
                dataChr3 <- percent_dat[percent_dat$chrom == 'chr3',]
                output$chr3 <- renderPlotly({ggplotly(ggplot(data=dataChr3, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr3$gene,"\n","area:",dataChr3$start,"-",dataChr3$end))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,198022430) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr3_chrom <- renderPlot(chromNumber("3"))
                dataChr4 <- percent_dat[percent_dat$chrom == 'chr4',]
                output$chr4 <- renderPlotly({ggplotly(ggplot(data=dataChr4, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr4$gene,"\n","area:",dataChr4$start,"-",dataChr4$end))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,191154276) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr4_chrom <- renderPlot(chromNumber("4"))
                dataChr5 <- percent_dat[percent_dat$chrom == 'chr5',]
                output$chr5 <- renderPlotly({ggplotly(ggplot(data=dataChr5, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr5$gene,"\n","area:",dataChr5$start,"-",dataChr5$end))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,180915260) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr5_chrom <- renderPlot(chromNumber("5"))
                dataChr6 <- percent_dat[percent_dat$chrom == 'chr6',]
                output$chr6 <- renderPlotly({ggplotly(ggplot(data=dataChr6, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr6$gene,"\n","area:",dataChr6$start,"-",dataChr6$end))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,171115067) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr6_chrom <- renderPlot(chromNumber("6"))
                dataChr7 <- percent_dat[percent_dat$chrom == 'chr7',]
                output$chr7 <- renderPlotly({ggplotly(ggplot(data=dataChr7, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr7$gene,"\n","area:",dataChr7$start,"-",dataChr7$end))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,159138663) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr7_chrom <- renderPlot(chromNumber("7"))
                dataChr8 <- percent_dat[percent_dat$chrom == 'chr8',]
                output$chr8 <- renderPlotly({ggplotly(ggplot(data=dataChr8, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr8$gene,"\n","area:",dataChr8$start,"-",dataChr8$end))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,146364022) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr8_chrom <- renderPlot(chromNumber("8"))
                dataChr9 <- percent_dat[percent_dat$chrom == 'chr9',]
                output$chr9 <- renderPlotly({ggplotly(ggplot(data=dataChr9, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr9$gene,"\n","area:",dataChr9$start,"-",dataChr9$end))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,141213431) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr9_chrom <- renderPlot(chromNumber("9"))
                dataChr10 <- percent_dat[percent_dat$chrom == 'chr10',]
                output$chr10 <- renderPlotly({ggplotly(ggplot(data=dataChr10, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr10$gene,"\n","area:",dataChr10$start,"-",dataChr10$end))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,135534747) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr10_chrom <- renderPlot(chromNumber("10"))
                dataChr11 <- percent_dat[percent_dat$chrom == 'chr11',]
                output$chr11 <- renderPlotly({ggplotly(ggplot(data=dataChr11, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr11$gene,"\n","area:",dataChr11$start,"-",dataChr11$end))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,135006516) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr11_chrom <- renderPlot(chromNumber("11"))
                dataChr12 <- percent_dat[percent_dat$chrom == 'chr12',]
                output$chr12 <- renderPlotly({ggplotly(ggplot(data=dataChr12, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr12$gene,"\n","area:",dataChr12$start,"-",dataChr12$end))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,133851895) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr12_chrom <- renderPlot(chromNumber("12"))
                dataChr13 <- percent_dat[percent_dat$chrom == 'chr13',]
                output$chr13 <- renderPlotly({ggplotly(ggplot(data=dataChr13, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr13$gene,"\n","area:",dataChr13$start,"-",dataChr13$end))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,115169878) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr13_chrom <- renderPlot(chromNumber("13"))
                dataChr14 <- percent_dat[percent_dat$chrom == 'chr14',]
                output$chr14 <- renderPlotly({ggplotly(ggplot(data=dataChr14, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr14$gene,"\n","area:",dataChr14$start,"-",dataChr14$end))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,107349540) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr14_chrom <- renderPlot(chromNumber("14"))
                dataChr15 <- percent_dat[percent_dat$chrom == 'chr15',]
                output$chr15 <- renderPlotly({ggplotly(ggplot(data=dataChr15, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr15$gene,"\n","area:",dataChr15$start,"-",dataChr15$end))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,102531392) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr15_chrom <- renderPlot(chromNumber("15"))
                dataChr16 <- percent_dat[percent_dat$chrom == 'chr16',]
                output$chr16 <- renderPlotly({ggplotly(ggplot(data=dataChr16, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr16$gene,"\n","area:",dataChr16$start,"-",dataChr16$end))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,90354753) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr16_chrom <- renderPlot(chromNumber("16"))
                dataChr17 <- percent_dat[percent_dat$chrom == 'chr17',]
                output$chr17 <- renderPlotly({ggplotly(ggplot(data=dataChr17, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr17$gene,"\n","area:",dataChr17$start,"-",dataChr17$end))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,81195210) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr17_chrom <- renderPlot(chromNumber("17"))
                dataChr18 <- percent_dat[percent_dat$chrom == 'chr18',]
                output$chr18 <- renderPlotly({ggplotly(ggplot(data=dataChr18, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr18$gene,"\n","area:",dataChr18$start,"-",dataChr18$end))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,78077248) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr18_chrom <- renderPlot(chromNumber("18"))
                dataChr19 <- percent_dat[percent_dat$chrom == 'chr19',]
                output$chr19 <- renderPlotly({ggplotly(ggplot(data=dataChr19, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr19$gene,"\n","area:",dataChr19$start,"-",dataChr19$end))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,59128983) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr19_chrom <- renderPlot(chromNumber("19"))
                dataChr20 <- percent_dat[percent_dat$chrom == 'chr20',]
                output$chr20 <- renderPlotly({ggplotly(ggplot(data=dataChr20, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr20$gene,"\n","area:",dataChr20$start,"-",dataChr20$end))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,63025520) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr20_chrom <- renderPlot(chromNumber("20"))
                dataChr21 <- percent_dat[percent_dat$chrom == 'chr21',]
                output$chr21 <- renderPlotly({ggplotly(ggplot(data=dataChr21, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr21$gene,"\n","area:",dataChr21$start,"-",dataChr21$end))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,48129895) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr21_chrom <- renderPlot(chromNumber("21"))
                dataChr22 <- percent_dat[percent_dat$chrom == 'chr22',]
                output$chr22 <- renderPlotly({ggplotly(ggplot(data=dataChr22, mapping=aes(x=start, y=t, text = paste0("gene:",dataChr22$gene,"\n","area:",dataChr22$start,"-",dataChr22$end))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,51304566) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chr22_chrom <- renderPlot(chromNumber("22"))
                dataChrX <- percent_dat[percent_dat$chrom == 'chrX',]
                output$chrX <- renderPlotly({ggplotly(ggplot(data=dataChrX, mapping=aes(x=start, y=t, text = paste0("gene:",dataChrX$gene,"\n","area:",dataChrX$start,"-",dataChrX$end))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,155270560) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chrX_chrom <- renderPlot(chromNumber("X"))
                dataChrY <- percent_dat[percent_dat$chrom == 'chrY',]
                output$chrY <- renderPlotly({ggplotly(ggplot(data=dataChrY, mapping=aes(x=start, y=t, text = paste0("gene:",dataChrY$gene,"\n","area:",dataChrY$start,"-",dataChrY$end))) +geom_point( color="red", size=0.1) +ylim(-15,15) + xlim(0,59373566) + labs(x='position',y='t score')+theme(legend.position="blank",plot.title = element_text(hjust = 0.5)))})
                output$chrY_chrom <- renderPlot(chromNumber("Y"))
                output$chrY_chrom <- renderPlot(chromNumber("Y"))

                ##output table
                DPPercentData <- data.table::fread(paste0(sampleDir,"baseLineCorDepthPercent.txt"),sep='\t',header=TRUE)
                MergeData <- merge(DPPercentData, corData, by=c("chrom","start","end","gene"))
                output$table <- DT::renderDataTable({MergeData})
                dbDisconnect(conn)

                outFileName = paste0(testSample,".depth.csv")
                output$downloadData = downloadHandler(
                  filename = function() {
                    outFileName
                   },
                  content = function(file) {
                    write.table(MergeData, file,sep=",", row.names= FALSE)
                   }
                )
              }
            }

    })
}
