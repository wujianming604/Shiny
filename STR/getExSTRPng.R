argv <- commandArgs(TRUE)
sampleName <- argv[1]
library(ggplot2)
library(plotly)
library(data.table)
library(exSTRa)

workDir <- "/share_data/wujm/project/run_exSTR/run_NUOHE/"

#all_locus <- c("DRPLA","HD","HDL2","SBMA","SCA1","SCA2","SCA3","SCA6","SCA7","SCA12","SCA17")
all_locus <- c("SCA8")
str_score <- read_score(
  file <- paste0(workDir,"samples_NUOHE.exSTR2.txt"), 
  #database <- "/share_data/wujm/software/exSTRa/inst/extdata/nuohe_repeat_expansion_disorders_hg19.txt",
  database <- "/share_data/wujm/software/exSTRa/inst/extdata/repeat_expansion_disorders_hg19.txt",
  groups.regex = c(control = "^WGSrpt_0[24]$", case = "")
)
colo <- c(sample = "red")
names(colo) <- c(sampleName)

for(locus in all_locus){
    outfile <- paste0(sampleName,".",locus,'.png')
    png(outfile, width = 1232,height = 665)
    plot(str_score, locus, sample_col = colo)
    dev.off()
    # Return a list containing the filename
    list(src = outfile,
        contentType = 'image/png',
        width = 1232,
        height = 665)
  }
