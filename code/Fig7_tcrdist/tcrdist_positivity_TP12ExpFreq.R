###Calculate frequency of SARS-Cov-2 reactive clones determined by tcrdist3 at each time point

#load libraryies
library(ggplot2)
library(extrafont)
loadfonts("win")
library(foreach)
library(doParallel)
library(stringr)
library(dplyr)
library(data.table)

#Input layer
name.input <- "tmp/result/Fig7/tcrdist_clones.csv"
name.hla <- "tmp/metadata/HLAtyping_result.csv"
dir.output <- "tmp/result/Fig7"
name.output <- "tcrdist_epitope_TP12Exp_Freq.csv"
cores <- 12
FirstSecondType <- c("Early", "Main", "Third")
epitope_query <-c("S269", "S448", "S919", "S1208")
TPs <- c("TP1", "TP2", "TP3", "TP4", "TP5", "TP6", "TP8", "TP9", "TP10", "TP11")

############################### Define functions ########################################
###Define functions
#Table.read.fast
tableread_fast = function(i, header=TRUE, quote="", sep=","){
  tmp = fread(i, header=header, sep=sep, quote=quote, nThread=32)
  tmp = as.data.frame(tmp)
  return(tmp)
}
createEmptyDf = function( nrow, ncol, colnames = c() ){
  data.frame( matrix( vector(), nrow, ncol, dimnames = list( c(), colnames ) ) )
}

############################### Processing ########################################
##load data
#load tcrdist3 data
dir.create(dir.output, recursive = TRUE)
d <- tableread_fast(name.input, header = TRUE, quote="\"", sep=",")
d <- dplyr::filter(d, motif == "Posi")

output.all <- data.frame()
for(epitope in epitope_query){
  #Extract table of each epitope
  d$motif <- d[,which(names(d) == str_c(epitope, "_motif"))]
  d_sub <- dplyr::filter(d, motif == "Posi")
  IDs <- unique(d_sub$sample)
  
  for(id in IDs){
    #Extract data of each participants
    d_sub2 <- dplyr::filter(d_sub, sample == id)
    #Prepare output data frame
    output <- createEmptyDf(length(FirstSecondType), length(TPs), colnames = TPs)
    #Check TP12 exp pattern in tcrdist3 hits
    for(j in 1:length(FirstSecondType)){
        fst <- FirstSecondType[j]
        d_sub3 <- dplyr::filter(d_sub2, Resp == fst)
        if(nrow(d_sub3) > 0){
          output[j, 1:length(TPs)] <- dplyr::summarize(d_sub3, dplyr::across(TPs, sum) )
        } else{
          output[j, 1:length(TPs)] <- rep(0, times = length(TPs))
        }
    }
    output$TP12Exp <- FirstSecondType
    output$ID <- id
    output$Epitope <- epitope
    output.all <- rbind(output.all, output)
  }
}

name.output <- str_c(dir.output, name.output, sep = "/")
write.csv(output.all, name.output, row.names = FALSE)
