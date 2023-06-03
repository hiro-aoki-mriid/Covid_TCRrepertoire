##Merge different abundance clone into tcr repertoire dataset

library(foreach)
library(doParallel)
library(stringr)
library(data.table)
library(dplyr)

dir.input <- "tmp/result/intermediate/2_AIM/JoinTP_DifAbund_AIM"

cores <- 24

tableread_fast = function(i, header=TRUE, sep="\t"){
  tmp = fread(i, header=header, sep=sep, nThread=32)
  tmp = as.data.frame(tmp)
  return(tmp)
}
createEmptyDf = function( nrow, ncol, colnames = c() ){
  data.frame( matrix( vector(), nrow, ncol, dimnames = list( c(), colnames ) ) )
}

############################### Functions ########################################
##Extract beta-binomial test result
Extract <- function(file.name, dir.input, FirstSecondType, ThirdType, condition){
  #load repetoire data
  name.input <- str_c(dir.input, file.name, sep = "/")
  d <- tableread_fast(name.input, header = TRUE, sep=",")
  TPs <- c("TP1", "TP2", "TP3", "TP4", "TP5", "TP6", "TP8", "TP9", "TP10", "TP11")
  
  ###Summarize
  ##Count clones
  #Create output data.frame
  output <- createEmptyDf(1, length(TPs), colnames = TPs)
  dsub <- dplyr::filter(d, class == "AIM")
  if(nrow(dsub) > 0){
    output[1, 1:length(TPs)] <- dplyr::summarize(dsub, dplyr::across(TPs, sum) )
  } else{
    output[1, 1:length(TPs)] <- rep(0, times = length(TPs))
  }

  output$sample <- file.name
  
  return(output)
}

############################### Processing ########################################
files  <- list.files(dir.input, pattern="NaraCOVID")
  
t<-proc.time()
cl <- makeCluster(cores)
registerDoParallel(cl)   
output <- foreach(file.name = files, .combine = rbind,
                 .packages=c("data.table", "stringr", "dplyr")) %dopar% {Extract(file.name, dir.input,
                                                                                 FirstSecondType, ThirdType, condition)}
stopCluster(cl)
proc.time()-t
  
name.output <- "result/Fig2/AIM.timepoint.csv"
write.csv(output, name.output, row.names = FALSE)
