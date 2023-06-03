##Merge different abundance clone into tcr repertoire dataset

library(foreach)
library(doParallel)
library(stringr)
library(data.table)
library(dplyr)

dir.input <- "tmp/result/intermediate/2_AIM/JoinTP_DifAbund_AIM"

#Output params

FirstSecondType <- c("1st", "1stInc", "2nd", "2ndInc", "Dual", "DualInc", "3rd", "3rdInc", "Others")
AIMType <- c("AIM", "NonNV")

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
Extract <- function(file.name, dir.input, FirstSecondType, AIMType, condition){
  #load repetoire data
  name.input <- str_c(dir.input, file.name, sep = "/")
  d <- tableread_fast(name.input, header = TRUE, sep=",")
  
  ###Summarize
  ##Count clones
  #Create output data.frame
  output <- createEmptyDf(length(FirstSecondType), length(AIMType), colnames = AIMType)
  for(i in 1:length(FirstSecondType)){
    fst <- FirstSecondType[i]
    dsub <- dplyr::filter(d, Resp == fst)
    #Count of clones
    for(j in 1:length(AIMType)){
      tt <- AIMType[j]
      output[i,j] <- nrow(dplyr::filter(dsub, class == tt))
    }
  }
  if(sum(output$AIM) > 9){
    output$AIM_prop <- output$AIM / sum(output$AIM)
  } else {
    output$AIM_prop <- "Non"
  }
  if(sum(output$NonNV) > 9){
    output$NonNV_prop <- output$NonNV / sum(output$NonNV)
  } else {
    output$NonNV_prop <- "Non"
  }
  
  
  output$Resp <- FirstSecondType
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
                                                                                 FirstSecondType, AIMType, condition)}
stopCluster(cl)
proc.time()-t
  
name.output <- "result/Fig2/FirstSecond.AIM.csv"
write.csv(output, name.output, row.names = FALSE)
