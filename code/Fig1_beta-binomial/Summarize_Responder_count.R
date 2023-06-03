##Merge different abundance clone into tcr repertoire dataset

library(foreach)
library(doParallel)
library(stringr)
library(data.table)
library(dplyr)

dir.input <- "tmp/result/intermediate/1_beta-binomial/JoinTP_DifAbund"

#Output params

FirstSecondType <- c("1st", "1stInc", "2nd", "2ndInc", "Dual", "DualInc", "3rd", "3rdInc", "Others")
RespPoint <- c("TP1TP2dif", "TP2TP3dif", "TP6TP8dif")
TP_anal <- c("TP2", "TP3", "TP8") 

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
Extract <- function(file.name, dir.input, FirstSecondType, RespPoint){
  #load repetoire data
  name.input <- str_c(dir.input, file.name, sep = "/")
  d <- tableread_fast(name.input, header = TRUE, sep=",")
  
  ###Summarize
  ##Count clones
  #Create output data.frame
  count.all <- nrow(d)
  output <- createEmptyDf(1, length(RespPoint), colnames = RespPoint)
  for(j in 1:length(RespPoint)){
    tt <- RespPoint[j]
    point_anal <- TP_anal[j]
    d_anal <- dplyr::select(d, c(tt, point_anal))
    names(d_anal) <- c("DifPoint", "AnalPoint")
    count.all <- sum(d_anal$AnalPoint > 0)
    #Count of clones
    output[1,j] <- nrow(dplyr::filter(d_anal, DifPoint == "Exp"))/count.all
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
                                                                                 FirstSecondType, RespPoint)}
stopCluster(cl)
proc.time()-t
  
name.output <- "result/Fig1/1D.Responder.countprop.csv"
write.csv(output, name.output, row.names = FALSE)
