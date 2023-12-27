##Count and calculate total freq. of expanded clones

library(foreach)
library(doParallel)
library(stringr)
library(data.table)
library(dplyr)

dir.input <- "tmp/result/intermediate/1_beta-binomial/JoinTP_DifAbund"
dir.output <- "tmp/result/Fig1"
file.output <- "1C.Responder.Freq.csv"

#Output parameters
Point_Dif <- c("TP2dif", "TP3dif", "TP4dif", "TP5dif", "TP6dif", "TP8dif", "TP9dif", "TP10dif", "TP11dif")

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
Extract <- function(file.name, dir.input){
  #load repetoire data
  name.input <- str_c(dir.input, file.name, sep = "/")
  d <- tableread_fast(name.input, header = TRUE, sep=",")
  
  ###Summarize
  ##Count clones
  #Create output data.frame
  Resp <- vector()
  output.freq <- createEmptyDf(1, length(Point_Dif), colnames = str_c(Point_Dif, "freq", sep = "."))
  output.count <- createEmptyDf(1, length(Point_Dif), colnames = str_c(Point_Dif, "count", sep = "."))
  for(i in 1:length(Point_Dif)){
    p_dif <- Point_Dif[i]
    #Extract freq and Resp timepoint
    dsub <- dplyr::select(d, contains(p_dif))
    names(dsub) <- c("resp", "p_dif")
    #Only Expanded clones
    dsub <- dplyr::filter(dsub, resp == "Exp")
    #Freq at each time point
    if(nrow(dsub) > 0){
      output.freq[1, i] <- sum(dsub$p_dif)
      output.count[1, i] <- nrow(dsub)/nrow(d)
    } else{
      output.freq[1, i] <- 0
      output.count[1, i] <- 0
    }
  }
  output <- cbind(output.freq, output.count)
  output$sample <- file.name

  return(output)
}

############################### Processing ########################################
files  <- list.files(dir.input, pattern="NaraCOVID")

t<-proc.time()
cl <- makeCluster(cores)
registerDoParallel(cl)   
output <- foreach(file.name = files, .combine = rbind,
                 .packages=c("data.table", "stringr", "dplyr")) %dopar% {Extract(file.name, dir.input)}
stopCluster(cl)
proc.time()-t

dir.create(dir.output, recursive = TRUE)
name.output <- str_c(dir.output, file.output, sep = "/")
write.csv(output, name.output, row.names = FALSE)