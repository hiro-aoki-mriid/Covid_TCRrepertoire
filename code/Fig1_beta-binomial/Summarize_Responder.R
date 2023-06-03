##Merge different abundance clone into tcr repertoire dataset

library(foreach)
library(doParallel)
library(stringr)
library(data.table)
library(dplyr)

dir.input <- "tmp/result/intermediate/1_beta-binomial/JoinTP_DifAbund"

#Output params
Point_Dif <- c("TP1_TP2", "TP2_TP3", "TP6_TP8")

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
  
  ###Summarize
  ##Count clones
  #Create output data.frame
  Resp <- vector()
  output <- createEmptyDf(length(Point_Dif), 1, colnames = "Extent.Expansion")
  for(i in 1:length(Point_Dif)){
    p_dif <- Point_Dif[i]
    resp <- str_remove(p_dif, "_") %>% str_c("dif")
    Resp <- c(Resp, resp)
    #Extract freq and Resp timepoint
    dsub <- dplyr::select(d, c(p_dif, resp))
    names(dsub) <- c("p_dif", "resp")
    #Only Expanded clones
    dsub <- dplyr::filter(dsub, resp == "Exp")
    #Freq at each time point
    if(nrow(dsub) > 0){
      output[i, 1] <- sum(dsub$p_dif)
    } else{
      output[i, 1] <- 0
    }
  }
  output$Resp <- Resp
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
  
name.output <- "result/Fig1/1C.Responder.Freq.csv"
write.csv(output, name.output, row.names = FALSE)