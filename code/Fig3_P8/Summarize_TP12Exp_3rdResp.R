##Merge different abundance clone into tcr repertoire dataset

library(foreach)
library(doParallel)
library(stringr)
library(data.table)
library(dplyr)

dir.input <- "tmp/result/intermediate/1_beta-binomial/JoinTP_DifAbund"
dir.output <- "tmp/result/Fig3"

#Output params

Types <- c("main") 
freq_params <- c("TP2", "TP3", "TP6", "TP8")
count_params <- c("TP3dif", "TP8dif")
clonality_params <- c("Pielou_TP3", "Pielou_TP8", "DE90_TP3", "DE90_TP8")

cores <- 24

############################### Functions ########################################
tableread_fast = function(i, header=TRUE, sep="\t"){
  tmp = fread(i, header=header, sep=sep, nThread=32)
  tmp = as.data.frame(tmp)
  return(tmp)
}
createEmptyDf = function( nrow, ncol, colnames = c() ){
  data.frame( matrix( vector(), nrow, ncol, dimnames = list( c(), colnames ) ) )
}
#Calculating 1-Pielou index
CalcPielou <- function(d_vec){
  #Normalize
  d_vec <- d_vec / sum(d_vec)
  #Calculation
  pi <- 1 + sum(d_vec * log(d_vec))/log(length(d_vec))
}
#Calculate DE90 diversity index
DE90 <- function(d_vec){
  #Normalize
  d_vec <- d_vec / sum(d_vec)
  #Sort
  d_vec <- sort(d_vec, decreasing = TRUE)
  #Calculate
  sum <- 0
  i <- 0
  while (sum <= 0.9) {
    i <- i+1
    sum <- sum + d_vec[i] 
  }
  DEX <- i/length(d_vec)
}

##Extract beta-binomial test result
Extract <- function(file.name, dir.input, FirstSecondType, ThirdType){
  #load repetoire data
  name.input <- str_c(dir.input, file.name, sep = "/")
  d <- tableread_fast(name.input, header = TRUE, sep=",")
  TPs <- c("TP1", "TP2", "TP3", "TP4", "TP5", "TP6", "TP8", "TP9", "TP10", "TP11")
  
  ###Summarize
  ##Count clones
  #Create output data.frame
  output.all <- data.frame()
  all_params <- c(freq_params, count_params, clonality_params)
  output <- createEmptyDf(1, length(all_params), colnames = all_params)

  #Only analyze main responder clones
  dsub <- dplyr::filter(d, Resp == "Main")
    
  #Calculate total freq. or freq. diff
  dsub_freq <- dplyr::select(dsub, freq_params)
  output[1, 1:(length(freq_params))] <- apply(dsub_freq, 2, sum)
    
  #Calculate count of expanding clones
  dsub_exp <- dplyr::select(dsub, count_params)
  count_exp <- function(vec){sum(vec == "Exp")}
  output[1, (1+length(freq_params)):length(c(freq_params, count_params))] <- apply(dsub_exp, 2, count_exp)
    
  #Calculate clonality
  #Only analyze clonality when count of clones at TP8 or TP8 is more than 10.
  TP3_count <- sum(dsub$TP3 > 0)
  TP8_count <- sum(dsub$TP8 > 0)
  if(TP3_count > 9 & TP8_count >9){
    dsub_TP3 <- dplyr::select(dsub, c("TP3")) %>% dplyr::filter(TP3 > 0)
    names(dsub_TP3) <- c("Freq")
    dsub_TP8 <- dplyr::select(dsub, c("TP8")) %>% dplyr::filter(TP8 > 0)
    names(dsub_TP8) <- c("Freq")
      
    #Calculate 1-Pielou index for Third-expanded 2nd responder
    pi_TP3 <- CalcPielou(as.vector(dsub_TP3$Freq))
    pi_TP8 <- CalcPielou(as.vector(dsub_TP8$Freq))
    de_TP3 <- DE90(as.vector(dsub_TP3$Freq))
    de_TP8 <- DE90(as.vector(dsub_TP8$Freq))
    output[1, (1+length(c(freq_params, count_params))):length(all_params)] <- c(pi_TP3, pi_TP8, de_TP3, de_TP8)
  } else {
    output[1, (1+length(c(freq_params, count_params))):length(all_params)] <- rep("NA", 4)
  }
  
  output$sample <- file.name
  return(output)
}

############################### Processing ########################################
dir.create(dir.output, recursive = TRUE)
files  <- list.files(dir.input, pattern="NaraCOVID")
  
t<-proc.time()
cl <- makeCluster(cores)
registerDoParallel(cl)   
output <- foreach(file.name = files, .combine = rbind,
                 .packages=c("data.table", "stringr", "dplyr")) %dopar% {Extract(file.name, dir.input,
                                                                                 FirstSecondType, ThirdType)}
stopCluster(cl)
proc.time()-t
  
name.output <- str_c(dir.output, "3A-D.MainResp_TP2vsTP8.csv", sep = "/")
write.csv(output, name.output, row.names = FALSE)
