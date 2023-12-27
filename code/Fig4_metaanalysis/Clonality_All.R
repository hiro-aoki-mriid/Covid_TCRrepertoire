##Merge different abundance clone into tcr repertoire dataset

library(foreach)
library(doParallel)
library(stringr)
library(data.table)
library(dplyr)

dir.input <- "tmp/data/original/raw.data"
dir.output <- "tmp/result/Fig4"

cores <- 24

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

############################### Functions ########################################
##Extract beta-binomial test result
Extract <- function(file.name, dir.input, FirstSecondType, ThirdType){
  #load repetoire data
  name.input <- str_c(dir.input, file.name, sep = "/")
  d <- tableread_fast(name.input, header = TRUE, sep="\t")
  
  #Calculate clonality and evenness
  pi <- CalcPielou(as.vector(d$freq))
  de <- DE90(as.vector(d$freq))

  output <- c(pi, de, file.name)
  return(output)
}

############################### Processing ########################################
dir.create(dir.output, recursive = TRUE)
files  <- list.files(dir.input, pattern="NaraCOVID_TP1_")
  
t<-proc.time()
cl <- makeCluster(cores)
registerDoParallel(cl)   
output <- foreach(file.name = files, .combine = rbind,
                 .packages=c("data.table", "stringr", "dplyr")) %dopar% {Extract(file.name, dir.input,
                                                                                 FirstSecondType, ThirdType)}
stopCluster(cl)
proc.time()-t

output <- as.data.frame(output)
names(output) <- c("Pielou", "DE90", "file.name")
  
name.output <- str_c(dir.output, "Clonality_TP1.csv", sep = "/")
write.csv(output, name.output, row.names = FALSE)
