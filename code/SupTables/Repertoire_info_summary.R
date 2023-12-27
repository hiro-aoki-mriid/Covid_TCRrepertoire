##Merge different abundance clone into tcr repertoire dataset

library(ggplot2)
library(foreach)
library(doParallel)
library(stringr)
library(data.table)
library(dplyr)

dir.input <- "tmp/data/original/raw.data"
dir.output <- "tmp/result/SupTables"
dir.input <- "data/original/raw.data"
dir.output <- "result/SupTables"
name.output <- "clone.count.table.csv"
cores <- 24

############################### Functions ########################################
tableread_fast = function(i, header=TRUE, quote="", sep="\t"){
  tmp = fread(i, header=header, sep=sep, quote=quote, nThread=32)
  tmp = as.data.frame(tmp)
  return(tmp)
}
createEmptyDf = function( nrow, ncol, colnames = c() ){
  data.frame( matrix( vector(), nrow, ncol, dimnames = list( c(), colnames ) ) )
}

##Extract beta-binomial test result
Extract <- function(name){
  #load repetoire data
  name.input <- str_c(dir.input, file.name, sep = "/")
  d <- tableread_fast(name.input, header = TRUE, quote="", sep="\t")
  
  out <- c(nrow(d), name.input)

  #Output
  return(out)
}

############################### Processing ########################################
dir.create(dir.output, recursive = TRUE)

files <- list.files(dir.input, pattern="COVID")
t<-proc.time()
cl <- makeCluster(cores)
registerDoParallel(cl)   
output <- foreach(file.name = files, .combine = rbind,
                    .packages=c("data.table", "stringr", "dplyr")) %dopar% {Extract(file.name)}
stopCluster(cl)
proc.time()-t

name.output <- str_c(dir.output, name.output, sep = "/")
write.csv(output, name.output, row.names = FALSE)
