###Fig 4I_5D: Extratract Top10 dominant clones among AIM+ clones at P2, P3 and P4

#load libraryies

library(ggplot2)
library(extrafont)
loadfonts("win")
library(foreach)
library(doParallel)
library(stringr)
library(dplyr)
library(data.table)
library(tidyverse)

#Input layer
dir.input <- "WGCNA_AIMcombined"
dir.output <- "4I_5D"
rank <- 10 #examine top X clones
classes <- c("AIM", "NonNV")
cores <- 12

#################################### Processing layer ####################################################
###Define functions
#Table.read.fast
tableread_fast = function(i, header=TRUE, quote="", sep="\t"){
  tmp = fread(i, header=header, sep=sep, quote=quote, nThread=32)
  tmp = as.data.frame(tmp)
  return(tmp)
}

createEmptyDf = function( nrow, ncol, colnames = c() ){
  data.frame( matrix( vector(), nrow, ncol, dimnames = list( c(), colnames ) ) )
}

#Main function
Main <- function(file.name, dir.input, class){
  ##load data
  #load data
  name.input <- str_c(dir.input, file.name, sep = "/")
  d <- tableread_fast(name.input, header = TRUE, quote="\"", sep=",")
  #Select only AIM+ clones
  d <- dplyr::filter(d, class.y == class)
  
  if(nrow(d) > rank){
    #re-order at P2
    d <- dplyr::arrange(d, -P2)
    d$P2dom <- c(rep(1, rank), rep(0, (nrow(d)-rank)))
    #re-order at P3
    d <- dplyr::arrange(d, -P3)
    d$P3dom <- c(rep(1, rank), rep(0, (nrow(d)-rank)))  
    #re-order at P4
    d <- dplyr::arrange(d, -P4)
    d$P4dom <- c(rep(1, rank), rep(0, (nrow(d)-rank))) 
    
    result <- createEmptyDf(3, 3, colnames = c("P2dom", "P3dom", "P4dom"))
    
    #Summarize for top10 clones at P2 or P3 or P4
    d_sub <- dplyr::filter(d, P2dom == 1)
    result[1,] <- summarise(d_sub, P2dom = sum(P2dom), P3dom = sum(P3dom), P4dom = sum(P4dom))
    d_sub <- dplyr::filter(d, P3dom == 1)
    result[2,] <- summarise(d_sub, P2dom = sum(P2dom), P3dom = sum(P3dom), P4dom = sum(P4dom))
    d_sub <- dplyr::filter(d, P4dom == 1)
    result[3,] <- summarise(d_sub, P2dom = sum(P2dom), P3dom = sum(P3dom), P4dom = sum(P4dom))
    result$Point <- c("P2dom", "P3dom", "P4dom")
    result$name <- file.name
    
    #Output
    return(result)
  }
}

###Main module
files  <- list.files(dir.input, pattern="csv")
dir.create(dir.output)
for(class in classes){
  t<-proc.time()
  cl <- makeCluster(cores)
  registerDoParallel(cl)    
  out <- foreach(file.name = files,
                 .combine = rbind, .packages=c("ggplot2", "extrafont", "stringr", "dplyr", "data.table")) %dopar% {Main(file.name, dir.input, class)}
  stopCluster(cl)
  proc.time()-t
  
  out <- as.data.frame(out)
  name.output <- str_c(dir.output, class, sep = "/") %>% str_c("Top10.clones.dominance.csv", sep = ".")
  write.csv(out, name.output, row.names = FALSE)
}

