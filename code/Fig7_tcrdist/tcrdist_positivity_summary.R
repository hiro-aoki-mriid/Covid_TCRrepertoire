###Extract SARS-Cov-2 reactive clones determined by tcrdist3

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
dir.input <- "tmp/result/intermediate/7_tcrdist/JoinTP_DifAbund_Tet_tcrdist"
dir.output <- "tmp/result/Fig7"
name.output <- "tmp/result/Fig7/tcrdist_clones.csv"
cores <- 12

#################################### Processing layer ####################################################
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

#Main function
Main <- function(file.name, dir.input, dir.output){
  ##load data
  #load AIM data
  name.input <- str_c(dir.input, file.name, sep = "/")
  d <- tableread_fast(name.input, header = TRUE, quote="\"", sep=",")
  
  d_sub <- dplyr::filter(d, tetramer != "Nega" |
                           motif != "Nega")
  if(nrow(d_sub) > 0){
    d_sub$sample <- file.name
  }
  return(d_sub)
}

############################### Processing ########################################
dir.create(dir.output, recursive = TRUE)
files  <- list.files(dir.input, pattern="NaraCOVID")
t<-proc.time()
cl <- makeCluster(cores)
registerDoParallel(cl)   
d_difabund <- foreach(file.name = files, .combine = rbind,
               .packages=c("ggplot2", "extrafont", "stringr", "dplyr", "data.table")) %dopar% {Main(file.name, dir.input, dir.output)}
stopCluster(cl)
proc.time()-t
d_difabund <- as.data.frame(d_difabund)

write.csv(d_difabund, name.output, row.names = FALSE)
