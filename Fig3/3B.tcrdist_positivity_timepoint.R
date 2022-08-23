###Calculate frequency of SARS-Cov-2 reactive clones determined by tcrdist3 at each time point

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
dir.input <- "WGCNA_tcrdist"
dir.output <- "3B"
cores <- 12
params <- c("P1", "P2", "P3", "P4", "count")
epitopes <- c("ORF1ab", "ORF10", "M", "ORF3a", "N", "ORF7a", "ORF8", "S", "ORF7b", "E", "ORF6")

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
Main <- function(file.name, dir.input, dir.output){
  ##load data
  #load AIM data
  name.input <- str_c(dir.input, file.name, sep = "/")
  d <- tableread_fast(name.input, header = TRUE, quote="\"", sep="\t")
  d$count <- 1
  #Extract clones that assigned their epitopes by tcrdist3
  d <- dplyr::filter(d, protein != "non") %>% dplyr::select(c(params, "protein"))
  
  #Prepare dataframe for output
  summary.wgcna <- createEmptyDf(length(epitopes), length(params), colnames = params)
  #Summarize total freq or clone count
  for(j in 1:length(epitopes)){
    k <- epitopes[j]
    d_sub <- dplyr::filter(d, protein==k) %>% dplyr::select(-"protein")
    if(nrow(d_sub) > 0){
      summary.wgcna[j,] <- apply(d_sub, 2, sum)
    } else
      summary.wgcna[j,] <- c(0, 0, 0, 0, 0)
  }
  summary.wgcna$epitope <- epitopes
  summary.wgcna$name <- file.name
  
  return(summary.wgcna)
}

###Main module
files  <- list.files(dir.input, pattern="wgcna.tcrdist3.tsv")
dir.create(dir.output)
t<-proc.time()
cl <- makeCluster(cores)
registerDoParallel(cl)   
out <- foreach(file.name = files,
        .combine = rbind, .packages=c("ggplot2", "extrafont", "stringr", "dplyr", "data.table")) %dopar% {Main(file.name, dir.input, dir.output)}
stopCluster(cl)
proc.time()-t

out <- as.data.frame(out)
names(out) <- c("freq.P1", "freq.P2", "freq.P3", "freq.P4", "clone.count", "epitope", "name")
name.output <- str_c(dir.output, "3B.epitope_frequency.csv", sep = "/")
write.csv(out, name.output, row.names = FALSE)