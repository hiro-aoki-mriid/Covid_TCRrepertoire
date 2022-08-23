###Summarize the stat of clones used for WGCNA analysis

#Load Libraries

library(foreach)
library(doParallel)
library(stringr)
library(data.table)
library(dplyr)

#Input layer
dir <- "Join"
dir.output <- "S1A"
cores <- 12
topX <- 150 #rank thredhold sed for WGCNA analysis
name.output <- "S1A"

#################################### Processing layer ####################################################
###Define functions
#Subfunction
tableread_fast = function(i, header=TRUE, quote="", sep="\t"){
  tmp = fread(i, header=header, sep=sep, quote=quote, nThread=32)
  tmp = as.data.frame(tmp)
  return(tmp)
}

#Main function
Basic <- function(file.name){
  name.input <- str_c(dir, file.name, sep = "/")
  d <- tableread_fast(name.input, header = TRUE, quote="", sep="\t")
  names(d) <- c("count", "freq",  "cdr3nt", "cdr3aa", "v", "d", "j", "VEnd", "DStart", "DEnd", "JStart",
                "peak", "occurences", "sampling.p","P1", "P2", "P3", "P4")
  
  #Extract top150 clones
  d$rank <- c(1:nrow(d))
  d_sub <- dplyr::select(d, c("P1", "P2", "P3", "P4"))
  d_4point <- dplyr::filter(d, occurences >= 4)
  d_Top <- dplyr::filter(d, occurences >= 4 & rank <= topX) 
  d_Top_sub <- dplyr::select(d_Top, c("P1", "P2", "P3", "P4"))
  
  #Count clones at each time point
  count <- function(x) {return (sum(x>0))}
  table.count <- apply(d_sub, 2, count)
  
  #Count clones detected at four timepoint
  count.4p <- nrow(d_4point)
  
  #Calculate freq of top150 clones at each time point
  table.freq <- apply(d_Top_sub, 2, sum)
  
  ##Output data
  out <- c(table.count, count.4p, table.freq, file.name)
  
  return(out)
}

###Main module
dir.create(dir.output)
files  <- list.files(dir, pattern=".txt")

t<-proc.time()
cl <- makeCluster(cores)
registerDoParallel(cl)   
out <- foreach(file.name = files, .combine = rbind, 
               .packages=c("data.table", "stringr", "dplyr")) %dopar% {Basic(file.name)}
stopCluster(cl)
proc.time()-t

#Format output table
out <- as.data.frame(out)
names(out) <- c("P1_count", "P2_count", "P3_count", "P4_count", "count_4p",
                "P1_TopFreq", "P2_TopFreq", "P3_TopFreq", "P4_TopFreq", "name")

#Output
name.out <- str_c(dir.output, name.output, sep = "/") %>% str_c(topX, "Clone.Information.csv", sep = ".")
write.csv(out, name.out, row.names = F)


