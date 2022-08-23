###Calculate frequency of AIM+ clones at each timepoint

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
dir.input.bulk <- "Fig1_Fig2/Join"
dir.input.AIM <- "Fig4_Fig5/AIM_Table"
dir.output <- "Fig4_Fig5/4D_S3DE"
cores <- 12

#################################### Processing layer ####################################################
###Define functions
#Table.read.fast
tableread_fast = function(i, header=TRUE, quote="", sep="\t"){
  tmp = fread(i, header=header, sep=sep, quote=quote, nThread=32)
  tmp = as.data.frame(tmp)
  return(tmp)
}

#Main function
Main <- function(file.name, dir.input.AIM, dir.input.bulk){
  ##load data
  #load AIM data
  name.input <- str_c(dir.input.AIM, file.name, sep = "/")
  AIM <- tableread_fast(name.input, header = TRUE, quote="\"", sep=",")
  #load bulk TCR repertoire data
  name.input <- str_c(dir.input.bulk, file.name, sep = "/") %>% 
    str_replace("COVIDAIM", "NaraCOVID") %>% str_replace("th.4.AIMposi.csv", "join.strict.table.txt")
  Bulk <- tableread_fast(name.input, header = TRUE, quote="\"", sep="\t")
  names(Bulk) <- c("count", "freq",  "cdr3nt", "cdr3aa", "v", "d", "j", "VEnd", "DStart", "DEnd", "JStart",
                "peak", "occurences", "sampling.p","P1", "P2", "P3", "P4")

  ##Extract AIM+ clones from Bulk dataset
  #Prepare dataset
  AIM <- dplyr::filter(AIM, class == "AIM")
  #Define clonotype by CDR3nt, V usage, and J usage
  AIM$query <- str_c(AIM$cdr3nt, AIM$v, AIM$j, sep = "_")
  Bulk$query <- str_c(Bulk$cdr3nt, Bulk$v, Bulk$j, sep = "_")
  #Search dataset
  Bulk_AIM <- dplyr::filter(Bulk, query %in% unique(AIM$query)) %>% dplyr::select(c("P1", "P2", "P3", "P4"))
  
  ##Calculate frequency of AIM+ clones at each timepoint
  sum_table <- apply(Bulk_AIM, 2, sum)
  sum_table <- c(file.name, as.vector(sum_table), nrow(Bulk_AIM))
  
  return(sum_table)
}

###Main module
setwd("../")
files  <- list.files(dir.input.AIM, pattern="AIMposi.csv")
dir.create(dir.output)
t<-proc.time()
cl <- makeCluster(cores)
registerDoParallel(cl)   
out <- foreach(file.name = files,
        .combine = rbind, .packages=c("ggplot2", "extrafont", "stringr", "dplyr", "data.table")) %dopar% {Main(file.name, dir.input.AIM, dir.input.bulk)}
stopCluster(cl)
proc.time()-t

out <- as.data.frame(out)
names(out) <- c("name", "freq.P1", "freq.P2", "freq.P3", "freq.P4", "clone.count")
name.output <- str_c(dir.output, "AIM_frequency.csv", sep = "/")
write.csv(out, name.output, row.names = FALSE)

setwd("Fig4_Fig5")