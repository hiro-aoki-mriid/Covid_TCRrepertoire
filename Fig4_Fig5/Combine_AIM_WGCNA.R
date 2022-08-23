###Determine kinetics of AIM+ clones based on the WGCNA results

#load libraryies
library(ggplot2)
library(foreach)
library(doParallel)
library(stringr)
library(dplyr)
library(data.table)

#Input layer
dir.input.WGCNA <- "Fig1_Fig2/WGCNA_output_150/Correlation.Table"
dir.input.AIM <- "Fig4_Fig5/AIM_Table"
dir.output <- "Fig4_Fig5/WGCNA_AIMcombined"
extract <- c("cdr3nt", "v", "j", #Definition of clones
             "P1", "P2", "P3", "P4", #Frequency at each point
             "class") #classification of clones
cores <- 12

#################################### Processing layer ####################################################
###Define functions
tableread_fast = function(i, header=TRUE, quote="", sep="\t"){
  tmp = fread(i, header=header, sep=sep, quote=quote, nThread=32)
  tmp = as.data.frame(tmp)
  return(tmp)
}

##Combine WGCNA and AIM informations
Combine <- function(file.name, extract, dir.input.AIM, dir.input.WGCNA){ 
  #load data
  name.input <- str_c(dir.input.WGCNA, file.name, sep = "/")
  WGCNA <- tableread_fast(name.input, header = TRUE, quote="\"", sep=",")
  name.input <- str_c(dir.input.AIM, file.name, sep = "/") %>% 
    str_replace("NaraCOVID", "COVIDAIM") %>% str_replace("wgcna.csv", "th.4.AIMposi.csv")
  if(file.exists(name.input) == TRUE){
    AIM <- tableread_fast(name.input, header = TRUE, quote="\"", sep=",")
    
    ##Extract AIM+ clones from WGCNA dataset
    #Define clonotype by CDR3nt, V usage, and J usage
    AIM$query <- str_c(AIM$cdr3nt, AIM$v, AIM$j, sep = "_")
    AIM <- dplyr::select(AIM, c("class", "query"))
    WGCNA$query <- str_c(WGCNA$cdr3nt, WGCNA$v, WGCNA$j, sep = "_")
    #Merge AIM information
    WGCNA_AIM <- merge(WGCNA, AIM, by = "query", all.x = TRUE)
    WGCNA_AIM[is.na(WGCNA_AIM)] <- "non"
    
    #Output
    name.output <- str_c(dir.output, file.name, sep = "/") %>% str_replace("wgcna", "wgcna.aim")
    write.csv(WGCNA_AIM, name.output, row.names = FALSE)
  }
}

###Main module
setwd("../")
dir.create(dir.output)
files  <- list.files(dir.input.WGCNA, pattern=".csv")

##Combine WGCNA and AIM informations
t<-proc.time()
cl <- makeCluster(cores)
registerDoParallel(cl)   
foreach(file.name = files,
               .packages=c("data.table", "stringr", "dplyr", "ggplot2"),
               .combine = rbind) %dopar% {Combine(file.name, extract, dir.input.AIM, dir.input.WGCNA)}
stopCluster(cl)
proc.time()-t

setwd("Fig4_Fig5")
