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
name.output <- "tmp/result/Fig1/responder_clones.csv"
dir.input <- "tmp/result/intermediate/1_beta-binomial/JoinTP_DifAbund"
cores <- 12
patterns <- c("2nd", "1st", "3rd", "Dual")

#################################### Processing layer ####################################################
###Define functions
#Table.read.fast
tableread_fast = function(i, header=TRUE, quote="", sep=","){
  tmp = fread(i, header=header, sep=sep, quote=quote, nThread=32)
  tmp = as.data.frame(tmp)
  return(tmp)
}

### Merge scTCR clone table ###
##Combine beta-binomial result and AIM information
Combine <- function(file.name, extract, dir.AIM, dir.input){ 
  #load data
  name.input <- str_c(dir.input, file.name, sep = "/")
  data <- tableread_fast(name.input, header = TRUE, quote="\"", sep=",")
  data <- dplyr::filter(data, Resp %in% patterns)
  data$name <- file.name
  
  return(data)
}

###Main module
files  <- list.files(dir.input, pattern="NaraCOVID")

##Combine beta-binomial result and AIM informations
t<-proc.time()
cl <- makeCluster(cores)
registerDoParallel(cl)   
out <- foreach(file.name = files,
        .packages=c("data.table", "stringr", "dplyr", "ggplot2"),
        .combine = rbind) %dopar% {Combine(file.name, extract, dir.AIM, dir.input)}
stopCluster(cl)
proc.time()-t

write.csv(out, name.output, row.names = FALSE)



