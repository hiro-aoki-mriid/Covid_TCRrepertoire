#Prepare input file for searching clones that match metaclonotype

#load libraryies
library(dplyr)
library(stringr)
library(foreach)
library(doParallel)
library(data.table)

#Input layer
dir.input <- "tmp/result/intermediate/7_tcrdist/Metaclonotype_query"
dir.output <- "tmp/result/intermediate/7_tcrdist/Metaclonotype_query/participants"
name.hla.table <- "tmp/metadata/HLAtyping_result.csv"
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
Main <- function(file.name, dir.input, dir.output){
  ##load data
  name.input <- str_c(dir.input, file.name, sep = "/")
  table <- tableread_fast(name.input, header = TRUE, quote="\"", sep="\t")
  hla_table <- read.csv(name.hla.table, header = TRUE, colClasses=c("character", "character"))
  
  for(i in 1:nrow(hla_table)){
    id <- hla_table[i,1]
    
    hla <- as.vector(str_split(hla_table[i,2], pattern = "_", simplify = TRUE))
    table_sub <- dplyr::filter(table, HLA %in% hla)
    name_output <- str_c(dir.output, id, sep = "/") %>% str_c(names(hla_table)[2], file.name, sep = ".")
    write.table(table_sub, name_output, sep = "\t", row.names = FALSE, quote = FALSE)
  }
}

###Main module
files  <- list.files(dir.input, pattern="tsv")
dir.create(dir.output)
t<-proc.time()
cl <- makeCluster(cores)
registerDoParallel(cl)   
out <- foreach(file.name = files,
               .combine = rbind, .packages=c("ggplot2", "extrafont", "stringr", "dplyr", "data.table")) %dopar% {Main(file.name, dir.input, dir.output)}
stopCluster(cl)
proc.time()-t

