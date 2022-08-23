###Extract clones which correlated with WGCNA modules

#Load Libraries
library(foreach)
library(doParallel)
library(dplyr)
library(WGCNA)
library(flashClust)
options(stringsAsFactors=FALSE)
enableWGCNAThreads(14)            #maximum number of CPU threads of your PC
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
library(stringr)
library(data.table)
library(hfunk) #https://github.com/hstojic/hfunk

#Input layer
dir.input <- "S1BCEF"
Tcells <- c("CD4", "CD8") #Whether CD4 or CD8 repertoire is analysed
dir.output <- "WGCNA_output"
#Parameters for WGCNA analysis
topX <- 150 #Count of clones used for WGCNA for each repertoire
cores <- 14
ppi <- 600
th.correlation <- 0.85 #threshold for clone detection which correlated with WGCNA modules

#################################### Processing layer ####################################################
dir.output.raw <- str_c(dir.output, topX, sep = "_")
dir.create(dir.output.raw)

##### Prepare functions #####
#Create empty data.frame
createEmptyDf = function( nrow, ncol, colnames = c() ){
  data.frame( matrix( vector(), nrow, ncol, dimnames = list( c(), colnames ) ) )
}

#Table.read.fast
tableread_fast = function(i, header=TRUE, quote="", sep="\t"){
  tmp = fread(i, header=header, sep=sep, quote=quote, nThread=32)
  tmp = as.data.frame(tmp)
  return(tmp)
}

#stop parallel computing
unregister_dopar <- function() {
  env <- foreach:::.foreachGlobals
  rm(list=ls(name=env), pos=env)
}

#Extract clones which correlated with WGCNA modules
WGCNA.extract <- function(file.name, wgcna.result, dir.output.raw, th.correlation){ 
  #load repertoire data
  name.input <- str_c(dir, file.name, sep = "/")
  d_raw <- tableread_fast(name.input, header = TRUE, quote="\"", sep="\t")
  names(d_raw) <- c("count", "freq",  "cdr3nt", "cdr3aa", "v", "d", "j", "VEnd", "DStart", "DEnd", "JStart",
                "peak", "occurences", "sampling.p","P1", "P2", "P3", "P4")
  
  # Calculate correlation coefficient with clones detected over 2 time points
  d_sub <- dplyr::filter(d_raw, occurences >= 2)
  d_TCR <- dplyr::select(d_sub, c("cdr3nt", "cdr3aa", "v", "d", "j"))
  d_sub <- dplyr::select(d_sub, c("P1", "P2", "P3", "P4"))
  modules <- names(wgcna.result)
  wgcna.cor = createEmptyDf(nrow(d_sub), length(modules), colnames = modules )
  for(i in 1:length(modules)){
    mod <- modules[i]
    mod.val <- wgcna.result[,i]
    calc.cor <- function(x) {return (cor(x, mod.val, method="pearson"))}
    wgcna.cor[,i] <- apply(d_sub, 1, calc.cor)
  }
  
  # Define module which each clone is most correlated
  classdef <- function(x) {
    if(max(x) < th.correlation){
      return(0)
    }else{
      idx <- idxMax(x)
      return(modules[idx])
    }}
  class <- apply(wgcna.cor, 1, classdef)
  d <- cbind(d_sub, class, d_TCR)
  
  #Merge data of singlet clones
  d_singlet <- dplyr::filter(d_raw, occurences == 1)
  d_sub <- dplyr::select(d_singlet, c("P1", "P2", "P3", "P4"))
  class <- rep("singlet", nrow(d_singlet))
  d_TCR <- dplyr::select(d_singlet, c("cdr3nt", "cdr3aa", "v", "d", "j"))
  d_singlet <- cbind(d_sub, class, d_TCR)
  d <- rbind(d, d_singlet)
  
  #Output correlation table
  dir.create(str_c(dir.output.raw, "Correlation.Table", sep = "/"))
  name.out <- str_c(dir.output.raw, "Correlation.Table", file.name, sep = "/") %>% str_replace("join.strict.table.txt", "wgcna.csv")
  fwrite(d, name.out, row.names = FALSE)
  
  ### Summarize
  #count of clones correlated with each module
  counts <- table(d$class)
  summary.table <- as.data.frame(counts)
  names(summary.table) <- c("class", "counts")
  summary.table$prop <- summary.table$counts / sum(summary.table$counts)
  #frequency of clones correlated with each module
  sum_TP1 <- as.vector(as.numeric(tapply(d$P1, d$class, sum)))
  sum_TP2 <- as.vector(as.numeric(tapply(d$P2, d$class, sum)))
  sum_TP3 <- as.vector(as.numeric(tapply(d$P3, d$class, sum)))
  sum_TP4 <- as.vector(as.numeric(tapply(d$P4, d$class, sum)))
  
  summary.table$sum_TP1 <- sum_TP1
  summary.table$sum_TP2 <- sum_TP2
  summary.table$sum_TP3 <- sum_TP3
  summary.table$sum_TP4 <- sum_TP4
  
  ###Output
  summary.table$name <- file.name
  return(summary.table)
}

##### Main module #####
for(Tcell in Tcells){
  #load input data
  name.input <- str_c(dir.input, "module_result", sep = "/") %>% str_c(Tcell, topX, "csv", sep = ".")
  meord_output <- read.csv(name.input, header = TRUE)
  
  #Generate modules that inversely correlated with original eigengenes
  wgcna.result <- dplyr::select(meord_output, -sample)
  wgcna.inv <- -(wgcna.result)
  names(wgcna.inv) <- str_c("inv", names(wgcna.result), sep = "")
  wgcna.result <- cbind(wgcna.result, wgcna.inv)
  names(wgcna.result) <- str_replace(names(wgcna.result), "ME", "")
  
  files  <- list.files("Join", pattern=Tcell)
  
  #Extract correlated clones
  t<-proc.time()
  cl <- makeCluster(cores)
  registerDoParallel(cl)   
  out <- foreach(file.name = files, .combine = rbind,
                 .packages=c("stringr", "dplyr", "hfunk", "data.table")) %dopar% {WGCNA.extract(file.name, wgcna.result, dir.output.raw, th.correlation)}
  stopCluster(cl)
  proc.time()-t
  
  #Output
  out <- as.data.frame(out)
  name.out <- str_c(dir.output.raw, "WGCNA_extract", sep = "/") %>% str_c(Tcell, "th.correlation", th.correlation, "csv", sep = ".")
  write.csv(out, name.out, row.names = F)
  
}


