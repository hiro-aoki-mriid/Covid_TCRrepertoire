### Summarize AIM positivity within each WGCNA module

#load libraries
library(foreach)
library(doParallel)
library(data.table)
library(stringr)
library(dplyr)
library(ggplot2)

#Input layer
dir.input <- 'WGCNA_AIMcombined'
dir.output <- "."
cores <- 12
Tcell <- "CD4" #CD4 or CD8
wgcna.params <- c("singlet", "0", "green", "invgreen", "brown", "invbrown", "yellow", "invyellow",
                  "turquoise", "invturquoise", "blue", "invblue") #For CD4

#################################### Processing layer ####################################################
###Define functions
createEmptyDf = function( nrow, ncol, colnames = c() ){
  data.frame( matrix( vector(), nrow, ncol, dimnames = list( c(), colnames ) ) )
}
tableread_fast = function(i, header=TRUE, quote="", sep="\t"){
  tmp = fread(i, header=header, sep=sep, quote=quote, nThread=32)
  tmp = as.data.frame(tmp)
  return(tmp)
}
##Main functions
WGCNA_AIM_summarize = function(file.name, dir.input, types, params, wgcna.params){
  #load
  name.input <- str_c(dir.input, file.name, sep = "/")
  d_wgcna <- tableread_fast(name.input, header = TRUE, quote="\"", sep=",")
  d_wgcna$unique <- rep(1, nrow(d_wgcna))
  
  ##Calculate stats for denominator of relative values
  #Calculate total freq. of AIM, NonNV, non at each point
  d_wgcna <- group_by(d_wgcna, class.y)
  total_params <- summarise(d_wgcna, unique = sum(unique), P1 = sum(P1), P2 = sum(P2), P3 = sum(P3), P4 = sum(P4))
  d_wgcna <- ungroup(d_wgcna)
  
  #Prepare dataframe for output
  summary.wgcna <- createEmptyDf(length(params)*length(wgcna.params), length(types), colnames = types)
  summary.wgcna.relative <- createEmptyDf(length(params)*length(wgcna.params), length(types), colnames = types)
  #Summarize total freq or clone count
  row.vecs <- vector()
  for(j in 1:length(wgcna.params)){
    k <- wgcna.params[j]
    d_wgcna_sub <- dplyr::filter(d_wgcna, class.x==k)
    for(i in 1:length(types)){
      l <- types[i]
      d_calc <- dplyr::filter(d_wgcna_sub, class.y==l)
      d_calc <- dplyr::select(d_calc, params)
      #Raw values
      sum <- apply(d_calc, 2, sum)
      summary.wgcna[((j-1)*5+1):(j*5),i] <- sum
      #Relative values
      denominators <- dplyr::filter(total_params, class.y == l)
      sum_relatives <- as.numeric(as.vector(sum) / as.vector(denominators[,2:6]))
      summary.wgcna.relative[((j-1)*5+1):(j*5),i] <- sum_relatives
      #recall parameter vectors 
      row.vec <- str_c(k, params, sep = "_")
    }
    row.vecs <- c(row.vecs, row.vec)
  }
  
  #Merge output dataframes
  names(summary.wgcna.relative) <- str_c(names(summary.wgcna.relative), "relative", sep = "_")
  summary.wgcna <- cbind(summary.wgcna, summary.wgcna.relative)
  summary.wgcna$param <- row.vecs
  summary.wgcna$sample <- file.name
  
  return(summary.wgcna)
}

###Main modules
dir.create(dir.output)
#Set parameters for summarize
types <- c("AIM", "NonNV", "non")
params <- c("unique", "P1", "P2", "P3", "P4")

files  <- list.files(dir.input, pattern=Tcell)
t<-proc.time()
cl <- makeCluster(cores)
registerDoParallel(cl)   
out <- foreach(file.name = files, .combine = "rbind",
               .packages=c("stringr", "data.table", "dplyr")) %dopar% {WGCNA_AIM_summarize(file.name, dir.input, types, params, wgcna.params)}
stopCluster(cl)
proc.time()-t

#Output
name.output <- str_c(dir.output, Tcell, sep="/") %>% str_c("AIM.wgcna.summary.csv", sep = ".")
write.csv(out, name.output, row.names = FALSE, quote = F)

