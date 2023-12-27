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
dir.input <- "tmp/result/intermediate/5_Tet/JoinTP_DifAbund_Tet"
name.output.table <- "tmp/result/Fig5/tetramer_table.csv"
name.output <- "tmp/result/Fig5/tetramer_frequency.csv"
cores <- 12
params <- c("TP1", "TP2", "TP3", "TP4", "TP5", "TP6", "TP8", "TP9", "TP10", "TP11")
epitope_query <-c("S269", "S448", "S919", "S1208")

#################################### Processing layer ####################################################
###Define functions
#Table.read.fast
tableread_fast = function(i, header=TRUE, quote="", sep=","){
  tmp = fread(i, header=header, sep=sep, nThread=32)
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
  d <- tableread_fast(name.input, header = TRUE, sep=",")
  
  d_sub <- dplyr::filter(d, tetramer != "Nega")
  if(nrow(d_sub) > 0){
    d_sub$sample <- file.name
  }
  return(d_sub)
}

#################################### Main module ####################################################
files  <- list.files(dir.input, pattern="NaraCOVID")
t<-proc.time()
cl <- makeCluster(cores)
registerDoParallel(cl)   
d <- foreach(file.name = files, .combine = rbind,
               .packages=c("ggplot2", "extrafont", "stringr", "dplyr", "data.table")) %dopar% {Main(file.name, dir.input, dir.output)}
stopCluster(cl)
proc.time()-t

#Output tetramer+ clone table
write.csv(d, name.output.table, row.names = FALSE)

#Summarize
participants <- unique(d$sample)

output.all <- data.frame()
for(epitope in epitope_query){
  #Extract table of each epitope
  TPs <- names(d)[9:18]
  d_sub <- dplyr::select(d, c(TPs, epitope, sample))
  names(d_sub) <- c(TPs, "tetramer", "sample")
  d_sub <- dplyr::filter(d_sub, tetramer == "Posi")
    
  #Prepare output data frame
  output <- createEmptyDf(length(participants), length(params), colnames = params )
  
  #Check tcrdist hit for each participant
  for(i in 1:length(participants)){
    id <- participants[i]
    d_sub2 <- dplyr::filter(d_sub, str_detect(sample, id))
    if(nrow(d_sub > 0)){
      output[i,] <- dplyr::select(d_sub2, params) %>% apply(2, sum)
    } else {
      output[i,] <- 0
    } 
  }
  #Only export donors with specific clones
  output$ID <- participants
  output$detect <- apply(dplyr::select(output, contains("TP")), 1, sum)
  output <- dplyr::filter(output, detect > 0) %>% dplyr::select(-c("detect"))
  #Judge whether donors were responded to the epitope or not
  responder_judge <- function(x){
    out <- "Non"
    if(max(x[2], x[3], x[7]) > max(0.0005, 5*x[1])){
      out <- "Exp"
    }
    return(out)
  }
  output$classify <- apply(dplyr::select(output, contains("TP")), 1, responder_judge)
  
  output$epitope <- epitope

  output.all <- rbind(output.all, output)
}

write.csv(output.all, name.output, row.names = FALSE)
