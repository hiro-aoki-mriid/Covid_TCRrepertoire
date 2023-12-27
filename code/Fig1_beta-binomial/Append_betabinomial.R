##Merge different abundance clone into tcr repertoire dataset

library(foreach)
library(doParallel)
library(stringr)
library(data.table)
library(dplyr)

dir.binomial <- "tmp/result/intermediate/1_beta-binomial/binomial_test"
dir.repertoire <- "tmp/result/intermediate/1_beta-binomial/JoinTP"
dir.output <- "tmp/result/intermediate/1_beta-binomial/JoinTP_DifAbund"

cores <- 24

tableread_fast = function(i, header=TRUE, quote="", sep="\t"){
  tmp = fread(i, header=header, sep=sep, quote=quote, nThread=32)
  tmp = as.data.frame(tmp)
  return(tmp)
}
createEmptyDf = function( nrow, ncol, colnames = c() ){
  data.frame( matrix( vector(), nrow, ncol, dimnames = list( c(), colnames ) ) )
}

d.all <- data.frame()

############################### Functions ########################################
##Extract beta-binomial test result
Extract <- function(file.name, dir.binomial, dir.output){
  #load repetoire data
  name.input <- str_c(dir.repertoire, file.name, sep = "/")
  d <- tableread_fast(name.input, header = TRUE, quote="", sep="\t")
  #Extract only information used in analysis
  d <- dplyr::select(d, -c("VEnd", "DStart", "DEnd", "JStart", "peak", "occurences", "sampling.p"))
  names(d) <- c("count", "freq", "cdr3nt", "cdr3aa", "v", "d", "j",
                "TP1", "TP2", "TP3", "TP4", "TP5", "TP6", "TP8", "TP9", "TP10", "TP11")
  TPs <- c("TP1", "TP2", "TP3", "TP4", "TP5", "TP6", "TP8", "TP9", "TP10", "TP11")
  d$ntvj <- str_c(d$cdr3nt, d$v, d$j, sep = "_")
  
  #load beta-binomial result
  name.query <- str_remove(file.name, ".join.strict.table.txt") %>% str_remove("NaraCOVID_")
  beta.files <- list.files(dir.binomial, pattern = name.query)

  for(beta in beta.files){
    #Extract time point information
    time1 <- str_split(beta, pattern = "_", simplify = TRUE)[3]
    time2 <- str_split(beta, pattern = "_", simplify = TRUE)[9]
    
    #Load beta-binomial result
    name.input <- str_c(dir.binomial, beta, beta, sep = "/") %>% str_c(".differentialAbundance.tsv")
    d_beta <- tableread_fast(name.input, header = TRUE, quote="", sep="\t")
    d_beta <- dplyr::filter(d_beta, significance == "B > A")
    if(nrow(d_beta) > 0){
      d_beta$significance <- "Exp"
      #Extract query
      d_beta <- dplyr::select(d_beta, c("significance", "sequence"))
      names(d_beta) <- c(str_c(time2, "dif"), "cdr3nt")
      
      #Integrate binomial result into join table
      d <- dplyr::left_join(d, d_beta, by = c("ntvj"="cdr3nt"))
    } else{
      names.old <- names(d)
      d$none <- NA
      names(d) <- c(names.old, str_c(time2, "dif"))
    }
  }
  
  #For data of ID18 (no TP11 data)
  if(ncol(d) < 27){
    d$TP11dif <- "UC"
  }
  
  ###Define expanded clone with more than 2 fold increase of their frequency
  d[is.na(d)] <- "UC"
  d$TP2dif[which(d$TP1*2 > d$TP2 | d$TP1 > 0.0001)] <- "UC"
  d$TP3dif[which(d$TP1*2 > d$TP3 | d$TP2*2 > d$TP3 | d$TP1 > 0.0001)] <- "UC"
  d$TP4dif[which(d$TP1*2 > d$TP4 | d$TP3*2 > d$TP4 | d$TP1 > 0.0001)] <- "UC"
  d$TP5dif[which(d$TP1*2 > d$TP5 | d$TP4*2 > d$TP5 | d$TP1 > 0.0001)] <- "UC"
  d$TP6dif[which(d$TP1*2 > d$TP6 | d$TP5*2 > d$TP6 | d$TP1 > 0.0001)] <- "UC"
  d$TP8dif[which(d$TP1*2 > d$TP8 | d$TP6*2 > d$TP8 | d$TP1 > 0.0001)] <- "UC"
  d$TP9dif[which(d$TP1*2 > d$TP9 | d$TP8*2 > d$TP9 | d$TP1 > 0.0001)] <- "UC"
  d$TP10dif[which(d$TP1*2 > d$TP10 | d$TP9*2 > d$TP10 | d$TP1 > 0.0001)] <- "UC"
  d$TP11dif[which(d$TP1*2 > d$TP11 | d$TP10*2 > d$TP11 | d$TP1 > 0.0001)] <- "UC"
  
  #Append changes in frequency after vaccination
  d$TP2difFreq <- d$TP2 - d$TP1
  d$TP3difFreq <- d$TP3 - d$TP2
  d$TP4difFreq <- d$TP4 - d$TP3
  d$TP5difFreq <- d$TP5 - d$TP4
  d$TP6difFreq <- d$TP6 - d$TP5
  d$TP8difFreq <- d$TP8 - d$TP6
  d$TP9difFreq <- d$TP9 - d$TP8
  d$TP10difFreq <- d$TP10 - d$TP9
  d$TP11difFreq <- d$TP11 - d$TP10
  
  ###Judge early/main/third responders
  d$Resp <- "Others"
  #Third responder: not detected until TP8, expanded at TP8
  d$Resp[which(d$TP8dif == "Exp" & d$TP1 == 0 & d$TP2 == 0 & d$TP3 == 0 &
                 d$TP4 == 0 & d$TP5 == 0 & d$TP6 == 0)] <- "Third"
  #Early responder: expanded at TP2
  d$Resp[which(d$TP2dif == "Exp")] <- "Early"
  #Main responder: not expanded at TP2, expanded at TP3
  d$Resp[which(d$TP2dif == "UC" & d$TP3dif == "Exp")] <- "Main"
  ##Trial
  #Late responder: not expanded at TP2, expanded at TP3
  d$Resp[which(d$TP2dif == "UC" & d$TP3dif == "UC" & d$TP4dif == "Exp")] <- "Late"
  
  #Output
  name.output <- str_c(dir.output, file.name, sep = "/") %>% str_replace("join.strict.table.txt", "DifAbund.csv")
  fwrite(d, name.output, sep=",", nThread=32)
}

############################### Processing ########################################
files <- list.files(dir.repertoire, pattern="strict.table")
  
dir.create(dir.output, recursive = TRUE)
  
  
t<-proc.time()
cl <- makeCluster(cores)
registerDoParallel(cl)   
output <- foreach(file.name = files, .combine = rbind,
                    .packages=c("data.table", "stringr", "dplyr")) %dopar% {Extract(file.name, dir.binomial, dir.output)}
stopCluster(cl)
proc.time()-t
