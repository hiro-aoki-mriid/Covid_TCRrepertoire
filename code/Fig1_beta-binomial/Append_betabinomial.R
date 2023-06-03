##Merge different abundance clone into tcr repertoire dataset

library(foreach)
library(doParallel)
library(stringr)
library(data.table)
library(dplyr)

dir.binomial <- "tmp/result/intermediate/1_beta-binomial/binomial_test"
dir.repertoire <- "tmp/result/intermediate/1_beta-binomial/JoinTP"
dir.output <- "tmp/result/intermediate/1_beta-binomial/JoinTP_DifAbund"

#Output params

FirstSecondType <- c("1st", "2nd", "Dual", "Others", "UD")


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
Extract <- function(file.name, dir.binomial, dir.output, FirstSecondType, ThirdType, condition){
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
      names(d_beta) <- c(str_c(time1, time2, "dif"), "cdr3nt")
      
      #Integrate binomial result into join table
      d <- dplyr::left_join(d, d_beta, by = c("ntvj"="cdr3nt"))
    } else{
      names.old <- names(d)
      d$none <- NA
      names(d) <- c(names.old, str_c(time1, time2, "dif"))
    }
  }
  
  ###Define expanded clone with more than 2 fold increase of their frequency
  d$TP1TP3dif[which(d$TP1TP3dif == "Exp" & d$TP1*2 > d$TP3)] <- "UC"
  d$TP1TP2dif[which(d$TP1TP2dif == "Exp" & d$TP1*2 > d$TP2)] <- "UC"
  d$TP2TP3dif[which(d$TP2TP3dif == "Exp" & d$TP2*2 > d$TP3)] <- "UC"
  d$TP6TP8dif[which(d$TP6TP8dif == "Exp" & d$TP6*2 > d$TP8)] <- "UC"
  ###Define expanded clone with more than 2 fold increase of their frequency at P1
  d$TP2TP3dif[which(d$TP2TP3dif == "Exp" & d$TP1*2 > d$TP3)] <- "UC"
  d$TP6TP8dif[which(d$TP6TP8dif == "Exp" & d$TP1*2 > d$TP8)] <- "UC"
  
  #Add increase/decrease between time points
  d$TP1TP2dif[is.na(d$TP1TP2dif)] <- "UC"
  d$TP1TP3dif[is.na(d$TP1TP2dif)] <- "UC"
  d$TP1TP2dif[which(d$TP1TP2dif == "UC" & d$TP1*2 < d$TP2 & d$TP2 > 0.000025)] <- "Inc"
  d$TP2TP3dif[is.na(d$TP2TP3dif)] <- "UC"
  d$TP2TP3dif[which(d$TP2TP3dif == "UC" & d$TP2 < d$TP3 & d$TP1*2 < d$TP3 & d$TP3 > 0.000025)] <- "Inc"
  d$TP6TP8dif[is.na(d$TP6TP8dif)] <- "UC"
  d$TP6TP8dif[which(d$TP6TP8dif == "UC" & d$TP6 < d$TP8 & d$TP1*2 < d$TP8 & d$TP8 > 0.000025)] <- "Inc"
  
  ###Judge first/second/dual/third responders
  d$Resp <- "Others"
  #3rd responder: not detected from TP1 to TP6, and expanded at TP8
  d$Resp[which(d$TP1 == 0 & d$TP2 == 0 & d$TP3 == 0 & d$TP4 == 0 & d$TP5 == 0 & d$TP6 == 0 &  d$TP8 != 0 &
                 d$TP6TP8dif == "Exp")] <- "3rd"
  #3rd responder: not detected from TP1 to TP6, and increased at TP8
  d$Resp[which(d$TP1 == 0 & d$TP2 == 0 & d$TP3 == 0 & d$TP4 == 0 & d$TP5 == 0 & d$TP6 == 0 &  d$TP8 != 0 &
                 d$TP6TP8dif == "Inc")] <- "3rdInc"
  ##1st responder: Expand at TP2 -> not increase at TP3
  d$Resp[which(d$TP1TP2dif == "Exp" & d$TP2TP3dif == "UC" & d$TP1 < 0.0001)] <- "1st"
  ##1st increased: Increased at TP2 -> not increase at TP3
  d$Resp[which(d$TP1TP2dif == "Inc" & d$TP2TP3dif == "UC" & d$TP1 < 0.0001)] <- "1stInc"
  ##2nd responder: not increase at TP2 -> Expand at TP3
  d$Resp[which(d$TP2TP3dif == "Exp" & d$TP1TP2dif == "UC" & d$TP1 < 0.0001)] <- "2nd"
  ##2nd increased: not increase at TP2 -> Increased at TP3
  d$Resp[which(d$TP2TP3dif == "Inc" & d$TP1TP2dif == "UC" & d$TP1 < 0.0001)] <- "2ndInc"
  ##Dual Increased: Inc at TP2 & TP3 
  d$Resp[which(d$TP1TP2dif != "UC" & d$TP2TP3dif != "UC" & d$TP1 < 0.0001)] <- "DualInc"
  ##Dual responder: Exp from TP1 to TP3 -> Inc at TP2 & TP3 
  d$Resp[which(d$TP1TP3dif == "Exp" & d$Resp == "DualInc" & d$TP1 < 0.0001)] <- "Dual"

  #Append changes in frequency after vaccination
  d$TP1_TP2 <- d$TP2 - d$TP1
  d$TP2_TP3 <- d$TP3 - d$TP2
  d$TP6_TP8 <- d$TP8 - d$TP6
  
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
                    .packages=c("data.table", "stringr", "dplyr")) %dopar% {Extract(file.name, dir.binomial, dir.output,
                                                                                    FirstSecondType, ThirdType, condition)}
stopCluster(cl)
proc.time()-t
