###Merge metaclonotype information from tcrdist3 to result of difabund and AIM analysis

#load libraryies

library(foreach)
library(doParallel)
library(data.table)
library(stringr)
library(dplyr)
library(ggplot2)

#Input layer
dir.input <- "tmp/result/intermediate/5_Tet/JoinTP_DifAbund_Tet"
dir.input.tcrdist <- "tmp/result/intermediate/7_tcrdist/tcrdist3_output"
dir.output <- "tmp/result/intermediate/7_tcrdist/JoinTP_DifAbund_Tet_tcrdist"
name.HLA <- "tmp/metadata/epitope_HLA_list.csv"
cores <- 12
metaclonotypes <- c("vdjdb", "iedb", "KmayerB", "Tet")
query_params <- c("protein_coordinate_search", "protein_search", "pathogen_search", "feature_search", "HLA_search")
epitopes <- c("S269", "S448", "S919", "S1208")

############################### Functions ########################################
###Define functions
#Table.read.fast
tableread_fast = function(i, header=TRUE, quote="", sep="\t"){
  tmp = fread(i, header=header, sep=sep, quote=quote, nThread=32)
  tmp = as.data.frame(tmp)
  return(tmp)
}
createEmptyDf = function( nrow, ncol, colnames = c() ){
  data.frame( matrix( vector(), nrow, ncol, dimnames = list( c(), colnames ) ) )
}

#Main module
tcrdist_merge = function(file.name, dir.input, dir.input.tcrdist, dir.output, metaclonotypes){
  #load data
  name_difabund <- str_c(dir.input, file.name, sep = "/")
  d_difabund <- tableread_fast(name_difabund, header = TRUE, quote="\"", sep=",")
  for(metaclonotype in metaclonotypes){
    name_tcrdist <- str_c(dir.input.tcrdist, metaclonotype, sep = "/") %>%
      str_c(file.name, sep = "") %>%
      str_replace("DifAbund.tetramer.csv", "metaclonotype.csv")
    
    if(file.exists(name_tcrdist)){
      d_tcrdist3 <- tableread_fast(name_tcrdist, header = TRUE, quote="\"", sep=",")
      
      ##Extract Metaclonotype-matched clones from difabund dataset
      #Define clonotype by CDR3aa, V usage, and J usage
      d_tcrdist3_q <- dplyr::select(d_tcrdist3, query_params)
      #Change names of source protein
      if(metaclonotype == "KmayerB"){
        d_tcrdist3_q$protein_search <- str_replace(d_tcrdist3_q$protein_search, "Spike", "S")
        d_tcrdist3_q$protein_coordinate_search <- str_c(lapply(d_tcrdist3_q$protein_coordinate_search,
                                                         function(x) {as.vector(d_HLA$Name[match(x, d_HLA$Epitope)])}))
      }
      if(metaclonotype != "KmayerB"){
        d_tcrdist3_q$protein_coordinate_search <- str_split(d_tcrdist3_q$protein_coordinate_search, pattern = "_", simplify = TRUE)[,2]
      }
      
      names(d_tcrdist3_q) <- str_c(names(d_tcrdist3_q), metaclonotype, sep = "_")
      d_tcrdist3$v <- str_sub(d_tcrdist3$v_b_gene_bulk, end = -4)
      d_tcrdist3$j <- str_sub(d_tcrdist3$j_b_gene_bulk, end = -4)
      d_tcrdist3_q$query <- str_c(d_tcrdist3$cdr3_b_nucseq_bulk, d_tcrdist3$v, d_tcrdist3$j, sep = "_")
      
    } else{
      d_tcrdist3_q <- createEmptyDf(1, (length(query_params)+1), colnames = c(str_c(query_params, metaclonotype, sep = "_"), "query"))
    }
    #Merge tcrdist3 information
    d_difabund <- left_join(d_difabund, d_tcrdist3_q, by = c("ntvj" = "query"))
    d_difabund[is.na(d_difabund)] <- "non"
    }
  
  #Define specificity of clones
  target.finder <- function(x, y){
    i <- 0
    out <- "Nega"
    while(i < length(x)){
      i <- i+1
      if(x[i] == y){
        out <- "Posi"
        break
      }
    }
    return(out)
  }
  d_specificity <- createEmptyDf(nrow(d_difabund), length(epitopes), colnames = str_c(epitopes, "motif", sep = "_"))
  for(i in 1:length(epitopes)){
    d_specificity[,i] <- apply(dplyr::select(d_difabund, contains("protein_coordinate_search")),
                               1, target.finder, y=epitopes[i])
  }
  d_difabund <- cbind(d_difabund, d_specificity)
  
  #Define clones detected in metaclonotypes
  tet.finder <- function(x){
    i <- 0
    out <- "Nega"
    while(i < length(x)){
      i <- i+1
      if(x[i] != "Nega"){
        out <- "Posi"
        break
      }
    }
    return(out)
  }
  d_difabund$motif <- apply(d_specificity, 1, tet.finder)
  
  #Output
  name.output <- str_c(dir.output, file.name, sep = "/") %>% str_replace("csv", "tcrdist3.csv")
  write.table(d_difabund, name.output, row.names = FALSE, quote = FALSE, sep = ",")
}

############################### Processing ########################################
dir.create(dir.output)
files <- list.files(dir.input, pattern="NaraCOVID")
d_HLA <- read.csv(name.HLA, header = TRUE)

##Combine difabund and tcrdist3 informations
t<-proc.time()
cl <- makeCluster(cores)
registerDoParallel(cl)   
foreach(file.name = files,
        .packages=c("data.table", "stringr", "dplyr", "ggplot2"),
        .combine = rbind) %dopar% {tcrdist_merge(file.name, dir.input, dir.input.tcrdist, dir.output, metaclonotypes)}
stopCluster(cl)
proc.time()-t

