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
dir.input <- "tmp/result/intermediate/7_tcrdist/JoinTP_DifAbund_AIM_Tet_tcrdist"
name.output <- "tmp/result/Fig7/tcrdist_clones.csv"
cores <- 12

#################################### Processing layer ####################################################
###Define functions
#Table.read.fast
tableread_fast = function(i, header=TRUE, quote="", sep=","){
  tmp = fread(i, header=header, sep=sep, quote=quote, nThread=32)
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
  d <- tableread_fast(name.input, header = TRUE, quote="\"", sep=",")
  
  d_sub <- dplyr::filter(d, class == "AIM" | tetramer != "Nega" |
                           protein_coordinate_search_vdjdb != "non" |
                           protein_coordinate_search_iedb != "non" |
                           protein_coordinate_search_KmayerB != "non" |
                           protein_coordinate_search_Tet != "non")
  if(nrow(d_sub) > 0){
    d_sub$sample <- file.name
  }
  return(d_sub)
}

###Main module
files  <- list.files(dir.input, pattern="NaraCOVID")
t<-proc.time()
cl <- makeCluster(cores)
registerDoParallel(cl)   
d_difabund <- foreach(file.name = files, .combine = rbind,
               .packages=c("ggplot2", "extrafont", "stringr", "dplyr", "data.table")) %dopar% {Main(file.name, dir.input, dir.output)}
stopCluster(cl)
proc.time()-t
d_difabund <- as.data.frame(d_difabund)

#Annotate protein source of peptide
Annotate <- function(vector){
  vec.len <- length(vector)
  i <- 1
  while(i <= vec.len){
    out <- vector[i]
    i <- i+1
    if(out != "non&non"){
      break
    }
  }
  return(out)
}
#Concatenate epitope, protein, pathogen and HLA allele
target_KmayerB <- str_c(d_difabund$protein_coordinate_search_KmayerB,
                        d_difabund$HLA_search_KmayerB,
                        sep = "&")
target_vdjdb <- str_c(str_remove(d_difabund$protein_coordinate_search_vdjdb, "vdjdb_"),
                      d_difabund$HLA_search_vdjdb,
                      sep = "&")
target_iedb <- str_c(str_remove(d_difabund$protein_coordinate_search_iedb, "IEDB_"),
                      d_difabund$HLA_search_iedb,
                      sep = "&")
target_Tet <- str_c(d_difabund$protein_coordinate_search_Tet,
                    d_difabund$HLA_search_Tet,
                    sep = "&")

#Annotate
target <- vector()
for(i in 1:nrow(d_difabund)){
  out <- Annotate(c(target_vdjdb[i], target_iedb[i], target_KmayerB[i], target_Tet[i]))
  target <- c(target, out)
}
d_difabund$target <- target

write.csv(d_difabund, name.output, row.names = FALSE)
