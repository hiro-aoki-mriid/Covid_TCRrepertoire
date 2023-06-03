###Calculate frequency of SARS-Cov-2 reactive clones determined by tcrdist3 at each time point

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
name.input <- "tmp/result/Fig7/tcrdist_clones.csv"
name.hla <- "tmp/metadata/HLAtyping_result.csv"
name.output <- "tmp/result/Fig7/tcrdist_epitope_frequency.csv"
cores <- 12
params <- c("TP1", "TP2", "TP3", "TP4", "TP5", "TP6", "TP8", "TP9", "TP10", "TP11")

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

##load data
#load tcrdist3 data
d <- tableread_fast(name.input, header = TRUE, quote="\"", sep=",")
d <- dplyr::filter(d, target != "non&non")
#load hla typing data
hla_table <- read.csv(name.hla, header = TRUE, colClasses=c("character", "character"))

#Call array of epitopes
epitope_query <- unique(d$target)

output.all <- data.frame()
for(epitope in epitope_query){
  #Extract table of each epitope
  d_sub <- dplyr::filter(d, target == epitope)
  
  #Extract target participant
  HLA <- str_split(epitope, pattern = "&", simplify = TRUE)[2]
  hla_table_query <- dplyr::filter(hla_table, str_detect(HLAI, HLA))
  
  #Prepare output data frame
  output <- createEmptyDf(nrow(hla_table_query), length(params), colnames = params )
  
  #Check tcrdist hit for each participant
  for(i in 1:nrow(hla_table_query)){
    id <- hla_table_query$ID[i]
    d_sub2 <- dplyr::filter(d_sub, str_detect(sample, id))
    if(nrow(d_sub > 0)){
      output[i,] <- dplyr::select(d_sub2, params) %>% apply(2, sum)
    } else {
      output[i,] <- 0
    } 
  }
  output$ID <- hla_table_query$ID
  output$epitope <- epitope

  output.all <- rbind(output.all, output)
}

output.all <- left_join(output.all, hla_table, by = c("ID" = "ID"))

write.csv(output.all, name.output, row.names = FALSE)
