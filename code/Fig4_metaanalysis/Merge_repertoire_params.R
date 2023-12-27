##Merge repertoire parameters based on correlation heatmap

#load packages
library(ComplexHeatmap)
library(circlize)
library(stringr)
library(dplyr)
library(qvalue)
library(data.table)

#load Fonts
library(extrafont)
loadfonts("win")

#Input layer
Cells <- c("CD4", "CD8")

dir.input <- "tmp/metadata"
dir.output <- "tmp/result/Fig4"

################### Define functions ######################
tableread_fast = function(i, header=TRUE, sep="\t"){
  tmp = fread(i, header=header, sep=sep, nThread=32)
  tmp = as.data.frame(tmp)
  return(tmp)
}

createEmptyDf = function( nrow, ncol, colnames = c() ){
  data.frame( matrix( vector(), nrow, ncol, dimnames = list( c(), colnames ) ) )
}

#################### Processing ##########################
dir.create(dir.output, recursive = TRUE)

#Load parameter merge setting
name.input <- str_c(dir.input, "Metaanalysis_repertoire_merge.csv", sep = "/")
merge <- tableread_fast(name.input, header = TRUE, sep=",")

for(cell in Cells){
  #load data
  name.input <- str_c(dir.input, "Metaanalysis_", sep = "/") %>% str_c(cell, ".txt") 
  d <- read.table(name.input, header = TRUE)
  
  #Normalization
  #Use Spearman coefficient -> averaging the rank of parameters in each clusterる
  d_rank <- createEmptyDf(nrow(d), ncol(d), colnames = names(d))
  for(i in 1:ncol(d)){
    d_rank[,i] <- rank(d[,i])
  }
  
  ##Take Average rank for meta-parameters
  merge_sub <- dplyr::filter(merge, Cell == cell)
  #Prepare output
  d_out <- createEmptyDf(nrow(d), nrow(merge_sub), colnames = merge_sub[,2])
  #take average
  for(j in 1:nrow(merge_sub)){
    #Call meta-parameter name
    name_param <- merge_sub[j,2]
    #Call individual parameter names
    rep_params <- str_split(merge_sub[j,3], pattern = "-", simplify = TRUE)
    rep_params <- rep_params[rep_params > 0]
    #Average
    d_rank_sub <- dplyr::select(d_rank, rep_params)
    d_out[,j] <- apply(d_rank_sub, 1, mean)
  }
  
  #Output
  name.output <- str_c(dir.output, "Merged.Metaanalysis_", sep = "/") %>% str_c(cell, ".csv")
  write.csv(d_out, name.output, row.names = FALSE)
}

