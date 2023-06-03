##Convert to Immunoseq format

library(foreach)
library(doParallel)
library(stringr)
library(data.table)
library(dplyr)

dir.input <- "tmp/data/original/raw.data"
dir.output <- "tmp/result/intermediate/1_beta-binomial/Immunoseq"

files  <- list.files(dir.input, pattern="strict.table.txt")
cores <- 24

tableread_fast = function(i, header=TRUE, quote="", sep="\t"){
  tmp = fread(i, header=header, sep=sep, quote=quote, nThread=32)
  tmp = as.data.frame(tmp)
  return(tmp)
}
createEmptyDf = function( nrow, ncol, colnames = c() ){
  data.frame( matrix( vector(), nrow, ncol, dimnames = list( c(), colnames ) ) )
}

#====================================================================
ConvImmu <- function(file.name){ 
  name.input <- str_c(dir.input, file.name, sep = "/")
  d <- tableread_fast(name.input, header = TRUE, quote="", sep="\t")
  
  #Select raws used for Immunoseq format
  d_sub <- dplyr::select(d, c("cdr3nt", "cdr3aa", "count", "freq"))
  d_sub$cdr3Length <- nchar(d_sub$cdr3nt)
  names(d_sub) <- c("nucleotide", "aminoAcid", "count", "frequencyCount")
  d_sub$frequencyCount <- 100*d_sub$frequencyCount
  #Adjust total read count to 100K (with zero padding)
  d_sub$count <- round(100000*d_sub$count/sum(d_sub$count), digits = 0) + 1
  
  #Add v/j information to nucleotide column
  d_sub$nucleotide <- str_c(d_sub$nucleotide, d$v, d$j, sep = "_")
  
  #Additional columns for Immunoseq
  add_column1 <- c("vFamilyName", "vGeneName", "vGeneAllele", "vFamilyTies", "vGeneNameTies", "vGeneAlleleTies",
                   "dFamilyName", "dGeneName", "dGeneAllele", "dFamilyTies", "dGeneNameTies", "dGeneAlleleTies",
                   "jFamilyName", "jGeneName", "jGeneAllele", "jFamilyTies", "jGeneNameTies", "jGeneAlleleTies",
                   "vDeletion", "d5Deletion", "d3Deletion", "jDeletion", "n2Insertion", "n1Insertion", "vIndex",
                   "n2Index", "dIndex", "n1Index", "jIndex", "vdNormalizationFactor", "jNormalizationFactor")
  d_add1 <- createEmptyDf( nrow(d_sub), length(add_column1), colnames = add_column1 )
  d_add1[,] <- "NA"
  d_sub <- cbind(d_sub, d_add1)
  d_sub$inputTemplateEstimate <- d_sub$count
  d_sub$frequencyInputTemplateEstimate <- d_sub$frequencyCount
  d_sub$sequenceStatus <- "In"
  add_column2 <- c("vMaxResolved", "d2MaxResolved", "dMaxResolved", "jMaxResolved", "cloneResolved", "d2FamilyName",
                   "d2GeneName", "d2GeneAllele", "d2FamilyTies", "d2GeneNameTies", "d2GeneAlleleTies", "vScore",
                   "dScore", "jScore", "vAlignLength", "vAlignSubstitutionCount", "dAlignLength", "dAlignSubstitutionCount",
                   "jAlignLength", "jAlignSubstitutionCount", "vOrphon", "dOrphon", "jOrphon", "vFunction", "dFunction",
                   "jFunction", "vAlignSubstitutionIndexes", "dAlignSubstitutionIndexes", "jAlignSubstitutionIndexes",
                   "vAlignSubstitutionGeneThreePrimeIndexes", "dAlignSubstitutionGeneThreePrimeIndexes",
                   "jAlignSubstitutionGeneThreePrimeIndexes")
  d_add2 <- createEmptyDf( nrow(d_sub), length(add_column2), colnames = add_column2 )
  d_add2[,] <- "NA"
  d_sub <- cbind(d_sub, d_add2)
  
  ##Output data
  name.output <- str_c(dir.output, file.name, sep = "/") %>% str_replace("pool.strict.table.txt", "tsv")
  fwrite(d_sub, file = name.output, row.names = FALSE, sep="\t", nThread=32)
}

#===============================================================
dir.create(dir.output, recursive = TRUE)
t<-proc.time()
cl <- makeCluster(cores)
registerDoParallel(cl)   
foreach(file.name = files, .combine = rbind,
        .packages=c("data.table", "stringr", "dplyr")) %dopar% {ConvImmu(file.name)}
stopCluster(cl)
proc.time()-t


