###Merge metaclonotype information from tcrdist3 to result of WGCNA and AIM analysis

#load libraryies

library(foreach)
library(doParallel)
library(data.table)
library(stringr)
library(dplyr)
library(ggplot2)

#Input layer
dir.input <- "Fig1_Fig2/WGCNA_output_150/Correlation.Table"
dir.input.tcrdist <- "Fig3/tcrdist3_output"
dir.output <- "Fig3/WGCNA_tcrdist"
cores <- 12
metaclonotypes <- c("minervina", "kmayerb", "francis")

#################################### Processing layer ####################################################
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
  name_wgcna <- str_c(dir.input, file.name, sep = "/")
  d_wgcna <- tableread_fast(name_wgcna, header = TRUE, quote="\"", sep=",")
  d_wgcna$query <- str_c(d_wgcna$cdr3aa, d_wgcna$v, d_wgcna$j, sep = "_")
  for(metaclonotype in metaclonotypes){
    name_tcrdist <- str_c(dir.input.tcrdist, file.name, sep = "/") %>% str_replace("wgcna.csv", metaclonotype) %>% str_c("metaclonotype.csv", sep = ".")
    d_tcrdist3 <- tableread_fast(name_tcrdist, header = TRUE, quote="\"", sep=",")
    
    ##Extract Metaclonotype-matched clones from WGCNA dataset
    #Define clonotype by CDR3aa, V usage, and J usage
    d_tcrdist3$v_b_gene_bulk <-str_sub(d_tcrdist3$v_b_gene_bulk, end=-4)
    d_tcrdist3$j_b_gene_bulk <-str_sub(d_tcrdist3$j_b_gene_bulk, end=-4)
    d_tcrdist3_q <- dplyr::select(d_tcrdist3, c("protein_coordinate_search", "protein_search", "feature_search"))
    #Change names of source protein
    if(metaclonotype == "minervina"){
      d_tcrdist3_q$protein_search <- str_replace(d_tcrdist3_q$protein_search, "NSP3", "ORF1ab") %>%
        str_replace("ORF3", "ORF3a") %>% 
        str_replace("RNP", "N") %>% 
        str_replace("Spike", "S")
    }
    if(metaclonotype == "francis"){
      d_tcrdist3_q$protein_search <- str_replace_all(d_tcrdist3_q$protein_search, "3A", "ORF3a") %>%
        str_replace_all("SPIKE", "S")
    }
    names(d_tcrdist3_q) <- str_c(names(d_tcrdist3_q), metaclonotype, sep = "_")
    d_tcrdist3_q$query <- str_c(d_tcrdist3$cdr3_b_aa_bulk, d_tcrdist3$v_b_gene_bulk, d_tcrdist3$j_b_gene_bulk, sep = "_")
    
    #Merge tcrdist3 information
    d_wgcna <- left_join(d_wgcna, d_tcrdist3_q, by = "query")
    d_wgcna[is.na(d_wgcna)] <- "non"
  }
  
  #Annotate protein source of peptide
  #Trust the results of pMHC result (minervina, francis) than MIRA assay (kmayerb)
  Annotate <- function(protein_search_kmayerb, protein_search_minervina, protein_search_francis){
    if(protein_search_minervina == "non" & protein_search_francis == "non"){
      out <- protein_search_kmayerb
    } else {
      if(protein_search_minervina != "non" & protein_search_francis == "non"){
        out <- protein_search_minervina
      } else {
        if(protein_search_minervina == "non" & protein_search_francis != "non"){
          out <- protein_search_francis
        } else{
          out <- "non"
        }
      }
    }
    return(out)
  }
  protein <- vector()
  for(i in 1:nrow(d_wgcna)){
    out <- Annotate(d_wgcna$protein_search_kmayerb[i], d_wgcna$protein_search_minervina[i], d_wgcna$protein_search_francis[i])
    protein <- c(protein, out)
  }
  d_wgcna$protein <- protein

  #Output
  name.output <- str_c(dir.output, file.name, sep = "/") %>% str_replace("wgcna.csv", "wgcna.tcrdist3.tsv")
  write.table(d_wgcna, name.output, row.names = FALSE, quote = FALSE, sep = "\t")
}

###Main module
setwd("../")
dir.create(dir.output)
files  <- list.files(dir.input, pattern="CD8")

##Combine WGCNA and AIM informations
t<-proc.time()
cl <- makeCluster(cores)
registerDoParallel(cl)   
foreach(file.name = files,
        .packages=c("data.table", "stringr", "dplyr", "ggplot2"),
        .combine = rbind) %dopar% {tcrdist_merge(file.name, dir.input, dir.input.tcrdist, dir.output, metaclonotypes)}
stopCluster(cl)
proc.time()-t

setwd("Fig3")
