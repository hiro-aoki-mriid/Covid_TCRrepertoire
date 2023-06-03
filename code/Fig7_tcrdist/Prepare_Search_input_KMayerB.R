#Prepare input file for searching clones that match metaclonotype

#load libraryies
library(dplyr)
library(stringr)
library(foreach)
library(doParallel)
library(data.table)

#Input layer
name1 <- "KMayerB_Mira_HLAassociated.csv"
name2 <- "KMayerB_Mira_CoV2Protein.csv"

dir.input <- "tmp/data/public/Public.Clone.Input"
dir.output <- "tmp/result/intermediate/7_tcrdist/Metaclonotype_query"
cores <- 12
name.HLA <- "tmp/metadata/epitope_HLA_list.csv"

#################################### Processing layer ####################################################
dir.create(dir.output)
tableread_fast = function(i, header=TRUE, quote="", sep="\t"){
  tmp = fread(i, header=header, sep=sep, quote=quote, nThread=32)
  tmp = as.data.frame(tmp)
  return(tmp)
}

#Download and loaddata
d_main <- tableread_fast(str_c(dir.input, name1,sep = "/"), header = TRUE, quote="\"", sep=",")
d_protein <- tableread_fast(str_c(dir.input, name2,sep = "/"), header = TRUE, quote="\"", sep=",")

#Append target protein to main data
d1_sub <- left_join(d_main, d_protein, by= c("set" = "Mira"))
d1_sub$pathogen <- "SARS-CoV-2"
d1_sub <- dplyr::select(d1_sub, c("set", "protein", "pathogen", "hla_predicted", "protein_coordinate",
                        "feature", "cdr3_b_aa", "v_b_gene", "j_b_gene", "pgen", "radius", "regex"))
#Epitope???̕ϊ?
d1_sub$protein_coordinate <- str_replace(d1_sub$protein_coordinate, "YLQPRTFL_YLQPRTFLL_YYVGYLQPRTF", "YLQPRTFLL")
d1_sub$protein_coordinate <- str_replace(d1_sub$protein_coordinate, "NYLYRLFRK_NYNYLYRLF", "NYNYLYRLF")

#HLA name unification
d_HLA <- read.csv(name.HLA, header = TRUE)
d_sub <- data.frame()
for(i in 1:nrow(d_HLA)){
  epi <- d_HLA$Epitope[i]
  d1_sub2 <- dplyr::filter(d1_sub, protein_coordinate == epi)
  if(nrow(d1_sub2) > 0){
    d1_sub2$hla_predicted <- d_HLA$HLA[i]
    d_sub <- rbind(d_sub, d1_sub2)
  }
}

names(d_sub) <- c("set", "protein", "pathogen", "HLA", "protein_coordinate",
              "feature", "cdr3_b_aa", "v_b_gene", "j_b_gene", "pgen", "radius", "regex")
name.output <- str_c(dir.output, "KmayerB_metaclonotypes.tsv", sep = "/")
write.table(d_sub, name.output, sep = "\t", row.names = FALSE, quote = FALSE)



