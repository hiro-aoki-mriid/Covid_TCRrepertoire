#Prepare input file from Minervina et al. (2022, https://www.nature.com/articles/s41590-022-01184-4)
#Download from github (https://github.com/pogorely/COVID_vax_CD8), and format for analysis using tcrdist3
#Epitope.Information.csv was manually generated from Fig1D of Minervina et al. (2022)

#load libraryies
library(foreach)
library(doParallel)
library(HelpersMG)
library(dplyr)
library(stringr)
library(data.table)

#Input layer
dir.original <- "tmp/data/original/raw.data"
dir.input <- "tmp/data/public/Public.Clone.Input"
dir.output <- "tmp/result/intermediate/7_tcrdist/Input.tcrdist3"
###VDJdb
name1 <- "SearchTable-2023-01-25 05_53_41.437.tsv"
###Immune Epitope DataBase
name2 <- "tcell_receptor_table_export_1678716394.csv"
name.HLA <- "tmp/metadata/epitope_HLA_list.csv"
cores <- 12

######################################################################
#Create directory
dir.create(dir.output, recursive = TRUE)
#define function
tableread_fast = function(i, header=TRUE, quote="", sep="\t"){
  tmp = fread(i, header=header, sep=sep, quote=quote, nThread=32)
  tmp = as.data.frame(tmp)
  return(tmp)
}

createEmptyDf = function( nrow, ncol, colnames = c() ){
  data.frame( matrix( vector(), nrow, ncol, dimnames = list( c(), colnames ) ) )
}


###For vdjdb dataset
#Download and loaddata
#wget(c("https://raw.githubusercontent.com/pogorely/COVID_vax_CD8/master/cd8_only_dextr_rev_clean.tsv"))
d1 <- tableread_fast(str_c(dir.input, name1,sep = "/"), header = TRUE, quote="\"", sep="\t")

###For vdjdb input data
#Extract columns used for tcrdist3 analysis
d1_sub <- str_split(d1$Meta, pattern = "subject.id", simplify = TRUE)
d1_sub <- str_split(d1_sub[,2], pattern = ",", simplify = TRUE)
d1$subject <- str_c(d1$Reference, str_remove(d1_sub[,1], ":"), sep = "_")
d1$cell_type <- "PBMC"
d1_sub <- dplyr::select(d1, c("subject", "cell_type", "V", "J", "CDR3", "MHC A", "MHC B", "Epitope", "Epitope gene", "Epitope species")) %>%
  dplyr::filter(str_detect(V, "TRBV"))
names(d1_sub) <- c("subject", "cell_type", "v_b_gene", "j_b_gene", "cdr3_b_aa", "MHC.A", "MHC.B", "Epitope", "Epitope.gene", "Epitope.species")

##Select SARS-CoV-2 Spike epitopes
d1_sub <- dplyr::filter(d1_sub, Epitope.species == "SARS-CoV-2" & Epitope.gene == "Spike")

#Select epitopes: having more than 20 TCRs
epi.counts <- as.data.frame(table(d1_sub$Epitope))
epi.analyze <- as.vector(dplyr::filter(epi.counts, Freq >= 20)[,1])
d1_sub <- dplyr::filter(d1_sub, Epitope %in% epi.analyze)

#HLA name unification
d_HLA <- read.csv(name.HLA, header = TRUE)
d <- data.frame()
for(i in 1:nrow(d_HLA)){
  epi <- d_HLA$Epitope[i]
  d1_sub2 <- dplyr::filter(d1_sub, Epitope == epi)
  if(nrow(d1_sub2) > 0){
    d1_sub2$MHC <- d_HLA$HLA[i]
    d <- rbind(d, d1_sub2)
  }
}
d$Epitope <- str_c("vdjdb_", d$Epitope)

##Generate HLA-epitope combination
d$EpiHLA <- str_c(d$MHC, d$Epitope, d$Epitope.gene, d$Epitope.species, sep = "+")

#Separate data.table by targeted epitopes
list_epitopes <- unique(d$EpiHLA)
for(i in list_epitopes){
  #Extract cells matching epitope specificity
  d_sub <- dplyr::filter(d, EpiHLA == i) %>% dplyr::select(-c("EpiHLA", "MHC", "Epitope", "Epitope.gene", "Epitope.species", "MHC.A", "MHC.B"))
  
  name_output <- str_c(dir.output, i, sep = "/") %>% str_c("tcrdist3.csv", sep = ".")
  write.csv(d_sub, name_output, row.names = FALSE)
}

###For IEDB dataset
#Download and loaddata
d1 <- tableread_fast(str_c(dir.input, name2,sep = "/"), header = TRUE, quote="\"", sep=",")

##Select SARS-CoV-2 Spike epitopes
d1_sub <- dplyr::filter(d1, Organism == "SARS-CoV2" & str_detect(`MHC Allele Names`, ":") & `Chain 2 CDR3 Curated` != "")

#Select epitopes: having more than 20 TCRs
epi.counts <- as.data.frame(table(d1_sub$Description))
epi.analyze <- as.vector(dplyr::filter(epi.counts, Freq >= 20)[,1])
d1 <- dplyr::filter(d1_sub, Description %in% epi.analyze)

###For vdjdb input data
#Extract columns used for tcrdist3 analysis
d1$Antigen <- "Spike"
#Separate TCRs with multiple assays
d1_assay_count <- sum(str_split(d1$`Assay IDs`, pattern = ", ", simplify = TRUE) > 0)
d2 <- createEmptyDf(d1_assay_count, ncol(d1), colnames = names(d1) )
j <- 0
d1_assays <- vector()
for(i in 1:nrow(d1)){
  d1_sub <- d1[i,]
  d1_assay <- as.vector(str_split(d1_sub$`Assay IDs`, pattern = ", ", simplify = TRUE))
  len_d1_sub <- length(d1_assay)
  d2[(j+1):(j+ len_d1_sub),] <- d1_sub
  j <- j + len_d1_sub
  d1_assays <- c(d1_assays, d1_assay)
}
d2$Assay.IDs <- d1_assays

d2$cell_type <- "PBMC"
d1_sub <- dplyr::select(d2, c("Assay.IDs", "cell_type", "Curated.Chain.2.V.Gene", "Curated.Chain.2.J.Gene", "Chain.2.CDR3.Curated",
                              "MHC.Allele.Names", "Description", "Antigen", "Organism"))
names(d1_sub) <- c("subject", "cell_type", "v_b_gene", "j_b_gene", "cdr3_b_aa", "MHC", "Epitope", "Epitope.gene", "Epitope.species")

#HLA name unification
d_HLA <- read.csv(name.HLA, header = TRUE)
d <- data.frame()
for(i in 1:nrow(d_HLA)){
  epi <- d_HLA$Epitope[i]
  d1_sub2 <- dplyr::filter(d1_sub, Epitope == epi)
  if(nrow(d1_sub2) > 0){
    d1_sub2$MHC <- d_HLA$HLA[i]
    d <- rbind(d, d1_sub2)
  }
}
d$Epitope <- str_c("IEDB_", d$Epitope)

##Generate HLA-epitope combination
d$EpiHLA <- str_c(d$MHC, d$Epitope, d$Epitope.gene, d$Epitope.species, sep = "+")

#Separate data.table by targeted epitopes
list_epitopes <- unique(d$EpiHLA)
for(i in list_epitopes){
  #Extract cells matching epitope specificity
  d_sub <- dplyr::filter(d, EpiHLA == i) %>% dplyr::select(-c("EpiHLA", "MHC", "Epitope", "Epitope.gene", "Epitope.species"))
  d_sub$v_b_gene <- str_c(d_sub$v_b_gene, "*01")
  d_sub$j_b_gene <- str_c(d_sub$j_b_gene, "*01")
  
  name_output <- str_c(dir.output, i, sep = "/") %>% str_c("tcrdist3.csv", sep = ".")
  write.csv(d_sub, name_output, row.names = FALSE)
}

###For our original dataset
#Main function
Main <- function(file.name, dir.original, dir.output){
  ##load data
  name.input <- str_c(dir.original, file.name, sep = "/")
  d <- tableread_fast(name.input, header = TRUE, quote="\"", sep="\t")
  names <- str_split(file.name, pattern = "_", simplify = TRUE)
  d$sample <- names[5]
  d$cell_type <- "PBMC"
  
  ##change data format
  d_sub <- select(d, c("sample", "cell_type", "v", "j", "cdr3aa"))
  names(d_sub) <- c("subject", "cell_type", "v_b_gene", "j_b_gene", "cdr3_b_aa")
  d_sub$v_b_gene <- str_c(d_sub$v_b_gene, "*01")
  d_sub$j_b_gene <- str_c(d_sub$j_b_gene, "*01")
  
  return(d_sub)
}

epitopes <- c("Tet_NF9", "Tet_QI9")

for(epitope in epitopes){
  #For NF9 sepcific clones
  files  <- list.files(dir.original, pattern=epitope)
  cl <- makeCluster(cores)
  registerDoParallel(cl)   
  d <- foreach(file.name = files, .combine = rbind,
               .packages=c("ggplot2", "extrafont", "stringr", "dplyr", "data.table")) %dopar% {Main(file.name, dir.original, dir.output)}
  stopCluster(cl)
  d$MHC <- "A24-02"
  d$Epitope <- epitope
  d$Epitope.gene <- "Spike"
  d$Epitope.species <- "SARS-CoV-2"
  
  EpiHLA <- str_c("A24-02", epitope, "Spike", "SARS-CoV-2", sep = "+")
  
  #Extract cells matching epitope specificity
  d_sub <- dplyr::select(d, -c("MHC", "Epitope", "Epitope.gene", "Epitope.species"))
  
  #Output
  name_output <- str_c(dir.output, EpiHLA, sep = "/") %>% str_c("tcrdist3.csv", sep = ".")
  write.csv(d_sub, name_output, row.names = FALSE)
}



