#Prepare input file from Minervina et al. (2022, https://www.nature.com/articles/s41590-022-01184-4)
#Download from github (https://github.com/pogorely/COVID_vax_CD8), and format for analysis using tcrdist3
#Epitope.Information.csv was manually generated from Fig1D of Minervina et al. (2022)

#load libraryies
library(HelpersMG)
library(dplyr)
library(stringr)

#Input layer
dir.output <- "Minervina"

#Create directory
dir.create(dir.output)

#Download data
wget(c("https://raw.githubusercontent.com/pogorely/COVID_vax_CD8/master/cd8_only_dextr_rev_clean.tsv"))
d <- read.table("cd8_only_dextr_rev_clean.tsv", header = TRUE, sep = "\t")
Epitopes <- read.csv("Epitope.Information.csv", header = TRUE)

#Extract columns used for tcrdist3 analysis
d <- dplyr::select(d, c("cdr3b", "vb", "jb", "donor", "epitope"))
#remove cells containing NA
d <- na.omit(d)

#Separate data.table by targeted epitopes
list_epitopes <- unique(d$epitope)
for(i in list_epitopes){
  #Collect information about epitope 
  epitope.inf <- dplyr::filter(Epitopes, Epitope == i)
  #Only analyze epitopes in the list
  if(nrow(epitope.inf) > 0){
    peptide <- epitope.inf$Peptide
    
    #Extract cells matching epitope specificity
    d_sub <- dplyr::filter(d, epitope == i) %>% dplyr::select(-epitope)
    #change column names
    d_sub$cell_type <- "PBMC"
    names(d_sub) <- c("cdr3_b_aa", "v_b_gene", "j_b_gene", "subject", "cell_type")
    d_sub$v_b_gene <- str_c(d_sub$v_b_gene, "*01", sep = "")
    d_sub$j_b_gene <- str_c(d_sub$j_b_gene, "*01", sep = "")
    d_sub <- dplyr::select(d_sub, c("subject", "cell_type", "v_b_gene", "j_b_gene", "cdr3_b_aa"))
    
    #Output
    name_output <- str_c(dir.output, "MinervinaNI", sep = "/") %>%
      str_c(i, epitope.inf$Peptide, epitope.inf$Source, sep = "+") %>% str_c("tcrdist3.csv", sep = ".")
    write.csv(d_sub, name_output, row.names = FALSE)
  }
}
