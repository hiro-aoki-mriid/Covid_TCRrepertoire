#Prepare input file from Francis et al. (2021, https://www.science.org/doi/10.1126/sciimmunol.abk3070)
#Supplementary data file S3 was downloaded, then converted to csv file 

#load libraryies
library(dplyr)
library(stringr)

#Input layer
dir.output <- "Francis"

#Create directory
dir.create(dir.output)

#Download data
d <- read.table("Francis.SciImmunol.2021.DataS3.csv", header = TRUE, sep = ",")

#Extract columns used for tcrdist3 analysis
CDR3 <- str_remove_all(d$clonotype, pattern = ";") %>% str_split(pattern = "_", simplify = TRUE)
d$cdr3b <- CDR3[,2]

#Remove "*" from HLA name
HLA_alleles <- c("A", "B")
d_out <- data.frame()
for(allele in HLA_alleles){
  d_sub <- dplyr::filter(d, str_detect(HLA, allele))
  d_sub$HLA <- str_sub(d_sub$HLA, start=3)
  d_sub$HLA <- str_c(allele, d_sub$HLA, sep = "")
  d_out <- rbind(d_out, d_sub)
}
d_out$HLA <- str_remove(d_out$HLA, ":")

#Collect data needed for analysis
d$epitope <- str_c(d_out$HLA, d_out$epitope, d_out$antigen, sep = "+")
d <- dplyr::filter(d, cdr3b != "No-beta-detected") %>% dplyr::filter(str_detect(epitope, "SARS2")) %>%
  dplyr::select(c("cdr3b", "beta1_vgene", "beta1_jgene", "patient_ID", "epitope", "nCells"))

#Separate data.table by targeted epitopes
list_epitopes <- unique(d$epitope)
for(i in list_epitopes){
  #Extract cells matching epitope specificity
  d_sub <- dplyr::filter(d, epitope == i) %>% dplyr::select(-epitope)
  #change column names
  d_sub$cell_type <- "PBMC"
  names(d_sub) <- c("cdr3_b_aa", "v_b_gene", "j_b_gene", "subject", "count", "cell_type")
  d_sub$v_b_gene <- str_c(d_sub$v_b_gene, "*01", sep = "")
  d_sub$j_b_gene <- str_c(d_sub$j_b_gene, "*01", sep = "")
  d_sub <- dplyr::select(d_sub, c("subject", "cell_type", "v_b_gene", "j_b_gene", "cdr3_b_aa", "count"))
    
  #Output (only epitope with more than 15 different TCRs)
  if(nrow(d_sub) > 14){
    name_output <- str_c(dir.output, "FrancisSciImmunol", sep = "/") %>%
      str_c(i, sep = "+") %>% str_c("tcrdist3.csv", sep = ".")
    write.csv(d_sub, name_output, row.names = FALSE)
  }
}
