#Prepare input file for searching clones that match metaclonotype

#load libraryies
library(dplyr)
library(stringr)
library(foreach)
library(doParallel)
library(data.table)

#Input layer
dir.input <- "Minervina_hla_restricted_meta_clonotypes"
dir.output <- "Metaclonotype_query"
cores <- 12

#################################### Processing layer ####################################################
###Define functions
#Table.read.fast
tableread_fast = function(i, header=TRUE, quote="", sep="\t"){
  tmp = fread(i, header=header, sep=sep, quote=quote, nThread=32)
  tmp = as.data.frame(tmp)
  return(tmp)
}

#Main function
Main <- function(file.name, dir.input, dir.output){
  ##load data
  name.input <- str_c(dir.input, file.name, sep = "/")
  table <- tableread_fast(name.input, header = TRUE, quote="\"", sep="\t")
  
  if(nrow(table) > 0){
    #Collect epitope information
    Epitope.inf <- str_split(file.name, pattern = "\\.", simplify = TRUE)
    Epitope.inf <- str_split(Epitope.inf[1], pattern = "\\+", simplify = TRUE)
    protein <- Epitope.inf[4]
    protein_coodinate <- Epitope.inf[3]
    
    #Collect information needed for metaclonotype search
    out <- dplyr::select(table, c("cdr3_b_aa", "v_b_gene", "j_b_gene", "pgen", "radius", "regex"))
    out$set <- "1"
    out$sars_cov_2_genomic_position <- protein
    out$protein_coordinate <- protein_coodinate
    out$protein <- protein
    out$feature <- str_c("1E6", out$v_b_gene, out$cdr3_b_aa, out$radius, out$regex, sep = "+")
    out <- dplyr::select(out, c("set", "sars_cov_2_genomic_position", "protein_coordinate", "protein",
                                "feature", "cdr3_b_aa", "v_b_gene", "j_b_gene", "pgen", "radius", "regex"))
    
    return(out)
  }
}

###Main module
files  <- list.files(dir.input, pattern="ranked_centers_bkgd_ctlr_1E6.tsv")
dir.create(dir.output)
t<-proc.time()
cl <- makeCluster(cores)
registerDoParallel(cl)   
out <- foreach(file.name = files,
               .combine = rbind, .packages=c("ggplot2", "extrafont", "stringr", "dplyr", "data.table")) %dopar% {Main(file.name, dir.input, dir.output)}
stopCluster(cl)
proc.time()-t

out <- as.data.frame(out)
names(out) <- c("set", "sars_cov_2_genomic_position", "protein_coordinate", "protein",
                "feature", "cdr3_b_aa", "v_b_gene", "j_b_gene", "pgen", "radius", "regex")
name.output <- str_c(dir.output, "Minervina_metaclonotypes.tsv", sep = "/")
write.table(out, name.output, sep = "\t", row.names = FALSE, quote = FALSE)
