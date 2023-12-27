#Prepare input file for searching clones that match metaclonotype

#load libraryies
library(dplyr)
library(stringr)
library(foreach)
library(doParallel)
library(data.table)

#Input layer
dir.input <- "tmp/result/intermediate/7_tcrdist/hla_restricted_meta_clonotypes"
dir.output <- "tmp/result/intermediate/7_tcrdist/Metaclonotype_query"
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

    #Collect information needed for metaclonotype search
    out <- dplyr::select(table, c("cdr3_b_aa", "v_b_gene", "j_b_gene", "pgen", "radius", "regex", "nsubject"))
    out$set <- "1"
    out$protein <- Epitope.inf[3]
    out$pathogen <- Epitope.inf[4]
    out$HLA <- Epitope.inf[1]
    out$epitope <- Epitope.inf[2]
    out$feature <- str_c("1E5", out$v_b_gene, out$cdr3_b_aa, out$radius, out$regex, sep = "+")
    out <- dplyr::select(out, c("set", "protein", "pathogen", "HLA", "epitope",
                                "feature", "cdr3_b_aa", "v_b_gene", "j_b_gene", "pgen", "radius", "regex", "nsubject"))
    
    return(out)
  }
}

###Main module
files  <- list.files(dir.input, pattern="ranked_centers_bkgd_ctlr_1E5.tsv")
dir.create(dir.output)
t<-proc.time()
cl <- makeCluster(cores)
registerDoParallel(cl)   
out <- foreach(file.name = files,
               .combine = rbind, .packages=c("ggplot2", "extrafont", "stringr", "dplyr", "data.table")) %dopar% {Main(file.name, dir.input, dir.output)}
stopCluster(cl)
proc.time()-t

out <- as.data.frame(out)
names(out) <- c("set", "protein", "pathogen", "HLA", "protein_coordinate",
                "feature", "cdr3_b_aa", "v_b_gene", "j_b_gene", "pgen", "radius", "regex", "nsubject")
name.output <- str_c(dir.output, "vdjdb_metaclonotypes.tsv", sep = "/")

out$protein <- str_replace(out$protein, "Spike", "S")

###Split meta-clonotype by their source
out_tet <- dplyr::filter(out, str_detect(protein_coordinate, "Tet_"))
out_iedb <- dplyr::filter(out, str_detect(protein_coordinate, "IEDB_"))
out_vdjdb <- dplyr::filter(out, str_detect(protein_coordinate, "vdjdb_"))

#Output
write.table(out_vdjdb, name.output, sep = "\t", row.names = FALSE, quote = FALSE)
name.output.tet <- str_replace(name.output, "vdjdb", "Tet")
write.table(out_tet, name.output.tet, sep = "\t", row.names = FALSE, quote = FALSE)
name.output.iedb <- str_replace(name.output, "vdjdb", "iedb")
write.table(out_iedb, name.output.iedb, sep = "\t", row.names = FALSE, quote = FALSE)