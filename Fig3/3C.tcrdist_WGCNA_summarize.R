### Summarize tcrdist Metaclonotype pattern within each WGCNA module

#load libraries
library(foreach)
library(doParallel)
library(data.table)
library(stringr)
library(dplyr)
library(ggplot2)

#Input layer
dir.input <- 'WGCNA_tcrdist'
dir.output <- "3C"
cores <- 12 #Number of cores used for Foreach
Tcell <- "CD8"
wgcna.params <- c("0", "blue", "invblue", "yellow", "invyellow", "green", "invgreen",
                  "brown", "invbrown", "red", "invred", "turquoise", "invturquoise") #For CD8
types <- c("S", "other")

#################################### Processing layer ####################################################
###Define functions
createEmptyDf = function( nrow, ncol, colnames = c() ){
  data.frame( matrix( vector(), nrow, ncol, dimnames = list( c(), colnames ) ) )
}
tableread_fast = function(i, header=TRUE, quote="", sep="\t"){
  tmp = fread(i, header=header, sep=sep, quote=quote, nThread=32)
  tmp = as.data.frame(tmp)
  return(tmp)
}
##Main functions
WGCNA_summarize = function(file.name, dir.input, dir.output){
  #load
  name.input <- str_c(dir.input, file.name, sep = "/")
  d <- tableread_fast(name.input, header = TRUE, quote="\"", sep="\t")
  
  #Extract only clones assigned metaclonotype
  d <- dplyr::filter(d, protein != "non")
  d$name <- file.name
  
  return(d)
}

###Main modules
dir.create(dir.output)

files  <- list.files(dir.input, pattern="tcrdist3")
t<-proc.time()
cl <- makeCluster(cores)
registerDoParallel(cl)   
out <- foreach(file.name = files, .combine = "rbind",
               .packages=c("stringr", "data.table", "dplyr")) %dopar% {WGCNA_summarize(file.name, dir.input, dir.output)}
stopCluster(cl)
proc.time()-t

d <- as.data.frame(out)

##Simplify the classification: signelt or doublet / S-specific or other
#Judge whether clones were detected once (singlet) or more than two times (double)
d$detected <- "singlet"
d$detected[which(d$class != "singlet")] <- "doublet"
#Judge whether clones were spike reactive or not
d$protein_S <- "other"
d$protein_S[which(d$protein == "S")] <- "S"
#output
name.output <- str_c(dir.output, "tcrdist3_detected.allID.tsv", sep="/")
write.table(d, name.output, row.names = FALSE, quote = FALSE, sep = "\t")

##Count singlet / double S pecific clones in each participant
d_S <- dplyr::filter(d, protein == "S")
count_table <- table(d_S$name, d_S$detected)
#Output
name.output <- str_c(dir.output, "Sspecific.singlet.proportion.csv", sep="/")
write.csv(count_table, name.output, row.names = TRUE, quote = F)

##Summarize wgcna group in doublet tcrdist3-assigned clones
d_double <- dplyr::filter(d, detected== "doublet")
#Prepare dataframe for output
summary.wgcna <- createEmptyDf(length(wgcna.params), length(types), colnames = types)
#Summarize clone count
for(j in 1:length(wgcna.params)){
  k <- wgcna.params[j]
  d_double_sub <- dplyr::filter(d_double, class==k)
  for(i in 1:length(types)){
    l <- types[i]
    d_calc <- dplyr::filter(d_double_sub, protein_S==l)
    summary.wgcna[j,i] <- nrow(d_calc)
  }
}
summary.wgcna$wgcna <- wgcna.params
#Output
name.output <- str_c(dir.output, "Sspecific.doublet.wgcna.csv", sep="/")
write.csv(summary.wgcna, name.output, row.names = TRUE, quote = F)