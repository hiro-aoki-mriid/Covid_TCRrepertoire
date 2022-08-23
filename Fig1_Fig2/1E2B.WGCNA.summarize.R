###Calculate clonality of TCR repretoire

#Load Libraries

library(foreach)
library(doParallel)
library(stringr)
library(data.table)
library(dplyr)

#Input layer
dir <- "WGCNA_output_150"
Tcells <- c("CD8", "CD4") #Whether CD4 or CD8 repertoire is analysed
dir.outputs <- c("1E_S2A-D", "2B_S2E-H")
cores <- 12
name.output <- "WGCNA.summarize.temporal"

#################################### Processing layer ####################################################
for(i in 1:length(Tcells)){
  
  Tcell <- Tcells[i]
  dir.output <- dir.outputs[i]
  dir.create(dir.output)
  
  #load data
  name.input <- str_c(dir, "WGCNA_extract",sep = "/") %>% str_c(Tcell, "th.correlation.0.85.csv", sep = ".")
  d <- read.csv(name.input, header=TRUE)
  sample.name <- str_split(d$name, pattern = "_", simplify = TRUE)
  d$ID <- sample.name[,3]
  d_freqs <- dplyr::select(d, c("sum_TP1", "sum_TP2", "sum_TP3", "sum_TP4"))
  
  #Calculate maximum changes in frequency
  delta <- function(x) {return (max(x) - min(x))}
  d$max_change <- apply(d_freqs, 1, delta)
  
  #Change table format
  point_array <- c("sum_TP1", "sum_TP2", "sum_TP3", "sum_TP4", "max_change")
  wide_all <- data.frame()
  for(point in point_array){
    wide <- reshape2::dcast(d, ID ~ class, value.var = point)  
    wide$Point <- point
    wide_all <- rbind(wide_all, wide)
  }
  
  #Output
  name.out <- str_c(dir.output, name.output, sep = "/") %>% str_c(Tcell, "csv", sep = ".")
  write.csv(wide_all, name.out, row.names = F)
}



