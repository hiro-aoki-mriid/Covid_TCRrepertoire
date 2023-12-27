###Create Venn diagram for the correspondence between tcrdist, tetramer-TCrseq and AIM assay

#load libraryies

library(ggplot2)
library(extrafont)
loadfonts("win")
library(dplyr)
library(ggvenn)
library(data.table)
library(stringr)

#Input layer
name.input <- "tmp/result/Fig7/tcrdist_clones.csv"
dir.output <- "tmp/result/Fig7/S7A"
epitopes <- c("S269", "S448", "S919", "S1208")

############################### Define functions ########################################
tableread_fast = function(i, header=TRUE, quote="", sep=","){
  tmp = fread(i, header=header, sep=sep, quote=quote, nThread=32)
  tmp = as.data.frame(tmp)
  return(tmp)
}

##FigS7B, NF9 tcrdist vs tetramer
Venn <- function(data, epitope, dir.output){
  ppi <- 600
  name.output <- str_c(dir.output, epitope, sep = "/") %>% str_c(".venn.tiff")
  tiff(name.output, width=1.4*ppi, height=1*ppi, res = ppi)
  p <- ggvenn(data, auto_scale = TRUE,
              fill_color = c("blue", "red"),
              show_percentage = FALSE,
              stroke_size = 0.3,
              text_size = 0,
              show_outside = "none",
              set_name_size = 0
  )
  print(p)
  dev.off()
}


############################### Processing ########################################
dir.create(dir.output, recursive = TRUE)
#load tcrdist3 data
d <- tableread_fast(name.input, header = TRUE, quote="\"", sep=",")
d$clone.name <- str_c(d$ntvj, d$sample)

output.all <- data.frame()
for(epitope in epitopes){
  #Extract necessary informations
  d_sub <- dplyr::select(d, contains(epitope))
  d_sub$clone.name <- d$clone.name
  names(d_sub) <- c("tet", "motif", "clone.name")
  #Prepare dataframe
  d_Tet <- dplyr::filter(d_sub, tet == "Posi")
  d_motif <- dplyr::filter(d_sub, motif == "Posi")
  data <- list(motif = d_motif$clone.name, tet = d_Tet$clone.name)
  ##Create Venn plot
  Venn(data, epitope, dir.output)
  ##Summarize
  out <- as.data.frame(table(d_sub$tet, d_sub$motif))
  names(out) <- c("tetramer", "motif", "count")
  out$epitope <- epitope
  output.all <- rbind(output.all, out)
}

#Output
name.output <- str_c(dir.output, "Venn.summary.csv", sep = "/")
write.csv(output.all, name.output, row.names = FALSE)
