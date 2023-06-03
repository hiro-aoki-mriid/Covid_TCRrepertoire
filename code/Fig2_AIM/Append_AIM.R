###Determine kinetics of AIM+ clones based on the responding pattern

#load libraryies
library(ggplot2)
library(foreach)
library(doParallel)
library(stringr)
library(dplyr)
library(data.table)

#Input layer
dir.input <- "tmp/result/intermediate/1_beta-binomial/JoinTP_DifAbund"
dir.AIM <- "tmp/result/intermediate/2_AIM/AIMJoin"
dir.output <- "tmp/result/intermediate/2_AIM/JoinTP_DifAbund_AIM"
dir.scatter <- "tmp/result/Fig2/2C"
threshold <- 16 #AIM/NonNV ratio
cores <- 12

#################################### Processing layer ####################################################
cell.count.table <- read.csv(name.csv, header = TRUE)

###Define functions
tableread_fast = function(i, header=TRUE, quote="", sep="\t"){
  tmp = fread(i, header=header, sep=sep, quote=quote, nThread=32)
  tmp = as.data.frame(tmp)
  return(tmp)
}

#Generate Scatter plot for enrichment in AIM+ repertoire
Scatter <- function(data, name_out){
  #Prepare data
  data_sub <- dplyr::filter(data, AIM > 0)
  ppi <- 600
  tiff(name_out, width=1*ppi, height=1*ppi, res=ppi)
  p <- ggplot(data_sub, aes(x=log10(100*NonNV), y=log10(100*AIM), colour = class)) + 
    geom_point(size=0.25, shape = 16) + 
    scale_colour_manual(values=c(AIM="red", NonNV="blue"))+
    theme_bw(base_size = 6) +
    theme(
      axis.title.x=element_blank(), axis.title.y=element_blank(),
      axis.text.x = element_text(family="Arial"),
      axis.text.y = element_text(family="Arial"),
      plot.margin= unit(c(0.25, 0.3, 0.25, 0.25), "lines")) +
    guides(colour= "none")
  print(p)
  dev.off()
}

##Combine beta-binomial result and AIM informations
Combine <- function(file.name, extract, dir.AIM, dir.input){ 
  #AIM assay, #vdjtools.output?̓ǂݍ???
  name.input <- str_c(dir.AIM, file.name, sep = "/") %>%
    str_replace(".DifAbund.csv", ".join.strict.table.txt") %>%
    str_replace("NaraCOVID", "COVIDAIM")
  d_aim <- read.table(name.input, header = TRUE)
  names(d_aim) <- c("count", "freq", "cdr3nt", "cdr3aa", "v", "d", "j", "VEnd", "DStart", "DEnd", "JStart",
                    "peak", "occurences", "sampling.p", "AIM", "NonNV")
  
  ##Define AIM+ clonotypes based on the enrichment in AIM+ repertoire
  d_aim$class <- "NonNV"
  d_aim$class[which((d_aim$AIM + 1/100000/2) / (d_aim$NonNV + 1/100000/2) > threshold)] <- "AIM"
  
  #Generate Scatter plot for enrichment in AIM+ repertoire
  name_out <- str_c(dir.scatter, file.name, sep = "/") %>% str_replace("DifAbund.csv", "Scatter.tiff")
  Scatter(d_aim, name_out)
  
  #Merge AIM information
  #Define clonotype by CDR3nt
  d_aim$ntvj <- str_c(d_aim$cdr3nt, d_aim$v, d_aim$j, sep = "_")
  AIM <- dplyr::select(d_aim, c("class", "ntvj"))
  #load repertoire table and merge
  name.input <- str_c(dir.input, file.name, sep = "/")
  data <- tableread_fast(name.input, header = TRUE, quote="\"", sep=",")
  data_AIM <- merge(data, AIM, by = "ntvj", all.x = TRUE)
  data_AIM[is.na(data_AIM)] <- "non"
  
  #Output
  name.output <- str_c(dir.output, file.name, sep = "/") %>% str_replace("csv", "aim.csv")
  write.csv(data_AIM, name.output, row.names = FALSE)
}

###Main module
dir.create(dir.output)
dir.create(dir.scatter)
files  <- list.files(dir.input, pattern="NaraCOVID")

##Combine beta-binomial result and AIM informations
t<-proc.time()
cl <- makeCluster(cores)
registerDoParallel(cl)   
foreach(file.name = files,
        .packages=c("data.table", "stringr", "dplyr", "ggplot2"),
        .combine = rbind) %dopar% {Combine(file.name, extract, dir.AIM, dir.input)}
stopCluster(cl)
proc.time()-t

