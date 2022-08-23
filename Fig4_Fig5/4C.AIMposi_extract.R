###Define AIM+ clones by their enrichment into AIM relative to NonNV

#load libraryies

library(ggplot2)
library(extrafont)
loadfonts("win")
library(foreach)
library(doParallel)
library(stringr)
library(dplyr)
library(reshape2)
library(scales)

#Input layer
dir <- "Join"
dir.scatter <- "4C"
dir.collapsed <- "S3BC"
dir.table <- "AIM_Table"
name.csv <- "cell_count.csv"
cores <- 12
threshold <- 4 #Threshold of AIM/NonNV ratio for AIM+ clones

#################################### Processing layer ####################################################
###Define functions
#Generate Scatter plot for AIM-NonNV overlap
Scatter <- function(data, dir.scatter, name_out){
  #Prepare data
  data_sub <- dplyr::filter(data, NonNV > 0 & AIM > 0)
  ppi <- 600
  image.file <- str_c(dir.scatter, name_out, sep = "/") %>% str_c("tiff", sep = ".")
  tiff(image.file, width=1*ppi, height=1*ppi, res=ppi)
  p <- ggplot(data_sub, aes(x=log2(NonNVcount), y=log2(AIMcount), colour = class)) + 
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

#Generate collapsed plot for AIM-NonNV overlap
Collapsed <- function(data, dir.collapsed, name_out){
  #Prepare data
  data_sub <- dplyr::select(data, c("AIM", "NonNV", "class"))
  data_sub$class[which(data_sub$AIM == 0)] <- "nonNVonly"
  data_sub$class[which(data_sub$NonNV == 0)] <- "AIMonly"
  data_sub <- group_by(data_sub, class)
  sum_table <- summarize(data_sub, Posi = sum(AIM), Nega = sum(NonNV))
  sum_table <- melt(sum_table, id.vars = "class", variable.name = "organ", value.name = "freq") 
  
  #Create collapsed plot
  ppi <- 600
  image.file <- str_c(dir.collapsed, name_out, sep = "/") %>% str_c("tiff", sep = ".")
  tiff(image.file, width=0.8*ppi, height=0.85*ppi, res=ppi)
  #Make plot
  p <- ggplot(sum_table, aes(x=organ, y=100*freq,
                             fill = factor(class, levels=c("nonNVonly", "NonNV", "AIM", "AIMonly")))) + 
    geom_bar(stat = "identity", width = 0.5) +
    scale_fill_manual(values=c(AIMonly="red", NonNV="skyblue", AIM="pink", nonNVonly="darkblue")) +
    scale_x_discrete(expand = c(0.2, 0.2)) +
    theme_bw(base_size = 6) +
    theme(
      axis.title.x=element_blank(), axis.title.y=element_blank(),
      axis.text.x = element_text(family="Arial"),
      axis.text.y = element_text(family="Arial"),
      plot.margin= unit(c(0.25, 0.65, 0.25, 0.25), "lines")) +
    guides(fill="none")
  print(p)
  dev.off() 
  
  #Output
  sum_table$name <- name_out
  return(sum_table)
}

#Main function
Main <- function(file.name, dir, dir.scatter, dir.collapsed, dir.table,
                 cell.count.table, threshold){
  #load data
  name.input <- str_c(dir, file.name, sep = "/")
  data <- read.table(name.input, header = TRUE)
  names(data) <- c("count", "freq", "cdr3nt", "cdr3aa", "v", "d", "j", "VEnd", "DStart", "DEnd", "JStart",
                   "peak", "occurences", "sampling.p", "AIM", "NonNV")
  
  #load data of cell count
  name_out <- str_replace(file.name, ".join.strict.table.txt", "")
  name_aim <- str_replace(name_out, "COVIDAIM", "COVID_AIM_AIM")
  name_nonnv <- str_replace(name_out, "COVIDAIM", "COVID_AIM_NonNV")
  count_aim <- dplyr::filter(cell.count.table, name == name_aim)$Cell_count
  count_nonnv <- dplyr::filter(cell.count.table, name == name_nonnv)$Cell_count
  int <- count_aim / count_nonnv #ratio of total cell count
  
  ##Define AIM+ clones based on AIM/nonAIM ratio
  #Calculate AIM/nonAIM ratio using the cell count for each clone
  data$AIMcount <- data$AIM * count_aim + 0.001 #Add to avoid dividing by zero
  data$NonNVcount <- data$NonNV * count_nonnv + 0.001 #Add to avoid dividing by zero
  data$ratio <- data$AIMcount / data$NonNVcount
  #Define and extract AIM+ clones
  data$class <- "NonNV"
  data$class[which(data$ratio > int*threshold)] <- "AIM"
  
  ##Generate Scatter plot for AIM-NonNV overlap
  Scatter(data, dir.scatter, name_out)

  ##Generate collapsed plot for AIM-NonNV overlap
  sum_table <- Collapsed(data, dir.collapsed, name_out)
  
  #Output table of AIM+ clones
  name.output <- str_c(dir.table, file.name, sep = "/") %>%
    str_remove(".join.strict.table.txt") %>% 
    str_c("th", threshold, "AIMposi.csv", sep = ".")
  write.csv(data, name.output, row.names = FALSE)
  
  return(sum_table)
}

###Main module
files  <- list.files(dir, pattern="join.strict.table.txt")
cell.count.table <- read.csv(name.csv, header = TRUE)
dir.create(dir.table)
dir.create(dir.scatter)
dir.create(dir.collapsed)
t<-proc.time()
cl <- makeCluster(cores)
registerDoParallel(cl)   
out <- foreach(file.name = files,
        .combine = rbind, .packages=c("ggplot2", "extrafont", "stringr", "dplyr", "reshape2")) %dopar% {Main(file.name, dir, dir.scatter, dir.collapsed, dir.table,
                                                                                                 cell.count.table, threshold)}
stopCluster(cl)
proc.time()-t

#Output AIM separation stats
name.output <- str_c(dir.collapsed, "AIM-NonNV.overlap.csv", sep ="/")
write.csv(out, name.output, row.names = FALSE)
