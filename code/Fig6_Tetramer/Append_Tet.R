###Determine kinetics of AIM+ clones based on the WGCNA results

#load libraryies
library(ggplot2)
library(foreach)
library(doParallel)
library(stringr)
library(dplyr)
library(data.table)

#Input layer
dir.input <- "tmp/result/intermediate/2_AIM/JoinTP_DifAbund_AIM"
dir.Tet <- "tmp/data/original/raw.data"
dir.output <- "tmp/result/intermediate/6_Tet/JoinTP_DifAbund_AIM_Tet"
cores <- 12

epitope_array <- c("NF9", "QI9")
point_array <- c("TP2", "TP3", "TP8")
points_all <- c("TP1", "TP2", "TP3", "TP4", "TP5", "TP6", "TP8", "TP9", "TP10", "TP11")
threshold <- 16 #Threshold of Tetposi/total ratio for Tet+ clones

#################################### Processing layer ####################################################
###Define functions
tableread_fast = function(i, header=TRUE, quote="", sep="\t"){
  tmp = fread(i, header=header, sep=sep, quote=quote, nThread=32)
  tmp = as.data.frame(tmp)
  return(tmp)
}

##Combine beta-binomial result and AIM information
Combine <- function(file.name, extract, dir.AIM, dir.input){ 
  #load data
  name.input <- str_c(dir.input, file.name, sep = "/")
  data <- tableread_fast(name.input, header = TRUE, quote="\"", sep=",")
  
  #Search data
  name.input <- str_remove(file.name, "NaraCOVID_CD8_") %>%
    str_remove(".DifAbund.aim.csv")
  
  #Check whether tetramer test is performed or not
    for(epitope in epitope_array){
      for(point in point_array){
        query <- str_c(epitope, point, sep = "_")
        search.input <- str_c(dir.Tet, "COVIDTet", sep = "/") %>%
          str_c(query, "CD8", name.input, sep = "_") %>% 
          str_c(".txt")
        if(file.exists(search.input) == TRUE){
          Tet <- tableread_fast(search.input, header = TRUE, quote="\"", sep="\t")
          Tet$ntvj <- str_c(Tet$cdr3nt, Tet$v, Tet$j, sep = "_")
          #Define clonotype by CDR3nt
          Tet <- dplyr::select(Tet, c("freq", "ntvj"))
          names(Tet) <- c("freq.Tet", "ntvj")
          #Merge AIM information
          data <- merge(data, Tet, by = "ntvj", all.x = TRUE)
          data[is.na(data)] <- "0"
          #Thresholding
          data_points <- dplyr::select(data, points_all)
          data_points$max <- apply(data_points, 1, max)
          data_sub <- as.data.frame(cbind(data_points$max, data$freq.Tet))
          names(data_sub) <- c("total", "freq.Tet") 
          data_sub$total <- as.numeric(data_sub$total)
          data_sub$freq.Tet <- as.numeric(data_sub$freq.Tet)
          data_sub$class <- "Nega"
          names(data_sub) <- c("total", "freq.Tet", "class")
          data_sub$class[which((data_sub$freq.Tet + 1/100000/2) / (data_sub$total + 1/100000/2) > threshold)] <- "Posi"
          #Return thresholding results
          data$class.tet <- data_sub$class
          data <- dplyr::select(data, -c("freq.Tet"))
          names(data) <- c(names(data)[1:(ncol(data)-1)], query)
        } else {
          data$class.Tet <- "Nega"
          names(data) <- c(names(data)[1:(ncol(data)-1)], query)
        }
      }
    }
    
    #Deterine NF9/QI9 positivity
    data$NF9 <- "Nega"
    data$NF9[which(data$NF9_TP2 == "Posi" | data$NF9_TP3 == "Posi" | data$NF9_TP8 == "Posi")] <- "Posi"
    data$QI9 <- "Nega"
    data$QI9[which(data$QI9_TP2 == "Posi" | data$QI9_TP3 == "Posi" | data$QI9_TP8 == "Posi")] <- "Posi"
    data$tetramer <- "Nega"
    data$tetramer[which(data$NF9 == "Posi")] <- "NF9"
    data$tetramer[which(data$QI9 == "Posi")] <- "QI9"
    data$tetramer[which(data$NF9 == "Posi" & data$QI9 == "Posi")] <- "Dual"
    
    #Output 
    name.output <- str_c(dir.output, file.name, sep = "/") %>% str_replace("csv", "tetramer.csv")
    write.csv(data, name.output, row.names = FALSE)
}

###Main module
dir.create(dir.output, recursive = TRUE)
files  <- list.files(dir.input, pattern="NaraCOVID_CD8")

##Combine beta-binomial result and AIM informations
t<-proc.time()
cl <- makeCluster(cores)
registerDoParallel(cl)   
foreach(file.name = files,
               .packages=c("data.table", "stringr", "dplyr", "ggplot2"),
               .combine = rbind) %dopar% {Combine(file.name, extract, dir.AIM, dir.input)}
stopCluster(cl)
proc.time()-t

