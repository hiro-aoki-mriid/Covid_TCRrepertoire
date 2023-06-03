###Extract SARS-Cov-2 reactive clones determined by tcrdist3

#load libraryies

library(ggplot2)
library(extrafont)
loadfonts("win")
library(foreach)
library(doParallel)
library(stringr)
library(dplyr)
library(data.table)

#Input layer
dir.output <- "tmp/result/Fig4"
dir.input <- "tmp/result/intermediate/2_AIM/JoinTP_DifAbund_AIM"
cores <- 12
params <- c("TP3", "TP3FC")

#################################### Processing layer ####################################################
###Define functions
#Table.read.fast
tableread_fast = function(i, header=TRUE, quote="", sep=","){
  tmp = fread(i, header=header, sep=sep, quote=quote, nThread=32)
  tmp = as.data.frame(tmp)
  return(tmp)
}

### Merge scTCR clone table ###
##Combine beta-binomial result and AIM information
Combine <- function(file.name, extract, dir.AIM, dir.input){ 
  #load data
  name.input <- str_c(dir.input, file.name, sep = "/")
  data <- tableread_fast(name.input, header = TRUE, quote="\"", sep=",")
  data <- dplyr::filter(data, Resp == "2nd")
  if(nrow(data) > 0){
    data$name <- file.name
  }
  
  return(data)
}

###Main module
files  <- list.files(dir.input, pattern="NaraCOVID")
files <- files[str_detect(files, "CD8")]

##Combine beta-binomial result and AIM informations
t<-proc.time()
cl <- makeCluster(cores)
registerDoParallel(cl)   
out <- foreach(file.name = files,
        .packages=c("data.table", "stringr", "dplyr", "ggplot2"),
        .combine = rbind) %dopar% {Combine(file.name, extract, dir.AIM, dir.input)}
stopCluster(cl)
proc.time()-t

##Classify clones based on TP6-TP8 Expansion
out$TP6TP8dif[which(out$TP8 == 0)] <- "UD"
out$TP6TP8dif[which(out$TP6TP8dif == "Inc")] <- "UC"

##Summarize statistcal parameters

data_plot <- dplyr::select(out, c("TP6TP8dif", param))
names(data_plot) <- c("class", "param")

medians <- dplyr::group_by(data_plot, class) %>% dplyr::summarise(Median = median(param))

#Create histogram
ppi <- 600
file.name <- str_c(dir.output, "S4C", sep = "/") %>% str_c(param, "vsTP8State.tiff", sep = ".")
tiff(file.name, width=1.5*ppi, height=1.3*ppi, res=ppi)
p <- ggplot(data_plot, aes(x=log10(100*param), fill=class)) +
  geom_density(alpha = 0.6, size = 0.25) +
  scale_fill_manual(values=c(Exp="red", UC = "green", UD="cyan"))+
  theme_bw(base_size = 7) +
  theme(
    axis.title.x=element_blank(),
    axis.title.y=element_blank(),
    axis.text.x = element_text(family="Arial"),
    axis.text.y = element_text(family="Arial")) +
  annotate("text", x=0.2, y=2, parse = TRUE, size = 2.5, colour = "red",
           label=str_c("Exp.:", round(100*as.numeric(medians[1,2]), 3))) +
  annotate("text", x=0.2, y=1.65, parse = TRUE, size = 2.5, colour = "green",
           label=str_c("UC.:", round(100*as.numeric(medians[2,2]), 3))) +
  annotate("text", x=0.2, y=1.3, parse = TRUE, size = 2.5, colour = "cyan",
           label=str_c("ND.:", round(100*as.numeric(medians[3,2]), 3))) +
  guides(fill = FALSE)
print(p)
dev.off()



