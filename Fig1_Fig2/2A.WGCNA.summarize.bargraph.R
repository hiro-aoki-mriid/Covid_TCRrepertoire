### Dummarize data for barplot of frequencies of each module

#load libraries
library(reshape2)
library(stringr)

#Input layer
dir <- "WGCNA_output_150"
Tcell <- "CD4" #Whether CD4 or CD8 repertoire is analysed
dir.output <- "2A"
name.output <- "WGCNA.summarize.barplot"
index_array <- c("sum_TP1", "sum_TP2", "sum_TP3", "sum_TP4") #Parameters for summarize 
class_order <- c("singlet", "0", "green", "invgreen", "brown", "invbrown", "yellow", "invyellow",
                 "turquoise", "invturquoise", "blue", "invblue")

############ Processing layer #####################################
dir.create(dir.output)

#load data
name.input <- str_c(dir, "WGCNA_extract",sep = "/") %>% str_c(Tcell, "th.correlation.0.85.csv", sep = ".")
d <- read.csv(name.input, header=TRUE)
sample.name <- str_split(d$name, pattern = "_", simplify = TRUE)
d$ID <- sample.name[,3]

#Change table format
wide_all <- data.frame()
for (m in 1:length(index_array)){
  value_i <- index_array[m]
  wide <- reshape2::dcast(d, class~ID, value.var = value_i)  
  wide$index <- value_i
  wide_all <- rbind(wide_all, wide)
}
wide_all[is.na(wide_all)] <- 0

wide_all$class <- factor(wide_all$class, levels = class_order)

wide_split <- split(wide_all, wide_all$class)
n <- length(wide_split)
temp <- wide_split[[1]]
for (i in 2:n) {
  temp <- merge(temp, wide_split[[i]], by = "index")
}
res <- temp

#Output
file.name <- str_c(dir.output, name.output, sep = "/") %>% str_c(Tcell, "csv", sep = ".") 
write.csv(res, file.name, row.names = F)
