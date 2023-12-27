#Load packages
library(ComplexHeatmap)
library(circlize)
library(stringr)
library(dplyr)
library(qvalue)


#Input layer
Cells <- c("CD4", "CD8")

dir.input <- "tmp/metadata"
dir.output <- "tmp/result/Fig4/4AB_Correlation"
merge_threshold <- 0.8

#################### Define functions ##########################
createEmptyDf = function( nrow, ncol, colnames = c() ){
  data.frame( matrix( vector(), nrow, ncol, dimnames = list( c(), colnames ) ) )
}

Correlation <- function(data1, data2, name1, name2, dir.name){
  ##Prepare output
  #directory
  dir.name.heatmap <- dir.name
  dir.create(dir.name.heatmap)
  #file
  csv.name <- str_c(dir.name, name1, sep = "/") %>%
    str_c("vs", name2, ".correlation.result.csv", sep = '')
  heatmapname <- str_c(dir.name.heatmap, name1, sep = "/") %>%
    str_c("vs", name2, ".heatmap.tiff", sep = '')
  
  ###Call parameter names
  param1 <- names(data1)
  param2 <- names(data2)
  
  ###Calculate correlation coefficient
  ##Prepare output
  result.csv <- createEmptyDf(length(param1)*length(param2), 3, colnames = c("rvalue", "r_abs", "pvalue") )
  result.heatmap <- createEmptyDf(length(param1), length(param2), colnames = param2 )
  x.param.names <- vector()
  y.param.names <- vector()
  ##Iterative process
  k <- 0
  for (i in 1:length(param1)){
    for (j in 1:length(param2)){
      k <- k+1
      #Call parameter name
      x.param_name <- param1[i]
      y.param_name <- param2[j]
      #Spearman test
      cor_t <- cor.test(data1[,i], data2[,j], use="pairwise", method="spearman")
      #Store result for csv output
      if(i > j){
        result.csv[k,1] <- as.vector(cor_t$estimate)
        result.csv[k,2] <- abs(result.csv[k,1])
        result.csv[k,3] <- cor_t$p.value
      }
      x.param.names <- c(x.param.names, x.param_name)
      y.param.names <- c(y.param.names, y.param_name)
      #Store result for heatmap output
      result.heatmap[i,j] <- cor_t$estimate
    }
  }
  
  ##Add information for correlation coefficient table
  result.csv$x.param <- x.param.names
  result.csv$y.param <- y.param.names
  #Calculate q-value
  result.csv <- na.omit(result.csv)
  qobj<- qvalue(p = result.csv$pvalue, lambda=0)
  result.csv$qvalue <- qobj$qvalues
  #Append significance asterisks
  result.csv$signif <- "n.s."
  result.csv$signif[which(result.csv$qvalue < 0.05)] <- "*"
  result.csv$signif[which(result.csv$qvalue < 0.01)] <- "**"
  result.csv$signif[which(result.csv$qvalue < 0.001)] <- "***"
  #Output
  write.csv(result.csv, csv.name, row.names = FALSE)
  
  ###Create heatmap
  row.names(result.heatmap) <- param1
  d.mat <- as.matrix(result.heatmap)
  #Prepare color function
  col_fun = colorRamp2(c(-1, -0.5, 0, 0.5, 1), c("blue","cyan","white","magenta", "red"))
  #Output
  ppi <- 300
  tiff(heatmapname, width=3.3*ppi, height=2.2*ppi, res = ppi)
  p <- Heatmap(d.mat, name = "name", col = col_fun,
               cluster_rows = TRUE,
               cluster_columns = TRUE,
               row_dend_width = unit(0.5, "cm"),
               column_dend_height = unit(0.5, "cm"),
               row_names_gp = gpar(fontsize = 6, family="Arial"),
               show_column_names = FALSE,
               heatmap_legend_param =list(labels_gp = gpar(fontsize=6, family="Arial"),
                                          title = "Z value", 
                                          title_gp = gpar(fontsize=6, family="Arial"),
                                          legend_direction = "vertical",
                                          legend_width = unit(3, "cm"),
                                          legend_height = unit(0.15, "cm")),
               )
  q <- draw(p, heatmap_legend_side = "left")
  print(q)
  dev.off()
}

#################### Processing ##########################
dir.create(dir.output, recursive = TRUE)

for(Cell in Cells){
  #load data
  name.input <- str_c(dir.input, "Metaanalysis_", sep = "/") %>% str_c(Cell, ".txt") 
  d <- read.table(name.input, header = TRUE)
  
  #Extract repertoire data
  a <- dplyr::select(d, -c("ID"))
  
  Correlation(a, a, Cell, Cell, dir.output)
}

