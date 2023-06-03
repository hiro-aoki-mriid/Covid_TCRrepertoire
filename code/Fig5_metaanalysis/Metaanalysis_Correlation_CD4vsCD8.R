# パッケージの読み込み
library(ComplexHeatmap)
library(circlize)
library(stringr)
library(dplyr)
library(qvalue)

#Fontの導入
library(extrafont)
loadfonts("win")

#Input layer

dir.input <- "tmp/result/Fig5"
dir.output <- "tmp/result/Fig5/5C"
merge_threshold <- 0.8

#################### Define functions ##########################
createEmptyDf = function( nrow, ncol, colnames = c() ){
  data.frame( matrix( vector(), nrow, ncol, dimnames = list( c(), colnames ) ) )
}

Correlation <- function(data1, data2, name1, name2, dir.name){
  ##Prepare output
  #file
  csv.name <- str_c(dir.name, name1, sep = "/") %>%
    str_c("vs", name2, ".correlation.result.csv", sep = '')
  heatmapname <- str_c(dir.name, name1, sep = "/") %>%
    str_c("vs", name2, ".heatmap.tiff", sep = '')
  
  ###Call parameter names
  param1 <- names(data1)
  param2 <- names(data2)
  
  ###Calculate correlation coefficient
  ##Prepare output
  result.csv <- createEmptyDf(length(param1)*length(param2), 3, colnames = c("rvalue", "r_abs", "pvalue") )
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
      result.csv[k,1] <- as.vector(cor_t$estimate)
      result.csv[k,2] <- abs(result.csv[k,1])
      result.csv[k,3] <- cor_t$p.value
      x.param.names <- c(x.param.names, x.param_name)
      y.param.names <- c(y.param.names, y.param_name)
    }
  }
  
  ##Add information for correlation coefficient table
  result.csv$cd4.param <- x.param.names
  result.csv$cd8.param <- y.param.names
  #Calculate q-value
  result.csv <- na.omit(result.csv)
  qobj<- qvalue(p = result.csv$pvalue)
  result.csv$qvalue <- qobj$qvalues
  #Append significance asterisks
  result.csv$signif <- ""
  result.csv$signif[which(result.csv$qvalue < 0.05)] <- "*"
  result.csv$signif[which(result.csv$qvalue < 0.01)] <- "**"
  result.csv$signif[which(result.csv$qvalue < 0.001)] <- "***"
  #Output
  write.csv(result.csv, csv.name, row.names = FALSE)
  
  ##Prepare matrix for heatmap
  result.heat <- matrix(result.csv$rvalue, ncol = length(param2), byrow = TRUE)
  colnames(result.heat) <- str_remove(param2, "X")
  row.names(result.heat) <- str_remove(param1, "X")
  result.heat.q <- matrix(result.csv$signif, ncol = length(param2), byrow = TRUE)
  colnames(result.heat.q) <- param2
  row.names(result.heat.q) <- param1
  
  ###Create heatmap
  #Prepare color function
  col_fun = colorRamp2(c(-1, -0.5, 0, 0.5, 1), c("blue","cyan","white","magenta", "red"))
  #Output
  ppi <- 600
  tiff(heatmapname, width=3.3*ppi, height=2.8*ppi, res = ppi)
  p <- Heatmap(result.heat, name = "name", col = col_fun,
               cell_fun = function(j, i, x, y, width, height, fill) {
                 grid.text(sprintf(result.heat.q[i, j]), x, y, gp = gpar(fontsize = 10))
               },
               cluster_rows = TRUE,
               cluster_columns = TRUE,
               row_dend_width = unit(0.5, "cm"),
               column_dend_height = unit(0.5, "cm"),
               row_names_gp = gpar(fontsize = 6, family="Arial"),
               column_names_gp = gpar(fontsize = 6, family="Arial"),
               heatmap_legend_param =list(labels_gp = gpar(fontsize=6, family="Arial"),
                                          title = "Z value", 
                                          title_gp = gpar(fontsize=6, family="Arial"),
                                          legend_direction = "vertical"))
  q <- draw(p, heatmap_legend_side = "left")
  print(q)
  dev.off()
}

#################### Processing ##########################
dir.create(dir.output, recursive = TRUE)

##load data
#CD4 data
name.input <- str_c(dir.input, "Merged.Metaanalysis_CD4.csv", sep = "/")
d_cd4 <- tableread_fast(name.input, header = TRUE, sep=",")
#CD8 data
name.input <- str_c(dir.input, "Merged.Metaanalysis_CD8.csv", sep = "/")
d_cd8 <- tableread_fast(name.input, header = TRUE, sep=",")
  
Correlation(d_cd4, d_cd8, "CD4", "CD8", dir.output)


