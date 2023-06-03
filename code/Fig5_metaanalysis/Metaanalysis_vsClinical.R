# パッケージの読み込み
library(ComplexHeatmap)
library(circlize)
library(stringr)
library(dplyr)
library(ggplot2)
library(ggsignif)
library(data.table)
library(qvalue)

#Fontの導入
library(extrafont)
loadfonts("win")

#Input layer
Cells <- c("CD4", "CD8")

dir.input <- "tmp/result/Fig5"
name.clinical <- "tmp/metadata/Metaanalysis_Clinical.txt"
dir.name.heatmap <- "tmp/result/Fig5/5AB"
dir.name.scatter <- "tmp/result/Fig5/S5C"
dir.name.bar <- "tmp/result/Fig5/Association_Factor"

#Type of parameters
param.fact <- c("Sex")
param.num <- c("Age", "Ancestral.P1", "Ancestral.P2", "Ancestral.P3", "Ancestral.P4", "Ancestral.P5",
               "Ancestral.P6", "Ancestral.P8", "Ancestral.P9", "Ancestral.P10", "Ancestral.P11",
               "Omicron.P3", "Omicron.P6", 
               "AE.1st", "AE.2nd", "AE.3rd")

#####################   Prepare functions   ###########################################
createEmptyDf = function( nrow, ncol, colnames = c() ){
  data.frame( matrix( vector(), nrow, ncol, dimnames = list( c(), colnames ) ) )
}

tableread_fast = function(i, header=TRUE, sep="\t"){
  tmp = fread(i, header=header, sep=sep, nThread=32)
  tmp = as.data.frame(tmp)
  return(tmp)
}

##Scatter plot
Scatter <- function(for.plot, x.param_name, y.param_name, dir.name.scatter, Cell){
  #Prepare output
  ppi <- 600
  image.file <- str_c(dir.name.scatter, Cell, sep = "/") %>%
    str_c(str_c(x.param_name, "-", y.param_name), "png", sep = '.') 
  png(image.file, width=1*ppi, height=0.9*ppi, res=ppi)
  
  names(for.plot) <- c("x.param", "y.param")
  for.plot <- na.omit(for.plot)
  
  p <- ggplot(for.plot, aes(x=log10(x.param), y=y.param)) + 
    geom_point(size = 0.5) +
    theme_bw(base_size = 6) +
    labs(x=x.param_name, y=y.param_name) +
    theme(axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.x = element_text(size=6, family="Arial"),
          axis.text.y = element_text(size=6, family="Arial")) +
    guides(colour="none")
  print(p)
  dev.off() 
}

###main Function
Correlation <- function(data1, data2, name1, name2,
                        dir.name.heatmap, dir.name.scatter){
  ##Prepare output
  #directory
  dir.create(dir.name.heatmap)
  dir.create(dir.name.scatter)
  #file
  csv.name <- str_c(dir.name.heatmap, name1, sep = "/") %>%
    str_c("vs", name2, ".correlation.result.csv", sep = '')
  heatmapname <- str_c(dir.name.heatmap, name1, sep = "/") %>%
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
      
      #Create Scatter plot when p-value is lower than 0.01
      if(as.numeric(cor_t$p.value) < 0.01){
        for.plot <- data.frame(data1[,i], data2[,j])
        Scatter(for.plot, x.param_name, y.param_name, dir.name.scatter, name1)
      } else {
      }
    }
  }
  
  ##Add information for correlation coefficient table
  result.csv$x.param <- x.param.names
  result.csv$y.param <- y.param.names
  #Calculate q-value
  result.csv <- na.omit(result.csv)
  qobj<- qvalue(p = result.csv$pvalue)
  result.csv$qvalue <- qobj$qvalues
  #Judge significance
  result.csv$significance <- ""
  result.csv$significance[which(result.csv$qvalue < 0.05)] <- "*"
  result.csv$significance[which(result.csv$qvalue < 0.01)] <- "**"
  result.csv$significance[which(result.csv$qvalue < 0.001)] <- "***"
  #Output
  write.csv(result.csv, csv.name, row.names = FALSE)
  
  ##Prepare matrix for heatmap
  result.heat <- matrix(result.csv$rvalue, ncol = length(param2), byrow = TRUE)
  colnames(result.heat) <- param2
  row.names(result.heat) <- str_remove(param1, "X")
  result.heat.q <- matrix(result.csv$significance, ncol = length(param2), byrow = TRUE)
  colnames(result.heat.q) <- param2
  row.names(result.heat.q) <- param1
  
  ###Create heatmap
  #Prepare color function
  col_fun = colorRamp2(c(-1, -0.5, 0, 0.5, 1), c("blue","cyan","white","magenta", "red"))
  #Output
  ppi <- 300
  tiff(heatmapname, width=5.8*ppi, height=2.8*ppi, res = ppi)
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
  q <- draw(p, heatmap_legend_side = "right")
  print(q)
  dev.off()
}

##Bar Plot
BarPlot <- function(d_sub, dir.name.bar, x.param_name, y.param_name, test, p.value){
  #Prepare output
  ppi <- 150
  image.file <- str_c(dir.name.bar, x.param_name, sep = "/") %>%
    str_c("-", y.param_name, ".png", sep = '') 
  png(image.file, width=3*ppi, height=2.4*ppi, res=ppi)
  
  ##Plotting graphs
  #p-valueを表示する位置を調整
  ymax <- max(log10(d_sub$rep))
  #pvalueを表示する位置
  if(ymax > 0){
    ylim <- ymax*1.5
    ypos <- ylim*0.85
  }else{
    ylim <- ymax*0.7
    ypos <- ylim*1.1
  }
  d_sub <- na.omit(d_sub)
  p <- ggplot(d_sub, mapping = aes(x=as.factor(clinical), y=log10(repertoire))) +
    geom_boxplot(outlier.shape = NA, size=0.2) + 
    geom_jitter(
      data=d_sub,
      aes(x= as.factor(clinical),
          y= log10(repertoire)),
      stat="identity",
      size=2,
      shape = 1,
      position=position_jitter(0.17)) +
    scale_y_continuous(limits = c(NA,ylim)) +
    theme_bw(base_size = 6) +
    theme(axis.title.x=element_blank(), axis.title.y=element_blank(),
          axis.text.x = element_text(family="Arial"),
          axis.text.y = element_text(family="Arial"))
  
  print(p)
  dev.off()
}

Factors <- function(data1, data2, name1, name2, dir.name.bar){
  ##Prepare output
  #file
  csv.name <- str_c(dir.name.bar, name1, sep = "/") %>%
    str_c("vs", name2, ".t-test.result.csv", sep = '')
  
  ###Call parameter names
  param1 <- names(data1)
  param2 <- names(data2)
  
  ###Calculate correlation coefficient
  ##Prepare output
  result.csv <- createEmptyDf(length(param1)*length(param2), 1, colnames = c("pvalue") )
  x.param.names <- vector()
  y.param.names <- vector()
  test.names <- vector()
  ##Iterative process
  k <- 0
  for (i in 1:length(param1)){
    for (j in 1:length(param2)){
      k <- k+1
      #Call parameter name
      x.param_name <- param1[i]
      y.param_name <- param2[j]
      x.param.names <- c(x.param.names, x.param_name)
      y.param.names <- c(y.param.names, y.param_name)
      
      #Create data.frame for analysis
      d_sub <- data.frame(data1[,i], data2[,j])
      names(d_sub) <- c("repertoire", "clinical")
      
      #classify based on clinical parameter
      clin.class <- length(unique(d_sub$clinical))
      #class.num = 2 -> Perform wilcox.test
      if(clin.class == 2){
        ##Perform test
        split.data <- split(d_sub, d_sub$clinical)
        test <- "wilcox.test"
        d_test <- wilcox.test(split.data[[1]]$repertoire, split.data[[2]]$repertoire)
        p.value <- as.numeric(d_test$p.value)
      }
      #class.num > 2 -> Perform Kruskal-Wallis test
      if(clin.class > 2){
        ##Perform test
        test <- "kruskal.test"
        d_test <- kruskal.test(repertoire~clinical, d_sub)
        p.value <- as.numeric(d_test$p.value)
      }
      #Keep records
      result.csv[k,1] <- p.value
      test.names <- c(test.names, test)
      
      ###Create Bar Plot
      dir.name.output <- dir.name.bar
      if(as.numeric(p.value) < 0.05){
        BarPlot(d_sub, dir.name.output, x.param_name, y.param_name, test, p.value)
      }
    } 
  }
  
  ##Add information for correlation coefficient table
  result.csv$test <- test.names
  result.csv$x.param <- x.param.names
  result.csv$y.param <- y.param.names
  #Calculate q-value
  result.csv <- na.omit(result.csv)
  qobj<- qvalue(p = result.csv$pvalue)
  result.csv$qvalue <- qobj$qvalues
  #Judge significance
  result.csv$signif <- ""
  result.csv$signif[which(result.csv$qvalue < 0.05)] <- "*"
  result.csv$signif[which(result.csv$qvalue < 0.01)] <- "**"
  result.csv$signif[which(result.csv$qvalue < 0.001)] <- "***"
  #Output
  write.csv(result.csv, csv.name, row.names = FALSE)
}

######################################################################################
dir.create(dir.name.heatmap, recursive = TRUE)
dir.create(dir.name.scatter, recursive = TRUE)

#Load dataset
for(Cell in Cells){
  #load data
  name.input <- str_c(dir.input, "Merged.Metaanalysis_", sep = "/") %>% str_c(Cell, ".csv")
  d <- tableread_fast(name.input, header = TRUE, sep=",")
  clinical <- read.table(name.clinical, header = TRUE)
  
  #HLA allele保有データへの変換
  HLAs <- as.vector(str_split(clinical$HLA.All, pattern = "_", simplify = TRUE))
  HLAs <- unique(HLAs[HLAs>0])
  table_HLA <- createEmptyDf( nrow(clinical), length(HLAs), colnames = HLAs )
  for(i in 1:length(HLAs)){
    hla <- HLAs[i]
    table_HLA[,i] <- str_detect(clinical$HLA.All, hla)
  }
  HLA_count <- apply(table_HLA, 2, sum)
  table_HLA_analyze <- dplyr::select(table_HLA, names(HLA_count[HLA_count>4]))
  
  #カテゴリーデータと数値データに分類
  clin.fact <- cbind(dplyr::select(clinical, c(param.fact)), table_HLA_analyze)
  clin.num <- dplyr::select(clinical, c(param.num))
  
  ###数値データ: 相関解析実施
  ##Perform analysis
  Correlation(d, clin.num, Cell, "Clinical",
              dir.name.heatmap, dir.name.scatter)
  
  #カテゴリーデータ: グラフ作成, One-way ANOVA実施
  ##Prepare output directory
  dir.create(dir.name.bar, recursive = TRUE)
  ##Perform analysis
  Factors(d, clin.fact, Cell, "Clinical", dir.name.bar)
}
