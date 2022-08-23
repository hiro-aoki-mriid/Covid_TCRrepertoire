###Perform WGCNA analysis for variation patterns in clone frequency

#Load Libraries
library(foreach)
library(doParallel)
library(dplyr)
library(WGCNA)
library(flashClust)
options(stringsAsFactors=FALSE)
enableWGCNAThreads(14)            #maximum number of CPU threads of your PC
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
library(stringr)
library(data.table)
library(hfunk) #https://github.com/hstojic/hfunk

#Input layer
dir <- "Join"
Tcells <- c("CD4", "CD8") #Whether CD4 or CD8 repertoire is analysed
dir.output <- "S1BCEF"
#Parameters for WGCNA analysis
topXs <- c(150, 200, 250) #Count of clones used for WGCNA for each repertoire
sft_limit <- 30 #upper limit for soft threshold for WGCNA analysis
merge_thres=0.35 #threshold for merging clusters identified by WGCNA
cores <- 14
ppi <- 600
th.correlation <- 0.85 #threshold for clone detection which correlated with WGCNA modules

#################################### Processing layer ####################################################
dir.create(dir.output)

##### Prepare functions #####
#Create empty data.frame
createEmptyDf = function( nrow, ncol, colnames = c() ){
  data.frame( matrix( vector(), nrow, ncol, dimnames = list( c(), colnames ) ) )
}

#Table.read.fast
tableread_fast = function(i, header=TRUE, quote="", sep="\t"){
  tmp = fread(i, header=header, sep=sep, quote=quote, nThread=32)
  tmp = as.data.frame(tmp)
  return(tmp)
}

#stop parallel computing
unregister_dopar <- function() {
  env <- foreach:::.foreachGlobals
  rm(list=ls(name=env), pos=env)
}

#load repertoire data and merge into one data table
Combine <- function(file.name, dir){ 
  #load repertoire data
  name.input <- str_c(dir, file.name, sep = "/")
  d <- tableread_fast(name.input, header = TRUE, quote="\"", sep="\t")
  names(d) <- c("count", "freq",  "cdr3nt", "cdr3aa", "v", "d", "j", "VEnd", "DStart", "DEnd", "JStart",
                "peak", "occurences", "sampling.p","P1", "P2", "P3", "P4")
  
  #Append sample information and clone ID
  name <- str_split(file.name, pattern = "_", simplify = TRUE)
  d$sample <- str_replace(name[3], ".join.strict.table.txt", "")
  d$rank <- c(1:nrow(d))
  return(d)
}

#Determine soft threshold for WGCNA by identifying elbow point of scale free topology 
elbow_finder <- function(thre, sft_limit) { 
  #Modify input
  x_values <- thre[,1]
  y_values <- thre[,2]-min(thre[,2])
  y_values <- abs(y_values - max(y_values))
  max_x_x <- max(x_values)
  max_x_y <- y_values[which.max(x_values)]
  max_y_y <- max(y_values)
  max_y_x <- x_values[which.max(y_values)]
  max_df <- data.frame(x = c(max_y_x, max_x_x), y = c(max_y_y, max_x_y))
  
  # Creating straight line between the max values
  fit <- stats::lm(max_df$y ~ max_df$x)
  # Distance from point to line
  distances <- c()
  for(i in 1:length(x_values)) {
    distances <- c(distances,
                   abs(coef(fit)[2]*x_values[i] - y_values[i] + coef(fit)[1]) / sqrt(coef(fit)[2]^2 + 1^2))
  }
  
  # Max distance point
  pw1_x <- x_values[which.max(distances)]
  
  #If Auto-detect soft threshold exceeds user-supplied limit (sft_limit), use user-supplied value"
  if (pw1_x>=sft_limit){ 
    message(paste("Auto-detect soft threshold is ", pw1_x,
                  " ,exceeds user-supplied limit. use user-supplied value", sep=""))
    pw1_x = sft_limit
  }
  
  #return y value when x = pw1
  pw1_y <- thre[pw1_x,2]
  
  #Set soft threshold to the value that over Scale Free Topology Model Fit > 0.9
  y_near_09 <-1
  y_values <- thre[,2]
  for (i in y_values) {
    if(i>0.9){
      y_near_09 <- i
      break
    }
  }
  x_y_near09 <- thre[,1][which(y_values==y_near_09)]
  if(pw1_y>0.9){
    pw1_x<-x_y_near09
    pw1_y<-y_near_09
  }
  
  return(c(pw1_x, pw1_y))
}

### Plot eigengene value for each module as heatmap
ModuleHeatmap <- function(dir.output, meord_output, topX, Tcell){

  #prepare data
  meord <- dplyr::select(meord_output, -sample)
  names(meord) <- str_remove(names(meord), "ME") %>% str_to_sentence
  temp = t(meord)
  
  #Output
  image.file <- str_c(dir.output, "ModuleEigengene", sep = "/") %>% str_c(Tcell, topX, "tiff", sep = ".")
  #出力ファイルの準備
  ppi <- 600
  png(image.file, width=2.6*ppi, height=1.3*ppi, res=ppi)
  
  p <- Heatmap(temp, name = "Eigengenes", 
               col = colorRamp2(c(-1, -0.5, 0, 0.5, 1), c("blue","cyan","white","magenta", "red")),
               row_names_gp = gpar(fontsize = 8, family="Arial"),
               show_column_names = FALSE,
               cluster_columns = FALSE,
               row_dend_width = unit(0.2, "inch"),
               heatmap_legend_param =list(labels_gp = gpar(fontsize=6, family="Arial"),
                                          title = "Value", 
                                          title_gp = gpar(fontsize=8, family="Arial"),
                                          legend_direction = "vertical",
                                          title_position = "topcenter"))
  q <- draw(p, show_annotation_legend = TRUE,
            heatmap_legend_side = "right")
  #q
  print(q)
  
  invisible(dev.off())
}


###Iterate for Tcells (CD4/CD8) or TopXs (150/200/250)
for(Tcell in Tcells){
  ##### load repertoire data and merge into one data table #####
  files  <- list.files(dir, pattern=Tcell)
  cl <- makeCluster(cores)
  registerDoParallel(cl)   
  out <- foreach(file.name = files, .combine = "rbind",
                 .packages=c("stringr", "data.table")) %dopar% {Combine(file.name, dir)}
  stopCluster(cl)
  out <- as.data.frame(out)
  unregister_dopar()

  ###### Perform WGCNA analysis for variation pattern of clone frequency ######
  for(topX in topXs){
    #Extract clones and prepare data for WGCNA analysis
    #clones detected at all 4 points and top 250 in their averaged frequency
    data <- dplyr::filter(out, occurences >= 4 & rank <= topX) 
    data_ex <- dplyr::select(data, c("P1", "P2", "P3", "P4"))
    data_ex <- t(as.matrix(data_ex))
    
    ##### Determine soft threshold to fulfill scalefreeness of the network #####
    powerseq=c(seq(1,50,by=1))   # sequence of the values that you want to test
    sft<- pickSoftThreshold(data_ex, powerVector=powerseq, verbose=5)
    thre=data.frame(power=sft$fitIndices[,1], fit=-sign(sft$fitIndices[,3])*sft$fitIndices[,2])
    thre1 = elbow_finder(thre, sft_limit)
    pw1_x=thre1[1]
    pw1_y=thre1[2]
    
    ##### Detect variation motif based on the changes in frequency of clone #####
    ##Generate adjacent matrix
    adj1=adjacency(data_ex,power=pw1_x)
    ##Calculate connectivity for each node
    convec=softConnectivity(data_ex,
                            corFnc = "cor", corOptions = "use = 'p'",type = "unsigned",
                            power =pw1_x,blockSize = 1500,
                            minNSamples = NULL,verbose = 2, indent = 0)
    #Convert into TOM (Topological Overlap  matrix)
    stom1=TOMsimilarity(adj1)
    
    #Clustering of genes
    dissim=1-stom1      # Distance matrix for clustering
    genetree=flashClust(as.dist(dissim), method="average")       # you can also use hclust 
    
    #initial detection of modules
    dmd1=cutreeDynamic(dendro=genetree, distM=dissim,
                       deepSplit=2, pamRespectsDendro=FALSE,minClusterSize=30)
    dc1=labels2colors(dmd1)
    melist=WGCNA::moduleEigengenes(data_ex, colors=dc1)
    
    ##Merge ckusters from initial clustering
    #Calculate eigengene for each module
    #Clustering modules based on the correlation between eigengenes
    mes=melist$eigengenes
    mediss=1-WGCNA::cor(mes)
    metree=flashClust::flashClust(as.dist(mediss), method="average")
    merge=WGCNA::mergeCloseModules(data_ex,dc1,cutHeight=merge_thres,
                                   relabel=TRUE, verbose=3)
    ##Reduce the number of clusters
    #If the count of modules is over 8, reduce modules by increasing the threshold for merging
    while (length(unique(merge$colors))>=8){
      message(paste("more than 15 modules are identified, reset merge_thres as ", merge_thres, ".", sep=""))
      merge_thres = merge_thres+0.01
      merge=WGCNA::mergeCloseModules(data_ex,dc1,cutHeight=merge_thres,
                                     relabel=TRUE, verbose=2)
      if (merge_thres >= 0.4) break
    }
    
    ##Determine the module names
    colorOrder=c("grey", standardColors(50))
    modcol=merge$colors
    modlab=match(modcol, colorOrder)-1
    modme=merge$newMEs
    
    ##### Return the results of WGCNA analysis #####
    ##### Output eigengenes for each module #####
    me1=moduleEigengenes(data_ex,modcol)$eigengenes
    meord=orderMEs(me1)
    meord_output <- meord
    meord_output$sample=rownames(data_ex)
    name.out <- str_c(dir.output, "module_result", sep = "/") %>% str_c(Tcell, topX, "csv", sep = ".")
    write.csv(meord_output, name.out, row.names = F)
    
    ##Plot eigengene value for each module as heatmap
    ModuleHeatmap(dir.output, meord_output, topX, Tcell)
  }
}
