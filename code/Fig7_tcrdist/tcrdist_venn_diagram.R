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

#################################### Processing layer ####################################################
tableread_fast = function(i, header=TRUE, quote="", sep=","){
  tmp = fread(i, header=header, sep=sep, quote=quote, nThread=32)
  tmp = as.data.frame(tmp)
  return(tmp)
}

#load tcrdist3 data
d <- tableread_fast(name.input, header = TRUE, quote="\"", sep=",")
d$clone.name <- str_c(d$ntvj, d$sample)

#Extract clonotypes
d_NF9_tet <- dplyr::filter(d, tetramer == "NF9")
d_QI9_tet <- dplyr::filter(d, tetramer == "QI9")
d_NF9_dist <- dplyr::filter(d, target == "NYNYLYRLF&A24-02")
d_QI9_dist <- dplyr::filter(d, target == "QYIKWPWYI&A24-02")
d_Spike_dist <- dplyr::filter(d, target != "non&non")
tcrdist_samples <- as.vector(unique(d_Spike_dist$sample))
d_AIM <- dplyr::filter(d, sample %in% tcrdist_samples) %>% dplyr::filter(class == "AIM")


#Create Venn diagram

##FigS7B, NF9 tcrdist vs tetramer
data <- list(tcrdist = d_NF9_dist$clone.name, Tetramer = d_NF9_tet$clone.name)
ppi <- 600
tiff("result/Fig7/S7B_NF9_venn.tiff", width=1.4*ppi, height=1*ppi, res = ppi)
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

##FigS7B, QI9 tcrdist vs tetramer
data <- list(tcrdist = d_QI9_dist$clone.name, Tetramer = d_QI9_tet$clone.name)
ppi <- 600
tiff("result/Fig7/S7B_QI9_venn.tiff", width=1.4*ppi, height=1*ppi, res = ppi)
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

##FigS7C, Atcrdist vs AIM
data <- list(tcrdist = d_Spike_dist$clone.name, AIM = d_AIM$clone.name)
ppi <- 600
tiff("result/Fig7/S7C_AIM-tcrdist_venn.tiff", width=1.4*ppi, height=1*ppi, res = ppi)
p <- ggvenn(data, auto_scale = TRUE,
            fill_color = c("blue", "green"),
            show_percentage = FALSE,
            stroke_size = 0.3,
            text_size = 0,
            show_outside = "none",
            set_name_size = 0
)
print(p)
dev.off()