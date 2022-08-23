######Covid mRNA vaccine, repertoire analysis
#####Figure 4, 5: Repertoire analysis on Spike reactive CD4 and CD8 T cells determined by AIM assay

###Run the following docker container, and perform analysis in the "tmp" directory
#haokimriid/singlecell
##Usage
#docker run --rm -p 8888:8888 -v /mnt/"Your working directory"/:/work haokimriid/singlecell jupyter notebook --allow-root
#(Ex) docker run --rm -p 8888:8888 -v /mnt/c/Users/Owner/Desktop/COVID_paper/Fig6_Fig7:/work haokimriid/singlecell jupyter notebook --allow-root

##Following analysis are performed on jupyter notebook on docker images
#Preprocessing single-cell data: CovidscTCR_preprocess.ipynb
#Create graph for CD8 dataset: Fig6.CovidscTCR_plotting_CD8.ipynb
#Create graph for CD4 dataset: Fig7.CovidscTCR_plotting_CD4.ipynb