#R source codes for generating the graphs in Figure 4, 5, S3, and S4.

#Figure 4C: Extract AIM+ clones and generate collapsed plot for AIM posi-nega overlap
#Collapsed plots are also used in Fig S3B, C
source("4C.AIMposi_extract.R")

#Figure 4D: Calculate frequency of AIM+ clones at each timepoint
#Also used for Figure S3D and S3E
source("4D.AIMpositivity_timepoint.R")

###Merge information about WGCNA and AIM
source("Combine_AIM_WGCNA.R")

#Figure 4E-H: Summarize AIM positivity within each WGCNA module For CD8 T cells
#Also used for Figure S4A-D
source("4E-H.CD8.AIM_WGCNA_summarize.R")

#Figure 4I, 5D: Extratract Top10 dominant clones among AIM+ clones at P2, P3 and P4
source("4I_5D.AIMposi_Major_extract.R")

#Figure 5A-C: Summarize AIM positivity within each WGCNA module For CD4 T cells
#Also used for Figure S4E-H
source("5A-C.CD4.AIM_WGCNA_summarize.R")

