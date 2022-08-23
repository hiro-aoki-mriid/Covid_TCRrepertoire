#R source codes for generating the graphs in Figure 1, 2, S1, and S2.

#Figure S1A: Calculate total frequency of clones used for WGCNA (top 150 clones,n=41).
source("S1A.Top150Clone.Information.R")

#Figure S1B-F: Generate heatmap of WCGNA modules
#"module_result.CD4.150.csv" and "module_result.CD8.150.csv" are also used as Fig S1D and Fig 1C.
source("S1BCEF_WGCNA.preprocess_rank.R")

#Grouping T-cell clones based on their kinetic profile using WGCNA
source("WGCNA.preprocess.R")

#Figure 1D: Total frequency of CD8+ clones classified into kinetic profile patterns at each time point.
source("1D.WGCNA.summarize.bargraph.R")
#Figure 2A: Total frequency of CD4+ clones classified into kinetic profile patterns at each time point.
source("2A.WGCNA.summarize.bargraph.R")

#Figure 1E and 2B: Comparison of CD8+ clonal response between kinetic profile patterns.
#Output was slso used for FIgure S2 
source("1E2B.WGCNA.summarize.R")
