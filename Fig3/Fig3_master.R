#R source codes for generating the graphs in Figure 3.

####Merge metaclonotype information from tcrdist3 to result of WGCNA analysis
source("Combine_tcrdist_WGCNA.R")

#Figure 3B: Total frequency of CD8+ SARS-CoV-2 reactive meta-clonotypes at each time point.
source("3B.tcrdist_positivity_timepoint.R")

#Figure 3C: Summarize tcrdist Metaclonotype pattern within each WGCNA module
source("3C.tcrdist_WGCNA_summarize.R")
