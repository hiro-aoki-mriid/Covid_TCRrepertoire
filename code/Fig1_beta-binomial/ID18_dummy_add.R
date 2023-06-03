###Add dummy data for ID18 TP11
d <- read.table("tmp/result/intermediate/1_beta-binomial/JoinTP/NaraCOVID_CD4_018_TCR.join.strict.table.txt", header = TRUE)
d$NaraCOVID_TP11_CD4_018_TCR.pool.strict.table <- 0
write.table(d, "tmp/result/intermediate/1_beta-binomial/JoinTP/NaraCOVID_CD4_018_TCR.join.strict.table.txt", row.names = FALSE, sep = "\t", quote = FALSE)

d <- read.table("tmp/result/intermediate/1_beta-binomial/JoinTP/NaraCOVID_CD8_018_TCR.join.strict.table.txt", header = TRUE)
d$NaraCOVID_TP11_CD8_018_TCR.pool.strict.table <- 0
write.table(d, "tmp/result/intermediate/1_beta-binomial/JoinTP/NaraCOVID_CD8_018_TCR.join.strict.table.txt", row.names = FALSE, sep = "\t", quote = FALSE)
