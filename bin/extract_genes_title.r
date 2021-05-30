#!/bin/R




argv <- commandArgs(trailingOnly=T)



x <- read.csv(argv[1], header = T, sep = "\t")

y <- x[order(x[,"log2FoldChange"], decreasing = T),]


z <- y[which(!is.na(y[,"log2FoldChange"])),]

write.table( z[,1], "../source/ranked_ensembl_ids.txt", col.names =T, row.names =F, quote =F, sep = "\t")



