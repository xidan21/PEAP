#!/bin/R




argv <- commandArgs(trailingOnly=T)
x <- read.csv(argv[1], header=F, sep = "\t") # gsea.txt

y <- na.omit(x)

a <- data.frame(y[,1],p.adjust(y[,2], method = "BH", length(y[,2])))

match_index <- match(x[,1], a[,1])

b <- data.frame(x[,1:2], a[match_index,2], x[,3])

write.table(b, sub(".txt","_with_fdr.txt",argv[1]), col.names = F, row.names = F, quote = F, sep = "\t")

annotations <- read.csv(paste0("../lib/", tolower(argv[2]),"/annotation_pathway.txt"), header = F, sep = "\t")

match_index_2 <- match(sub("path:","", annotations[,1]), b[,1])

c <- data.frame(annotations[which(!is.na(match_index_2)),], b[match_index_2[which(!is.na(match_index_2))],c(2:4)])

web_links <- read.csv(sub("gsea.txt", "web_link.txt", argv[1]), header = F, sep ="\t")

match_index_3 <- match(sub("path:","", c[,1]), web_links[,1])

d <- data.frame(c[which(!is.na(match_index_3)),], web_links[match_index_3[which(!is.na(match_index_3))],])

colnames(d) <- c("Pathway ids", "Descriptions", "P values", "FDR", "MES", "Ids", "Web links")

write.table(d, sub("gsea.txt", "pathway_analysis.txt", argv[1]), col.names = T, row.names = F, quote = F, sep = "\t")
