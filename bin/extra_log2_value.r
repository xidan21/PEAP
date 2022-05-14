#!/bin/R



#require("biomaRt")

argv <- commandArgs(trailingOnly=T)


y <- read.csv(argv[1], header = T, sep = "\t")

y[(which(y[,"log2FoldChange"] > 0)),"log2FoldChange"]/max(abs(y[(which(y[,"log2FoldChange"] > 0)),"log2FoldChange"])) -> y[(which(y[,"log2FoldChange"] > 0)),"log2FoldChange"] 
y[(which(y[,"log2FoldChange"] < 0)),"log2FoldChange"]/max(abs(y[(which(y[,"log2FoldChange"] < 0)),"log2FoldChange"])) -> y[(which(y[,"log2FoldChange"] < 0)),"log2FoldChange"] 


annotations <- read.csv(paste0("../lib/", tolower(argv[2]),"/", tolower(argv[2]), "_annotations.txt"), header = T, sep = "\t")

files <- list.files(paste0("../lib/",tolower(argv[2]),"/"))


if (argv[2] == "Human"){ files_index <- grep ("hsa*", files) 
}else if (argv[2] == "Mouse"){ files_index <- grep ("mmu*", files) 
}else if (argv[2] == "Rat"){ files_index <- grep ("mmu*", files) }

pathway_files <- files[files_index]
for (i in 1:length(pathway_files)){

	x <- read.csv(paste0("../lib/",argv[2],"/",pathway_files[i]), header = T, sep = "\t")
	match_index_1 <- match (x[,1], annotations[,2])
	annotations_1 <- annotations[match_index_1[which(!is.na(match_index_1))],]
	print (head(match_index_1))
	match_index <- match(annotations_1[,1], y[,1])
	output <- data.frame(annotations_1[which(!is.na(match_index)),2], annotations_1[which(!is.na(match_index)),1], y[match_index[which(!is.na(match_index))],"log2FoldChange"])

	print(match_index)	
	
	write.table(output, paste("../source/",gsub("\\.txt","",pathway_files[i]),"_entrez_id_with_log2_value.txt", sep=""), row.names = F, col.names = F, quote = F)

}

