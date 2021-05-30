#!/bin/bash






#require("biomaRt")
argv <- commandArgs(trailingOnly=T)

w <- read.csv(argv[1], header = T, sep = "\t")
#w <- read.csv("scr_vs_i_keep_dup_comparsion_ranked_ensembl_ids.txt", header = T, sep = "\t")


# nPerm have to be less than the number of total inputList

gsea <- function(inputSet, inputList, nPerm, setName){

	######################### clean input NA ##########################

	geneSet <- inputSet[which(!is.na (inputSet))]
	geneList <- inputList[which(!is.na (inputList))]	

	per_matrix <- data.frame(geneList,replicate(nPerm-1, sample(geneList, length(geneList))))
	
	######################### clean index_matrix NA ##########################

	rank_geneSet_match_index <- match(geneSet, geneList)
	
	if (sum(is.na(rank_geneSet_match_index)) == length(geneSet))
	{
		print (paste0("there is no gene matched !!!","/",setName))
	
	}else if (sum(!is.na(rank_geneSet_match_index)) == 1)
	{
	
		print (paste0("there is only one gene matched !!!","/", setName))

	}else{
	
#		print (paste("there are ", eval(sum(!is.na(rank_geneSet_match_index))), " of ", length(rank_geneSet_match_index)," gene matched.", sep = ""))
		rank_geneSet <- data.frame(geneSet, rank_geneSet_match_index)












		#print(rank_geneSet)

		index_matrix_primitive <- data.frame(apply(per_matrix, 2, function(per_matrix) {match (geneSet, per_matrix)}))
		index_matrix <- apply(index_matrix_primitive, 2, function(index_matrix_primitive) {index_matrix_primitive [which(!is.na(index_matrix_primitive))]})


		######################### permutation by matrix ##########################

		decrease_value <- -((nrow(index_matrix) /(length(geneList)-nrow(index_matrix)))^0.5)

		increase_value <- ((length(geneList)- nrow(index_matrix))/nrow(index_matrix))^0.5

		######################### input increase/decrease values  ##########################
	
		gsea_matrix <- matrix(decrease_value, nrow=length(geneList), ncol=nPerm)

		for (i in 1: ncol(index_matrix))
		{
			gsea_matrix[index_matrix[,i],i] <- increase_value
		}

		######################### extra max cumsum both positive and negative ##########################

		max_cumsum_index <- apply(gsea_matrix, 2, function(gsea_matrix){which(abs(cumsum(gsea_matrix)) == max(abs(cumsum(gsea_matrix))))})
	
		max_cumsum <- c()
		for (i in 1: ncol(gsea_matrix))
		{
		
			max_cumsum <- c(max_cumsum, (cumsum(gsea_matrix[,i]))[as.numeric(unlist(max_cumsum_index[i]))])
		}
	
		######################### calculate p value of permutation test ##########################

		gsea_score <- max_cumsum[1]
		MES <- (nrow(index_matrix))*increase_value
	
		if (gsea_score >0)		
		{

			pvalue = sum(max_cumsum >= gsea_score) / length(max_cumsum)
	
		}	else
		{
			pvalue = sum(max_cumsum <= gsea_score) / length(max_cumsum)
		}


		######################### plot enrichment bar ##########################

		total_rank = length(geneList)

		pdf(file= paste("../source/",setName, '.pdf', sep=""))
		layout(matrix(c(1,1,1,1,2,2), nrow = 6, ncol = 1, byrow = TRUE))

		plot(1:total_rank, cumsum(gsea_matrix[,1]),  main = setName, type = "l", col = "green", xaxt='n', ann=FALSE,bty="n")
		title(main = setName); lines(1:total_rank,replicate(total_rank,0), col = "red")
		plot(1:total_rank,replicate(total_rank,1),  type = "l",bty="n", axes=FALSE, xaxt='n', ann=FALSE, col = "white")
		segments(c(1:total_rank),1.0,c(1:total_rank),5.5, col=2)

		query_index <- match(geneSet, geneList)


		query_index_clean <- query_index[which(!is.na(query_index))]
		segments(query_index_clean,1.0,query_index_clean,5.5, col=3)

		dev.off()
		######################### print results ##########################

		output <- c(pvalue, gsea_score/MES)
		sprintf("the pvalue is %f", pvalue)

		sprintf ("the enrichment socres of gene sets is %f", gsea_score/MES)

		######################### return results ##########################

		bla <- paste(eval(pvalue), eval(gsea_score/MES), sep="/")

		print (paste("there are ", eval(sum(!is.na(rank_geneSet_match_index))), " of ", length(rank_geneSet_match_index)," gene matched.", "\t", setName,"/",eval(bla),sep = ""))
	}
	
}




















pathways <- read.csv(paste0("../lib/",tolower(argv[2]),"/annotation_pathway_formatted.txt"), header = F, sep = "\t")

annotations <- read.csv(paste0("../lib/", tolower(argv[2]),"/", tolower(argv[2]), "_annotations.txt"), header = T, sep = "\t")

match_index <- match(w[,1], annotations[,1])

annotations_2 <- annotations[match_index[which(!is.na(match_index))],2]



files <- list.files(paste0("../lib/",tolower(argv[2]),"/"))

if (argv[2] == "Human"){ files_index <- grep ("hsa*", files) 
}else if (argv[2] == "Mouse"){ files_index <- grep ("mmu*", files)
}else if (argv[2] == "Rat"){ files_index <- grep ("mmu*", files) }

pathway_files <- files[files_index]
for (i in 1:length(pathway_files))
{
	x <- read.csv(paste0("../lib/",argv[2],"/",pathway_files[i]), header = T, sep = "\t")
	inputSet = x[,1]	
	inputList = annotations_2

	setName = pathways[i,1]
	
	nPerm = 1000

	gsea(inputSet, inputList, nPerm, setName)	
}














