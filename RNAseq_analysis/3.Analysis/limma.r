
args <- commandArgs(trailingOnly = T)

library(limma)

library(multtest)

Expr <- read.table(file = args[1], sep = "\t", header = T, row.names = 1)

dim(Expr)

Factor <- read.table(file = args[2], header = T, stringsAsFactors = T)

row.names(Factor) <- paste0(Factor$Sample, "_count.txt_norm")

Factor[["Sample"]] <- paste0(Factor$Sample, "_count.txt_norm")

type <- args[3]

Focus <- args[4]

Cases <- unlist(strsplit(args[5], split = "-"))

#Factor <- Factor[Factor[[Focus]] %in% Cases, ]

Expr <- Expr[, intersect(rownames(Factor), names(Expr))]

Factor <- Factor[intersect(rownames(Factor), names(Expr)), ]

#if(args[6]){

#	Diff <- read.table(file = args[6], sep = "\t", header = T)

#}

dim(Expr)

Expr[1:5, ]

#table(Factor[[Focus]])

#Cases <- unlist(strsplit(args[5], split = "-"))

print("Starting limma")

if(type=="paired"){

	pair <- factor(Factor[[Focus]], levels = c(Cases[2], Cases[1]))

	stage <- Factor$Identity

	#stage <- as.factor(substr(Factor$Sample, 1, 10))

	design<-model.matrix(~pair+stage)

	fit<-lmFit(Expr, design)

	fit <- eBayes(fit)

        result <- topTable(fit, coef=grep("pair", colnames(design), value = TRUE), adjust="BH", number = 20000)

        result <- result[order(result$adj.P.Val), ]

        names(result) <- c("Log2Rat", "AveExpr", "t", "limma.p", "limma.padj", "B")

        write.table(result, file = paste(args[4], args[5], "paired_diff.all.xls", sep="_"), sep = "\t", quote = F, col.names = NA)

	result <- result[abs(result$Log2Rat) >= as.numeric(args[6]) & result$limma.padj <= as.numeric(args[7]), ]

	#Ann <- read.table("/sc/orga/projects/zhangw09a/PANDA/db_ZS/Ensembl/Ann_ensembl_refseq_gene.xls", sep = "\t", header = TRUE, stringsAsFactors = FALSE, fill = TRUE, comment.char = "")

	Ann <- read.delim(file = "/sc/orga/projects/zhangw09a/PANDA/db_ZS/Ensembl/Ann_ensembl_refseq_gene.xls", stringsAsFactors = F)

	result <- merge(result, Ann, by.x = 0, by.y = 2)

	result <- merge(result, Expr, by.x = 1, by.y = 0)

        write.table(result, file = paste(args[4], args[5], "paired_diff.ann.xls", sep="_"), sep = "\t", quote = F, row.names = FALSE)

	#Result <- topTable(fit, sort.by="p")

	#multadjust <- mt.rawp2adjp(fit$p.value[, grep("stage", colnames(fit$p.value))], proc=c("BH"))

	#eBpvalues <- multadjust$adjp[order(multadjust$index),]

}else{

	design <- model.matrix(~ -1+Factor[[Focus]])

	colnames(design) <- levels(Factor[[Focus]])

	fit <-lmFit(Expr, design)

	if(length(Cases) == 2){

		Contrast <- args[5]

	}else{

		Contrast <- c(paste0(Cases[1], "-", Cases[2]), paste0(Cases[1], "-", Cases[3]), paste0(Cases[2], "-", Cases[3]))

	}

	contrast.matrix <- makeContrasts(contrasts=Contrast, levels=design)

	fit <- contrasts.fit(fit, contrast.matrix)

	fit <- eBayes(fit)

	for(i in 1:length(Contrast)){

		result <- topTable(fit, coef=i, adjust="BH", number = 20000)

		result <- result[order(result$adj.P.Val), ]

		names(result) <- c("Log2Rat", "AveExpr", "t", "limma.p", "limma.padj", "B")		

		write.table(result, file = paste(args[4], args[5], Contrast[i], "diff.xls", sep="_"), sep = "\t", quote = F, col.names = NA)

	        result <- result[abs(result$Log2Rat) >= as.numeric(args[6]) & result$limma.p <= as.numeric(args[7]), ]

	        Ann <- read.delim(file = "/sc/orga/projects/zhangw09a/PANDA/db_ZS/Ensembl/Ann_ensembl_refseq_gene.xls", stringsAsFactors = F)

	        result <- merge(result, Ann, by.x = 0, by.y = 2)

	        result <- merge(result, Expr, by.x = 1, by.y = 0)

                result <- result[order(abs(result$Log2Rat), decreasing = TRUE), ]

       		write.table(result, file = paste(args[4], args[5], Contrast[i], "diff.ann.xls", sep="_"), sep = "\t", quote = F, row.names = FALSE)

	}

}

	#design <- model.matrix(~ 0 + Gender + Age + Cigarette + Alcohol + Family + Location, data = Factor)

#	design <- model.matrix(~ Gender + Age + Cigarette + Alcohol + Family + Location, data = Factor)
#
#	fit <-lmFit(Expr, design)
#
#	fit2 <- eBayes(fit)
#
#	 <- apply(fit2$p.value, 2, function(x) rownames(fit2$p.value)[x<0.01])
#
#	tt[[1]] <- NULL

