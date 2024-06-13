
#R/4.0.3

library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(RColorBrewer)
library(ggraph)
library(GOSemSim)
library(enrichplot)
library(pheatmap)

Args <- commandArgs(trailingOnly = TRUE)

GeneList <- Args[1]

IDType <- Args[2]


## Subfunction

heatplot2 <- function(Obj, showCategory = 100, foldChange=NULL){

	Obj <- head(Obj, showCategory)

	Genes <- lapply(Obj$geneID, function(x) unlist(strsplit(x, split = "/")))

	names(Genes) <- Obj$ID

	GeneAll <- unique(unlist(Genes))

	Matr <- matrix(0, nrow = length(Genes), ncol = length(GeneAll), dimnames = list(names(Genes), GeneAll))

	for(nn in names(Genes)){

		if(length(foldChange) == 0){

			Matr[nn, Genes[[nn]]] <- 1

		}else{

			Matr[nn, Genes[[nn]]]  <- foldChange[Genes[[nn]]]

		}

	}

	rownames(Matr) <- Obj[rownames(Matr), "Description"]

	Matr

}




Colors <- brewer.pal(12, "Set3")

GeneL <- read.delim(file = GeneList, header = FALSE)

GeneL <- GeneL[!duplicated(GeneL[[1]]), ]

GeneL <- GeneL[!is.na(GeneL[[1]]), ]

rownames(GeneL) <- GeneL[[1]]

GeneL <- GeneL[, -1, drop = FALSE]

Tmp <- unlist(GeneL)

names(Tmp) <- rownames(GeneL)

GeneL <- Tmp

Sym <- names(GeneL)

if(IDType == "ENSEMBL"){

	Sym <- mapIds(org.Hs.eg.db, keys = sapply(Sym, function(x) unlist(strsplit(x, split = "\\."))[1]), keytype = "ENSEMBL", column="SYMBOL")

        names(GeneL) <- Sym[names(GeneL)]

}else if(IDType == "Entrez"){

	Sym <- mapIds(org.Hs.eg.db, keys = sapply(Sym, function(x) unlist(strsplit(x, split = "\\."))[1]), keytype = "ENTREZID", column="SYMBOL")

	names(GeneL) <- Sym[names(GeneL)]

}

Sym_map <- bitr(Sym, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

GeneL <- GeneL[Sym_map$SYMBOL]

GeneL.sym <- GeneL

names(GeneL) <- Sym_map$ENTREZID

gg <- enrichGO(gene = Sym_map$ENTREZID, keyType = "ENTREZID", OrgDb = org.Hs.eg.db, ont = "BP", pAdjustMethod = "BH", pvalueCutoff  = 2, qvalueCutoff  = 2, readable = TRUE)

d <- godata('org.Hs.eg.db', ont="BP")

gg2 <- pairwise_termsim(gg, method="Wang", semData = d)


ggx <- setReadable(gg, 'org.Hs.eg.db', 'ENTREZID')

#p1 <- cnetplot(ggx, foldChange=GeneL, showCategory = 5)
#
#ggsave(paste0(GeneList, "_conceptNet.pdf"), p1)


#p1 <- heatplot(ggx, showCategory = 100, foldChange=GeneL) 
#
#p1 <- p1 + theme(axis.text.x = element_blank())
#
#ggsave(paste0(GeneList, "_Heat.pdf"), p1, height = 30, width = 15)


Matr <- heatplot2(ggx, showCategory = 100, foldChange=GeneL.sym)

pdf(paste0(GeneList, "_Heat2.pdf"), width = 12)

	pheatmap(Matr, fontsize_row = 3, fontsize_col = 0.5)

dev.off()


p1 <- emapplot(gg2, layout = "kk", showCategory = 100, cex_label_category=0.5, cex_line = 0.2)

#p1 <- p1 + geom_node_text(aes(label = name), size = 0.5)

#p1 <- p1 + theme(text = element_text(size = 10), cex_label_category=2)

ggsave(paste0(GeneList, "_Over.pdf"), p1)


gg_ori <- gg

gg <- data.frame(gg)

gg <- gg[gg$pvalue < 0.1, ]

write.table(gg_ori, file = paste0(GeneList, "_GO_enrich.xls"), sep = "\t", row.names = FALSE, quote = FALSE)


gg <- gg[1:20, ]

gg$Description <- factor(gg$Description, levels = rev(gg$Description))

p <- ggplot(gg, aes(x=1, y=Description, color=-log10(pvalue))) +

        geom_point(aes(size=Count)) +

        #scale_color_gradient(low=Colors[5], high=Colors[4]) +

	scale_color_gradient(high = brewer.pal(9, "YlOrRd")[9], low = brewer.pal(9, "YlOrRd")[4]) +

	theme_minimal() +

        #scale_color_gradient2(midpoint=mid, low=Colors[5], mid="white", high=Colors[4], space ="Lab" ) +

        theme(text = element_text(size=15), axis.text.x = element_blank(), axis.title.x = element_blank())

        #theme(panel.background = element_blank())

ggsave(file = paste0(GeneList, "_GO_enrich.pdf"), p, width = 9)

