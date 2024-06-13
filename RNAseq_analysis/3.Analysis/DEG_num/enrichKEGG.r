
#!/usr/bin/R

library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(RColorBrewer)
#library(enrichplot)

Args <- commandArgs(trailingOnly = TRUE)

GeneList <- Args[1]

IDType <- Args[2]

Colors <- brewer.pal(12, "Set3")

Sym <- read.table(file = GeneList, sep = "\t", stringsAsFactors = FALSE, header = FALSE)

Sym <- Sym[[1]]

if(IDType == "ENSEMBL"){

	Sym <- mapIds(org.Hs.eg.db, keys = sapply(Sym, function(x) unlist(strsplit(x, split = "\\."))[1]), keytype = "ENSEMBL", column="SYMBOL")

}

Sym_map <- bitr(Sym, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

#write.table(Sym_map, file = paste0(Args[1], "_Sym_map.xls"), sep = "\t", quote = FALSE, col.names = NA)

#kk <- enrichKEGG(gene = Sym_map$ENTREZID, organism = 'hsa', keyType = "ncbi-geneid", pvalueCutoff = 2, qvalueCutoff = 2, universe = all_genes$ENTREZID)


kk <- enrichKEGG(gene = Sym_map$ENTREZID, organism = 'hsa', pvalueCutoff = 2, qvalueCutoff = 2)

#kk_ori <- kk

kk <- data.frame(kk)

#kk <- kk[kk$pvalue < 0.1, ]

kk$geneID <- sapply(kk$geneID, function(x) {

	IDs <- unlist(strsplit(x, split = "/"))

	paste(Sym_map[match(IDs, Sym_map$ENTREZID), "SYMBOL"], collapse = "/")

})

write.table(kk, file = paste0(GeneList, "_KEGG_enrich.xls"), sep = "\t", row.names = FALSE, quote = FALSE)

#p1 <- dotplot(kk, showCategory=30)

if(nrow(kk) >= 20){

	kk <- kk[1:20, ]

}

kk$Description <- factor(kk$Description, levels = rev(kk$Description))

p <- ggplot(kk, aes(x=1, y=Description, color=-log10(pvalue))) + 

	geom_point(aes(size=Count)) +

	#scale_color_gradient(low=Colors[5], high=Colors[4]) +

	#scale_color_gradient2(midpoint=mid, low=Colors[5], mid="white", high=Colors[4], space ="Lab" ) +

        scale_color_gradient(high = brewer.pal(9, "YlOrRd")[9], low = brewer.pal(9, "YlOrRd")[4]) +

	theme_minimal() +

	theme(text = element_text(size=15), axis.text.x = element_blank(), axis.title.x = element_blank())

	#theme(panel.background = element_blank())

ggsave(file = paste0(GeneList, "_KEGG_enrich.pdf"), p, width = 9)

