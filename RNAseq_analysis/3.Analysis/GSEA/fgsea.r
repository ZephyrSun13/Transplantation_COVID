
Args <- commandArgs(trailingOnly = TRUE)

library(fgsea)
library(ggplot2)
library(data.table)
library(RColorBrewer)

plotEnrichmentSun <- function (pathway, stats, gseaParam = 1, ticksSize = 0.2) 
{
    rnk <- rank(-stats)
    ord <- order(rnk)
    statsAdj <- stats[ord]
    statsAdj <- sign(statsAdj) * (abs(statsAdj)^gseaParam)
    statsAdj <- statsAdj/max(abs(statsAdj))
    pathway <- unname(as.vector(na.omit(match(pathway, names(statsAdj)))))
    pathway <- sort(pathway)
    gseaRes <- calcGseaStat(statsAdj, selectedStats = pathway, 
        returnAllExtremes = TRUE)
    bottoms <- gseaRes$bottoms
    tops <- gseaRes$tops
    n <- length(statsAdj)
    xs <- as.vector(rbind(pathway - 1, pathway))
    ys <- as.vector(rbind(bottoms, tops))
    toPlot <- data.frame(x = c(0, xs, n + 1), y = c(0, ys, 0))
    diff <- (max(tops) - min(bottoms))/8
    x = y = NULL
    Colors <- brewer.pal(12, "Set3")
    g <- ggplot(toPlot, aes(x = x, y = y)) + 
	geom_point(color = Colors[1], size = 0.1) + 
	geom_hline(yintercept = max(tops), colour = Colors[4], linetype = "dashed") + 
	geom_hline(yintercept = min(bottoms), colour = Colors[4], linetype = "dashed") + 
	geom_hline(yintercept = 0, colour = "black") + 
	geom_area(fill = Colors[1]) + 
	theme_bw() + 
        geom_segment(data = data.frame(x = pathway), mapping = aes(x = x, y = -diff/2, xend = x, yend = diff/2), size = ticksSize) + 
        theme(panel.border = element_blank(), panel.grid.minor = element_blank(), text = element_text(size = 23)) + 
        labs(x = "rank", y = "enrichment score")
    g
}

#rnk.file <- system.file("extdata", "naive.vs.th1.rnk", package="fgsea")
#gmt.file <- system.file("extdata", "mouse.reactome.gmt", package="fgsea")

ranks <- read.table(Args[1], header=TRUE, colClasses = c("character", "numeric"))

names(ranks) <- c("ID", "t")

ranks <- setNames(ranks$t, ranks$ID)

str(ranks)

pathways <- gmtPathways(Args[2])

str(head(pathways))

fgseaRes <- fgsea(pathways = pathways, stats = ranks, minSize  = 15, maxSize  = 500, nperm = 10000)

saveRDS(fgseaRes, file = "fgseaRes.rds")

head(fgseaRes[order(pval), ])

#for(cc in fgseaRes$pathway[fgseaRes$pval <= as.numeric(Args[3])]){
#
#	pp <- plotEnrichmentSun(pathways[[cc]], ranks) +  ggtitle(paste0(cc, "    P = ", round(fgseaRes[fgseaRes$pathway == cc, "pval"], digits = 3))) + theme(plot.title = element_text(hjust = 0.5, size = 23), axis.text.x = element_text(angle = 45, hjust = 1))
#
#	ggsave(paste0(cc, ".pdf"), pp, height = 4, width = 4)
#
#}

topPathwaysUp <- fgseaRes[ES > 0][head(order(pval), n=10), pathway]

topPathwaysDown <- fgseaRes[ES < 0][head(order(pval), n=10), pathway]

topPathways <- c(topPathwaysUp, rev(topPathwaysDown))

pp <- plotGseaTable(pathways[topPathways], ranks, fgseaRes, gseaParam=0.5)

ggsave("GSEATab.pdf", pp)

collapsedPathways <- collapsePathways(fgseaRes[order(pval)][padj < 0.01], pathways, ranks)

mainPathways <- fgseaRes[pathway %in% collapsedPathways$mainPathways][order(-NES), pathway]

pp <- plotGseaTable(pathways[mainPathways], ranks, fgseaRes, gseaParam = 0.5)

ggsave("GSEATabMain.pdf", pp)

fwrite(fgseaRes, file="fgseaRes.txt", sep="\t", sep2=c("", " ", ""))


#library(org.Mm.eg.db)

#fgseaResMain <- fgseaRes[match(mainPathways, pathway)]

#fgseaResMain[, leadingEdge := mapIdsList(x=org.Mm.eg.db, keys=leadingEdge,keytype="ENTREZID", column="SYMBOL")]

#fwrite(fgseaResMain, file="fgseaResMain.txt", sep="\t", sep2=c("", " ", ""))


# ReactomePathways

#pathways <- reactomePathways(names(exampleRanks))
#
#fgseaRes <- fgsea(pathways, exampleRanks, maxSize=500)
#
#head(fgseaRes)

