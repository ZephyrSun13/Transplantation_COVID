
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
        labs(x = "Rank", y = "Enrichment score")
    g
}


fgseaRes <- readRDS(file = Args[1])

ranks <- read.table(Args[3], header=TRUE, colClasses = c("character", "numeric"))

names(ranks) <- c("ID", "t")

ranks <- setNames(ranks$t, ranks$ID)

pathways <- gmtPathways(Args[4])


Plots <- unlist(read.delim(file = Args[2], as.is = T, header = F))

for(cc in Plots){

	pp <- plotEnrichmentSun(pathways[[cc]], ranks) +  ggtitle(paste0(cc, " \n P = ", round(fgseaRes[fgseaRes$pathway == cc, "pval"], digits = 3))) + theme(plot.title = element_text(hjust = 0.5, size = 15), axis.text.x = element_text(angle = 45, hjust = 1))

	ggsave(paste0(cc, ".pdf"), pp, height = 5, width = 7)

}

