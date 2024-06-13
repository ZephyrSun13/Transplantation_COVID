
library(RColorBrewer)
library(VennDiagram)
library(grDevices)
library(ggplot2)

FileLis <- read.table(file = "EnrichRLis", as.is = T, row.names = 2)

TabLis <- lapply(FileLis[[1]], function(x) read.delim(x, as.is = T))

names(TabLis) <- rownames(FileLis)

TabLis <- lapply(TabLis, function(x) x[x$P.value <= 0.05, ])


DBs <- c("GO_Biological_Process_2021", "KEGG_2021_Human", "COVID-19_Related_Gene_Sets_2021", "WikiPathway_2021_Human")

Over <- lapply(DBs, function(x) Reduce(intersect, lapply(TabLis, function(y) y[y$DB == x, "Term"])))

OverCount <- lapply(DBs, function(x) c(unlist(lapply(TabLis, function(y) nrow(y[y$DB == x, ]))), length(Reduce(intersect, lapply(TabLis, function(y) y[y$DB == x, "Term"])))))

names(OverCount) <- DBs

BgCount <- read.delim(file = "EnrichR.bg.count", row.names = 1, header = F)

lapply(names(OverCount), function(x) c(OverCount[[x]][3], OverCount[[x]][1] - OverCount[[x]][3], OverCount[[x]][2] - OverCount[[x]][3], BgCount[x, ], phyper(OverCount[[x]][3] - 1, OverCount[[x]][1], BgCount[x, ] - OverCount[[x]][1], OverCount[[x]][2], lower.tail = F)))


OverTab <- Reduce(cbind, lapply(TabLis, function(x) x[match(unlist(Over), x$Term), ]))

write.table(OverTab, file = "OverTabEnrichR.xls", sep = "\t", quote = F, row.names = F)


#Over <- unlist(read.delim(file = "EnrichR.over", header = F, as.is = T))

Over <- read.delim(file = "EnrichR.over", header = F, as.is = T)

TabLis <- lapply(names(TabLis), function(x) {TabLis[[x]]$Group = x; TabLis[[x]]})

Tab <- Reduce(rbind, lapply(TabLis, function(x) x[x$Term %in% Over[[1]], ]))

Tab$Term <- Over[match(Tab$Term, Over[[1]]), 2]

Tab$Term <- factor(as.character(Tab$Term), levels = rev(Over[[2]]))

pp <- ggplot(Tab, aes(x = Group, y = Term, color = -log10(P.value))) +

        geom_point(aes(size = Count)) +

        scale_size(range = c(5, 10)) +

        scale_color_gradient(high = brewer.pal(9, "YlOrRd")[9], low = brewer.pal(9, "YlOrRd")[4]) +

	scale_x_discrete(breaks=c("Acute", "Post_acute-Acute"), labels=c("Acute", "Post-acute vs. acute")) +

        theme_minimal() +

        xlab("") +

        ylab("") +

        theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 20, face = "bold"), axis.text.y = element_text(size = 24, face = "bold"), legend.text = element_text(size = 18), legend.title = element_text(size = 18))

ggsave(file = paste0("EnrichROver", "_GOPlot.pdf"), pp, height = 8, width = 15)

