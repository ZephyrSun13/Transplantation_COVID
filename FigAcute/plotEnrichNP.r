
Args <- commandArgs(trailingOnly = T)

library(ggplot2)
library(RColorBrewer)

Enrich <- read.delim(file = Args[1], as.is = T)

DEG <- read.delim(file = Args[2], row.names = 1, as.is = T, check.names = F)

DEG <- DEG[, 1:7]

DEG <- DEG[!duplicated(DEG$Symbol), ]

rownames(DEG) <- DEG$Symbol


Enrich$Pos <- sapply(Enrich$Genes, function(x) {Genes <- unlist(strsplit(x, split = ";")); length(which(DEG[Genes, "Log2Rat"]>0))})

Enrich$Neg <- sapply(Enrich$Genes, function(x) {Genes <- unlist(strsplit(x, split = ";")); length(which(DEG[Genes, "Log2Rat"]<0))})

write.table(Enrich, file = paste0(Args[1], ".np.xls"), sep = "\t", quote = F, row.names = F)

ForPlot <- unlist(read.delim(file = Args[3], as.is = T))

Enrich <- Enrich[Enrich$Term %in% ForPlot, ]

Enrich$Count <- Enrich$Pos + Enrich$Neg

Tab <- rbind(data.frame(Term = Enrich$Term, P = -log10(Enrich$P.value), Rate = -log10(Enrich$P.value)/Enrich$Count*Enrich$Pos, Group = "Pos"), data.frame(Term = Enrich$Term, P = -log10(Enrich$P.value), Rate = -log10(Enrich$P.value)/Enrich$Count*Enrich$Neg, Group = "Neg"))

Tab$Term <- factor(as.character(Tab$Term), levels = rev(ForPlot))

pp <- ggplot(Tab, aes(x = Term, y = Rate, fill = Group)) +

	geom_bar(stat="identity") +

	scale_fill_manual(values = brewer.pal(8, "Set2")[2:3])+

	theme_minimal() +

	ylab("-log(P)") +

	xlab("") +

	coord_flip() 

ggsave(paste0(Args[1], ".plot.pdf"), pp, width = 12, height = 5)

