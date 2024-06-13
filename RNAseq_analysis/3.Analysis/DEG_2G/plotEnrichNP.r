
Args <- commandArgs(trailingOnly = T)

library(ggplot2)
library(RColorBrewer)

Colors <- read.delim(file = "/sc/arion/projects/zhangw09a/PANDA/ext_ZS/bin/code/ColorLibrary", row.names = 1, header = F, as.is = T)

Enrich <- read.delim(file = Args[1], as.is = T)

#DEG <- read.delim(file = Args[2], row.names = 1, as.is = T, check.names = F)

#DEG <- DEG[, 1:7]

#DEG <- DEG[!duplicated(DEG$Symbol), ]

#rownames(DEG) <- DEG$Symbol


#Enrich$Pos <- sapply(Enrich$Genes, function(x) {Genes <- unlist(strsplit(x, split = ";")); length(which(DEG[Genes, "Log2Rat"]>0))})

#Enrich$Neg <- sapply(Enrich$Genes, function(x) {Genes <- unlist(strsplit(x, split = ";")); length(which(DEG[Genes, "Log2Rat"]<0))})

#write.table(Enrich, file = paste0(Args[1], ".np.xls"), sep = "\t", quote = F, row.names = F)

#ForPlot <- unlist(read.delim(file = Args[2], as.is = T, header = F))

ForPlot <- read.delim(file = Args[2], as.is = T, header = F)

ForPlot <- ForPlot[ForPlot[[1]] %in% Enrich$Term, ]

Enrich <- Enrich[Enrich$Term %in% ForPlot[[1]], ]

#Enrich$Count <- Enrich$Pos + Enrich$Neg

Tab <- rbind(data.frame(Term = ForPlot[match(Enrich$Term, ForPlot[[1]]), 2], P = -log10(Enrich$P.value), Rate = -log10(Enrich$P.value)/Enrich$Count*Enrich$Pos, Group = "Pos"), data.frame(Term = ForPlot[match(Enrich$Term, ForPlot[[1]]), 2], P = -log10(Enrich$P.value), Rate = -log10(Enrich$P.value)/Enrich$Count*Enrich$Neg, Group = "Neg"))

Tab$Term <- factor(as.character(Tab$Term), levels = rev(ForPlot[[2]]))

Tab$Group <- factor(Tab$Group, levels = c("Neg", "Pos"))

pp <- ggplot(Tab, aes(x = Term, y = Rate, fill = Group)) +

	geom_bar(stat="identity") +

	#scale_fill_manual(values = brewer.pal(8, "Set2")[c(3, 2)]) +

	scale_fill_manual(values = Colors[c("Blue", "Red"), ]) +

	theme_minimal() +

	ylab("-log10(P value)") +

	xlab("") +

	coord_flip() +

	theme(legend.title = element_blank(), axis.text.y = element_text(face = "bold"))

ggsave(paste0(Args[1], ".plot.pdf"), pp, width = 8, height = 4)

