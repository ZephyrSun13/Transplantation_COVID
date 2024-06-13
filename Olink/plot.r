
library(pheatmap)
library(ggplot2)
library(RColorBrewer)
library(ggrepel)

Colors <- read.delim(file = "/sc/arion/projects/zhangw09a/PANDA/ext_ZS/bin/code/ColorLibrary", row.names = 1, as.is = T, header = F)

ThresP <- 0.05

DEGs <- read.delim(file = "DEG_diff.xls", row.names = 1, as.is = T)


Tab <- DEGs

Tab$Symbol <- rownames(Tab)

Tab$Group <- "No"

Tab$Group[Tab$Log2Rat > 0 & Tab$limma.p <= ThresP] <- "Up"

Tab$Group[Tab$Log2Rat < 0 & Tab$limma.p <= ThresP] <- "Down"


TabSig <- Tab[Tab$limma.p <= ThresP, ]

#TabSig <- TabSig[c(head(order(TabSig$Log2Rat, decreasing = T), 20), head(order(TabSig$Log2Rat, decreasing = F), 20)), ]

pp <- ggplot(Tab, aes(y=-log10(limma.p), x=Log2Rat, color = Group)) +

        geom_point(alpha = 1) +

        scale_color_manual(values = brewer.pal(8, "Set2")[c(3, 8, 2)]) +

        geom_text_repel(data = TabSig, aes(label = Symbol), size = 3) +

        xlab("log2(Fold change)") +

        ylab("-log10(P value)") +

        geom_hline(yintercept=-log10(ThresP), linetype="dashed", size = 0.5) +

	ggtitle("Severity") +

        theme_classic() +

        theme(legend.position = "none", plot.title = element_text(hjust = 0.5))

ggsave(paste0("Severity", "_Volcano.pdf"), width = 3, height = 2.5)


DEGs <- DEGs[DEGs$limma.p <= 0.05, ]

DEGs <- DEGs[order(DEGs$Log2Rat, decreasing = T), ]

Expr <- read.delim(file = "../Expr.xls", row.names = 1, as.is = T)

Info <- read.delim(file = "../Info.xls", row.names = 1, as.is = T)

Expr.z <- t(apply(Expr, 1, function(x) (x - mean(x, na.rm = T))/sd(x, na.rm = T)))

DEGs.P <- DEGs[, c("Log2Rat", "limma.p"), drop = F]

DEGs.P$limma.p <- -log10(DEGs.P$limma.p)

colnames(DEGs.P) <- c("Log2FC", "Pvalue")

Expr.z <- Expr.z[rownames(DEGs), ]

Expr.z <- rbind(Expr.z, Mean = colMeans(Expr.z), Score = Info[colnames(Expr.z), "SeverityScore"])

Expr.z <- Expr.z[, order(Expr.z["Score", ], Expr.z["Mean", ], decreasing = F)]

pdf("Heatmap.sig.pdf", height = 6, width = 7)

        #pheatmap(Expr.z[rownames(DEGs), rownames(Info)[order(Info$SeverityScore, decreasing = T)]], annotation_col = Info, annotation_row = DEGs.P, cluster_cols = F, cluster_rows = F, fontsize = 5)

        pheatmap(Expr.z[rownames(DEGs), ], annotation_col = Info[, "SeverityScore", drop = F], annotation_row = DEGs.P, cluster_cols = F, cluster_rows = F, fontsize_row = 10, fontsize_col = 3, show_colnames = F, color = colorRampPalette(Colors[paste0("HeatRd", 9:1), ])(100))

dev.off()

