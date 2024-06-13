
#Args <- commandArgs(trailingOnly = T)

Args <- c("Expr.sig.tpm.log2.xls", "FactorTmp.xls", "KeyGenes", "Group")

library(ggplot2)
library(RColorBrewer)
library(gridExtra)

Expr <- read.delim(file = Args[1], row.names = 1, as.is = T, check.names = F)

Factor <- read.delim(file = Args[2], row.names = 1, as.is = T)

Genes <- unlist(read.table(file = Args[3], as.is = T))

Factor[[Args[4]]] <- as.factor(Factor[[Args[4]]])

Plots <- list()

for(gg in Genes){

        Factor$Expr <- unlist(Expr[gg, rownames(Factor)])

        pp <- ggplot(Factor, aes_string(x = Args[4], y = "Expr")) +

		geom_boxplot(outlier.shape = NA, fill = brewer.pal(8, "Set2")[3]) +

		geom_jitter(shape = 16, position = position_jitter(0.2)) +

                #geom_point(color = brewer.pal(8, "Set2")[2]) +

                ylab(gg) +

		xlab("Severity Score") +

                #geom_boxplot(outlier.shape = NA) +

                theme_minimal()

        Plots[[gg]] <- pp

}

ggsave(paste0(Args[4], "_TopGene.pdf"), arrangeGrob(grobs = Plots, ncol = 6), height = 2, width = 12)


