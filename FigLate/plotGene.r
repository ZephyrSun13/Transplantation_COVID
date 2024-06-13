
#Args <- commandArgs(trailingOnly = T)

Args <- c("Expr.sig.tpm.log2.xls", "FactorTmp.xls", "KeyGenes", "Group")

#Args <- c("Expr.sig.tpm.log2.xls", "FactorTmp.xls", "KeyGenes.LILR", "Group")

library(ggplot2)
library(RColorBrewer)
library(gridExtra)

Expr <- read.delim(file = Args[1], row.names = 1, as.is = T, check.names = F)

Factor <- read.delim(file = Args[2], row.names = 1, as.is = T)

Genes <- unlist(read.table(file = Args[3], as.is = T))

Factor[[Args[4]]] <- as.factor(Factor[[Args[4]]])

levels(Factor[[Args[4]]]) <- c("Acute", "Post_Acute")

Plots <- list()

for(gg in Genes){

        Factor$Expr <- unlist(Expr[gg, rownames(Factor)])

        pp <- ggplot(Factor, aes_string(x = Args[4], y = "Expr", fill = Args[4])) +

		geom_boxplot(outlier.shape = NA) +

		scale_fill_manual(values = brewer.pal(8, "Set2")[c(3, 3)]) +

		geom_jitter(shape = 16, position = position_jitter(0.2)) +

	        scale_x_discrete(breaks=c("Acute", "Post_Acute"), labels=c("Acute", "Post-acute")) +

                #geom_point(color = brewer.pal(8, "Set2")[2]) +

                ylab(gg) +

		xlab("") +

                #geom_boxplot(outlier.shape = NA) +

                theme_minimal() +

		theme(legend.position = "none") +

		theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "bold"), axis.title.y = element_text(face = "bold"))

        Plots[[gg]] <- pp

}

ggsave(paste0(Args[4], "_TopGene.pdf"), arrangeGrob(grobs = Plots, ncol = 6), height = 3, width = 14)


