
Expr <- read.table(file = "Expr.xls", row.names = 1, header = TRUE, sep = "\t", check.names = F)

Ann <- read.table(file = "/sc/arion/projects/zhangw09a/PANDA/db_ZS/Ensembl/gene_symbol.txt", header = TRUE, sep = "\t", fill = TRUE, stringsAsFactors = FALSE)

Expr <- merge(Ann, Expr, by.x =1 , by.y = 0)

Expr <- Expr[which(Expr$HGNC.symbol != ""), ]

Expr <- Expr[!duplicated(Expr$HGNC.symbol), ]

rownames(Expr) <- Expr$HGNC.symbol

#row.names(Expr) <- make.names(Expr[[2]], unique = TRUE)

Expr <- Expr[, c(-1, -2)]

write.table(Expr, file = "Expr.xls.sym", sep = "\t", col.names = NA, quote = FALSE)

library(edgeR)

dge <- DGEList(counts = Expr)

dge <- calcNormFactors(dge)

Expr.norm <- cpm(dge, log = F, normalized.lib.sizes = F)

write.table(Expr.norm, file = "Expr.xls.sym.cpm.xls", sep = "\t", col.names = NA, quote = FALSE)

