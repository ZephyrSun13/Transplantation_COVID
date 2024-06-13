
Expr <- read.delim(file = "/sc/arion/projects/zhangw09a/PANDA/db_ZS/Xena/LUSC/TCGA-LUSC.htseq_counts.tsv.sym", row.names = 1, check.names = FALSE)

Expr <- Expr[, grepl("01A", names(Expr))]

names(Expr) <- sapply(names(Expr), function(x) paste(unlist(strsplit(x, split = "-"))[1:3], collapse = "-"))

write.table(Expr, file = "Expr.xls", sep = "\t", quote = FALSE, col.names = NA)

