
library(limma)

ExprFile <- "Expr.xls"

FactorFile <- "Info.xls"

Focus <- "SeverityScore"

Expr <- read.table(file = ExprFile, sep = "\t", header = T, row.names = 1, check.names = FALSE)

dim(Expr)

Factor <- read.delim(file = FactorFile, header = T, stringsAsFactors = T, row.names = 1)


                design <- model.matrix(~ Factor[[Focus]])

                #colnames(design) <- levels(Factor[[Focus]])

                colnames(design)[2] <- Focus

        fit <- lmFit(Expr, design)

        fit <- eBayes(fit, trend = T)

                result <- topTable(fit, coef=Focus, adjust="BH", number = 20000)

                result <- result[order(result$P.Val), ]

                names(result) <- c("Log2Rat", "AveExpr", "t", "limma.p", "limma.padj", "B")        
                write.table(result, file = paste0("DEG", "_diff.xls"), sep = "\t", quote = F, col.names = NA)

