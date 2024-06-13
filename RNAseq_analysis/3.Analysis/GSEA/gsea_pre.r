
Res <- read.delim(file = "../Group_Group_Score_Age_Sex_Batch_Lymph_Group_diff.ann.xls", row.names = 1)

length(which(duplicated(Res$Symbol)))

Res <- Res[!duplicated(Res$Symbol), ]

Res <- data.frame(Gene = Res$Symbol, T = Res$t)

Res <- Res[order(Res$T, decreasing = TRUE), ]

write.table(Res, file = "Ranks", row.names = FALSE, sep = "\t", quote = FALSE)

