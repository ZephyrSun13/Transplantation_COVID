
#Factor <- read.table(file = "../FactorAll2Batch.txt", sep = "\t", header = TRUE, row.names = 1, as.is = T)

Factor <- read.table(file = "../Factor.txt", sep = "\t", header = TRUE, row.names = 1, as.is = T)

Factor$Sex <- factor(as.character(Factor$Sex), levels = c("M", "F"))

Factor <- Factor[Factor$AcutevsChronic == "Acute", ]

Factor$Group <- factor(Factor$SeverityScore)

#levels(Factor$Group) <- c(1, 1, 2, 2, 3, 3)

table(Factor$Group)

#Factor$Group[Factor$SeverityScore == 7 | Factor$SeverityScore == 6] <- "S67"

write.table(Factor, file = "Factor.xls", sep = "\t", quote = F, col.names = NA)

Stat <- list()

Stat[["Age"]] <- tapply(Factor$Age, Factor$Group, function(x) paste0(round(mean(x), 2), "+-", round(sd(x), 2)))

Stat[["Lymph"]] <- tapply(Factor$Lymph, Factor$Group, function(x) paste0(round(mean(x, na.rm = T), 2), "+-", round(sd(x, na.rm = T), 2)))

Stat[["Sex"]] <- tapply(Factor$Sex, Factor$Group, function(x) paste(table(x), collapse = "/"))

Stat[["Batch"]] <- tapply(Factor$Batch, Factor$Group, function(x) paste(table(x), collapse = "/"))

Stat <- Reduce(cbind, Stat)

colnames(Stat) <- c("Age", "Lymph", "Sex", "Batch")

write.table(Stat, file = "StatClin.xls", sep = "\t", quote = F, col.names = NA)

