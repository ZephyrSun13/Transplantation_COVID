
multiCat <- function(tt){

        Tmp <- table(tt)

        Resu <- paste0(Tmp, "/", sum(Tmp), " (", round(Tmp/sum(Tmp)*100, 1), ")")

        names(Resu) <- names(Tmp)

        Resu

}

TabAll <- read.delim(file = "Demography_merged.full_MM.txt", row.names = 1, as.is = T)

Variates <- c("FKT", "MMFBefore", "Steroids")

Tab <- TabAll[TabAll$Condition == "Acute", ]

Factor <- Tab$SeverityScore

StatMulti <- Reduce(rbind, lapply(Variates, function(x) {cbind(multiCat(factor(Tab[[x]])), Reduce(cbind, tapply(factor(Tab[[x]]), Factor, multiCat)), round(fisher.test(Reduce(cbind, tapply(factor(Tab[[x]]), Factor, table)))$p.value, 2))}))

colnames(StatMulti) <- c("All", "1", "3", "4", "5", "6", "7", "P")

write.table(StatMulti, file = "ImmuSupr_severity.acute.xls", sep = "\t", quote = F, col.names = NA)


Tab <- TabAll[TabAll$Condition == "Chronic", ]

Factor <- Tab$SeverityScore

StatMulti <- Reduce(rbind, lapply(Variates, function(x) {cbind(multiCat(factor(Tab[[x]])), Reduce(cbind, tapply(factor(Tab[[x]]), Factor, multiCat)), round(fisher.test(Reduce(cbind, tapply(factor(Tab[[x]]), Factor, table)))$p.value, 2))}))

colnames(StatMulti) <- c("All", "1", "2", "3", "4", "5", "6", "P")

write.table(StatMulti, file = "ImmuSupr_severity.chronic.xls", sep = "\t", quote = F, col.names = NA)

