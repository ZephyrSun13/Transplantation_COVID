
tableF <- function(tt){

        tmp <- table(tt)[c("0", "1")]

        tmp[is.na(tmp)] <- 0

        tmp

}

biCat <- function(tt){

        Tmp <- table(tt)

        paste0(Tmp["1"], "/", sum(Tmp), " (", round(Tmp["1"]/(sum(Tmp))*100, 1), ")")

}

multiCat <- function(tt){

        Tmp <- table(tt)

        Resu <- paste0(Tmp, "/", sum(Tmp), " (", round(Tmp/sum(Tmp)*100, 1), ")")

        names(Resu) <- names(Tmp)

        Resu

}

numSum <- function(tt){

        paste0(round(mean(tt, na.rm = T), 2),  "+-", round(sd(tt, na.rm = T), 2))

}

Tab <- read.delim(file = "Demography_merged.full_90.txt", row.names = 1, as.is = T)

Tab$transplant_vintage <- Tab$transplant_vintage/365

rownames(Tab) <- gsub(" ", "", rownames(Tab))

Tab.seq <- read.delim(file = "Demography_merged.full_MM.txt", row.names = 1, as.is = T)

rownames(Tab.seq) <- gsub(" ", "", Tab.seq$Code)

Tab$Seq <- "No"

Tab$Seq[rownames(Tab) %in% rownames(Tab.seq)] <- "Seq"

Tab.1 <- Tab

Tab.1$Seq <- "All"

Tab.2 <- Tab[Tab$Seq == "Seq", ]

Tab <- rbind(Tab.1, Tab.2)

Bivariate <- c("Race", "Smoking", "Diabetes", "HTN", "DonorStatus", "Sex", "Death", "Graftloss_Within_1_year", "Induction")

MutiVariate <- c("SeverityScore")

Numeric <- c("Age", "DaysAfPCR", "transplant_vintage")

Factor <- factor(Tab$Seq, levels = c("All", "Seq"))

StatBi <- t(apply(Tab[, Bivariate, drop = F], 2, function(x) {c(biCat(x), tapply(x, Factor, biCat), round(fisher.test(Reduce(cbind, tapply(x, Factor, tableF)))$p.value, 2))}))

colnames(StatBi) <- c("CC", "All", "Seq", "Pval")

StatMulti <- Reduce(rbind, lapply(MutiVariate, function(x) {cbind(multiCat(Tab[[x]]), Reduce(cbind, tapply(Tab[[x]], Factor, multiCat)), round(chisq.test(Reduce(cbind, tapply(Tab[[x]], Factor, table)))$p.value, 2))}))

colnames(StatMulti) <- c("CC", "All", "Seq", "Pval") 

StatNum <- t(apply(Tab[, Numeric, drop = F], 2, function(x) {Tmp <- data.frame(Val = x, Fac = Factor); c(numSum(x), tapply(x, Factor, numSum), round(t.test(Val ~ Fac, data = Tmp)$p.value, 2))}))

colnames(StatNum) <- c("CC", "All", "Seq", "Pval")

Stat <- rbind(StatBi, StatMulti, StatNum)

write.table(Stat, file = "Stat_All89.xls", sep = "\t", quote = F, col.names = NA)

