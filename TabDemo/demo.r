
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


Tab <- read.delim(file = "Demography_merged.full_MM.txt", row.names = 1, as.is = T)

Tab["Monte006SA", "AKI"] <- 0

Tab["Monte012RN", "AKI"] <- 1

Tab$MMFBefore <- factor(Tab$MMFBefore)

Tab$Steroids <- factor(Tab$Steroids)

Tab$FKT <- factor(Tab$FKT)

Tab$SeverityScore <- factor(Tab$SeverityScore)

Tab$transplant_vintage <- Tab$transplant_vintage/365.25


Bivariate <- c("Induction", "CNI", "Race", "Maintenance.Steroids", "MMF", "Smoking", "Diabetes", "HTN", "ACEI.ARB", "Donor.Status", "Respiratory.symptoms", "Fever", "GI.symptoms", "Myalgia", "Sex", "AKI", "Graftloss_Within_1_year", "Death")

BivariateCat <- c("Batch", "Cohort")

MutiVariate <- c("Native.kidney.disease", "FKT", "MMFBefore", "Steroids", "SeverityScore")

Numeric <- c("BMI", "Age", "Duration", "Steroids_during_admission", "Lymph", "Baseline_SCr", "Peak_SCr_COVID", "transplant_vintage")

## Condition

Factor <- Tab$Condition

StatBi <- t(apply(Tab[, Bivariate, drop = F], 2, function(x) {c(biCat(x), tapply(x, Factor, biCat), round(fisher.test(Reduce(cbind, tapply(x, Factor, tableF)))$p.value, 2))}))

colnames(StatBi) <- c("All", "Acute", "Chronic", "Pval")

StatBiCat <- Reduce(rbind, lapply(BivariateCat, function(x) {cbind(multiCat(factor(Tab[[x]])), Reduce(cbind, tapply(factor(Tab[[x]]), Factor, multiCat)), round(fisher.test(Reduce(cbind, tapply(factor(Tab[[x]]), Factor, table)))$p.value, 2))}))

colnames(StatBiCat) <- c("All", "Acute", "Chronic", "Pval")

StatMulti <- Reduce(rbind, lapply(MutiVariate, function(x) {cbind(multiCat(Tab[[x]]), Reduce(cbind, tapply(Tab[[x]], Factor, multiCat)), round(chisq.test(Reduce(cbind, tapply(Tab[[x]], Factor, table)))$p.value, 2))}))

colnames(StatMulti) <- c("All", "Acute", "Chronic", "Pval")

StatNum <- t(apply(Tab[, Numeric, drop = F], 2, function(x) {Tmp <- data.frame(Val = x, Fac = Factor); c(numSum(x), tapply(x, Factor, numSum), round(t.test(Val ~ Fac, data = Tmp)$p.value, 2))}))

colnames(StatNum) <- c("All", "Acute", "Chronic", "Pval")

Stat <- rbind(StatBi, StatBiCat, StatMulti, StatNum)

write.table(Stat, file = "Stat_condition.xls", sep = "\t", quote = F, col.names = NA)

#StatMulti <- apply(Tab[, MutiVariate, drop = F], 2, function(x) {cbind(multiCat(x), Reduce(cbind, tapply(x, Factor, multiCat)), round(chisq.test(Reduce(cbind, tapply(x, Factor, table)))$p.value, 2))})


## Viromia

ExprViro <- read.delim(file = "/sc/arion/projects/GOCAR/Sun/5.DRWork/8.COVID/3.Ana_COVIDAll/Expr.xls.sym", row.names = 1, check.names = F)

Viremia.lst <- unlist(read.table(file = "Viremia.lst", as.is = T))

Factor <- rep("No", ncol(ExprViro))

Factor[colnames(ExprViro) %in% Viremia.lst] <- "Yes"

#Factor <- factor(ExprViro["S", rownames(Tab)] > 0)

levels(Factor) <- c("No", "Yes")

StatBi <- t(apply(Tab[, Bivariate, drop = F], 2, function(x) {c(biCat(x), tapply(x, Factor, biCat), round(fisher.test(Reduce(cbind, tapply(x, Factor, tableF)))$p.value, 2))}))

colnames(StatBi) <- c("All", "No", "Yes", "Pval")

#StatBiCat <- Reduce(rbind, lapply(BivariateCat, function(x) {cbind(multiCat(Tab[[x]]), Reduce(cbind, tapply(Tab[[x]], Factor, multiCat)), round(fisher.test(Reduce(cbind, tapply(Tab[[x]], Factor, table)))$p.value, 2))}))

StatBiCat <- Reduce(rbind, lapply(BivariateCat, function(x) {cbind(multiCat(factor(Tab[[x]])), Reduce(cbind, tapply(factor(Tab[[x]]), Factor, multiCat)), round(fisher.test(Reduce(cbind, tapply(factor(Tab[[x]]), Factor, table)))$p.value, 2))}))

colnames(StatBiCat) <- c("All", "No", "Yes", "Pval")

#StatMulti <- Reduce(rbind, lapply(MutiVariate, function(x) {cbind(multiCat(Tab[[x]]), Reduce(cbind, tapply(Tab[[x]], Factor, multiCat)), round(chisq.test(Reduce(cbind, tapply(Tab[[x]], Factor, table)))$p.value, 2))}))

StatMulti <- Reduce(rbind, lapply(MutiVariate, function(x) {cbind(multiCat(factor(Tab[[x]])), Reduce(cbind, tapply(factor(Tab[[x]]), Factor, multiCat)), round(chisq.test(Reduce(cbind, tapply(factor(Tab[[x]]), Factor, table)))$p.value, 2))}))

colnames(StatMulti) <- c("All", "No", "Yes", "Pval")

StatNum <- t(apply(Tab[, Numeric, drop = F], 2, function(x) {Tmp <- data.frame(Val = x, Fac = Factor); c(numSum(x), tapply(x, Factor, numSum), round(t.test(Val ~ Fac, data = Tmp)$p.value, 2))}))

colnames(StatNum) <- c("All", "No", "Yes", "Pval")

Stat <- rbind(StatBi, StatBiCat, StatMulti, StatNum)

write.table(Stat, file = "Stat_Viro.xls", sep = "\t", quote = F, col.names = NA)


TabAll <- Tab

## Acute

Tab <- TabAll[TabAll$Condition == "Acute", ]

Tab$SeverityScore <- factor(Tab$SeverityScore)

levels(Tab$SeverityScore) <- c("Low", "Low", "Median", "Median", "High", "High")


Bivariate <- c("Induction", "CNI", "Race", "Maintenance.Steroids", "MMF", "Smoking", "Diabetes", "HTN", "ACEI.ARB", "Donor.Status", "Respiratory.symptoms", "Fever", "GI.symptoms", "Myalgia", "Sex", "AKI", "Graftloss_Within_1_year", "Death")

BivariateCat <- c("Batch", "Cohort")

MutiVariate <- c("Native.kidney.disease", "FKT", "MMFBefore", "Steroids")

Numeric <- c("BMI", "Age", "Duration", "Steroids_during_admission", "Lymph", "Baseline_SCr", "Peak_SCr_COVID", "transplant_vintage")


Factor <- Tab$SeverityScore

StatBi <- t(apply(Tab[, Bivariate, drop = F], 2, function(x) {c(biCat(x), tapply(x, Factor, biCat), round(fisher.test(Reduce(cbind, tapply(x, Factor, tableF)))$p.value, 2))}))

colnames(StatBi) <- c("All", "Low", "Median", "High", "Pval")

StatBiCat <- Reduce(rbind, lapply(BivariateCat, function(x) {cbind(multiCat(Tab[[x]]), Reduce(cbind, tapply(Tab[[x]], Factor, multiCat)), round(fisher.test(Reduce(cbind, tapply(Tab[[x]], Factor, table)))$p.value, 2))}))

colnames(StatBiCat) <- c("All", "Low", "Median", "High", "Pval")

StatMulti <- Reduce(rbind, lapply(MutiVariate, function(x) {cbind(multiCat(Tab[[x]]), Reduce(cbind, tapply(Tab[[x]], Factor, multiCat)), round(chisq.test(Reduce(cbind, tapply(Tab[[x]], Factor, table)))$p.value, 2))}))

colnames(StatMulti) <- c("All", "Low", "Median", "High", "Pval")

StatNum <- t(apply(Tab[, Numeric, drop = F], 2, function(x) {Tmp <- data.frame(Val = x, Fac = Factor); c(numSum(x), tapply(x, Factor, numSum), round(summary(aov(Val ~ Fac, data = Tmp))[[1]][["Pr(>F)"]][1], 2))}))

colnames(StatNum) <- c("All", "Low", "Median", "High", "Pval")

Stat <- rbind(StatBi, StatBiCat, StatMulti, StatNum)

write.table(Stat, file = "Stat_severity_acute.xls", sep = "\t", quote = F, col.names = NA)

## Chronic

#Tab <- TabAll[TabAll$Condition == "Chronic", ]

#Tab$SeverityScore <- factor(Tab$SeverityScore)

#levels(Tab$SeverityScore) <- c("Low", "Low", "Median", "Median", "High", "High")

