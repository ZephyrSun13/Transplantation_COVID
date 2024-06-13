
Tab <- read.delim(file = "Demography_merged.txt", row.names = 2)

Fac <- read.table(file = "../../3.AnaAll/Factor.txt", row.names = 1)

Tab$Condition <- Fac[rownames(Tab), "AcutevsChronic"]

Tab$Batch <- Fac[rownames(Tab), "Batch"]

Tab$SeverityScore <- Fac[rownames(Tab), "SeverityScore"]

Tab$Duration <- Fac[rownames(Tab), "Duration"]

Tab$Race <- Tab$Race - 1

Tab$Native.kidney.disease <- as.character(Tab$Native.kidney.disease)

Tab$Native.kidney.disease[Tab$Native.kidney.disease == " E"] <- "E"

Tab$Lymph <- Fac[rownames(Tab), "Lymph"]

Tab$Date.of.PCR.confirmation <- Fac[rownames(Tab), "Date.of.PCR.confirmation"]

Tab$EnrollmentDate <- Fac[rownames(Tab), "EnrollmentDate"]

write.table(Tab, file = "Demography_merged.full.xls", sep = "\t", quote = F, col.names = NA)

