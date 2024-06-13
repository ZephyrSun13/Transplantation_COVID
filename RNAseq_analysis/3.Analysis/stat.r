
#!/usr/bin/R

Args <- commandArgs(trailingOnly = TRUE)

Expr <- read.table(file = Args[1], sep = "\t", header = TRUE, row.names = 1, check.names = F)

#names(Expr) <- gsub("_star", "", names(Expr))

## PCA analysis

library(ggplot2)

Factor <- read.table(file = "Factor.xls", sep = "\t", header = TRUE, row.names = 1, as.is = T)

rownames(Factor) <- gsub(" sample 1", "", rownames(Factor))

Expr <- Expr[, names(Expr) %in% Factor$SeqID]

names(Expr) <- rownames(Factor)[match(names(Expr), Factor$SeqID)]

write.table(Expr, file = "Expr.match.xls", sep = "\t", quote = F, col.names = NA)

Factor <- Factor[names(Expr), ]

int_PCA <- prcomp(t(Expr))

#Factor <- data.frame(Name = paste0(Factor$Row.names, "_count.txt_norm"), Factor)

#Factor[["Row.names"]] <- paste0(Factor$Sample, "_count.txt_norm")

int_PCA_x <- merge(Factor, int_PCA$x[, 1:8], by.x = 0, by.y = 0)

checks <- colnames(Factor)[-1]

for(check in checks){

        #condition_col <- as.factor(Factor[[check]])

        #levels(condition_col) <- brewer.pal(nlevels(condition_col), "Set3")[1:nlevels(condition_col)]

        #plotMDS(int, labels = gsub("_norm", "", gsub("Sample_", "", colnames(int))), col = as.character(condition_col), dim = c(1, 2))

        p<-ggplot(int_PCA_x,aes_string(x="PC1",y="PC2",color=check, label="Row.names"))

        p<-p+geom_point() + geom_text(size=3) + theme_minimal()

        ggsave(paste0(check, "_MDS_dim_1_2.pdf"), p)

        p<-ggplot(int_PCA_x,aes_string(x="PC3",y="PC4",color=check, label="Row.names"))

        p<-p+geom_point() + geom_text(size=3) + theme_minimal()

        ggsave(paste0(check, "_MDS_dim_3_4.pdf"), p)

        p<-ggplot(int_PCA_x,aes_string(x="PC5",y="PC6",color=check, label="Row.names"))

        p<-p+geom_point() + geom_text(size=3) + theme_minimal()

        ggsave(paste0(check, "_MDS_dim_5_6.pdf"), p)

        p<-ggplot(int_PCA_x,aes_string(x="PC7",y="PC8",color=check, label="Row.names"))

        p<-p+geom_point() + geom_text(size=3) + theme_minimal()

        ggsave(paste0(check, "_MDS_dim_7_8.pdf"), p)

}


Expr <- read.table(file = "Expr.xls", row.names = 1, header = TRUE, sep = "\t")

Ann <- read.table(file = "/sc/arion/projects/zhangw09a/PANDA/db_ZS/Ensembl/gene_symbol.txt", header = TRUE, sep = "\t", fill = TRUE, stringsAsFactors = FALSE)

Expr <- merge(Ann, Expr, by.x =1 , by.y = 0)

row.names(Expr) <- make.names(Expr[[2]], unique = TRUE)

Expr <- Expr[, c(-1, -2)]

write.table(Expr, file = "Expr.xls.sym", sep = "\t", col.names = NA, quote = FALSE)

