
library(ggplot2)
library(RColorBrewer)
library(gridExtra)

Stat <- read.delim(file = "map.stat", stringsAsFactors = FALSE)

Stat[1, 1] <- "Batch"

rownames(Stat) <- Stat[[1]]

Stat <- Stat[, -1]

Stat <- data.frame(t(Stat))


Check <- c("X..........................Number.of.input.reads..", "X........................Uniquely.mapped.reads....", "X...................Uniquely.mapped.reads.number..", "X...............of.reads.mapped.to.multiple.loci..", "X...................of.reads.unmapped..too.short..")

Map <- c("Total Number of Reads", "Percentage of Unique Mapped", "Number of Unique Mapped", "Percentage of Multiple Mapped", "Percentage of Unmapped")

names(Map) <- Check

Plots <- list()

for(nn in Check){

        Tmp <- Stat

        Tmp[[nn]] <- as.numeric(gsub("%" ,"" , Stat[[nn]]))

        Tmp <- Tmp[order(Tmp[[nn]], decreasing = T), ]

        Tmp$Sam <- factor(rownames(Tmp), levels = rownames(Tmp))

        pp <- ggplot(Tmp, aes_string(x = "Sam", y = nn, fill = "Batch")) +

                geom_bar(stat = "identity") +

                scale_fill_manual(values = brewer.pal(8, "Set2")) +

                xlab("") +

                ylab(Map[nn]) +

                theme_minimal() +

                theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6), legend.title = element_blank())

        Plots[[nn]] <- pp

}

ggsave("MapedStat.pdf", arrangeGrob(grobs = Plots, ncol = 1), height = 3*length(Plots), width = 10)

Stat <- Stat[, c("Batch", Check)]

names(Stat) <- c("Batch", gsub(" ", "_", Map))

write.table(Stat, file = "map.stat.clean.xls", sep = "\t", quote = F, col.names = NA)

