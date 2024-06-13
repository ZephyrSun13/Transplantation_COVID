
library(ggplot2)
library(RColorBrewer)

Stat <- read.delim(file = "map.stat", stringsAsFactors = FALSE, row.names = 1)

Stat <- data.frame(t(Stat))

Stat <- Stat[, 9:10]

names(Stat) <- c("ReadCount", "MappingRate")

Stat$ReadCount <- as.numeric(as.character(Stat$ReadCount))

Stat$MappingRate <- as.numeric(gsub("%", "", as.character(Stat$MappingRate)))

Tab <- Stat

        Tab <- Tab[order(Tab[,"ReadCount"], Tab[, "MappingRate"], decreasing = FALSE), ]

        Tab$Sam<- factor(rownames(Tab), levels = rownames(Tab))

        pp <- ggplot(data=Tab, aes_string(x="Sam", y="ReadCount")) +

                geom_bar(stat="identity", width=0.5, fill = brewer.pal(12, "Set3")[4]) +

		theme_minimal() +

		xlab("") +
	
                theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x=element_blank())

	ggsave("ReadCount.pdf", pp, width = 12)

       	Tab <- Tab[order(Tab[,"MappingRate"], Tab[, "ReadCount"], decreasing = FALSE), ]

        Tab$Sam<- factor(rownames(Tab), levels = rownames(Tab))

        pp <- ggplot(data=Tab, aes_string(x="Sam", y="MappingRate")) +

                geom_bar(stat="identity", width=0.5, fill = brewer.pal(12, "Set3")[4]) +

                theme_minimal() +

                xlab("") +

                ylim(0, NA) +

                theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x=element_blank()) 

        ggsave("MappingRate.pdf", pp, width = 12)

