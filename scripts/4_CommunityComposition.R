# Teresita M. Porter, Aug. 4, 2022

library(data.table) # setDT
library(tidyr) # gather

# ESV table ----

# ESV table, under-sequenced samples removed, rare ESVs removed
tab <- read.csv("outfiles/ESVtable.csv", header=TRUE, row.names=1, stringsAsFactors = FALSE)
# transpose to get ESVs in columns and set rownames
tab.t <- data.frame(t(tab))
# move rownames to first col
setDT(tab.t, keep.rownames = TRUE)[]
names(tab.t)[1] <- "GlobalESV"

# convert from wide to long format
tab.long <- gather(tab.t, sample, reads, -GlobalESV)

# Taxonomy ----
# read in filtered taxonomy
tax <- read.csv("outfiles/taxonomy.csv", header=TRUE, row.names = 1, stringsAsFactors = FALSE)

# map taxonomy on to ESV table, focus on family level diversity
tab2 <- merge(tab.long, tax, by="GlobalESV", all=TRUE)

# if Phylum "" then "Unknown"
tab2$Phylum <- ifelse(tab2$Phylum=="", "Unknown", tab2$Phylum)

# Phylum ----

# summarize community composition at the phylum rank
gg <- data.frame(tab2 %>% group_by(sample,Phylum) %>% dplyr::summarize(sum(reads)))

# split up sample name to get separate fields
gg <- data.frame(gg, do.call(rbind, strsplit(gg$sample,'_')), stringsAsFactors = FALSE)
names(gg)[4:7] <- c("time", "site", "rep", "preservative")

names(gg)[3] <- "Reads"

# create factors
gg$time <- factor(gg$time, levels=c("T0", "T1", "T2"))
gg$preservative <- factor(gg$preservative, 
                            levels=c("B", "A", "N", "R", "Z", "Q"),
                            labels=c("Bas", "Ant", "Nor", "RNA", "Zym", "Qia"))
gg$site <- factor(gg$site, levels=c("S1", "S2", "S3"))

# visualize as stacked bar plot
p <- ggplot(gg, aes(fill=Phylum, y=Reads, x=preservative)) + 
  geom_bar(position="stack", stat="identity") +
  facet_wrap(~time) +
  theme_classic()
p

ggsave("outfiles/composition.pdf", p, width = 8, height = 4)
# phylum has a manageable number of categories for a stacked bar plot
# for more categories (family, genus, species) a heat map would be better
