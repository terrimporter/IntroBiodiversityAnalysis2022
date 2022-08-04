# Teresita M. Porter, Aug. 4/22

library(vegan) # specnum
library(stringr) # str_split

# ESV.table3 under-sequenced samples removed, rare ESVs removed
t <- read.csv("outfiles/ESVtable.csv", header=TRUE, row.names=1, stringsAsFactors = FALSE)

# ESV Richness ----
r <- data.frame(sample=rownames(t), richness=specnumber(t))

# split up sample into fields
# Parse out metadata from sample
r.df <- data.frame(r, do.call(rbind, str_split(r$sample,"_")), stringsAsFactors = FALSE)
names(r.df)[3:6]<-c("time","site","replicate","preservative")

# create factors
r.df$time <- factor(r.df$time, levels=c("T0", "T1", "T2"))
r.df$preservative <- factor(r.df$preservative, 
                               levels=c("B", "A", "N", "R", "Z", "Q"),
                               labels=c("Bas", "Ant", "Nor", "RNA", "Zym", "Qia"))

# plot
p <- ggplot(r.df, aes(x=preservative, y=richness)) + 
  geom_boxplot() +
  labs(x="Preservative", y="ESV Richness") +
  facet_wrap(~time) +
  theme_classic()

ggsave("outfiles/richness.pdf", p, height = 4, width = 8)


