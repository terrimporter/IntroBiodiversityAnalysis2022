# Teresita M. Porter, Aug. 3, 2022

library(reshape2) # dcast, pivot table
library(stringr) # str_split
library(vegan) # rarecurve
library(ggplot2)
library(purrr) # map_dfr
library(scales) # comma

# Read in ESV
a <- read.csv("infiles/Brennan_RNA1_results.csv", header=TRUE, stringsAsFactors = FALSE)

# Data cleaning ----
# remove fields that aren't important
a$SampleName <- gsub("Brennan_RNAPres1_", "", a$SampleName)

# dim(a)
# [1] 21390    39

# PCR negative controls (NC)
ctrls <- unique(a$SampleName[grepl("NC", a$SampleName)])

# ctrls
# time_ctrltype_ctrl_replicate_illuminaSample_primer
# [1] "T0_PCR1_NC_R3_S88_BR5" 
# time_ctrltype_ctrl_replicate_illuminaSample_primer ATM/AJM? RM?
# "T2_ATM_RM_NC_S96_MLJG"
# [3] "T0_ExNC_R3_B_S86_MLJG" "T2_AJM_RM_NC_S94_MLJG"
# [5] "T1_PCR1_R1_NC_S92_BR5" "T2_BS_RS_NC_S93_MLJG" 
# [7] "T0_ExNC_R3_B_S86_BR5"  "T1_ATM_RM_NC_S89_BR5" 
# [9] "T1_BS_RS_NC_S90_BR5"   "T2_R1_R1_NC_S95_MLJG" 
# [11] "T2_R1_R1_NC_S95_BR5" 

# exclude controls for now
a2 <- a[!(a$SampleName %in% ctrls),]

# file already has separated columns for amplicon(x), primer, site, replicate, pres? 
# remove these and start fresh
a2 <- a2[,-c(33:39)]

# Split up SampleName with pkg 'stringr' to get amplicon & replicate fields
a3 <- data.frame(a2, do.call(rbind, str_split(a2$SampleName,"_")), stringsAsFactors = FALSE)
names(a3)[33:39] <- c("time", "site", "replicate", "preservative", "marker", "IlluminaSample", "amplicon")

# ensure data in columns contain the expected type of data
# time
unique(a3$time)
# [1] "T0" "T2" "T1"
# site
unique(a3$site)
# [1] "S2" "S3" "S1"
# replicate
unique(a3$replicate)
# [1] "R3" "R2" "R1"
unique(a3$preservative)
# [1] "B" "Z" "N" "A" "Q" "R"
unique(a3$marker) 
# [1] "COI"
unique(a3$amplicon)
# [1] "MLJG" "BR5" 

# ESV table ----

# create new sample name with just the essential fields

# extract ESV.table, pivot table
ESV.table <- reshape2::dcast(a3, time+site+replicate+preservative ~ GlobalESV, value.var = "ESVsize", fun.aggregate = sum)

# set rownames
rownames(ESV.table) <- paste(ESV.table$time, ESV.table$site, ESV.table$replicate, ESV.table$preservative, sep="_")
ESV.table <- ESV.table[,-c(1:4)]

#remove columns with only zeros
ESV.table <- ESV.table[,colSums(ESV.table) !=0]

#remove rows with only zeros & edit rownames
ESV.table <- ESV.table[rowSums(ESV.table) !=0,]

# Remove under-sequenced samples ----

# visualize read depth distribution
hist(as.numeric(rowSums(ESV.table)))

# remove samples with < 10000 reads
ESV.table2 <- ESV.table[!rowSums(ESV.table) < 10000,]
# removed 10 samples
setdiff(rownames(ESV.table), rownames(ESV.table2))
#  [1] "T0_S1_R1_B" "T0_S2_R2_B" "T1_S1_R1_Q" "T1_S1_R2_A" "T1_S1_R2_N"
# [6] "T1_S2_R1_Z" "T1_S2_R2_A" "T1_S2_R2_B" "T2_S1_R1_Q" "T2_S1_R2_Q"
# dropped samples come from across different preservative types so hopefully ok

# Remove rare ESVs ----

# remove ESVs that represent less than 0.0001 or 0.01% of all reads
cutoff <- sum(colSums(ESV.table2)) * 0.0001
# using this cutoff will remove ESVs with < ~ 216 reads
# [1] 216.1664

ESV.table3 <- ESV.table2
ESV.table3[colSums(ESV.table2) < cutoff] <- NULL
length(names(ESV.table2)) - length(names(ESV.table3))
# 6438 rare ESVs removed after this step
rare <- setdiff(names(ESV.table2), names(ESV.table3))

# Check sequencing depth ----

# visualize curves
rc <- rarecurve(ESV.table3, step=100, label=FALSE)

# extract fields for ggplot for a nicer figure that's easier to customize
# Reformat vegan list as df (cols OTU, raw.read)
rare.df <- lapply(rc, function(x){
  z <- as.data.frame(x)
  z <- data.frame(OTU = z[,1], raw.read = rownames(z))
  z$raw.read <- as.numeric(gsub("N", "",  z$raw.read))
  return(z)
})

# Add sample names to vegan output (df) (rownames)
sample_names <- rownames(ESV.table2)
names(rare.df) <- sample_names

# Map rownames to vegan output (df)
rare.df <- map_dfr(rare.df, function(x){
  z <- data.frame(x)
  return(z)
}, .id = "sample")

# Parse out metadata from sample
rare.df <- data.frame(rare.df, do.call(rbind, str_split(rare.df$sample,"_")), stringsAsFactors = FALSE)
names(rare.df)[4:7]<-c("time","site","replicate","preservative")

# create factors
rare.df$time <- factor(rare.df$time, levels=c("T0", "T1", "T2"))
rare.df$preservative <- factor(rare.df$preservative, 
                               levels=c("B", "A", "N", "R", "Z", "Q"),
                               labels=c("Baseline", "Antifreeze", "Norgen", "RNALater", "Zymo", "Qiagen"))

# color by amplicon
p1 <- ggplot(data = rare.df) +
  ggtitle("") +
  labs(x="Reads", y="ESVs") +
  geom_point(aes(x = raw.read, y = OTU, color = preservative), size=0.05) +
  scale_x_continuous(label = comma) +
  theme_bw() + 
  theme(text = element_text(size=8),
        plot.title = element_text(size=8, hjust=0.01),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(),
        axis.title.x = element_text(vjust = -0.75),
        legend.title = element_blank(),
        legend.text = element_text(size=8),
        legend.position = "bottom",
        legend.key.size = unit(0.4, "cm"),
        strip.background = element_blank(),
        strip.text.x = element_text(angle = 0, hjust = 0.5)) +
  guides(color = guide_legend(override.aes = list(size=2)))
p1

ggsave("outfiles/rarefaction.pdf", p1, width = 4, height = 4)

# looks good, from now on use ESV.table3 for all further analyses
write.csv(ESV.table3, "outfiles/ESVtable.csv", quote=FALSE)

# Filter taxonomic assignments ----

# See https://github.com/terrimporter/CO1Classifier for cutoffs
# Ensure 95% correct at species rank, 99%+ correct at genus+ ranks for a ~ 200bp frag of COI
a4 <- a3
a4$Species <- ifelse(a4$sBP >= 0.7, a4$Species, "") 
a4$Genus <- ifelse(a4$gBP >= 0.3, a4$Genus, "") 
a4$Family <- ifelse(a4$fBP >= 0.2, a4$Family, "")

# only keep essential fields needed to map taxonomy
# drop sequences to reduce filesize
a4 <- a4[,c(1,9,12,15,18,21,24,27,30)]

# drop the ESVs that were filtered out of the ESV table to reduce file size
a4 <- unique(a4[!a4$GlobalESV %in% rare,])

# from now on use a4 for all further analyses that need mapped taxonomy
write.csv(a4, "outfiles/taxonomy.csv", quote=FALSE)

