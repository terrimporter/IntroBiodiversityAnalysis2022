# Teresita M. Porter, Aug. 4, 2022

library(goeveg) # scree
library(plyr) # ddply

# ESV table, under-sequenced samples removed, rare ESVs removed
t <- read.csv("outfiles/ESVtable.csv", header=TRUE, row.names = 1, stringsAsFactors = FALSE)

# normalize ----
# convert to proportions (reads / total reads per sample) if using metric that needs quantitative data
# t2 <- t / rowSums(t)
# sanity check, row sums should equal 1
# rowSums(t2)

# convert to proportions if using a metric that needs presence-absence data
# t[t>0] <- 1

# Scree ----
# Scree plots to determine number of dimensions to use for NMDS
pdf("outfiles/scree.pdf")
# check dims
dimcheckMDS(t)
dev.off()
# use k=3

# NMDS ----

# binary Bray Curtis = Sorensen dissimilarity (presence-absence)
m <- vegdist(t, method="bray", binary=TRUE)

# Do 3 dimensional NMDS
nmds3 <- metaMDS(m, k=3, trymax=100)
# stress = 0.07397908

# Goodness of fit ----
# Stressplot Shephards curve to assess goodness of fit between observed and ordination distances
pdf("outfiles/stressplot.pdf")
stressplot(nmds3)
gof <-goodness(nmds3)
gof
plot(nmds3, display = "sites", type="n", main="")
points(nmds3, display="sites",cex=2*gof/mean(gof))
dev.off()
# linear R2 = 0.976

# Create grouping matrix for samples by grabbing row names from above matrix
names <- data.frame(row.names(t), stringsAsFactors = FALSE)
# Rename the column
names(names) <- "sample"
# Copy column to row names
row.names(names) <- names$sample
# Split first column into their own fields
names2 <- data.frame(names, do.call(rbind, strsplit(names$sample,'_')), stringsAsFactors = FALSE)
names(names2)[2:5]<-c("time", "site", "rep", "preservative")
# Remove first column
names2 <- names2[,-1]
# Grab sites/species scores from NMDS output
df <- data.frame(scores(nmds3, display = "sites"))
# Put it all in one df for ggplot
gg <- merge(df, names2, by="row.names")

# create factors
gg$time <- factor(gg$time, levels=c("T0", "T1", "T2"))
gg$preservative <- factor(gg$preservative, 
                            levels=c("B", "A", "N", "R", "Z", "Q"),
                            labels=c("Bas", "Ant", "Nor", "RNA", "Zym", "Qia"))
gg$site <- factor(gg$site, levels=c("S1", "S2", "S3"))

# # enclose each site with a polygon (hulls)
# chulls.site <- ddply(gg, .(site), function(gg) gg[chull(gg$NMDS1, gg$NMDS2), ])

# NMDS plot, add sig env vars
p <- ggplot(data=gg, aes(x=NMDS1, y=NMDS2)) + 
  geom_point(data=gg, size = 2.5, aes(color = site, shape=preservative)) +
  facet_wrap(~time) +
  theme_classic() +
  theme(legend.position="bottom") 
p

ggsave("outfiles/NMDS.pdf", p, width = 8, height = 4)


# Beta dispersion ----

# If significant beta dispersion, PERMANOVA will reflect variation due to within and among-sample variation
bd.time <- betadisper(m, as.factor(gg$time))
bd.preservative <- betadisper(m, as.factor(gg$preservative))
bd.site <- betadisper(m, as.factor(gg$site))

set.seed(1234)
anova(bd.time) # n/s
anova(bd.preservative) # 0.0003303 ***
anova(bd.site) # 2.688e-06 ***

pdf("outfiles/betadisp.pdf")
par(mfrow = c(2,2))
boxplot(bd.time, main="Time")
boxplot(bd.preservative, main="Preservative")
boxplot(bd.site, main="Site")
dev.off()

# So there is significant variation within groups of preservatives and sites

# PERMANOVA ----

# See which groupings explain a significant amount of variance in beta diversity across groups
adonis2(m ~ time, data=gg, permutations=999)
#           Df SumOfSqs      R2      F Pr(>F)
# time      2   0.3627 0.02269 0.9985  0.381
# Residual 86  15.6213 0.97731              
# Total    88  15.9840 1.00000 

# Time does not explain a significant amount of variation in beta diversity

adonis2(m ~ preservative, data=gg, permutations=999)
#           Df SumOfSqs      R2      F Pr(>F)
# preservative  5   1.3693 0.08567 1.5554   0.04 *
# Residual     83  14.6147 0.91433                
# Total        88  15.9840 1.00000 

# preservative explains 8.5% of the variation in beta diversity (p=0.04) 
# due to both variation within and among groups

adonis2(m ~ site, data=gg, permutations=999)
#           Df SumOfSqs      R2      F Pr(>F)
# site      2   8.0168 0.50155 43.267  0.001 ***
# Residual 86   7.9672 0.49845                  
# Total    88  15.9840 1.00000

# site explains 50%% of the variation in beta diversity (p=0.001) 
# due to both variation within and among groups

