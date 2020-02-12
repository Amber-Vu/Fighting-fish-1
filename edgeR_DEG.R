

# TMM normalize
library(edgeR)

dat <- read.delim("input file", row.names=1)

head(dat)


# A: before 
# B: 20m
# C: 60m

grp <- c(rep("A", 5), rep("B", 10), rep("C", 10))

D <- DGEList(dat, group=grp)
D$samples

D <- calcNormFactors(D, method="TMM")    # TMM normalization

# dump normalized count data
D.cpm.tmm <- cpm(D, normalized.lib.size=T)
write.table(D.cpm.tmm, file="Nor_TMM.txt", sep="\t", quote=F)

D <- estimateCommonDisp(D)
D$common.dispersion

# Before vs 60m (A vs C)
de.tagwise.AC <- exactTest(D, pair=c("A", "C"))
tmp.AC <- topTags(de.tagwise.AC, n=nrow(de.tagwise.AC$table))
write.table(tmp.AC$table, "de.tagwise_Before_vs_D60.txt", quote=F)

# Before vs 20m (A vs B)
de.tagwise.AB <- exactTest(D, pair=c("A", "B"))
tmp.AB <- topTags(de.tagwise.AB, n=nrow(de.tagwise.AB$table))
write.table(tmp.AB$table, "de.tagwise_Before_vs_D20.txt", quote=F)

# 20m vs 60m (B vs C)
de.tagwise.BC <- exactTest(D, pair=c("B", "C"))
tmp.BC <- topTags(de.tagwise.BC, n=nrow(de.tagwise.BC$table))
write.table(tmp.BC$table, "de.tagwise_20m-60m.txt", quote=F)










