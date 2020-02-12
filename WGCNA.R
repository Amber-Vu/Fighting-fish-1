# WGCNA analysis

install.packages("BiocManager")
BiocManager::install("WGCNA")
library(WGCNA)
options(stringsAsFactors = FALSE)


femData <- read.delim("input file",row.names=1)
expData <- data.frame(femData)
dim(femData)
names(femData)
head(femData)
datExpr0 = as.data.frame(t(expData[c(1:25)]))
names(datExpr0) = expData$Geneid
rownames(datExpr0) = names(expData)[c(1:25)]


gsg = goodSamplesGenes(datExpr0, verbose = 3);
gsg$allOK
  if (!gsg$allOK)
  {
    if (sum(!gsg$goodGenes)>0) 
      printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
    if (sum(!gsg$goodSamples)>0) 
      printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
    datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
  }


sampleTree = hclust(dist(datExpr0), method = "average");
sizeGrWindow(12,9)
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 2)


abline(h = 19000, col = "red");# cut the height from 15 up
clust = cutreeStatic(sampleTree, cutHeight = 19000, minSize = 10)
table(clust)
keepSamples = (clust==1)
datExpr = datExpr0[keepSamples, ]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)



traitData <- read.delim("Trait file",row.names=1)
dim(traitData)
names(traitData)
allTraits = traitData[, c(1,3:5) ]

dim(allTraits)
names(allTraits)


femaleSamples = rownames(datExpr)
traitRows = match(femaleSamples, allTraits$Fish)
datTraits = allTraits[traitRows, -1]
rownames(datTraits) = allTraits[traitRows, 1]

collectGarbage()



sampleTree2 = hclust(dist(datExpr), method = "average")
traitColors = numbers2colors(datTraits, signed = FALSE);
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits), 
                    main = "Sample dendrogram and trait heatmap")
save(datExpr, datTraits, file = "FemaleLiver-01-dataInput.RData") 

  powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
abline(h=0.90,col="red")

plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")


net = blockwiseModules(datExpr, power = 6,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, 
                      # mergeCutHeight = 0.355, # ?????????????????????chao-Li ?????????module????????????
                       mergeCutHeight = 0.25, # amber ???????????????threshold
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "femaleMouseTOM", 
                       verbose = 3)


table(net$colors)

sizeGrWindow(12, 9)

mergedColors = labels2colors(net$colors)

plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];
save(MEs, moduleLabels, moduleColors, geneTree, 
    file = "FemaleLiver-02-networkConstruction-auto.RData")  ### ???????????????



  nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);


sizeGrWindow(10,6)
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));

labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.4,
               cex.axis = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))




