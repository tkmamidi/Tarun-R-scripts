source("http://bioconductor.org/biocLite.R")
biocLite("rgl")
biocLite("ArrayExpress")
vignette(package="Biobase", "ExpressionSetIntroduction")
library("Biobase")
library(limma)
library(pd.hugene.2.0.st)
library(hugene20sttranscriptcluster.db)
library(annotate)
library(genefilter)
library(lattice)
install.packages("rgl",dependencies = TRUE)
install.packages("digest", type="source")
library("rgl")

library(gplots)
library(NMF)

setwd("Z:/Tarun_projects/Tanja_project/limma")
pheno <- read.csv(file = "All_data_pheno.csv")
hdat <- read.csv(file = "all_data.csv", row.names = 1)


rownames(pheno) <- pheno$sample
colnames(hdat) <- rownames(pheno)
rownames(pheno) == colnames(hdat)

##### PCA
sample.groups <- c(rep("red", 12), rep("blue",5), rep("green", 6), rep("purple", 7), rep("orange", 7), rep("black", 7), rep("yellow", 7), rep("brown", 12))
pc <- princomp(hdat)
plot(pc$loadings, col=sample.groups, xlab= "Component 1", ylab="Component 2", main = "PCA")
plot3d(pc$loadings, col=sample.groups, xlab= "Component 1", ylab="Component 2", zlab="Component 3", type = "s", size = 1.5)

pc$loadings

samples <- factor(pheno$time)
design <- model.matrix(~0 + samples)
colnames(design) <- levels(samples)
design

fit <- lmFit(hdat, design)
names(fit)

fit1 <- contrasts.fit(fit, c(-1,0,0,-1,0,0,1,1))



fit1 <- eBayes(fit1)
toptable1 <- toptable(fit1)
nrow(topTable(fit1, coef=1, number=10000, lfc=2))
probeset.list1 <- topTable(fit1, coef=1, number=10000, lfc=2)
head(probeset.list1)

write.table(probeset.list1, file = "Day-0,3-vs-12,14.txt", sep="\t", quote = F, row.names=T)
volcanoplot(fit1,  main = "Day-0,3 VS Day-12,14")
fit1
