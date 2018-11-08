source("http://bioconductor.org/biocLite.R")
biocLite("rgl")
library(limma)
library("rgl")

setwd("Path/to/directory/")
pheno <- read.csv(file = "pheno.csv")
hdat <- read.csv(file = "data.csv", row.names = 1)


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

fit1 <- contrasts.fit(fit, c(-1,0,0,1,0,0,0,0))



fit1 <- eBayes(fit1)
toptable1 <- toptable(fit1)
nrow(topTable(fit1, coef=1, number=10000, lfc=2))
probeset.list1 <- topTable(fit1, coef=1, number=10000, lfc=2)
head(probeset.list1)

write.table(probeset.list1, file = "Output.txt", sep="\t", quote = F, row.names=T)
volcanoplot(fit1,  main = "Day-0 VS Day-3")
fit1
