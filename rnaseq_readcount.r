library(limma)
library(edgeR)
library(Glimma)
library(gplots)
source("https://bioconductor.org/biocLite.R")
biocLite("org.Hs.eg.db")
library(org.Hs.eg.db)
library(RColorBrewer)
library(annotate)


setwd("Z:/Tarun_projects/prostate_cancer_new/DE_in_R/DEseq")

# Read the data into R
countdata <- read.delim("data.txt",row.names = 1,check.names=F)

head(countdata)

# Read the sample information into R
sampleinfo <- read.delim("pheno1.txt")

#geneinfo <- read.delim("mart_export.txt")
#sampleinfo
#colnames(countdata)=sampleinfo$PRAD
#colnames(countdata)

#check samples names are same in data and pheno files
colnames(countdata)==sampleinfo$PRAD
table(colnames(countdata)==sampleinfo$PRAD)

# Obtain CPMs
myCPM <- cpm(countdata)
head(myCPM)

selr <- rowSums(cpm(countdata)>0.5)>=165
counts.keep <- countdata[selr,]


#extra filters (optional)

#selr<-rowSums(countdata) >= 200
selc <- colSums(countdata)>= 10000000 
counts.keep <- countdata[selr,selc]
#which(selc==FALSE)
selc <- colSums(countdata) < 100000000
pheno <- sampleinfo[selc,]
dim(pheno)
#keep<-rowSums(countdata) >= 225

# Which values in myCPM are greater than 0.5?
#thresh <- myCPM > 0.5

# This produces a logical matrix with TRUEs and FALSEs
#head(thresh)

# Summary of how many TRUEs there are in each row
# There are 11406 genes that have TRUEs in all 547 samples.
#table(rowSums(thresh))

# we would like to keep genes that have at least 2 TRUES in each row of thresh
#keep <- rowSums(thresh) >= 2
# Subset the rows of countdata to keep the more highly expressed genes
#counts.keep <- countdata[keep,]
#summary(keep)
#Mode   FALSE    TRUE 
#logical   33186   27266 

dim(counts.keep)

# Let's have a look and see whether our threshold of 0.5 does indeed correspond to a count of about 10-15
# We will look at the first sample
plot(myCPM[,1],countdata[,1])
# Let us limit the x and y-axis so we can actually look to see what is happening at the smaller counts
plot(myCPM[,1],countdata[,1],ylim=c(0,50),xlim=c(0,3))
# Add a vertical line at 0.5 CPM
abline(v=0.5)

#Convert counts to DGEList object
y <- DGEList(counts.keep)
y
names(y)

# Library size information is stored in the samples slot
y$samples

#Library sizes and distribution plots
y$samples$lib.size

# The names argument tells the barplot to use the sample names on the x-axis
# The las argument rotates the axis names
barplot(y$samples$lib.size,names=colnames(y),las=2)
# Add a title to the plot
title("Barplot of library sizes")


## Get log2 counts per million
logcounts <- cpm(y,log=TRUE)
# Check distributions of samples using boxplots
boxplot(logcounts, xlab="Samples", ylab="Log2 counts per million",las=2)
# Let's add a blue horizontal line that corresponds to the median logCPM
abline(h=median(logcounts),col="blue")
title("Boxplots of logCPMs (unnormalised)")

#Multidimensional scaling plots
plotMDS(y)
# How many cell types and in what order are they stored?
levels(sampleinfo$Class)
## Let's choose purple for Aggressive and orange for Indolent and Green for Normal
col.cell <- c("purple","orange","dark green")[sampleinfo$Class]
data.frame(sampleinfo$Class,col.cell)

# Redo the MDS with cell type colouring
plotMDS(y,dim=c(7,8),col=col.cell,pch = 16)
# Let's add a legend to the plot so we know which colours correspond to which cell type
legend("topleft",fill=c("purple","orange","dark green"),legend=levels(sampleinfo$Class))
title("Cell type")

#Another alternative is to generate an interactive MDS plot using the Glimma package. This allows the user to interactively explore the different dimensions.
labels <- paste(sampleinfo$PRAD , sampleinfo$Class )
group <- factor(sampleinfo$Class)
group <- factor(pheno$Class)
glMDSPlot(y, labels=labels, groups=group, folder="mds")

#Hierarchical clustering with heatmaps
# We estimate the variance for each row in the logcounts matrix
var_genes <- apply(logcounts, 1, var)
head(var_genes)

# Get the gene names for the top 500 most variable genes
select_var <- names(sort(var_genes, decreasing=TRUE))[1:500]
head(select_var)

# Subset logcounts matrix
highly_variable_lcpm <- logcounts[select_var,]
dim(highly_variable_lcpm)
head(highly_variable_lcpm)

## Get some nicer colours
mypalette <- brewer.pal(11,"RdYlBu")
morecols <- colorRampPalette(mypalette)
# Set up colour vector for celltype variable
col.cell <- c("blue","red")[sampleinfo$Class]

# Plot the heatmap
# Save the heatmap
png(file="High_var_genes.heatmap.png")
heatmap.2(highly_variable_lcpm,col=rev(morecols(50)),trace="none", main="Top 500 most variable g
enes across samples",ColSideColors=col.cell,scale="row" ,Colv="NA")
dev.off()

#Normalisation for composition bias
# Apply normalisation to DGEList object
y <- calcNormFactors(y, method = "TMM")
y$samples


#Plot the biased and unbiased MD plots side by side for the same sample to see the before and after TMM normalisation effect.
par(mfrow=c(1,2))
plotMD(logcounts,column = 11)
plotMD(y,column = 11)
#save a few data objects
save(group,y,logcounts,sampleinfo,file="day1objects.Rdata")

#Create the design matrix
design <- model.matrix(~ 0 + group)
## Make the column names of the design matrix a bit nicer
colnames(design) <- levels(group)
head(design)
dim(design)
#Voom transform the data
par(mfrow=c(1,1))
v <- voom(y,design,plot = TRUE)
v

#Boxplot the data
par(mfrow=c(1,2))
boxplot(logcounts, xlab="", ylab="Log2 counts per million",las=2,main="Unnormalised logCPM")
## Let's add a blue horizontal line that corresponds to the median logCPM
abline(h=median(logcounts),col="blue")
boxplot(v$E, xlab="", ylab="Log2 counts per million",las=2,main="Voom transformed logCPM")
## Let's add a blue horizontal line that corresponds to the median logCPM
abline(h=median(v$E),col="blue")


write.csv(v$E,file="voom_normalized.csv",row.names=TRUE)

#Testing for differential expression
# Fit the linear model
fit <- lmFit(v)
names(fit)
fit$coefficients
#prepare contrast matrix depending on the class used for D.E. analysis
cont.matrix <- makeContrasts(B.TumorVsNormal=Tumor - Normal,levels=design)
cont.matrix
#comparision
fit.cont <- contrasts.fit(fit, c(1,-1,0))

#performs empirical Bayes shrinkage on the variances, and estimates moderated t-statistics and the associated p-values.
fit.cont <- eBayes(fit.cont)
names(fit.cont)
#We can use the limma decideTests function to generate a quick summary of DE genes for the contrasts.
summa.fit <- decideTests(fit.cont)
summary(summa.fit)
#vennDiagram(summa.fit,include="both")

#The limma topTable function summarises the output in a table format.
topTable(fit.cont,coef=1,sort.by="p")

volcanoplot(fit.cont,coef=1, main="Tumor vs Normal")
nrow(topTable(fit.cont, coef=1,lfc = 1.5,n="Inf"))
probeset.list1 <- topTable(fit.cont, coef=1,n="Inf")
write.table(probeset.list1, file = "AggressiveVsIndolent.txt", sep="\t", quote = F, row.names=T)
probeset.list1 <- topTable(fit.cont, coef=1,n="Inf",lfc = 1.2)
threshold.high <- sort(probeset.list1$logFC, decreasing = TRUE)
threshold.low <- sort(probeset.list1$logFC, decreasing = FALSE)
with(subset(probeset.list1, logFC > threshold.high),
     points(logFC, -log10(P.Value), pch=20, col="red"))

with(subset(probeset.list1, logFC < threshold.high),
     points(logFC, -log10(P.Value), pch=20, col="red"))

View(probeset.list1)
genes <- row.names(probeset.list1)
#genes <- probeset.list1[,1]
norm.data <- v$E[genes,(1:495)]
head(norm.data)
dim(norm.data)
write.csv(norm.data,file="IndvsAgg_lfc_1.2.csv",row.names=TRUE)


#check any genes are duplicated
any(duplicated(rownames(fit.cont)))

#Adding annotation and saving the results using Human annotation database
suppressMessages(library("AnnotationDbi"))
suppressMessages(library("org.Hs.eg.db"))

#check D.E. genes and annotate them
fit.cont$SYMBOL <- mapIds(org.Hs.eg.db,keys=gsub("\\..*","",row.names(fit.cont)),column="SYMBOL",keytype="ENSEMBL",multiVals="first")
names(fit.cont)
nrow(topTable(fit.cont,coef=1,sort.by="p",n="Inf",lfc = 2))
limma.res <- topTable(fit.cont,coef=1,sort.by="p",n="Inf")
#add gene symbols
limma.res$SYMBOL <- mapIds(org.Hs.eg.db,keys=gsub("\\..*","",row.names(limma.res)),column="SYMBOL",keytype="ENSEMBL",multiVals="first")
names(limma.res)

#write the output
write.csv(limma.res,file="tumorvsnormal_results.csv",row.names=TRUE)

#Plots after testing for DE
par(mfrow=c(1,2))
plotMD(fit.cont,coef=1,status=summa.fit[,"B.TumorVsNormal"], values = c(-1, 1))
# For the volcano plot we have to specify how many of the top genes to highlight.
# We can also specify that we want to plot the gene symbol for the highlighted genes.
# let's highlight the top 100 most DE genes
volcanoplot(fit.cont,coef=1,highlight=100,names=limma.res$SYMBOL)

#An interactive version of the volcano plot above that includes the raw per sample values in a separate panel is
#possible via the glXYPlot function in the Glimma package.

levels(group)

glXYPlot(x=fit.cont$coefficients[,1], y=fit.cont$lods[,1],
         xlab="logFC", ylab="-log10(P-value)", main="B.TumorVsNormal",
         counts=y$counts, groups=group, status=summa.fit[,1],
         anno=fit.cont$genes, folder="volcano")

#Testing relative to a threshold (TREAT)
# Let's decide that we are only interested in genes that have a absolute logFC of 1.
# This corresponds to a fold change of 2, or 0.5 (i.e. double or half).
# We can perform a treat analysis which ranks our genes according to p-value AND logFC.
# This is easy to do after our analysis, we just give the treat function the fit.cont object and specify our cut-off.
fit.treat <- treat(fit.cont,lfc=2)
row.names(fit.treat)<-fit.treat$SYMBOL
res.treat <- decideTests(fit.treat)
summary(res.treat)
#B.TumorVsNormal
#Down              7159
#NotSig            8344
#Up               11763
topTable(fit.treat,coef=1,sort.by="p")

# Notice that much fewer genes are highlighted in the MAplot
par(mfrow=c(1,1))
plotMD(fit.treat,coef=1,status=res.treat[,"B.TumorVsNormal"], values=c(-1,1))
abline(h=0,col="grey")

#An interactive version of the mean-difference plots is possible via the glMDPlot function in the Glimma package.
glMDPlot(fit.treat, coef=1, counts=y$counts, groups=group,
         status=res.treat,  main="TumorVsNormal",
         folder="md")


sample.groups <- c(rep("red", 190), rep("blue",304))
sample.groups <- c(rep("red", 190), rep("blue",304), rep("green",52))
sample.groups <- c(rep("red", 495), rep("blue",52))
pc <- princomp(v)
plot(pc$loadings, col=sample.groups, xlab= "Component 1", ylab="Component 2", main = "PCA")
plot3d(pc$loadings, col=sample.groups, xlab= "Component 1", ylab="Component 2", zlab="Component 3", type = "s", size = 1.5)

