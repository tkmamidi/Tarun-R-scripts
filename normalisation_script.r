install.packages("combinat")
library("combinat")
input = read.delim("matrix_GE.txt", header = TRUE)
ge_matrix = as.matrix(hdat)

library(preprocessCore)
ge_nml = normalize.quantiles(ge_matrix, copy=FALSE )
max(ge_nml, na.rm=TRUE)

 tum.ge<-data.matrix(ge_nml[1:17813, 1:155])
 num.ge<-data.matrix(ge_nml[1:17813, 156:173])
 wil.ge <- matrix(nrow = 17814, ncol = 1)
rownames(wil.ge) = rownames(ge_matrix)
 for(i in 1:17813)
 {
   r1=as.numeric(tum.ge[i,])
   r2=as.numeric(num.ge[i,])
   wil.ge[i,1] = wilcox.test(r1, r2)$p.value
 }
 write.table(wil.ge, file="wilcox_GE_output.txt", quote=F, sep="\t")
 
 ge.fdr = p.adjust(wil.ge,method="fdr")
 gene.expression_fdr = as.matrix(ge.fdr)
 rownames(gene.expression_fdr)= rownames(wil.ge)
 write.table(gene.expression_fdr, file="fdr_GE_output.txt", quote=F, sep="\t")

 
 dim(ge_matrix)
 dim(ge_nml)
 
 
 
 edata=log2(hdat+1)
 edata=edata[rowMeans(edata)>3,]
 
 
 
 colramp = colorRampPalette(c(3,"white",2))(20)
 plot(density(ge_nml[,1]),col=colramp[1],lwd=3,ylim=c(0,.20))
 for(i in 2:20){lines(density(ge_nml[,i]),lwd=3,col=colramp[i])}
 
 plot(ge_nml[1,],col=as.numeric(pheno$class))
 
 
 svd1 = svd(ge_nml - rowMeans(ge_nml))
 plot(svd1$v[,1],svd1$v[,2],xlab="PC1",ylab="PC2",
      col=as.numeric(pheno$class))
 
 
 