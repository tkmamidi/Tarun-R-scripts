v <- read.csv(file = "voom_normalized.csv", row.names = 1,check.names=FALSE)
sample.groups <- c(rep("red", 45), rep("blue",204))
pc <- princomp(training)
plot(pc$loadings, col=sample.groups, xlab= "Component 1", ylab="Component 2", main = "PCA")
plot3d(pc$loadings, col=sample.groups, xlab= "Component 1", ylab="Component 2", zlab="Component 3", type = "s", size = 1.5)

testing <- read.csv(file = "testing_pca.csv", row.names = 1,check.names=FALSE)
sample.groups <- c(rep("red", 145), rep("blue",101))
pc1 <- princomp(training)
plot(pc1$loadings, col=sample.groups, xlab= "Component 1", ylab="Component 2", main = "PCA")
plot3d(pc1$loadings, col=sample.groups, xlab= "Component 1", ylab="Component 2", zlab="Component 3", type = "s", size = 1.5)
legend("topleft",fill=c("red","blue"),legend=c("Indolent","Aggressive"))
pc1$loadings

select_var <-row.names(probeset.list1)[1:550]
top_500<- hdat[select_var,]
dim(top_500)
head(select_var)
View(top_500)

sample.groups <- c(rep("red", 190), rep("blue",304), rep("green",52))
pc <- princomp(top_500)
plot(pc$loadings, col=sample.groups, xlab= "Component 1", ylab="Component 2", main = "PCA")
plot3d(pc$loadings, col=sample.groups, xlab= "Component 1", ylab="Component 2", zlab="Component 3", type = "s", size = 1.5)
write.csv(top_500,file="top_550_DE.csv",row.names=TRUE)
