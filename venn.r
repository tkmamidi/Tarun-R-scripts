install.packages('VennDiagram')
library(VennDiagram)
library(limma)


head(venn)
nrow(subset(venn, Somatic >1))

grid.newpage()
draw.triple.venn(area1 =214, area2 = 12879, area3 = 12527, n12 = 169, n23 =6897, n13 = 153, n123 = 125,
category = c("GWAS", "Somatic", "RNA-seq"), lty = "blank", 
fill = c("blue", "red", "green"),cex = 2,cat.cex = 2,
cat.col = c("blue", "red", "green"), scaled = FALSE)

?draw.triple.venn
