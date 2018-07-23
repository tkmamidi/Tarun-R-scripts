rname=sign_exp[,1]
sign_matrix= as.data.frame(sign_exp[,-1])
row.names(sign_matrix)=rname
View(sign_matrix)
ln2=log2(sign_matrix)
View(ln2)
row.names(as.data.frame(ln2))=rname
write.csv(ln2,file = "C:/Users/tmamid/Desktop/sign_exp.csv")
