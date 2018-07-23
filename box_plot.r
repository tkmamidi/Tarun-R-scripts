library(readr)

sign_exp <- read_csv("O:/Tanja_project/Day_0_3_14/sign_exp_ln2_0_3_14.csv")
rname=sign_exp[,1]
sign_matrix= as.data.frame(sign_exp[,-1])


boxplot(sign_matrix)
stripchart(sign_matrix, vertical = TRUE, 
           method = "jitter", add = TRUE, pch = 20, col = 'blue')

scatter.smooth( sign_matrix$`Cd3-1`~ sign_matrix$`Cd14-1` )
