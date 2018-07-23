library(readr)
whole_data <- read_csv("Z:/prostate_cancer_new/whole_data/7.csv")
t(whole_data[22,2])
sorted=whole_data[order(t(whole_data[22,]))]
View(sorted)
write.csv(sorted,"Z:/prostate_cancer_new/whole_data/sorted_7.csv")
