
setwd("D:/Desktop/R包开发/LXPCoA")
library(openxlsx)

PCoA_data_example <- read.xlsx("PCoA_example.xlsx")

usethis::use_data(PCoA_data_example,overwrite = T)

rm(list=ls())

data(PCoA_data_example)

