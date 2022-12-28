
setwd("D:/Desktop/R包开发/LXPCoA")
library(openxlsx)

LXPCoA_data_example <- read.xlsx("LXPCoA_data_example.xlsx")

usethis::use_data(LXPCoA_data_example,overwrite = T)

rm(list=ls())

data(LXPCoA_data_example)

