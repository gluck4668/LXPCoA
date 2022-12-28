\name{LXPCoA::LXPCoA}
\alias{LXPCoA}
\title{LXPCoA::principal co-ordinates analysis}
\usage{
LXPCoA (PCoA_data)
}
\description{
LXPCoA can be used for the principal co-ordinates analysis.
}
\examples{

if(!requireNamespace("devtools"))
  install.packages("devtools")
library(devtools)

install_github("gluck4668/LXPCoA")

library(LXPCoA)
??LXPCoA

#-----------------------------
data(LXPCoA_data_example)

#----------------------------


rm(list=ls())

setwd("D:/Desktop/R_example/LXPCoA")

PCoA_data = "PCoA_example.xlsx"


LXPCoA (PCoA_data)
}
