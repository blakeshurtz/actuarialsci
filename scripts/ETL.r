###is there a better way to import these libraries? Packrat?
library(readr)
library(lubridate)
library(tidyverse)
library(data.table)
library(ChainLadder)
library(runjags)
library(coda)

setwd("data")

#models are run on a per-insurer basis. Define variables for industry and group.
#See Appendix A for groups included in model
insurer.data="comauto_pos.csv"
grpcode="353"

#import data
a=read.csv(insurer.data)

setwd("..")

