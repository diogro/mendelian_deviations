library(gdata)
library(plyr)
library(dplyr)
library(reshape2)

raw_mouse_gen = llply(paste0("./data/genotypes/chrom", 1:19, ".csv"), read.csv, as.is = TRUE)
names(raw_mouse_gen) = paste0("chrom", 1:19)
mouse_gen = raw_mouse_gen

rm(list = ls(pattern='raw'))
