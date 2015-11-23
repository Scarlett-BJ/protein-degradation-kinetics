library("dplyr")
library("ggplot2")
library("reshape2")

prots <- read.table("halflife/data/protwise.tsv", header=TRUE)
df <- melt(prots)
