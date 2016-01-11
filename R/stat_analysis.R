library("dplyr")
library("ggplot2")

prots <- read.table("halflife/data/protwise.tsv", header=TRUE)

mean_score <- function(comp){return (mean(filter(prots, struc==comp)$s))}
unq_chains <- function(comp){return (mean(filter(prots, struc==comp)$unq))}
tot_chains <- function(comp){return (mean(filter(prots, struc==comp)$tot))}

df <- data.frame(levels(factor(prots$struc)))
df$V2 <- apply(df, 1, mean_score)
colnames(df) <- c("struc", "meanscore")
df$unq <- apply(df, 1, unq_chains)
df$tot <- apply(df, 1, tot_chains)

