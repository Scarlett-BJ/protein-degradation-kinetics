library("ggplot2")
library("dplyr")

df_corum <- read.table("halflife/halflife/data/coexpressdb_corum_combined.tsv")
colnames(df_corum) <- c("complex.id", "entrez.id", "uniprot.id",
                        "avg.coexpression", "decay.class", "species")

density_plotter <- function(df){
  df <- filter(df, !is.na(decay.class))
  df$decay.class <- factor(df$decay.class, levels = c("UN", "ED", "NED"))
  plt <- ggplot(df, aes(x = avg.coexpression, fill = decay.class)) +
    geom_density(alpha = 0.7) +
    scale_fill_manual(values = c("cadetblue4", "darkgoldenrod1", "darkred")) +
    xlab("Average coexpression (Pearson's r)") +
    ylab("Density") +
    labs(fill = "Decay class")
  return(plt)
}
density_plotter(df_corum)

bin_plotter <- function(df, interval){
  df_counts <- count(df, complex.id)
  df_counts <- mutate(df_counts,
                      intervals = findInterval(df_counts$n,
                                               seq(1, max(df_counts$n), 3)))
  comps <- filter(df_counts, intervals == interval)$complex.id
  df_bin <- filter(df, complex.id %in% comps)
  plt <- density_plotter(df_bin)
  return(plt)
}
bin_plotter(df_corum, 10)

box_plotter <- function(df){
  # df <- filter(df, !is.na(decay.class))
  df$decay.class <- factor(df$decay.class, levels = c("NED", "UN", "ED"))
  plt <- ggplot(df, aes(x = decay.class, y = avg.coexpression,
                        fill = decay.class)) +
    geom_boxplot(notch = TRUE, outlier.size = 0.75) +
    scale_fill_manual(values = c("darkred", "cadetblue4", "darkgoldenrod1")) +
    xlab("Decay class") +
    ylab("Average coexpression (Pearson's r)") +
    theme(legend.position = "none")
  return(plt)
}
box_plotter(df_corum)
