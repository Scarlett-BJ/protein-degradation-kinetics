library("ggplot2")
library("dplyr")

df_corum <- read.table("halflife/halflife/data/coexpressdb_corum_combined.tsv")
colnames(df_corum) <- c("complex.id", "entrez.id", "uniprot.id",
                        "avg.coexpression", "decay.class", "species")

density_plotter <- function(df){
  df <- df %>% 
    group_by(complex.id) %>%
    mutate(unq = n())
  df$bin <- cut(df$unq, c(0, 5, 10, 15, 20, 25, 30, 40, 80, 120))
  levels(df$bin) <- c("0-5", "5-10", "10-15", "15-20", "20-25", "25-30",
                      "30-40", "40-80", "80-120")
  df <- na.omit(df)
  df$decay.class <- factor(df$decay.class, levels = c("UN", "ED", "NED"))
  plt <- ggplot(df, aes(x = avg.coexpression, fill = decay.class)) +
    geom_density(alpha = 0.7) +
    scale_fill_manual(values = c("cadetblue4", "darkgoldenrod1", "darkred")) +
    labs(x = "Average subunit coexpression (Pearson's r)",
         y = "Density",
         fill = "Decay") +
    theme(text = element_text(size = 10),
          legend.key.size = unit(0.5, "cm"))
  return(plt)
}
density_plotter(df_corum) + facet_wrap(~ bin)


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
  df <- filter(df, !is.na(decay.class))
  df$decay.class <- factor(df$decay.class, levels = c("NED", "UN", "ED"))
  plt <- ggplot(df, aes(x = decay.class, y = avg.coexpression,
                        fill = decay.class)) +
    geom_boxplot(notch = TRUE, outlier.size = 0.75) +
    scale_fill_manual(values = c("darkred", "cadetblue4", "darkgoldenrod1")) +
    labs(x = "", y = "Average subunit coexpression (Pearson's r)") +
    theme(text = element_text(size = 10),
          legend.position = "none")
  return(plt)
}
box_plotter(df_corum)
