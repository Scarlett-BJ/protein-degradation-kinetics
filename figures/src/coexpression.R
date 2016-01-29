library("ggplot2")
library("dplyr")

df_corum <- read.table("halflife/halflife/data/coexpressdb_corum_mouse_homologs.tsv",
                       header = TRUE, sep = "\t")
colnames(df_corum) <- c("complex.id", "usubs", "uniprot.id",
                        "avg.coexpression", "decay.class", "species")

density_plotter <- function(df){
  df$bin <- cut(df$usubs, c(1, 5, 10, 15, 20, 25, 30, 40, 80, 200))
  levels(df$bin) <- c("2-5", "6-10", "11-15", "16-20", "21-25", "26-30",
                      "31-40", "41-80", "81+")
  df <- na.omit(df)
  df$decay.class <- factor(df$decay.class, levels = c("UN", "ED", "NED"))
  levels(df$decay.class) <- c("Undefined", "ED", "NED")
  plt <- ggplot(df, aes(x = avg.coexpression, fill = decay.class)) +
    geom_density(alpha = 0.7) +
    scale_fill_manual(values = c("grey", "#f22a1b", "#2f868b")) +
    labs(x = "Average subunit coexpression (Pearson's r)",
         y = "Density",
         fill = "Decay") +
    theme(text = element_text(size = 10),
          legend.key.size = unit(0.5, "cm"))
  return(plt)
}

density_plotter(df_corum) + facet_wrap(~ bin)


# Deprecated, Correct but not for publication
box_plotter <- function(df){
  df <- filter(df, !is.na(decay.class))
  df$decay.class <- factor(df$decay.class, levels = c("NED", "UN", "ED"))
  levels(df$decay.class) <- c("NED", "Undefined", "ED")
  plt <- ggplot(df, aes(x = decay.class, y = avg.coexpression,
                        fill = decay.class)) +
    geom_boxplot(notch = TRUE, outlier.size = 0.75) +
    scale_fill_manual(values = c("#2f868b", "grey", "#f22a1b")) +
    labs(x = "", y = "Average subunit coexpression (Pearson's r)") +
    theme(text = element_text(size = 10),
          legend.position = "none")
  return(plt)
}
box_plotter(df_corum)
