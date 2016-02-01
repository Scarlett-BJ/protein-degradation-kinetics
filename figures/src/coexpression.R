library("ggplot2")
library("dplyr")

mouse_file <- "halflife/halflife/data/coexpressdb_corum_mouse_homologs.tsv"
human_file <- "halflife/halflife/data/coexpressdb_corum_human.tsv"
header <- c("cid", "usubs", "uniprot.id", "avg.coex", "decay.class", "species")
df_corum_mouse <- read.table(mouse_file, header = TRUE, sep = "\t")
df_corum_human <- read.table(human_file, header = TRUE, sep = "\t")
colnames(df_corum_mouse) <- header
colnames(df_corum_human) <- header

density_plotter <- function(df){
  df$bin <- cut(df$usubs, c(2, 5, 10, 15, 20, 25, 30, 40, 80, 200))
  levels(df$bin) <- c("3-5", "6-10", "11-15", "16-20", "21-25", "26-30",
                      "31-40", "41-80", "81+")
  df <- na.omit(df)
  df$decay.class <- factor(df$decay.class, levels = c("UN", "ED", "NED"))
  levels(df$decay.class) <- c("Undefined", "ED", "NED")
  plt <- ggplot(df, aes(x = avg.coex, fill = decay.class)) +
    geom_density(alpha = 0.7, lwd = 0.4) +
    scale_fill_manual(values = c("grey", "#f22a1b", "#2f868b")) +
    labs(x = "Average subunit coexpression (Pearson's r)",
         y = "Density",
         fill = "Decay") +
    theme(text = element_text(size = 10),
          legend.key.size = unit(0.5, "cm"))
  return(plt)
}

density_plotter(df_corum_mouse) + facet_wrap(~ bin)
density_plotter(df_corum_human) + facet_wrap(~ bin)

# # Deprecated, Correct but not for publication
# box_plotter <- function(df){
#   df <- filter(df, !is.na(decay.class))
#   df$decay.class <- factor(df$decay.class, levels = c("NED", "UN", "ED"))
#   levels(df$decay.class) <- c("NED", "Undefined", "ED")
#   plt <- ggplot(df, aes(x = decay.class, y = avg.coexpression,
#                         fill = decay.class)) +
#     geom_boxplot(notch = TRUE, outlier.size = 0.75) +
#     scale_fill_manual(values = c("#2f868b", "grey", "#f22a1b")) +
#     labs(x = "", y = "Average subunit coexpression (Pearson's r)") +
#     theme(text = element_text(size = 10),
#           legend.position = "none")
#   return(plt)
# }
# box_plotter(df_corum)
