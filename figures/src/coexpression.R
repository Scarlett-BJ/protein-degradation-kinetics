library("ggplot2")
library("dplyr")

mouse_file <- "halflife/halflife/data/coexpression/coexpressdb_corum_mouse_homologs.tsv"
human_file <- "halflife/halflife/data/coexpression/coexpressdb_corum_human.tsv"
header <- c("cid", "usubs", "uniprot.id", "avg.coex", "decay.class", "species")
df_corum_mouse <- read.table(mouse_file, header = TRUE, sep = "\t")
df_corum_human <- read.table(human_file, header = TRUE, sep = "\t")
colnames(df_corum_mouse) <- header
colnames(df_corum_human) <- header
df_corum_mouse$decay.class <- factor(df_corum_mouse$decay.class, 
                                     levels = c("UN", "ED", "NED"))
levels(df_corum_mouse$decay.class) <- c("Undefined", "ED", "NED")

df_struc <- read.table("halflife/halflife/data/coexpression/struc_coex_filt.txt",
                       header = TRUE, sep = "\t")
colnames(df_struc) <- c("gene", "decay.class", "cid", "avg.coex", "seqid", "usubs")
df_struc$decay.class <- factor(df_struc$decay.class, 
                               levels = c("U", "E", "N"))
levels(df_struc$decay.class) <- c("Undefined", "ED", "NED")


density_plotter <- function(df){
  df$bin <- cut(df$usubs, c(2, 5, 10, 15, 20, 25, 30, 40, 80, 200))
  levels(df$bin) <- c("3-5", "6-10", "11-15", "16-20", "21-25", "26-30",
                      "31-40", "41-80", "81+")
#   df$bin <- cut(df$usubs, c(2, 4, 6,  10, 100))
#   levels(df$bin) <- c("3-4", "5-6", "7-10", "11+")
  df <- na.omit(df)
  # df$decay.class <- factor(df$decay.class, levels = c("UN", "ED", "NED"))
  # levels(df$decay.class) <- c("Undefined", "ED", "NED")
  plt <- ggplot(df, aes(x = avg.coex)) +
    geom_density(aes(fill = decay.class), alpha = 0.7, lwd = 0.4) +
    scale_fill_manual(values = c("grey", "#f22a1b", "#2f868b")) +
    labs(x = "Average subunit coexpression (Pearson's r)",
         y = "Density",
         fill = "") +
    scale_x_continuous(limits = c(-0.25, 0.9),
                       breaks = c(-0.2, 0, 0.2, 0.4, 0.6, 0.8),
                       labels = c("-0.2", "0.0", "0.2", "0.4", "0.6", "0.8")) +
    theme(text = element_text(size = 10),
          legend.key.size = unit(0.5, "cm"),
          legend.title = element_blank(),
          legend.margin = unit(0, "cm")) +
    guides(fill=guide_legend(reverse=TRUE))
  return(plt)
}
paneld <- density_plotter(filter(df_struc, usubs >= 3)) +
  geom_text(data = data.frame(x = 0.5, y = 2.25, lab = "P = 1.2e-08"), 
            aes(x, y, label = lab), size = 2.5, fontface = "italic")
density_plotter(filter(df_corum_mouse, usubs >= 3)) + 
  geom_text(data = data.frame(x = 0.45, y = 2.5, lab = "P < 2.2e-16"), 
            aes(x, y, label = lab), size = 2.5)
density_plotter(df_corum_mouse) + facet_wrap(~ bin)
density_plotter(df_corum_human) + facet_wrap(~ bin)
density_plotter(df_struc) + facet_wrap(~ bin)
