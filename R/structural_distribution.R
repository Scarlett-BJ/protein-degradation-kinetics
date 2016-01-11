library("ggplot2")
library("dplyr")
library("scales")
library("gridExtra")
library("RColorBrewer")

humanf <- "halflife/halflife/data/ned_structural_files/NED_quaternary_human.txt"
mousef <- "halflife/halflife/data/ned_structural_files/NED_quaternary_mouse.txt"
df_human <- read.table(humanf, header = TRUE)
df_mouse <- read.table(mousef, header = TRUE)

df_mouse <- mutate(df_mouse, species = "mouse")
df_human <- mutate(df_human, species = "human")
df_combined <- merge(df_human, df_mouse, all = TRUE)


stacked_plotter <- function(df){
  df$decay.class <- factor(df$decay.class, levels = c("N", "U", "E"))
  levels(df$decay.class) <- c("NED", "UN", "ED")
  df$qtype <- factor(df$qtype, levels = c("mon", "hom", "het"))
  levels(df$qtype) <- c("monomer", "homomer", "heteromer")
  cols <- colorRampPalette(brewer.pal(8, "RdYlBu"))(8)
  plt <- ggplot(df, aes(x = df$decay.class, fill = df$qtype)) +
    geom_bar(position = "fill") +
    scale_y_continuous(labels = percent) +
    scale_fill_manual(values = cols) +
    xlab("Decay class") +
    ylab("Genes") +
    labs(fill = "") +
    facet_grid(. ~ species)
  return(plt)
}
stacked_plotter(df_combined)
grid.arrange(stacked_plotter(df_mouse), stacked_plotter(df_human), ncol = 2)
