library("ggplot2")
library("dplyr")
library("scales")
library("gridExtra")
library("RColorBrewer")
library("dgof")

humanf <- "halflife/halflife/data/structural/NED_quaternary_human.txt"
mousef <- "halflife/halflife/data/structural/NED_quaternary_mouse.txt"
df_human <- read.table(humanf, header = TRUE)
df_mouse <- read.table(mousef, header = TRUE)
df_mouse <- mutate(df_mouse, species = "mouse")
df_human <- mutate(df_human, species = "human")
df_combined <- merge(df_human, df_mouse, all = TRUE)

# Returns fisher.test pval for matrix containing chosen classes
get_fisher <- function(df, matrix){
  mvals <- c()
  for (i in matrix[,1]){
    for (j in matrix[,2]){
      len <- length(filter(df, qtype == i, decay.class == j)$gene)
      print(c(i, j, len))
      mvals <- append(mvals, len) 
    }
  }
  f <- fisher.test(matrix(mvals, nrow = 2, ncol = 2))
  return(signif(f$p.value, 2))
}

stacked_plotter <- function(df){
  df <- filter(df, decay.class != "X")
  df$decay.class <- factor(df$decay.class, levels = c("N", "U", "E"))
  levels(df$decay.class) <- c("NED", "UN", "ED")
  df$qtype <- factor(df$qtype, levels = c("mon", "hom", "het"))
  # levels(df$qtype) <- c("monomer", "homomer", "heteromer")
  cols <- colorRampPalette(brewer.pal(8, "RdYlBu"))(8)
  # Must define fill in geom_bar (i.e. locally) for geom_text to work later
  plt <- ggplot(df, aes(x = decay.class)) +
    geom_bar(aes(fill = qtype), position = "fill", 
             colour = c("black"), lwd = 0.4) +
    scale_y_continuous(labels = percent) +
    scale_fill_manual(values = cols) +
    labs(x = "", y = "Genes", fill = "") +
    theme(text = element_text(size = 10),
          legend.key.size = unit(0.5, "cm")) +
    facet_grid(. ~ species) 
    # ggplot ninja...
  return(plt)
}

## Annotate plot
lab1 <- c(get_fisher(df_human, matrix(c("mon", "het", "N", "U"), nrow = 2)),
          get_fisher(df_mouse, matrix(c("mon", "het", "N", "U"), nrow = 2)))
lab2 <- c(get_fisher(df_human, matrix(c("mon", "het", "N", "E"), nrow = 2)),
          get_fisher(df_mouse, matrix(c("mon", "het", "N", "E"), nrow = 2)))

stacked <- stacked_plotter(df_combined) +
geom_text(aes(x, y, label = lab),
          data=data.frame(x = 2, y = 0.95, lab = lab1,
                          species = c("human", "mouse")), size = 2) +
  geom_text(aes(x, y, label = lab),
            data=data.frame(x = 3, y = 0.95, lab = lab2,
                            species = c("human", "mouse")), size = 2)
stacked

density_plotter <- function(df){
  df <- filter(df, qtype == "het", decay.class != "X")
  df$decay.class <- factor(df$decay.class, levels = c("U", "E", "N"))
  levels(df$decay.class) <- c("UN", "ED", "NED")
  plt <- ggplot(df, aes(x = unq)) +
    geom_density(aes(fill = decay.class), alpha = 0.7, adjust = 1, lwd = 0.4) +
    scale_fill_manual(values =  c("grey", "#f22a1b", "#2f868b")) +
    scale_x_log10(lim = c(2, 128), 
                  breaks = c(2, 4, 8, 16, 32, 64, 128)) +
    labs(x = "Number of unique subunits", y = "Density", fill = "") +
    theme(text = element_text(size = 10), 
          legend.key.size = unit(0.5, "cm")) +
    facet_grid(. ~ species)
  return(plt)
}

dgof_ks <- function(df, c1, c2){
  c1v <- factor(filter(df, qtype == "het", decay.class == c1)$unq)
  c2v <- factor(filter(df, qtype == "het", decay.class == c2)$unq)
  print(ks.test(c1v, ecdf(c2v), simulate.p.value = TRUE, exact = FALSE))
}
dgof_ks(df_human, "N", "E")

stacked 
density_plotter(df_combined)
