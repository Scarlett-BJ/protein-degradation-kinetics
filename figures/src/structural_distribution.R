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
df_combined <- merge(df_human, df_mouse, all = TRUE)

## Statistical tests
# Returns fisher.test pval for matrix containing chosen classes
get_fisher <- function(df, matrix){
  mvals <- c()
  for (i in matrix[,1]){
    for (j in matrix[,2]){
      len <- length(filter(df, qtype == i, decay.class == j)$gene)
      mvals <- append(mvals, len) 
    }
  }
  f <- fisher.test(matrix(mvals, nrow = 2, ncol = 2))
  return(signif(f$p.value, 1))
}

# Kolmogorov-Smirnov test, adjusted to account for discrete distribution of subs
dgof_ks <- function(df, c1, c2){
  c1v <- factor(filter(df, qtype == "het", decay.class == c1)$unq)
  c2v <- factor(filter(df, qtype == "het", decay.class == c2)$unq)
  print(ks.test(c1v, ecdf(c2v), simulate.p.value = TRUE, exact = FALSE))
}


## Plot functions
stacked_plotter <- function(df){
  options(scipen = 10000)
  df <- filter(df, decay.class != "X")
  df$decay.class <- factor(df$decay.class, levels = c("E", "U", "N"))
  levels(df$decay.class) <- c("ED", "Undefined", "NED")
  df$qtype <- factor(df$qtype, levels = c("het", "hom", "mon"))
  # get_fisher dataframe must be specified manually for some reason?
  pvals <- c(get_fisher(df_mouse, matrix(c("mon", "het", "N", "U"), nrow = 2)),
             get_fisher(df_mouse, matrix(c("mon", "het", "N", "E"), nrow = 2)))
  cols <- c("#003366", "#6699CC", "#C0C0C0")
  # Must define fill in geom_bar (i.e. locally) for geom_text to work later
  plt <- ggplot(df, aes(x = decay.class)) +
    geom_bar(aes(fill = qtype), position = "fill", 
             colour = c("black"), lwd = 0.4) +
    scale_y_continuous(labels = percent) +
    scale_fill_manual(values = cols) +
    labs(x = "", y = "Subunits", fill = "") +
    theme(text = element_text(size = 10),
          legend.key.size = unit(0.5, "cm")) +
    geom_text(data = data.frame(x = c(2, 1), y = 0.95, lab = pvals),
              aes(x, y, label = lab), size=2)
    # facet_grid(. ~ species) 
    # ggplot ninja...
  return(plt)
}
stacked_plotter(df_mouse)

## Distribution of heteromer subunit counts
density_plotter <- function(df){
  df <- filter(df, qtype == "het", decay.class != "X")
  df$decay.class <- factor(df$decay.class, levels = c("U", "E", "N"))
  levels(df$decay.class) <- c("Undefined", "ED", "NED")
  plt <- ggplot(df, aes(x = unq)) +
    geom_density(aes(fill = decay.class), alpha = 0.7, adjust = 1, lwd = 0.4) +
    scale_fill_manual(values =  c("grey", "#f22a1b", "#2f868b")) +
    scale_x_log10(lim = c(2, 128), 
                  breaks = c(2, 4, 8, 16, 32, 64, 128)) +
    labs(x = "Number of unique subunits", y = "Density", fill = "Decay") +
    theme(text = element_text(size = 10), 
          legend.key.size = unit(0.5, "cm"))
  return(plt)
}

density_plotter(df_human)

