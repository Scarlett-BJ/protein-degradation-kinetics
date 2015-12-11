library("ggplot2")
library("dplyr")
library("reshape2")

################################################################################

# Build data frames
mouse_file <- "halflife/halflife/data/mouse_paxdb_tissue_expression.txt"
df_mouse <- read.table(mouse_file, header = TRUE)
df_mouse <- melt(df_mouse, id.vars = c("prot", "tcount", "def"),
                 variable.name = "tissue", value.name = "abund")
df_mouse$bin <- cut(log(df_mouse$abund), c(-7, -3, 1, 5, 9, 13))
df_mouse <- na.omit(df_mouse)
df_mouse$tcount <- df_mouse$tcount - 1
df_mouse_tissues <- filter(df_mouse, tissue != "WHOLE_ORGANISM")
df_mouse_tissues <- droplevels(df_mouse_tissues)
df_mouse_whole <- filter(df_mouse, tissue == "WHOLE_ORGANISM")

human_file <- "halflife/halflife/data/human_paxdb_tissue_expression.txt"
df_human <- read.table(human_file, header = TRUE)
df_human <- melt(df_human, id.vars = c("prot", "tcount", "def"),
                 variable.name = "tissue", value.name = "abund")
df_human$bin <- cut(log(df_human$abund), c(-7, -3, 1, 5, 9, 13))
df_human <- na.omit(df_human)
df_human$tcount <- df_human$tcount - 1
df_human_tissues <- filter(df_human, tissue != "WHOLE_ORGANISM" &
                             tissue != "SALIVA" & tissue != "URINE" &
                             tissue != "PLASMA")
df_human_tissues <- droplevels(df_human_tissues)
df_human_whole <- filter(df_human, tissue == "WHOLE_ORGANISM")

human_prodb <- "halflife/halflife/data/human_proteomicsdb_tissue_expression.txt"
df_protdb_human <- read.table(human_prodb, header = TRUE, sep = "\t")
df_protdb_human <- melt(df_protdb_human, id.vars = c("prot", "tcount", "def"),
                        variable.name = "tissue", value.name = "abund")
df_protdb_human <- filter(df_protdb_human, tissue != 'Unknown')
df_protdb_human <- droplevels(df_protdb_human)
tissues <- c()
for (t in levels(df_protdb_human$tissue)){
  corr <- cor.test(filter(df_protdb_human, tissue == t)$abund,
                   filter(df_protdb_human, tissue == t)$tcount)
  if (corr$estimate >= 0.6){
    tissues <- c(tissues, t)
  }
}
df_protdb_human <- filter(df_protdb_human, tissue %in% tissues)
df_protdb_human <- droplevels(df_protdb_human)
df_protdb_human <- na.omit(df_protdb_human)
df_protdb_human$bin <- cut(df_protdb_human$abund, c(1, 2, 3, 4, 5, 6, 7, 8, 9))


# Linear model plots
plot_linear_models <- function(df){
  # df[df$def == "UN",]$def <- "NED"
  plt <- ggplot(df, aes(x = abund, y = tcount))
  cat_colours <- c("dodgerblue", "firebrick", "darkgoldenrod1")
  line_colours <- c("dodgerblue", "firebrick4", "darkorange3")
  index <- 1
  plt_lm <- plt
  for (category in c("UN", "NED", "ED")){
    plt_lm <- plt_lm +
      geom_jitter(data = filter(df, def == category),
                  alpha = 0.35,
                  colour = c(cat_colours[index]))
    index <- index + 1
  }
  index <- 1
  for (category in c("UN", "NED", "ED")){
    plt_lm <- plt_lm +
      geom_smooth(data = filter(df, def == category), 
                  method = "lm",
                  # formula = y ~ s(x, bs = "cs"),
                  colour = c(line_colours[index]),
                  fill = c(line_colours[index]),
                  lwd = 1)
    index <- index + 1
  }
  return (plt_lm)
}
plot_linear_models(filter(df_protdb_human, tissue == "ovary")) + ylim(0, 65)

# Display linear model plots
plot_linear_models(df_protdb_human) + facet_wrap(~ tissue)


plot_linear_models(df_human_tissues) + 
  facet_wrap(~ tissue) +
  ylim(0.5, 19.5)

plot_linear_models(df_mouse_tissues) + 
  facet_wrap(~ tissue) +
  ylim(0.5, 8.5)
  
plot_linear_models(df_human_whole) + 
  ylim(0.5, 19.5)

plot_linear_models(df_mouse_whole) + 
  ylim(0.5, 8.5)

plot_linear_models(df_human_tissues)

plot_linear_models(df_mouse_tissues) +
  ylim(0.5, 8.5)


# Binned boxplots
plot_binned_boxplots <- function(df){
  df[df$def == "UN",]$def <- "NED"
  df <- na.omit(df)
  plt <- ggplot(df, aes(x = bin, y = tcount, fill = def))
  plt_bbox <- plt + 
    geom_boxplot(position = position_dodge(), outlier.size = 0.1,
                 notch = TRUE) +
    scale_fill_manual(values =  c("darkgoldenrod1", "darkred"))
  return(plt_bbox)
}

plot_binned_boxplots(df_protdb_human) + facet_wrap(~ tissue)
plot_binned_boxplots(filter(df_protdb_human, tissue == "ovary"))

# Median boxplot
plot_median_boxplot <- function(df){
  df[df$def == "UN",]$def <- "NED"
  plt <- ggplot(filter(df, abund >= summary(df$abund)[2] &
                         abund <= summary(df$abund)[3]),
                aes(x = def, y = tcount, fill = def)) +
    # geom_violin(adjust = 2) +
    geom_boxplot(notch = T) +
    scale_fill_manual(values =  c("darkgoldenrod1", "darkred")) +
    theme(legend.position = "none")
  print(wilcox.test(filter(df, def == "ED")$abund, 
                    filter(df, def == "NED")$abund))
  print(median(filter(df, def == "NED")$tcount))
  print(median(filter(df, def == "ED")$tcount))
  return(plt)
}
plot_median_boxplot(df_protdb_human) + facet_wrap(~ tissues)

