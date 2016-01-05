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
tissues <- c()
for (t in levels(df_protdb_human$tissue)){
  if (var(filter(df_protdb_human, tissue == t)$tcount) > 150){
    tissues <- c(tissues, t)   
  }
}
df_protdb_human <- filter(df_protdb_human, tissue %in% tissues)
df_protdb_human <- droplevels(df_protdb_human)
df_protdb_human <- na.omit(df_protdb_human)
df_protdb_human$bin <- cut(df_protdb_human$abund, c(1, 2, 3, 4, 5, 6, 7, 8, 9))

# Binned boxplots
plot_binned_boxplots <- function(df){
  # df[df$def == "UN",]$def <- "NED"
  # df <- filter(df, def != "UN")
  df <- na.omit(df)
  df <- mutate(df, tbin = cut(df$tcount, c(1, 10, 20, 30, 40, 50, 60)),
               abin = cut(df$abund, c(1, 2, 3, 4, 5, 6, 7, 8, 9)))
  plt <- ggplot(df, aes(y = tcount, x = abin, fill = def))
  plt_bbox <- plt + 
    geom_boxplot(position = position_dodge(), outlier.size = 0.1,
                 notch = FALSE) +
#     scale_x_discrete(labels = c("2 - 3", "3 - 4", "4 - 5", 
#                                 "5 - 6", "6 - 7", "7 - 8")) +
    scale_fill_manual(values =  c("darkgoldenrod1", "darkred", "cadetblue4")) +
#     ylab("Tissues in which protein is expressed") +
#     xlab("Binned protein abundance, (log10 iBAQ intensity)") +
    labs(fill = "Decay class")
  return(plt_bbox)
}

plot_binned_boxplots(df_protdb_human)
plot_binned_boxplots(filter(df_protdb_human, tissue == "ovary"))

# Median boxplot
plot_median_boxplot <- function(df){
  df[df$def == "UN",]$def <- "NED"
  df <- filter(df, def != "UN")
  quartile1 <- summary(df$abund)[2]
  median <- summary(df$abund)[3]
  quartile3 <- summary(df$abund)[5]
  df <- filter(df, abund >= quartile1 & abund <= quartile3)
  plt <- ggplot(df,aes(x = def, y = tcount, fill = def)) +
    # geom_violin(adjust = 2) +
    geom_boxplot(notch = T) +
    scale_fill_manual(values =  c("darkgoldenrod1", "darkred")) +
    theme(legend.position = "none")
  print(wilcox.test(filter(df, def == "ED")$tcount, 
                    filter(df, def == "NED")$tcount))
  print(mean(filter(df, def == "NED")$tcount))
  print(mean(filter(df, def == "ED")$tcount))
  return(plt)
}
plot_median_boxplot(filter(df_protdb_human, tissue == "ovary"))

# Linear model
plot_linear <- function(df){
  df[df$def == "UN",]$def <- "NED"
  df <- droplevels(df)
  df$def <- factor(df$def, levels = c("ED", "NED"))
  df
  plt <- ggplot(df, aes(x = abund, y = tcount, colour = def, fill = def)) +
    geom_jitter(alpha = 0.35, size = 1.5) +
    geom_smooth(lwd = 1.2, method = "lm") +
    scale_colour_manual(values = c("darkorange", "firebrick"),
                        labels = c("ED", "NED or UN")) +
    scale_fill_manual(values = c("darkgoldenrod1", "firebrick"),
                      labels = c("ED", "NED or UN")) +
    # ylim(c(0, 65)) +
#     ylab("Tissues in which protein is expressed") +
#     xlab("Protein abundance, (log10 iBAQ intensity)") +
    labs(fill = "Decay class", colour = "Decay class")
  return(plt)
}
plot_linear(filter(df_protdb_human, tissue == "ovary"))

# Simple Density plot
tissue_count_density <- function(df){
  df$def <- factor(df$def, levels = c("NED", "UN", "ED"))
  plt <- ggplot(df, aes(x = tcount, fill = def)) +
    geom_density(alpha = 0.7, adjust = 2) +
    scale_fill_manual(values = c("darkred", "cadetblue4", "darkgoldenrod1")) +
    xlab("Number of tissues in which protein is expressed") +
    ylab("Density") +
    labs(fill = "Decay class")
  return(plt)
}
tissue_count_density(df_protdb_human)