library("ggplot2")
library("dplyr")
library("reshape2")
library("dgof")
library(gridExtra)
################################################################################

## Build data frames

# PaxDB data
paxdb_dataframes <- function(df){
  df <- melt(df, id.vars = c("prot", "tcount", "def"),
                   variable.name = "tissue", value.name = "abund")
  df$bin <- cut(log(df$abund), c(-7, -3, 1, 5, 9, 13))
  df <- na.omit(df)
  df$tcount <- df$tcount - 1
  df_tissues <- droplevels(filter(df, tissue != "WHOLE_ORGANISM"))
  df_whole <- droplevels(filter(df, tissue == "WHOLE_ORGANISM"))
  return(list(df_tissues, df_whole))
}

mouse_file <- "halflife/halflife/data/abundance/mouse_paxdb_tissue_expression.txt"
human_file <- "halflife/halflife/data/abundance/human_paxdb_tissue_expression.txt"
dfv <- paxdb_dataframes(read.table(human_file, header = TRUE))
df_human_tissues <- dfv[[1]]
df_human_whole <- dfv[[2]]
dfv <- paxdb_dataframes(read.table(mouse_file, header = TRUE))
df_mouse_tissues <- dfv[[1]]
df_mouse_whole <- dfv[[2]]
remove(dfv)

# Human ProteomicsDB data
human_prodb <- "halflife/halflife/data/abundance/human_proteomicsdb_tissue_expression.txt"
df_protdb_human <- read.table(human_prodb, header = TRUE, sep = "\t")
df_protdb_human <- melt(df_protdb_human, id.vars = c("prot", "tcount", "def"),
                        variable.name = "tissue", value.name = "abund")
df_protdb_human <- filter(df_protdb_human, tissue != 'Unknown')
df_protdb_human <- droplevels(df_protdb_human)
# tissues <- c()
# for (t in levels(df_protdb_human$tissue)){
#   corr <- cor.test(filter(df_protdb_human, tissue == t)$abund,
#                    filter(df_protdb_human, tissue == t)$tcount)
#   if (corr$estimate >= 0.6){
#     tissues <- c(tissues, t)
#   }
# }
# df_protdb_human <- filter(df_protdb_human, tissue %in% tissues)
# df_protdb_human <- droplevels(df_protdb_human)
# tissues <- c()
# for (t in levels(df_protdb_human$tissue)){
#   if (var(filter(df_protdb_human, tissue == t)$tcount) > 150){
#     tissues <- c(tissues, t)   
#   }
# }
# df_protdb_human <- filter(df_protdb_human, tissue %in% tissues)
df_protdb_human <- droplevels(df_protdb_human)
df_protdb_human <- na.omit(df_protdb_human)

# Tissue lists ranked by quantity of available data
protdb_ranked <- arrange(count(df_protdb_human, tissue), desc(n))$tissue
paxdb_human_ranked <- arrange(count(df_human_tissues, tissue), desc(n))$tissue
paxdb_mouse_ranked <- arrange(count(df_mouse_tissues, tissue), desc(n))$tissue

################################################################################


# Tissue count density plot
tissue_count_density <- function(df){
  # Count each protein once, even if present in many tissues
  df <- df[match(unique(df$prot), df$prot),]
  df$def <- factor(df$def, levels = c("UN", "ED", "NED"))
  levels(df$def) <- c("Undefined", "ED", "NED")
  plt <- ggplot(df, aes(x = tcount, fill = def)) +
    geom_density(alpha = 0.7, lwd = 0.4) +
    scale_fill_manual(values = c("grey", "#f22a1b", "#2f868b")) +
    labs(x = "Number of tissues in which protein is detected",
         y = "Density",
         fill = "Decay") +
    theme(text = element_text(size = 10),
          legend.key.size = unit(0.5, "cm"))
  return(plt)
}
# tissue_count_density(filter(df_protdb_human, tissue == "ovary"))
tissue_count_density(df_protdb_human)
################################################################################

## Linear Model Plots

# PaxDB plots
paxdb_linear <- function(df){
  df <- na.omit(df)
  df <- filter(df, def != "UN")
#   quartile1 <- summary(df$abund)[2]
#   quartile3 <- summary(df$abund)[5]
#   df <- filter(df, abund >= quartile1 & abund <= quartile3)
  plt <- ggplot(df, aes(x = log(abund), y = tcount, colour = def, fill = def)) +
    geom_jitter(alpha = 0.35, size = 1.5) +
    geom_smooth(lwd = 1.2, method = "lm") +
    scale_colour_manual(values = c("darkorange", "firebrick", "cadetblue4"),
                        labels = c("ED", "NED",  "UN")) +
    scale_fill_manual(values = c("darkgoldenrod1", "firebrick", "cadetblue4"),
                      labels = c("ED", "NED", "UN"))
  return(plt)
}

# ProteomicsDB plots
protdb_linear <- function(df){
  #   df[df$def == "UN",]$def <- "NED"
  #   df <- droplevels(df)
  # df <- filter(df, def != "UN")
  stdev <- sqrt(var(df$abund))
  amean <- mean(df$abund)
  minval <- amean - 2 * stdev
  maxval <- amean + 2 * stdev
  df <- filter(df, abund >= minval & abund <= maxval)
  plt <- ggplot(df, aes(x = abund, y = tcount, colour = def, fill = def)) +
    geom_jitter(alpha = 0.35, size = 0.5) +
    geom_smooth(lwd = 1) +
    scale_colour_manual(values = c("darkorange", "firebrick", "cadetblue4"),
                        labels = c("ED", "NED",  "UN")) +
    scale_fill_manual(values = c("darkgoldenrod1", "firebrick", "cadetblue4"),
                      labels = c("ED", "NED", "UN")) +
    ylim(c(0, 65)) +
    labs(x = expression(paste("Protein abundance ", 
                              log[10], "(iBAQ intensity)")),
         y = "Tissues in which protein is expressed",
         fill = "Decay", colour = "Decay")+
    theme(text = element_text(size = 10),
          legend.key.size = unit(0.5, "cm"))
  return(plt)
}
protdb_linear(df)
# Display Plots
# paxdb_linear(filter(df_mouse_tissues, tissue %in% paxdb_mouse_ranked[1:8])) +
#   facet_wrap(~ tissue)
# paxdb_linear(filter(df_human_tissues, tissue %in% paxdb_human_ranked[1:4])) +
#   facet_wrap(~ tissue)
protdb_linear(filter(df_protdb_human, tissue %in% protdb_ranked[1:6])) +
  facet_wrap(~ tissue)

################################################################################

## Boxplots

# ProteomicsDB binned boxplots - Not so nice, can't see distribution properly


bin_wilcoxon <- function(df){
  df <- na.omit(df)
  df <- mutate(df, tbin = cut(df$tcount, c(1, 20, 30, 40, 50, 70)))
  pvals <- c()
  for (bin in levels(df$tbin)){
    ned <- filter(df, tbin == bin, def == "NED")$abund
    ed <- filter(df, tbin == bin, def == "ED")$abund
    p <- wilcox.test(ned, ed)$p.value
    pvals <- append(pvals, as.character(signif(p, 2)))
  }
  return(pvals)
}

protdb_binned_boxplots <- function(df){
  df <- na.omit(df)
  df$def <- factor(df$def, levels = c("UN", "ED", "NED"))
  levels(df$def) <- c("Undefined", "ED", "NED")
  df <- mutate(df, tbin = cut(df$tcount, c(1, 20, 30, 40, 50, 70)))
               # abin = cut(df$abund, c(1, 3, 4, 5, 6, 7, 9)))
  pvals <- bin_wilcoxon(df)
  text_data <- data.frame(x = c(c(1:5)), y = 3, label = pvals)
  df <- na.omit(df)
  plt <- ggplot(df, aes(y = abund, x = tbin))
  plt_bbox <- plt + 
    geom_boxplot(aes(fill = def), position = position_dodge(), 
                 outlier.size = 0.1, notch = FALSE, lwd=0.3) +
    scale_x_discrete(labels = c("1-20", "21-30", "31-40", "41-50", "51+")) +
    scale_fill_manual(values =  c("grey", "#f22a1b", "#2f868b")) +
    labs(y = expression(paste("Protein abundance ", 
                              log[10], "(iBAQ intensity)")),
         x = "Number of tissues in which protein is detected",
         fill = "Decay") +
    geom_text(data = text_data, aes(x, y, label = label), size = 2) +
    theme(text = element_text(size = 10),
          legend.key.size = unit(0.5, "cm"),
          legend.position = "none")
  return(plt_bbox)
}

figa <- protdb_binned_boxplots(df)

protdb_binned_boxplots(filter(df_protdb_human, tissue == "ovary"))
# bin_dgof_ks(filter(df_protdb_human, tissue == "ovary"))

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

# Display plots


protdb_binned_boxplots(df_protdb_human)
bin_wilcoxons(df_protdb_human)
# plot_median_boxplot(filter(df_protdb_human, tissue == "ovary"))



bin_dgof_ks <- function(df){
  df <- na.omit(df)
  df <- mutate(df, abin = cut(df$abund, c(-10, -6, -4, -2, 0, 4)))
  pvals <- c()
  for (bin in levels(df$abin)){
    ned <- filter(df, abin == bin, def == "NED")$tcount
    ed <- filter(df, abin == bin, def == "ED")$tcount
    p <- wilcox.test(ned, ed, simulate.p.value = TRUE, exact = FALSE)$p.value
    pvals <- append(pvals, as.character(signif(p, 2)))
  }
  return(pvals)
}
bin_dgof_ks(df)


protdb_binned_boxplots <- function(df){
  df <- na.omit(df)
  df$def <- factor(df$def, levels = c("UN", "ED", "NED"))
  levels(df$def) <- c("Undefined", "ED", "NED")
  df <- mutate(df, abin = cut(df$abund, c(-10, -6, -4, -2, 0, 4)))
  pvals <- bin_dgof_ks(df)
  text_data <- data.frame(x = c(c(1:5)), y = 65, label = pvals)
  df <- na.omit(df)
  plt <- ggplot(df, aes(y = tcount, x = abin))
  plt_bbox <- plt + 
    geom_boxplot(aes(fill = def), position = position_dodge(), 
                 outlier.size = 0.1, notch = FALSE, lwd=0.3) +
        scale_x_discrete(labels = c("-10 - -6", "-6 - -4", "-4 - -2", 
                                    "-2 - 0", "0 - 4")) +
    scale_fill_manual(values =  c("grey", "#f22a1b", "#2f868b")) +
    labs(x = expression(paste("log(protein abundance)")),
         y = "Number of tissues in which protein is detected",
         fill = "Decay") +
    geom_text(data = text_data, aes(x, y, label = label), size = 2) +
    theme(text = element_text(size = 10),
          legend.key.size = unit(0.5, "cm"))
  return(plt_bbox)
}
figb <- protdb_binned_boxplots(df)






