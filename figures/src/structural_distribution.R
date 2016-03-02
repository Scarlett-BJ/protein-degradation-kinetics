library("ggplot2")
library("dplyr")
library("scales")
library("gridExtra")
library("grid")
library("lattice")
library("RColorBrewer")
library("dgof")
library("boot")

humanf <- "halflife/halflife/data/structural/NED_quaternary_human.txt"
mousef <- "halflife/halflife/data/structural/NED_quaternary_mouse_ribo.txt"
df_human <- read.table(humanf, header = TRUE)
df_mouse <- read.table(mousef, header = TRUE)

df_isize <- read.table("halflife/halflife/data/structural/isize.txt",
                       header = TRUE, sep = "\t")

df_assembly <- read.table("halflife/halflife/data/structural/assembly.txt",
                          header = TRUE, sep = "\t")

################################################################################

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

# Get absolute size of each bar in stacked barplot
get_bar_sizes <- function(df){
  barvals <- c()
  for (i in c("mon", "hom", "het")){
    for (j in c("E", "U", "N")){
      bsize <- length(filter(df, qtype == i, decay.class == j)$gene)
      barvals <- append(barvals, bsize)
    }
  }
  return(matrix(barvals, nrow = 3))
}

# Plot proportions of mon, ned etc.
stacked_plotter <- function(df){
  dfinternal <- filter(df, decay.class != "X")
  dfinternal$decay.class <- factor(dfinternal$decay.class, 
                                   levels = c("E", "U", "N"))
  levels(dfinternal$decay.class) <- c("ED", "Undefined", "NED")
  dfinternal$qtype <- factor(dfinternal$qtype, 
                             levels = c("het", "hom", "mon"))
  levels(dfinternal$qtype) <- c("heteromer", "homomer", "monomer")
  barvals <- get_bar_sizes(df)
  cols <- c("#7BC342", "#6699CC", "#C0C0C0")
  # Must define fill in geom_bar (i.e. locally) for geom_text to work later
  plt <- ggplot(dfinternal, aes(x = decay.class)) +
    geom_bar(aes(fill = qtype), position = "fill", 
             colour = c("black"), lwd = 9) +
    scale_y_continuous(labels = percent) +
    scale_fill_manual(values = cols) +
    labs(x = "", y = "Subunits", fill = "") +
    theme(text = element_text(size = 10),
          legend.key.size = unit(0.3, "cm"),
          legend.position = "top",
          axis.title.x = element_blank(),
          legend.margin = unit(0, "cm")) +
#     RIBSOMES INCLUDED
    geom_text(data = data.frame(x = 1, y = c(0.56, 0.38, 0.05), 
                                lab = barvals[1,]),
              aes(x, y, label = lab), size = 2.5, alpha = 0.75) +
    geom_text(data = data.frame(x = 2, y = c(0.625, 0.44, 0.05), 
                                lab = barvals[2,]),
              aes(x, y, label = lab), size = 2.5, alpha = 0.75) +
    geom_text(data = data.frame(x = 3, y = c(0.815, 0.63, 0.05), 
                                lab = barvals[3,]),
              aes(x, y, label = lab), size = 2.5, alpha = 0.75) +
#     RIBOSOMES FILTERED
#     geom_text(data = data.frame(x = 1, y = c(0.55, 0.37, 0.05), 
#                                 lab = barvals[1,]),
#               aes(x, y, label = lab), size=2, alpha = 0.75) +
#     geom_text(data = data.frame(x = 2, y = c(0.62, 0.42, 0.05), 
#                                 lab = barvals[2,]),
#               aes(x, y, label = lab), size=2, alpha = 0.75) +
#     geom_text(data = data.frame(x = 3, y = c(0.75, 0.52, 0.05), 
#                                 lab = barvals[3,]),
#               aes(x, y, label = lab), size=2, alpha = 0.75) +
    guides(fill=guide_legend(reverse=TRUE))
  return(plt)
}
panela <-stacked_plotter(df_mouse) 
panela
## Distribution of heteromer subunit counts
density_plotter <- function(df){
  df <- filter(df, qtype == "het", decay.class != "X")
  df$decay.class <- factor(df$decay.class, levels = c("U", "E", "N"))
  levels(df$decay.class) <- c("Undefined", "ED", "NED")
  plt <- ggplot(df, aes(x = unq)) +
    geom_density(aes(fill = decay.class), alpha = 0.7, adjust = 1, lwd = 0.4) +
    scale_fill_manual(values =  c("grey", "#f22a1b", "#2f868b")) +
    scale_x_log10(lim = c(2, 512), 
                  breaks = c(2, 4, 8, 16, 32, 64, 128, 256, 512)) +
    labs(x = "Number of unique subunits", y = "Density", fill = "") +
    theme(text = element_text(size = 9), 
          legend.key.size = unit(0.5, "cm"))
  return(plt)
}
density_plotter(df_mouse)


################################################################################

## Interface size plots

# Calculate wilcoxon rank sum tests comparing NED vs ED interface size
calc_pvals <- function(df){
  df_internal <- df
  df_internal[df_internal$Unique.subunits >= 4,]$Unique.subunits <- 4
  pvals <- c()
  for (i in c(1:4)){
    neds <- filter(df_internal, Unique.subunits == i, Class == "N")
    eds <- filter(df_internal, Unique.subunits == i, Class == "E")
    p <- wilcox.test(neds$Interface.size, eds$Interface.size)$p.value
    pvals <- append(pvals, signif(p, 1))
  }
  return(data.frame(x = c(1:4), pvals = pvals))
}

# Build table summarising interface size data
build_data <- function(df){
  # Calculate bootstrap confidence intervals
  bootstrap_ci <- function(x){
    set.seed(1)
    gmean <- function(dat) exp(mean(log(dat)))
    boots.out <- boot(data = x, 
                      statistic = function(d, i) gmean(d[i]), R = 10000)
    bootresults <- boot.ci(boots.out, conf = 0.65, type = "norm")$normal
    results <- data.frame(mean = c(gmean(x)), 
                          ymin = c(bootresults[2]), 
                          ymax = c(bootresults[3]))
    return(results)
  }
  
  df_internal <- df
#   df_internal$Class <- factor(df_internal$Class, levels = c("N", "U", "E"))
#   levels(df_internal$Class) <- c("NED", "Undefined", "ED")
  df_internal$Class <- factor(df_internal$Class, levels = c("E", "U", "N"))
  levels(df_internal$Class) <- c("ED", "Undefined", "NED")
  df_internal[df_internal$Unique.subunits >= 4,]$Unique.subunits <- "4+"
  df_summary <- df_internal %>% 
    group_by(Class, Unique.subunits) %>% 
    do(bootstrap_ci(.$Interface.size))
  return(df_summary)
}

bar_plotter <- function(df){
  df_summary <- build_data(df)
  annots <- mutate(calc_pvals(df), y = c(2500, 1900, 2600, 3600))
  dodge <- position_dodge(width = 0.7)
  plt <- ggplot(df_summary, aes(x = Unique.subunits)) +
    geom_bar(aes(y = mean,  fill = Class), stat = "identity", width = 0.7,
             position = dodge, colour = "black", lwd = 0.4) +
    geom_errorbar(aes(ymin = ymin, ymax = ymax, fill = Class),
                  position = dodge, width = 0.25, lwd = 0.4) +
    # scale_fill_manual(values = c( "#2f868b", "grey", "#f22a1b")) +
    scale_fill_manual(values = c( "#f22a1b", "grey", "#2f868b")) +
    labs(x = "Unique subunits",
         y = expression(paste("Mean interface size (\uc5"^"2",")")),
         fill = "Decay") +
    theme(text = element_text(size=10),
          legend.key.size = unit(0.5, "cm"),
          axis.text=element_text(color='black')) +
    geom_text(data = annots, aes(x, y, label = pvals), 
              size = 2.5) +
    guides(fill=guide_legend(reverse=TRUE))
  return(plt)
}
bar_plotter(df_isize)

box_plotter <- function(df){
  annots <- mutate(calc_pvals(df), y = 13500)
  df_internal <- df
  df_internal$Class <- factor(df_internal$Class, levels = c("E", "U", "N"))
  levels(df_internal$Class) <- c("ED", "Undefined", "NED")
  df_internal[df_internal$Unique.subunits >= 4,]$Unique.subunits <- 4
  plt <- ggplot(df_internal, aes(x=factor(Unique.subunits))) + 
    geom_boxplot(aes(y = Interface.size, fill = Class), position = "dodge",
                 outlier.size = 0.5, lwd = 0.4) +
    scale_fill_manual(values = c( "#f22a1b", "grey", "#2f868b")) +
    scale_y_log10(breaks = c(200, 800, 3200, 12800),
                  labels = c("200", "800", 
                             "3,200", "12,800")) +
    scale_x_discrete(labels = c("1", "2", "3", "4+")) +
    labs(x = "Unique subunits",
         y = expression(paste("Interface size (\uc5"^"2",")")),
         fill = "") +
    theme(text = element_text(size=10),
          legend.key.size = unit(0.5, "cm"),
          legend.title = element_blank(),
          axis.title.y = element_text(margin = margin(r = -8)),
          legend.margin = unit(0, "cm")) +
    # guides(fill=guide_legend(reverse=TRUE)) + 
    geom_text(data = annots, aes(x, y, label = pvals), 
              size = 2.5)
  
  return(plt)
}
panelb <- box_plotter(df_isize)
panelb

################################################################################

## Assembly jitter plots

jitter_plotter <- function(df){
  df_plt <- df
  df_plt$Class <- factor(df_plt$Class, levels = c("N", "U", "E"))
  levels(df_plt$Class) <- c("NED", "Undefined", "ED")
  plt <- ggplot(df_plt, aes(x = Class)) +
    # geom_boxplot(aes(y = Normalised.assembly.order), lwd = 0.6) +
    geom_boxplot(aes(y = Normalised.assembly.order, fill = Class), 
                 alpha = 1, lwd = 0.4) +
#     geom_jitter(aes(y = Normalised.assembly.order, fill = Class), 
#                 size = 1.8, width = 0.6, colour = "black", pch=21) +
    scale_fill_manual(values = c("#2f868b", "grey", "#f22a1b")) +
    scale_colour_manual(values = c("blue4", "gray8", "red4")) +
    labs(x = "", y = "Normalised assembly order (early to late)") +
    theme(text = element_text(size=10),
          legend.key.size = unit(0.5, "cm"),
          legend.position = "none",
          axis.title.y = element_blank())
  return(plt)
}
panelc <-jitter_plotter(df_assembly) + coord_flip()
panelc

violin_plotter <- function(df){
  df_plt <- df
  df_plt$Class <- factor(df_plt$Class, levels = c("N", "U", "E"))
  levels(df_plt$Class) <- c("NED", "Undefined", "ED")
  plt <- ggplot(df_plt, aes(x = Class)) +
#     geom_violin(aes(y = Normalised.assembly.order),
#                 lwd = 0.5, alpha = 0.0) +
    geom_violin(aes(y = Normalised.assembly.order, fill = Class),
                lwd = 0.5, alpha = 0.75) +
    geom_boxplot(aes(y = Normalised.assembly.order, fill = Class),
                 lwd = 0.6, width = 0.15) +
    scale_fill_manual(values = c("#2f868b", "grey", "#f22a1b")) +
    labs(x = "", y = "Normalised assembly order") +
    theme(text = element_text(size=10),
          legend.key.size = unit(0.5, "cm"),
          legend.position = "none",
          axis.text=element_text(color='black'))
  return(plt)
}
violin_plotter(df_assembly)
