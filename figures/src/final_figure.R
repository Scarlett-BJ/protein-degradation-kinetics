library("ggplot2")
library("dplyr")
library("scales")
library("gridExtra")
library("dgof")

# Panel A data (with ribosomes)
df_mouse <- read.table("halflife/halflife/data/revised_data/NED_quaternary_mouse.txt", 
                       header = TRUE)

# Panel B data
df_isize <- read.table("halflife/halflife/data/structural/isize.txt",
                       header = TRUE, sep = "\t")

# Panel C data
df_assembly <- read.table("halflife/halflife/data/structural/assembly.txt",
                          header = TRUE, sep = "\t")
df_assembly$Class <- factor(df_assembly$Class, levels = c("N", "U", "E"))
levels(df_assembly$Class) <- c("NED", "Undefined", "ED")
# Panel D data
df_struc <- read.table("halflife/halflife/data/coexpression/struc_coex.txt",
                       header = TRUE, sep = "\t")
colnames(df_struc) <- c("gene", "decay.class", "cid", "avg.coex", "seqid", "usubs")
df_struc$decay.class <- factor(df_struc$decay.class, 
                               levels = c("U", "E", "N"))
levels(df_struc$decay.class) <- c("Undefined", "ED", "NED")


## Panel A - Stacked barplot showing proportion of ED in diff quat. strucs.

# Returns fisher.test pval for matrix containing chosen classes
get_fisher <- function(dfx, dclass){
  ned_hets <- length(filter(dfx, qtype == "het", 
                            decay.class == "N")$gene)
  ned_non_hets <- length(filter(dfx, qtype != "het", 
                                decay.class == "N")$gene)
  other_hets <- length(filter(dfx, qtype == "het", 
                              decay.class == dclass)$gene)
  other_non_hets <- length(filter(dfx, qtype != "het", 
                                  decay.class == dclass)$gene)
  fmatrix <- matrix(c(ned_hets, ned_non_hets, other_hets, other_non_hets),
                    nrow = 2, ncol = 2)
  f <- fisher.test(fmatrix)
  return(signif(f$p.value, 1))
}

# Get absolute size of each bar in stacked barplot
get_bar_sizes <- function(dfx){
  barvals <- c()
  for (i in c("mon", "hom", "het")){
    for (j in c("E", "U", "N")){
      bsize <- length(filter(dfx, qtype == i, decay.class == j)$gene)
      barvals <- append(barvals, bsize)
    }
  }
  return(matrix(barvals, nrow = 3))
}

# Return panel A plot
stacked_plotter <- function(df_plt){
  barvals <- get_bar_sizes(df_plt)
  pvals <- c(get_fisher(df_plt, "E"), get_fisher(df_plt, "U"))
  dfinternal <- filter(df_plt, decay.class != "X")
  dfinternal$decay.class <- factor(dfinternal$decay.class, 
                                   levels = c("E", "U", "N"))
  levels(dfinternal$decay.class) <- c("ED", "Undefined", "NED")
  dfinternal$qtype <- factor(dfinternal$qtype, 
                             levels = c("het", "hom", "mon"))
  levels(dfinternal$qtype) <- c("heteromer", "homomer", "monomer")
  cols <- c("#7BC342", "#6699CC", "#C0C0C0")
  # Must define fill in geom_bar (i.e. locally) for geom_text to work later
  plt <- ggplot(dfinternal, aes(x = decay.class)) +
    geom_bar(aes(fill = qtype), position = "fill", 
             colour = c("black"), lwd = 0.4) +
    scale_y_continuous(labels = percent,
                       breaks = c(0.0, 0.25, 0.50, 0.75, 1.00)) +
    scale_fill_manual(values = cols) +
    labs(x = "", y = "Subunits", fill = "") +
    theme(text = element_text(size = 9),
          legend.key.size = unit(0.3, "cm"),
          legend.position = "top",
          axis.title.x = element_blank(),
          legend.margin = unit(0, "cm")) +
    # annotations is for ribosomes-sincluded
    geom_text(data = data.frame(x = 1, y = c(0.72, 0.46, 0.05), 
                                lab = barvals[1,]),
              aes(x, y, label = lab), size = 2.5, alpha = 0.6) +
    geom_text(data = data.frame(x = 2, y = c(0.79, 0.53, 0.05), 
                                lab = barvals[2,]),
              aes(x, y, label = lab), size = 2.5, alpha = 0.6) +
    geom_text(data = data.frame(x = 3, y = c(0.90, 0.68, 0.05), 
                                lab = barvals[3,]),
              aes(x, y, label = lab), size = 2.5, alpha = 0.6) +
    geom_text(data = data.frame(x = c(2, 2.5), y = c(1.135, 1.075), 
                                lab = pvals),
              aes(x, y, label = lab), size = 2.5, fontface = "italic") +
    guides(fill=guide_legend(reverse=TRUE))
  return(plt)
}

panela <-stacked_plotter(df_mouse) 
panela

## Panel B - Dodged boxplots showing increase in interface size with NEDs
# Cacluate wilcoxon p-values
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

# Return panel B plot
box_plotter <- function(df){
  pval_annots <- mutate(calc_pvals(df), y = 15000)
  df_internal <- df
  df_internal$Class <- factor(df_internal$Class, levels = c("E", "U", "N"))
  levels(df_internal$Class) <- c("ED", "Undefined", "NED")
  df_internal[df_internal$Unique.subunits >= 4,]$Unique.subunits <- 4
  rcounts <- df_internal %>% group_by(Unique.subunits) %>% summarise(counts = n())
  plt <- ggplot(df_internal, aes(x=factor(Unique.subunits))) + 
    geom_boxplot(aes(y = Interface.size, fill = Class), position = "dodge",
                 outlier.size = 0.5, lwd = 0.4) +
    scale_fill_manual(values = c( "#f22a1b", "grey", "#2f868b")) +
    scale_y_log10(limits = c(140, 18000),
                  breaks = c(200, 800, 3200, 12800),
                  labels = c("200", "800", 
                             "3,200", "12,800")) +
    scale_x_discrete(labels = c("1", "2", "3", "4+")) +
    labs(x = "Unique subunits",
         y = expression(paste("Interface size (\uc5"^"2",")")),
         fill = "") +
    theme(text = element_text(size=9),
          legend.key.size = unit(0.3, "cm"),
          legend.title = element_blank(),
          axis.title.y = element_text(margin = margin(r = -3)),
          legend.margin = unit(0, "cm")) +
    # guides(fill=guide_legend(reverse=TRUE)) + 
    geom_text(data = pval_annots, aes(x, y, label = pvals), 
              size = 2.5, fontface = "italic") +
    geom_text(data = rcounts, aes(x = c(1, 2, 3, 4), y = 150, label = counts), 
              size = 2.5, alpha = 0.6)
  return(plt)
}

panelb <- box_plotter(df_isize)
panelb

## Panel C - boxplots showing NED proteins assemble earlier
get_raw_counts <- function(df){
  ned <- length(filter(df, Class == "NED")$Class)
  un <- length(filter(df, Class == "Undefined")$Class)
  ed <- length(filter(df, Class == "ED")$Class)
  return(c(ed, un, ned))
}

get_wilcoxon_pvals <- function(df){
  ned_ed <- wilcox.test(filter(df, Class == "NED")$Normalised.assembly.order,
                        filter(df, Class == "ED")$Normalised.assembly.order)
  ned_un <- wilcox.test(filter(df, Class == "NED")$Normalised.assembly.order,
                        filter(df, Class == "Undefined")$Normalised.assembly.order)
  ed_un <- wilcox.test(filter(df, Class == "ED")$Normalised.assembly.order,
                       filter(df, Class == "Undefined")$Normalised.assembly.order)
  return(c(signif(ned_un$p.value, 1), 
           signif(ned_ed$p.value, 1),
           signif(ed_un$p.value, 1)))
}

jitter_plotter <- function(df){
  df_plt <- df
  rcounts <- get_raw_counts(df_plt)
  pvals <- get_wilcoxon_pvals(df_plt)
  plt <- ggplot(df_plt, aes(x = Class)) +
    geom_boxplot(aes(y = Normalised.assembly.order, fill = Class), 
                 alpha = 1, lwd = 0.4) +
    #     geom_jitter(aes(y = Normalised.assembly.order, fill = Class), 
    #                 size = 1.8, width = 0.6, colour = "black", pch=21) +
    scale_fill_manual(values = c("#2f868b", "grey", "#f22a1b")) +
    scale_colour_manual(values = c("blue4", "gray8", "red4")) +
    scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1),
                   limits = c(0, 1.18),
                   labels = c("0.00", "0.25", "0.50", "0.75", "1.00")) +
    labs(x = "", y = "Normalised assembly order, early to late") +
    theme(text = element_text(size = 9),
          legend.key.size = unit(0.3, "cm"),
          legend.position = "none",
          axis.title.y = element_blank()) +
    geom_text(data = data.frame(x = c(1, 2, 3), y = c(0.25, 0.3, 0.58), 
                                lab = rcounts),
              aes(x, y, label = lab), size = 2.5, alpha = 0.6) +
    geom_text(data = data.frame(x = c(1.5, 2, 2.5), y = c(1.05, 1.15, 1.05), 
                                lab = pvals),
              aes(x, y, label = lab), size = 2.5, fontface = "italic")
  return(plt)
}

panelc <-jitter_plotter(df_assembly) + coord_flip()
panelc

## Panel D - NED subunits are more highly coepressed
density_plotter <- function(df, bintype = "struc"){
  if (bintype == "struc"){ 
    df$bin <- cut(df$usubs, c(2, 4, 6,  10, 100))
    levels(df$bin) <- c("3-4", "5-6", "7-10", "11+")
  } else if (bintype == "corum"){
    df$bin <- cut(df$usubs, c(2, 5, 10, 15, 20, 25, 30, 40, 80, 200))
    levels(df$bin) <- c("3-5", "6-10", "11-15", "16-20", "21-25", "26-30",
                        "30-40", "40-80", "80+")
  }
  df <- na.omit(df)
  plt <- ggplot(df, aes(x = avg.coex)) +
    geom_density(aes(fill = decay.class), alpha = 0.7, lwd = 0.4) +
    scale_fill_manual(values = c("grey", "#f22a1b", "#2f868b")) +
    labs(x = "Average subunit coexpression (Pearson's r)",
         y = "Density",
         fill = "") +
    scale_x_continuous(limits = c(-0.25, 0.9),
                       breaks = c(-0.2, 0, 0.2, 0.4, 0.6, 0.8),
                       labels = c("-0.2", "0.0", "0.2", "0.4", "0.6", "0.8")) +
    theme(text = element_text(size = 9),
          legend.key.size = unit(0.3, "cm"),
          legend.title = element_blank(),
          legend.margin = unit(0, "cm")) +
    guides(fill=guide_legend(reverse=TRUE))
  return(plt)
}

paneld <- density_plotter(filter(df_struc, usubs >= 3)) + 
  geom_text(data = data.frame(x = 0.55, y = 2.5, lab = "P = 1.2e-08"), 
            aes(x, y, label = lab), size = 2.5, fontface = "italic")

grid.arrange(panela, panelb, panelc, paneld, 
             widths = c(5, 7), heights = c(1, 1), nrow = 2)

################################################################################

## Supplementary figures

# Panel A with ribosmomes filtered out
mouse_no_ribo_file <- "halflife/halflife/data/structural/NED_quaternary_mouse.txt"
df_mouse_no_ribosomes <- read.table(mouse_no_ribo_file, header = TRUE)

spanel_1 <- stacked_plotter(df_mouse_no_ribosomes)
spanel_1

# Structural data, coexpression density filtered by number of unique subunits
spanel_2a <- density_plotter(filter(df_struc, usubs >= 3)) +
  facet_wrap(~ bin, ncol=2)
spanel_2a

mouse_file <- "halflife/halflife/data/coexpression/coexpressdb_corum_mouse_homologs.tsv"
header <- c("cid", "usubs", "uniprot.id", "avg.coex", "decay.class", "species")
df_corum_mouse <- read.table(mouse_file, header = TRUE, sep = "\t")
colnames(df_corum_mouse) <- header
df_corum_mouse$decay.class <- factor(df_corum_mouse$decay.class, 
                                     levels = c("UN", "ED", "NED"))
levels(df_corum_mouse$decay.class) <- c("Undefined", "ED", "NED")

spanel_2bi <-  density_plotter(filter(df_corum_mouse, usubs >= 3), 
                              bintype = "corum") + facet_wrap(~ bin, ncol=3)
spanel_2bi

spanel_2c <-  density_plotter(filter(df_corum_mouse, usubs >= 3), 
                              bintype = "corum") +
  geom_text(data = data.frame(x = 0.55, y = 2.5, lab = "P = 3.8e-55"), 
            aes(x, y, label = lab), size = 2.5, fontface = "italic")
spanel_2c



dgof_ks <- function(df, c1, c2){
  c1v <- factor(filter(df, qtype == "het", decay.class == c1)$unq)
  c2v <- factor(filter(df, qtype == "het", decay.class == c2)$unq)
  print(length(c1v))
  print(length(c2v))
  result = ks.test(c1v, ecdf(c2v), simulate.p.value = TRUE, exact = FALSE)
  print(result)
  return(result$p.value)
}

heteromer_density_plotter <- function(df){
  df <- filter(df, qtype == "het", decay.class != "X")
  df$decay.class <- factor(df$decay.class, levels = c("U", "E", "N"))
  levels(df$decay.class) <- c("Undefined", "ED", "NED")
  plt <- ggplot(df, aes(x = unq)) +
    geom_density(aes(fill = decay.class), alpha = 0.7, adjust = 1, lwd = 0.4) +
    scale_fill_manual(values =  c("grey", "#f22a1b", "#2f868b")) +
    scale_x_log10(lim = c(2, 512), 
                  breaks = c(2, 4, 8, 16, 32, 64, 128, 256, 512)) +
    labs(x = "Number of unique subunits", y = "Density", fill = "") +
    theme(text = element_text(size=9),
          legend.key.size = unit(0.3, "cm"),
          legend.title = element_blank(),
          axis.title.y = element_text(margin = margin(r = -3)),
          legend.margin = unit(0, "cm"))
  return(plt)
}


spanel_3 <- grid.arrange(heteromer_density_plotter(df_mouse) + 
                           geom_text(data = data.frame(x = 8, y = 2, 
                                                       lab = "P < 2.2e-16"), 
                                     aes(x, y, label = lab), size = 2.5, 
                                     fontface = "italic"),
                         heteromer_density_plotter(df_mouse_no_ribosomes) +
                           geom_text(data = data.frame(x = 8, y = 2, 
                                                       lab = "P < 2.2e-16"), 
                                     aes(x, y, label = lab), size = 2.5, 
                                     fontface = "italic"))

