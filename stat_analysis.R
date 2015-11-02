library("dplyr")
library("ggplot2")

df <- read.table("protein_halflife/data/ortholog_abundances.txt")
colnames(df) <- c("ENSMUSP", "ENSP", "nedscore", "relabund", "orthabund")
df <- mutate(df, class = "U")
df[df$nedscore < 0.1,]$class <- "ED"
df[df$nedscore > 0.9,]$class <- "NED"

corplot <- ggplot(df, aes(relabund, orthabund, fill=class)) +
  geom_point() +
  scale_y_log10() +
  scale_x_log10() +
  ylab("Human ortholog abundance") +
  xlab("Relative mouse abundance")

corplot

ortho_box <- ggplot(df, aes(class, orthabund, fill=class)) +
  geom_boxplot() +
  scale_y_log10() +
  ylab("Human ortholog abundance")

ortho_box

rel_box <- ggplot(df, aes(class, relabund, fill=class)) +
  geom_boxplot() +
  scale_y_log10() +
  ylab("Mouse relative abundance")

rel_box

df <- read.table("protein_halflife/data/pairs.tmp", sep="\t")
colnames(df) <- c("comp", "ENSMUSP1", "ENSMUSP2", "score1", "score2", "class1",
                  "class2", "tcount1", "tcount2", "interface")
df <- mutate(df, class1 = "U", class2 = "U")
df[df$score1 > 0.9,]$class1 <- "NED"
df[df$score2 > 0.9,]$class2 <- "NED"
df[df$score2 < 0.1,]$class2 <- "ED"
df[df$score1 < 0.1,]$class1 <- "ED"
write.table(df, "protein_halflife/data/pairs.tmp", sep="\t", quote=F, row.names=F, col.names=F)

ED_ex <- length(filter(df, class == "ED", !is.na(orthabund))$class)
ED_nex <- length(filter(df, class == "ED", is.na(orthabund))$class)
NED_ex <- length(filter(df, class == "NED", !is.na(orthabund))$class)
NED_nex <-length(filter(df, class == "NED", is.na(orthabund))$class)

NEDoverED <- filter(df, class1 == "NED" & class2 != "NED" | class1 != "NED" & class2 == "NED")
length(filter(NEDoverED, class2 == "NED" & tcount2 > tcount1, interface > 200)$comp) + 
  length(filter(NEDoverED, class1 == "NED" & tcount1 > tcount2, interface > 200)$comp)

length(NEDoverED$comp) - length(filter(NEDoverED, class1 == "NED" & tcount1 > tcount2 | 
                                        class2 == "NED" & tcount2 > tcount1)$comp)
