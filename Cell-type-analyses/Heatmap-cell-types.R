## Script to plot gene expression data across different cell types as a heatmap
## Use for a short list of genes identified to be of interest (cell type-specific, responsive to cue, gene family)

# k.kajala at uu dot nl
# github: kaisakajala

# last update 2018.08.13 - KK



library(ggplot2)
library(reshape2)

setwd("~/Rwork/1808-Atlas-paper/Heatmaps/TF-heatmap/")

#ggplot method
data <- read.csv("TFs.csv")
head(data)
#data <- data[,1:15] ##weird NA rows appear, remove them
data.m <- melt(data)
data.m


pdf("EXO-MCO-specific-TF-Heatmap.pdf")
ggplot(data.m, aes(variable, X)) + geom_tile(aes(fill = value),
     colour = "white") + scale_fill_gradient(low = "white", high = "magenta") +
     theme(text = element_text(size=8))
dev.off()



data <- read.csv("TFs-noMCO.csv")
data.m <- melt(data)

pdf("EXO-specific-TF-Heatmap.pdf")
ggplot(data.m, aes(variable, Symbol)) + geom_tile(aes(fill = value),
  colour = "white") + scale_fill_gradient(low = "white", high = "magenta") +
  theme(text = element_text(size=8))
dev.off()





#ggplot method
data <- read.csv("TFs-by-venn-group.csv")
head(data)
#data <- data[,1:15] ##weird NA rows appear, remove them
data.m <- melt(data)
data.m


pdf("EXO-MCO-specific-TF-Heatmap-by-venn-group.pdf")
ggplot(data.m, aes(variable, Symbol)) + geom_tile(aes(fill = value),
        colour = "white") + scale_fill_gradient(low = "white", high = "magenta") +
  theme(text = element_text(size=8))
dev.off()