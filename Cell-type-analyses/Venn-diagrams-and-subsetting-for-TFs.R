## Script to create Venn diagrams across different DE lists
## First subsetting for transcription factors of interest.

# k.kajala at uu dot nl
# github: kaisakajala

# scripts updated on script by Mauricio Reynoso
# Julia Bailey-Serres lab

# last update 2018.08.10 - KK



library(systemPipeR)

removeDotGene <- function(dotGeneName,ColumnName="",verbose=F){
  if (verbose) cat ("removeDotGene function called\n")
  tmp <- gsub("\\.[0-9].*$","",unlist(dotGeneName))
  return(tmp)
}



# Set working directory
setwd("~/Rwork/1808-Atlas-paper/Exo-specific lists/Venn and Subsetting-for-TFs")


### Load gene lists
TFs <- read.table("inputs/TF_combined-180813.txt", header=F)
TFs <- as.vector(TFs[,1])

MCOroku <- read.csv("inputs/MCO-specific-by-roku.csv", header=T, row.names=1)
MCOroku <- removeDotGene(as.vector(row.names(MCOroku)))

EXOroku <- read.csv("inputs/EXO-specific-by-roku.csv", header=T, row.names=1)
EXOroku <- removeDotGene(as.vector(row.names(EXOroku)))

EXObrady <- read.csv("inputs/DEG_AllContrasts_180813-Atlas-FILTERED.csv", header=T,row.names=1)
EXObrady <- removeDotGene(as.vector(row.names(EXObrady)))

EXOrokuOld <- read.csv("inputs/180730-exo-specific-roku-quant-norm-means.csv", header=T, row.names=1)
EXOrokuOld <- removeDotGene(as.vector(row.names(EXOrokuOld)))

EXObradyOld <- read.csv("inputs/DEG_AllContrasts_180730-Atlas-filtered-annotated.csv", header=T,row.names=1)
EXObradyOld <- removeDotGene(as.vector(row.names(EXObradyOld)))

EXOstar.roku <- read.csv("inputs/STAR-eXpress-EXO-sp-genes.csv",header=T,row.names=1)
EXOstar.roku <- removeDotGene(as.vector(row.names(EXOstar.roku)))


MCOrokuTF <- MCOroku[MCOroku %in% TFs]
EXOrokuTF <- EXOroku[EXOroku %in% TFs]
EXObradyTF <- EXObrady[EXObrady %in% TFs]
EXOrokuOldTF <- EXOrokuOld[EXOrokuOld %in% TFs]
EXObradyOldTF <- EXObradyOld[EXObradyOld %in% TFs]
EXOstar.rokuTF <- EXOstar.roku[EXOstar.roku %in% TFs]



####### List of vectors
setlist<-list(EXOroku=EXOroku, EXObrady=EXObrady, EXOrokuOld=EXOrokuOld, EXObradyOld=EXObradyOld, EXOstar.roku=EXOstar.roku)
setlistTF<-list(EXO.roku.TF=EXOrokuTF, EXO.Brady.TF=EXObradyTF, EXO.roku.old.TF=EXOrokuOldTF, EXO.Brady.old.TF=EXObradyOldTF, EXO.star.roku.TF=EXOstar.rokuTF)

setlist2<-list(EXOroku=EXOroku,MCOroku=MCOroku)
setlistTF2<-list(EXO.roku.TF=EXOrokuTF,MCO.roku.TF=MCOrokuTF)



pdf("results/Exodermis-specific-lists.pdf")
## set a venn object using the list ## See the order can be modified here
vennset_up <- overLapper(setlist[c(1,2,3,4,5)],type="vennsets")
##  Plot set up and finalization
vennPlot(list(vennset_up), mymain="Exodermis-specific genes by different methods", colmode=2, ccol=c("red"))
dev.off()

#output data in csv
n=list
n=vennlist(vennset_up)
dl_all <- data.frame(ID = rep(names(n), sapply(n, length)),
                     Obs = unlist(n))
head(dl_all)
write.csv(dl_all,"results/Exolist_Venn_data.csv")

 
pdf("results/Exodermis-specific-TFs.pdf")
## set a venn object using the list ## See the order can be modified here
vennset_up <- overLapper(setlistTF[c(1,2,3,4,5)],type="vennsets")
##  Plot set up and finalization
vennPlot(list(vennset_up), mymain="Exodermis-specific TFs by different methods", colmode=2, ccol=c("red"))
dev.off()

#output data in csv
n=list
n=vennlist(vennset_up)
dl_all <- data.frame(ID = rep(names(n), sapply(n, length)),
                     Obs = unlist(n))
head(dl_all)
write.csv(dl_all,"results/Exo_TF_list_Venn_data.csv")




pdf("results/MCO-vs-EXO-roku-lists.pdf")
## set a venn object using the list ## See the order can be modified here
vennset_up <- overLapper(setlist2[c(1,2)],type="vennsets")
##  Plot set up and finalization
vennPlot(list(vennset_up), mymain="EXO vs MCO-specific genes by roku", colmode=2, ccol=c("red"))
dev.off()

#output data in csv
n=list
n=vennlist(vennset_up)
dl_all <- data.frame(ID = rep(names(n), sapply(n, length)),
                     Obs = unlist(n))
head(dl_all)
write.csv(dl_all,"results/EXO-MCO_Venn_data.csv")


pdf("results/MCO-vs-EXO-roku-TF-lists.pdf")
## set a venn object using the list ## See the order can be modified here
vennset_up <- overLapper(setlistTF2[c(1,2)],type="vennsets")
##  Plot set up and finalization
vennPlot(list(vennset_up), mymain="EXO vs MCO-specific TFs by roku", colmode=2, ccol=c("red"))
dev.off()

#output data in csv
n=list
n=vennlist(vennset_up)
dl_all <- data.frame(ID = rep(names(n), sapply(n, length)),
                     Obs = unlist(n))
head(dl_all)
write.csv(dl_all,"results/EXO-MCO_TF_Venn_data.csv")