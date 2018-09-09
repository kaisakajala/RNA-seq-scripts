## RNA seq analysis using limma-voom

# k.kajala at uu dot nl
# github: kaisakajala

# scripts updated on limma-voom script by:
# jrodriguezm at ucdavis dot edu
# github: rodriguezmDNA
# brady lab

# last update 2018.08.13 - KK

######################## User defined options ########################
######################################################################
# Main directory:
setwd("~/Rwork/1808-Atlas-paper/limma-voom/") #Full path of the working directory. 
#It must contain: 
# ** A directory named "Counts" with a delimited file for raw counts. Rows are genes and columns samples.
# ** A directory named "meta" with a metadata file with information about the samples. Ideally the number of rows in the metadata is the same as in the raw counts.
# ** A directory named Scripts with this script and the 'functions.R' script.


## Metadata options
metaFile <- "20180811_atlas-SL-meta_final2.csv" #Name of metadata file in csv format.
doFilter <- T #F
whichFilter <- c("PH.T","XY.TQ","WOX.Ta","XY.I","WOX.Tb","MCO.S")

## Counts file name (with extension) 
countsFile <- "20180811_ATLAS_ITAG3.2_Kallisto_raw_counts.txt" #Name of counts file


shortName <- "180813-Atlas" #Short name to append at the end of the filenames. If missing it will append the name of the folder where the scripts where run.

## Filter genes with low expression using CPM counts.
filterByCPM <- T #T
CPMcutoff <- 0.5 #2


## pValue (default = 0.05 ) and absolute logFC (default = 2) to color genes on volcano plots
pValCut=0.05  #0.05  
logCut=2 #2

########################
########################

###

library(edgeR)
library(reshape)
library(gplots)
library(RColorBrewer)
library(calibrate)
library(Glimma) 

## Output
outDir = "Atlas-DiffExprs-ALL-withoutOutliers/"
dir.create(outDir, showWarnings=T)

sink('Atlas-DiffExprs-ALL-withoutOutliers/Atlas-drops-analysis-output.txt')

geneListsDir = "Atlas-DiffExprs-ALL-withoutOutliers/GeneLists"
dir.create(geneListsDir, showWarnings=T)
#
imgDir = "Atlas-DiffExprs-ALL-withoutOutliers/images/"
dir.create(imgDir, showWarnings=T)
## --

if (is.na(shortName)){
  shortName <- basename(getwd())
}

# Load functions
source("Scripts/functions.R")
######## --- --- --- 


## Start of analysis
####################################################################################
####################################################################################
cat("Reading metadata file \n")

meta <- metaDataProcessing(metaFile,doFilter,whichFilter)
head(meta)

#
cat("Reading counts file:",countsFile,"\n")

GeneCounts <- read.delim(paste0("Counts/",countsFile),row.names = 1)
dim(GeneCounts)

## Check that samples in both counts and metadata are the same.
## Use function filterCounts(counts,meta)
tmp <- filterCounts(GeneCounts,meta)
GeneCounts <- tmp[["counts"]]
meta <- tmp[["meta"]]
rm(tmp)
## --


###### Design matrix
## Convert experimental metadata to factors for the design
experimentFactors <- lapply(apply(meta,2,split,""),unlist)
experimentFactors <- as.data.frame(lapply(experimentFactors,as.factor))

cat ("Create the design with these factors:\n")
print(head(experimentFactors))

###  User modified:
####Simplest design taking into account all possible interactions
Groups <- as.factor(paste0(experimentFactors$Promoter))
design <- model.matrix(~0+Groups) 
# Example of an interaction
#design <- model.matrix(~0+experimentFactors$Sample*experimentFactors$Treatment) #Sample*Treatment interaction

## Ensures column names are optimal for the contrast design
fixCols <- paste(c("Groups","experimentFactors","\\$","\\:","\\-",
                  colnames(experimentFactors)),sep="",collapse = "|")

colnames(design) <- gsub(fixCols,"",colnames(design))
head(design)


####################################################################################
cat("Removing genes with 0 counts on all conditions \n")
cat("Initial number of genes:",nrow(GeneCounts),"\n")
rmIDX <- which(rowSums(GeneCounts) == 0)
cat("Removing",length(rmIDX),"genes \n")
GeneCounts <- GeneCounts[-rmIDX,]
cat("Remaining number of genes:",nrow(GeneCounts),"\n")

### Use cpms to uncover lowly expressed genes
dge <- DGEList(counts=GeneCounts,remove.zeros = F)


# Filter genes with low CPMs accross replicates 
cat("Replicates of samples range between:", range(table(Groups)),"\n")

#
if (filterByCPM){
  
  sampleMin <- min(table(Groups))
  cat("Filtering reads with low CPMs ( <",CPMcutoff,") in at least",sampleMin,"replicates \n")
  #
  cpm <- cpm(dge)
  keep.exprs <- rowSums(cpm>CPMcutoff)>=sampleMin
  table(keep.exprs)
  

  cat("Removing",table(keep.exprs)[1],"genes \n")
  cat("Remaining number of genes:",table(keep.exprs)[2],"\n")
  
  #
  
  y <- dge[keep.exprs, , keep.lib.size = FALSE]
} else {
  cat("Not doing CPM filtering")
}

## Normalization with TMM
#y <- calcNormFactors(y, method = "TMM")
# A paper from Gordon K Smyth, who developed limma uses TMM before voom:
# http://dx.doi.org/10.12688/f1000research.9005.1
#####


normalizedExpression <- cpm(y)
#Save a file with cpm normalized expression
tmpSave <- paste(outDir,"normalizedExpression","_",shortName,".csv",sep="")
#cat("Saving normalized data to:", tmpSave, "\n")
#save(normalizedExpression,file = "cpm_normalizedExpression.RData")
write.csv(x=normalizedExpression,tmpSave,quote = F,row.names = T)


## Easier visualization for MDS plots

cat("Using glimma for MDS plot visualization - Normalized data \n")
glMDSPlot(y, labels=rownames(y$samples),
          groups=meta, launch=T)



#### Start PDF
tmpSave <- paste(imgDir,"DEG_Analysis_",shortName,".pdf",sep="")
pdf(tmpSave,paper = "USr")


### Use voom on the dge object.
v <- voom(y, design, plot = TRUE,normalize.method ="quantile")
# OR:
#v <- voomWithQualityWeights(y,design, normalization="quantile",plot = T)
##Where sample-level variation is evident from earlier inspections of the MDS plot, the voomWithQualityWeights
#function can be used to simultaneously incorporate sample-level weights together with the abundance dependent
#weights estimated by voom
###
cat("Analyzing",nrow(v),"with",ncol(v),"libraries \n")

######## Visualization and quality control

##################
## Correlation between replicates of samples belonging to same group
corrSamples <- cor(v$E)
## --
tmpSave <- paste(imgDir,"CorrelationBetweenReplicates_",shortName,".pdf",sep="")
pdf(tmpSave,paper = "USr")#width = 8,height = 6)
#colors <- colorRampPalette(c("darkgoldenrod4","darkgoldenrod1","white","white","steelblue1","steelblue4"))

for (each in (levels(Groups))){
  hmData <- corrSamples[grep(each, rownames(corrSamples)),grep(each, colnames(corrSamples))]
  #hmData <- corrSamples[,grep(each, colnames(corrSamples))]
  hm <- T
  if(!hm){
    cat("Heatmaps with NMF \n")
    NMF::aheatmap(hmData,col=colors(125),
                  txt = ifelse(hmData<0.85,"<",NA),#Rowv = F,Colv = F,
                  main=paste0("Correlation between samples of group ",each))
    
  } else {
    cat("Heatmaps with heatmap.2 \n")
    heatmap.2(hmData,col=cm.colors(125), keysize = 0.75,
              cellnote = ifelse(hmData<0.8,"*",NA), notecol = "black",
              
              #margins = c(16,16),
              
              dendrogram = "none", trace = "none",density.info='none',
              cexCol  = 0.8 ,cexRow = 0.8,
              lmat=rbind(c(4, 3, 9),
                         c(2, 1, 6),
                         c(8, 5, 7)),
              lhei=c(0.3, 0.6,0.8),
              lwid=c(0.25, 0.4,0.2),
              main=paste0("Correlation\n",each))
    legend("bottomleft",legend = "* means correlation < 0.85",bty = "n")
  }
  
}
dev.off()
################## --


## Assign colors to each of the experimental factors. 
ColorTable <- assignColorsToFactors(experimentFactors)

## Boxplot of normalized counts ordered by Groups
boxplot(v$E[,order(Groups)], range=0,col=customColors[Groups[order(Groups)]], 
        ylab="log2[counts]", xlab="sample", main="Quantile normalized Counts",
        cex.axis=0.5,las=2)



### Do contrasts
Groups
cont.matrix= makeContrasts(
    "EXO-EP"=EXO-EP,
    "EXO-COR"=EXO-COR,
    "EXO-MCO"=EXO-MCO,
    "EXO-EN"=EXO-EN,
    "EXO-V"=EXO-V,
    "EXO-XY"=EXO-XY,
    "EXO-PH"=EXO-PH,
    "EXO-WOX"=EXO-WOX,
    "EXO-MZ"=EXO-MZ,
    "EXO-ACT"=EXO-ACT,
    "EXO-35S"=EXO-X35S,
    
    "EP-EXO"=EP-EXO,  
    "EP-COR"=EP-COR,
    "EP-MCO"=EP-MCO,
    "EP-EN"=EP-EN,
    "EP-V"=EP-V,
    "EP-XY"=EP-XY,
    "EP-PH"=EP-PH,
    "EP-WOX"=EP-WOX,
    "EP-MZ"=EP-MZ,
    "EP-ACT"=EP-ACT,
    "EP-35S"=EP-X35S,
    
    "COR-EP"=COR-EP,
    "COR-EXO"=COR-EXO,  
    "COR-MCO"=COR-MCO,
    "COR-EN"=COR-EN,
    "COR-V"=COR-V,
    "COR-XY"=COR-XY,
    "COR-PH"=COR-PH,
    "COR-WOX"=COR-WOX,
    "COR-MZ"=COR-MZ,
    "COR-ACT"=COR-ACT,
    "COR-35S"=COR-X35S,
    
    "MCO-EP"=MCO-EP,
    "MCO-EXO"=MCO-EXO,  
    "MCO-COR"=MCO-COR,
    "MCO-EN"=MCO-EN,
    "MCO-V"=MCO-V,
    "MCO-XY"=MCO-XY,
    "MCO-PH"=MCO-PH,
    "MCO-WOX"=MCO-WOX,
    "MCO-MZ"=MCO-MZ,
    "MCO-ACT"=MCO-ACT,
    "MCO-35S"=MCO-X35S,
    
    "EN-EP"=EN-EP,
    "EN-EXO"=EN-EXO,  
    "EN-COR"=EN-COR,
    "EN-MCO"=EN-MCO,
    "EN-V"=EN-V,
    "EN-XY"=EN-XY,
    "EN-PH"=EN-PH,
    "EN-WOX"=EN-WOX,
    "EN-MZ"=EN-MZ,
    "EN-ACT"=EN-ACT,
    "EN-35S"=EN-X35S,
    
    "V-EP"=V-EP,
    "V-EXO"=V-EXO,  
    "V-COR"=V-COR,
    "V-MCO"=V-MCO,
    "V-EN"=V-EN,
    "V-XY"=V-XY,
    "V-PH"=V-PH,
    "V-WOX"=V-WOX,
    "V-MZ"=V-MZ,
    "V-ACT"=V-ACT,
    "V-35S"=V-X35S,
    
    "XY-EP"=XY-EP,
    "XY-EXO"=XY-EXO,  
    "XY-COR"=XY-COR,
    "XY-MCO"=XY-MCO,
    "XY-EN"=XY-EN,
    "XY-V"=XY-V,
    "XY-PH"=XY-PH,
    "XY-WOX"=XY-WOX,
    "XY-MZ"=XY-MZ,
    "XY-ACT"=XY-ACT,
    "XY-35S"=XY-X35S,

    "PH-EP"=PH-EP,
    "PH-EXO"=PH-EXO,  
    "PH-COR"=PH-COR,
    "PH-MCO"=PH-MCO,
    "PH-EN"=PH-EN,
    "PH-V"=PH-V,
    "PH-XY"=PH-XY,
    "PH-WOX"=PH-WOX,
    "PH-MZ"=PH-MZ,
    "PH-ACT"=PH-ACT,
    "PH-35S"=PH-X35S,
    
    "WOX-EP"=WOX-EP,
    "WOX-EXO"=WOX-EXO,  
    "WOX-COR"=WOX-COR,
    "WOX-MCO"=WOX-MCO,
    "WOX-EN"=WOX-EN,
    "WOX-V"=WOX-V,
    "WOX-XY"=WOX-XY,
    "WOX-PH"=WOX-PH,
    "WOX-MZ"=WOX-MZ,
    "WOX-ACT"=WOX-ACT,
    "WOX-35S"=WOX-X35S,
    
    "MZ-EP"=MZ-EP,
    "MZ-EXO"=MZ-EXO,  
    "MZ-COR"=MZ-COR,
    "MZ-MCO"=MZ-MCO,
    "MZ-EN"=MZ-EN,
    "MZ-V"=MZ-V,
    "MZ-XY"=MZ-XY,
    "MZ-PH"=MZ-PH,
    "MZ-WOX"=MZ-WOX,
    "MZ-ACT"=MZ-ACT,
    "MZ-35S"=MZ-X35S,
    
    "ACT-EP"=ACT-EP,
    "ACT-EXO"=ACT-EXO,  
    "ACT-COR"=ACT-COR,
    "ACT-MCO"=ACT-MCO,
    "ACT-EN"=ACT-EN,
    "ACT-V"=ACT-V,
    "ACT-XY"=ACT-XY,
    "ACT-PH"=ACT-PH,
    "ACT-WOX"=ACT-WOX,
    "ACT-MZ"=ACT-MZ,
    "ACT-35S"=ACT-X35S,
    
    "X35S-EP"=X35S-EP,
    "X35S-EXO"=X35S-EXO,  
    "X35S-COR"=X35S-COR,
    "X35S-MCO"=X35S-MCO,
    "X35S-EN"=X35S-EN,
    "X35S-V"=X35S-V,
    "X35S-XY"=X35S-XY,
    "X35S-PH"=X35S-PH,
    "X35S-WOX"=X35S-WOX,
    "X35S-MZ"=X35S-MZ,
    "X35S-ACT"=X35S-ACT,
    levels=design)


#### Fit and do Diff Expression
fit <- lmFit(v, design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)

## -- Summary and Venn diagrams , only good for up to 5 comparisons.
results <- decideTests(fit2)
summary(results)
if (ncol(results) <= 5){
  cat ("Doing Venn Diagrams \n")
  vennDiagram(results,include = c("up","down"), main="DE")
} else {
  cat ("More than 5 comparisons, skipping Venn Diagrams  \n")
}
DESummary <- t(summary(decideTests(fit2)))[,-2]
colnames(DESummary) = c("Downregulated","Upregulated")

# Save as csv
tmpSave <- paste(outDir,"DESummaryInteractions_",shortName,".csv",sep="")
write.csv(x=DESummary,tmpSave,quote = F,row.names = T)

# Write to PDF
plotData <- t(DESummary)
yMax <- max(colSums(plotData))
rownames(plotData) <- c("Down","Up")

barplot(plotData,legend.text = rownames(plotData),col=c("orange","steelblue4"),
        xlab = "Contrast", ylab = "Number of genes",
        beside = T,
        ylim = c(0,yMax*1.2),
        las=2,
        cex.names = 0.6, border = T, bty="n",
        main="DE genes per contrast")


## Prepare for gene annotation
annotationAvail <- F
if (annotationAvail){
  cat("Reading annotation file \n")
  genealiasfile <- "gene_aliases.txt"
  ID2Symbol <- getGeneSymbols(genealiasfile)
} else{cat("Annotation file unavailable \n")}


## --

DEList <- list()
for (contrast in colnames(cont.matrix)){
  print(contrast)
  ## Sorting by none ensures all contrasts will be in the same order
  tmp <- topTable(fit2, coef=contrast,number = Inf,sort.by = "none")
  #
  pValpassed <- table(tmp$adj.P.Val < 0.05)[2]
  cat ("Number of genes with pVal < 0.05 on ",contrast,":",pValpassed,"\n")
  
  
  ## Write genes that are up or downregulated (logFC > 0; logFC < 0)
  upGenes <- as.data.frame(rownames(tmp[tmp$adj.P.Val < 0.05 & tmp$logFC > 0,]))
  tmpSave <- paste(geneListsDir,"/",contrast,"_up",".csv",sep="")
  write.csv(x=upGenes,tmpSave,quote = F,row.names = T)
  #
  downGenes <- as.data.frame(rownames(tmp[tmp$adj.P.Val < 0.05 & tmp$logFC < 0,]))
  tmpSave <- paste(geneListsDir,"/",contrast,"_down",".csv",sep="")
  write.csv(x=downGenes,tmpSave,quote = F,row.names = T)
  #####
  
  #-- Add gene symbols if available
  tmp[,"Symbol"] <- rownames(tmp)
  if (annotationAvail){
    cat("Adding annotation \n")
    Genes <- rownames(tmp)
    idx <- intersect(names(AGI2Symbol),Genes)
    tmp[idx,"Symbol"] <- AGI2Symbol[idx]
    Genes
  }
  #--
  
  ## Add contrast name to the column names, in case of multiple contrasts.
  colnames(tmp) <- paste(colnames(tmp),contrast,sep = ".")
  
  # Write each contrast to file
  tmpSave <- paste(outDir,contrast,"_",shortName,".csv",sep="")
  write.csv(x=tmp,tmpSave,quote = F,row.names = T)

  # Save result to list
  DEList[[contrast]] <- tmp 
}

tmpSave <- paste(outDir,"DEList_",shortName,".RData",sep="")
save(DEList,file = tmpSave)
### ------------


tmpSave <- paste(imgDir,"VolcanoPlots_",shortName,".pdf",sep="")
pdf(tmpSave,paper = "USr")
makeVolcanoPlots(DEList,pValCut=0.01,logCut=2,plotGenes=F) #plotGenes=T to print genes in the plot
dev.off()

### Condense into a single list
DE_All <- condenseListTables(DEList) ## Use a custom function
DE_All <- DE_All[,-grep("t.|B.|P.Value|AveExpr|Symbol",colnames(DE_All))] #Remove unwanted columns

meanExp <- meanNormalizedExpression(normalizedExpression,levels(Groups)) 
#By default uses the unique names on 'Groups' to calculate mean between replicates of the 
# 'normalizedExpression' data frame.

DE_All <- cbind(DE_All,meanExp[rownames(DE_All),])

tmpSave <- paste(outDir,"DEG_AllContrasts_",shortName,".csv",sep="")
write.csv(x = DE_All,file = tmpSave,quote = F,row.names = T)


## Close main img
dev.off()
sink()

