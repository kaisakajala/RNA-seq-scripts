## Script to create a Boxplot after RNASeq to visualize the changes in gene expression and variation
## across different sample types (treatment, tissue, etc).

# k.kajala at uu dot nl
# github: kaisakajala

# scripts based on script by:
# jrodriguezm at ucdavis dot edu
# github: rodriguezmDNA
# brady lab

# last update 2017.10.17 - KK


setwd("~/Rwork/1707-SUB-paper/CoreHyp-boxplots-SL/")

#read in normalized data from limma
load("cpm_normalizedExpression_SL_SUB.RData")
head(normalizedExpression)

#read in a list of genes for visualization
CoreHypoxia <- read.csv("SolycCoreHypoxia2.csv")
head(CoreHypoxia)
symbols <- as.character(CoreHypoxia$SolycID)

######### User variables.
testSample <- c("TOT","TRA")
controlSample <- "INRNA"
controlTissue <- "HAIRY" 
testTreatment <- "SUB"
controlTreatment <- "CON"
#########

## Pull up meta data and subset so that only lines of interest are being looked at
meta1 <- read.csv("SUB-SL-metadata-drops-orgsub.csv")
meta1 <- subset(meta, meta$Tissue == "HAIRY")
meta1 <- subset(meta1, meta1$Genotype == "SL")
meta1$Group <- paste0(meta1$Tissue,meta1$Treatment,meta1$Sample)
Names <- as.character(meta1$Name)
Names
meta1

# Subset the expression data to have only these SL ROOT/SHOOT TOT tissues
plotdata1 <- subset(normalizedExpression, select=Names)
head(plotdata1)
str(plotdata1)

#Subset the expression data to only have the wanted genes
plotdata2 <- t(plotdata1)
plotdata2 <- subset(plotdata2, select=symbols)
plotdata2 <- t(plotdata2)
geneNames <- read.csv("SolycCoreHypoxia3.csv")
rownames(plotdata2) <- geneNames$Symbol


## Transform table.  
temp1 <- data.frame(meta1[,c("Name","Tissue","Genotype","Sample","Treatment","Replicate","Group")],t(plotdata2))
row.names(temp1) <- Names
str(temp1)
row.names(temp1)
head(colnames(temp1))

resultsList <- list()
### Operate over every sample AND the control
for (x in testTissue){
  each <- gsub("\\/","",x)
  each <- paste(each,controlSample,sep= "|")
  print(each)
  
  ## Get the factor levels
  #
  Treatment <- temp1$Treatment
  Treatment <- relevel(factor(Treatment),ref=controlTreatment)
  #
  Sample <- temp1$Sample
  Sample <- relevel(factor(Sample),ref=controlSample)
  #
  Group <- temp1$Group
  Group <- relevel(factor(Group),ref=paste0(controlTissue,controlTreatment,controlSample))
  
  #########
  # Save the boxplots as pdf
  
  titulo = ("SL-HAIRY-core-hypoxia-geneID.pdf")
  pdf(titulo,paper = "a4r",width = 12,height = 15)
  par(mfrow=c(3,4))
  
  lmResults <- list()
  #for (i in seq(6,10)){
  for (i in seq(8,ncol(temp1))){  
    
    # Get the gene name. Since the table is transposed, Gene Names are now in the columns
    genNam <- colnames(temp1)[i]
    print(genNam)
    
    #options("scipen"=-100, "digits"=4)
    
    #########
    ####
    #summary(groupModel)
    boxplot(as.numeric(temp1[,i]) ~ Treatment*Sample ,main=genNam,las=2,cex.axis=0.8)  
  }
  dev.off()
}



