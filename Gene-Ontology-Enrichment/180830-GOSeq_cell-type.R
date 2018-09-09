## Gene Ontology Enrichment analysis using GO-seq

# k.kajala at uu dot nl
# github: kaisakajala

# scripts built upon GoSeq scripts set up by:
# jrodriguezm at ucdavis dot edu
# github: rodriguezmDNA
# brady lab

# last update 2018.08.30 - KK
## https://pathwaycommons.github.io/guide/primers/statistics/fishers_exact_test/


setwd("~/Rwork/1808-Atlas-paper/Cell type lists/")

## Load libraries and Functions
  library("goseq")
  library(rtracklayer)
  library(GenomicRanges)
  library(Rsamtools)
  library(purrrlyr) 

source("~/Rwork/JoelsRcode/GOSeq_jrmScript/omicFunctions/metaFunctions_forNetworkAnalysis.R")
source("~/Rwork/JoelsRcode/GOSeq_jrmScript/omicFunctions/generic_functions.R")
source("~/Rwork/JoelsRcode/GOSeq_jrmScript/omicFunctions/orthologueFindingFunctions.R")
source("~/Rwork/JoelsRcode/GOSeq_jrmScript/omicFunctions/GOseq_functions.R")

######

## Load genes to test for enrichment
## Data input from exodermis-specific list

  de.genes <- read.csv("~/Rwork/1808-Atlas-paper/Cell type lists/EXO-filtered.csv", header=T, row.names = 1)
  de.genes <- removeDotGene(rownames(de.genes))
  de.genes
  
  
  
  ## load all DE genes
  load ("~/Rwork/1808-Atlas-paper/limma-voom/Atlas-DiffExprs-ALL/DEList_180813-Atlas.RData")
  
  assayed.genes <- removeDotGene(rownames(DEList$`EXO-EP`))  
  head(assayed.genes,100)
  
  gene.vector=as.integer(assayed.genes%in%de.genes)
  names(gene.vector)=removeDotGene(assayed.genes)
  head(gene.vector)
  
  

  ######### Prepare GFF data
  ####################################
  GFFfile = "~/Rwork/JoelsRcode/GOSeq_jrmScript/meta/ITAG3.2_gene_models.gff"
  GFF <- import.gff(GFFfile,version="3",feature.type="gene")
  grl <- IRanges::reduce(split(GFF, mcols(GFF)$Name))
  reducedGTF <- unlist(grl, use.names=T)
  mcols(reducedGTF)$Name <- rep(names(grl), elementNROWS(grl))
  reducedGTF
  AllLengths <- width(reducedGTF)
  names(AllLengths) <- mcols(reducedGTF)$Name
  names(AllLengths) <- removeDotGene(names(AllLengths))
  head(AllLengths)
  
  ## Get the gene lengths only of the genes assayed 
  # (the AllLengths contains information for ALL the genes in the genome but not all of them are *normally* used)
  GeneLengths <- AllLengths[names(gene.vector)]
  
  
#########


######### Prepare list of GO terms
####################################

{ # If reading the 3.2 version
  go <- read.delim("~/Rwork/JoelsRcode/GOSeq_jrmScript/meta/ITAG3.2_protein_go_Split_jrm.tsv",header = F,stringsAsFactors = F)
  colnames(go) <- c("ITAG",
                    "GOID")
  
  go.goseq <- go[,c("ITAG", "GOID")]
  go.goseq$ITAG <- removeDotGene (go.goseq$ITAG)
  head(go.goseq)
}
#########

######### GO terms
####################################################################################################
  

  ########## GOseq functions
  ##################################################
  ## nullp creates a weight probability of a gene based DE based on it's length
  pdf("GO-results/EXO-Probability-plot.pdf")
  pwf=nullp(gene.vector, bias.data=GeneLengths)
  dev.off()
  
  head(pwf)

  
  go.all <- goseq(pwf, "ITAG3.2", gene2cat=go.goseq)
  ################################################################################
  
  ### Filter significant categories
  table(go.all$over_represented_pvalue < 0.05)
  go.sign <- go.all[go.all$over_represented_pvalue < 0.05,]
  
  # Save to a list
  write.csv(go.sign,"GO-results/EXO-enrichment.csv")


