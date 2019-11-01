
###########################################################*
# *****************     GDCRNATools     ***************** #
###########################################################*
#https://github.com/Jialab-UCR/Jialab-UCR.github.io/blob/master/GDCRNATools.workflow.R#L1
#http://bioconductor.org/packages/devel/bioc/vignettes/GDCRNATools/inst/doc/GDCRNATools.html
# GDCRNATools is an R package which provides a standard, 
# easy-to-use and comprehensive pipeline for downloading, 
# organizing, and integrative analyzing RNA expression data 
# in the GDC portal with an emphasis on deciphering the 
# lncRNA-mRNA related ceRNAs regulatory network in cancer.


# Here we provide code of the basic steps for data analysis
# by GDCRNATools. Detailed instructions can be found here:
# http://htmlpreview.github.io/?https://github.com/Jialab-UCR/Jialab-UCR.github.io/blob/master/GDCRNATools_manual.html




#=========================================================#
#          1. GDCRNATools package installation         ####
#=========================================================#

# Get the current working directory, make sure that it is 
# writable, otherwise, change to a new directory
getwd()
#setwd(workingDirectory)

# installation of GDCRNATools from Bioconductor
# source("https://bioconductor.org/biocLite.R")
# biocLite("GDCRNATools")
install.packages("remotes")
library(remotes)
remotes::install_github("Jialab-UCR/GDCRNATools")
library(GDCRNATools)






#=========================================================#
#                      2. Quick start                 ####
#=========================================================#

# A small internal dataset is used here to show the most basic 
# steps for ceRNAs network analysis in GDCRNATools




###########################################################*
##         2.1 Normalization of HTSeq-Counts data     ####


### load RNA counts data
data(rnaCounts)
rnaCounts[1:5,1:5]

### load miRNAs counts data
data(mirCounts)
mirCounts[1:5,1:5]

### Normalization of RNAseq data
rnaExpr <- gdcVoomNormalization(counts = rnaCounts, filter = FALSE)
rnaExpr[1:5,1:5]

### Normalization of miRNAs data
mirExpr <- gdcVoomNormalization(counts = mirCounts, filter = FALSE)
mirExpr[1:5,1:5]


###########################################################*





###########################################################*
##          2.2 Parse and filter RNAseq metadata       ####

metaMatrix.RNA <- gdcParseMetadata(project.id = 'TCGA-HNSC',
                                   data.type  = 'RNAseq', 
                                   write.meta = FALSE)
metaMatrix.RNA <- gdcFilterDuplicate(metaMatrix.RNA)
metaMatrix.RNA <- gdcFilterSampleType(metaMatrix.RNA)
metaMatrix.RNA[1:5,]


###########################################################*





###########################################################*
##             2.3 ceRNAs network analysis             ####


### Identification of differentially expressed genes ###

DEGAll <- gdcDEAnalysis(counts     = rnaCounts, 
                        group      = metaMatrix.RNA$sample_type, 
                        comparison = 'PrimaryTumor-SolidTissueNormal', 
                        method     = 'limma')
DEGAll[1:5,]

# All DEGs
deALL <- gdcDEReport(deg = DEGAll, gene.type = 'all')
deALL[1:5,]

# DE long-noncoding genes
deLNC <- gdcDEReport(deg = DEGAll, gene.type = 'long_non_coding')
deLNC[1:5,]

# DE protein coding genes
dePC <- gdcDEReport(deg = DEGAll, gene.type = 'protein_coding')
dePC[1:5,]



########### ceRNAs network analysis of DEGs ############
ceOutput <- gdcCEAnalysis(lnc         = rownames(deLNC), 
                          pc          = rownames(dePC), 
                          lnc.targets = 'starBase', 
                          pc.targets  = 'starBase', 
                          rna.expr    = rnaExpr, 
                          mir.expr    = mirExpr)

ceOutput[1:5,]


### Export ceRNAs network to Cytoscape
ceOutput2 <- ceOutput[ceOutput$hyperPValue<0.01 
                      & ceOutput$corPValue<0.01 & ceOutput$regSim != 0,]

# Export edges
edges <- gdcExportNetwork(ceNetwork = ceOutput2, net = 'edges')
edges[1:5,]

# Export nodes
nodes <- gdcExportNetwork(ceNetwork = ceOutput2, net = 'nodes')
nodes[1:5,]


###########################################################*






#=========================================================#
#               3. START: TCGA-HNSC               ####
#=========================================================#


###########################################################*
##                   3.1 Download data                   ####

library(GDCRNATools)
# set up directories for downloaded data
project <- 'TCGA-HNSC'
rnadir <- paste(project, 'RNAseq', sep='/')
mirdir <- paste(project, 'miRNAs', sep='/')

# install gdc-client tool
# macOS$ curl -O https://gdc.cancer.gov/system/files/authenticated%20user/0/gdc-client_v1.3.0_OSX_x64.zip
# i4$ wget https://gdc.cancer.gov/system/files/authenticated%20user/0/gdc-client_v1.4.0_Ubuntu_x64.zip --no-check-certificate 
# https://docs.gdc.cancer.gov/Data_Transfer_Tool/Users_Guide/Data_Download_and_Upload/
# 
### Download RNAseq data, n=546
gdcRNADownload(project.id     = 'TCGA-HNSC', 
               data.type      = 'RNAseq', 
               write.manifest = FALSE,
               method = 'gdc-client', ## use gdc-client tool to download data
               directory      = rnadir)

### Download miRNAs data, n=569
gdcRNADownload(project.id     = 'TCGA-HNSC', 
               data.type      = 'miRNAs', 
               write.manifest = FALSE,
               method = 'gdc-client', ## use gdc-client tool to download data
               directory      = mirdir)

###########################################################*






###########################################################*
##                 3.2 Data organization                 ####


### Parse RNAseq metadata
### there is part of clinicopathological data: tumor stage and survival
metaMatrix.RNA <- gdcParseMetadata(project.id = 'TCGA-HNSC',
                                   data.type  = 'RNAseq', 
                                   write.meta = FALSE)

# Filter duplicated samples in RNAseq metadata
metaMatrix.RNA <- gdcFilterDuplicate(metaMatrix.RNA)
# Filter non-Primary Tumor and non-Solid Tissue Normal samples in RNAseq metadata
metaMatrix.RNA <- gdcFilterSampleType(metaMatrix.RNA)



### Parse miRNAs metadata
metaMatrix.MIR <- gdcParseMetadata(project.id = 'TCGA-HNSC',
                                   data.type  = 'miRNAs', 
                                   write.meta = FALSE)

# Filter duplicated samples in miRNAs metadata
metaMatrix.MIR <- gdcFilterDuplicate(metaMatrix.MIR)
# Filter non-Primary Tumor and non-Solid Tissue Normal samples in miRNAs metadata
metaMatrix.MIR <- gdcFilterSampleType(metaMatrix.MIR)




### Merge raw counts data
# Merge RNAseq data
rnaCounts <- gdcRNAMerge(metadata  = metaMatrix.RNA, 
                         path      = rnadir, 
                         organized = FALSE, ## if target data are in folders
                         data.type = 'RNAseq')
#Number of samples: 544
#Number of genes: 60483


# Merge miRNAs data
mirCounts <- gdcRNAMerge(metadata  = metaMatrix.MIR,
                         path      = mirdir,
                         organized = FALSE, ## if target data are in folders
                         data.type = 'miRNAs')


# doc http://bioconductor.org/packages/release/bioc/manuals/GDCRNATools/man/GDCRNATools.pdf
### TMM normalization and voom transformation
# Normalization of RNAseq data (rnaCounts)
rnaExpr <- gdcVoomNormalization(counts = rnaCounts, filter = FALSE)

# Normalization of miRNAs data (mirCounts)
mirExpr <- gdcVoomNormalization(counts = mirCounts, filter = FALSE)


### Differential gene expression analysis
# R4> colnames(DEGAll)
#[1] "symbol"  "group"   "logFC"   "AveExpr" "t"       "PValue" 
#[7] "FDR"     "B" 
# https://www.statisticshowto.datasciencecentral.com/false-discovery-rate/
# #The FDR approach is used as an alternative to the Bonferroni correction and controls for a low proportion of false positives
DEGAll <- gdcDEAnalysis(counts     = rnaCounts, 
                        group      = metaMatrix.RNA$sample_type, 
                        comparison = 'PrimaryTumor-SolidTissueNormal', 
                        method     = 'limma')
#data(DEGAll)

# All DEGs
deALL <- gdcDEReport(deg = DEGAll, gene.type = 'all')
# n=2194

# DE long-noncoding
deLNC <- gdcDEReport(deg = DEGAll, gene.type = 'long_non_coding')
# n=108

# DE protein coding genes
dePC <- gdcDEReport(deg = DEGAll, gene.type = 'protein_coding')
# n=2043

# table(deALL$group)

# IG long_non_coding           ncRNA  protein_coding 
# 4             108               2            2043 
# pseudogene             TEC              TR 
# 34               3               0 
###########################################################*










###########################################################*
##    3.3 Competing endogenous RNAs network analysis    ####


### The 3 steps of ceRNAs network analysis:
# Hypergeometric test
# Pearson correlation analysis
# Regulation pattern analysis

### All of the 3 steps can be performed in a single function


### ceRNAs network analysis using internal databases
ceOutput <- gdcCEAnalysis(lnc         = rownames(deLNC), 
                          pc          = rownames(dePC), 
                          lnc.targets = 'starBase', 
                          pc.targets  = 'starBase', 
                          rna.expr    = rnaExpr, 
                          mir.expr    = mirExpr)
# Error in cor.test.default(pcDa, lncDa, alternative = "greater") : 
# not enough finite observations


### ceRNAs network analysis using user-provided datasets
# load miRNA-lncRNA interactions
data(lncTarget)
lncTarget[1:3]

# load miRNA-mRNA interactions
data(pcTarget)
pcTarget[1:3]

ceOutput <- gdcCEAnalysis(lnc         = rownames(deLNC), 
                          pc          = rownames(dePC), 
                          lnc.targets = lncTarget, 
                          pc.targets  = pcTarget, 
                          rna.expr    = rnaExpr, 
                          mir.expr    = mirExpr)


### Network visulization in Cytoscape

# Filter potential ceRNA interactions
ceOutput2 <- ceOutput[ceOutput$hyperPValue<0.01 & 
                        ceOutput$corPValue<0.01 & ceOutput$regSim != 0,]


# Edges and nodes can be simply imported into Cytoscape 
# for network visualization
edges <- gdcExportNetwork(ceNetwork = ceOutput2, net = 'edges')
edges[1:5,]

nodes <- gdcExportNetwork(ceNetwork = ceOutput2, net = 'nodes')
nodes[1:5,]

write.table(edges, file='edges.txt', sep='\t', quote=F) ### Network of Cytoscape
write.table(nodes, file='nodes.txt', sep='\t', quote=F) ### Table of Cytoscape


### Correlation plot on a local webpage

shinyCorPlot(gene1    = rownames(deLNC), 
             gene2    = rownames(dePC), 
             rna.expr = rnaExpr, 
             metadata = metaMatrix.RNA)


###########################################################*

# an short example
# https://rdrr.io/github/Jialab-UCR/GDCRNATools/man/gdcSurvivalAnalysis.html
# In Jialab-UCR/GDCRNATools: GDCRNATools: an R/Bioconductor package for integrative analysis of lncRNA, mRNA, and miRNA data in GDC
# https://rdrr.io/github/Jialab-UCR/GDCRNATools/man/
#source("https://bioconductor.org/biocLite.R")
#biocLite("GDCRNATools") 
chooseCRANmirror()
install.packages("BiocManager")
library("BiocManager")
BiocManager::install("GDCRNATools")
library(GDCRNATools) # citation("pathview") for GSEA


genes <- c('ENSG00000000938','ENSG00000000971','ENSG00000001036',
           'ENSG00000001084','ENSG00000001167','ENSG00000001460')

samples <- c('TCGA-2F-A9KO-01', 'TCGA-2F-A9KP-01',
             'TCGA-2F-A9KQ-01', 'TCGA-2F-A9KR-01',
             'TCGA-2F-A9KT-01', 'TCGA-2F-A9KW-01')

metaMatrix <- data.frame(sample_type=rep('PrimaryTumor',6),
                         sample=samples,
                         days_to_death=seq(100,600,100),
                         days_to_last_follow_up=rep(NA,6))
rnaExpr <- matrix(c(2.7,7.0,4.9,6.9,4.6,2.5,
                    0.5,2.5,5.7,6.5,4.9,3.8,
                    2.1,2.9,5.9,5.7,4.5,3.5,
                    2.7,5.9,4.5,5.8,5.2,3.0,
                    2.5,2.2,5.3,4.4,4.4,2.9,
                    2.4,3.8,6.2,3.8,3.8,4.2),6,6)
rownames(rnaExpr) <- genes
colnames(rnaExpr) <- samples
survOutput <- gdcSurvivalAnalysis(gene=genes,
                                  rna.expr=rnaExpr, metadata=metaMatrix)
# hazard ratio, 95% confidence interval, P-value, and FDR
# # => where is FDR? 
# gdcDEAnalysis -> DEGAll has FDR.
# to screen out the differentially expressed genes (DEGs) with a fold change >2, and P value was defined as .05 to be statistically significant. 
# Volcano plot was drafted in RStudio and genes whose fold-change >2 along with false discovery rate (FDR) <0.1 were marked with red (upregulated) and green (downregulated). 
# 



###########################################################*
##             3.4 Other downstream analyses          ####


############### Univariate survival analysis          ####
# # how about: grouping by expression high/low with a optimized cutoff point?
# CoxPH analysis
survOutput_cox <- gdcSurvivalAnalysis(gene     = rownames(deALL), 
                                  method   = 'coxph', 
                                  rna.expr = rnaExpr, 
                                  metadata = metaMatrix.RNA)


# KM analysis
# # hazard ratio, 95% confidence interval, P-value.
survOutput_km <- gdcSurvivalAnalysis(gene     = rownames(deALL), 
                                  method   = 'KM', 
                                  rna.expr = rnaExpr, 
                                  metadata = metaMatrix.RNA, 
                                  sep      = 'median')
# P-value < 0.05
alpha_HNSCC <- 0.05; zcut <- 1.5
attach(survOutput_km)
# removal of as.factor
#survOutput_km$pValue <- formatC(as.numeric(as.character(pValue)), format = "e", digits = 2)
survOutput_km$pValue <- signif(as.numeric(as.character(pValue)), digits = 3)
survOutput_km$HR <- signif(as.numeric(as.character(HR)), digits=3)
survOutput_km$lower95 <- signif(as.numeric(as.character(lower95)), digits=3)
survOutput_km$upper95 <- signif(as.numeric(as.character(upper95)), digits=3)
survOutput_km <- survOutput_km[order(pValue, -HR), ] #sorting by order(ascending)
# FDR calculation: BH, Benjamini and Hochberg; http://www.omicshare.com/forum/thread-173-1-1.html
# p.adjust.methods # c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")
pValue_adj <- p.adjust(survOutput_km$pValue, method="fdr", n = nrow(survOutput_km))
survOutput_km <- cbind(survOutput_km, pValue_adj)
detach(survOutput_km)
#  & HR>=zcut, FDR <0.05
survOutput_km005 <- survOutput_km[which(survOutput_km$pValue_adj<=alpha_HNSCC), 1:6]
#x n=317; badguy 153, goodguy 5
# FDR <0.05 adjustment n=只有7個
survOutput_km005_bad <- survOutput_km[which(survOutput_km$pValue<=alpha_HNSCC & survOutput_km$HR>=zcut), 1:5]
survOutput_km005_good <- survOutput_km[which(survOutput_km$pValue<=alpha_HNSCC & survOutput_km$HR<(0.8)), 1:5]
# badguy 3(HR>2.2), goodguy 4(HR<0.77) => FDR 很棒

# KM plot on a local webpage by shinyKMPlot, a dynamic plot
# shiny => Listening on http://127.0.0.1:3118

shinyKMPlot(gene = rownames(deALL), rna.expr = rnaExpr, 
            metadata = metaMatrix.RNA)




############## Functional enrichment analysis          ####

### All the functional enrichment analyses can be 
### performed in a single function, including:
# Gene Ontology (BP, CC, MF) analysis
# KEGG pathway analysis
# Disease Ontology analysis

enrichOutput <- gdcEnrichAnalysis(gene = rownames(deALL), simplify = TRUE)

#data(enrichOutput)

# Barplot
gdcEnrichPlot(enrichOutput, type = 'bar', category = 'GO', num.terms = 10)

# Bubble plot
gdcEnrichPlot(enrichOutput, type='bubble', category='GO', num.terms = 10)


# View pathway maps on a local webpage
library(pathview)

deg <- deALL$logFC
names(deg) <- rownames(deALL)
pathways <- as.character(enrichOutput$Terms[enrichOutput$Category=='KEGG'])

shinyPathview(deg, pathways = pathways, directory = 'pathview')


###########################################################*

