# install JRE on BASH
# >tar zxvf jre-9_osx-x64_bin.tar.gz
# 
# 
# https://cran.r-project.org/doc/manuals/r-patched/R-admin.html
# sudo R CMD javareconf
# 
# install packages in batch on BASH
# > Rscript -e 'install.packages("rms", repos="http://cran.csie.ntu.edu.tw/")' # or "https://ftp.yzu.edu.tw/"
# packages: PairedData, psych, survival, reshape, data.table, R.utils, compositions, graphics, ggplot2, rms, xlsx, r2excel, rJava

# https://cran.r-project.org/web/packages/ggplot2/
# packages: stringi, stringr, reshape2
download.file("https://cran.r-project.org/src/contrib/ggplot2_2.2.1.tar.gz", destfile="~/Downloads/ggplot2_2.2.1.tar.gz")
# # >install.packages("~/Downloads/reshape2_1.4.3.tar.gz", repos = NULL, type = "source")
install.packages("~/Downloads/ggplot2_2.2.1.tar.gz", repos = NULL, type = "source")
library("ggplot2")



# 
# # brew cask install java
# ?sudo ln -f -s $(/usr/libexec/java_home)/jre/lib/server/libjvm.dylib /usr/local/lib
# # JDK 9 for Mac => http://www.oracle.com/technetwork/java/javase/downloads/jdk9-downloads-3848520.html
# # JAVA_HOME        : /Library/Java/JavaVirtualMachines/jdk-9.0.1.jdk/Contents/Home
# Java library path: $(JAVA_HOME)/lib/server
# 
# 
# # DO NOT RUN this {https://www.r-bloggers.com/upgrade-r-without-losing-your-packages/
# tmp <- installed.packages()
# installedpkgs <- as.vector(tmp[is.na(tmp[,"Priority"]), 1])
# save(installedpkgs, file="installed_old.rda")
# # 
# # curl -#ROL https://cran.rstudio.com/bin/macosx/R-3.4.3.pkg
# # open R-3.4.3.pkg
# #
# # >brew install r --with-java
# # >brew upgrade r --with-java # To upgrade to 3.4.3
# # automatic update all installed packages
# #  updateR ?
# tmp <- installed.packages()
# installedpkgs.new <- as.vector(tmp[is.na(tmp[,"Priority"]), 1])
# missing <- setdiff(installedpkgs, installedpkgs.new)
# install.packages(missing)
# update.packages()
# 
# # any packages from BioConductor
# # 
# chooseBioCmirror()
# biocLite() 
# load("installed_old.rda")
# tmp <- installed.packages()
# installedpkgs.new <- as.vector(tmp[is.na(tmp[,"Priority"]), 1])
# missing <- setdiff(installedpkgs, installedpkgs.new)
# for (i in 1:length(missing)) biocLite(missing[i])
# 
# # } 






# 
# data.frame prepare (from HNSCC TMA core results, Qupath) ####


# 
# 
# # path_ta51 <- "/Users/Apple/Documents/Tex_proposal_MHlab/TissueMicroarray_HNSCC_TA51/"
# path_ta51 <- "/Users/apple/Documents/Tex_proposal_MHlab/TissueMicroarray_HNSCC_TA51/TCGA_HNSCC_surgical_margin/"
# setwd(path_ta51) # .R in iCloud; new ID of Macbook Air is apple not Apple
# TA51B_map <- read.csv("/Users/Apple/Documents/Tex_proposal_MHlab/TissueMicroarray_HNSCC_TA51/TA51_OSCC_clinicopathological_dataset_PMM1_IHC/TMA51C_map-Table 1.csv", header=F)
# #TA51B_core <- read.csv("/Users/Apple/Documents/Tex_proposal_MHlab/TissueMicroarray_HNSCC_TA51/TA51_OSCC_clinicopathological_dataset_PMM1_IHC/clinical_clean-Table 1.csv", header=T)
# 
# TA1 <- unlist(TA51B_map[1:10, 2:17], byrow=T)
# TA2 <- unlist(TA51B_core$Unique.ID)
# View(TA51B_map[1:10, 2:17])
# TA2
# #transform (melting) the data.frame
# TA <- unlist(TA51B)
# TA51B <- matrix(TA, nrow = dim(TA51B)[1]/4, ncol = 4, byrow=T)
# #
#cgdsr R package (http://www.cbioportal.org/cgds_r.jsp)
#
## install.packages("TCGA2STAT") # since 2015
library("TCGA2STAT")
LUAD <- getTCGA(disease="LUAD", data.type="RNASeq2", type="RSEM", clinical=TRUE, cvars=c("yearstobirth","gender","vitalstatus","daystodeath","daystolastfollowup","daystolastknownalive","pathologicstage","pathologyTstage","pathologyNstage","pathologyMstage","residualtumor"))
# http://www.liuzlab.org/TCGA2STAT/ClinicalVariables.pdf with data schema: only 11 features



# http://www.learn-r-the-easy-way.tw/chapters/16 for dplyr and more and apply()
# dependencies ‘knitr’, ‘rmarkdown’ 
install.packages("/Users/apple/Downloads/knitr_1.18.tar.gz", repos = NULL, type="source")
install.packages("/Users/apple/Downloads/rmarkdown_1.8.tar.gz", repos = NULL, type="source")
install.packages("/Users/apple/Downloads/reprex_0.1.1.tar.gz", repos = NULL, type="source")
install.packages("tidyverse") # with reprex package; for all 3 packages:
library(tidyverse)
# ── Attaching packages ─────────────────────── tidyverse 1.2.1 ──
# ✔ ggplot2 2.2.1     ✔ purrr   0.2.4
# ✔ tibble  1.4.1     ✔ dplyr   0.7.4
# ✔ tidyr   0.7.2     ✔ stringr 1.2.0
# ✔ readr   1.1.1     ✔ forcats 0.2.0
# http://marcoghislanzoni.com/blog/2014/09/01/pivot-tables-r-dplyr/
# 1) magrittr	Stefan Milton Bache	能夠使用 %>% operators, pipeline
# beyond_start <- 1969
# beyond_yr <- Sys.Date() %>%
#   +     format(., format = "%Y") %>% # format(Sys.Date(), format="%Y"); "." is the insert point
#   +     as.numeric() %>%
#   +     `-` (beyond_start) # `` is tilt
# beyond_yr

# 2) tidyr	Hadley Wickham	能夠進行長寬表格的轉換 => pivot table, by gather(), spread()
# in native r: melt, cast, and aggregate.
# or package: rpivotTable
# 
# 3) dplyr	Hadley Wickham	更有效率地作資料處理 data.frame manipulation in SQL manner
# and pivot table as well
library(dplyr) # drop = FALSE => to become a vector result
# filter()	篩選符合條件的觀測值 age > 65
# select()	選擇變數 newID = name_ID, new_age = age
# mutate()	新增變數 with calculation
# arrange()	依照變數排序觀測值, desc(age) or age
# summarise()	summarize, mean(age)
# > group_by(straw_hat_df, gender) %>%
# +     summarise(mean(age)) %>%
#   +     as.data.frame() # convert back from tibble.data.frame (or table.data.frame) to R data.frame
# 
# data_tdf <- tbl_df(data) # convet to table.data.frame



# [xLoad from TCGA: x TCGA2STAT] ##
# [TCGA2STAT] download from TCGA cBioportal, prepared a LIST is ready for survival analysis
# citation("CNTools") # from bioconductor
source("https://bioconductor.org/biocLite.R")
biocLite("CNTools") # to process the segmented CNV or CNA data into gene-level data.


# https://github.com/IARCbioinfo/awesome-TCGA # list of all usefull TCGA tools
# or https://gdc.cancer.gov/access-data/gdc-community-tools, such as GDCtools
# FirebrowseR - Paper describing the R FirebrowseR package.
# GenomicDataCommons - Paper describing the R GenomicDataCommons package.
## [1.Directly] start from beginning  get Broad Institute GDAC: TCGA/Firhose data into R (Retrieve TCGA CDEs verbatim) ####
# FirebrowseR - An R package to download directly the results of the analyses performed by Firehose in R.
#go to FireBrowse ( http://gdac.broadinstitute.org/ ):
install.packages("devtools")
library("devtools")
devtools::install_github("mariodeng/FirebrowseR") # with more features (81): such as residual_tumor, vital_status, days_to_last_followup, "smoking duration"

library(FirebrowseR)
#
#manuscript and TCGA survival analysis Interpretation ##
#https://www.biostars.org/p/153013/ Tutorial: Survival analysis of TCGA patients
#integrating gene expression (RNASeq) data Why would you want to do survival
#analysis based on gene expression data? Well, let's say you have a number of
#genes that you are interested in and they are differentially expressed between
#tumor and normal samples, it would be very powerful to show that alteration in
#gene expression correlates with worse survival or earlier tumor recurrence.
#https://groups.google.com/forum/m/#!msg/ucsc-cancer-genomics-browser/YvKnWZSsw1Q/3IAkkEMyFa4J
#from Mary Goldman and Jing Zhu, UCSC Cancer Browser
#https://genome-cancer.ucsc.edu/ 


# set path on google drive
path_LUAD <- "/Users/apple/Google\ Drive/2018_PhD_dissertation/LUAD_Peter_survival/"
# path_LUAD <- "/Users/apple/Documents/My\ Tableau\ Repository/Workbooks/"
setwd(path_LUAD) # change working directory to the google drive

TCGA_cohort <- "LUAD" # first define the cancer type of cohort
LUAD.clinical.Fire <- Samples.Clinical(cohort = "LUAD", format="csv") # csv or tsv? n=150? yes page one
# [LUAD.....Fire is our source .Rda file]
# 
#colnames(LUAD.clinical.Fire)
#TCGA_05_5420 <- LUAD.clinical.Fire[LUAD.clinical.Fire$tcga_participant_barcode == "TCGA-05-5420",]

#{ clinical features CDEs by FirebrowseR; n=522, 81 features (CDEs)
# obtaining all CDEs
all.Found <- F
page.Counter <- 1
LUAD.clinical.Fire <- list()
page.Size = 600 # using a bigger page size is faster
while(all.Found == F){
  LUAD.clinical.Fire[[page.Counter]] = Samples.Clinical(format = "csv",
                                                       cohort = "LUAD",
                                                       page_size = page.Size,
                                                       page = page.Counter)
  if(nrow(LUAD.clinical.Fire[[page.Counter]]) < page.Size)
    all.Found = T
  else
    page.Counter = page.Counter + 1
}
LUAD.clinical.Fire <- do.call(rbind, LUAD.clinical.Fire)
#}
save(LUAD.clinical.Fire, file="LUAD.clinical.Fire.Rda")







# skip
# x inner join: merging 2D tables # NOT by "patient_id"
# LUAD.clinical.Fire.Rda and LUAD.mRNA.Exp.Fire.Rda
# # set path on google drive
path_LUAD <- "/Users/apple/Google\ Drive/2018_PhD_dissertation/LUAD_Peter_survival/"
# path_LUAD <- "/Users/apple/Documents/My\ Tableau\ Repository/Workbooks/"
setwd(path_LUAD) # change working directory to the google drive

load(file="LUAD.clinical.Fire.Rda")
load(file="LUAD.mRNA.Exp.Fire.Rda")
LUAD.clinico_mRNA.Fire <- merge(LUAD.clinical.Fire, LUAD.mRNA.Exp.Fire, by="tcga_participant_barcode") #n=576 with sample "TP" or "TR", excluding "NT normal tissue"
# without duplicated ID
save(LUAD.clinico_mRNA.Fire, file=paste("LUAD.clinico_mRNA.", diff.Exp.Genes, ".Fire.Rda", sep=""))

# its name will be "LUAD.clinico_mRNA.ZZZ3.Fire.Rda"


save(LUAD, file="LUAD_TCGA2STAT.Rda")
load(file="LUAD_TCGA2STAT.Rda")
# 20501 genes have been imported (RNASeqV2 data)! Clinical data (n=515) will be imported.
# The returned object rnaseq.ov will be a list containing three elements:
#   
#   dat - The omics-profiles data matrix with genes in rows and patients in columns; in this example, this is a matrix of RPKM values.
# clinical - Clinical data matrix; in this example, NULL is returned as the clinical data is not specified (default).
# merged.dat - If clinical data is imported, a data matrix with both omics-profiles and clinical covariates (overall survival by default); patients are in the rows with the sample ID in the first column, followed by clinical covariates, and then the omics-profiles.
View(LUAD)
dim(LUAD$merged.dat) # n=515, genes=20501, and clinical features=10
# clinico.LUAD <- LUAD[["clinical"]] # LUAD$merged.dat[,1:21]
head(LUAD$merged.dat[1:3, c(1:3, 4:5, 6,7:9, 10, 20510:20511)])
# c(11:20509) for RNASeq2 by RSEM; gene profile start from col 11
# colnames(LUAD[["merged.dat"]])[1:21]
# [1] "bcr"             "YEARSTOBIRTH"    "GENDER"          "VITALSTATUS"    
# [5] "DAYSTODEATH"     "PATHOLOGICSTAGE" "PATHOLOGYTSTAGE" "PATHOLOGYNSTAGE"
# [9] "PATHOLOGYMSTAGE" "RESIDUALTUMOR"   "A1BG"            "A1CF"           
# ....
# [20508] "ZZEF1"     "ZZZ3"      "psiTPTE22" "tAKR"

head(LUAD$merged.dat[1:3, c(1:3, 4:5, 6,7:9, 10, 20510:20511)])
# rnaseq.LUAD <- getTCGA(disease="LUAD", data.type="RNASeq", type="RPKM", clinical=TRUE, 
#                       cvars=c("yearstobirth","gender","vitalstatus","daystodeath","pathologicstage","pathologyTstage","pathologyNstage","pathologyMstage"))
# list of clinical features: http://www.liuzlab.org/TCGA2STAT/ClinicalVariables.pdf
# colnames(clinico.LUAD[,1:21]) # full features
# I want to pick up c(11,2:4,7:10,19)
# [1] "Composite Element REF"           
# [2] "yearstobirth"                    
# [3] "vitalstatus"                     
# [4] "daystodeath"                     
# [5] "daystolastfollowup"              
# [6] "tumortissuesite"                 
# [7] "pathologicstage"                 
# [8] "pathologyTstage"                 
# [9] "pathologyNstage"                 
# [10] "pathologyMstage"                 
# [11] "gender"                          
# [12] "dateofinitialpathologicdiagnosis"
# [13] "daystolastknownalive"            
# [14] "radiationtherapy"                
# [15] "karnofskyperformancescore"       
# [16] "histologicaltype"                
# [17] "numberpackyearssmoked"           
# [18] "yearoftobaccosmokingonset"       
# [19] "residualtumor"                   
# [20] "race"                            
# [21] "ethnicity" 

# col1: bcr = barcode "TCGA-05-4244"
clinico.LUAD <- LUAD$merged.dat[, c(1:3, 4:5, 6,7:9, 10, 20510:20511)] #n=515; no RNASeq2
save(clinico.LUAD, file="clinico_LUAD.Rda") # clinical features only
#
RNASeq2.LUAD <- LUAD$merged.dat[, c(1,11:20509)] # [, 11:] for RNASeq2 by RSEM (20499 genes) only
save(RNASeq2.LUAD, file="RNASeq2_LUAD.Rda")

load(file="RNASeq2_LUAD.Rda") # dim 515 x 20500
# > colnames(RNASeq2.LUAD)[1:10]
# [1] "bcr"    "A1BG"   "A1CF"   "A2BP1"  "A2LD1"  "A2ML1"  "A2M"    "A4GALT" "A4GNT" 
# [10] "AAA1" ...
whole_genome <- colnames(RNASeq2.LUAD)[(1+1):(20499+1)]
save(whole_genome, file="whole_genome.Rda")
## paused ###
## 
## 
## 


# split whole cohort in the clinical and expresions of each gene, by gene Symbol, from Peter ###

# whole_genome <- unique(as.character(TA51_HSCORE$Symbol))
whole_genome <- as.character(TA51_HSCORE$Symbol[1:20499])
# save(whole_genome, file="whole_genome.Rda")
# or
load(file="whole_genome.Rda")
LUAD_n <- length(whole_genome) # n=16907, such as "ZZZ3"   "ZZEF1"  "ZYX"    "ZYG11B" "ZYG11A" "ZXDC"  .....


# nodata <- as.data.frame(setNames(replicate(5,numeric(0), simplify = F), letters[1:5]))
# a=data.frame(matrix(NA, nrow=100, ncol=2)) ; names(a) <- c("SNP","P_value")

beginning <- 1
# keep going from "RASL11A" at 5640
# beginning <- 5640
for (i in beginning:LUAD_n) {
  #  LUAD <- data.frame(matrix(NA, nrow=1, ncol= ncol((TA51_HSCORE)))) # reset
  #  names(LUAD) <- colnames(TA51_HSCORE) # reset
  LUAD <- data.frame()
  attach(TA51_HSCORE)
  # https://www.statmethods.net/management/subset.html
  #  LUAD <- subset(TA51_HSCORE, select=c("Unique.ID", "H.score_N", "H.score_T", "Gender")) # subset by variables
  LUAD <- TA51_HSCORE[which(Symbol == whole_genome[i]), ] # filtering by observations
  detach(TA51_HSCORE)
  
  duplicated(LUAD$patient_ID) # passed [all should be FALSE], n = 100
  
  save(LUAD, file = paste("TCGA_LUAD_100_clinical_fullFeatures10_", whole_genome[i], ".Rda", sep = ""))
}




#
# [x Load HSCORE from Tableau] ####
# transfer from A1-A4 B1-B4 ... to individual case N-T-T-T
# 
install.packages("googledrive")
library("googledrive")
library(readr)
install.packages("xml2")
library("xml2")

#TA51_HSCORE <- read_xml("https://drive.google.com/open?id=0B16FrnckS21aUEZKNDExSUFrck0")

# direct access google drive
#drive_auth(cache=T)
#drive_get(id="0B16FrnckS21aLThqZkRaRkZqXzg") # Enter authorization code: 4/zEhM0LaXwO8HjfT4QHNbJ438GmYs6tH_SDgoqWMYrMw

##
# split them into 3Mb chunk in BASH command: split
# https://eikhart.com/blog/autosplit-csv
# split -a a -l 50722 LUAD_survival_Tablau_raw.csv ALPS # split by line number
# -p ZZZ3
# retained their header by script https://gist.github.com/steezeburger/98114746b2e4c5fa1ad1

# reCsvEdit is working to fix the hearder error from Tableau: LUAD_survival_Tablau_unicode.csv
# http://recsveditor.sourceforge.net/

# (OK) TA51_HSCORE_ALPS001 <- read.csv(file.choose(), header=T, nrows = -1, sep=",") # fileEncoding = "UTF-8", encoding= "UTF-8") # stdin()
TA51_HSCORE <- read.csv(file.choose(), header=T, skip =0 , nrows = -1, sep="") # "tabs" seperated
save(TA51_HSCORE, file="TCGA_LUAD_100_clinical_fullFeatures10_gene16907.Rda")


#
# split whole cohort in the clinical and expresions of each gene, by gene Symbol ###
load(file= "TCGA_LUAD_100_clinical_fullFeatures10_gene16907.Rda")

# whole_genome <- unique(as.character(TA51_HSCORE$Symbol)) # n= 42652 ?? why ??
whole_genome <- as.character(TA51_HSCORE$Symbol[1:16907])
# save(whole_genome, file="whole_genome.Rda")
# or
load(file="whole_genome.Rda")
LUAD_n <- length(whole_genome) # n=16907, such as "ZZZ3"   "ZZEF1"  "ZYX"    "ZYG11B" "ZYG11A" "ZXDC"  .....


# nodata <- as.data.frame(setNames(replicate(5,numeric(0), simplify = F), letters[1:5]))
# a=data.frame(matrix(NA, nrow=100, ncol=2)) ; names(a) <- c("SNP","P_value")

beginning <- 1
# keep going from "RASL11A" at 5640
# beginning <- 5640
for (i in beginning:LUAD_n) {
#  LUAD <- data.frame(matrix(NA, nrow=1, ncol= ncol((TA51_HSCORE)))) # reset
#  names(LUAD) <- colnames(TA51_HSCORE) # reset
  LUAD <- data.frame()
  attach(TA51_HSCORE)
# https://www.statmethods.net/management/subset.html
#  LUAD <- subset(TA51_HSCORE, select=c("Unique.ID", "H.score_N", "H.score_T", "Gender")) # subset by variables
  LUAD <- TA51_HSCORE[which(Symbol == whole_genome[i]), ] # filtering by observations
  detach(TA51_HSCORE)
  
  duplicated(LUAD$patient_ID) # passed [all should be FALSE], n = 100
  
  save(LUAD, file = paste("TCGA_LUAD_100_clinical_fullFeatures10_", whole_genome[i], ".Rda", sep = ""))
}



#
H <- unlist(TA51_HSCORE)
HB <- as.data.frame(matrix(H, ncol = 4, byrow=T))
colnames(HB) <- c("N","T_a","T_b","T_c")
#write.csv(HB, stdout(), row.names=FALSE) # => paste into clinicopathological tables
#
#
#
# Read table from .csv ##
## TMU HNSCC TA51BCDE dataset, n=82 (Jan 2017)
# check Uniqueness of ID
#ID_HSCORE <- read.csv(stdin(),header=T)
# or
# # > colnames(ID_HSCORE_clinico) # N-TTT 4 cores as one case
# [1] "TMA" (core name)                   "Unique.ID"              
# [3] "H.score_N"               "H.score_Ta"             
# [5] "H.score_Tb"              "H.score_Tc"             
# [7] "H.score_T"               "Patho_Diagnosis"        
# [9] "primary.site"            "patho"                  
# [11] "Differ"                  "margin"                 
# [13] "LVI"                     "PNI"                    
# [15] "ECM"                     "X"                      
# [17] "Gender"                  "age.at.diagnosis"       
# [19] "T"                       "N"                      
# [21] "M"                       "stage"                  
# [23] "stage_2"                 "surgery"                
# [25] "procedure"               "R.T"                    
# [27] "C.T"                     "dead"                   
# [29] "death.related.to.SCCHN." "death.date"             
# [31] "latest.F.U"              "cencored.date"          
# [33] "OS..months._from.biopsy" "OS.from.op"             
# [35] "RFS..months._from.op"    "recurrence"             
# [37] "date.of.recurrence"      "recurrent.site"         
# [39] "X_Censored"              "Last_FU_date"           
# [41] "surgery.date"            "X1st.R.T.date"          
# [43] "X1st.R.T.total.dose"     "X1st.C.T.date"          
# [45] "date.of.recurrence.1"    "date.of.death"          
# [47] "cause.of.death"





#==
#x Process LUAD cohort into individual gene ####
load(file="LUAD_TCGA2STAT.Rda")   # as LUAD
load(file="whole_genome.Rda") # the name list of protein coding genome

ID_clinico <- LUAD$merged.dat[, c(1:3, 4:5, 6,7:9, 10, 20510:20511)] # dim 515 x 12
colnames(ID_clinico)[1] <- "Unique.ID" # rename it
len_ID_clinico <- length(ID_clinico) # e.x. 12

# LUAD$merged.dat[, c(1,11:20509)] # RNAseqV2

i <- beginning <- 11 # gene: "A1BG"
LUAD_n <- 20509 # ending, last gene: "ZZZ3"

while (i<= LUAD_n)  # for (i in beginning:LUAD_n)
  {
  geneName <- colnames(LUAD$merged.dat)[i] # e.g. "A1BG" or "PMM1"
  ID_HSCORE_clinico <- cbind(ID_clinico, LUAD$merged.dat[, i]) # dim 515 x 13
  colnames(ID_HSCORE_clinico)[len_ID_clinico+1] <- geneName
  save(ID_HSCORE_clinico, file=paste("TCGA_LUAD_515_clinical_fullFeatures11_", geneName, ".Rda", sep=""))
  i <- i + 1
  }
# done







## Peter's .xlsx ## archived
# colname processing
# > colnames(ID_HSCORE_clinico)
# [1] "X_INTEGRATION"    "Unique.ID"       "Symbol"           "X_OS"            
# [5] "X_OS_IND"         "X_RFS"            "X_RFS_IND"        "Expression"      
# [9] "gender.M.F."      "N.T"   [1:T;0:N]           "pathologic_M"     "pathologic_N"    
# [13] "pathologic_stage" "pathologic_T" 
# 
 # as.data.frame(table(ID_HSCORE_clinico$N.T)) # frequency table
 # Var1 Freq
 # 1    0    7 => paired normal
 # 2    1   60 => tumour :-)
 # 3    2    6
 # 4    3    5
# 
# ID_HSCORE_clinico <- read.csv("/Users/Apple/Documents/Tex_proposal_MHlab/TissueMicroarray_HNSCC_TA51/TA51_OSCC_clinicopathological_dataset_PMM1_IHC/clinical_clean-Table 1.csv", header=T) 
duplicated(ID_HSCORE_clinico$Unique.ID) # passed [all should be FALSE]

# NT_B <- ID_HSCORE_clinico[complete.cases(ID_HSCORE_clinico[, c(2,3,7)]), c(2,3,7)] # Id and HSCORE NTTT
# # n=70 (TA51BCDE), there is NO censored cases in TA51 cohort

#osccT <- ID_HSCORE_clinico[complete.cases(ID_HSCORE_clinico[, c(2,7)]), c(2,3,7,17,18,19:21,23,26,27,33,35,36)]
# or
# # %%% # includes [12] surgical_margin (+) as c(0.5, 1, 2) to DO survival curve again !
osccT <- ID_HSCORE_clinico[complete.cases(ID_HSCORE_clinico[, c(2,8)]), c(2,9,8,4,5,6,7,14,12,11,13)]
osccT[, 1] <- sapply(osccT[, 1], as.character) # convert from factor to character

save(osccT, file="TMU_TA51BCDE_T70_clinical_fullFeatures13.Rda")
# or
save(osccT, file="TMU_TA51BCDE_T70_clinical_fullFeatures14_margin.Rda")


#
#
#
#
#
#
# [matched N-T pair, n=35] Analysis of Matched Tumor and Normal tissues
# ex. "For PMM1 gene, we found expression to be higher in adjacent normal samples than
# tumors in most cases (p-value???=???0.007, paired t-test)"
TA51B_HSCORE <- read.csv(stdin(), header=T)
TA51B_HSCORE <- cbind(TA51B_HSCORE, read.csv(stdin(), header=T))
TA51B_HSCORE <- cbind(read.csv(stdin(), header=T), TA51B_HSCORE)

NT_B <- TA51B_HSCORE[complete.cases(TA51B_HSCORE), ]

# to identify normal and tumor samples. this is done using the TCGA barcode (https://wiki.nci.nih.gov/display/TCGA/TCGA+barcode). 
# The two digits at position 14-15 of the barcode will indicate teh sample type, from the link:
# Sample barcode: TCGA-02-0001-01
#   "Tumor types range from 01 - 09, normal types from 10 - 19 and control samples from 20 - 29."



# [2.(optional)Prepareing NT pair]  ####
# NT pair plot ---------------------Using graphs to display matched pairs data ####
### osccNT: NT pair ###
## matched N T statistics: Using the paired t-test
#library(MASS)
t.test(NT_B$H.score_N, NT_B$H.score_T, paired=T) 
# Result: p-value = 0.6024,
# N minus T by mean of the differences as -3.257143 with 95%CI [-15.843704, 9.329419]
#  # scatter plot
plot(NT_B$H.score_N, NT_B$H.score_T, xlab="Adjacent Normal", ylab="Tumour")
abline(lm(NT_B$H.score_N ~ NT_B$H.score_T))
## The expression levels of PMM1 mRNAs in OSCC patients and controls were compared using Mann-hitney U-test
# The Wilcoxon-Matt-Whitney test (or Wilcoxon rank sum test, or Mann-Whitney U-test) is used when is asked to compare the means of two groups that do not follow a normal distribution: it is a non-parametrical test. 
# is the equivalent of the t test, applied for independent samples.
# wilcox.test(mtcars$mpg, mtcars$am, correct=FALSE)
# wilcox.test(mpg ~ am, data=mtcars) 
#we conclude by accepting the hypothesis H0 of equality of means, if p-value > 0.05
#
#
install.packages("PairedData")
library(PairedData)
#data(PrisonStress)
#paired.plotProfiles(PrisonStress,"PSSbefore","PSSafter", subjects="Subject",groups="Group")
#paired.plotMcNeil(PrisonStress,"PSSbefore","PSSafter", subjects="Subject")

paired.plotProfiles(NT_B, "H.score_N", "H.score_T")
paired.plotMcNeil(NT_B, "H.score_N", "H.score_T", subjects="Unique.ID")


# TCGABiolinks - A R/Bioconductor package to search, download and prepare relevant data for analysis in R. Very powerful and well documented.





# [3.prepare T] OS and RFS
# by CDE: https://gdc.cancer.gov/clinical-data-elements
# or https://docs.gdc.cancer.gov/Data_Dictionary/viewer/#?view=table-definition-view&id=follow_up
# [3.prepare T] Binary/binomial of clinicalMatrix (clinicopathological features) ####
# or categorize, or factorize or dichotomize by recoding
# osccT_fileName <- load(file.choose()) # LUAD, "TCGA_LUAD_515_clinical_fullFeatures11_ZZZ3.Rda" with residual tumour status
# or transfer from LUAD

load(file.choose()) # LUAD.clinical.Fire.Rda
# or
load(file="LUAD.clinical.Fire.Rda")
osccT <- LUAD.clinical.Fire # n=522; #LUAD.clinico_mRNA.Fire

# ignore patient_ID (it is NOT true)
# https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/ data schema
# https://wiki.nci.nih.gov/display/TCGA/TCGA+barcode
# ID_table <- as.data.frame(table(LUAD.clinico_mRNA.Fire[,c(1,54, 86)]))
# TCGA LUAD n=576
# 
# #https://www.biostars.org/p/153013/ a step-by-step tutorial,
# Survival analysis of TCGA patients integrating gene expression (RNASeq) data 
source("https://bioconductor.org/biocLite.R")
biocLite("curatedTCGAData") # utf8 and pillar and tibble
library("curatedTCGAData")
biocLite("TCGAWorkflow") #https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5302158/
# GenomicDataCommons - A R/Bioconductor package for querying, accessing, and mining genomic datasets available from the GDC.
options(repos="http://cran.csie.ntu.edu.tw/") # or http://cran.ism.ac.jp/; yzu is working :-)
biocLite("GenomicDataCommons") # tibble, https://github.com/Bioconductor/GenomicDataCommons

# UCSC Xena: https://genome-cancer.ucsc.edu/proj/site/composite/datapages/?dataset=TCGA.LUAD.sampleMap/LUAD_clinicalMatrix&host=https://tcga.xenahubs.net
# list of all Phenotypes CDE of the LUAD dataset=>
# https://genome-cancer.ucsc.edu/proj/site/composite/datapages/?host=https%3A%2F%2Ftcga.xenahubs.net&dataset=TCGA.LUAD.sampleMap%2FLUAD_clinicalMatrix&label=Phenotypes&allIdentifiers=true
# _OS
# _OS_IND
# _OS_UNIT
#_RFS
# _RFS_IND
# _RFS_UNIT
#primary_therapy_outcome_success
# progression_determined_by
# person_neoplasm_cancer_status => yes we have it in LUAD of TCGA
# new_neoplasm_event_type
# new_tumor_event_after_initial_treatment
# days_to_new_tumor_event_after_initial_treatment
# days_to_additional_surgery_locoregional_procedure
# days_to_additional_surgery_metastatic_procedure

# days from the end of observation to ...
# [18] "days_to_index"  ??                            
# [19] "days_to_initial_pathologic_diagnosis" 
# (baseline: the date of initial pathologic diagnosis)    
# [20] "days_to_last_followup"                      
# [21] "days_to_last_known_alive"  
# [32] "followup_treatment_success"
# [56] "person_neoplasm_cancer_status" # : "tumor free" or "with tumor"
# [61] "primary_therapy_outcome_success" 


#
# https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/sample-type-codes; GDC Sample Type Codes (data schema)
# # > colnames(LUAD.clinico_mRNA.Fire)
# https://www.bioconductor.org/packages/devel/data/experiment/manuals/curatedOvarianData/man/curatedOvarianData.pdf (data schema)
# In the "Available sample meta-data" sections:
# [1] "tcga_participant_barcode"                   
# [2] "age_at_initial_pathologic_diagnosis"        
#...                    
# [13] "date"                                       
# [14] "day_of_dcc_upload"                          
# [15] "day_of_form_completion"                     
# [16] "days_to_birth"                              
# [17] "days_to_death"                              
#...
# [20] "days_to_last_followup"         # alive      
# => Days To Recurrence = RFS time        
# [21] "days_to_last_known_alive"                   
#..     
# [29] "ethnicity"                                  
# [30] "file_uuid"                                  
# ..
# [32] "followup_treatment_success"                 
# [33] "gender"                                     
# [34] "histological_type"                          
# [35] "history_of_neoadjuvant_treatment"           
# [36] "icd_10"                                     
# ...                
# [45] "lost_follow_up"                             
# [46] "month_of_dcc_upload"                        
# [47] "month_of_form_completion"                   
# [48] "number_pack_years_smoked"                   
# [49] "other_dx"                # co-morbility                   
# [50] "pathologic_m"                               
# [51] "pathologic_n"                               
# [52] "pathologic_stage"                           
# [53] "pathologic_t"                               
# [54] "patient_id" #???      # four digital     , n=516                     
# [55] "performance_status_scale_timing"            
# [56] "person_neoplasm_cancer_status"         # tumour free or with tumour  (RFS index)   
# ...         
# [61] "primary_therapy_outcome_success": completeresponse|partialresponse|progressivedisease|stabledisease:
#response to any kind of therapy (including radiation only).        
# [62] "progression_determined_by.3"        # positive biomarker(s)        
# [63] "pulmonary_function_test_performed"          
# [64] "race"                                       
# [65] "radiation_therapy"                          
# [66] "residual_tumor"           # residual tumor, or residual tumour     
# # (in detail:  http://www.ilincs.org/GenomicsPortals/differentialExpressionSetup.do?data_set=TCGA_LUAD_RNASeqV2&db=tcgaIlluminaHiSeq)
             
# [67] "stopped_smoking_year"                       
# ...                         
# [73] "tobacco_smoking_history"                    
# ...                         
# [76] "vital_status"      # 363 alive or 213 dead  # > table(LUAD.clinico_mRNA.Fire[,76]) # for this summary                       
# ...                         
# [79] "year_of_form_completion"                    
# [80] "year_of_initial_pathologic_diagnosis"       
# [81] "year_of_tobacco_smoking_onset"   
# 
# belowing is mRNA.Fire:           
# [82] "gene"     # ZZZ3                                  
# [83] "expression_log2"                            
# [84] "z.score"                                    
# [85] "cohort.y"           #= LUAD                        
# [86] "sample_type"          # 59 NT(Solid Tissue Normal) 515 TP(Primary Solid Tumor) or 2 TR(Recurrent Solid Tumor)
# [87] "protocol"                       # RSEM                                
# [88] "geneID"   # 26009 for ZZZ3
# 
# 




# [categorizing the survival data from clinical features] ####
# from osccT to oscc
# i= from 20499(ZZZ3) to 1(A1BG)
oscc <- data.frame()
oscc <- subset(osccT, select=c("tcga_participant_barcode", "gender")) # gender: male as 1, female as 2
# H_score as ZZZ3 RNAseq level
oscc$gender[oscc$gender == "male"] <- 1 # male as 1, recoding
oscc$gender[oscc$gender == "female"] <- 2 # female as 2
# oscc <- cbind(oscc, subset(osccT, select= c("OS..months._from.biopsy", "OS..months._from.biopsy_IND", "X_RFS", "X_RFS_IND")))
#
oscc$ageDx[osccT$age_at_initial_pathologic_diagnosis <= 65] <- 1 # younger
oscc$ageDx[osccT$age_at_initial_pathologic_diagnosis > 65] <- 2 # older#
# or create 2 age categories  by
# mydata$agecat <- ifelse(mydata$age > 65, 
#                        c("older"), c("younger")) 
# mydata$Agecat4<-cut(mydata$Age, seq(0,30,5), right=FALSE, labels=c(1:6))



# keep grouping and Recoding variable by Factor
# https://stackoverflow.com/questions/3875608/grouping-recoding-factors-in-the-same-data-frame
# > levels(as.factor(osccT$PATHOLOGYTSTAGE))
# [1] "t1"  "t1a" "t1b" "t2"  "t2a" "t2b" "t3"  "t4"  "tx"
# %in% is "within"
oscc$pathologic_T[osccT$pathologic_t %in% c("t1",  "t1a", "t1b", "t2",  "t2a", "t2b")] <- 1 # T1 T2
oscc$pathologic_T[osccT$pathologic_t %in% c("t3",  "t4")] <- 2 # T3 T4
oscc$pathologic_T[osccT$pathologic_t %in% c("tx", NA)] <- NA # unknown T
#
#> levels(as.factor(osccT$PATHOLOGYNSTAGE))
# [1] "n0" "n1" "n2" "n3" "nx"
oscc$pathologic_N[osccT$pathologic_n == "n0"] <- 0 # N0
oscc$pathologic_N[osccT$pathologic_n %in% c("n1", "n2", "n3")] <- 1 # N+ (N1,2,3)
oscc$pathologic_N[osccT$pathologic_n %in% c("nx", NA)] <- NA # unknown N

#> levels(as.factor(osccT$PATHOLOGYMSTAGE))
# [1] "m0"  "m1"  "m1a" "m1b" "mx" in LUAD
# M as 0 or 2 in TMU TA51 (why M2 instead of M1?)
oscc$pathologic_M[osccT$pathologic_m == "m0"] <- 0 # M0
oscc$pathologic_M[osccT$pathologic_m %in% c("m1",  "m1a", "m1b")] <- 1 # M+
oscc$pathologic_M[osccT$pathologic_m %in% c("mx", NA)] <- NA # unknown M
#
#> levels(as.factor(osccT$PATHOLOGICSTAGE))
# [1] "stage i"    "stage ia"   "stage ib"   "stage ii"   "stage iia"  "stage iib" 
# [7] "stage iiia" "stage iiib" "stage iv"
oscc$stage[osccT$pathologic_stage %in% c("stage i",    "stage ia",   "stage ib",   "stage ii",   "stage iia",  "stage iib")] <- 1 # stage 1 2
oscc$stage[osccT$pathologic_stage %in% c("stage iiia", "stage iiib", "stage iv")] <- 2 # stage 3 4
#


# ** margin issue: 0 or 1 or NaN ####
# with r1 or r2 (residual tumor) => person_neoplasm_cancer_status (with tumor) is NOT recurrence
# while clincal "tumor free" cases with r1 or r2 => Amazing ! selflimited :-)
# residual_tumor in LUAD: The status of a tissue margin following surgical resection:
#=> r0, 無殘留，tumor free from margin (就是有 「開乾淨」的意思)
#=> r1 (micro) r2 (macro) 有殘留，
#=> rx 是 unknown status。
# r0  r1  r2  rx 
# 371  15   4  27 
# [12] surgical_margin (+) in HNSCC; definition as c(0, 0.5, 1) [0.5 is close margin]; there is NO NaN.
oscc$margin[osccT$residual_tumor %in% c("r0")] <- 0 # margin free from tumour 0, n=348 
# NO close margin 0.5
oscc$margin[osccT$residual_tumor %in% c("r1", "r2")] <- 1 # positive surgical margin 1, n=18
# rx => NaN, n=156



# X NaN was impuated by oscc[17,9]<- 1; oscc[22,9] <- "1" in TA51 HNSCC

# censoring data: Since we want to do censored analysis, we need to have
# something to censor the data with. For example, if a patient has no death data
# BUT there is a date to last followup it means that after that day we know
# nothing about the patient, therefore after that day it cannot be used for
# calculations/Kaplan Meier plot anymore, therefore we censor it. so now we need
# to create vectors for both 'time to new tumor' and 'time to death' that
# contain also the data from censored individuals.
# https://www.biostars.org/p/153013/


# [OS] overall survival ####
# # OS and RFS 
# (at birth [16] -> healthy -> cancer diagnosed [19] -> therapies -> [recurrence event] -> death event)
# GDC data dictionary viewer for CDE: https://cdebrowser.nci.nih.gov/cdebrowserClient/cdeBrowser.html#/search
# rule and definition by colnames number [1:81] of LUAD.clinical.Fire `: 
# OS: OS_IND Event <- death or alive [76] "vital_status"
# "vital_status" == 1 as death; == 0 as alive
oscc$OS_IND[osccT$vital_status == "dead"] <- 1 # 1==death event (dead), n=213
oscc$OS_IND[osccT$vital_status == "alive"] <- 0 # 0==no event (alive), 363 (censored)

# [19] "days_to_initial_pathologic_diagnosis" 
# (baseline: the date of initial pathologic diagnosis) => all are NA or zero :-)
# OS Time to event <- days to death [17] or Max of (last followup [20]/ last known alive [21]); 
oscc <- cbind(oscc, subset(osccT, select=c("days_to_death")))
# max of case c(15,412,455)
osccT_2021 <- osccT[c(20,21)]
#OS_live_time_ind <- apply(osccT_2021, 1, which.max) # 1: in a row-wise matter
#osccT_2021[OS_live_time_ind] 
# ifelse function: ifelse(test, yes, no)
# => try to pick up max(osccT[c(20,21)], na.rm=T) by row-wise matter
OS_live_days <- apply(osccT_2021, 1, function(x) ifelse(all(is.na(x)), NA, max(x, na.rm=T)))
OS_live_days <- as.data.frame(OS_live_days) # R help from Jorge I Velez

# at data cleaning stage: removal of 3-(NAs or 0) # OStime <- (LUAD.clinical.Fire[c(76,17,19, 20,21)])
# to calculate overall survival time (months) of TCGA LUAD? months until death
oscc$OS_months <- oscc$days_to_death / 30
oscc$OS_live_months <- OS_live_days /30 # derived from last followup [20]/ last known alive [21])
#=> if there is no event (alive) => censored (status 0)
# "OS_live_months" with label.var as "OS_live_days"
names(oscc$OS_live_months) <- ""

# OS Censoring: create vector time_to_death containing values to censor (oscc$OS_live_months) for death
# oscc$censored_OS[!is.na(osccT$OS..months._from.biopsy)] <- 0 # censored -
# oscc$censored_OS[is.na(osccT$OS..months._from.biopsy)] <- 1 # censored +
# ifelse(test, yes, no):
for (i in (1:nrow(oscc))){ # for (i in beginning:LUAD_n)
oscc[i,11] <- ifelse (is.na(oscc[i,11]),
                      oscc[i,12], oscc[i,11])
# 11: oscc$OS_months
# 12: oscc$OS_live_months
# as.numeric(as.character())
}

# data("mgus2") # https://cran.r-project.org/web/packages/survival/survival.pdf
# etime <- with(mgus2, ifelse(pstat==0, futime, ptime)) #pstat=PCM: 0=no (OStime) or 1=yes (RFStime)
# event <- with(mgus2, ifelse(pstat==0, 2*death, 1)) #0:PCM 0+2*alive(0); 2:PCM 0+2*death(1); 1:PCM 1
# event <- factor(event, 0:2, labels=c("censor", "pcm", "death")) # level 0:censor 2:death 1:PCM
# View(Surv(etime, event)) # censoring status= 0 or "censor"; while 1 as "event"
# # ptime=RFStime
# time until progression to a plasma cell myeloma (PCM) or last contact, in months
# pstat=tumor event
# occurrence of PCM: 0=no, 1=yes
# futime=OStime
# time until death or last contact, in months
# death=death event
# occurrence of death: 0=no, 1=yes
# 
# within function: calculate a new variable
# mgus2 <- within(mgus2, {delta <- abs(death - 1) }) # In the “mgus2” data, we want to define a new censoring variable
# "delta" which takes the values 1 for a censored variable and 0 for an death event
# 
# # assiging labels to variable name:
names(oscc$OS_live_months) <- ""
library(Hmisc)
label(oscc$OS_live_months) <- ""
# >OS_live_days 
# set.seed(22)
# data <- data.frame(age = floor(rnorm(6,25,10)), 
#                    sex = gl(2,1,6, labels = c("f","m")))
# var.labels <- c(age = "Age in Years", 
#                 sex = "Sex of the participant")
# dplyr::as.tbl(data)
# data <- Hmisc::upData(data, labels = var.labels)
# Hmisc::label(data)
# Hmisc::contents(data)
## end of OS ##





#
# [RFS]: ####
# # (os vs rfs) there is no RFS time in LUAD [2018/01/27]
# however, we must to duplicate OS to RFS for convinence of code running :-)
oscc$RFS_IND <- oscc$OS_IND
oscc$RFS_months <- oscc$OS_months
# jumpt to [preparedT]

# or There is true RFS calculation below:
# impute missing RFS by OS time; OS RFS Interpretation
# https://www.biostars.org/p/153013/
# https://groups.google.com/forum/m/#!msg/ucsc-cancer-genomics-browser/YvKnWZSsw1Q/3IAkkEMyFa4J

# => A notice that this is a competitive risk problem, where, although a patient can recur and then die, if a patient is dead, it will not recur, therefore is more accurate to censor for death events.
# {
RFStime <- (LUAD.clinical.Fire[c(76,17,20,21,66,56)])
# months of tumour free and alive by using "person_neoplasm_cancer_status": n=330 "tumor free" or n=191 "with tumor"
# https://www.researchgate.net/post/Are_there_differences_between_progression-free_survival_relapse-free_survival_and_recurrence-free_survival2
# Recurrence-free survival (RFS) includes (1) any recurrence (local or regional
# [including invasive ipsilateral tumor and invasive locoregional tumor], or
# distant) and (2) death due to any cause (both LUAD and non-LUAD causes of death).
# to consider death to be one of your outcomes. So your event becomes "Event OR All-cause Mortality". 

# (tumour free or new tumor event = residual[66]/progression, recurrence or new/second primary malignacies)
# RFS: RFS_IND Event by "tumor free" or "with tumor" [56]
# RFS Time to event = "Days To Recurrence", Time interval from the date of disease recurrence to the date of initial pathologic diagnosis, represented as a calculated number of days.
# (https://cdebrowser.nci.nih.gov/cdebrowserClient/cdeBrowser.html#/search?publicId=3008295&version=1.0)
# = Max of (known a tumour even alive / days to tumour recurrence); no event => censored as OS time 

#Calculation of RFS (recurrence): between 6 months and 5 years in HNSC
# DFS is not the same with RFS
#Residual tumor [66] => in HNSC tumors patients
#re-diagnosed with a tumor ("with tumor") within 6 months of treatment completion it is
#considered part of the original tumor, not a new tumor, so it is residual tumor [66] and not a recurrence.
# (or pathologic report: margin is NOT free) [66]
#
#For HNSC at least a re-diagnosed patient is considered to have a recurrence of the original tumor
#if it's between 6 months and 5 years, and 5 years post-treatment it is considered a new primary tumor.
# }
# 


# [19] "days_to_initial_pathologic_diagnosis" 
# (baseline: the date of initial pathologic diagnosis)  => all are NA or zero :-)
# > View(as.data.frame(osccT[,c(16:17,20, 56, 76)])) 
# View(oscc[(oscc$OS_IND == 1 & oscc$months_alive >0), c(1,11,13)])
oscc$RFS_months <- osccT$days_to_last_followup / 30 # days_to_last_followup > 0 => all are alive => RFS time

oscc$RFS_IND[osccT$person_neoplasm_cancer_status == "with tumor"] <- 1 # 1==tumour recurrency, event
oscc$RFS_IND[osccT$person_neoplasm_cancer_status == "tumor free"] <- 0 # 0==NO recurrency
# 
#In R, the operators "|" and "&" indicate the logical operations OR and AND.
#
# RFS Censoring: create vector with time to new tumor containing data to censor for new_tumor
# oscc$censored_RFS[!is.na(osccT$RFS..months._from.op)] <- 0 # censored -
# oscc$censored_RFS[is.na(osccT$RFS..months._from.op)] <- 1 # censored +
# end of RFS
# 
# 

# to [preparedT] here ####
# oscc <- cbind(oscc, subset(osccT, select=c("expression_log2", "z.score")))
# save(oscc, file="LUAD_T576_clinical_fullFeatures13_dichotomized.Rda") # with "margin"
# 
oscc <- oscc[,!(colnames(oscc) %in% c("days_to_death", "OS_live_months"))] # removal of junk columns

save(oscc, file="LUAD_T522_clinical_fullFeatures11_dichotomized.Rda") # with "margin" at column 8




# [4] get all RNAseq by FirebrowseR ####
#LUAD.mRNA.Exp.Fire <- Samples.mRNASeq()
## # set path on google drive
library(FirebrowseR)
path_LUAD <- "/Users/apple/Google\ Drive/2018_PhD_dissertation/LUAD_Peter_survival/"
# path_LUAD <- "/Users/apple/Documents/My\ Tableau\ Repository/Workbooks/"
setwd(path_LUAD) # change working directory to the google drive
load(file="LUAD.clinical.Fire.Rda")
load(file="whole_genome.Rda") # the name list of protein coding genome

#http://mazamascience.com/WorkingWithData/?p=912 =>
# warning(…) — generates warnings
# stop(…) — generates errors
# suppressWarnings(expr) — evaluates expression and ignores any warnings
# tryCatch(…) — evaluates code and assigns exception handlers
#suppressWarnings(
#  { #my expr, non-stop run and ignores any warnings

    
# tryCatch.Rscript -- experiments with tryCatch ### 
#for error handling during Samples.mRNASeq (missing genes)
    
    # Get any arguments
    arguments <- commandArgs(trailingOnly=TRUE)
    a <- arguments[1]
    
    # Define a division function that can issue warnings and errors
    # dummy function :-)
    myDivide <- function(d, a) {
      if (a == 'warning') {
        return_value <- 'myDivide warning result'
        warning("myDivide warning message")
      } else if (a == 'error') {
        return_value <- 'myDivide error result'
        stop("myDivide error message")
      } else {
        return_value = d / as.numeric(a)
      }
      return(return_value)
    }

####
# save them in each .Rda
#while (i>= beginning)  # for (i in beginning:LUAD_n) , start from last one ZZZ3
LUAD_n <- length(whole_genome)
beginning <- 1 # = 1 first gene: "A1BG"
#i <- LUAD_n # ending, = 20499 last gene: "ZZZ3"
    
for (i in seq(LUAD_n,beginning, by=-1) )  
    {
      # Evalute the desired series of expressions inside of tryCatch
      result <- tryCatch({
        
        diff.Exp.Genes <- whole_genome[i] # c("ZZZ3") #, "ESR1", "GATA3", "XBP1", "FOXA1", "ERBB2", "GRB7", "EGFR", "FOXC1", "MYC")
        all.Found = F
        page.Counter = 1
        LUAD.mRNA.Exp.Fire = list()
        page.Size = 600 # using a bigger page size is faster
        while(all.Found == F)
        {
          LUAD.mRNA.Exp.Fire[[page.Counter]] <- Samples.mRNASeq(format = "csv",
                                                                gene = diff.Exp.Genes,
                                                                cohort = "LUAD",
                                                                tcga_participant_barcode =
                                                                  LUAD.clinical.Fire$tcga_participant_barcode,
                                                                page_size = page.Size,
                                                                page = page.Counter)
          
          #return code of Retrieve mRNASeq data: 
          #Error in download.Data(url, format,
          #page) : No samples matching your query; there is no 19423, 18757 ... (i)
          #try block: This skips over the error-causing non-numeric input with an error message
          #(you can suppress the error message with the silent=T argument to try), and
          #continues on with the rest of the input.
          
          if(nrow(LUAD.mRNA.Exp.Fire[[page.Counter]]) < page.Size)
            all.Found = T
          else
            page.Counter = page.Counter + 1
          # Warning: In strptime(x, fmt, tz = "GMT") :
          #   unknown timezone 'zone/tz/2017c.1.0/zoneinfo/Asia/Taipei'
          
        }
        LUAD.mRNA.Exp.Fire <- do.call(rbind, LUAD.mRNA.Exp.Fire)
        # dim(LUAD.mRNA.Exp.Fire) # LUAD.mRNA.Exp.Fire$expression_log2 or LUAD.mRNA.Exp.Fire$z.score
        # 
        save(LUAD.mRNA.Exp.Fire, file=paste("LUAD.mRNA.Exp.", diff.Exp.Genes, ".Fire.Rda", sep=""))
        #  } # end of real work
        #  
        #  
        b <- 2
        c <- b^2
        d <- c+2
        # if (a == 'suppress-warnings') {
        #   e <- suppressWarnings(myDivide(d,a))
        # } else {
        #   e <- myDivide(d,a) # 6/a
        # }
        e <- d
        f <- e + 100
        
      }, warning = function(war) {
        
        # warning handler picks up where error was generated
        print(paste("MY_WARNING:  ",war))
        b <- "changing 'b' inside the warning handler has no effect"
        e <- myDivide(d,0.1) # =60
        f <- e + 100
        print(b)
        return(f)
        
      }, error = function(err) {
        
        # error handler picks up where error was generated
        print(paste("MY_ERROR:  ",err))
        b <- "changing 'b' inside the error handler has no effect"
        e <- myDivide(d,0.01) # =600
        f <- e + 100
        print(b)
        return(f)
        
      }, finally = {
        
        print(paste("a =",a))
        print(paste("b =",b))
        print(paste("c =",c))
        print(paste("d =",d))
        # NOTE:  Finally is evaluated in the context of of the inital
        # NOTE:  tryCatch block and 'e' will not exist if a warning
        # NOTE:  or error occurred.
        #print(paste("e =",e))
        
      }) # END tryCatch
      
      print(paste("i =", whole_genome[i]))
      
} # end of for loop
    
#  } # end of my expr
#  ) # end of suppressWarnings
# done
# skip .....
# # https://www.biostars.org/p/153013/
# # addiontally, try to find the genes whose expression is == 0 in more than 50% of the samples:
# i<-1 # "A1BG"
# load(file=paste("LUAD.mRNA.Exp.", whole_genome[i], ".Fire.Rda", sep=""))
# rna <- LUAD.mRNA.Exp.Fire[,1:4]
# trash_check <- function(x){
#   x <- as.matrix(x)
#   x <- t(apply(x,1,as.numeric))
#   r <- as.numeric(apply(x,1,function(i) sum(i == 0)))
#   ind <- which(r > dim(x)[2]*0.5)
#   return(ind)
# }
# trash_ind <- trash_check(rna)
# trash_rna <- rna[trash_ind,]
# rna <- rna[-trash_ind,]
# #




#ggplot2 package to plot the expression
library(ggplot2)
p = ggplot(mRNA.Exp, aes(factor(gene), z.score))
p +
  geom_boxplot(aes(fill = factor(sample_type))) +
  # we drop some outlier, so plot looks nicer, this also causes the warning
  scale_y_continuous(limits = c(-1, 5)) +
  scale_fill_discrete(name = "Tissue")



##
#x [4.5 .Load all RNAseq; LUAD.clinico_mRNA.ZZZ3.Fire.Rda] #
while (i>= beginning)  # for (i in beginning:LUAD_n) , start from last one ZZZ3
{
  load(file = paste("LUAD.clinico_mRNA.", whole_genome[i], ".Fire.Rda", sep = "")) # as LUAD.clinico_mRNA.Fire
  #  assign(paste("LUAD.clinico_mRNA.", diff.Exp.Genes, ".Fire", sep=""), LUAD.clinico_mRNA.Fire) # rename it with gene ZZZ3
  
  LUAD <- ID_HSCORE_clinico <- LUAD.clinico_mRNA.Fire # as LUAD
  # [1] "tcga_participant_barcode" "gene"                     "expression_log2"         
  # [4] "z.score"                  "cohort"                   "sample_type"             
  # [7] "protocol"                 "geneID"      
  # jump to [survival analysis]
  # 
  # 
  # 
  # 
  i <- i - 1
}
## Preparation finished (for each cancer type or cohort) ###





# [2018/03/06] with Peter Lai birthday
# Resume:[5a.START] ZZZ3 and TMSB4X ####
# Start the survival analysis for each individual gene
## # set path on google drive
library(FirebrowseR)
TCGA_cohort <- "LUAD" # cancer type
path_LUAD <- "/Users/apple/Google\ Drive/2018_PhD_dissertation/LUAD_Peter_survival/"
# path_LUAD <- "/Users/apple/Documents/My\ Tableau\ Repository/Workbooks/"
setwd(path_LUAD) # set the working directory to the google drive


#load(file="TMU_TA51BCDE_T70_clinical_fullFeatures13_dichotomized.Rda") # as oscc (without maring feature)
# or
load(file="LUAD_T522_clinical_fullFeatures11_dichotomized.Rda") # as oscc with "margin", negative n=348 while positive n=18 and 156:NaN)
#
load(file="whole_genome.Rda") # the name list of protein coding genome
LUAD_n <- length(whole_genome)
geneName <- whole_genome[LUAD_n] # aka. "ZZZ3"


# load and prepare $$.clinico_mRNA.Fire for survival analysis
load(file=paste("LUAD.mRNA.Exp.", geneName, ".Fire.Rda", sep="")) # as LUAD.mRNA.Exp.Fire
# "sample_type"          # 59 NT(Solid Tissue Normal) 515 TP(Primary Solid Tumor) or 2 TR(Recurrent Solid Tumor)

LUAD_T_mRNA.Fire <- LUAD.mRNA.Exp.Fire[LUAD.mRNA.Exp.Fire$sample_type %in% c("TP"), c("tcga_participant_barcode","z.score")] # n=517 -2 = 515
# removing duplicated ID %in%  c("TCGA-50-5066","TCGA-50-5946") => TR recurrent sample

# inner join by merge
LUAD.clinico_mRNA.Fire <- merge(oscc, LUAD_T_mRNA.Fire, by="tcga_participant_barcode") #n=515 with sample_type "TP", excluding "NT normal tissue" or "TR"

# [colnames correction] LUAD(LUAD.clinico_mRNA.Fire) => HNSCC(osccT) "TMU_TA51BCDE_T70_clinical_fullFeatures14_margin.Rda"
# # > colnames(osccT) of HNSCC
# [1] "Unique.ID"               "H.score_N"               "H.score_T"              
# [4] "margin"                  "Gender"                  "age.at.diagnosis"       
# [7] "T"                       "N"                       "M"                      
# [10] "stage_2"                 "R.T"                     "C.T"                    
# [13] "OS..months._from.biopsy" "RFS..months._from.op"    "recurrence" (as RFS_IND)

# > colnames(LUAD.clinico_mRNA.Fire)
# [1] "tcga_participant_barcode" "gender"                  
# [3] "ageDx"                                     
# [4] "pathologic_T"             "pathologic_N"            
# [6] "pathologic_M"             "stage"     "margin"             
# [9] "OS_IND"                   "OS_months"               
# [11] "RFS_IND"                  "RFS_months"                
# [13] expression "z.score" of RNAseq (tumour)

# rename all of colnames as HNSCC's osccT
colnames(LUAD.clinico_mRNA.Fire) <- c("Unique.ID","Gender","age.at.diagnosis",
                                      "T","N",
                                      "M","stage_2","margin",
                                      "OS_IND", "OS..months._from.biopsy",
                                      "RFS_IND", "RFS..months._from.op", "H.score_T")
# n=515 in LUAD
oscc <- LUAD.clinico_mRNA.Fire # starting analysis with "oscc"
# oscc$H.score_T as LUAD.mRNA.Exp.Fire$z.score; expression level: H.score_T as RNAseq z.score

## a dummy universal variable for binomial (high/low) oscc$geneName_median, all are zero
oscc$PMM1_median <-(oscc$H.score_T >= median(oscc$H.score_T, na.rm=T)) +0 # higher 1 or lower 0
osccM_pos <- which(colnames(oscc) == "PMM1_median") # at column 14

#x oscc$geneM <- 0
#x colnames(oscc)[colnames(oscc) == "geneM"] <- paste(geneName, "_median", sep="") # rename as oscc$geneName_median
#x osccM_pos <- which(colnames(oscc) == paste(geneName, "_median", sep="")) # at column 19 :-)




##
### $setup global variables...(including "margin") ####
library("psych") # for describe()
library(survival)
# Create the Sheet title and subtitle; underline as 0 (default value, no underline), 1 (underline with one line), 2 (underline with two lines)
# Table 3. Univariate/Multivariate Cox proportional hazards regression analyses on OS time
# with 8 items:
featuresUni <- c("Gender",
                 "Age at diagnosis",
#                 "Primary site",
#                 "Clinical T Status",
#                 "Clinical N Status",
#                 "Clinical Stage",
                 "Pathologic T status",
                 "Pathologic N status",
                 "Pathologic M status",
                 "Pathologic Stage",
                 "Surgical Margin status",
#                 "Lymphovascular Invasion",
#                 "Perineural Invasion",
#                 "Extranodal spreading of neck LN",
#                 "SCC Histologic Grade",
#                 "Radiotherapy",
#                 "Chemotherapy",
                 paste(geneName, "z-score", sep = " ") # "IHC score" => "z-score"
)
#
# z-score calculation = [(value gene X in tumor Y)-(mean gene X in normal)]/(standard deviation X in normal)
#
colUni_1 <- c("Male",
              ">65y",
#              "higher risk",
#              "T3+T4",
#              "N1-3",
#              "Stage III+IV",
              "T3+T4",
              "N1-3",
              "M1",
              "Stage III+IV",
              "Positive", # positive safety margin
#              "Yes",# RT
#              "Yes",# CT
#              "G3+G4",
              "High") # PMM1

colUni_0 <- c("Female",
              "<=65y",
#              "lower risk",
#              "T1+T2",
#              "N0",
#              "Stage I+II",
              "T1+T2",
              "N0",
              "M0",
              "Stage I+II",
              "Negative", # negative safety margin
#              "no", # RT
#              "no", # CT
#              "G1+G2",
              "Low") # PMM1

#
# Define global functions
# contingency function, Tex's design (a confusion matrix with TP TN FP FN)



# #x for DEBUG only {
# # x osccCleanNA <- oscc
# # osccCleanNA has n=245 in LUAD (with margin free)
# # geneName
# #?? cutoff1 <- i <- round(nrow(oscc)/2) #213
# osccCleanNA_pos <- which(colnames(osccCleanNA) == "H.score_T")
# exp_geneName <- t(osccCleanNA[, osccCleanNA_pos ]) # z-score of ZZZ3
# cutoff1 <- quantile(exp_geneName)[3] # at 50%

# # contingency P-value
# contiT <- contingencyTCGA(osccCleanNA, geneName, cutoff1) # calling this function (OSCC cohort, PMM1, cutoff);
# # "margin" at column 8 of oscc
# #  oscc <- contiT[[1]] # updating PMM1_median
# chiT <- contiT[[2]] # extrac it from list by using [[]]; chiT$X2 is the P-value
# freq <- contiT[[3]] # well DONE
# }


contingencyTCGA <- function(osccCleanNA, geneName, cutoff1) {
  
  # create correlation table 1 with P-value by chisq.test
  library(reshape)
  library(data.table)
  #library(ca) # for Simple correspondence analysis
  
  #!!!!!# L <- 2 ("Gender"); R <- 8 ("margin") # in LUAD
  ## boundary of features column from first to last # 
  L <- which(colnames(osccCleanNA) == "Gender"); R <- which(colnames(osccCleanNA) == "margin")
  chiT <- data.frame(matrix(data = NA, nrow = R, ncol = 2)) # create a empty data.frame, 2D matrix as R*2
  #rown = 8
  freq <- data.frame(matrix(data = NA, nrow = 1, ncol = 4)) #, dimnames = list(c(1:rown+1), c("Var2", "L", "H"))))
  colnames(freq) <- c("Var2", 0, 1, "Features") # column 4 is for mapi
  #x freq <- array(NA, dim = c(5, 3, R)) # row, col, and R as 3D array

  #!!! PMM1 score## at col 13, 14
  osccCleanNA_pos <- which(colnames(osccCleanNA) == "H.score_T") #paste(geneName, sep="")) # position of PMM1 IHC score
  osccCleanNAM_pos <- which(colnames(osccCleanNA) == paste("PMM1", "_median", sep=""))
#  exp_geneName <- t(osccCleanNA[, osccCleanNA_pos ])
  

  # then [binomial of gene_median] resume the correlation tables
  osccCleanNA[osccCleanNAM_pos] <-(osccCleanNA[, osccCleanNA_pos ] >= cutoff1) +0 # binomial after osccCleanNA 
  # ***** addNA for counting all NA (e.g. there is 0, no 1) in "M" "stage" "margin"
  osccCleanNA$margin <- addNA(osccCleanNA$margin, ifany=F) # "always" add NA as a factor
  osccCleanNA$M <- addNA(osccCleanNA$M, ifany=F) # "always" add NA as a factor
  osccCleanNA$stage_2 <- addNA(osccCleanNA$stage_2, ifany=F) # "always" add NA as a factor
  
  
  # chisq.test(matrix(c(22, 5, 38, 21), ncol = 2), correct=F)$p.value # a example
  for (ii in L:R){ # generate freq table
    # Chi-square, Contingent table correlation, binary variables

    #  build a contingency table https://www.rdocumentation.org/packages/base/versions/3.4.3/topics/table
    # Powered by DataCamp, it might be running online
    # t, table from col of PMM1_median vs col ii (L "Gender" to R "margin"); deparse the argument (colnames) by deparse.level 2
#  the “pathological” case of two different kinds of NAs which are treated differently: exclude = if (useNA == "no") c(NA, NaN)
#  "unusual NA comes from addNA() as factor
    t<- NULL
    t <- table(osccCleanNA[,osccCleanNAM_pos], osccCleanNA[,ii], useNA = "ifany") #dnn=c(colnames(osccCleanNA[osccCleanNAM_pos])))  
               # useNA = "ifany"; "always" is for "margin" (with n=0 count)
    chiT[ii,1] <- colnames(osccCleanNA[ii]) # name list of feature variables from L to R
    chiT[ii,2] <- chisq.test(t)$p.value # retrieved in chiT$X2
    #print(ii) # debug
    
    obs <- as.data.frame(chisq.test(t)$observed)
    mdata <- melt(obs, id=c("Var2", "Var1")) # using the colnames of table t
    cdata <- as.data.frame(cast(mdata, Var2 ~ Var1)) 
    cdata <- cdata[, c(1:3)] # cdata should has 3 columns
    # names(dimnames(cdata)) = c("expression", "severity")
    # NA is necessary in column Var2
    if (!anyNA(cdata$Var2)) {
      cdata <- rbind(cdata, c(NA, 0, 0))
    }
    if (length(cdata$Var2) < 3) {
      cdata$ind <- seq_len(nrow(cdata)) # index as 1 2 3...
      crow <- data.frame("Var2"=2, "0"=0, "1"=0, "ind"=1.5)
      colnames(crow) <- colnames(cdata)
      cdata <- rbind(cdata, crow)
      cdata <- cdata[order(cdata$ind),]
      cdata <- cdata[,-4] # removal of index
    }
    cdata <- cbind(cdata, chiT[ii,1]) # add features name => # cdata has 4 columns
    # print(cdata)
    #freq[1+(ii-1)*rown,] <- cdata
    freq <- rbindlist(list(freq, cdata))
    #plot(ca(as.integer(cdata[-3,])))
  }
  freq <- freq[-1,] #removal of 1st row: NA
  name_freq <- colnames(osccCleanNA[L:R])
  name_freq <- t(name_freq) # "Gender" "age.at.diagnosis" "T"  "N"  "M"  "stage_2" and "margin"
  
  # array indexing https://cran.r-project.org/doc/manuals/r-release/R-intro.html#Array-indexing
  results <- list (osccCleanNA, chiT, freq) # it should to be returned/updated the osccCleanNA
  return(results)
} # end of contingencyTCGA function
#




# Define function introdcue Matthews Correlation Coefficient,
# \cite{Sanz-Pamplona2012} The Matthews Correlation Coefficient (MCC) [57] was
# chosen as measure of classification accuracy [58]. This index combines test
# sensitivity and specificity. It ranges from ???1 to 1 and its interpretation is
# similar to the Pearson???s correlation coefficient. In the context of a
# classification problem it is expected that MCC ranges 
# from 0 (no prediction ability at all) to +1 (perfect prediction) with negative values near zero
# possibly occurring in random classifiers due to sample variability. MCC values
# higher than 0.3 can be considered as indicative of high predictive value as they
# correspond to more than 65% accuracy in balanced data.

contingencyBin <- function (osccCleanNA, chiT, freq) { 
  # enough, geneName is not necessary a parameter
  # generate Table 2 by processing chiT and processing freq; place a "remark" and Matthews
  # for binary and integer convertion
  library(R.utils) # intToBin()
  library(compositions) # unbinary()
  # sigContig <- c(NA, "*") # remarkable when (p<0.05 || scoreContig == c(5,10))
  # for table 2
  contLowHighN <- c("Features",	"Low", "(%)",	"High",	"(%)", "Case no",	"P-value", "Remark", "Matthews")
  
  # rows need to be selected and reordered
  # *** mapping featuresUni <- osccCleanNA ; it needs to be corrected.
  indTbChi <- data.frame(cbind(featuresUni[-length(featuresUni)], c("Gender","age.at.diagnosis", "T", "N", "M", "stage_2","margin" )) )
#  indTbChi <- data.frame(cbind(featuresUni[-length(featuresUni)], c("Gender","ageDx", "pathologic_T", "pathologic_N", "pathologic_M", "stage","margin" )) )
  colnames(indTbChi) <- c("featuresUni", "osccCleanNA")
  
  #
  # rownames(tableChi1) <- featuresUni
  tableChi1 <- as.data.frame(setNames(replicate(length(contLowHighN), numeric(0), simplify = F), contLowHighN)) # declare a xlsx output data.frame (dynamic table)
  # colnames(tableChi1) <- contLowHighN
  for (i in 1:(length(featuresUni)-1)) { # without geneName_expression in last row
    # generate table by every two rows
    mapi <- as.character(freq$Features) == as.character(indTbChi$osccCleanNA[i]) # "gender" in [47:49] of freq
    mapi_pos <- which(mapi == T) #[47:49]
    var2U <- as.data.frame(freq[mapi_pos[1], ])
    var2L <- as.data.frame(freq[mapi_pos[2], ])
    subtotal <- sum(var2U[2:3]) + sum(var2L[2:3])
    U2 <- var2U[2];  U3 <- var2U[3]
    L2 <- var2L[2];  L3 <- var2L[3]
    pContig <- round(chiT[which(chiT$X1 == indTbChi$osccCleanNA[i]), 2], 4) # P-value retrieving from chiT$X2
    abin <- as.numeric(c(U3>U2, L3>L2, L2>U2, L3>U3)) # pattern of Sn=L3/(L2+L3), Sp=U2/(U2+U3)
    #    abin <- c(1,0,1,0) ; binary 1010 == decimal 10; while binary 0101 == decimal 5.
    achar <- paste(abin, collapse="")
    scoreContig <- unbinary(achar) # a integer score, 10 or 5 and more is significant, which will be marked as "*" in remark
    # Matthews correlation coefficient: 0-1 perfect prediction
    matthews <- (U2*L3-U3*L2)/sqrt((L3+U3)*(L3+L2)*(U2+U3)*(U2+L2))
    row0T2 <- data.frame(t(c(colUni_0[i], as.numeric(c(U2, round(U2/subtotal*100, 1), U3, round(U3/subtotal*100, 1), subtotal, pContig, (as.numeric((!(pContig>0.05) || scoreContig == c(5,10,1,2,4,7,8,11,13,14)))+0), matthews ))))) # remark 0 or 1
    row1T2 <- data.frame(t(c(colUni_1[i], as.numeric(c(L2, round(L2/subtotal*100, 1), L3, round(L3/subtotal*100, 1), NA, NA, NA, NA)))))
    #  colnames(row0) <- contLowHighN
    #  row0[,2:ncol(tableChi1)] <- as.numeric(as.character(row0[,2:ncol(tableChi1)]))
    tableChi1 <- rbind(tableChi1, row0T2, row1T2)
  }
  # tableChi1 <- cbind(tableChi1, chiT[-(1:5), 2]) # can not direct paste the p-value vector
  colnames(tableChi1) <- contLowHighN
  
  return(tableChi1)
} # end of contingencyBin define +++++++++++++++++++++++++++++ 
# ++++++++= end of Function defined ++++++++++## [reset for Repeat] end +++++++++++++++++++

#
#
#
#
#
#
#                   [run to here]
#
#
#
#
#
#
#
# [Data cleaning or osccClean], with margin feature ####
# clinicopathologic table [2018/03/09]
osccClean <- 0 # reset it
#osccClean <- oscc[!is.na(oscc$RFS..months._from.op), ] # optional ! # (removal of NA of X_RFS) => n=372 in OS and RFS cohort ??
osccClean <- oscc

# commonFeatures <- c(4:9) # common features: gender, age, margin and TNM :-)
commonFeatures <- which(colnames(oscc) %in% c("Gender", "age.at.diagnosis", "T", "N", "M", "margin")) #  6 essential features
osccCleanNA <- osccClean[complete.cases(osccClean[, commonFeatures]), ] # copy all features but remove Na or NaN entities in 6 essential features for their "completeness"
oscc_n256 <- osccCleanNA # n=256; removal of NA cases
osccClean1 <- osccClean # original cohort, with positive margin cohort, n=256
#+++ end of data cleaning ++++
#
# {
# surgical margin status: keeping 0 and excluding + margin (as 1; n=11)
osccCleanNA_freeMargin <- osccCleanNA[osccCleanNA$margin == 0, ] # margin==0
osccCleanNA <- osccCleanNA_freeMargin # n=245, LUAD s/p OP with margin free
# n=11, margin involved; 11/256 = 5.3% (how about it's survival impact on each individual genes in LUAD?)
# }
# 
# 
# margin free cohort (n=245): #### 
# osccCleanNA_freeMargin
# 
# margin positive cohort (n=256) #### 
# osccClean1



## 5b. (repeat) Cutoff finding [osccCleanNA] ####
oscc <- oscc0 <- osccCleanNA

# checking the completeness
which(complete.cases(oscc$H.score_T)==F) # no NaN -> 0
which(complete.cases(oscc$OS_IND)==F) # no Nan -> 0

# column 9 should be OS_IND
which(complete.cases(oscc[oscc$OS_IND==1,9])==F) #OS_IND ==1, death event (dead) => no NaN



# (skipped, if LUAD RFS is copied from OS)
which(complete.cases(oscc$RFS_IND)==F) # RFS_IND
# which(complete.cases(oscc$RFS..months._from.op)==F) #n=103; it may be imputed from OS time (oscc$OS..months._from.biopsy)
osccCleanNA_RFS <- oscc[which(complete.cases(oscc$RFS..months._from.op)==F), ]
osccCleanNA_RFS$RFS..months._from.op <- osccCleanNA_RFS[osccCleanNA_RFS$RFS_IND==1, 9] # imputed from OS..months._from.biopsy





# start 100 round
# oscc <- cbind(oscc0, osccCleanNA$H.score_T) # expression(IHC score) of PMM1
#exp_geneName519 <- exp_geneName 
find_repeat <- 100 # searching 100 slices in the interval
oscc0_pos <- which(colnames(oscc0) == "H.score_T") # oscc0$H.score_T as colunm 13
exp_geneName <- t(oscc0[, oscc0_pos]) # it needs horizontal matrix to store expression level
p_OS <- 0
cases_OS <- 0
p_RFS <- 0
cases_RFS <- 0
j <- 0
cut_featuresRemark <- 0 # finder results
# oscc_pos <- which(colnames(oscc) == geneName)
# oscc[oscc_pos] <- exp_geneName

cutoff <- quantile(exp_geneName, c(0.30,0.70)) # 30 percentile, 70 percentile

if (!require(pkg)){ 
  install.packages(pkg) 
} # Install package automatically if not there

# https://cran.r-project.org/web/packages/survival/survival.pdf
for (i in seq(cutoff[1], cutoff[2], length.out = find_repeat)){
  oscc[osccM_pos] <- (oscc0[oscc0_pos] >= i) +0 # oscc$PMM1_median <- (oscc0$PMM1 >= i) +0
  j <- j + 1
  # Surv {survival} run survival analysis: OS 
  # 1) event: For right censored data, the status indicator (event, OS_IND), normally 0=alive, 1=dead. 
  #   Other choices are TRUE/FALSE (TRUE =death) or 1/2 (2=death). 
  # 2) (For interval censored data, the status indicator is
  # 0=right censored, 1=event at time, 2=left censored, 3=interval censored.)
  # 3) Although unusual, the event indicator can be omitted, in which case all
  # subjects are assumed to have an event.
  # Usage: Surv(time, time2, event,
  # type=c('right', 'left', 'interval', 'counting', 'interval2', 'mstate'),
  # origin=0)
  mysurv <- Surv(oscc$OS..months._from.biopsy, oscc$OS_IND==1) #1==death event
  # Test for difference (log-rank test)
  surv_OS <- survdiff(mysurv ~ as.vector(unlist(oscc[osccM_pos]), mode="numeric"), data=oscc) # grouping by PMM1_median
  p_OS[j] <- format(pchisq(surv_OS$chisq, length(surv_OS$n)-1, lower.tail = FALSE), digits=3)
  cases_OS[j] <- surv_OS$n[1]
  #[coxph] - fits a Cox proportional hazards regression model
  #OS.km <- survfit(mysurv ~ oscc[osccM_pos], data=oscc, conf.type = "log-log")
  #  jpeg(file=paste("KMplot_OS_", geneName, i, ".jpg", sep = ""))
  #  plot(OS.km, lty=1, col=c("blue","red"), sub=paste("p-value =", p_OS[j]), main=paste("OS in OSCC(n=", surv_OS$n[1]+surv_OS$n[2],")/", geneName, "cutoff=", i), ylab="Percent Survival", xlab="Days")
  #  legend("topright", legend=c(paste("low(",surv_OS$n[1], ")"), paste("high(",surv_OS$n[2], ")")), lty=1:2, col=c("blue","red"))
  #  dev.off()
    
  
  # run survival analysis: RFS
  mysurv <- Surv(oscc$RFS..months._from.op, oscc$RFS_IND==1) #1==tumour recurrency
  # Test for difference (log-rank test)
  surv_RFS <- survdiff(mysurv ~ as.vector(unlist(oscc[osccM_pos]), mode="numeric"), data=oscc)
  p_RFS[j] <- format(pchisq(surv_RFS$chisq, length(surv_RFS$n)-1, lower.tail = FALSE), digits=3)
  cases_RFS[j] <- surv_RFS$n[1]
  #RFS.km <- survfit(mysurv ~ oscc[osccM_pos], data=oscc, conf.type = "log-log")
  #  jpeg(file=paste("KMplot_RFS_", geneName, i, ".jpg", sep = ""))
  #  plot(RFS.km, lty=1, col=c("blue","red"), sub=paste("p-value =", p_RFS[j]), main=paste("RFS in OSCC(n=", surv_RFS$n[1]+surv_RFS$n[2],")/", geneName, "cutoff=", i), ylab="Percent Survival", xlab="Days")
  #  legend("topright", legend=c(paste("low(",surv_RFS$n[1], ")"), paste("high(",surv_RFS$n[2], ")")), lty=1:2, col=c("blue","red"))
  #  dev.off()
  
  
  
  #  ++++ Beside KM P-value, contingency P-value is also important ! ++++++
  # Calling for contingency table to find significant clinicopathological features
  # contingency P-value
  contiT <- contingencyTCGA(oscc, geneName, i) # calling this function (OSCC cohort, PMM1); i is cutoff
  # "margin" at column 8 of oscc
    #  oscc <- contiT[[1]] # updating PMM1_median
  chiT <- contiT[[2]] # extrac it from list by using [[]]; chiT$X2 is the P-value
  freq <- contiT[[3]] # well DONE
  
  # then go through Table 2 construction by 
  # Calling contignecyBin
  tableChi1 <- contingencyBin (oscc, chiT, freq) # calculating the P-value
  # "margin" at column 8
  # to save "*" of column "remark" in 300 tableChi1 ###
  cut_featuresRemark <- cbind (cut_featuresRemark, tableChi1[7], tableChi1$Remark) # retrieve the P-value and remark
  #
} # end of 100 cut/slice, defined by [find_repeat]

# rownames(cut_featuresRemark) <- tableChi1$Features # there is duplicated row names :-)
# summary(cut_featuresRemark) to see result

# to generate the 201 first odd numbers; c(T,F) or seq.int(1, by = 2, length.out = find_repeat + 1)
# Sequence Generation: seq() or seq.int()
# https://www.rdocumentation.org/packages/base/versions/3.4.1/topics/seq

#describe(cut_featuresRemark[, c(1:(find_repeat * 2 + 1))[c(T,F)]])

###
library(graphics)
# cumulative P value curves
library(ggplot2) # http://www.cookbook-r.com/Graphs/Legends_(ggplot2)/
OS <- data.frame(cases_OS, as.numeric(p_OS))
# removal of duplicated item
OS <- OS[!duplicated(OS),]
colnames(OS)[2] <- "p_OS"
ggplot(OS, aes(x=cases_OS, y=p_OS)) + geom_point(size=2) +
  scale_x_discrete(limit = seq(0,500,50), labels = paste(seq(0,500,50), sep=",")) +
  scale_y_discrete(breaks = c(0, 0.05, 0.10, 0.50, 0.99)) + #, labels = c("0", "0.05", "0.10", "0.50", "0.99")) +
  coord_trans(y = "log10") +
  ggtitle(paste("Cumulative P value curves for OS under", geneName, "expression")) +
  xlab("# of patients") + ylab("P-value(log10)") +
  geom_hline(yintercept=0.05, linetype="dashed", color = "red", size=1, show.legend = T)
#  geom_vline(xintercept=300, color = "green", size=2, show.legend = T)
View(OS[OS$p_OS <= 0.05,1:2])
#
#plot(cases_OS, p_OS, type="p", log="y", main=paste("Cumulative P value curves for", geneName), xlab="# of patients", ylab="P-value")
#axis(side=4, at=c(0.01, 0.05, 0.2, 0.5, 0.99))
#abline(h=0.05, lty=2, col="red")
#abline(v=320, col="green")
#lines(x=c(cases_OS[1], cases_OS[length(cases_OS)]), y=c(0.05, 0.05), col = "red")

# plot(cases_RFS, p_RFS, type="p")
library(graphics)
# cumulative P value curves
library(ggplot2)
RFS <- data.frame(cases_RFS, as.numeric(p_RFS))
# removal of duplicated item
RFS <- RFS[!duplicated(RFS),]
colnames(RFS)[2] <- "p_RFS"
ggplot(RFS, aes(x=cases_RFS, y=p_RFS)) + geom_point(size=2) +
  scale_x_discrete(limit = seq(0,500,50), labels = paste(seq(0,500,50), sep=",")) +
  scale_y_discrete(breaks = c(0, 0.05, 0.10, 0.50, 0.99)) + #, labels = c("0", "0.05", "0.10", "0.50", "0.99")) +
  coord_trans(y = "log10") +
  ggtitle(paste("Cumulative P value curves for RFS under", geneName, "expression")) +
  xlab("# of patients") + ylab("P-value(log10)") +
  geom_hline(yintercept=0.05, linetype="dashed", color = "red", size=1, show.legend = T)
geom_vline(xintercept=150, color = "green", size=2, show.legend = T)
View(RFS[RFS$p_RFS <= 0.05,1:2])

# describe(cut_featuresRemark$Remark)

# [done cutoff finder] ++++++++++++++++++++++++

#
#
#
#
#
#
#
#{ PAUSE here.....
#
#
#
# {5c.Human}: decision cutoff=> please choice ONE from OS or RFS ####
case50_n <- 245 * 0.5
cutoff1 <- quantile(exp_geneName, c(case50_n/(surv_OS$n[1]+surv_OS$n[2])))  # 56 out of missingless LUAD cases
# or 
cutoff1 <- quantile(exp_geneName, c(case50_n/(surv_RFS$n[1]+surv_RFS$n[2])))  # 56 out of missingless LUAD cases

# for PMM1; cutoff at case n=56 (OS) or case n=? (RFS),  at lowerest P-value
# Kaplan-Meier survival of OS, P-value = 0.0198; of RFS, P-value = ? in LUAD



# 5d.[Statistic procedure] [Option-Cmd + E] for running to the end ####
# The features from column L to R; running to the end of R2Excel export +++++++++++
# The correlation of gene expression and clinical features  ### contingency tables
if (!require(pkg)){ 
  install.packages(pkg) 
} # Install package automatically if not there



# # call function for chiT, freq data.frame
# # x debug {contingencyTCGA)
# contiT <- contingencyTCGA(osccCleanNA, geneName, cutoff1) # calling this function (OSCC cohort, PMM1, cutoff1)
# osccCleanNA <- contiT[[1]] # updated
# chiT <- contiT[[2]] # extrac it from list by using [[]]
# freq <- contiT[[3]] # well DONE
# View(chiT) # for p-value (=0.0389 in pathologic_T)
# View(freq) # contingency table
# # }
#### ++++ Table 2 calculated well ++++++++++++++++




## KM survival curves in a specific cohort ####
# right censored? yes
# Kaplan-Meier curve + Log-rank test, + Cox proportional regression
# textbook: David Kleinbaum???Survival analysis: A self-learning text
library(survival)
library(rms)

osccCleanNA1 <- osccCleanNA
# osccCleanNA_pos <- which(colnames(osccCleanNA) == paste(geneName, sep=""))
osccCleanNAM_pos <- which(colnames(osccCleanNA) == paste("PMM1", "_median", sep=""))
# 3650 days = 10 years; cut it into 3650 days (do not need to remove 6 cases)
# over3650 <- osccCleanNA[osccCleanNA$X_RFS > 3650, 1]
# [1] TCGA-CV-5432-01 TCGA-CV-7183-01 TCGA-CV-7435-01 TCGA-CV-A45Q-01 TCGA-CV-A45R-01
# [6] TCGA-CV-A45T-01
# TCGA-CV-5430-01; removal one more
#osccCleanNA <- osccCleanNA[!osccCleanNA$sampleID %in% over3650,] # n=366
#osccCleanNA <- osccCleanNA[!osccCleanNA$sampleID == "TCGA-CV-5430-01",] # n=356

#OS
# cancer type shouble be defined at TCGA_cohort
mysurv <- Surv(osccCleanNA$OS..months._from.biopsy, osccCleanNA$OS_IND==1) #1==dead
# Test for difference (log-rank test)
surv_OS1 <- survdiff(mysurv ~ as.vector(osccCleanNA[, osccCleanNAM_pos], mode="numeric"), data=osccCleanNA) # don't need "unlist"
p_OS1 <- format(pchisq(surv_OS1$chisq, length(surv_OS1$n)-1, lower.tail = FALSE), digits=3)
cases_OS1 <- surv_OS1$n[1]
#[coxph] - fits a Cox proportional hazards regression model
# or [survfit] - Kaplan-Meier curve
# confidence intervals as log hazard or log(-log(survival))
OS.km <- survfit(mysurv ~ as.vector(unlist(osccCleanNA[, osccCleanNAM_pos]), mode="numeric"), data=osccCleanNA, type= "kaplan-meier", conf.type = "log-log")
# 365.25 days in a year for xscale => 3650 days for 10 years
# 12 months per year for 5 years => 60 months

plot(OS.km, lty=1, xscale=12, xmax=60, col=c("blue","red"), sub=paste("P-value =", p_OS1), main=paste("OS in TCGA", TCGA_cohort, "(n=", surv_OS1$n[1]+surv_OS1$n[2],")/", geneName), ylab="Percent Survival", xlab="Years")
legend("topright", legend=c(paste("low(",surv_OS1$n[1], ")"), paste("high(",surv_OS1$n[2], ")")), lty=1:2, col=c("blue","red"))
#plot(OS.km, lty=1, col=c("blue","red"), sub="p= 0.816", main="OS in osccCleanNA(n=505)/gene level", ylab="Percent Survival", xlab="Days")
#legend("topright", legend=c('Low', 'High'), lty=1:2, col=c("blue","red"))
summary(OS.km, times = seq(0, 3000, 100))




#RFS
mysurv <- Surv(osccCleanNA$RFS..months._from.op, osccCleanNA$RFS_IND==1) #1==tumour recurrency
# Test for difference (log-rank test)
surv_RFS1 <- survdiff(mysurv ~ as.vector(unlist(osccCleanNA[, osccCleanNAM_pos]), mode="numeric"), data=osccCleanNA)
p_RFS1 <- format(pchisq(surv_RFS1$chisq, length(surv_RFS1$n)-1, lower.tail = FALSE), digits=3)
cases_RFS1 <- surv_RFS1$n[1]
#[coxph] - fits a Cox proportional hazards regression model
#https://rstudio-pubs-static.s3.amazonaws.com/5588_72eb65bfbe0a4cb7b655d2eee0751584.html
RFS.km <- survfit(mysurv ~ as.vector(unlist(osccCleanNA[, osccCleanNAM_pos]), mode="numeric"), data=osccCleanNA, conf.type = "log-log")
# 12 months per year for 5 years => 60 months

plot(RFS.km, lty=1, xscale=12, xmax=60, col=c("blue","red"), sub=paste("p-value =", p_RFS1), main=paste("RFS in TCGA", TCGA_cohort, "(n=", surv_RFS1$n[1]+surv_RFS1$n[2],")/", geneName), ylab="Percent Survival", xlab="Years")
legend("topright", legend=c(paste("low(",surv_RFS1$n[1], ")"), paste("high(",surv_RFS1$n[2], ")")), lty=1:2, col=c("blue","red"))
#plot(RFS.km, lty=1, col=c("blue","red"), sub="p= 0.0487", main="RFS in osccCleanNA(n=372)/gene level", ylab="Percent Survival", xlab="Days")
#legend("topright", legend=c('Low', 'High'), lty=1:2, col=c("blue","red"))
summary(RFS.km, times = seq(0, 3000, 100))
# [DONE] ##  back to Cutoff finder (n=185) or set 50:50 on section [5c]






## x6a.to connect tableau with R studio for visualization and export ####
#https://community.tableau.com/thread/236068 a TCP/IP server that allows other
#programs to use facilities of R 
#
# http://www.rforge.net/Rserve/doc.html 
# in High Sierra: there is no :-) => R must have been configured with --enable-R-shlib in order to use Rserve
#If the compilation fails, please check that R shared library exists and is properly
#installed. It is located in $RHOME/bin and is named libR.so or libR.dylib. If
#it is missing, get an R distribution with shared library included or compile R
#with --enable-R-shlib configure flag.
# 
# the location of the R home directory, which is the root of the installed R tree
# >R.home()
# [1] "/usr/local/Cellar/r/3.4.3/lib/R" # my RHOME
# at /usr/local/Cellar/r/3.4.3/lib/R/lib/libR.dylib
# configure the path to the R library by:
# $ export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/Cellar/r/3.4.3/lib/R/bin
# $ export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/Cellar/r/3.4.3/lib/R/lib
#x $ export R_HOME=$R_HOME:/usr/local/Cellar/r/3.4.3/lib/R
# git clone --recursive https://github.com/s-u/Rserve.git Rserve
# compiling from the source
install.packages("/Users/apple/Downloads/Rserve_1.8-5.tar.gz", repos = NULL, type="source")
install.packages("Rserve", repos="http://rforge.net/", type="source")
install.packages("Rserve", repos="http://rforge.net/")
# or >R CMD INSTALL Rserve_1.8-5.tar.gz

library(Rserve)
Rserve(port=6311, args = "--no-save")
# Rserv started in daemon mode
# $R CMD Rserve --help
# 
# http://www.snaq.net/software/rjava-macos.php
# install.packages("JRI", repos="https://www.rforge.net/", type="source") # JRI is now part of rJava
# > system.file("jri",package="rJava")
# [1] "/usr/local/lib/R/3.4/site-library/rJava/jri"
# 
# http://www.clearlyandsimply.com/clearly_and_simply/2016/06/writing-and-reading-tableau-views-to-and-from-databases-and-text-files-part-1.html
install.packages("/Users/apple/Downloads/RODBC_1.3-15.tar.gz", repos = NULL, type="source")
library("RODBC") #configure: error: "ODBC headers sql.h and sqlext.h not found"


# Open tableau
# Help->Settings and performance-> Manage external server connection
# Server : localhost # or "127.0.0.1"
# Port: 6311
# then test the connection ...ok

# Start using the R scripts in Tableau
# SCRIPT_BOOL , SCRIPT_INT , SCRIPT_REAL , SCRIPT_STR

# http://www.simafore.com/blog/bid/120209/Integrating-Tableau-and-R-for-data-analytics-in-four-simple-steps
# flowchart Use the calculated Correlation as another measure to create a chart
# and dashboard
# https://www.tableau.com/learn/whitepapers/tableau-and-r-faq?__src=liftigniter&__widget=blog-widget&li_source=LI&li_medium=blog-widget
# e.x. “k-mean clustering." This could easily be done in a calculated field with
# a single line of R code:
# SCRIPT_INT(‘kmeans(data.frame(.arg1,.arg2,.arg3,.arg4),3)$cluster;',
# SUM([Petal length]), SUM([Petal width]),SUM([Sepal length]),SUM([Sepal
# width]))
            
            



##4$$##
# [for updating required packages if R upgraded]

# install.packages("r2excel") # dependence: ##
# Loading required package: xlsx
# Loading required package: xlsxjars


# Loading required package: rJava [loading with error]
# rJava installation => it DOESN'T need to load library("rJava")?

# x guide: https://github.com/MTFA/CohortEx/wiki/Run-rJava-with-RStudio-under-OSX-10.10,-10.11-(El-Capitan)-or-10.12-(Sierra)
# A common error when trying to load rJava in RStudio due to the known limitation with RStudio on OS X.
# JDK 8 for Mac => http://download.oracle.com/otn-pub/java/jdk/8u144-b01/090f390dda5b47b9b721c7dfaa008135/jdk-8u144-macosx-x64.dmg




# 2017 new guide for rJava: https://github.com/s-u/rJava/issues/119 #
# JNI issue resolved ; Getting rJava working on macOS
# a good reference for JAVA_HOME finding http://www.snaq.net/software/rjava-macos.php
# >/usr/libexec/java_home -V
# ==> it shows ""Java SE 9.0.1"	/Library/Java/JavaVirtualMachines/jdk-9.0.1.jdk/Contents/Home"
# >java -version ==> java version "9.0.1"
# >which java
# >ls -l /usr/bin/java
# ==> it shows "/System/Library/Frameworks/JavaVM.framework/Versions/Current/Commands/java"
# # >ls -l /usr/libexec/java_home 
# ==> /System/Library/Frameworks/JavaVM.framework/Versions/Current/Commands/java_home
#
# > export JAVA_HOME="$(/usr/libexec/java_home)"
# JAVA_HOME should at /Library/Java/JavaVirtualMachines/jdk-9.0.1.jdk/Contents/Home

# >export LIBJVM=$(find "${JAVA_HOME}" -name 'libjvm.dylib') # Java library path: $(JAVA_HOME)/lib/server
# ==> # /Library/Java/JavaVirtualMachines/jdk-9.0.1.jdk/Contents/Home/lib/server/libjvm.dylib
# # >echo $JAVA_HOME
 
#?? >exec /System/Library/Frameworks/JavaVM.framework/Versions/1.7/Commands/java "$@"

# https://www.ntu.edu.sg/home/ehchua/programming/howto/Environment_Variables.html
# >printenv
# JAVA_HOME and JRE_HOME: maintain the locations of JDK and JRE installed directory, respectively
#x?? in "/etc/profile" # to set JAVA_HOME
#xx setenv JAVA_HOME /System/Library/Frameworks/JavaVM.framework/Home  #the older csh (C-shell) command
# >"set varname=value" # local variable
# >"unset varname"
# # Bash shell uses colon (:) as the path separator
#x >export JAVA_HOME="/usr/libexec/java_home -v"
# >export JAVA_HOME="/System/Library/Frameworks/JavaVM.framework/Versions/Current/Commands/java_home" # as the global environment
# >source /etc/profile # to refresh it
# >echo $JAVA_HOME
# # >R CMD javareconf

# # # % [2017/12/17 not any more] Resolution with this R script (CohortEx R script) =>
# 
# if (Sys.info()['sysname'] == 'Darwin') {
#   libjvm <- paste0(system2('/usr/libexec/java_home', stdout = TRUE)[1], '/lib/server/libjvm.dylib')
#   message (paste0('Load libjvm.dylib from: ', libjvm))
#   dyn.load(libjvm)
# }

# install.packages("rJava") # not available in R 3.4.2
# # or
# >wget https://www.rforge.net/rJava/snapshot/rJava_0.9-9.tar.gz
# >R CMD INSTALL rJava_0.9-9.tar.gz
# # [install rJava from source]
install.packages("rJava", repos="http://rforge.net", type="source") # checking JNI data types... ok
# ** testing if installed package can be loaded
# * DONE (rJava)
# >sudo Rscript -e 'install.packages("rJava", repos="http://rforge.net", type="source")'


library(rJava) # PASSED

library("xlsx")
# solved ERROR: compilation failed for package ‘git2r’; ERROR: dependencies ‘httr’, ‘git2r’ are not available for package ‘devtools’
# brew install openssl
library("openssl")
install.packages('git2r', repos="https://github.com/ropensci/git2r", type='source')
# or clone then build (make) in BASH
# $ git clone https://github.com/ropensci/git2r.git
# $ cd git2r
# $ make install
# it response as * DONE (git2r) for R version 3.4.3)

#xx no .c file $ git clone https://github.com/r-lib/httr
install.packages("~/Downloads/httr_1.3.1.tar.gz", repos = NULL, type = "source")

install.packages("devtools")
library(devtools)
install_github("kassambara/r2excel") # the only way OK
# xx install.packages("r2excel", repos="https://api.github.com/repos/kassambara/r2excel/zipball/master", type="source") # 
library(r2excel) # (package ‘r2excel’ is available for R version 3.4.1 below)



## [Export] Begin
# *** [cont] all features needs to be adjusted one by one; updated on [2018/03/09] for LUAD osccCleanNA
## (7.)[coxph]  Export part I, univariate survival: ####
## => coxph, for gene high, low and other features #
## multivariant survival: coxph, Cox proportional Hazard regression Model
## covariates (group, sex...)
# http://www.uni-kiel.de/psychologie/rexrepos/posts/survivalCoxPH.html

#[coxph] - fits a Cox proportional hazards regression model
#<<vdxfitcox>> to fit a Cox regression model with the 7 genes.
#Multivariate Cox proportional hazards regression analyses on recurrence-free survival time in gene expression datasets.
#The status indicator, normally 0=alive, 1=dead. Other choices are TRUE/FALSE (TRUE = death) or 1/2 (2=death). For interval censored data, the status indicator is 0=right censored, 1=event at time, 2=left censored, 3=interval censored. A
# time (months) is OVERALL.SURVIVAL OS or RECURRENCE.FREE.SURVIVAL RFS
# status = overall.survival.indicator or recurrence.free.survival.indicator
library(survival)
# [Multivariate for OS] Table 3
osHRmulti <- 0
# selectedFeatures <- colnames(osccCleanNA[commonFeatures])
# colname of oscc = c("Unique.ID","Gender","age.at.diagnosis",
# "margin","T","N",
# "M","stage_2","OS..months._from.biopsy",
# "OS_IND","RFS..months._from.op","RFS_IND","H.score_T")

# rename all of colnames (features) to osccT
colnames(osccCleanNA) <- c("Unique.ID","Gender","ageDx",
                                      "pathologic_T","pathologic_N",
                                      "pathologic_M","stage", "margin",
                                      "OS_IND","OS..months._from.biopsy",
                                      "RFS_IND", "RFS..months._from.op", "H.score_T", paste("PMM1", "_median", sep=""))
# 
oscox <- coxph(Surv(osccCleanNA$OS..months._from.biopsy, osccCleanNA$OS_IND==1) ~ +
  Gender +
  ageDx +
#                 primary_site +
#                 clinical_T +
#                 clinical_N +
#                 clinical_stage +
  pathologic_T +
  pathologic_N +
  pathologic_M +
  stage +
  margin +
#  R.T +
#  C.T +
#                 presence_of_pathological_nodal_extracapsular_spread +
#                 neoplasm_histologic_grade +
                 as.vector(osccCleanNA[, osccCleanNAM_pos], mode="numeric"), data=osccCleanNA) # PMM1_median: IHC score (H.score_T) is [low, high]

summary(oscox)
osHRmulti <- cbind(data.frame(summary(oscox)$conf.int)[-2], data.frame(summary(oscox)$coefficients)[5])
rownames(osHRmulti)[nrow(osHRmulti)] <- paste("PMM1", "_median", sep="")
osHRmulti <- round(osHRmulti, 3) # rounding them by 3 decimal
View(osHRmulti)



# [Univariate for OS] copy and past one by one
osHR <- 0; oscox <- 0
# features 13 -> 14 (with margin status)
# Gender
oscox <- coxph(Surv(osccCleanNA$OS..months._from.biopsy, osccCleanNA$OS_IND==1) ~ +
                 Gender, data=osccCleanNA)
osHR <- rbind(osHR, cbind(data.frame(summary(oscox)$conf.int)[-2], data.frame(summary(oscox)$coefficients)[5]))
# ageDx
oscox <- coxph(Surv(osccCleanNA$OS..months._from.biopsy, osccCleanNA$OS_IND==1) ~ +
                 ageDx, data=osccCleanNA)
osHR <- rbind(osHR, cbind(data.frame(summary(oscox)$conf.int)[-2], data.frame(summary(oscox)$coefficients)[5]))

#oscox <- coxph(Surv(osccCleanNA$OS..months._from.biopsy, osccCleanNA$OS_IND==1) ~ +
#                 primary_site, data=osccCleanNA)
#osHR <- rbind(osHR, cbind(data.frame(summary(oscox)$conf.int)[-2], data.frame(summary(oscox)$coefficients)[5]))

# oscox <- coxph(Surv(osccCleanNA$OS..months._from.biopsy, osccCleanNA$OS_IND==1) ~ +
#                  clinical_T, data=osccCleanNA)
# osHR <- rbind(osHR, cbind(data.frame(summary(oscox)$conf.int)[-2], data.frame(summary(oscox)$coefficients)[5]))
# 
# oscox <- coxph(Surv(osccCleanNA$OS..months._from.biopsy, osccCleanNA$OS_IND==1) ~ +
#                  clinical_N, data=osccCleanNA)
# osHR <- rbind(osHR, cbind(data.frame(summary(oscox)$conf.int)[-2], data.frame(summary(oscox)$coefficients)[5]))
# 
# oscox <- coxph(Surv(osccCleanNA$OS..months._from.biopsy, osccCleanNA$OS_IND==1) ~ +
#                  clinical_stage, data=osccCleanNA)
# osHR <- rbind(osHR, cbind(data.frame(summary(oscox)$conf.int)[-2], data.frame(summary(oscox)$coefficients)[5]))

# pathologic_T
oscox <- coxph(Surv(osccCleanNA$OS..months._from.biopsy, osccCleanNA$OS_IND==1) ~ +
                 pathologic_T, data=osccCleanNA)
osHR <- rbind(osHR, cbind(data.frame(summary(oscox)$conf.int)[-2], data.frame(summary(oscox)$coefficients)[5]))
# pathologic_N
oscox <- coxph(Surv(osccCleanNA$OS..months._from.biopsy, osccCleanNA$OS_IND==1) ~ +
                 pathologic_N, data=osccCleanNA)
osHR <- rbind(osHR, cbind(data.frame(summary(oscox)$conf.int)[-2], data.frame(summary(oscox)$coefficients)[5]))
# pathologic_M
oscox <- coxph(Surv(osccCleanNA$OS..months._from.biopsy, osccCleanNA$OS_IND==1) ~ +
                 pathologic_M, data=osccCleanNA)
osHR <- rbind(osHR, cbind(data.frame(summary(oscox)$conf.int)[-2], data.frame(summary(oscox)$coefficients)[5]))
# stage
oscox <- coxph(Surv(osccCleanNA$OS..months._from.biopsy, osccCleanNA$OS_IND==1) ~ +
                 stage, data=osccCleanNA)
osHR <- rbind(osHR, cbind(data.frame(summary(oscox)$conf.int)[-2], data.frame(summary(oscox)$coefficients)[5]))

# margin
oscox <- coxph(Surv(osccCleanNA$OS..months._from.biopsy, osccCleanNA$OS_IND==1) ~ +
                 margin, data=osccCleanNA)
osHR <- rbind(osHR, cbind(data.frame(summary(oscox)$conf.int)[-2], data.frame(summary(oscox)$coefficients)[5]))


# oscox <- coxph(Surv(osccCleanNA$OS..months._from.biopsy, osccCleanNA$OS_IND==1) ~ +
#                  R.T, data=osccCleanNA)
# osHR <- rbind(osHR, cbind(data.frame(summary(oscox)$conf.int)[-2], data.frame(summary(oscox)$coefficients)[5]))
# 
# oscox <- coxph(Surv(osccCleanNA$OS..months._from.biopsy, osccCleanNA$OS_IND==1) ~ +
#                  C.T, data=osccCleanNA)
# osHR <- rbind(osHR, cbind(data.frame(summary(oscox)$conf.int)[-2], data.frame(summary(oscox)$coefficients)[5]))

# oscox <- coxph(Surv(osccCleanNA$OS..months._from.biopsy, osccCleanNA$OS_IND==1) ~ +
#                  presence_of_pathological_nodal_extracapsular_spread, data=osccCleanNA)
# osHR <- rbind(osHR, cbind(data.frame(summary(oscox)$conf.int)[-2], data.frame(summary(oscox)$coefficients)[5]))
# 
# oscox <- coxph(Surv(osccCleanNA$OS..months._from.biopsy, osccCleanNA$OS_IND==1) ~ +
#                  neoplasm_histologic_grade, data=osccCleanNA)
# osHR <- rbind(osHR, cbind(data.frame(summary(oscox)$conf.int)[-2], data.frame(summary(oscox)$coefficients)[5]))

# IHC score
oscox <- coxph(Surv(osccCleanNA$OS..months._from.biopsy, osccCleanNA$OS_IND==1) ~ +
                 as.vector(osccCleanNA[, osccCleanNAM_pos], mode="numeric"), data=osccCleanNA)
osHR <- rbind(osHR, cbind(data.frame(summary(oscox)$conf.int)[-2], data.frame(summary(oscox)$coefficients)[5]))

#
osHR <- round(osHR[-1,], 3) # skip first row 
View(osHR) #-> univariate table


#
# test for the assumption of proportionality of hazards
cox.zph(oscox)
plot(cox.zph(oscox)) # flat curve is good, not violate asumption
##


## [Multivariate for RFS]
rfsHRmulti <- 0
#Significant codes:  0 ???***??? 0.001 ???**??? 0.01 ???*??? 0.05 ???.??? 0.1 ??? ??? 1
rfscox <- coxph(Surv(osccCleanNA$RFS..months._from.op, osccCleanNA$RFS_IND==1) ~ +
                  Gender +
                  ageDx +
                  #                 primary_site +
                  #                 clinical_T +
                  #                 clinical_N +
                  #                 clinical_stage +
                  pathologic_T +
                  pathologic_N +
                  pathologic_M +
                  stage +
                  margin +
                  # R.T +
                  # C.T +
                  #                 presence_of_pathological_nodal_extracapsular_spread +
                  #                 neoplasm_histologic_grade +
                  as.vector(unlist(osccCleanNA[, osccCleanNAM_pos]), mode="numeric"), data=osccCleanNA)

summary(rfscox)
rfsHRmulti <- cbind(data.frame(summary(rfscox)$conf.int)[-2], data.frame(summary(rfscox)$coefficients)[5])
rownames(rfsHRmulti)[nrow(rfsHRmulti)] <- paste("PMM1", "_median", sep="")
rfsHRmulti <- round(rfsHRmulti, 3)
View(rfsHRmulti)

# test for the assumption of proportionality of hazards
cox.zph(rfscox)
plot(cox.zph(rfscox))
#
# [Univariate for RFS]
rfsHR <- 0; rfscox <- 0
#
# Gender
rfscox <- coxph(Surv(osccCleanNA$RFS..months._from.op, osccCleanNA$RFS_IND==1) ~ +
                  Gender, data=osccCleanNA)
rfsHR <- rbind(rfsHR, cbind(data.frame(summary(rfscox)$conf.int)[-2], data.frame(summary(rfscox)$coefficients)[5]))
# ageDx
rfscox <- coxph(Surv(osccCleanNA$RFS..months._from.op, osccCleanNA$RFS_IND==1) ~ +
                  ageDx, data=osccCleanNA)
rfsHR <- rbind(rfsHR, cbind(data.frame(summary(rfscox)$conf.int)[-2], data.frame(summary(rfscox)$coefficients)[5]))

# rfscox <- coxph(Surv(osccCleanNA$RFS..months._from.op, osccCleanNA$RFS_IND==1) ~ +
#                   primary_site, data=osccCleanNA)
# rfsHR <- rbind(rfsHR, cbind(data.frame(summary(rfscox)$conf.int)[-2], data.frame(summary(rfscox)$coefficients)[5]))

# rfscox <- coxph(Surv(osccCleanNA$RFS..months._from.op, osccCleanNA$RFS_IND==1) ~ +
#                   clinical_T, data=osccCleanNA)
# rfsHR <- rbind(rfsHR, cbind(data.frame(summary(rfscox)$conf.int)[-2], data.frame(summary(rfscox)$coefficients)[5]))
# 
# rfscox <- coxph(Surv(osccCleanNA$RFS..months._from.op, osccCleanNA$RFS_IND==1) ~ +
#                   clinical_N, data=osccCleanNA)
# rfsHR <- rbind(rfsHR, cbind(data.frame(summary(rfscox)$conf.int)[-2], data.frame(summary(rfscox)$coefficients)[5]))
# 
# rfscox <- coxph(Surv(osccCleanNA$RFS..months._from.op, osccCleanNA$RFS_IND==1) ~ +
#                   clinical_stage, data=osccCleanNA)
# rfsHR <- rbind(rfsHR, cbind(data.frame(summary(rfscox)$conf.int)[-2], data.frame(summary(rfscox)$coefficients)[5]))

# pathologic_T
rfscox <- coxph(Surv(osccCleanNA$RFS..months._from.op, osccCleanNA$RFS_IND==1) ~ +
                  pathologic_T, data=osccCleanNA)
rfsHR <- rbind(rfsHR, cbind(data.frame(summary(rfscox)$conf.int)[-2], data.frame(summary(rfscox)$coefficients)[5]))
# pathologic_N
rfscox <- coxph(Surv(osccCleanNA$RFS..months._from.op, osccCleanNA$RFS_IND==1) ~ +
                  pathologic_N, data=osccCleanNA)
rfsHR <- rbind(rfsHR, cbind(data.frame(summary(rfscox)$conf.int)[-2], data.frame(summary(rfscox)$coefficients)[5]))
# pathologic_M
rfscox <- coxph(Surv(osccCleanNA$RFS..months._from.op, osccCleanNA$RFS_IND==1) ~ +
                  pathologic_M, data=osccCleanNA)
rfsHR <- rbind(rfsHR, cbind(data.frame(summary(rfscox)$conf.int)[-2], data.frame(summary(rfscox)$coefficients)[5]))
# stage
rfscox <- coxph(Surv(osccCleanNA$RFS..months._from.op, osccCleanNA$RFS_IND==1) ~ +
                  stage, data=osccCleanNA)
rfsHR <- rbind(rfsHR, cbind(data.frame(summary(rfscox)$conf.int)[-2], data.frame(summary(rfscox)$coefficients)[5]))

# margin
rfscox <- coxph(Surv(osccCleanNA$RFS..months._from.op, osccCleanNA$RFS_IND==1) ~ +
                  margin, data=osccCleanNA)
rfsHR <- rbind(rfsHR, cbind(data.frame(summary(rfscox)$conf.int)[-2], data.frame(summary(rfscox)$coefficients)[5]))

# rfscox <- coxph(Surv(osccCleanNA$RFS..months._from.op, osccCleanNA$RFS_IND==1) ~ +
#                   R.T, data=osccCleanNA)
# rfsHR <- rbind(rfsHR, cbind(data.frame(summary(rfscox)$conf.int)[-2], data.frame(summary(rfscox)$coefficients)[5]))
# 
# rfscox <- coxph(Surv(osccCleanNA$RFS..months._from.op, osccCleanNA$RFS_IND==1) ~ +
#                   C.T, data=osccCleanNA)
# rfsHR <- rbind(rfsHR, cbind(data.frame(summary(rfscox)$conf.int)[-2], data.frame(summary(rfscox)$coefficients)[5]))

# rfscox <- coxph(Surv(osccCleanNA$RFS..months._from.op, osccCleanNA$RFS_IND==1) ~ +
#                   presence_of_pathological_nodal_extracapsular_spread, data=osccCleanNA)
# rfsHR <- rbind(rfsHR, cbind(data.frame(summary(rfscox)$conf.int)[-2], data.frame(summary(rfscox)$coefficients)[5]))
# 
# rfscox <- coxph(Surv(osccCleanNA$RFS..months._from.op, osccCleanNA$RFS_IND==1) ~ +
#                   neoplasm_histologic_grade, data=osccCleanNA)
# rfsHR <- rbind(rfsHR, cbind(data.frame(summary(rfscox)$conf.int)[-2], data.frame(summary(rfscox)$coefficients)[5]))

# IHC score
rfscox <- coxph(Surv(osccCleanNA$RFS..months._from.op, osccCleanNA$RFS_IND==1) ~ +
                  as.vector(osccCleanNA[, osccCleanNAM_pos], mode="numeric"), data=osccCleanNA)
rfsHR <- rbind(rfsHR, cbind(data.frame(summary(rfscox)$conf.int)[-2], data.frame(summary(rfscox)$coefficients)[5]))
#
rfsHR <- round(rfsHR[-1,], 3) # skip first row
View(rfsHR) # univariate table
## DONE and export part II ##




## R2Excel export ####
## # http://www.sthda.com/english/wiki/r2excel-read-write-and-format-easily-excel-files-using-r-software
library("xlsx")
library("r2excel")
# Create an Excel workbook. Both .xls and .xlsx file formats can be used.
path_LUAD <- "/Users/apple/Google\ Drive/2018_PhD_dissertation/LUAD_Peter_survival/"
# path_LUAD <- "/Users/apple/Documents/My\ Tableau\ Repository/Workbooks/"
setwd(path_LUAD) # change working directory to the google drive

filenamex <- paste(TCGA_cohort, "_survivalAnalysis_", geneName, ".xlsx", sep = "")
wb <- createWorkbook(type="xlsx")
# Create a sheet in that workbook
sheet <- xlsx::createSheet(wb, sheetName = paste(geneName, "_multivariate"))
# [add data row by row, start from column 2]
#+++++++++++++++++++++++++++++++
## Add paragraph : Author
library("tis") # by Brian Salzer
# today(), arg must be ti, tis, ts, tif, or tifName
author <- paste("Reported by Tex Li-Hsing Chi. \n",
                "tex@gate.sinica.edu.tw \n", Sys.Date(), sep="")
xlsx.addParagraph(wb, sheet, value=author, isItalic=TRUE, colSpan=5, 
                  rowSpan=4, fontColor="darkgray", fontSize=14)
xlsx.addLineBreak(sheet, 3)


# Part I: table 2 (calling function contingencyBin2)
contingencyBin2 <- function (osccCleanNA, chiT, freq) { 
  # enough, geneName is not necessary a parameter
  # generate Table 2 by processing chiT and processing freq; place a "remark" and Matthews
  # for binary and integer convertion
  library(R.utils) # intToBin()
  library(compositions) # unbinary()
  # sigContig <- c(NA, "*") # remarkable when (p<0.05 || scoreContig == c(5,10))
  # for table 2
  contLowHighN <- c("Features",	"Low", "(%)",	"High",	"(%)", "Case no",	"P-value", "Remark", "Matthews")
  
  # rows need to be selected and reordered
  # *** mapping featuresUni <- osccCleanNA ; it needs to be corrected.
  #  indTbChi <- data.frame(cbind(featuresUni[-length(featuresUni)], c("Gender","age.at.diagnosis", "T", "N", "M", "stage_2","margin" )) )
  indTbChi <- data.frame(cbind(featuresUni[-length(featuresUni)], c("Gender","ageDx", "pathologic_T", "pathologic_N", "pathologic_M", "stage","margin" )) )
  colnames(indTbChi) <- c("featuresUni", "osccCleanNA")
  
  #
  # rownames(tableChi1) <- featuresUni
  tableChi1 <- as.data.frame(setNames(replicate(length(contLowHighN), numeric(0), simplify = F), contLowHighN)) # declare a xlsx output data.frame (dynamic table)
  # colnames(tableChi1) <- contLowHighN
  for (i in 1:(length(featuresUni)-1)) { # without geneName_expression in last row
    # generate table by every two rows
    mapi <- as.character(freq$Features) == as.character(indTbChi$osccCleanNA[i]) # "gender" in [47:49] of freq
    mapi_pos <- which(mapi == T) #[47:49]
    var2U <- as.data.frame(freq[mapi_pos[1], ])
    var2L <- as.data.frame(freq[mapi_pos[2], ])
    subtotal <- sum(var2U[2:3]) + sum(var2L[2:3])
    U2 <- var2U[2];  U3 <- var2U[3]
    L2 <- var2L[2];  L3 <- var2L[3]
    pContig <- round(chiT[which(chiT$X1 == indTbChi$osccCleanNA[i]), 2], 4) # P-value retrieving from chiT$X2
    abin <- as.numeric(c(U3>U2, L3>L2, L2>U2, L3>U3)) # pattern of Sn=L3/(L2+L3), Sp=U2/(U2+U3)
    #    abin <- c(1,0,1,0) ; binary 1010 == decimal 10; while binary 0101 == decimal 5.
    achar <- paste(abin, collapse="")
    scoreContig <- unbinary(achar) # a integer score, 10 or 5 and more is significant, which will be marked as "*" in remark
    # Matthews correlation coefficient: 0-1 perfect prediction
    matthews <- (U2*L3-U3*L2)/sqrt((L3+U3)*(L3+L2)*(U2+U3)*(U2+L2))
    row0T2 <- data.frame(t(c(colUni_0[i], as.numeric(c(U2, round(U2/subtotal*100, 1), U3, round(U3/subtotal*100, 1), subtotal, pContig, (as.numeric((!(pContig>0.05) || scoreContig == c(5,10,1,2,4,7,8,11,13,14)))+0), matthews ))))) # remark 0 or 1
    row1T2 <- data.frame(t(c(colUni_1[i], as.numeric(c(L2, round(L2/subtotal*100, 1), L3, round(L3/subtotal*100, 1), NA, NA, NA, NA)))))
    #  colnames(row0) <- contLowHighN
    #  row0[,2:ncol(tableChi1)] <- as.numeric(as.character(row0[,2:ncol(tableChi1)]))
    tableChi1 <- rbind(tableChi1, row0T2, row1T2)
  }
  # tableChi1 <- cbind(tableChi1, chiT[-(1:5), 2]) # can not direct paste the p-value vector
  colnames(tableChi1) <- contLowHighN
  
  return(tableChi1)
} # end of contingencyBin2 define +++++++++++++++++++++++++++++ 
tableChi1 <- contingencyBin2 (osccCleanNA, chiT, freq)

# Export Table 2
xlsx.addHeader(wb, sheet, value=paste("Table 2. The correlation of", geneName,  "expression and clinical features."),
               level=5, color="black", underline=0)
# add TCGA analysis tables
xlsx.addLineBreak(sheet, 1)
xlsx.addTable(wb, sheet, data = t(data.frame(c(paste(geneName, "expression"), "", paste(geneName, "expression"), "", "(Optimised)"))), fontSize=12, startCol=4,
              fontColor="darkblue", row.names = F, col.names = F) #, colSpan=1, rowSpan=1)

xlsx.addTable(wb, sheet, data= tableChi1, startCol=2,
              fontColor="darkblue", fontSize=12,
              rowFill=c("white", "white")
)
xlsx.addLineBreak(sheet, 2) # export the real row names of table 2, for copy and paste
xlsx.addTable(wb, sheet, data =featuresUni[-length(featuresUni)], fontSize=12, startCol=1,
              row.names = F, col.names = F)

#xlsx.addHeader(wb, sheet, value=paste("p-value of", geneName, "gene expression in OSCC"), level=5)

#xlsx.addTable(wb, sheet, data= chiT, startCol=2,
#              fontColor="darkblue", fontSize=14,
#              rowFill=c("white", "white")
#)
#xlsx.addLineBreak(sheet, 3)
# 
# 
# 
### Part II: table 3 and table 4 ####
# # {processing Table 3 OS: cbind(osHR, osHRmulti)
ciUniMti <- c("Features",	"HR",	"CI95%(L)",	"CI95%(H)",	"P-value",	"HR",	"CI95%(L)",	"CI95%(H)",	"P-value")
tableOS <- cbind(as.data.frame(colUni_1), osHR, osHRmulti)
rownames(tableOS) <- featuresUni
colnames(tableOS) <- ciUniMti
tableOS1 <- as.data.frame(setNames(replicate(length(ciUniMti), numeric(0), simplify = F), ciUniMti)) # for xlsx output (dynamic table)
colnames(tableOS1) <- ciUniMti
for (i in 1:nrow(tableOS)) {
  row0 <- data.frame(t(c(colUni_0[i], as.numeric(c(1, NA, NA, NA, 1, NA, NA, NA)))))
  colnames(row0) <- ciUniMti
  row0[,2:ncol(tableOS1)] <-as.numeric(as.character(row0[,2:ncol(tableOS1)]))
  tableOS1 <- rbind(tableOS1, row0, tableOS[i, ])
}
# }
# 
# 
# # {processing Table 4 RFS: cbind(rfsHR, rfsHRmulti)
#
ciUniMti <- c("Features",	"HR",	"CI95%(L)",	"CI95%(H)",	"P-value",	"HR",	"CI95%(L)",	"CI95%(H)",	"P-value")
tableRFS <- cbind(as.data.frame(colUni_1), rfsHR, rfsHRmulti) # it can't be rouned here ??
rownames(tableRFS) <- featuresUni
colnames(tableRFS) <- ciUniMti
tableRFS1 <- as.data.frame(setNames(replicate(length(ciUniMti), numeric(0), simplify = F), ciUniMti)) # for xlsx output (dynamic table)
colnames(tableRFS1) <- ciUniMti
for (i in 1:nrow(tableRFS)) {
  row0 <- data.frame(t(c(colUni_0[i], as.numeric(c(1, NA, NA, NA, 1, NA, NA, NA)))))
  colnames(row0) <- ciUniMti
  row0[,2:ncol(tableRFS1)] <- as.numeric(as.character(row0[,2:ncol(tableRFS1)]))
  tableRFS1 <- rbind(tableRFS1, row0, tableRFS[i, ])
}
# }
# 
# Export to workbook: using colSpan=5, rowSpan=4 to merge cells ; output row by row
# # fontColor="darkgray"

xlsx.addLineBreak(sheet, 5)
xlsx.addHeader(wb, sheet, value=paste("Table 3. Univariate/Multivariate Cox proportional hazards regression analyses on OS time of", geneName, "gene expression in ", TCGA_cohort), level=5)
xlsx.addLineBreak(sheet, 1)
# xlsx.addParagraph(wb, sheet, value = "Overall Survival", fontSize=12, isItalic=F, startCol=4,
#                  colSpan=8, rowSpan=1)
xlsx.addTable(wb, sheet, data = t(data.frame(c("Univariate", "\t", "\t", "\t", "Multivariate"))), fontSize=12, startCol=4,
              row.names = F, col.names = F) #, colSpan=4, rowSpan=1)
#
xlsx.addTable(wb, sheet, data = tableOS1, startCol=2,
              fontColor="darkblue", fontSize=12,
              rowFill=c("white", "lightblue")
)
xlsx.addLineBreak(sheet, 5)
xlsx.addHeader(wb, sheet, value=paste("Table 4. Univariate/Multivariate Cox proportional hazards regression analyses on RFS time of", geneName, "gene expression in ", TCGA_cohort), level=5)
xlsx.addLineBreak(sheet, 1)
xlsx.addTable(wb, sheet, data = t(data.frame(c("Univariate", "\t", "\t", "\t", "Multivariate"))), fontSize=12, startCol=4,
              row.names = F, col.names = F)
#
xlsx.addTable(wb, sheet, data= tableRFS1, startCol=2,
              fontColor="darkblue", fontSize=12,
              rowFill=c("white", "lightblue")
)
xlsx.addLineBreak(sheet, 5)

# Add paragraph
#+++++++++++++++++++++++++++++
#xlsx.addHeader(wb, sheet, "Add paragraph", level=2)
#xlsx.addLineBreak(sheet, 2)
#paragraph="Lorem Ipsum is simply dummy text of the printing and typesetting industry. Lorem Ipsum has been the industry's standard dummy text ever since the 1500s, when an unknown printer took a galley of type and scrambled it to make a type specimen book. It has survived not only five centuries, but also the leap into electronic typesetting, remaining essentially unchanged."
#xlsx.addParagraph(wb, sheet, paragraph, fontSize=14, isItalic=TRUE, 
#                  fontColor="darkred", backGroundColor="gray")
#xlsx.addLineBreak(sheet, 2)
# Add Hyperlink
#++++++++++++++ show cutoff +++++++++++++++
# xlsx.addParagraph(wb, sheet, value = paste("Kaplan-Meier survival estimate: \n Cutoff at", round(cutoff1[2], 3), "(", round(cutoff1[1], 2), ")", sep = ""), level=2, colSpan = 5, rowSpan = 2)
km <- paste("Kaplan-Meier survival estimate: Cutoff at ", round(cutoff1, 3), " (", names(cutoff1[1]), ")", sep = "")
xlsx.addParagraph(wb, sheet, value=km, isItalic=TRUE, colSpan=5, 
                  rowSpan=4, fontSize=14) # fontColor="darkgray"

xlsx.addLineBreak(sheet, 1)
xlsx.addHyperlink(wb, sheet, "http://www.sthda.com/english/wiki/r2excel-read-write-and-format-easily-excel-files-using-r-software", "Click-me!!", fontSize=12)
xlsx.addLineBreak(sheet, 2)

# library(forestplot) # is also fine; https://rd.springer.com/content/pdf/bbm%3A978-3-319-31245-3%2F1.pdf

# Add KM plots
#+++++++++++++++++++++++++++++
# xlsx.addHeader(wb, sheet, "Kaplan-Meier survival estimate: ", level=3)
xlsx.addLineBreak(sheet, 1)
# plotFunction <- function(){boxplot(len ~ dose, data = ToothGrowth, col = 1:3)}
plotFunction <- function() { # OS
  plot(OS.km, lty=1, xscale=12, xmax=60, col=c("blue","red"), sub=paste("P-value =", p_OS1), main=paste("OS in TCGA ", TCGA_cohort, "(n=", surv_OS1$n[1]+surv_OS1$n[2],")/", geneName), ylab="Percent Survival", xlab="Years")
  legend("topright", legend=c(paste("low(",surv_OS1$n[1], ")"), paste("high(",surv_OS1$n[2], ")")), lty=1:2, col=c("blue","red"))
}
xlsx.addPlot(wb, sheet, plotFunction())
plotFunction <- function() { # RFS
  plot(RFS.km, lty=1, xscale=12, xmax=60, col=c("blue","red"), sub=paste("p-value =", p_RFS1), main=paste("RFS in TCGA ", TCGA_cohort, "(n=", surv_RFS1$n[1]+surv_RFS1$n[2],")/", geneName), ylab="Percent Survival", xlab="Years")
  legend("topright", legend=c(paste("low(",surv_RFS1$n[1], ")"), paste("high(",surv_RFS1$n[2], ")")), lty=1:2, col=c("blue","red"))
}
xlsx.addPlot(wb, sheet, plotFunction())
#
xlsx.addLineBreak(sheet, 8)
# plot cumulative p-value from cutoff finder
plotFunction <- function() {
  pOS <- ggplot(OS, aes(x=cases_OS, y=p_OS)) + geom_point(size=2) +
    scale_x_discrete(limit = seq(0,500,50), labels = paste(seq(0,500,50), sep=",")) +
    scale_y_discrete(breaks = c(0, 0.05, 0.10, 0.50, 0.99)) + #, labels = c("0", "0.05", "0.10", "0.50", "0.99")) +
    coord_trans(y = "log10") +
    ggtitle(paste("Cutoff Finder: Cumulative P value curves for OS under ", geneName, " expression \n", "(cutoff at ", surv_OS1$n[1], " )", sep = "")) +
    xlab("# of patients") + ylab("P-value(log10)") +
    geom_hline(yintercept=0.05, linetype="dashed", color = "red", size=1, show.legend = T)
  print(pOS) # export to xlsx
}
xlsx.addPlot(wb, sheet, plotFunction())
#
plotFunction <- function() {
  pRFS <- ggplot(RFS, aes(x=cases_RFS, y=p_RFS)) + geom_point(size=2) +
    scale_x_discrete(limit = seq(0,500,50), labels = paste(seq(0,500,50), sep=",")) +
    scale_y_discrete(breaks = c(0, 0.05, 0.10, 0.50, 0.99)) + #, labels = c("0", "0.05", "0.10", "0.50", "0.99")) +
    coord_trans(y = "log10") +
    ggtitle(paste("Cutoff Finder: Cumulative P value curves for RFS under ", geneName, " expression \n", "(cutoff at ", surv_RFS1$n[1], " )", sep = "")) +
    xlab("# of patients") + ylab("P-value(log10)") +
    geom_hline(yintercept=0.05, linetype="dashed", color = "red", size=1, show.legend = T)
  geom_vline(xintercept=150, color = "green", size=2, show.legend = T)
  print(pRFS) # export to xlsx
  #  ggsave("plot.png", width=4, height=4, dpi=100) # save as .png
  #  dev.off() # "R, plot is done, please print it instead of on-screen device."
}
xlsx.addPlot(wb, sheet, plotFunction())
#
# create a row name table for copy/paste :-( redundent one row
# #xlsx.addParagraph
#for (i in 1:(length(featuresUni)-1)) {
#  xlsx.addParagraph(wb, sheet, value=featuresUni[i], isItalic=TRUE, startCol=12, colSpan=1,
#                      rowSpan=2, fontSize=14, fontColor="darkblue")
#}

xlsx.addLineBreak(sheet, 5) # export the p-value from Cutoff Finding
xlsx.addTable(wb, sheet, data =OS[OS$p_OS <= 0.05,1:2], fontSize=12, startCol=4,
              row.names = F, col.names = T)
#View(OS[OS$p_OS <= 0.05,1:2])
xlsx.addTable(wb, sheet, data =RFS[RFS$p_RFS <= 0.05,1:2], fontSize=12, startCol=5,
              row.names = F, col.names = T)
#View(RFS[RFS$p_RFS <= 0.05,1:2])
# save the workbook to an Excel file and write the file to disk.
xlsx::saveWorkbook(wb, filenamex)
# xlsx.openFile(filenamex) # open file to review
# +++++
# the END of R2Excel ###
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ R2Excel +++++++++

## (8.) google spreadsheet (with adds-on)  ####



################====the final END=======#################















##--------spared code------------####
# [method for copy and paste]
# 
#
# data <- read.csv(stdin())
# dput(data) # a command to recreate its data structure
# 
# # paste it in this script
# dat <- read.table(header=TRUE, text='
# H_score
# 135.9
#                   NaN
#                   NaN
#                   246.29')

#.........
# # % it doesn't work on Mac or Linux OS
# dat <- read.csv('clipboard',sep="\t",header=T)
# 
# read.Numbers <- function(header=TRUE,...) {
#   read.table('clipboard',sep="\t",header=header,...)
# }
# 
# dat <- read.Numbers()

# #
# write.Numbers <- function(x,row.names=FALSE,col.names=TRUE,...) {
#   write.table(x,"clipboard",sep="\t",row.names=row.names,col.names=col.names,...)
# }
# write.Numbers(dat)
# Well done.
