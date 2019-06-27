# debug: geneName <- whole_genome[main_i]
survival_marginFree <- function(geneName) { 
#  source(paste(path_LUAD, "cutofFinder_func.R", sep="")) # cutofFinder_func <- function(geneName) {} in cutofFinder_func.R
# for HNSCC: already defined on main, source(file=file.path(path_cohort, "cutofFinder_func_HNSCC.R")) # cutofFinder_func <- function(geneName) {} in cutofFinder_func.R
  
#  system.time( # as a timer
  #save(ZSWIM2, file=desti_ZSWIM2) # save it globally
  #global declaration was moved onto main.R
  #
  # for debug
  # sink(stdout(), type="message")
  # options(error=traceback, showWarnCalls=T)
  
# Declare as a function to be called
# 
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

# # https://cran.r-project.org/web/packages/ggplot2/
# # packages: stringi, stringr, reshape2
# download.file("https://cran.r-project.org/src/contrib/ggplot2_2.2.1.tar.gz", destfile="~/Downloads/ggplot2_2.2.1.tar.gz")
# # # >install.packages("~/Downloads/reshape2_1.4.3.tar.gz", repos = NULL, type = "source")
# install.packages("~/Downloads/ggplot2_2.2.1.tar.gz", repos = NULL, type = "source")
# library("ggplot2")



# 
# # brew cask install java
# ?sudo ln -f -s $(/usr/libexec/java_home)/jre/lib/server/libjvm.dylib /usr/local/lib
# # JDK 9 for Mac => http://www.oracle.com/technetwork/java/javase/downloads/jdk9-downloads-3848520.html
# # JAVA_HOME        : /Library/Java/JavaVirtualMachines/jdk-9.0.1.jdk/Contents/Home
# Java library path: $(JAVA_HOME)/lib/server
# 
# 
# x# DO NOT RUN this {https://www.r-bloggers.com/upgrade-r-without-losing-your-packages/
# tmp <- installed.packages()
# installedpkgs <- as.vector(tmp[is.na(tmp[,"Priority"]), 1])
# save(installedpkgs, file="installed_old.rda")
# # 
# # curl -#ROL https://cran.rstudio.com/bin/macosx/R-3.4.3.pkg
# # open R-3.4.3.pkg
# 
# brew reinstall R --with-openblas --with-java
# https://cran.rstudio.com/src/base/R-3/R-3.4.2.tar.gz
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




# [2018/03/06] with Peter Lai birthday
# 
# [Option-Cmd + E] for running to the end
# ***geneName<-"SLC2A9" # for debug
# Resume:[5.START] ZZZ3 and TMSB4X ###
# Start the survival analysis for each individual gene
# TCGA_cohort <- "HNSCC" # cancer type

# 6) call cutofFinder_func_HNSCC.R, run100; source defined on main ####
  marginTag <- "_marginFree_" # with margin 0 alone
  #browser()
  cutoffReturn <- cutofFinder_func(geneName, marginTag) # with return cutoff1 at # of patient
if (length(cutoffReturn) == 8) {
  
  osccCleanNA <- cutoffReturn[[1]] # taken for freq, chiT as well as tableChi1
  case50_n <- cutoffReturn[[2]] # of patients to cut
  cutoff1 <- cutoffReturn[[3]] # cutoff1 at RNAseq level
  surv_OS1 <- cutoffReturn[[4]]
  p_OS0 <- p_OS1 <- as.numeric(cutoffReturn[[5]]) # optimized P-Value
  OS <- cutoffReturn[[6]] # cases_OS and p_OS
  RFS <- cutoffReturn[[7]] # cases_RFS and p_RFS
  OS_pvalue <- cutoffReturn[[8]]
  osccCleanNAM_pos <- which(colnames(osccCleanNA) == paste("PMM1", "_median", sep=""))
  
# => find cutoff1, generate OS P-Value plot and KM plot.
# It was masked _by_ .GlobalEnv:  cases_OS, p_OS
# 
#**spared code:  cutoffReturn <- cutofFinder_func(geneName) # with return cutoff1 at # of patients
  #p_RFS0 <- cutoffReturn[[2]]
  # => find cutoff1, generate RFS P-Value plot and KM plot.
 
} else if (length(cutoffReturn) == 1) {
  return(cutoffReturn) # return with error code
  }


## 7)[COXPH modelling] ####
library(rJava) # PASSED
library("xlsx")
# # solved ERROR: compilation failed for package ‘git2r’; ERROR: dependencies ‘httr’, ‘git2r’ are not available for package ‘devtools’
# # $ brew install openssl # macOS
# $ sudo apt-get install openssl # linux
# $ sudo apt-get install libssl-dev
library("openssl")
# install.packages('git2r', repos="https://github.com/ropensci/git2r", type='source')
# # or clone then build (make) in BASH
# # $ git clone https://github.com/ropensci/git2r.git
# # $ cd git2r
# # $ make install
# # it response as * DONE (git2r) for R version 3.4.3)
# 
# #xx no .c file $ git clone https://github.com/r-lib/httr
# install.packages("~/Downloads/httr_1.3.1.tar.gz", repos = NULL, type = "source")
# 
# install.packages("devtools")
# library(devtools)
#install_github("kassambara/r2excel", force=T) # the only way OK; reinstall it
library(r2excel) # (package ‘r2excel’ is available for R version 3.4.1 below)
# x installing curl directly then r2excel
# install r2excel by install_github()



# *** [cont] all features needs to be adjusted one by one; updated on [2018/03/09] for LUAD osccCleanNA

#[coxph/ageDx]  Export part I, univariate survival: ##
## => coxph, for gene high, low and other features #
## multivariant survival: coxph, Cox proportional Hazard regression Model
## covariates (group, sex...)
# http://www.uni-kiel.de/psychologie/rexrepos/posts/survivalCoxPH.html

#[coxph] - fits a Cox proportional hazards regression model
#https://rstudio-pubs-static.s3.amazonaws.com/5588_72eb65bfbe0a4cb7b655d2eee0751584.html

#<<vdxfitcox>> to fit a Cox regression model with the 7 genes.
#Multivariate Cox proportional hazards regression analyses on recurrence-free survival time in gene expression datasets.
#The status indicator, normally 0=alive, 1=dead. Other choices are TRUE/FALSE (TRUE = death) or 1/2 (2=death). For interval censored data, the status indicator is 0=right censored, 1=event at time, 2=left censored, 3=interval censored. A
# time (months) is OVERALL.SURVIVAL OS or RECURRENCE.FREE.SURVIVAL RFS
# status = overall.survival.indicator or recurrence.free.survival.indicator
# library(survival)


#*** [Multivariate for OS] Table 3 right panle

# selectedFeatures <- colnames(osccCleanNA[commonFeatures])
# colname of oscc = c("Unique.ID","Gender","age.at.diagnosis",
# "margin","T","N",
# "M","stage_2","OS..months._from.biopsy",
# "OS_IND","RFS..months._from.op","RFS_IND","H.score_T")

# ***rename all of colnames (features) as osccT [old coding]
colnames(osccCleanNA) <- coln_osccT # colnames (features) of osccT

osHRmulti <- 0
# 8 + OS/RFS features in HNSCC
# #warning() X matrix deemed to be singular (margin) in coxph
# https://stackoverflow.com/questions/20977401/coxph-x-matrix-deemed-to-be-singular
oscox <- coxph(Surv(osccCleanNA$OS..months._from.biopsy, osccCleanNA$OS_IND==1) ~ +
  Gender +
  ageDx +
#                 primary_site +
                 clinical_T +
                 clinical_N +
                 clinical_M +
#  pathologic_T +
#  pathologic_N +
#  pathologic_M +
  stage +
  margin +
#  R.T +
#  C.T +
#                 presence_of_pathological_nodal_extracapsular_spread +
#                 neoplasm_histologic_grade +
  as.vector(osccCleanNA[, osccCleanNAM_pos], mode="numeric"), 
  data=osccCleanNA) # PMM1_median == IHC score == (H.score_T) in [low, high]
 # warning() X matrix deemed to be singular; variable 8:margin
# summary(oscox)

osHRmulti <- cbind(data.frame(summary(oscox)$conf.int)[-2], data.frame(summary(oscox)$coefficients)[5])

skipNA <- rownames(osHRmulti) %in% c(1, "marginNA")
osHRmulti <- osHRmulti[which(!skipNA), ]

# correction of rownames
rownames(osHRmulti) <- featuresUni
#rownames(osHRmulti)[nrow(osHRmulti)] <- paste("PMM1", "_median", sep="")
osHRmulti <- round(osHRmulti, 3) # rounding them by 3 decimal

#View(osHRmulti) # table 3 right panel



#*** [Univariate for OS]  table 3 left panel
# # 8 features in LUAD, put in coxph one by one
# number of features 13 -> 14 (add a feature: margin status)
oscox <- 0
features_os <- colnames(osccCleanNA)[c(2:8,14)] # features selection ["Gender" to "margin", "PMM1"] 8 out of 14
osHR <- data.frame(matrix(0, nrow=length(features_os), ncol=4))
## *** looping the cox regression model over several features
## https://stackoverflow.com/questions/13092923/looping-cox-regression-model-over-several-predictor-variables
## as.formula: text to code (class: forumla list)
coxph_func <- function(x) as.formula(paste("Surv(osccCleanNA$OS..months._from.biopsy, osccCleanNA$OS_IND==1)", x, sep="~"))
formlist <- lapply(features_os, coxph_func)
#coxph_func2 <- function(x) as.formula(paste(x, ',', "data=osccCleanNA", sep=""))
#formlist2 <- lapply(formlist, coxph_func2)
#oscox <- coxph(coxph_func(x), data = osccCleanNA, ties="efron")
oscox <- lapply(formlist, coxph, data = osccCleanNA) # run coxph, "data=" as additional arguments to FUN. 
#-> model oscox (list 8)
# warnings() In FUN(X[[i]], ...) : X matrix deemed to be singular; variable 2

coef_func <- function(i, oscox1) {
  x1 <- data.frame(summary(oscox1[[i]])$conf.int)[-2]
  x2 <- data.frame(summary(oscox1[[i]])$coefficients)[5]
  return(list(x1, x2))
}
#osHR <- rbind(osHR, cbind(unlist(uni_CI)))
uni_CI <- lapply(c(1:length(oscox)), coef_func, oscox) # for (i = 1:8)
for (ii in c(1:length(uni_CI))) {
  # ii=7 has 8 values in two rows <- marginNA issue => pick up first row by index [1, ]
  osHR[ii, ] <- t(c(unlist(uni_CI[[ii]][[1]][1,]), unlist(uni_CI[[ii]][[2]][1,])))
}
# > ii<-7; print(unlist(uni_CI[[ii]]))
# exp.coef.1   exp.coef.2   lower..951   lower..952 
# 3.4057062147           NA 1.7000673239           NA 
# upper..951   upper..952    Pr...z..1    Pr...z..2 
# 6.8225738228           NA 0.0005463051           NA
# => marginNA with NA

#skipNA <- rownames(osHR) %in% c(1) #, "marginNA")
#osHR <- osHR[which(!skipNA), ]
# correction of rownames
rownames(osHR) <- featuresUni
# round
osHR <- round(osHR, 3)
# P-value notes:
# Significant codes:  0 ???***??? 0.001 ???**??? 0.01 ???*??? 0.05 ???.??? 0.1 ??? ??? 1
# if <0.001 => mark as "***"
osHR$X4[osHR$X4<0.001] <- "***"
#View(osHR) #-> (univariate) table 3 left panel



# # [skip]
# # test for the assumption of proportionality of hazards
# cox.zph(oscox)
# plot(cox.zph(oscox)) # flat curve is good, not violate asumption
# ##


## *** [Multivariate for RFS] table 4 right
## # 8 features in HNSCC
rfsHRmulti <- 0
# P-value notes:
#Significant codes:  0 ???***??? 0.001 ???**??? 0.01 ???*??? 0.05 ???.??? 0.1 ??? ??? 1
rfscox <- coxph(Surv(osccCleanNA$RFS..months._from.op, osccCleanNA$RFS_IND==1) ~ +
                  Gender +
                  ageDx +
                  #                 primary_site +
                                   clinical_T +
                                   clinical_N +
                                   clinical_M +
                  #pathologic_T +
                  #pathologic_N +
                  #pathologic_M +
                  stage +
                  margin +
                  # R.T +
                  # C.T +
                  #                 presence_of_pathological_nodal_extracapsular_spread +
                  #                 neoplasm_histologic_grade +
                  as.vector(unlist(osccCleanNA[, osccCleanNAM_pos]), mode="numeric"), data=osccCleanNA)

# summary(rfscox)

rfsHRmulti <- cbind(data.frame(summary(rfscox)$conf.int)[-2], data.frame(summary(rfscox)$coefficients)[5])
#rownames(rfsHRmulti)[nrow(rfsHRmulti)] <- paste("PMM1", "_median", sep="")
skipNA <- rownames(rfsHRmulti) %in% c("marginNA")
rfsHRmulti <- rfsHRmulti[which(!skipNA), ]
# correction of rownames
rownames(rfsHRmulti) <- featuresUni

rfsHRmulti <- round(rfsHRmulti, 3)
#View(rfsHRmulti)

# # [skip]
# test for the assumption of proportionality of hazards
# cox.zph(rfscox)
# plot(cox.zph(rfscox))



#*** [Univariate for RFS]  table 4 left panel
# # 8 features in HNSCC, put in coxph one by one
# number of features 13 -> 14 (add a feature: margin status)
rfscox <- 0
features_os <- colnames(osccCleanNA)[c(2:8,14)] # features selection ["Gender" to "margin", "PMM1"] clincical TNM, 8 out of 14
rfsHR <- data.frame(matrix(0, nrow=length(features_os), ncol=4))
## *** looping the cox regression model over several features
## as.formula: text to code (class: forumla list)
coxph_func <- function(x) as.formula(paste("Surv(osccCleanNA$OS..months._from.biopsy, osccCleanNA$OS_IND==1)", x, sep="~"))
formlist <- lapply(features_os, coxph_func)
#coxph_func2 <- function(x) as.formula(paste(x, ',', "data=osccCleanNA", sep=""))
#formlist2 <- lapply(formlist, coxph_func2)
#rfscox <- coxph(coxph_func(x), data = osccCleanNA, ties="efron")
rfscox <- lapply(formlist, coxph, data = osccCleanNA) # run coxph, "data=" as additional arguments to FUN. 
#-> model rfscox (list 8)
# warnings() In FUN(X[[i]], ...) : X matrix deemed to be singular; variable 2

coef_func <- function(i, rfscox1) {
  x1 <- data.frame(summary(rfscox1[[i]])$conf.int)[-2]
  x2 <- data.frame(summary(rfscox1[[i]])$coefficients)[5]
  return(list(x1, x2))
}
#rfsHR <- rbind(rfsHR, cbind(unlist(uni_CI)))
uni_CI <- lapply(c(1:length(rfscox)), coef_func, rfscox) # for (i = 1:8)
for (ii in c(1:length(uni_CI))) {
  # ii=7 has 8 values in two rows <- marginNA issue => pick up first row by index [1, ]
  rfsHR[ii, ] <- t(c(unlist(uni_CI[[ii]][[1]][1,]), unlist(uni_CI[[ii]][[2]][1,])))
}
# > ii<-7; print(unlist(uni_CI[[ii]]))
# exp.coef.1   exp.coef.2   lower..951   lower..952 
# 3.4057062147           NA 1.7000673239           NA 
# upper..951   upper..952    Pr...z..1    Pr...z..2 
# 6.8225738228           NA 0.0005463051           NA
# => marginNA with NA

#skipNA <- rownames(rfsHR) %in% c(1) #, "marginNA")
#rfsHR <- rfsHR[which(!skipNA), ]
# correction of rownames
rownames(rfsHR) <- featuresUni
# round
rfsHR <- round(rfsHR, 3)
# P-value notes:
# Significant codes:  0 ???***??? 0.001 ???**??? 0.01 ???*??? 0.05 ???.??? 0.1 ??? ??? 1
# if <0.001 => mark as "***"
rfsHR$X4[rfsHR$X4<0.001] <- "***"
#View(rfsHR) # univariate table
## DONE and then export part II ##




##8) R2Excel export ####
## # http://www.sthda.com/english/wiki/r2excel-read-write-and-format-easily-excel-files-using-r-software
library("xlsx")
library("r2excel")
# Create an Excel workbook. Both .xls and .xlsx file formats can be used.
#TCGA_cohort <- "HNSCC" # cancer type
#path_LUAD <- "~/R/LUAD_Peter_survival/" # under rstudio-server on GCP
#path_LUAD <- "/Users/apple/Google\ Drive/2018_PhD_dissertation/LUAD_Peter_survival/"
# path_LUAD <- "/Users/apple/Documents/My\ Tableau\ Repository/Workbooks/"
setwd(path_cohort) # change working directory to the HNSCC, GCP

# filename for all margins cases, defined as HNSCC
filenamex <- paste("xlsx/", TCGA_cohort, "_survivalAnalysis_marginFree_", geneName, ".xlsx", sep = "")
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



# under optimized cutoff1
# Part I: # Table 2 construction by Calling contingencyTCGA and contignecyBin ####
# ??using tableChi2 by contigencyBin2 with c("Gender","ageDx", "pathologic_T", "pathologic_N", "pathologic_M", "stage","margin" )
# under beset cutoff value (auto choosen)
# Table 2: clinical correlation tables (tableChi1)
contiT <- contingencyTCGA(osccCleanNA, geneName) # calling this function (OSCC cohort, PMM1) 2 parameters, no more cutoff1; 
# "margin" at column 8 of oscc
#  oscc <- contiT[[1]] # updating PMM1_median
chiT <- contiT[[2]] # extrac it from list by using [[]]; chiT$X2 is the P-value
freq <- contiT[[3]] # well DONE
tableChi1 <- contingencyBin (osccCleanNA, chiT, freq) # calculating the P-value
# ***add rownames
nrow_featuresUni <- length(featuresUni) # aka. 8
nrow_tableChi1 <- as.character(seq(1, 2*(nrow_featuresUni-1))) # aka. 7 x2 NA (a vector)
nrow_tableChi1[c(T,F)] <- featuresUni[-nrow_featuresUni] # skip last one row (z-score)
# > nrow_tableChi1
# [1] "Gender"                 NA                      
# [3] "Age at diagnosis"       NA                      
# [5] "clinical T status"    NA                      
# [7] "clinical N status"    NA                      
# [9] "clinical M status"    NA                      
# [11] "clinical Stage"       NA                      
# [13] "Surgical Margin status" NA     
rownames(tableChi1) <- nrow_tableChi1 # duplicate 'row.names' are not allowed
# "margin" at column 8 of osccCleanNA

# to save "*" of column "remark" in 300 tableChi1 ###
#?? cut_featuresRemark <- cbind (cut_featuresRemark, tableChi1[7], tableChi1$Remark) # retrieve the P-value and remark
#

# header
xlsx.addHeader(wb, sheet, value=paste("Table 2. The correlation of", geneName,  "expression and clinical features."),
               level=5, color="black", underline=0)
xlsx.addHeader(wb, sheet, value=paste("Cutoff at ", round(cutoff1, 3), " (", percent(surv_OS1$n[1]/(surv_OS1$n[1]+surv_OS1$n[2])), ")", sep = ""),
               level=5, color="red", underline=0) # total n is taken from surv_OS1$n

#xlsx.addLineBreak(sheet, 1) # add one blank line

xlsx.addTable(wb, sheet, data = t(data.frame(c(paste(geneName, "expression"), "", paste(geneName, "expression"), "", "(Optimised)"))), fontSize=12, startCol=4,
              fontColor="darkblue", row.names = F, col.names = F) #, colSpan=1, rowSpan=1)
# tableChi1
xlsx.addTable(wb, sheet, data= tableChi1, startCol=2,
              fontColor="darkblue", fontSize=12,
              rowFill=c("white", "white"), row.names = TRUE)

xlsx.addLineBreak(sheet, 2)  # add two blank lines

# export the real row names of table 2, for copy and paste
# xlsx.addTable(wb, sheet, data = featuresUni[-length(featuresUni)], fontSize=12, fontColor="black",
#               row.names = F, col.names = F) # 8 features => -8 as removal of last feature z-score
# #[2018/03/25] why ? (Error in 1:col.n : argument of length 0) ####

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
# # {processing Table 3 OS: cbind(osHR, osHRmulti) => tableOS1
#ciUniMti <- c("Features",	"HR",	"CI95%(L)",	"CI95%(H)",	"P-value",	"HR",	"CI95%(L)",	"CI95%(H)",	"P-value")
tableOS <- cbind(as.data.frame(colUni_1), osHR, osHRmulti) # colUni_1: values of feature [male/female...]
rownames(tableOS) <- featuresUni # 8 features of rownames [gender, age....]
colnames(tableOS) <- ciUniMti # colnames

# for xlsx output (dynamic tableOS)
# setNames, This is a convenience function that sets the names on an object and returns the object.
# a.k.a. create a empty data.frame: tableOS1; however, "()" will be assigned as "."
tableOS1 <- as.data.frame(setNames(replicate(length(ciUniMti), numeric(0), simplify = F), ciUniMti)) 
colnames(tableOS1) <- ciUniMti
for (i in 1:nrow(tableOS)) { # eg. 1~8
  # row0
  row0 <- data.frame(t(c(colUni_0[i], as.numeric(c(1, NA, NA, NA, 1, NA, NA, NA)))))
  colnames(row0) <- ciUniMti
  row0[,2:ncol(tableOS1)] <-as.numeric(as.character(row0[,2:ncol(tableOS1)]))
  # row1 == tableOS
  tableOS1 <- rbind(tableOS1, row0, tableOS[i, ]) # rbind every 2 rows
  # (OK) warnings() In as.numeric(as.character(row0[, 2:ncol(tableOS1)])) :
  # NAs introduced by coercion
}
# }
# 
# 
# # {processing Table 4 RFS: cbind(rfsHR, rfsHRmulti) => tableRFS1
#
#ciUniMti <- c("Features",	"HR",	"CI95%(L)",	"CI95%(H)",	"P-value",	"HR",	"CI95%(L)",	"CI95%(H)",	"P-value")
tableRFS <- cbind(as.data.frame(colUni_1), rfsHR, rfsHRmulti) # it can't be rouned here ??
rownames(tableRFS) <- featuresUni
colnames(tableRFS) <- ciUniMti
tableRFS1 <- as.data.frame(setNames(replicate(length(ciUniMti), numeric(0), simplify = F), ciUniMti)) # for xlsx output (dynamic table)
colnames(tableRFS1) <- ciUniMti
for (i in 1:nrow(tableRFS)) {
  # row0
  row0 <- data.frame(t(c(colUni_0[i], as.numeric(c(1, NA, NA, NA, 1, NA, NA, NA)))))
  colnames(row0) <- ciUniMti
  row0[,2:ncol(tableRFS1)] <- as.numeric(as.character(row0[,2:ncol(tableRFS1)]))
  # row1 == tableRFS
  tableRFS1 <- rbind(tableRFS1, row0, tableRFS[i, ])
  # (OK) warnings() In as.numeric(as.character(row0[, 2:ncol(tableOS1)])) :
  # NAs introduced by coercion
}
# }



# Export to tables workbook: using colSpan=5, rowSpan=4 to merge cells ####
# output row by row (tableOS1, tableRFS1)
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
              rowFill=c("white", "lightblue"), row.names = TRUE
)
xlsx.addLineBreak(sheet, 5)
xlsx.addHeader(wb, sheet, value=paste("Table 4. Univariate/Multivariate Cox proportional hazards regression analyses on RFS time of", geneName, "gene expression in ", TCGA_cohort), level=5)
xlsx.addLineBreak(sheet, 1)
xlsx.addTable(wb, sheet, data = t(data.frame(c("Univariate", "\t", "\t", "\t", "Multivariate"))), fontSize=12, startCol=4,
              row.names = F, col.names = F)
#
xlsx.addTable(wb, sheet, data= tableRFS1, startCol=2,
              fontColor="darkblue", fontSize=12,
              rowFill=c("white", "lightblue"), row.names = TRUE
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
km <- paste("Kaplan-Meier survival estimate: Cutoff at ", round(cutoff1, 3), " (", percent(surv_OS1$n[1]/(surv_OS1$n[1]+surv_OS1$n[2])), ")", sep = "")
xlsx.addParagraph(wb, sheet, value=km, isItalic=TRUE, colSpan=5, 
                  rowSpan=4, fontSize=14) # fontColor="darkgray"

# xlsx.addLineBreak(sheet, 1)
# xlsx.addHyperlink(wb, sheet, "http://www.sthda.com/english/wiki/r2excel-read-write-and-format-easily-excel-files-using-r-software", "Click-me!!", fontSize=12)
# xlsx.addLineBreak(sheet, 2)
# # a data.frame to write to the workbook
# xlsx.writeFile(wb, file=filenamex, 
#                sheetName=sheet, append=FALSE)
#xlsx::saveWorkbook(wb, filenamex)

# library(forestplot) # is also fine; https://rd.springer.com/content/pdf/bbm%3A978-3-319-31245-3%2F1.pdf



# Export KM plots by plotFunction() ####
#+++++++++++++++++++++++++++++
# xlsx.addHeader(wb, sheet, "Kaplan-Meier survival estimate: ", level=3)
xlsx.addLineBreak(sheet, 1)

# # # debug
# plot(OS.km, lty=1, xscale=12, xmax=60, col=c("blue","red"), sub=paste("P-value =", p_OS1), main=paste("OS in TCGA", TCGA_cohort, "(n=", surv_OS1$n[1]+surv_OS1$n[2],")/", geneName), ylab="Percent Survival", xlab="Years")
# legend("topright", legend=c(paste("low(",surv_OS1$n[1], ")"), paste("high(",surv_OS1$n[2], ")")), lty=1:2, col=c("blue","red"))
# # # debug

# below code was copied from cutofFinder_func.R
mysurv <- Surv(osccCleanNA$OS..months._from.biopsy, osccCleanNA$OS_IND==1) #1==dead
# Test for difference (log-rank test) between groups (by PMM1_median 0 vs 1)
tryCatch(surv_OS1 <- survdiff(mysurv ~ as.vector(osccCleanNA[, osccCleanNAM_pos], mode="numeric"), data=osccCleanNA), error = function(e) return(NA)) # PMM1 high or low
# pchisq gives the distribution function
#p_OS1 <- format(pchisq(surv_OS1$chisq, length(surv_OS1$n)-1, lower.tail = FALSE), digits=3)
# (p_OS1 == p_OS0) is TRUE
#cases_OS1 <- surv_OS1$n[1]
OS.km <- survfit(mysurv ~ as.vector(unlist(osccCleanNA[, osccCleanNAM_pos]), mode="numeric"), data=osccCleanNA, type= "kaplan-meier", conf.type = "log-log")

xlsx.addPlot.OS<-function(OS.km, surv_OS1, wb, sheet, startRow=NULL, startCol=2,
                        width=480, height=480,... )
{ # OS KM plot
  #library("xlsx")
  
  png(filename = "plot.png", width = width, height = height,...)
  plot(OS.km, lty=1, xscale=12, xmax=60, col=c("blue","red"), sub=paste("optimized P-Value =", p_OS1), main=paste("OS in TCGA ", TCGA_cohort, "(n=", surv_OS1$n[1]+surv_OS1$n[2],")/", geneName), ylab="Percent Survival", xlab="Years")
  legend("topright", legend=c(paste("low(",surv_OS1$n[1], ")"), paste("high(",surv_OS1$n[2], ")")), lty=1:1, col=c("blue","red"))
  
  dev.off() 
  #Append plot to the sheet
  if(is.null(startRow)){
    rows<- getRows(sheet) #list of row object
    startRow=length(rows)+1
  } 
  # Add the file created previously
  addPicture("plot.png", sheet=sheet,  startRow = startRow, startColumn = startCol) 
  xlsx.addLineBreak(sheet, round(width/20)+1)
  res<-file.remove("plot.png")
}
# calling
xlsx.addPlot.OS(OS.km, surv_OS1, wb, sheet)

#
# below code was copied from cutofFinder_func.R
#RFS.km <- survfit(mysurv ~ as.vector(unlist(osccCleanNA[, osccCleanNAM_pos]), mode="numeric"), data=osccCleanNA, type= "kaplan-meier", conf.type = "log-log")

xlsx.addPlot.RFS<-function(RFS.km, surv_RFS1, wb, sheet, startRow=NULL, startCol=2,
                           width=480, height=480,... )
{  # RFS KM plot
  #library("xlsx")
  png(filename = "plot.png", width = width, height = height,...)
  # plot fuction here
  plot(RFS.km, lty=1, xscale=12, xmax=60, col=c("blue","red"), sub=paste("optimized P-Value =", p_RFS1), main=paste("RFS in TCGA ", TCGA_cohort, "(n=", surv_RFS1$n[1]+surv_RFS1$n[2],")/", geneName), ylab="Percent Survival", xlab="Years")
  legend("topright", legend=c(paste("low(",surv_RFS1$n[1], ")"), paste("high(",surv_RFS1$n[2], ")")), lty=1:1, col=c("blue","red"))
  
  dev.off() 
  #Append plot to the sheet
  if(is.null(startRow)){
    rows<- getRows(sheet) #list of row object
    startRow=length(rows)+1
  } 
  # Add the file created previously
  addPicture("plot.png", sheet=sheet,  startRow = startRow, startColumn = startCol) 
  xlsx.addLineBreak(sheet, round(width/20)+1)
  res<-file.remove("plot.png")
}
# # NOT calling
# xlsx.addPlot.RFS(RFS.km, surv_RFS1, wb, sheet)
# #

xlsx.addLineBreak(sheet, 8)



# plot cumulative p-value from cutoff finder
xlsx.addPlot.OSpval<-function(OS, case50_n, p_OS0, wb, sheet, startRow=NULL,startCol=2,
                            width=480, height=480,... )
{   # OS P-values dots plot
  #library("xlsx")
  png(filename = "plot.png", width = width, height = height,...)
  # plot fuction here
  g1<- subset(OS, cases_OS==case50_n) # optimized cutoff1 at (x=case50_n, y=p_OS0)
  
  if (g1$p_OS<=0.05) {
    pOS <- ggplot(OS, aes(x=cases_OS, y=p_OS)) + geom_point(size=2) +
      #  xlim(70, cutoff_n[2]) +
      scale_x_discrete(limit = seq(80, 180, 20)) + # cutoff_n[2])) + #, labels = paste(seq(0,500,50), sep=",")) +
      #  scale_x_discrete(limit = seq(0,500,50), labels = paste(seq(0,500,50), sep=",")) +
      scale_y_log10(breaks=c(p_OS0, 1e-05, 1e-03, 0.05), labels=c(p_OS0, "0.00001", 0.001, 0.05)) + #(limits=break_y, labels=c("0", "0.05", "0.10", "0.50")) +
      #  scale_y_discrete(breaks = break_y, labels = break_y) + #c("0", "0.05", "0.10", "0.50", "0.99")) +
      #  coord_trans(y = "log10") +
      ggtitle(paste("Cumulative P-Value plot for OS under", geneName, "expression")) +
      xlab("# of patients in 'low exp' group") + ylab("P-Value(log10)") +
      geom_hline(yintercept=0.05, linetype="dashed", color = "red", size=1, show.legend = T) +
      geom_point(data=g1, color="red") +  # this adds a red point
      geom_text(data=g1, label="optimized", color="red", hjust=-0.5) # this adds a label for the red point
    #geom_vline(xintercept=case50_n, linetype="dashed", color = "green", size=2, show.legend = T)
  } else { #(g1$p_OS > 0.05)
    pOS <- ggplot(OS, aes(x=cases_OS, y=p_OS)) + geom_point(size=2) +
      #  xlim(70, cutoff_n[2]) +
      scale_x_discrete(limit = seq(80, 180, 20)) + # cutoff_n[2])) + #, labels = paste(seq(0,500,50), sep=",")) +
      #  scale_x_discrete(limit = seq(0,500,50), labels = paste(seq(0,500,50), sep=",")) +
      #scale_y_log10(breaks=c(p_OS0, 1e-05, 1e-03, 0.05), labels=c(p_OS0, "0.00001", 0.001, 0.05)) + #(limits=break_y, labels=c("0", "0.05", "0.10", "0.50")) +
      #  scale_y_discrete(breaks = break_y, labels = break_y) + #c("0", "0.05", "0.10", "0.50", "0.99")) +
      #  coord_trans(y = "log10") +
      ggtitle(paste("Cumulative P-value plot for OS under", geneName, "expression")) +
      xlab("# of patients in 'low exp' group") + ylab("P-Value(log10)") +
      geom_hline(yintercept=0.05, linetype="dashed", color = "red", size=1, show.legend = T) +
      geom_point(data=g1, color="yellow") +  # this adds a red point
      geom_text(data=g1, label="50:50", color="red", hjust=-0.5, vjust=1.0) # this adds a label for the red point
    #geom_vline(xintercept=case50_n, linetype="dashed", color = "green", size=2, show.legend = T)
    
  }
  
  print(pOS) # export to xlsx
  dev.off() 
  #Append plot to the sheet
  if(is.null(startRow)){
    rows<- getRows(sheet) #list of row object
    startRow=length(rows)+1
  } 
  # Add the file created previously
  addPicture("plot.png", sheet=sheet,  startRow = startRow, startColumn = startCol) 
  xlsx.addLineBreak(sheet, round(width/20)+1)
  res<-file.remove("plot.png")
}
# calling
xlsx.addPlot.OSpval(OS, case50_n, p_OS0, wb, sheet)
#
#
#
#
xlsx.addPlot.RFSpval<-function(RFS, case50_n, p_RFS0, wb, sheet, startRow=NULL,startCol=2,
                               width=480, height=480,... )
{   # RFS P-values dots plot
  #library("xlsx")
  png(filename = "plot.png", width = width, height = height,...)
  # plot fuction here
  g1<- subset(RFS, cases_RFS==case50_n) # optimized cutoff1 at (x=case50_n, y=p_RFS0)
  
  if (g1$p_RFS<=0.05) {
    pRFS <- ggplot(RFS, aes(x=cases_RFS, y=p_RFS)) + geom_point(size=2) +
      #  xlim(70, cutoff_n[2]) +
      scale_x_discrete(limit = seq(80, 180, 20)) + # cutoff_n[2])) + #, labels = paste(seq(0,500,50), sep=",")) +
      #  scale_x_discrete(limit = seq(0,500,50), labels = paste(seq(0,500,50), sep=",")) +
      scale_y_log10(breaks=c(p_RFS0, 1e-05, 1e-03, 0.05), labels=c(p_RFS0, "0.00001", 0.001, 0.05)) + #(limits=break_y, labels=c("0", "0.05", "0.10", "0.50")) +
      #  scale_y_discrete(breaks = break_y, labels = break_y) + #c("0", "0.05", "0.10", "0.50", "0.99")) +
      #  coord_trans(y = "log10") +
      ggtitle(paste("Cumulative P-Value plot for RFS under", geneName, "expression")) +
      xlab("# of patients in 'low exp' group") + ylab("P-Value(log10)") +
      geom_hline(yintercept=0.05, linetype="dashed", color = "red", size=1, show.legend = T) +
      geom_point(data=g1, color="red") +  # this adds a red point
      geom_text(data=g1, label="optimized", color="red", hjust=-0.5) # this adds a label for the red point
    #geom_vline(xintercept=case50_n, linetype="dashed", color = "green", size=2, show.legend = T)
  } else { #(g1$p_OS > 0.05)
    pRFS <- ggplot(RFS, aes(x=cases_RFS, y=p_RFS)) + geom_point(size=2) +
      #  xlim(70, cutoff_n[2]) +
      scale_x_discrete(limit = seq(80, 180, 20)) + # cutoff_n[2])) + #, labels = paste(seq(0,500,50), sep=",")) +

      ggtitle(paste("Cumulative P-value plot for RFS under", geneName, "expression")) +
      xlab("# of patients in 'low exp' group") + ylab("P-Value(log10)") +
      geom_hline(yintercept=0.05, linetype="dashed", color = "red", size=1, show.legend = T) +
      geom_point(data=g1, color="yellow") +  # this adds a red point
      geom_text(data=g1, label="50:50", color="red", hjust=-0.5, vjust=1.0) # this adds a label for the red point
    #geom_vline(xintercept=case50_n, linetype="dashed", color = "green", size=2, show.legend = T)
    
  }
  
  print(pRFS) # export to xlsx
  #  ggsave("plot.png", width=4, height=4, dpi=100) # save as .png
  #  dev.off() # "R, plot is done, please print it instead of on-screen device."
  dev.off() 
  #Append plot to the sheet
  if(is.null(startRow)){
    rows<- getRows(sheet) #list of row object
    startRow=length(rows)+1
  } 
  # Add the file created previously
  addPicture("plot.png", sheet=sheet,  startRow = startRow, startColumn = startCol) 
  xlsx.addLineBreak(sheet, round(width/20)+1)
  res<-file.remove("plot.png")
}
# # NOT calling
# xlsx.addPlot.RFSpval(RFS, case50_n, p_RFS0, wb, sheet)
# #


##
#
#
# create a row name table for copy/paste :-( redundent one row
# #xlsx.addParagraph
#for (i in 1:(length(featuresUni)-1)) {
#  xlsx.addParagraph(wb, sheet, value=featuresUni[i], isItalic=TRUE, startCol=12, colSpan=1,
#                      rowSpan=2, fontSize=14, fontColor="darkblue")
#}

xlsx.addLineBreak(sheet, 5) 
# export the p-value from Cutoff Finding
xlsx.addTable(wb, sheet, data = subset(OS, p_OS <= 0.05)[order(subset(OS, p_OS <= 0.05)$p_OS),], fontSize=12, startCol=8,
              row.names = F, col.names = T)
#View(OS[OS$p_OS <= 0.05,1:2])
xlsx.addTable(wb, sheet, data = subset(RFS, p_RFS <= 0.05)[order(subset(RFS, p_RFS <= 0.05)$p_RFS),], fontSize=12, startCol=9,
              row.names = F, col.names = T)
#View(RFS[RFS$p_RFS <= 0.05,1:2])
#

# save the workbook to an Excel file and write the file to disk.
xlsx::saveWorkbook(wb, filenamex)

print(paste("case50_n=",case50_n,";", filenamex, "successfully."))
# xlsx.openFile(filenamex) # open file to review
# +++++
# the END of R2Excel ###


## 
## save to  individual gene to an R data file.####

# Added code for table 2, table 3, table 4 and survival p-value save as .Rda
# save(list = c("tableChi1", "tableOS1", "tableRFS1"), file=paste("LUAD_survivalAnalysis_marginS_", geneName, ".Rda"))

#tryCatch(
  RFS_pvalue <- OS_pvalue #:-) for HNSCC only
  save(list = c("tableChi1", "tableOS1", "tableRFS1", "OS_pvalue", "RFS_pvalue"), file=paste("xlsx/HNSCC_survivalAnalysis_marginFree_", geneName, ".Rda", sep=""))
#, error = function(e) return(NA))
#print(paste("Create", paste("LUAD_survivalAnalysis_marginFree_", geneName, ".Rda", sep=""), "successfully."))

 
## 
# Final Return ####
  if (nrow(OS_pvalue) > 0)  { # we hit one gene with P-value < 0.05 in KM plot
    return(which.min(OS_pvalue$p_OS))} else {return(0)}

#  ) #system.time end
} # function END (survival_marginFree)
##-- -- -- -- --


### spared coded (rJava) installation
# [for updating required packages if R upgraded]
# brew upgrade r
# 
# x install.packages("r2excel") # dependence: ##
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
# 
# brew cask list
# java

#***# install java JDK on ubuntu
# $ apt-get install default-jdk
# then configuration: https://thishosting.rocks/install-java-ubuntu/
# $ date-alternatives --config java
# There is only one alternative in link group java (providing /usr/bin/java): /usr/lib/jvm/java-8-openjdk-amd64/jre/bin/java
#interpreter : '/usr/bin/java'
#archiver    : '/usr/bin/jar'
#compiler    : '/usr/bin/javac'
#header prep.: '/usr/bin/javah'


# Low-level R to Java interface
# # install.packages("rJava") # not available in R 3.4.2
# # # or
# # $ wget https://www.rforge.net/rJava/snapshot/rJava_0.9-9.tar.gz # 2018/01
# # $ wget https://www.rforge.net/rJava/snapshot/rJava_0.9-10.tar.gz # 2018/04/22
# # $ curl -O https://www.rforge.net/rJava/snapshot/rJava_0.9-9.tar.gz
# $ ~$ sudo R CMD javareconf # JNI issue
# # $ R CMD INSTALL rJava_0.9-10.tar.gz # [install rJava from source]
# >install.packages("rJava", repos="http://rforge.net", type="source") # checking JNI data types... ok
# or
# $ Rscript -e 'install.packages("rJava", repos="http://rforge.net", type="source")'
# # ** testing if installed package can be loaded
# however, libjvm.so, java.so can't be opened=>
#(ok) https://stackoverflow.com/questions/28462302/libjvm-so-cannot-open-shared-object-file-no-such-file-or-directory
# answered Nov 8 '15 at 10:44, by minhas23
#1) sudo rstudio-server stop
#2) export LD_LIBRARY_PATH=/usr/lib/jvm/jre/lib/amd64:/usr/lib/jvm/jre/lib/amd64/default
#3) sudo rstudio-server start
## library(rJava) # * DONE (rJava) on GCP(Linux) and macOS(MBA)


# https://docs.oracle.com/javase/7/docs/webnotes/install/mac/mac-jdk.html
# http://www.oracle.com/technetwork/java/javase/downloads/index.html # JDK 7 8 9
# https://github.com/MTFA/CohortEx/wiki/Run-rJava-with-RStudio-under-OSX-10.10,-10.11-(El-Capitan)-or-10.12-(Sierra)
# https://stackoverflow.com/questions/26755013/install-xlsx-and-rjava-on-macos-mavericks-10-9-5/28886808#28886808
# error in MBP (Apple LLVM)
# If you do plan using RStudio, you can use the original Mac Java 1.7
# (In Xcode, the LLVM compiler uses the Clang front end (a C-based languages project on LLVM.org) to parse source code and turn it into an interim format. )
# >clang -o libjri.jnilib Rengine.o jri.o Rcallbacks.o Rinit.o globals.o rjava.o  -dynamiclib -framework JavaVM -F/opt/local/Library/Frameworks/R.framework/.. -framework R -llzma -lm -liconv -licuuc -licui18n
# ld: library not found for -licuuc
# clang: error: linker command failed with exit code 1 (use -v to see invocation)
# => https://github.com/Homebrew/homebrew-science/issues/5968
## [2018/03/15] good solution:
# sudo vim /opt/local/Library/Frameworks/R.framework/Resources/etc/Makeconf and change the line
# => LIBS =  -llzma -lm -liconv 
# # -licuuc -licui18n # removal of -licuuc
# # install the original 1.9 Mac Java
# >sudo curl -#ROL -b "oraclelicense=a" http://download.oracle.com/otn-pub/java/jdk/9.0.4+11/c2514751926b4512b076cc82f959763f/jdk-9.0.4_osx-x64_bin.dmg
# >open jdk-9.0.4_osx-x64_bin.dmg
# >export JAVA_HOME=/Library/Java/JavaVirtualMachines/jdk-9.0.4.jdk/Contents/Home
# X >export JAVA_HOME=/System/Library/Java/JavaVirtualMachines/1.7.0.jdk/Contents/Home
# to tell R to use our Java 9 as it's JAVA_HOME
# >sudo R CMD javareconf
# 
# error(OK) (conftest.c:1:10: fatal error: 'jni.h' file not found)
#  => brew reinstall R --with-openblas --with-java
#  
#  >sudo R CMD javareconf \ JAVA_CPPFLAGS='-I/System/Library/Frameworks/JavaVM.framework/Headers -I/Library/Java/JavaVirtualMachines/jdk-9.0.4.jdk/'
#  
# >sudo Rscript -e 'install.packages("rJava", repos="http://rforge.net", type="source")'


