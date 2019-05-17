survival_marginFree <- function(geneName) {
#  system.time( # as a timer
  #save(ZSWIM2, file=desti_ZSWIM2) # save it globally
  #global decalration on main.R
  #
  # for debug
  # sink(stdout(), type="message")
  # options(error=traceback, showWarnCalls=T)
  
# delare as a function to be called
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
# # DO NOT RUN this {https://www.r-bloggers.com/upgrade-r-without-losing-your-packages/
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



# https://github.com/IARCbioinfo/awesome-TCGA # list of all usefull TCGA tools
# or https://gdc.cancer.gov/access-data/gdc-community-tools, such as GDCtools
# FirebrowseR - Paper describing the R FirebrowseR package.
# GenomicDataCommons - Paper describing the R GenomicDataCommons package.
# 
## [1.Directly] start from beginning  get Broad Institute GDAC: TCGA/Firhose data into R (Retrieve TCGA CDEs verbatim) #
# FirebrowseR - An R package to download directly the results of the analyses performed by Firehose in R.
#go to FireBrowse ( http://gdac.broadinstitute.org/ ):
# install.packages("devtools")
# library("devtools")
# devtools::install_github("mariodeng/FirebrowseR") # with more features (81): such as residual_tumor, vital_status, days_to_last_followup, "smoking duration"

# library(FirebrowseR)
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




# [2018/03/06] with Peter Lai birthday
# 
# [Option-Cmd + E] for running to the end
# ***geneName<-"SLC2A9" # for debug
# Resume:[5a.START] ZZZ3 and TMSB4X ####
# Start the survival analysis for each individual gene
## # set path on google drive
#library(FirebrowseR)
path_cohort <- "~/R/HNSCC_Tex_survival/hnscc_github"
TCGA_cohort <- "HNSCC" # cancer type
# path_LUAD <- "~/R/LUAD_Peter_survival" # under rstudio-server on GCP
#path_LUAD <- "~/R/LUAD_Peter_survival/mount" # mount google drive to ubuntu@instances-4 at GCP
#path_LUAD <- "/Users/apple/Google\ Drive/2018_PhD_dissertation/LUAD_Peter_survival/"
# path_LUAD <- "/Users/apple/Documents/My\ Tableau\ Repository/Workbooks/"
setwd(path_cohort) # set the working directory to the google drive => GCP


#load(file="TMU_TA51BCDE_T70_clinical_fullFeatures13_dichotomized.Rda") # as oscc (without margin feature)
# or
load(file="LUAD_T522_clinical_fullFeatures11_dichotomized.Rda") # as oscc with "margin", negative n=348 while positive n=18 and 156:NaN)
#
#load(file="whole_genome.Rda") # the name list of protein coding genome: whole_genome
# x LUAD_n <- length(whole_genome) # we do not need this anymore


# # debug; run inside the function()
# geneName <- whole_genome[20483] # aka. ZSWIM3
# # debug
#
# geneName <- whole_genome[LUAD_n] # aka. "ZZZ3"
# print(paste("Gene name = ", geneName, sep=""))

# load and prepare HNSCC.clinico_mRNA.Fire for survival analysis
# check file permissions: modes 4 test for read permission; 2 write; 1 excutable; 0 exist.
#if (file.access(".", 4))
# 
if (is.na(tryCatch(load(file=paste("LUAD.mRNA.Exp.", geneName, ".Fire.Rda", sep="")), error = function(e) return(NA)))) {return(3)} # there is NO such file
                   # as LUAD.mRNA.Exp.Fire
# there is NO such file: "LUAD.mRNA.Exp.ZFP91.CNTF.Fire.Rda" %in% dir() => False
# 
# "sample_type"          # 59 NT(Solid Tissue Normal) 515 TP(Primary Solid Tumor) or 2 TR(Recurrent Solid Tumor)

LUAD_T_mRNA.Fire <- LUAD.mRNA.Exp.Fire[LUAD.mRNA.Exp.Fire$sample_type %in% c("TP"), c("tcga_participant_barcode","z.score")] # n=517 -2 = 515
# removing duplicated ID %in%  c("TCGA-50-5066","TCGA-50-5946") => TR recurrent sample

# inner join by merge
if (is.na(tryCatch(LUAD.clinico_mRNA.Fire <- merge(oscc, LUAD_T_mRNA.Fire, by="tcga_participant_barcode") , error = function(e) return(NA)))) {return(5)} # merge error
                   #n=515 with sample_type "TP", excluding "NT normal tissue" or "TR"

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
# n=? in HNSCC
oscc <- LUAD.clinico_mRNA.Fire # starting analysis with "oscc" HNSCC
# oscc$H.score_T as LUAD.mRNA.Exp.Fire$z.score; expression level: H.score_T as RNAseq z.score

## a dummy universal variable for binomial (high/low) oscc$geneName_median, all are zero
oscc$PMM1_median <-(oscc$H.score_T >= median(oscc$H.score_T, na.rm=T)) +0 # higher 1 or lower 0
osccM_pos <- which(colnames(oscc) == "PMM1_median") # at column 14

#x oscc$geneM <- 0
#x colnames(oscc)[colnames(oscc) == "geneM"] <- paste(geneName, "_median", sep="") # rename as oscc$geneName_median
#x osccM_pos <- which(colnames(oscc) == paste(geneName, "_median", sep="")) # at column 19 :-)




##
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

# *** addNA for counting all NA (e.g. there is 0, no 1, in margin-free cohort) in "M" "stage_2" or "margin"
# the “pathological” case of two different kinds of NAs which are treated differently: exclude = if (useNA == "no") c(NA, NaN)
# "unusual NA comes from addNA() as factor
# [2018/03/14] finallized debug
osccCleanNA$margin <- addNA(osccCleanNA$margin, ifany=F) # ifany=="always" add NA as a factor
#   0    1 <NA> => 3 levels of factor in "margin"
# 245   11    0 
#  osccCleanNA$M <- addNA(osccCleanNA$M, ifany=F) # "always" add NA as a factor
#  osccCleanNA$stage_2 <- addNA(osccCleanNA$stage_2, ifany=F) # "always" add NA as a factor

oscc_n256 <- osccCleanNA # n=256; removal of NA cases
#osccClean1 <- osccClean # original cohort n=515

#+++ end of data cleaning ++++
#
# {
# surgical margin status: keeping 0 and excluding + margin (as 1; n=11)
osccCleanNA_freeMargin <- osccCleanNA[osccCleanNA$margin == 0, ] # margin==0
#osccCleanNA <- osccCleanNA_freeMargin # n=245, LUAD s/p OP with margin free##
# n=11, margin involved; 11/256 = 5.3% (how about it's survival impact on each individual genes in LUAD?)
# 
# option1) marginFree
# margin free cohort (n=245): #### 
osccCleanNA <- osccCleanNA_freeMargin
# 
# option2) marginS
# margin positive and negative cohort (n=256) #### 
#osccCleanNA <- oscc_n256
#}


## 5b. (repeat100) Cutoff finder [osccCleanNA] ####
oscc <- oscc0 <- osccCleanNA # syncrhonize them

# checking the completeness
which(complete.cases(oscc$H.score_T)==F) # no NaN -> 0
which(complete.cases(oscc$OS_IND)==F) # no Nan -> 0

# *** column 9 should be OS_IND
which(complete.cases(oscc[oscc$OS_IND==1,9])==F) #OS_IND ==1, death event (dead) => no NaN



# # (skipped, if LUAD RFS is copied from OS)
# which(complete.cases(oscc$RFS_IND)==F) # RFS_IND
# # which(complete.cases(oscc$RFS..months._from.op)==F) #n=103; it may be imputed from OS time (oscc$OS..months._from.biopsy)
# osccCleanNA_RFS <- oscc[which(complete.cases(oscc$RFS..months._from.op)==F), ]
# osccCleanNA_RFS$RFS..months._from.op <- osccCleanNA_RFS[osccCleanNA_RFS$RFS_IND==1, 9] # imputed from OS..months._from.biopsy
# 




# start 100 round ####
# oscc <- cbind(oscc0, osccCleanNA$H.score_T) # expression(IHC score) of PMM1
# debug
# 100 slicing is NOT right !
#find_repeat <- 100 # searching 100 slices in the interval
cohort_n <- nrow(oscc)
# debug
oscc0_pos <- which(colnames(oscc0) == "H.score_T") # oscc0$H.score_T as colunm 13


exp_geneName <- t(oscc0[, oscc0_pos]) # horizontal vector of RNAseq, why t?
#exp_geneName <- oscc0[, oscc0_pos] # vertical vector
# $x is RNAseq be sorted, $x is its original position
# sorting issue ##
exp_geneName_sorted <- data.frame(sort(exp_geneName, decreasing = F, method="radix", index.return=T))
#asc <- order(exp_geneName, method="radix")


p_OS <-  data.frame(matrix(data = NA, nrow = cohort_n, ncol = 1))
cases_OS <-  data.frame(matrix(data = NA, nrow = cohort_n, ncol = 1))
p_RFS <-  data.frame(matrix(data = NA, nrow = cohort_n, ncol = 1))
cases_RFS <-  data.frame(matrix(data = NA, nrow = cohort_n, ncol = 1))
#run100 <- 0
#cut_featuresRemark <- 0 # finder results

# oscc_pos <- which(colnames(oscc) == geneName)
# oscc[oscc_pos] <- exp_geneName
#if (is.na(tryCatch(cutoff <- quantile(exp_geneName, c(0.30,0.70)), error = function(e) return(NA)))) {return(4)} # return to main by "XKRY"(19642)
#xxx cutoff <- quantile(exp_geneName, c(0.30,0.70)) # by RNAseq value
# by cases
cutoff_n <- round(quantile(c(1:cohort_n), c(0.30,0.70))) # 30 percentile, 70 percentile
cutoff_n[2] <- cutoff_n[2] -1 # 78~179 in cases 256
# *debug, error and stop on "XKRY"(19642)#
# (ok)Error in quantile.default(exp_geneName, c(0.3, 0.7)) : 
#   missing values and NaN's not allowed if 'na.rm' is FALSE
# (OK)Error during wrapup: names() applied to a non-vector
#debug

# if (!require(pkg)){ 
#   install.packages(pkg) 
# } # Install package automatically if not there
# 



# https://cran.r-project.org/web/packages/survival/survival.pdf
#for (i in seq(cutoff[1], cutoff[2], length.out = find_repeat)){
# P-value according to KM survival analysis (alone)
for (run100 in seq(cutoff_n[1], cutoff_n[2])){ 
  # sorted (by RNAseq) no.78~179 in cases 256 of LUAD
  # use oscc0 for "positioning"; 
  # oscc is the dataset to be analysed here.
  
  # Binominal H.Score_T (RNAseq) by exp_geneName_sorted$x[i]) #RNAseq cutoff
  # exp_geneName_sorted$x[i], cutoff100 is the cutoff value of RNAseq in rank i
  cutoff100 <- exp_geneName_sorted$x[run100]
  oscc[osccM_pos] <- (oscc0[oscc0_pos] >= cutoff100) +0 # oscc$PMM1_median <- (oscc0$PMM1 >= i) +0
  #run100 <- i - as.numeric(cutoff_n[1]) +1 # *** incremental (start from 1) by i
#print(paste(run100, cutoff100))
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
  # 
  # run survival analysis: OS
  mysurv <- Surv(oscc$OS..months._from.biopsy, oscc$OS_IND==1) #1==death event
  # Test for difference (log-rank test) with P-value
  if (is.na(tryCatch(surv_OS <- survdiff(mysurv ~ as.vector(unlist(oscc[osccM_pos]), mode="numeric"), data=oscc), error = function(e) return(NA)))) {return(2)} # grouping by PMM1_median
  # OS, error 2 due to ZSWIM2? (one group only) ####
  # extract P-value from surv_OS, put in "original" position (unsorted)
  p_OS[exp_geneName_sorted$ix[run100], 1] <- format(pchisq(surv_OS$chisq, length(surv_OS$n)-1, lower.tail = FALSE), digits=3)
#  cases_OS[j] <- surv_OS$n[1] #cutoffs by cases, remaping sorted ####
  #[coxph] - fits a Cox proportional hazards regression model
  #OS.km <- survfit(mysurv ~ oscc[osccM_pos], data=oscc, conf.type = "log-log")
  #  jpeg(file=paste("KMplot_OS_", geneName, i, ".jpg", sep = ""))
  #  plot(OS.km, lty=1, col=c("blue","red"), sub=paste("p-value =", p_OS[j]), main=paste("OS in OSCC(n=", surv_OS$n[1]+surv_OS$n[2],")/", geneName, "cutoff=", i), ylab="Percent Survival", xlab="Days")
  #  legend("topright", legend=c(paste("low(",surv_OS$n[1], ")"), paste("high(",surv_OS$n[2], ")")), lty=1:2, col=c("blue","red"))
  #  dev.off()
    
  
  # run survival analysis: RFS
  mysurv <- Surv(oscc$RFS..months._from.op, oscc$RFS_IND==1) #1==tumour recurrency
  # Test for difference (log-rank test) with P-value
  # RFS, error 2 = function(e) return(NA)) for "Survdiff.fit error on there is only one group"   # due to ZSWIM2 (one group only) ####
  if (is.na(tryCatch(surv_RFS <- survdiff(mysurv ~ as.vector(unlist(oscc[osccM_pos]), mode="numeric"), data=oscc), error = function(e) return(NA)))) {return(2)}

  p_RFS[exp_geneName_sorted$ix[run100], 1] <- format(pchisq(surv_RFS$chisq, length(surv_RFS$n)-1, lower.tail = FALSE), digits=3)
#  cases_RFS[j] <- surv_RFS$n[1] #cutoffs
  #RFS.km <- survfit(mysurv ~ oscc[osccM_pos], data=oscc, conf.type = "log-log")
  #  jpeg(file=paste("KMplot_RFS_", geneName, i, ".jpg", sep = ""))
  #  plot(RFS.km, lty=1, col=c("blue","red"), sub=paste("p-value =", p_RFS[j]), main=paste("RFS in OSCC(n=", surv_RFS$n[1]+surv_RFS$n[2],")/", geneName, "cutoff=", i), ylab="Percent Survival", xlab="Days")
  #  legend("topright", legend=c(paste("low(",surv_RFS$n[1], ")"), paste("high(",surv_RFS$n[2], ")")), lty=1:2, col=c("blue","red"))
  #  dev.off()
  
  #  ++++ Beside KM P-value, contingency P-value is also important ! ++++++
  # Calling for contingency table to find significant clinicopathological features
  # contingency P-value
  # **debug(contingencyTCGA) ####
  contiT <- contingencyTCGA(oscc, geneName) # calling this function (OSCC cohort, PMM1) 2 parameters, no more cutoff1; 
  # There were 12 warnings (use warnings() to see them)
  # exp_geneName_sorted$x[i] is cutoff value of RNAseq
  if (contiT[[1]]=="skip") {return(1)} # ***ZSWIM2 skip "error prone" genes, error 1
  #(OK)Error in if (contiT[[1]] == "skip") :NA, missing value where TRUE/FALSE needed

  # "margin" at column 8 of oscc
    #  oscc <- contiT[[1]] # updating PMM1_median
  chiT <- contiT[[2]] # extrac it from list by using [[]]; chiT$X2 is the P-value
  freq <- contiT[[3]] # well DONE
  
  # no more table 2 here
  
} # end of 100 cut/slice, defined by [find_repeat]
# cases_OS[, 1] <- exp_geneName_sorted$ix[seq(cutoff_n[1],cutoff_n[2])] #remaping the index after sorting
# cases_RFS[, 1] <-exp_geneName_sorted$ix[seq(cutoff_n[1],cutoff_n[2])]
cases_OS[, 1] <- exp_geneName_sorted$ix # unsorted rank
cases_RFS[, 1] <-exp_geneName_sorted$ix
print(paste("Cutoff finding is done, for", geneName, "(", which(whole_genome==geneName), ")"))
# no more use of interim surv_OS, surv_RFS







# Generate OS and RFS "table" with cutoffs and its P-value in KM ####
#{
OS <- data.frame(cases_OS, as.numeric(unlist(p_OS)), c(1:cohort_n), as.numeric(unlist(exp_geneName)))
# removal of duplicated item
OS <- OS[!duplicated(OS),]
colnames(OS) <- c("cases_OS","p_OS", "rank", "exp")
OS <- OS[complete.cases(OS$p_OS), ] # removal of NA, keep n=102

RFS <- data.frame(cases_RFS, as.numeric(unlist(p_RFS)), c(1:cohort_n),  as.numeric(unlist(exp_geneName)))
# removal of duplicated item and NA
RFS <- RFS[!duplicated(RFS),]
colnames(RFS) <- c("cases_RFS","p_RFS", "rank", "exp")
RFS <- RFS[complete.cases(RFS$p_RFS), ] # removal of NA

#}

# # rownames(cut_featuresRemark) <- tableChi1$Features # there is duplicated row names :-)
# # summary(cut_featuresRemark) to see result
# 
# # to generate the 201 first odd numbers; c(T,F) or seq.int(1, by = 2, length.out = find_repeat + 1)
# # Sequence Generation: seq() or seq.int()
# # https://www.rdocumentation.org/packages/base/versions/3.4.1/topics/seq
# 
# #describe(cut_featuresRemark[, c(1:(find_repeat * 2 + 1))[c(T,F)]])
# 

library(graphics) # plot HOLD ####
# #{
# cumulative P value curves
library(ggplot2) # http://www.cookbook-r.com/Graphs/Legends_(ggplot2)/

# ggplot(OS, aes(x=cases_OS, y=p_OS)) + geom_point(size=2) +
#   scale_x_discrete(limit = seq(0,500,50), labels = paste(seq(0,500,50), sep=",")) +
#   scale_y_discrete(breaks = c(0, 0.05, 0.10, 0.50, 0.99)) + #, labels = c("0", "0.05", "0.10", "0.50", "0.99")) +
#   coord_trans(y = "log10") +
#   ggtitle(paste("Cumulative P value curves for OS under", geneName, "expression")) +
#   xlab("# of patients") + ylab("P-value(log10)") +
#   geom_hline(yintercept=0.05, linetype="dashed", color = "red", size=1, show.legend = T)
# #  geom_vline(xintercept=300, color = "green", size=2, show.legend = T)
# View(OS[OS$p_OS <= 0.05,1:2])
# #
# #plot(cases_OS, p_OS, type="p", log="y", main=paste("Cumulative P value curves for", geneName), xlab="# of patients", ylab="P-value")
# #axis(side=4, at=c(0.01, 0.05, 0.2, 0.5, 0.99))
# #abline(h=0.05, lty=2, col="red")
# #abline(v=320, col="green")
# #lines(x=c(cases_OS[1], cases_OS[length(cases_OS)]), y=c(0.05, 0.05), col = "red")
# 
# plot(cases_RFS, p_RFS, type="p")
library(graphics)
# cumulative P value curves
library(ggplot2)

# ggplot(RFS, aes(x=cases_RFS, y=p_RFS)) + geom_point(size=2) +
#   scale_x_discrete(limit = seq(0,500,50), labels = paste(seq(0,500,50), sep=",")) +
#   scale_y_discrete(breaks = c(0, 0.05, 0.10, 0.50, 0.99)) + #, labels = c("0", "0.05", "0.10", "0.50", "0.99")) +
#   coord_trans(y = "log10") +
#   ggtitle(paste("Cumulative P value curves for RFS under", geneName, "expression")) +
#   xlab("# of patients") + ylab("P-value(log10)") +
#   geom_hline(yintercept=0.05, linetype="dashed", color = "red", size=1, show.legend = T)
# geom_vline(xintercept=150, color = "green", size=2, show.legend = T)
# View(RFS[RFS$p_RFS <= 0.05,1:2])
# #}
# describe(cut_featuresRemark$Remark)

# [done cutoff finder] ++++++++++++++++++++++++
#{ PAUSE here.....
#
#
#
#
#
#

# {5c.Auto}: select a cutoff=> please choice ONE from OS (in LUAD) or RFS ####
# auto pickup one cutoff from OS/RFS tables, which has lowest P-value, to run directly then export
OS_pvalue <- OS[OS$p_OS <= 0.05, ] # with rank registration
cohort_n <- nrow(oscc) # length(oscc) is columne number (variable number)
#cohort_n <- (surv_OS$n[1]+surv_OS$n[2])
# case50_n is the cutoff1 by cases number in unsorted.
if (nrow(OS_pvalue) > 0) {
  case50_n <- OS_pvalue[which.min(OS_pvalue$p_OS), 3] # there is a hit; original rank
#  cutoff1 <- exp_geneName[case50_n] # unsorted, original rank####
  cutoff1 <- OS_pvalue[which.min(OS_pvalue$p_OS), 4]
  #cohort_rank <- OS_pvalue[which(OS_pvalue$cases_OS==case50_n), 3] # for KM plot legend
  } else {
    # sad, no hit, using default 50% cutoff in SORTED manner at n=128 around
    # [if multiple candidate, it will pickup the first hit]
    default50 <- 0.5
    case50_n <- cohort_n * default50
    cutoff1 <- quantile(exp_geneName_sorted$x, c(default50))  # ?half of missingless LUAD cases
    } # cutoff1 submit for statistics export


# (OK)Error in quantile.default(exp_geneName, c(case50_n/cohort_n)) : 
# 'probs' outside [0,1] => cohort_n
# or 
##cutoff1 <- quantile(exp_geneName, c(case50_n/(surv_RFS$n[1]+surv_RFS$n[2])))  # 56 out of missingless LUAD cases

# for PMM1; cutoff at case n=56 (OS) or case n=? (RFS),  at lowerest P-value
# Kaplan-Meier survival of OS, P-value = 0.0198; of RFS, P-value = ? in LUAD



# 5d.[Statistic procedure] [Option-Cmd + E] for running to the end ###
# The features from column L to R; running to the end of R2Excel export +++++++++++
# The correlation of gene expression and clinical features  ### contingency tables
# if (!require(pkg)){ 
#   install.packages(pkg) 
# } # Install package automatically if not there



# # x debug {contingencyTCGA)
# # # call function for chiT, freq data.frame
# contiT <- contingencyTCGA(osccCleanNA, geneName, cutoff1) # calling this function (OSCC cohort, PMM1, cutoff1)
# osccCleanNA <- contiT[[1]] # updated
# chiT <- contiT[[2]] # extrac it from list by using [[]]
# freq <- contiT[[3]] # well DONE
# View(chiT) # for p-value (=0.0389 in pathologic_T)
# View(freq) # contingency table
# # } # debug



# Statistics for osccCleanNA by cutoff1
## (6) KM survival curves in a specific cohort ####
# right censored? yes
# Kaplan-Meier curve + Log-rank test, + Cox proportional regression
# textbook: David Kleinbaum???Survival analysis: A self-learning text
library(survival)
library(rms)

osccCleanNA <- oscc # synchronize; be sorted (osccCleanNA) by RNAseq ##
#osccCleanNA1 <- osccCleanNA
osccCleanNAM_pos <- which(colnames(osccCleanNA) == paste("PMM1", "_median", sep=""))
osccCleanNA_pos <- which(colnames(osccCleanNA) == "H.score_T")


# Refresh binominal by a definite cutoff value by AUTO selection: cutoff1
# [**binomial of PMM1_median] resume the correlation table 2 ####
osccCleanNA[osccCleanNAM_pos] <- (osccCleanNA[, osccCleanNA_pos] >= cutoff1) +0 # binomial after osccCleanNA 
# vertical vector is corrected

# osccCleanNA[osccCleanNAM_pos] <- (osccCleanNA[osccCleanNA_pos] >= cutoff1) +0
# it is also ok # horizontal vector
#all.equal(as.numeric(osccCleanNA$H.score_T), as.numeric(exp_geneName), check.attributes = F)
#> TRUE

# 3650 days = 10 years; cut it into 3650 days (do not need to remove 6 cases)
# over3650 <- osccCleanNA[osccCleanNA$X_RFS > 3650, 1]
# [1] TCGA-CV-5432-01 TCGA-CV-7183-01 TCGA-CV-7435-01 TCGA-CV-A45Q-01 TCGA-CV-A45R-01
# [6] TCGA-CV-A45T-01
# TCGA-CV-5430-01; removal one more
#osccCleanNA <- osccCleanNA[!osccCleanNA$sampleID %in% over3650,] # n=366
#osccCleanNA <- osccCleanNA[!osccCleanNA$sampleID == "TCGA-CV-5430-01",] # n=356


# KM survival curves
#OS1
# cancer type shouble be defined at TCGA_cohort <- "LUAD" "HNSC"
mysurv <- Surv(osccCleanNA$OS..months._from.biopsy, osccCleanNA$OS_IND==1) #1==dead
# Test for difference (log-rank test)
tryCatch(surv_OS1 <- survdiff(mysurv ~ as.vector(osccCleanNA[, osccCleanNAM_pos], mode="numeric"), data=osccCleanNA), error = function(e) return(NA)) # PMM1 high or low
# pchisq gives the distribution function
p_OS1 <- format(pchisq(surv_OS1$chisq, length(surv_OS1$n)-1, lower.tail = FALSE), digits=3)
#cases_OS1 <- surv_OS1$n[1]
# #Value of survdiff => a list with components: help("survdiff")
# n => the number of subjects in each group.
# obs
# the weighted observed number of events in each group. If there are strata, this will be a matrix with one column per stratum.
# exp
# the weighted expected number of events in each group. If there are strata, this will be a matrix with one column per stratum.
# chisq
# the chisquare statistic for a test of equality.
# var
# the variance matrix of the test.
# strata
# optionally, the number of subjects contained in each stratum.

#

# [survfit] - Kaplan-Meier curve

# confidence intervals as log hazard or log(-log(survival))
OS.km <- survfit(mysurv ~ as.vector(unlist(osccCleanNA[, osccCleanNAM_pos]), mode="numeric"), data=osccCleanNA, type= "kaplan-meier", conf.type = "log-log")
# 365.25 days in a year for xscale => 3650 days for 10 years
# 12 months per year for 5 years => 60 months

plot(OS.km, lty=1, xscale=12, xmax=60, col=c("blue","red"), sub=paste("P-value =", p_OS1), main=paste("OS in TCGA", TCGA_cohort, "(n=", surv_OS1$n[1]+surv_OS1$n[2],")/", geneName), ylab="Percent Survival", xlab="Years")
legend("topright", legend=c(paste("low(",surv_OS1$n[1], ")"), paste("high(",surv_OS1$n[2], ")")), lty=1:2, col=c("blue","red"))
#plot(OS.km, lty=1, col=c("blue","red"), sub="p= 0.816", main="OS in osccCleanNA(n=505)/gene level", ylab="Percent Survival", xlab="Days")
#legend("topright", legend=c('Low', 'High'), lty=1:2, col=c("blue","red"))
# summary(OS.km, times = seq(0, 3000, 100))




#RFS1
mysurv <- Surv(osccCleanNA$RFS..months._from.op, osccCleanNA$RFS_IND==1) #1==tumour recurrency
# Test for difference (log-rank test)
tryCatch(surv_RFS1 <- survdiff(mysurv ~ as.vector(unlist(osccCleanNA[, osccCleanNAM_pos]), mode="numeric"), data=osccCleanNA), error = function(e) return(NA))
p_RFS1 <- format(pchisq(surv_RFS1$chisq, length(surv_RFS1$n)-1, lower.tail = FALSE), digits=3)
#cases_RFS1 <- surv_RFS1$n[1]

RFS.km <- survfit(mysurv ~ as.vector(unlist(osccCleanNA[, osccCleanNAM_pos]), mode="numeric"), data=osccCleanNA, conf.type = "log-log")
# 12 months per year for 5 years => 60 months

plot(RFS.km, lty=1, xscale=12, xmax=60, col=c("blue","red"), sub=paste("P-value =", p_RFS1), main=paste("RFS in TCGA", TCGA_cohort, "(n=", surv_RFS1$n[1]+surv_RFS1$n[2],")/", geneName), ylab="Percent Survival", xlab="Years")
legend("topright", legend=c(paste("low(",surv_RFS1$n[1], ")"), paste("high(",surv_RFS1$n[2], ")")), lty=1:2, col=c("blue","red"))
#plot(RFS.km, lty=1, col=c("blue","red"), sub="p= 0.0487", main="RFS in osccCleanNA(n=372)/gene level", ylab="Percent Survival", xlab="Days")
#legend("topright", legend=c('Low', 'High'), lty=1:2, col=c("blue","red"))
# summary(RFS.km, times = seq(0, 3000, 100))
# [DONE] ##  back to Cutoff finder (n=185) or set 50:50 on section [5c]





## ## 7)[COXPH modelling] ####
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

# rename all of colnames (features) as osccT
colnames(osccCleanNA) <- coln_osccT # colnames (features) of osccT

osHRmulti <- 0
# 8 features in LUAD
# #warning() X matrix deemed to be singular (margin) in coxph
# https://stackoverflow.com/questions/20977401/coxph-x-matrix-deemed-to-be-singular
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
## # 8 features in LUAD
rfsHRmulti <- 0
# P-value notes:
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
# # 8 features in LUAD, put in coxph one by one
# number of features 13 -> 14 (add a feature: margin status)
rfscox <- 0
features_os <- colnames(osccCleanNA)[c(2:8,14)] # features selection ["Gender" to "margin", "PMM1"] 8 out of 14
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
#TCGA_cohort <- "LUAD" # cancer type
#path_LUAD <- "~/R/LUAD_Peter_survival/" # under rstudio-server on GCP
#path_LUAD <- "/Users/apple/Google\ Drive/2018_PhD_dissertation/LUAD_Peter_survival/"
# path_LUAD <- "/Users/apple/Documents/My\ Tableau\ Repository/Workbooks/"
path_cohort <- "~/R/HNSCC_Tex_survival/hnscc_github"
setwd(path_cohort) # change working directory to the google drive

# filename for all margins cases, eighter _marginFree_ or _marginS_
filenamex <- paste(TCGA_cohort, "_survivalAnalysis_marginFree_", geneName, ".xlsx", sep = "")
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




# Part I: # Table 2 construction by Calling contignecyBin ####
# ??using tableChi2 by contigencyBin2 with c("Gender","ageDx", "pathologic_T", "pathologic_N", "pathologic_M", "stage","margin" )
# under beset cutoff value (auto choosen)
# Table 2: clinical correlation tables (tableChi1)
tableChi1 <- contingencyBin (osccCleanNA, chiT, freq) # calculating the P-value
# ***add rownames
nrow_featuresUni <- length(featuresUni) # aka. 8
nrow_tableChi1 <- as.character(seq(1, 2*(nrow_featuresUni-1))) # aka. 7 x2 NA (a vector)
nrow_tableChi1[c(T,F)] <- featuresUni[-nrow_featuresUni] # skip last one row (z-score)
# > nrow_tableChi1
# [1] "Gender"                 NA                      
# [3] "Age at diagnosis"       NA                      
# [5] "Pathologic T status"    NA                      
# [7] "Pathologic N status"    NA                      
# [9] "Pathologic M status"    NA                      
# [11] "Pathologic Stage"       NA                      
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

xlsx.addPlot.OS<-function( wb, sheet, startRow=NULL,startCol=2,
                        width=480, height=480,... )
{ # OS KM plot
  #library("xlsx")
  png(filename = "plot.png", width = width, height = height,...)
  plot(OS.km, lty=1, xscale=12, xmax=60, col=c("blue","red"), sub=paste("P-value =", p_OS1), main=paste("OS in TCGA ", TCGA_cohort, "(n=", surv_OS1$n[1]+surv_OS1$n[2],")/", geneName), ylab="Percent Survival", xlab="Years")
  legend("topright", legend=c(paste("low(",surv_OS1$n[1], ")"), paste("high(",surv_OS1$n[2], ")")), lty=1:2, col=c("blue","red"))
  
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
xlsx.addPlot.OS(wb, sheet)

#
xlsx.addPlot.RFS<-function( wb, sheet, startRow=NULL,startCol=2,
                           width=480, height=480,... )
{  # RFS KM plot
  #library("xlsx")
  png(filename = "plot.png", width = width, height = height,...)
  # plot fuction here
  plot(RFS.km, lty=1, xscale=12, xmax=60, col=c("blue","red"), sub=paste("p-value =", p_RFS1), main=paste("RFS in TCGA ", TCGA_cohort, "(n=", surv_RFS1$n[1]+surv_RFS1$n[2],")/", geneName), ylab="Percent Survival", xlab="Years")
  legend("topright", legend=c(paste("low(",surv_RFS1$n[1], ")"), paste("high(",surv_RFS1$n[2], ")")), lty=1:2, col=c("blue","red"))
  
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
xlsx.addPlot.RFS(wb, sheet)
#

xlsx.addLineBreak(sheet, 8)

# plot cumulative p-value from cutoff finder
xlsx.addPlot.OSpval<-function( wb, sheet, startRow=NULL,startCol=2,
                            width=480, height=480,... )
{   # OS P-values dots plot
  #library("xlsx")
  png(filename = "plot.png", width = width, height = height,...)
  # plot fuction here
  pOS <- ggplot(OS, aes(x=cases_OS, y=p_OS)) + geom_point(size=2) +
    scale_x_discrete(limit = seq(0,500,50), labels = paste(seq(0,500,50), sep=",")) +
    scale_y_discrete(breaks = c(0, 0.05, 0.10, 0.50, 0.99)) + #, labels = c("0", "0.05", "0.10", "0.50", "0.99")) +
    coord_trans(y = "log10") +
    ggtitle(paste("Cutoff Finder: Cumulative P value curves for OS under ", geneName, " expression \n", "(cutoff at ", surv_OS1$n[1], " )", sep = "")) +
    xlab("# of patients") + ylab("P-value(log10)") +
    geom_hline(yintercept=0.05, linetype="dashed", color = "red", size=1, show.legend = T)
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
xlsx.addPlot.OSpval(wb, sheet)
#
#
#
#
xlsx.addPlot.RFSpval<-function( wb, sheet, startRow=NULL,startCol=2,
                               width=480, height=480,... )
{   # RFS P-values dots plot
  #library("xlsx")
  png(filename = "plot.png", width = width, height = height,...)
  # plot fuction here
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
xlsx.addPlot.RFSpval(wb, sheet)
#


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
xlsx.addTable(wb, sheet, data =OS[OS$p_OS <= 0.05, ], fontSize=12, startCol=8,
              row.names = F, col.names = T)
#View(OS[OS$p_OS <= 0.05,1:2])
xlsx.addTable(wb, sheet, data =RFS[RFS$p_RFS <= 0.05, ], fontSize=12, startCol=9,
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
## save to append an object to an R data file.####
# append.Rda <- function(x, file) {
#   old.objects <- load(file, new.env())
#   save(list = c(old.objects, deparse(substitute(x))), file = file)
# }
#
# # re-do ?in OS
# OS <- data.frame(cases_OS, as.numeric(p_OS))
# # removal of duplicated item
# OS <- OS[!duplicated(OS),]
# colnames(OS)[2] <- "p_OS"
# OS_pvalue <- OS[OS$p_OS <= 0.05,1:2]
# # *** appending survival
# #append.Rda(OS[OS$p_OS <= 0.05,1:2], paste("LUAD_survivalAnalysis_marginS_", geneName, ".Rda"))
# 
# # re-do ? in RFS
# RFS <- data.frame(cases_RFS, as.numeric(p_RFS))
# # removal of duplicated item
# RFS <- RFS[!duplicated(RFS),]
# colnames(RFS)[2] <- "p_RFS"
# RFS_pvalue <- RFS[RFS$p_RFS <= 0.05,1:2]
# # *** appending survival
# #append.Rda(RFS[RFS$p_RFS <= 0.05,1:2], paste("LUAD_survivalAnalysis_marginS_", geneName, ".Rda"))

# Added code for table 2, table 3, table 4 and survival p-value save as .Rda
# save(list = c("tableChi1", "tableOS1", "tableRFS1"), file=paste("LUAD_survivalAnalysis_marginS_", geneName, ".Rda"))

#tryCatch(
  RFS_pvalue <- OS_pvalue #:-) for LUAD only
  save(list = c("tableChi1", "tableOS1", "tableRFS1", "OS_pvalue", "RFS_pvalue"), file=paste("LUAD_survivalAnalysis_marginFree_", geneName, ".Rda", sep=""))
  # save _marginFree_ or _marginS_
  #, error = function(e) return(NA))
#print(paste("Create", paste("LUAD_survivalAnalysis_marginS_", geneName, ".Rda", sep=""), "successfully."))


## (9.) Adding code for google spreadsheet export (with adds-on into a single file)  ####
## 
## 
## 
# Final Return
  if (nrow(OS_pvalue) > 0)  { # we hit one gene with P-value < 0.05 in KM plot
    return(which.min(OS_pvalue$p_OS))} else {return(0)}

#  ) #system.time end
} # function END (survival_marginS)
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


