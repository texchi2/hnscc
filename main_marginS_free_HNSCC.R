# main_marginS_free.R [R script on tex@35.201.169.0:/R/LUAD_Peter_survival/]
# or http://10.140.0.3:8787 (internal IP) through $ ssh -D 2001 tex@35.201.169.0
# ssh has more security mechanism (better than http:// alone)
# 這是最主要的 R code (marginS and marginFree) since [2018/05/20]-[2018/07/08] instance-4
# porting to HNSCC_Tex_survival since [2018/10/24] {3 levels:main, TCGA_HNSCC_[marginS/marginFree], cutofFinder_func.R}

# [2018/10/26] sync: pull and push with remote github (well done); Screens connect from iPad pro (ok)
#  git config remote.origin.url git@github.com:texchi2/hnscc.git  # for pull/push (sync by ssh without id/pw)
# Happy Birthday [2018/11/07] resume from project of hnscc_git
# 

# Tutorial: Survival analysis of TCGA patients integrating gene expression (RNASeq) data
# https://www.biostars.org/p/153013/
# HNSCC or LUAD, there is many web tools: https://github.com/mdozmorov/TCGAsurvival/blob/master/Cancer_DB.Rmd
# # https://www.r-bloggers.com/storing-a-function-in-a-separate-file-in-r/ #Rstudio
# times <- dget("times.R") # function {} embraced the times.R
# survival <- dget("TCGA_LUAD_marginS.R", keep.source = FALSE)
# times(-4:4, 2)
# # dput(“function {}”, file = "times.R", control = c("keepNA", "keepInteger", "showAttributes"))
# # or
# source("fun.R") # multi <- function() {} in fun.R
# mult(-4:4, 2)

# install packages
#{
# make from the source "curl" and its libcurl, compiling under shell
# $ wget https://github.com/curl/curl/releases/download/curl-7_59_0/curl-7.59.0.tar.gz
# $ tar -xzvf curl-7.59.0.tar.gz
# $ cd curl-7.59.0/
# $ ./configure # make  # sudo make install
install.packages(c("git2r", "curl", "httr"), repos="https://cran.r-project.org", type="source")
# write something here; try jjj by :imap jjj <Esc>
install.packages(c("R.utils", "compositions", "openssl"))
install.packages(c("psych", "survival", "reshape", "data.table"))
install.packages(c("scales", "dplyr", "magrittr"))
install.packages(c("plyr")) #ddply()
install.packages("ca")
# followings from github
install.packages("devtools")
library(devtools)
install_github("jsugarelli/debugr") # package requires R >= 3.5.0
# dwatch() debugr_switchOff() debugr_switchOn()
devtools::install_github("hoxo-m/pforeach")
devtools::install_github("xvrdm/ggrough") # ggrough converts my ggplot2 plots to rough/sketchy charts, using the excellent javascript roughjs library.
# dependencies gdtools, ‘svglite’, ‘xml2’; https://xvrdm.github.io/ggrough/
# usage:
# ggplot() -> p 
# library(ggrough)  
# options <- list(
#   Background=list(roughness=8),
#   GeomCol=list(fill_style="zigzag", angle_noise=0.5, fill_weight=2))
# get_rough_chart(p, options)
#{
# $ sudo apt-get install libcairo2-dev
#system("sudo apt-get install libcairo2-dev")
devtools::install_github('davidgohel/gdtools')
devtools::install_github("r-lib/svglite")
#$ sudo apt-get install libxml2
devtools::install_github("r-lib/xml2")
# }
devtools::install_github("ismayc/rticles") # R Markdown "Reed Senior Thesis" template, https://www.r-bloggers.com/r-markdown-senior-thesis-template/
# https://www.r-bloggers.com/r-markdown-senior-thesis-template/

install.packages(c("gmailr", "graphics", "ggplot2", "rms", "xlsx", "r2excel", "tis"))
# $ wget https://www.rforge.net/rJava/snapshot/rJava_0.9-10.tar.gz
# install.packages("rJava", repos="~/R/rJava_0.9-10.tar.gz") # it is tuft to install

install.packages(c("xlsx")) # after rJava OK...
# installing curl directly then r2excel
install_github("kassambara/r2excel", force=T) # the only way OK; reinstall it
# or
install.packages("githubinstall")
library("githubinstall")
githubinstall("r2excel")
#} all libraries are installed.

install.packages("minpack.lm")
install.packages("gplots") # also venn
install.packages("VennDiagram")
install.packages("venn")
install.packages("binaryLogic") # as.binary

# for sweave and knitr
install.packages("knitr")


## START: set path on google drive ####
#library(FirebrowseR)
TCGA_cohort <- "LUAD" # cancer type
#path_LUAD <- "/Users/apple/Google\ Drive/2018_PhD_dissertation/LUAD_Peter_survival/"
path_LUAD <- "~/R/LUAD_Peter_survival" # under rstudio-server on GCP
#path_LUAD <- "~/R/LUAD_Peter_survival/mount" # mount google drive to ubuntu@instances-4 at GCP
# it can’t read “mounted” directory from R
#path_LUAD <- "/Users/apple/Documents/My\ Tableau\ Repository/Workbooks/"
setwd(path_LUAD) # set the working directory to the google drive

load(file="whole_genome.Rda") # the name list of protein coding genome
# below 2 lines: it needs to be modified for a function cal
LUAD_n <- length(whole_genome) #last one 20499: "ZZZ3"

# # source our functions; automatically source all of the functions in a directory, "/tmp", 
# # which makes it easy to run a long script with a single run:
# code.dir <- "/tmp"
# code.files = dir(code.dir, pattern = "[.r]")
# for (file in code.files){
#   source(file = file.path(code.dir,file))
# }

# by "Source on Save" checked
#source(paste(path_LUAD, "TCGA_LUAD_marginS.R", sep="")) # survival_marginS <- function() {} in TCGA_LUAD_marginS.R
#source(paste(path_LUAD, "TCGA_LUAD_marginFree.R", sep="")) # survival_marginFree <- function() {} in TCGA_LUAD_marginFree.R
source(file=file.path(path_LUAD, "TCGA_LUAD_marginS.R")) # survival_marginS <- function() {} in TCGA_LUAD_marginS.R
source(file=file.path(path_LUAD, "TCGA_LUAD_marginFree.R")) # survival_marginFree <- function() {} in TCGA_LUAD_marginFree.R
source(file=file.path(path_LUAD, "cutofFinder_func.R")) # cutofFinder_func <- function(geneName) {} in cutofFinder_func.R
# https://support.rstudio.com/hc/en-us/articles/200484448-Editing-and-Executing-Code

#install.packages("gmailr") # dependencies ‘curl’, ‘openssl’ are not available for package ‘httr’
#  'libcurl'
# make from the source "curl" and its libcurl 
# $ wget https://github.com/curl/curl/releases/download/curl-7_59_0/curl-7.59.0.tar.gz
# $ tar -xzvf curl-7.59.0.tar.gz 
# $ cd curl-7.59.0/
# $ ./configure # make  # sudo make install
# then install curl (done)
# > install.packages("curl", repos="https://cran.r-project.org", type="source") 

# "openssl"(done), "httr"(done), "git2r"(done)
# install.packages("git2r", repos="https://cran.r-project.org", type="source")
#install.packages("devtools")
#devtools::install_github("hoxo-m/pforeach")
library(gmailr) # notify me the script progress by email...https://developers.google.com/gmail/api/
library(scales)
# library(doParallel) # core n=2 in my CPU
# library(foreach)
# #library(pforeach)
library(plyr); library(dplyr) #ddply()
library(magrittr)
recipient <- "texchi2@gmail.com"
sender <- "texchi2@gmail.com"


### $setup global variables...(including "margin") ####
library("psych") # for describe()
library(survival)

# ***colnames of osccT
coln_osccT <- c("Unique.ID","Gender","ageDx",
                "pathologic_T","pathologic_N",
                "pathologic_M","stage", "margin",
                "OS_IND","OS..months._from.biopsy",
                "RFS_IND", "RFS..months._from.op", "H.score_T", paste("PMM1", "_median", sep=""))


# for tableChi1 (table 2)
contLowHighN <- c("Features",	"Low", "(%)",	"High",	"(%)", "Case no",	"P-value", "Remark", "Matthews")

# Create the Sheet title and subtitle; underline as 0 (default value, no underline), 1 (underline with one line), 2 (underline with two lines)
# Table 3. Univariate/Multivariate Cox proportional hazards regression analyses on OS time
ciUniMti <- c("Features",	"HR",	"CI95%(L)",	"CI95%(H)",	"P-value",	"HR",	"CI95%(L)",	"CI95%(H)",	"P-value")

# with 8 features: in table 3 and table 4
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
                 paste("RNAseq(z-score)") # "IHC score" == "z-score"
)
# z-score calculation = [(value gene X in tumor Y)-(mean gene X in normal)]/(standard deviation X in normal)

# lower rowname of table 3 and table 4
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
# upper rowname of table 3 and table 4
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

####
# Define global functions
# contingency function, Tex's design (a confusion matrix with TP TN FP FN)



# # #x for DEBUG only {
# osccCleanNA <- oscc # T,N,M and stage_2
# # osccCleanNA has n=245 in LUAD (with margin free) or n=256 in LUAD (margin 0 or 1)
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




## function declare
# create correlation table 1 with P-value by chisq.test
library(reshape)
library(data.table)
#library(ca) # for Simple correspondence analysis
contingencyTCGA <- function(osccCleanNA_conTCGA, geneName) { # no more "run100"; do not need cutoff1 here
  
  #!!!!!# L <- 2 ("Gender"); R <- 8 ("margin") # in LUAD
  ## boundary of features column from first to last # 
  L <- which(colnames(osccCleanNA_conTCGA) == "Gender")
  R <- which(colnames(osccCleanNA_conTCGA) == "margin")
  chiT <- data.frame(matrix(data = NA, nrow = R, ncol = 2)) # create a empty data.frame, 2D matrix as R*2
  #rown = 8
  freq <- data.frame(matrix(data = NA, nrow = 1, ncol = 4)) #, dimnames = list(c(1:rown+1), c("Var2", "L", "H"))))
  colnames(freq) <- c("Var2", 0, 1, "Features") # column 4 is for mapi
  #x freq <- array(NA, dim = c(5, 3, R)) # row, col, and R as 3D array
  
  #!!! PMM1 score## at col 13, 14
  #osccCleanNA_conTCGA_pos <- which(colnames(osccCleanNA_conTCGA) == "H.score_T") # position of PMM1 IHC score
  osccCleanNA_conTCGAM_pos <- which(colnames(osccCleanNA_conTCGA) == paste("PMM1", "_median", sep=""))
  #  exp_geneName <- t(osccCleanNA[, osccCleanNA_pos ])
  
  # # it should be done at marginS.R
  # # then [binomial of gene_median] resume the correlation tables
  # osccCleanNA[osccCleanNAM_pos] <-(osccCleanNA[, osccCleanNA_pos ] >= cutoff1) +0 # binomial after osccCleanNA 
  # # ***** addNA for counting all NA (e.g. there is 0, no 1, in margin-free cohort) in "M" "stage_2" or "margin"
  # # the “pathological” case of two different kinds of NAs which are treated differently: exclude = if (useNA == "no") c(NA, NaN)
  # # "unusual NA comes from addNA() as factor
  # # [2018/03/14] finallized debug
  # osccCleanNA$margin <- addNA(osccCleanNA$margin, ifany=F) # ifany=="always" add NA as a factor
  # #   0    1 <NA> => 3 levels of factor in "margin"
  # # 245   11    0 
  # #  osccCleanNA$M <- addNA(osccCleanNA$M, ifany=F) # "always" add NA as a factor
  # #  osccCleanNA$stage_2 <- addNA(osccCleanNA$stage_2, ifany=F) # "always" add NA as a factor
  # 
  
  # chisq.test(matrix(c(22, 5, 38, 21), ncol = 2), correct=F)$p.value # a example
  for (ii in L:R){ # generate freq table
    # Chi-square, Contingent table correlation, binary variables
    
    # build a contingency table https://www.rdocumentation.org/packages/base/versions/3.4.3/topics/table
    # Powered by DataCamp, it might be running online
    # t, table from col of PMM1_median vs col ii (L "Gender" to R "margin"); deparse the argument (colnames) by deparse.level 2
    # # {DEBUG 
    # ii<-L-1
    #     ii<-ii+1
    # # DEBUG}
    t<- NULL
    t <- table(osccCleanNA_conTCGA[,osccCleanNA_conTCGAM_pos], osccCleanNA_conTCGA[,ii], useNA = "ifany") #dnn=c(colnames(osccCleanNA_conTCGA[osccCleanNA_conTCGAM_pos])))  
    # useNA = "ifany"; "always" is for "margin" (with n=0 count)
    chiT[ii,1] <- colnames(osccCleanNA_conTCGA[ii]) # name list of feature variables from L to R
    check_p <- chisq.test(t)$p.value # retrieved in chiT$X2
    if (is.na(check_p)==T) {check_p <- chisq.test(t[1:2,1:2])$p.value}
    chiT[ii,2] <- check_p
    # chisq is sum( (o-e)^2/e ); the Na should be removed, You have zero frequencies in 2 counts.
    #warnings(): In chisq.test(t) : Chi-squared approximation may be incorrect
    #=> fisher.test(a) # Fisher exact test is better in small counts. 
    #table(t); #print(c("ii=",paste(ii)))
    # t
    # 0   5   6 122 123 
    # 2   1   1   1   1 
    #print(ii) # debug
    
    obs <- as.data.frame(chisq.test(t)$observed)
    
    
    # >table(t)
    # 50 57 62 76 
    # 1  1  1  1 
    # [1] "-0.0765952029885 \n c(127, 129)"
    #   Var1 Var2 Freq
    # 1    0    0  127
    # 2    1    0   129
    # 3    0    1    0
    # 4    1    1    0
    # debug
    # with error # mdata <- melt(obs, id=c("Var2", "Var1")) # using the colnames of table t
    # id variables not found in data:" Var2, Var1 [2018/03/27] , only error on "Gene name = ZSWIM2"(20482)
    
    #print(paste("Run", run100, geneName, "(", which(whole_genome==geneName)," ): obs", obs, t))
    # under cutoff finding 100
    
    # ***ZSWIM2 skip (?) 
    # Error when "melt": id variables not found in data: Var2, Var1; (360: mdata <- melt(obs, id=c("Var2", "Var1")) # using the colnames of table t)
    #if (is.na(tryCatch([expression], error = function(e) return(NA)))) {return(list("skip", chiT, freq))}
    #if (obs$Freq==c(127,129,"NA","NA")) {print(paste("skip", geneName)); return(list("skip", chiT, freq))} # back to marginS.R
    # warnings():  the condition has length > 1 and only the first element will be used
    # [3] "0.238555836961 \n c(170, 75, 9, 2, 0, 0)"
    #   Var1 Var2 Freq
    # 1    0    0  170
    # 2    1    0   75
    # 3    0    1    9
    # 4    1    1    2
    # 5    0 <NA>    0
    # 6    1 <NA>    0
    # debug
    if (is.na(tryCatch(mdata <- melt(obs, id.vars=c("Var2", "Var1")), error = function(e) return(NA)))) {return(list("skip", chiT, freq))} # back to marginS.R
    # mdata <- melt(obs, id.vars=c("Var2", "Var1")) # error when there is only one group
    # debug for "ZSWIM2"(20482)
    
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
    # DEBUG
  }
  freq <- freq[-1,] #removal of 1st row: NA
  name_freq <- colnames(osccCleanNA_conTCGA[L:R])
  name_freq <- t(name_freq) # "Gender" "age.at.diagnosis" "T"  "N"  "M"  "stage_2" and "margin"
  
  # array indexing https://cran.r-project.org/doc/manuals/r-release/R-intro.html#Array-indexing
  results <- list (TRUE, chiT, freq) # it x should to be returned/updated the osccCleanNA_conTCGA
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




##
library(R.utils) # intToBin()
library(compositions) # unbinary()
contingencyBin <- function (osccCleanNA_conBin, chiT, freq) {
  # to generate tableChi1 (table 2)
  # enough, geneName is not necessary a parameter (osccCleanNA_conBin <- osccCleanNA)
  # processing chiT and processing freq; place a "remark" and Matthews
  # for binary and integer convertion
  # library(R.utils) # intToBin()
  # library(compositions) # unbinary()
  # sigContig <- c(NA, "*") # remarkable when (p<0.05 || scoreContig == c(5,10))
  # 
  freq_features <- as.character(freq$Features[seq(1, nrow(freq), 3)]) # every 3 rows, pick one
  #c("Gender","age.at.diagnosis", "T", "N", "M", "stage_2","margin")
  # rows need to be selected and reordered
  # ***match freq and osccCleanNA_conBin of colnames
  colnames(osccCleanNA_conBin)[2:8] <- freq_features
  # ***a feature-name mapping (indTbChi) of " freq$Features  <- osccCleanNA" ; it needs to be updated.
  indTbChi <- data.frame(cbind(featuresUni[-length(featuresUni)],
                               freq_features,
                               colnames(osccCleanNA_conBin)[2:8])) # from Gender to margin
  colnames(indTbChi) <- c("featuresUni", "freq$Features", "osccCleanNA_conBin")
  # > indTbChi as a data.frame (mapping table)
  # #               featuresUni-1   freq$Features    osccCleanNA_conBin(osccCleanNA)
  # 1                 Gender           Gender       Gender
  # 2       Age at diagnosis age.at.diagnosis        age.at.diagnosis
  # 3    Pathologic T status                T         T
  # 4    Pathologic N status                N         N
  # 5    Pathologic M status                M         M
  # 6       Pathologic Stage          stage_2        stage_2
  # 7 Surgical Margin status           margin       margin
  # #
  # rownames(tableChi1) <- featuresUni
  
  
  # reset tableChi1
  tableChi1 <- as.data.frame(setNames(replicate(length(contLowHighN), numeric(0), simplify = F), contLowHighN)) # declare a xlsx output data.frame (dynamic table)
  # colnames(tableChi1) <- contLowHighN
  for (i in 1:(length(featuresUni)-1)) { # without PMM1 on the last row (i in 1:7)
    # generate table 2/tableChi1 by every two rows
    mapi <- as.character(freq$Features) == as.character(indTbChi$osccCleanNA_conBin[i]) # "Gender"... in [47:49] of freq
    mapi_pos <- which(mapi == T) #[47:49], every 3
    var2U <- as.data.frame(freq[mapi_pos[1], ])
    var2L <- as.data.frame(freq[mapi_pos[2], ])
    subtotal <- sum(var2U[2:3]) + sum(var2L[2:3])
    U2 <- var2U[2];  U3 <- var2U[3]
    L2 <- var2L[2];  L3 <- var2L[3]
    pContig <- round(chiT[which(chiT$X1 == as.character(indTbChi$osccCleanNA_conBin[i])), 2], 4) # P-value retrieving from chiT$X2
    
    abin <- as.numeric(c(U3>U2, L3>L2, L2>U2, L3>U3)) # pattern of Sn=L3/(L2+L3), Sp=U2/(U2+U3)
    #    abin <- c(1,0,1,0) ; binary 1010 == decimal 10; while binary 0101 == decimal 5.
    achar <- paste(abin, collapse="") # ex. "0001"
    #  # a integer score, 10 or 5 and more is significant, which will be marked as "*" in remark
    scoreContig <- unbinary(achar) # 'structure(list(), *) warning() :-)
    # Matthews correlation coefficient: 0-1 perfect prediction
    #[OFF it]:    matthews <- (U2*L3-U3*L2)/sqrt((L3+U3)*(L3+L2)*(U2+U3)*(U2+L2))
    matthews <- 0
    
    # colUni_0 and colUni_1: value of feature => "male"/"female"....
    row0T2 <- data.frame(t(c(colUni_0[i], as.numeric(c(U2, round(U2/subtotal*100, 1), U3, round(U3/subtotal*100, 1), subtotal, pContig, (as.numeric((!(pContig>0.05) || scoreContig == c(5,10,1,2,4,7,8,11,13,14)))+0), matthews ))))) # remark 0 or 1
    row1T2 <- data.frame(t(c(colUni_1[i], as.numeric(c(L2, round(L2/subtotal*100, 1), L3, round(L3/subtotal*100, 1), NA, NA, NA, NA)))))
    #  colnames(row0) <- contLowHighN
    #  row0[,2:ncol(tableChi1)] <- as.numeric(as.character(row0[,2:ncol(tableChi1)]))
    tableChi1 <- rbind(tableChi1, row0T2, row1T2)
  }
  # tableChi1 <- cbind(tableChi1, chiT[-(1:5), 2]) # can not direct paste the p-value vector
  colnames(tableChi1) <- contLowHighN
  
  return(tableChi1)
} # end of contingencyBin define ++++++

# contingencyBin2  (spared) ####
# library(R.utils) # intToBin()
# library(compositions) # unbinary()
# contingencyBin2 <- function (osccCleanNA, chiT, freq) { 
#   # enough, geneName is not necessary a parameter
#   # generate Table 2 by processing chiT and processing freq; place a "remark" and Matthews
#   # for binary and integer convertion
#   # library(R.utils) # intToBin()
#   # library(compositions) # unbinary()
#   # sigContig <- c(NA, "*") # remarkable when (p<0.05 || scoreContig == c(5,10))
#   # for table 2
#   contLowHighN <- c("Features",	"Low", "(%)",	"High",	"(%)", "Case no",	"P-value", "Remark", "Matthews")
#   
#   # rows need to be selected and reordered
#   # *** mapping featuresUni <- osccCleanNA ; it needs to be corrected.
#   #  indTbChi <- data.frame(cbind(featuresUni[-length(featuresUni)], c("Gender","age.at.diagnosis", "T", "N", "M", "stage_2","margin" )) )
#   indTbChi <- data.frame(cbind(featuresUni[-length(featuresUni)], c("Gender","ageDx", "pathologic_T", "pathologic_N", "pathologic_M", "stage","margin" )) )
#   colnames(indTbChi) <- c("featuresUni", "osccCleanNA")
#   
#   #
#   # rownames(tableChi2) <- featuresUni
#   # reset tableChi2
#   tableChi2 <- as.data.frame(setNames(replicate(length(contLowHighN), numeric(0), simplify = F), contLowHighN)) # declare a xlsx output data.frame (dynamic table)
#   # colnames(tableChi2) <- contLowHighN
#   for (i in 1:(length(featuresUni)-1)) { # without geneName_expression in last row
#     # generate table by every two rows
#     mapi <- as.character(freq$Features) == as.character(indTbChi$osccCleanNA[i]) # "gender" in [47:49] of freq
#     mapi_pos <- which(mapi == T) #[47:49]
#     var2U <- as.data.frame(freq[mapi_pos[1], ])
#     var2L <- as.data.frame(freq[mapi_pos[2], ])
#     subtotal <- sum(var2U[2:3]) + sum(var2L[2:3])
#     U2 <- var2U[2];  U3 <- var2U[3]
#     L2 <- var2L[2];  L3 <- var2L[3]
#     pContig <- round(chiT[which(chiT$X1 == indTbChi$osccCleanNA[i]), 2], 4) # P-value retrieving from chiT$X2
#     abin <- as.numeric(c(U3>U2, L3>L2, L2>U2, L3>U3)) # pattern of Sn=L3/(L2+L3), Sp=U2/(U2+U3)
#     #    abin <- c(1,0,1,0) ; binary 1010 == decimal 10; while binary 0101 == decimal 5.
#     achar <- paste(abin, collapse="")
#     scoreContig <- unbinary(achar) # a integer score, 10 or 5 and more is significant, which will be marked as "*" in remark
#     # Matthews correlation coefficient: 0-1 perfect prediction
#     matthews <- (U2*L3-U3*L2)/sqrt((L3+U3)*(L3+L2)*(U2+U3)*(U2+L2))
#     row0T2 <- data.frame(t(c(colUni_0[i], as.numeric(c(U2, round(U2/subtotal*100, 1), U3, round(U3/subtotal*100, 1), subtotal, pContig, (as.numeric((!(pContig>0.05) || scoreContig == c(5,10,1,2,4,7,8,11,13,14)))+0), matthews ))))) # remark 0 or 1
#     row1T2 <- data.frame(t(c(colUni_1[i], as.numeric(c(L2, round(L2/subtotal*100, 1), L3, round(L3/subtotal*100, 1), NA, NA, NA, NA)))))
#     #  colnames(row0) <- contLowHighN
#     #  row0[,2:ncol(tableChi2)] <- as.numeric(as.character(row0[,2:ncol(tableChi2)]))
#     tableChi2 <- rbind(tableChi2, row0T2, row1T2)
#   }
#   # tableChi2 <- cbind(tableChi2, chiT[-(1:5), 2]) # can not direct paste the p-value vector
#   colnames(tableChi2) <- contLowHighN
#   
#   return(tableChi2)
# } # end of contingencyBin2 define +++++++++++++++++++++++++++++

# ++++++++= end of Function defination ++++++++++
## [reset for Repeat] end +++++++++++++++++++

#
#
#
#
#
# [run] mainA or mainB ####
# after functions defined




# # #> table(ZSWIM2$X2) # "skip XK" XKRY(19642)
# 
# 0   1   2   3 
# 2 211  10   1 
# 
# 543:  quantile(exp_geneName, c(0.3, 0.7)) 
# Error in quantile.default(exp_geneName, c(0.3, 0.7)) : 
# missing values and NaN's not allowed if 'na.rm' is FALSE
# 
# # D17904D <- load(file="ZSWIM2_1.Rda") # saved ZSWIM2 data.frame from save("D17904D", envir=.rs.CachedDataEnv, file="ZSWIM2_1")
# D17904D$X3 <- addNA(D17904D$X2, ifany=F) # always factorize including NA
# # $ X3: Factor w/ 4 levels "0","2","3",NA
# > table(D17904D$X3)
# Summary of 800 genes run:
#    0     2     3  <NA>
#   361    18     1 20119
#   return(0) is OK
#   error return(1)"ZSWIM2 skip", return(2)"survdiff", return(3)"no such file: LUAD.mRNA.Exp.Fire"
#   error return(4)"XKRY"quantile.default(exp_geneName, c(0.3, 0.7)); 
#   


#### [mainA process #part A] { ####
# genome-wide scan for margin 0 and 1 cohort ###
# added tryCatch for Survdiff.fit error

start_time <- Sys.time() # counted in minutes
# file.exists or list.files()
desti_ZSWIM2 <- "ZSWIM2_archive.Rda" # appending summary data of error codes or P-value
if (file.exists(desti_ZSWIM2)) {load(file=desti_ZSWIM2)} else
{ZSWIM2 <- data.frame(matrix(data = NA, nrow = LUAD_n, ncol = 2))}
aa <- LUAD_n; bb<- 1

## [2018/04/05] 00:38-08:48, 8 hours, run ZZZ3(20499) to ZNF268(20023)

# debug
# (OK) ZSWIM2_1 <-  20482 #skip "ZSWIM2"(20482), ZSCAN10(20466), ZPBP(20460), "ZFP91.CNTF"(19866) or ZFP91, 
#  "XK"(19643), "skip?":"XKRY"(19642), "TSNAX.DISC1"(18757), "TRPA1"(18700), TRIL...STX10(17406),
# ***  hit good P-value: ZZEF1(20498), ZYG11B(20496), 18700 TRPA1, TRPC1, TRO..., TRIP13(best)...STX10 (for 16 hours),
# STRAP, SLC8A1(16396) for 9 hours; SLC30A3(16243), SLC2A9(16239), 
# SLC2A4 ( 16234 ) "GLUT4",  [2018/04/10] start to fully run: SLC26A2 ( 16204 )", SLC26A7,
#  VANGL2(19297) (Rstudio down once) :-)
# (:))Error on "TSNAX.DISC1", "TRIL"(18602) (P-value=0.0281) "hit" but rt=26????, and "TRILB3" rt=7, #

#main_i loops ####
# #> start_time;end_time
# [1] "2018-04-11 00:19:10 CST"
# [1] "2018-04-17 02:56:45 CST"
# aa <- 16209; bb<- 1
#aa <- 19508-1; 
#bb<- 19508 #"WDR63"
#xx ddply(whole_genome[main_i], , survival_marginS, failwith(NA, fun(x){...}, quiet = TRUE))
aa <- 18758-1 # [2018/06/19] run04 done
#bb<- 13666


# # Best good: by mclapply,
# mclapply(whole_genome[aa:bb], survival_marginFree)
# #error on : [2018-05-27 08:59:51] [error] handle_read_frame error: websocketpp.transport:7 (End of File)# #save(ZSWIM2, file=desti_ZSWIM2)
# save(ZSWIM2, file="ZSWIM2_archive.Rda")
# apply(), MARGIN=c(1,2) # applys over rows and over columns
# 
# x#pforeach(main_i = aa:bb, .cores=2, .errorhandling = c("stop"))({
# Second good: by for loop, for save the ZSWIM2 data
for (main_i in aa:bb) {
  ZSWIM2[main_i, 2] <- survival_marginS(whole_genome[main_i]) # codes at source("TCGA_LUAD_marginS.R")
  # gene scan; return() at X2; for loop, we need ZSWIM2 data to be saved
  save(ZSWIM2, file=desti_ZSWIM2)
}
##}


## ## email for P-value analysisi of "ZSWIM2_archive.Rda (main_i)"####
end_time <- Sys.time()
body_text <- data.frame(matrix(data = NA, nrow = LUAD_n, ncol = 2)) # aa>bb == True ; aa-bb+1
colnames(body_text) <- c("P-value") #, "body_text")

#debug
#aa<- 20499; bb<- 16349
main_j <- which(!is.na(ZSWIM2$X2)) # error !=0
for (main_i in main_j[length(main_j)]:main_j[1]) {
  #main_i <- 0 # it needs to be delared before pforeach usage
  
  #print(paste("main_i", main_i))
  geneName <- whole_genome[main_i] # ok on "ZSWIM3"(20483), then error on "ZSWIM2"(20482)
  #debug call (survival_marginS)
  rt <- ZSWIM2[main_i,2] #<- survival_marginS(geneName) # return either 0(ok) or 1(***ZSWIM2 "error prone") or 2 # traceback() if error
  rt <- rt +1 # return()+1 = rt, as (P-value)+1 or range 1-6
  
  if (is.na(rt)) {next()}
  # multiple conditions: #rt as 1 or 2~6
  #with(ZSWIM2[main_i], 0<X2 & X2<=0.05), or 1< rt <= 1.05} # we hit one gene with P-value < 0.05 in KM plot
  # body_text$X1[main_i]<- (rt-1) # this kind of indexing is a vector/matrix operation
  if ((0+1)<rt & rt<=(0.05+1)) {body_text[main_i, 1]<- (rt-1) ; body_text[main_i, 2] <- paste(geneName, ", we hit one gene with P-value as ", rt-1," in KM plot", sep="")
  } else if (rt<=6) {
    body_text[main_i, 2] <- switch(rt, "OK", paste(geneName, "ZSWIM2-skip from contingencyTCGA()."), paste(geneName, "has only one group in survdiff."), paste("There is no", geneName, " RNAseq data."), paste(geneName, "quantile.default with NA."), paste(geneName, "with merge error."))
  } 
  # rt from 0~5 => 1~6     #body_text <- "OK" # return(0) rt==1
  # rt==0 is "OK"
  # if (rt==2) {body_text <- paste(geneName, "ZSWIM2 skip from contingencyTCGA()"); print(body_text)}
  # # if return 1 => one group issue in Melt() in contingencyTCGA()
  # if (rt==3) {body_text <- paste(geneName, "has only one group in survdiff."); print(body_text)}
  # # if return 2 => one group issue in survdiff.
  # if (rt==4) {body_text <- paste(geneName, "there is NO such RNAseq."); print(body_text)}
  # # return 3: there is NO such RNAseq file: "LUAD.mRNA.Exp.ZFP91.CNTF.Fire.Rda"
  # if (rt==5) {body_text <- paste(geneName, "quantile.default with NA"); print(body_text)}
  # # quantile.default(exp_geneName, c(0.3, 0.7)) by "XKRY"19642
  # if (rt==6) {body_text <- paste(geneName, "merge error"); print(body_text)}
  # # merge(oscc, LUAD_T_mRNA.Fire...) applied to a non-vector by "TSNAX.DISC1"(18757)
  if (rt>6){
    body_text[main_i, 2] <- paste(geneName, ": rt=", rt, ", I don't know what's happened XD.")
  }
  #  print(paste("(", main_i, ")", body_text[main_i, 2]))
  #ZSWIM2$X3 <- addNA(ZSWIM2$X2, ifany=F) # factorized for summary
  # save(ZSWIM2, file="ZSWIM2_archive.Rda")
  # dt <- strptime(Sys.time(), format="%Y-%m-%d %H:%M:%S")
  # if ((rt==1) & ((as.numeric(dt$min) %% 10)!=0)) {next()} # (rt=1, return(0)) go to next iteration, no email
  # x# sending email when error OR every 10 mins
}

### sending cluster email
# the progress, with error code {
# authentication of gmail account in browser, once for all? yes
# pipe %>%
#formatC(body_text[which((0< body_text$`P-value`) & (body_text$`P-value` <=0.05)), ])
body_pvalue <- data.frame(matrix(data = NA, nrow = LUAD_n, ncol = 2))
# body_pvalue[1,]<-matrix(c(0.04,"SLC2A4"))
body_pvalue <- body_text[which((0< body_text$`P-value`) & (body_text$`P-value` <=0.05)), ]# %>%
#body_pvalue <- body_text 
#save(body_pvalue, file="gmailr_pvalue.Rda")
write.csv(body_text, file="gmailr_pvalue.csv") #debug body_pvalue .csv is empty?
#  paste0 # as.character then concatenating
#body_pvalue <- paste(body_pvalue, sep=" ")
#formatC # as character
# %>% gettextf
if (length(body_pvalue) > 0) {
  texchi2 <- mime() %>%
    to(recipient) %>%
    from(sender) %>%
    subject(paste("Rstudio!!! Significant finding:", "There is", length(body_pvalue), "genes in LUAD")) %>%
    html_body(paste("Moreover, there is P-value < 0.05, see attached file.")) %>%
    attach_file(file="gmailr_pvalue.csv")
  # email notifying only significant P-value genes (X1)
} else {
  texchi2 <- mime() %>%
    to(recipient) %>%
    from(sender) %>%
    subject(paste("Rstudio:", whole_genome[main_j[length(main_j)]], "...", whole_genome[main_j[1]], ",", length(main_j), "genes in LUAD"))
  #html_body("There is P-value < 0.05, see attached file.")
  # email notifying scanne genes (X1) with non-zero return
}

#create_draft(texchi2)
print("Sending you a email report...")
if (is.na(tryCatch( send_message(texchi2), error = function(e) return(NA)))) {
} else {print("Done !")}

## in case internet is lost=>
# (Ok)Error in curl, Could not resolve host: www.googleapis.com
#} Auto-refreshing stale OAuth token.Id: 162a34f5a88e06f1


#} # the end of whole genome scanning
print(paste("Expended duration:", end_time - start_time, "hours"))
# run to here




# [Results] ####
#====== Analysis of output .Rda of _marginS_ or _marginFree_
# [2018/06/20] they are stored at ./run04_marginS_
path_ZSWIM2 <- file.path(path_LUAD, "run04_marginS_")

# ZSWIM2_archive1000_20180408_0042_0933.Rda; 9 hours for 1,000 genes to be scanned
# 3 hours for 1,000 genes to be scanned under GCP Rstudio server
# STRAP to SLC8A1
#which(whole_genome==geneName)
#
#_marginS_
# merge ZSWIM2a ZSWIM2b and ZSWIM2c => ZSWIM2abc_archive.Rda [2018/06/19]
# para_seg <- c(20499, 13528, 6833, 1) # ZZZ3-> PLCE1, PNMT-> GARNL3, GAR1-> A1BG
load(file=file.path(path_ZSWIM2, "ZSWIM2c_archive.Rda")) #
ZSWIM2 <- ZSWIM2abc #(done, merged)
save(data=ZSWIM2, file=file.path(path_ZSWIM2, "ZSWIM2_archive.Rda")) # or ZSWIM2abc_archive 
# OR
# _marginFree_
load(file="ZSWIM2_free_archive.Rda") # as well as ZSWIM2



#*** Dissection of ZSWIM2 #
#1) KM plot is not at best cutoff, why??
#2) ZSWIM2 should be saved by appending.
#3)
load(file=file.path(path_ZSWIM2, "ZSWIM2_archive.Rda")) # merged abc
ZSWIM2$X3 <- addNA(ZSWIM2$X2, ifany=F)
ZSWIM2_2 <- table(ZSWIM2$X3) # try to tryCatch (1, 2)
# 0     1     2     3     4     5     6     7     8     9 
# 682    64   187    26    23    21     9    18    13    14 
# 10    11    12    13    14    15    16    17    18    19 
# 14    10    12    10     5     6     5     6     2     5 
# 20    21    22    23    24    25    26    27    28    29 
# 3     5     1     4     2     3     1     2     5     2 
# 30    31    34    35    36    37    40    41    44    49 
# 3     3     3     2     2     2     1     1     1     1 
# 51    52    57    60    66  <NA> 
#   1     2     1     1     1 19314 
ZSWIM2_2 <- data.frame(ZSWIM2_2)
#ZSWIM2_2$Freq[46] <- ZSWIM2_2$Freq[46] - (17406+(LUAD_n-18602)) #number of NA:11, from TRIL(18602) to "STX10"(17406)
#print(run16h <- (sum(ZSWIM2_2$Freq)-(LUAD_n-(18602-17406))))
# running 18602-17406=1196 genes in 16 hours
plot((ZSWIM2_2))
#> sum(ZSWIM2_2$Freq)
#[1] 20499  # scan completely

# reture(0): ok for analysis
text_pie <- paste("Workable genes \n in ", ZSWIM2_2$Freq[1]/sum(ZSWIM2_2$Freq)*100, " %") # _marginS_ about 50.96% usable RNAseq data
pie(ZSWIM2_2$Freq, labels=ZSWIM2_2$Var1, main=text_pie)
# OR _marginFree_ with 52%(?) usable RNAseq data
# deal with return(1), return(2) and return(3) errors, try to solve it
# _marginS_
n_percent_Bonferroni <- ZSWIM2_2$Freq[1]/sum(ZSWIM2_2$Freq) # 0.5095858
error01_sample <- which(ZSWIM2$X3==1) #  4.6% error01: "ZSWIM2 skip from contingencyTCGA()"): one group issue in Melt()
error02_sample <- which(ZSWIM2$X3==2) # 17.4% error02: There has only one group in survdiff.
error03_sample <- which(ZSWIM2$X3==3) # 0.42% error03: There has only one group in survdiff.
# OR _marginFree_
n_percent_Bonferroni <- ZSWIM2_2$Freq[1]/sum(ZSWIM2_2$Freq) # 0.52
error01_sample <- which(ZSWIM2$X3==1) #  5.0% error01: "ZSWIM2 skip from contingencyTCGA()"): one group issue in Melt()
error02_sample <- which(ZSWIM2$X3==2) # 18.4% error02: There has only one group in survdiff.
error03_sample <- which(ZSWIM2$X3==3) # 3.4% error03: There has only one group in survdiff.
# done #






# (Both) Retrieving the summary table to form z-score in LUAD (marginS), from all LUAD_survival*.Rda ####
# _marginFree_ or _marginS_ from .Rda
# # [choice ONE]: _marginFree_ or _marginS_ loading from .Rda
SFree <- "_marginS_"
#OR _marginFree_
SFree <- "_marginFree_"
# get_Rda_pvalue <- function(geneName) {
#   load(file=paste("LUAD_survivalAnalysis_marginS_", geneName, ".Rda", sep=""))
#   # load list = c("tableChi1", "tableOS1", "tableRFS1", "OS_pvalue", "RFS_pvalue")
#   # a example: geneName <- "TRIP13"
#   #OS_pvalue$p_OS[which.min(OS_pvalue$p_OS)]
#   
#   if (nrow(OS_pvalue) > 0)  { # we hit this gene with P-value < 0.05 in KM plot
#     return(min(OS_pvalue$p_OS))} else {return(NA)}
# }

# candidate_sample(KM) and candidate_cox (Cox list), two parts
# candidate_sample: Kaplan–Meier plots for genes significantly associated with survival.
# Table1: list genes, which it's KM plot with P-value < 0.05, and ranking by P-value (ascending) and z-score

aa <- LUAD_n; bb<- 1
#aa <- ZSWIM2_2$Freq[1]  # n=9445
# however, $ ls LUAD_sur*.Rda | wc -l
# 16194

#{ declare empty data.frame and list
#*
candidate_sample <- data.frame(matrix(data = NA, nrow = aa, ncol = 3)) # for num, gene ID and it's P-value
colnames(candidate_sample) <- c("number", "gene_id", "p_value") # OS P-value only (in LUAD); 
#"number" position in ZSWIM2
#candidate_sample$number <- data.frame(which(ZSWIM2$X3==0)) # gene "number"
candidate_sample$gene_id <- whole_genome[bb:aa] # retrieving gene name from whole_genome; return() at X2
#

#* a list for all cox survival datas
# http://www.cookbook-r.com/Manipulating_data/Converting_between_data_frames_and_contingency_tables/
# column name might be created by variable: e.x. df.sex[,"per"] <- df.sex$count1/sum(df.sex$count1); # df.sex$per
candidate_cox <- replicate(aa, list()) # it doesn't need colnames
#data.table(matrix(data = NA, nrow = aa, ncol = 3)) # for num, gene ID and it's P-value
# colnames(candidate_cox) <- c( "Features", "HR",       "P_value_uni",  "sig",      "Features", "HR",      
#                               "P_value_multi",  "sig",      "Features", "P_value_KM",  "sig" )
# #}


setwd(path_ZSWIM2) # set for a while (for ip loop)
for (ip in (bb:aa)) {
  geneName <- candidate_sample$gene_id[ip]
  #print(paste("At", path_ZSWIM2, "=> (", ip, ")", geneName), sep="")
  #candidate_sample$p_value <- lapply(unlist(candidate_sample$gene_id[bb:aa]), get_Rda_pvalue) # retrieving P-value by gene name from .Rda; return() at X3
  #[choice ONE]: _marginFree_ or _marginS_ loading from .Rda
  #  _marginS_
  
  load_filename <- file.path(paste("LUAD_survivalAnalysis", SFree, geneName, ".Rda", sep=""))
  #OR (automatically defined)
  #_marginFree_
  #load_filename <- paste("LUAD_survivalAnalysis", SFree, geneName, ".Rda", sep="")
  #
  #
  if (load_filename %in% dir()) {
    
    load(file=load_filename)
    #      if (is.na(tryCatch(load(file=load_filename), error = function(e) return(NA)))) {} # load file with error free :-)
    #  load(file=paste("LUAD_survivalAnalysis_marginS_", geneName, ".Rda", sep=""))
    # load list = c("tableChi1", "tableOS1", "tableRFS1", "OS_pvalue", "RFS_pvalue")
    # a example: geneName <- "TRIP13"
    print(paste(SFree, "with how many OS P-values: ", nrow(OS_pvalue)))
    if (nrow(OS_pvalue) > 0) {
      candidate_sample$p_value[ip] <- min(OS_pvalue$p_OS)
      candidate_sample$number[ip] <- nrow(OS_pvalue) #number of OS P-values in this gene
    }
    #candidate_sample$p_value[ip] <- get_Rda_pvalue(candidate_sample$gene_id[ip]) # retrieving P-value by gene name from .Rda; return() at X3
    #  print(paste(ip, geneName, ": ", candidate_sample$p_value[ip], sep=" "))
    
    #***
    # Significant features from table 2, table 3 and table 4 of output .Rda files ##
    #{
    # [uni_cox_pvalue and uni_HR]
    # [multi_cox_pvalue and multi_HR]
    # [exp_pvalue]
    # from tableChi1 (Table 2), tableOS1 (Table 3) /[tableRFS1]
    
    # focusing on tableOS1, Cox proportiopnal hazard model: (both univariate and multivariate) P-value <= 0.05,
    # and sorted by HR > 1 vs HR < 1
    # P-value notes:
    # Significant codes:  0 ???***??? 0.001 ???**??? 0.01 ???*??? 0.05 ???.??? 0.1 ??? ??? 1
    # if <0.001 => mark as "***"
    # osHR$X4[osHR$X4<0.001] <- "***"
    
    # uni_cox_pvalue
    uni_cox_pvalue <- tableOS1[c(FALSE, TRUE), c(1,2,5)] # get significant P-value, hazard ratio at column 2
    uni_cox_pvalue$`P-value` <- as.character(uni_cox_pvalue$`P-value`)
    uni_cox_pvalue$`P-value`[(uni_cox_pvalue$`P-value` == "***")] <- "0.001"
    uni_cox_pvalue$`P-value` <- as.numeric(uni_cox_pvalue$`P-value`)
    uni_cox_pvalue$sig <- uni_cox_pvalue$`P-value` <= 0.05 # "sig" marking for significant
    
    # multi_cox_pvalue
    multi_cox_pvalue <- tableOS1[c(FALSE, TRUE), c(1,6,9)] # get significant P-value, hazard ratio at column 6
    multi_cox_pvalue$`P-value` <- as.character(multi_cox_pvalue$`P-value`)
    multi_cox_pvalue$`P-value`[(multi_cox_pvalue$`P-value` == "***")] <- "0.001"
    multi_cox_pvalue$`P-value` <- as.numeric(multi_cox_pvalue$`P-value`)
    multi_cox_pvalue$sig <- multi_cox_pvalue$`P-value` <= 0.05 # "sig" marking for significant
    
    # *** we don't have RFS in LUAD cohort :-)
    
    # and tableChi1: P-value <= 0.05 # exp_pvalue
    # "sig" marking the significant "features": check "odd" position
    exp_pvalue <- tableChi1[c(TRUE, FALSE), c(1,7)] # get significant P-value from column 7, name from column 1; remark at column 8
    exp_pvalue$`P-value` <- as.numeric(as.character(exp_pvalue$`P-value`))
    exp_pvalue$sig <- exp_pvalue$`P-value` <= 0.05 # "sig" marking for significant
    
    # # merging them as one by common row names
    # > colnames(exp_pvalue)
    exp1 <- data.frame("X1"= geneName, "X2" = NA,  "X3"=NA) # append one row, with this geneName
    colnames(exp1) <- colnames(exp_pvalue) # precisely matched
    candidate_cox_ip <- cbind(uni_cox_pvalue, multi_cox_pvalue, rbind(exp_pvalue, exp1)) # colname is duplicated, bind 3 tables together
    colnames(candidate_cox_ip) <- c( "uni_Features", "uni_HR",       "uni_P_value",  "uni_sig",      "multi_Features", "multi_HR",
                                     "multi_P_value",  "multi_sig",      "KM_Features", "KM_P_value",  "KM_sig") # KM_sig is Remark at table 2(tableChi1)
    
    # it must be [[]] for assign the content !!! https://stat.ethz.ch/R-manual/R-devel/library/base/html/Extract.html.
    candidate_cox[[ip]] <- candidate_cox_ip #http://cran.r-project.org/doc/manuals/R-lang.html#Indexing
    print(paste(ip, "cox features saved; ncol = ", length(candidate_cox[[ip]]))) # a vector in [] is ok for indexing
    # length == 11 => that is correct !
    #}
  }
} # end of ip for loop

#_marginS_ or _marginFree_ by SFree; saving on ./run04_marginS_, files => 17030
save(candidate_sample, candidate_cox, n_percent_Bonferroni, file=file.path(path_ZSWIM2, paste("LUAD_OS", SFree, "pvalueKM_candidate_cox.Rda", sep=""))) #ok; with KM_sig and Remark, and Cox HR
setwd(path_LUAD) 
#_marginS_ by SFree
# saved file="LUAD_OS_marginS_pvalueKM_candidate_cox.Rda" above
#_marginFree_ by SFree
#save(candidate_sample, candidate_cox, n_percent_Bonferroni, file="LUAD_OS_marginFree_pvalueKM_candidate_cox.Rda") #ok; with KM_sig and Remark
#[2018/05/27]

# devtools::install_github(c('jeroenooms/jsonlite', 'rstudio/shiny', 'ramnathv/htmlwidgets', 'timelyportfolio/listviewer'))
# library(listviewer)
# jsonedit( candidate_cox )
# # finished (both) processes




## Post1 process _marginS_ ####
##  table1 (KM): candidate_sample(KM) ##
# #sort by mpg (ascending) and cyl (descending)
# 楊老師: Bonferroni correction (adjustment) sets the significance cut-off at α/n. of KM P-value; Cox P-value <=0.05 as well.
# https://www.stat.berkeley.edu/~mgoldman/Section0402.pdf
# newdata <- mtcars[order(mpg, -cyl),]
load(file=file.path(path_ZSWIM2, "LUAD_OS_marginS_pvalueKM_candidate_cox.Rda")) # as candidate_sample with candidate_cox and n_percent_Bonferroni 
attach(candidate_sample)
LUAD_OS_marginS_pvalue_sorted <- candidate_sample[order(p_value, -number),] # sorting by order(ascending)
detach(candidate_sample)

attach(LUAD_OS_marginS_pvalue_sorted)
# Bonferroni correction (adjustment) sets the significance cut-off at α/n 
# {
alpha_LUAD <- 0.05
Bonferroni_cutoff <- alpha_LUAD / (LUAD_n * n_percent_Bonferroni) # as 5.302486e-06
# => 4.786521e-06 (run04)
# if no Bonferroni: 
#Bonferroni_cutoff <- alpha_LUAD / 1
# } Bonferroni end

LUAD_OS_marginS_pvalue005_sorted <- LUAD_OS_marginS_pvalue_sorted[which(p_value<=alpha_LUAD & !is.na(p_value)), 1:3]
detach(LUAD_OS_marginS_pvalue_sorted)  # keeping drawing
# no correction: n=6260 in _marginS_;  and n=6345 in _marginFree_ of LUAD_OS_marginS_pvalue005_sorted


#  plot(OS.km, lty=1, col=c("blue","red"), sub=paste("p-value =", p_OS[j]), main=paste("OS in OSCC(n=", surv_OS$n[1]+surv_OS$n[2],")/", geneName, "cutoff=", i), ylab="Percent Survival", xlab="Days")
#  legend("topright", legend=c(paste("low(",surv_OS$n[1], ")"), paste("high(",surv_OS$n[2], ")")), lty=1:2, col=c("blue","red"))
library(stats)
library(scales)
library(minpack.lm)

attach(LUAD_OS_marginS_pvalue005_sorted)
# reg <- lm(number ~ p_value, data = LUAD_OS_marginS_pvalue005_sorted)
# abline(reg, col="blue")
# or
LUAD_OS_marginS_pvalue005_sorted$z_score <- scale(number, scale = T) # "number" frequency to z-scores ("Z" because the normal distribution is also known as the "Z distribution").
LUAD_OS_marginS_pvalue005_sorted$number_01 <- scales:::rescale(z_score, to = c(0, 1)) # rescaled to range minnew to maxnew (aka. 0 to 1 for binomial glm)

plot(p_value, number_01, type="p", ylab="Z-score", xlab="P-value from KM survival analysis \n (cutoff by Bonferroni correction)", log="x") # log scale x or y
#axis(side=3, at=c(1e-7, 1e-6, 1e-5, 0.01, 0.05)) #1=below, 2=left, 3=above and 4=right
abline(h=0.6, lty=2, col="red")
abline(v=Bonferroni_cutoff, lty=2, col="red") # *** as 4.693073e-06
# run a logistic regression model (categorical 0 vs 1)
g <- glm(number_01 ~ p_value, family=poisson, data=LUAD_OS_marginS_pvalue005_sorted)
# (in this case, generalized linear model with log link)(link = "log"), poisson distribution
curve(predict(g, data.frame(p_value = x), type="response"), add=TRUE, col="blue") # draws a curve based on prediction from logistic regression model

# x# run a non-linear regression: non-linear least squares approach (function nls in R) ; A nice feature of non-linear regression in an applied context is that the estimated parameters have a clear interpretation (Vmax in a Michaelis-Menten model is the maximum rate) which would be harder to get using linear models on transformed data for example.
# # the Levenberg-Marquardt algorithm for nonlinear regression
# # "eyeball" plot to set approximate starting values
# # https://stats.stackexchange.com/questions/160552/why-is-nls-giving-me-singular-gradient-matrix-at-initial-parameter-estimates
# # nls "singular gradient matrix at initial parameter estimates" error, using nlsLM instead
# a_start <- 0.9 #param a is the y value when x=0
# b_start <- 2*log(2)/a_start #b is the decay rate
# m <- nlsLM(number_01 ~ p_value, data=LUAD_OS_marginS_pvalue005_sorted, start=list(a=a_start, b=b_start))
# curve(predict(m, data.frame(p_value = x), type="response", col="green"), add=TRUE, col="blue") # draws a curve based on prediction from logistic regression model
# #lines(p_value, predict(m), lty=2, col="green", lwd=3)

legend("topright", legend=c(paste("LR"), paste("Cutoff")), lty=1:2, col=c("blue","red"), cex=0.7) # box and font size
# figure legend: logistic regression, LR, by Generalized linear model, glm
#detach(LUAD_OS_marginS_pvalue005_sorted)

# then...
# a candidate table, with "number of P-value" under cutoff finding (becoming z-score: bigger, much more significant cutoff "sites")
# => a kind of local minimal or global minimal of curve fitting (Levenberg-Marquardt Optimization)?
# subsetting the candidated genes table, according P-value (Bonferroni_cutoff) and Z-score
# ***after Bonferroni correction => n=34 in _marginS_
# ***after Bonferroni correction => n=33 in _marginFree_
#attach(LUAD_OS_marginS_pvalue005_sorted)
LUAD_OS_marginS_pvalue1e_6_zscore0_6 <- LUAD_OS_marginS_pvalue005_sorted[which(p_value<=Bonferroni_cutoff & z_score>=0.6), 2:5]
LUAD_OS_marginS_pvalue1e_6_zscore0_6[, 2] <- signif(LUAD_OS_marginS_pvalue1e_6_zscore0_6[, 2], 3)
detach(LUAD_OS_marginS_pvalue005_sorted) # n=17
# > colnames(LUAD_OS_marginS_pvalue1e_6_zscore0_6)
# [1] "gene_id"   "p_value"   "z_score"   "number_01"
# 
# _marginS_ (ok)
save(LUAD_OS_marginS_pvalue1e_6_zscore0_6, Bonferroni_cutoff, file=file.path(path_ZSWIM2, "LUAD_OS_marginS_pvalue1e_6_zscore0_6.Rda")) # cutoff by Bonferroni_cutoff
save(LUAD_OS_marginS_pvalue005_sorted, file=file.path(path_ZSWIM2, "LUAD_OS_marginS_pvalue005_sorted.Rda"))
# x # # *** OR _marginFree_ (x not here)
# LUAD_OS_marginFree_pvalue1e_6_zscore0_6 <- LUAD_OS_marginS_pvalue1e_6_zscore0_6
# save(LUAD_OS_marginFree_pvalue1e_6_zscore0_6, file="LUAD_OS_marginFree_pvalue1e_6_zscore0_6.Rda")
# LUAD_OS_marginFree_pvalue005_sorted <- LUAD_OS_marginS_pvalue005_sorted
# save(LUAD_OS_marginFree_pvalue005_sorted, file="LUAD_OS_marginFree_pvalue005_sorted.Rda")





## _marginS_ ( )
## Post2 process and add table2 (cox): candidate_cox ####
# "sig" marking for significant P-value (<=0.05)
# [uni_cox_pvalue, uni_HR, uni_sig]
# [multi_cox_pvalue, multi_HR, multi_sig]
# [exp_pvalue]
# to generate LUAD_OS_marginS_pvalue005KM_sorted_pvalueCox_HR.Rda; n=??

#LUAD_OS_marginS_pvalue005_sorted # n=6584
load(file=file.path(path_ZSWIM2, "LUAD_OS_marginS_pvalue005_sorted.Rda")) # as LUAD_OS_marginS_pvalue005_sorted
load(file=file.path(path_ZSWIM2, "LUAD_OS_marginS_pvalueKM_candidate_cox.Rda")) # as candidate_cox (a list, Bonferroni_cutoff), since 2018/05/16
# candidate_sample, candidate_cox, n_percent_Bonferroni

# new variable: KM + Cox, n=6261
LUAD_OS_marginS_pvalue005KM_sorted_pvalueCox_HR <- LUAD_OS_marginS_pvalue005_sorted
#dataframe[,"newName"] <- NA # add more named columns (NA)
LUAD_OS_marginS_pvalue005KM_sorted_pvalueCox_HR[,c("uni_HR", "uni_P_value", "uni_sig","multi_HR", "multi_P_value", "multi_sig")] <- NA
#LUAD_OS_marginS_pvalue005_sorted[,c(colnames(candidate_cox[[1]][8, c(2:4, 6:8)]))] <- NA
#colnames(LUAD_OS_marginS_pvalue005KM_sorted_pvalueCox_HR) <- c(...."uni_HR", "uni_P_value", "uni_sig",
#                                                "multi_HR", "multi_P_value", "multi_sig")

#attach(LUAD_OS_marginS_pvalue005KM_sorted_pvalueCox_HR)
#ipp <- 0
for (ip in 1:nrow(LUAD_OS_marginS_pvalue005KM_sorted_pvalueCox_HR)) {
  pos_gene <- which(whole_genome==LUAD_OS_marginS_pvalue005KM_sorted_pvalueCox_HR$gene_id[ip]) # LUAD_OS_marginS_pvalue005KM_sorted_pvalueCox_HR$gene_id
  # by geneid ?
  if (length(candidate_cox[[pos_gene]])==11) {
    # reference: candidate_cox_ip, which listing as candidate_cox
    LUAD_OS_marginS_pvalue005KM_sorted_pvalueCox_HR[ip, 6:11]  <- candidate_cox[[pos_gene]][8, c(2:4, 6:8)] # taking RNAseg: sig marked and P-values, HRs
    #   ipp <- ipp+1
    print(paste(ip, " out of ", nrow(LUAD_OS_marginS_pvalue005KM_sorted_pvalueCox_HR), " (", whole_genome[pos_gene], ")", sep=""))
    print(paste("..added...", LUAD_OS_marginS_pvalue005KM_sorted_pvalueCox_HR[ip, c(2,6)], sep=";"))
  }
  #?? error on : [2018-05-18 07:26:42] [error] handle_read_frame error: websocketpp.transport:7 (End of File)
  #"exp_pvalue" is a correlation of gene expression vs features (TNM....): 
  # <- candidate_cox[which(gene_id[ip]==whole_genome)][1:7, c(11)] # correlation
  # colnames <-  c("KM_Features", "KM_P_value",  "KM_sig")
  
}
#detach(LUAD_OS_marginS_pvalue005KM_sorted_pvalueCox_HR)
# c( "uni_Features", "uni_HR",       "uni_P_value",  "uni_sig",      "multi_Features", "multi_HR",
#"multi_P_value",  "multi_sig",      "KM_Features", "KM_P_value",  "KM_sig")
# add colname as HRs, "uni_cox"/"multi_cox", sig ... for HRs, P-values, sig of RNAseq(z-score)

##
# save table1 + table2 (ok) [2018/05/18]
save(LUAD_OS_marginS_pvalue005KM_sorted_pvalueCox_HR, file=file.path(path_ZSWIM2, "LUAD_OS_marginS_pvalue005KM_sorted_pvalueCox_HR.Rda"))

#save(list(LUAD_OS_marginS_pvalue1e_6_zscore0_6, ???))
# x[i], or might be x[i:j]
# x[i, j]
# x[[i]]; x[[expr]]; it can NOT be x[[i:j]]
# x[[i, j]]
# x$a
# x$"a"




# ** Pickup all significant genes list -> LUAD_OS_marginS_THREE_pvalue005 ####
# _marginS_
# [2018/05/18][2018/05/28]
# > colnames(LUAD_OS_marginS_pvalue005KM_sorted_pvalueCox_HR)
# [1] "number"        "gene_id"       "p_value"      
# [4] "z_score"       "number_01"     "uni_HR"       
# [7] "uni_P_value"   "uni_sig"       "multi_HR"     
# [10] "multi_P_value" "multi_sig"
LUAD_OS_marginS_THREE_pvalue005 <- subset(LUAD_OS_marginS_pvalue005KM_sorted_pvalueCox_HR, (uni_P_value <= 0.05) & (multi_P_value <= 0.05), 
                                          select=c(number, gene_id, p_value, uni_HR, uni_P_value, multi_HR, multi_P_value))
# n=4337, KM P-value only (NOT by Bonferroni_cutoff)
save(LUAD_OS_marginS_THREE_pvalue005, Bonferroni_cutoff, file=file.path(path_ZSWIM2, "LUAD_OS_marginS_THREE_pvalue005.Rda")) 
# as LUAD_OS_marginS_THREE_pvalue005, Bonferroni_cutoff <- 5.302486e-06

#*** RR>reproducible research resume: ####
load(file=file.path(path_ZSWIM2,"LUAD_OS_marginS_THREE_pvalue005.Rda"))
# <<<

# plot uni_HR, n=4131
attach(LUAD_OS_marginS_THREE_pvalue005)
plot(p_value, uni_HR, type="p", ylab="Cox Uni_HR", xlab="P-value from KM survival", log="x")
#axis(side=3, at=c(1e-7, 1e-6, 1e-5, 0.01, 0.05)) #1=below, 2=left, 3=above and 4=right
abline(h=1, lty=2, col="red")
abline(v=Bonferroni_cutoff, lty=2, col="red") #Bonferroni_cutoff
# run a logistic regression model (categorical 0 vs 1)
#g <- glm(uni_HR ~ p_value, family=poisson, data=LUAD_OS_marginS_THREE_pvalue005)
# (in this case, generalized linear model with log link)(link = "log"), poisson distribution
#curve(predict(g, data.frame(p_value = x), type="response"), add=TRUE, col="blue") # draws a curve based on prediction from logistic regression model
#legend("topright", legend=c(paste("LR"), paste("Cutoff")), lty=1:2, col=c("blue","red"), cex=0.7) # box and font size
# figure legend: logistic regression, LR, by Generalized linear model, glm
detach(LUAD_OS_marginS_THREE_pvalue005)


#plot(LUAD_OS_marginS_THREE_pvalue005$p_value,  LUAD_OS_marginS_THREE_pvalue005$uni_HR)
# plot multi_HR, n=4132
attach(LUAD_OS_marginS_THREE_pvalue005)
plot(p_value, multi_HR, type="p", ylab="Cox Multi_HR", xlab="P-value from KM survival", log="x")
#axis(side=3, at=c(1e-7, 1e-6, 1e-5, 0.01, 0.05)) #1=below, 2=left, 3=above and 4=right
abline(h=1, lty=2, col="red")
abline(v=Bonferroni_cutoff, lty=2, col="red")
# run a logistic regression model (categorical 0 vs 1)
#g <- glm(uni_HR ~ p_value, family=poisson, data=LUAD_OS_marginS_THREE_pvalue005)
# (in this case, generalized linear model with log link)(link = "log"), poisson distribution
#curve(predict(g, data.frame(p_value = x), type="response"), add=TRUE, col="blue") # draws a curve based on prediction from logistic regression model
#legend("topright", legend=c(paste("LR"), paste("Cutoff")), lty=1:2, col=c("blue","red"), cex=0.7) # box and font size
# figure legend: logistic regression, LR, by Generalized linear model, glm
detach(LUAD_OS_marginS_THREE_pvalue005)


## Venn_marginS_: cross matching HR>1 of Uni+Multi, HR<1 of Uni+Multi ####
#library(venn) #https://cran.r-project.org/web/packages/venn/venn.pdf
# venn(x, snames = "", ilabels = FALSE, counts = FALSE, ellipse = FALSE,
#      zcolor = "bw", opacity = 0.3, size = 15, cexil = 0.6, cexsn = 0.85,
#      borders = TRUE, ...)
#library(gplots)
# To get the list of gene present in each Venn compartment we can use the gplots package
#library(gplots) # capture the list of genes from venn

#{ pickup1 from LUAD_OS_marginS_THREE_pvalue005; Bonferroni_cutoff
#* Cox HR (>1 or) >=2.5 (bad guy genes) # & (uni_P_value <= 0.05) & (multi_P_value <= 0.05)
LUAD_OS_marginS_uni_CoxHR2p5 <- subset(LUAD_OS_marginS_THREE_pvalue005, (p_value <= Bonferroni_cutoff) & (uni_HR >=2.5), 
                                       select=c(number, gene_id, p_value, uni_HR, uni_P_value, multi_HR, multi_P_value))
# n=10
# {
#   $`uni_Cox HR >=2.5`
#   [1] "TFF1"    "UCK2"    "RRM2"    "EIF5AL1" "DKK1"   
#   [6] "FAM111B" "ECT2"    "KIF20A"  "KIF23"   "DSCC1" 
# # > LUAD_OS_marginS_uni_CoxHR2p5$gene_id #manuscript: with Bonferroni_cutoff, n=10; however, genes list is not the same :-)
# [1] "PLCD3"   "CYTSB"   "TFF1"    "TFAP2A"  "UCK2"    "RRM2"   
# [7] "KIF14"   "EIF5AL1" "DKK1"    "FAM111B"
# }


#.. # & (uni_P_value <= 0.05) & (multi_P_value <= 0.05) 
LUAD_OS_marginS_multi_CoxHR2p5 <- subset(LUAD_OS_marginS_THREE_pvalue005,  (p_value <= Bonferroni_cutoff) & (multi_HR >=2.5), 
                                         select=c(number, gene_id, p_value, uni_HR, uni_P_value, multi_HR, multi_P_value))
# n=5
# {
#   $`multi_Cox HR >=2.5`
#   [1] "DUSP5"  "SETD8"  "PLSCR1" "KLHDC5"
# > LUAD_OS_marginS_multi_CoxHR2p5$gene_id #manuscript: with Bonferroni_cutoff, n=5 and genes list is not the same (only DUSP5 remained) :-)
# [1] "PLCD3"  "CYTSB"  "TFAP2A" "KIF14"  "DUSP5"
# }
#...

# venn1 diagram of HR>=2.5 of Uni & Multi ###
#for list of genes by grouping; library(gplots)
venn_HR2p5 <- list(LUAD_OS_marginS_uni_CoxHR2p5$gene_id, LUAD_OS_marginS_multi_CoxHR2p5$gene_id)
# cutoff by Bonferroni_cutoff
names_HR2p5 <- c("uni_Cox HR >=2.5", "multi_Cox HR >=2.5")
library(gplots)
tmp <- venn(venn_HR2p5, names=names_HR2p5, show.plot=F) #library(gplots); the group count matrix alone
isect_HR2p5 <- attr(tmp, "intersections")
#isect_HR2p5
detach(package:gplots)

library(venn)
venn(venn_HR2p5, snames=names_HR2p5,
     ilabels = T, counts = T, ellipse = FALSE, zcolor = "red, deeppink", opacity = 0.6, size = 15, cexil = 0.6, cexsn = 0.85, borders = TRUE)
#      predefined colors if "style"
#http://www.stat.columbia.edu/~tzheng/files/Rcolor.pdf
# meta-language 1 0 or -
title <- paste(c("LUAD survival analysis", "KM P-Value <= ", signif(Bonferroni_cutoff, 3)), sep = "", collapse = "\n")
#coords <- unlist(getCentroid(getZones(venn_HR2p5, snames="uni_CoxHR>=2p5, multi_CoxHR>=2p5")))
# coords[1], coords[2], 
text(500,900, labels = title, cex = 0.85) # (0,0) on bottom_left corner
# n=6
#{
# # $`uni_Cox HR >=2.5:multi_Cox HR >=2.5`
# [1] "PLCD3"    "CYTSB"    "TFAP2A"   "KIF14"    "C9orf140"
# [6] "RG9MTD2" 
# $`uni_Cox HR >=2.5:multi_Cox HR >=2.5`; 
# n=4 #manuscript: with Bonferroni_cutoff
# [1] "PLCD3"  "CYTSB"  "TFAP2A" "KIF14" 
# 
#}

#https://stackoverflow.com/questions/43324180/adding-legend-to-venn-diagram
#legend("top", legend=c("B:multi_Cox HR >=2.5", "A:uni_Cox HR >=2.5"), lty=1:2, col=c("blue","red"), cex=0.7) # box and font size
# figure legend: 



#===
#{ pickup2 from LUAD_OS_marginS_THREE_pvalue005
#* Cox HR <0.4 (good guy genes)

#...
LUAD_OS_marginS_uni_CoxHR0p5 <- subset(LUAD_OS_marginS_THREE_pvalue005, (p_value <= Bonferroni_cutoff) & (uni_P_value <= 0.05) & (multi_P_value <= 0.05) & (uni_HR <0.4), 
                                       select=c(number, gene_id, p_value, uni_HR, uni_P_value, multi_HR, multi_P_value))
# n=8
# {
#  # $`uni_Cox HR <0.4`
# [1] "DBP"     "TXNDC11" "ZNF709"  "PDK2"    "PLEKHB1"
# # > LUAD_OS_marginS_uni_CoxHR0p5$gene_id #manuscript Bonferroni_cutoff => n=10
# [1] "DBP"       "CRHR2"     "MYLIP"     "ZNF682"    "TXNDC11"  
# [6] "SLC11A2"   "ZNF709"    "NUP210L"   "PDK2"      "LOC284440"
# }

# ...
LUAD_OS_marginS_multi_CoxHR0p5 <- subset(LUAD_OS_marginS_THREE_pvalue005, (p_value <= Bonferroni_cutoff) & (uni_P_value <= 0.05) & (multi_P_value <= 0.05) & (multi_HR <0.4), 
                                         select=c(number, gene_id, p_value, uni_HR, uni_P_value, multi_HR, multi_P_value))
# n=7
# {
#   # $`multi_Cox HR <0.4`
#   # [1] "PXMP4"     "FBP1"      "PLD3"      "ACAD8"     "ZNF14"    
#   # [6] "LIFR"      "LOC151174" "NFATC2"    "MEX3A"    
# #  LUAD_OS_marginS_multi_CoxHR0p5$gene_id #manuscript Bonferroni_cutoff => n=7
# [1] "CRHR2"     "MYLIP"     "ZNF682"    "SLC11A2"   "NUP210L"  
# [6] "PXMP4"     "LOC284440"
# }


# venn2 diagram of HR < 0.4 of Uni & Multi ###
#for list of genes by grouping; library(gplots)
venn_HR0p5 <- list(LUAD_OS_marginS_uni_CoxHR0p5$gene_id, LUAD_OS_marginS_multi_CoxHR0p5$gene_id)
names_HR0p5 <- c("uni_Cox HR < 0.4", "multi_Cox HR < 0.4")
library(gplots)
tmp <- venn(venn_HR0p5, names=names_HR0p5, show.plot=F) #library(gplots); the group count matrix alone
isect_HR0p5 <- attr(tmp, "intersections")
isect_HR0p5
detach(package:gplots)

library(venn)
venn(venn_HR0p5, snames=names_HR0p5,
     ilabels = T, counts = T, ellipse = FALSE,  zcolor = "forestgreen, darkolivegreen1", opacity = 0.6, size = 15, cexil = 0.6, cexsn = 0.85, borders = TRUE)
#      predefined colors if "style"
# meta-language 1 0 or -, https://cran.r-project.org/web/packages/gplots/vignettes/venn.pdf
title <- paste(c("LUAD survival analysis", "KM P-Value <= ", signif(Bonferroni_cutoff, 3)), sep = "", collapse = "\n")
text(500,900, labels = title, cex = 0.85) # (0,0) on bottom_left corner

#{ intersection of them (HR <0.4), n=8
# $`uni_Cox HR <0.4:multi_Cox HR <0.4`
# [1] "CRHR2"        "MYLIP"        "ZNF682"       "SLC11A2"
# [5] "NUP210L"      "LOC284440"    "LOC100130093" "ZKSCAN4"
# 
#   # $`uni_Cox HR < 0.4:multi_Cox HR < 0.4` #manuscript Bonferroni_cutoff => n=6
# [1] "CRHR2"     "MYLIP"     "ZNF682"    "SLC11A2"   "NUP210L"  
# [6] "LOC284440"
# }

#x #library(VennDiagram)
# VENN.LIST <- list(LUAD_OS_marginS_uni_CoxHR0p5$gene_id, LUAD_OS_marginS_multi_CoxHR0p5$gene_id)
# #venn.plot <- venn.diagram(VENN.LIST , NULL, fill=c("darkmagenta", "darkblue"), alpha=c(0.5,0.5), cex = 2, cat.fontface=4, category.names=c("A", "B"), main="LUAD: Cox HR <0.5")
# # To plot the venn diagram we will use the grid.draw() function to plot the venn diagram
# #grid.draw(venn.plot)
# # We can summarize the contents of each venn compartment, as follows:
# # in 1) ConditionA only, 2) ConditionB only, 3) ConditionA & ConditionB
# lapply(isect_HR0p5, head)
#}


#} venn the end




# Excluding the LUAD cancer driver genes list??? ##
#BioXpress*.csv # data from https://hive.biochemistry.gwu.edu/cgi-bin/prd/bioxpress/servlet.cgi




## Export r2excel and .Rda ####
# _marginS_ [2018/06/27]
# sink() for .xlsx export as well :-) https://stackoverflow.com/questions/34038041/how-to-merge-multiple-data-frame-into-one-table-and-export-to-excel
# [2018/05/29] => refinement of LUAD_OS_marginS_candidates_Venn.xlsx: show up Bonferroni_cutoff and 5e-6 (expression style).
# 
load(file=file.path(path_ZSWIM2, "LUAD_OS_marginS_pvalue1e_6_zscore0_6.Rda")) # in LUAD_OS_marginS_pvalue1e_6_zscore0_6, cutoff by Bonferroni_cutoff
load(file=file.path(path_ZSWIM2, "LUAD_OS_marginS_pvalueKM_candidate_cox.Rda")) # in candidate_cox (a list), candidate_sample, candidate_cox, n_percent_Bonferroni


library("xlsx")
library("r2excel")
# Create an Excel workbook. Both .xls and .xlsx file formats can be used.
filenamex <- paste("LUAD_OS", SFree, "candidates_Venn", ".xlsx", sep = "") # "LUAD_OS_marginS_candidates_Venn.xlsx"
wb <- createWorkbook(type="xlsx")

# Create a sheet in that workbook
sheet <- xlsx::createSheet(wb, sheetName = paste("Survival_candidates"))
# [add data row by row, start from column 2]
#+++++++++++++++++++++++++++++++
## Add paragraph : Author
library("tis") # by Brian Salzer
# today(), arg must be ti, tis, ts, tif, or tifName
author <- paste("Reported by Tex Li-Hsing Chi. \n",
                "tex@gate.sinica.edu.tw \n", SFree, "\n", Sys.Date(), sep="")
xlsx.addParagraph(wb, sheet, value=author, isItalic=TRUE, colSpan=5, 
                  rowSpan=4, fontColor="darkgray", fontSize=24)
xlsx.addLineBreak(sheet, 3)
# header
xlsx.addHeader(wb, sheet, value=paste("Table 1. The candiate genes expressed in ", TCGA_cohort,
                                      " (ranking by KM P-value, selected by Bonferroni cutoff, ", signif(Bonferroni_cutoff, 3), ") ", "\n", "n= ", nrow(LUAD_OS_marginS_pvalue1e_6_zscore0_6), sep=""),
               level=5, color="black", underline=0)
#xlsx.addHeader(wb, sheet, value=paste("Cutoff at ", round(cutoff1, 3), " (", percent(surv_OS1$n[1]/(surv_OS1$n[1]+surv_OS1$n[2])), ")", sep = ""),
#               level=5, color="red", underline=0) # total n is taken from surv_OS1$n

xlsx.addLineBreak(sheet, 1) # add one blank line

#xlsx.addTable(wb, sheet, data = t(data.frame(c(paste(geneName, "expression"), "", paste(geneName, "expression"), "", "(Optimised)"))), fontSize=12, startCol=4,
#              fontColor="darkblue", row.names = F, col.names = F) #, colSpan=1, rowSpan=1)

# a candidate genes table
# > colnames(LUAD_OS_marginS_pvalue1e_6_zscore0_6)
# [1] "gene_id"   "p_value"   "z_score"   "number_01"
colnames(LUAD_OS_marginS_pvalue1e_6_zscore0_6) <- c("Gene_id",   "P_value",   "z_score_raw",   "Z_score") # Z_score is rescaled as 0-1
# Bonferroni_cutoff; [, c(1,2,4)]
# in scientific notation: formatC of [, c(2)]
attach(LUAD_OS_marginS_pvalue1e_6_zscore0_6)
candidates_Bonferroni_pvalue <- LUAD_OS_marginS_pvalue1e_6_zscore0_6[which(P_value<=Bonferroni_cutoff), c(1,2,4)]
detach(LUAD_OS_marginS_pvalue1e_6_zscore0_6)

attach(candidates_Bonferroni_pvalue) # removal of as.factor (P_value)
candidates_Bonferroni_pvalue$P_value <- formatC(as.numeric(as.character(P_value)), format = "e", digits = 2)
candidates_Bonferroni_pvalue$Z_score <- signif(Z_score, digits=4)
candidates_Bonferroni_pvalue <- candidates_Bonferroni_pvalue[order(P_value, -Z_score), ] #sorting by order(ascending)
detach(candidates_Bonferroni_pvalue)

xlsx.addTable(wb, sheet, data = candidates_Bonferroni_pvalue, startCol=2,
              fontColor="darkblue", fontSize=12,
              rowFill=c("white", "white"), row.names = TRUE)

xlsx.addLineBreak(sheet, 5)  # add two blank lines

# Export z-score/P-Value plot by plotFunction() ####
xlsx.addLineBreak(sheet, 1)

# Define function
xlsx.addPlot.candidates<-function( wb, sheet, startRow=NULL,startCol=2,
                                   width=480, height=480,... )
{ # plot of z-score vs p-value digram from "the summary"
  
  png(filename = "plot.png", width = width, height = height,...)
  #{
  # plot(OS.km, lty=1, xscale=12, xmax=60, col=c("blue","red"), sub=paste("P-value =", p_OS1), main=paste("OS in TCGA ", TCGA_cohort, "(n=", surv_OS1$n[1]+surv_OS1$n[2],")/", geneName), ylab="Percent Survival", xlab="Years")
  # legend("topright", legend=c(paste("low(",surv_OS1$n[1], ")"), paste("high(",surv_OS1$n[2], ")")), lty=1:2, col=c("blue","red"))
  load(file=file.path(path_ZSWIM2, "LUAD_OS_marginS_pvalue005_sorted.Rda")) # as LUAD_OS_marginFree_pvalue005_sorted
  # "swap data" for following code running
  
  attach(LUAD_OS_marginS_pvalue005_sorted)
  plot(p_value, number_01, type="p", ylab="Z-score", xlab="P-value from KM survival analysis \n (cutoff by Bonferroni correction)", log="x") # log scale x or y
  #axis(side=3, at=c(1e-7, 1e-6, 1e-5, 0.01, 0.05)) #1=below, 2=left, 3=above and 4=right
  abline(h=0.6, lty=2, col="green")
  abline(v=Bonferroni_cutoff, lty=2, col="red") # *** as 4.693073e-06
  # run a logistic regression model (categorical 0 vs 1)
  g <- glm(number_01 ~ p_value, family=poisson, data=LUAD_OS_marginS_pvalue005_sorted)
  # (in this case, generalized linear model with log link)(link = "log"), poisson distribution
  curve(predict(g, data.frame(p_value = x), type="response"), add=TRUE, col="blue") # draws a curve based on prediction from logistic regression model
  legend("topright", legend=c(paste("LR"), paste("Cutoff Z"), paste("Cutoff B")), lty=1:2, col=c("blue","green","red"), cex=0.9) # box and font size
  # figure legend: logistic regression, LR, by Generalized linear model, glm
  detach(LUAD_OS_marginS_pvalue005_sorted)
  #}
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
} # Define function

# calling
currentRow <- xlsx.addPlot.candidates(wb, sheet) # startRow + a plot
print(paste("The z-score summary plot: ", filenamex, " was exported successfully."))



#..
# Export Venn1 by plotFunction() ####
xlsx.addLineBreak(sheet, 1)

# Define function
xlsx.addPlot.venn<-function( wb, sheet, startRow=NULL,startCol=2,
                             width=480, height=480,... )
{ # plot of venn digram"
  
  png(filename = "plot1.png", width = width, height = height,...)
  #{
  venn(venn_HR2p5, snames=names_HR2p5,
       ilabels = T, counts = T, ellipse = FALSE,  zcolor = "red, deeppink", opacity = 0.6, size = 15, cexil = 3, cexsn = 0.85, borders = TRUE)
  #      predefined colors if "style"; cexil = 0.6 (default text size of counts)
  title <- paste(c("LUAD survival analysis", "KM P-Value <= ", signif(Bonferroni_cutoff, 3)), sep = "", collapse = "\n")
  text(500,900, labels = title, cex = 0.85) # (0,0) on bottom_left corner
  
  #}
  dev.off()
  
  
  #Append plot1 to the sheet
  if(is.null(startRow)){
    rows<- getRows(sheet) #list of row object
    startRow=length(rows)+1
  } 
  # Add the file created previously
  addPicture("plot1.png", sheet=sheet,  startRow = startRow, startColumn = startCol) 
  xlsx.addLineBreak(sheet, round(width/20)+1)
  res<-file.remove("plot1.png")
  
} # end of Define function

# calling
#detach(package:gplots)
library(venn)
xlsx.addPlot.venn(wb, sheet)
print(paste("The Venn diagram (bad guy) and the candidate genes on", filenamex, ", which was exported successfully."))




# *generate (bad) prognostic features of those genes on lists in TCGA LUAD cohort ####
# store KM_sig remrark as a Byte; converted as base10 in list_KM_sigBin
load(file=file.path(path_ZSWIM2, "LUAD_OS_marginS_pvalueKM_candidate_cox.Rda")) # as candidate_cox (a list), since 2018/05/16; with KM_sig, Remark
# candidate_sample, candidate_cox, n_percent_Bonferroni
load(file=file.path(path_ZSWIM2, "LUAD_OS_marginS_THREE_pvalue005.Rda") )# as as LUAD_OS_marginS_THREE_pvalue005, Bonferroni_cutoff

#.. bad guy gene candidate
# => see "remark=1"; "which" or NOT "which", the order is wrong.
# isect_HR2p5 is a list
#geneid_bad_uni_HR2p5 <- c(isect_HR2p5$`A:B`, isect_HR2p5$`A`) # from venn_HR2p5: A + A:B = group uni_
geneid_bad_uni_HR2p5 <- c(isect_HR2p5[[3]], isect_HR2p5[[1]])
# whole_genome[which(whole_genome %in% geneid_bad_uni_HR2p5)]
# [1] "DKK1"    "DSCC1"   "ECT2"    "EIF5AL1" "FAM111B" "KIF20A"  "KIF23"  
# [8] "RRM2"    "TFF1"    "UCK2" 
candidate_bad_uni_HR2p5 <- cbind(data.frame(gene_id=geneid_bad_uni_HR2p5), LUAD_OS_marginS_THREE_pvalue005[which(LUAD_OS_marginS_THREE_pvalue005$gene_id %in% geneid_bad_uni_HR2p5), 3:7])
#as.binary(candidate_cox[which(whole_genome %in% geneid_bad_uni_HR2p5)]$KM_sig, n=7, logic=T) # store KM_sig remrark as a Byte
# can NOT use "whole_genome" any more
#x list_KM_sig <- candidate_cox[which(whole_genome %in% geneid_bad_uni_HR2p5)] 
list_KM_sig <- candidate_cox[as.numeric(rownames(candidate_bad_uni_HR2p5))] # a list (rownames has whole_genome's right position)
# km <- function(x) {print(list_KM_sig[[x]][8,9])}; lapply(c(1:10), km)
list_KM_sigBin <- data.frame(matrix(NA, nrow = length(geneid_bad_uni_HR2p5), ncol = nrow(candidate_cox[[1]])-1))
colnames(list_KM_sigBin) <- c(rownames(candidate_cox[[1]])[-8])
for (ikm in 1:length(list_KM_sig)) 
{
  #list_KM_sigBin[ikm, 2:nrow(candidate_cox[[1]])] <- 
  list_KM_sigBin[ikm,] <- as.numeric(list_KM_sig[[ikm]]$KM_sig[-8]) # or as .factor
}
candidate_bad_uni_HR2p5 <- cbind(candidate_bad_uni_HR2p5, list_KM_sigBin)
# ex. "RRM2" has good impact on TN and pathology when expression is high

# => see "remark=1"                                                # number, gene_id, p_value, uni_HR, uni_P_value, multi_HR, multi_P_value])
# # isect_HR2p5 is a list
geneid_bad_multi_HR2p5 <- c(isect_HR2p5[[3]], isect_HR2p5[[2]]) # from venn_HR2p5
candidate_bad_multi_HR2p5 <- cbind(data.frame(gene_id=geneid_bad_multi_HR2p5), LUAD_OS_marginS_THREE_pvalue005[(LUAD_OS_marginS_THREE_pvalue005$gene_id %in% geneid_bad_multi_HR2p5), 3:7])
# as.binary
list_KM_sig <- candidate_cox[as.numeric(rownames(candidate_bad_multi_HR2p5))] # a list
list_KM_sigBin <- data.frame(matrix(NA, nrow = length(geneid_bad_multi_HR2p5), ncol = nrow(candidate_cox[[1]])-1))
colnames(list_KM_sigBin) <- c(rownames(candidate_cox[[1]])[-8])
for (ikm in 1:length(list_KM_sig)) 
{
  #list_KM_sigBin[ikm, 2:nrow(candidate_cox[[1]])] <- 
  list_KM_sigBin[ikm,] <- as.numeric(list_KM_sig[[ikm]]$KM_sig[-8]) # or as .factor
}
candidate_bad_multi_HR2p5 <- cbind(candidate_bad_multi_HR2p5, list_KM_sigBin)
# ex. "DUSP5" has good impact on M stage

# => do NOT need to see "remark=1"
geneid_bad_unimulti_HR2p5 <- isect_HR2p5[[3]] #$`A:B`
candidate_bad_unimulti_HR2p5 <- cbind(data.frame(gene_id=geneid_bad_unimulti_HR2p5), LUAD_OS_marginS_THREE_pvalue005[(LUAD_OS_marginS_THREE_pvalue005$gene_id %in% geneid_bad_unimulti_HR2p5), 3:7])


# export tables
xlsx.addLineBreak(sheet, 4) # 
# header of candidate_bad_unimulti_HR2p5
xlsx.addHeader(wb, sheet, value=paste("Table 2A. The ", nrow(candidate_bad_unimulti_HR2p5), " consensus candiate genes (uni_ & multi_CoxHR>2.5) in ", TCGA_cohort, "(bad guy)"),
               level=5, color="black", underline=0)
xlsx.addLineBreak(sheet, 1) # add one blank line
candidate_bad_unimulti_HR2p5[,2] <- formatC(candidate_bad_unimulti_HR2p5[,2], format = "e", digits = 2) #as.numeric(as.character
xlsx.addTable(wb, sheet, data = candidate_bad_unimulti_HR2p5[, c(1:6)], startCol=2,
              fontColor="darkblue", fontSize=12,
              rowFill=c("white", "white"), row.names = F)
xlsx.addLineBreak(sheet, 2)  # add two blank lines

xlsx.addLineBreak(sheet, 3)
# header of candidate_bad_uni_HR2p5
xlsx.addHeader(wb, sheet, value=paste("Table 2B. The candiate genes (uni_CoxHR>2.5) in ", TCGA_cohort, "(bad guy)"),
               level=5, color="black", underline=0)
xlsx.addLineBreak(sheet, 1) # add one blank line
candidate_bad_uni_HR2p5[,2] <- formatC(candidate_bad_uni_HR2p5[,2], format = "e", digits = 2) #as.numeric(as.character
xlsx.addTable(wb, sheet, data = candidate_bad_uni_HR2p5[, c(1:6, 9:13)], startCol=2,
              fontColor="darkblue", fontSize=12,
              rowFill=c("white", "white"), row.names = F)
xlsx.addLineBreak(sheet, 2)  # add two blank lines

xlsx.addLineBreak(sheet, 3)
# header of candidate_bad_multi_HR2p5
xlsx.addHeader(wb, sheet, value=paste("Table 2C. The candiate genes (multi_CoxHR>2.5) in ", TCGA_cohort, "(bad guy)"),
               level=5, color="black", underline=0)
xlsx.addLineBreak(sheet, 1) # add one blank line
candidate_bad_multi_HR2p5[,2] <- formatC(candidate_bad_multi_HR2p5[,2], format = "e", digits = 2) #as.numeric(as.character
xlsx.addTable(wb, sheet, data = candidate_bad_multi_HR2p5[, c(1:6, 9:13)], startCol=2,
              fontColor="darkblue", fontSize=12,
              rowFill=c("white", "white"), row.names = F)
xlsx.addLineBreak(sheet, 2)  # add two blank lines


# Export Venn2 by plotFunction() ####
xlsx.addLineBreak(sheet, 1)

# Define function
xlsx.addPlot.venn<-function( wb, sheet, startRow=NULL,startCol=2,
                             width=480, height=480,... )
{ # plot of venn digram"
  
  
  png(filename = "plot2.png", width = width, height = height,...)
  #{
  venn(venn_HR0p5, snames=names_HR0p5,
       ilabels = T, counts = T, ellipse = FALSE, zcolor = "forestgreen, darkolivegreen1", opacity = 0.6, size = 15, cexil = 3, cexsn = 0.85, borders = TRUE)
  #      predefined colors if "style"
  title <- paste(c("LUAD survival analysis", "KM P-Value <= ", signif(Bonferroni_cutoff, 3)), sep = "", collapse = "\n")
  text(500,900, labels = title, cex = 0.85) # (0,0) on bottom_left corner
  
  #}
  dev.off()
  
  # append plot2
  if(is.null(startRow)){
    rows<- getRows(sheet) #list of row object
    startRow=length(rows)+1
  } 
  addPicture("plot2.png", sheet=sheet,  startRow = startRow , startColumn = startCol) 
  # jump to next venn ?
  xlsx.addLineBreak(sheet, round(width/20)+1)
  res<-file.remove("plot2.png")
} # end of Define function

# calling
#detach(package:gplots)
library(venn)
xlsx.addPlot.venn(wb, sheet)
print(paste("The Venn diagram2 and the candidate genes: ", filenamex, " were exported successfully."))


# *generate (good) prognostic features of those genes on lists in TCGA LUAD cohort ####
# => see "remark=1" 
geneid_good_uni_HR0p5 <- c(isect_HR0p5[[3]], isect_HR0p5[[1]]) #$`A:B`
candidate_good_uni_HR0p5 <- cbind(data.frame(gene_id=geneid_good_uni_HR0p5), LUAD_OS_marginS_THREE_pvalue005[(LUAD_OS_marginS_THREE_pvalue005$gene_id %in% geneid_good_uni_HR0p5), 3:7])
# number, gene_id, p_value, uni_HR, uni_P_value, multi_HR, multi_P_value])
# the list (order) of gene_id has error ???
list_KM_sig <- candidate_cox[as.numeric(rownames(candidate_good_uni_HR0p5))] # a list
list_KM_sigBin <- data.frame(matrix(NA, nrow = length(geneid_good_uni_HR0p5), ncol = nrow(candidate_cox[[1]])-1))
colnames(list_KM_sigBin) <- c(rownames(candidate_cox[[1]])[-8])
for (ikm in 1:length(list_KM_sig)) 
{
  #list_KM_sigBin[ikm, 2:nrow(candidate_cox[[1]])] <- 
  list_KM_sigBin[ikm,] <- as.numeric(list_KM_sig[[ikm]]$KM_sig[-8]) # or as .factor
}
candidate_good_uni_HR0p5 <- cbind(candidate_good_uni_HR0p5, list_KM_sigBin)
# ex."TXNDC11" higher expression has better prognosis impact as well as smaller T and less N and lower stage

# => see "remark=1" 
geneid_good_multi_HR0p5 <- c(isect_HR0p5[[3]], isect_HR0p5[[2]])
candidate_good_multi_HR0p5 <- cbind(data.frame(gene_id=geneid_good_multi_HR0p5), LUAD_OS_marginS_THREE_pvalue005[(LUAD_OS_marginS_THREE_pvalue005$gene_id %in% geneid_good_multi_HR0p5), 3:7])
# as.binary
list_KM_sig <- candidate_cox[as.numeric(rownames(candidate_good_multi_HR0p5))] # a list
list_KM_sigBin <- data.frame(matrix(NA, nrow = length(geneid_good_multi_HR0p5), ncol = nrow(candidate_cox[[1]])-1))
colnames(list_KM_sigBin) <- c(rownames(candidate_cox[[1]])[-8])
for (ikm in 1:length(list_KM_sig)) 
{
  #list_KM_sigBin[ikm, 2:nrow(candidate_cox[[1]])] <- 
  list_KM_sigBin[ikm,] <- as.numeric(list_KM_sig[[ikm]]$KM_sig[-8]) # or as .factor
}
candidate_good_multi_HR0p5 <- cbind(candidate_good_multi_HR0p5, list_KM_sigBin)
# ex."PLD3" has less N


# => do NOT need to see "remark=1"
geneid_good_unimulti_HR0p5 <- isect_HR0p5[[3]] #$`A:B`
candidate_good_unimulti_HR0p5 <- cbind(data.frame(gene_id=geneid_good_unimulti_HR0p5), LUAD_OS_marginS_THREE_pvalue005[(LUAD_OS_marginS_THREE_pvalue005$gene_id %in% geneid_good_unimulti_HR0p5), 3:7])
# *** using table 2 to interpretate the impact on survival by DEG of this candidate gene list 
# c("CRHR2", "MYLIP", "NUP210L", "ZKSCAN4")
# => see "remark=1" on LUAD_survivalAnalysis_marginS_CRHR2.xlsx => higher CRHR2 expression is associated with less LN metastaese
# => see "remark=1" on LUAD_survivalAnalysis_marginS_CRHR2.xlsx => higher MYLIP expression is associated with smaller tumor size (T)
# => see "remark=1" on LUAD_survivalAnalysis_marginS_CRHR2.xlsx => higher NUP210L expression is associated with less LN metastaese
# => see "remark=1" on LUAD_survivalAnalysis_marginS_CRHR2.xlsx => higher ZKSCAN4 expression is associated with smaller tumor size (T)

# export tables
xlsx.addLineBreak(sheet, 4) # 
# header of candidate_good_unimulti_HR0p5
xlsx.addHeader(wb, sheet, value=paste("Table 3A. The ", nrow(candidate_good_unimulti_HR0p5), " consensus candiate genes (uni_ & multi_CoxHR<0.5) in ", TCGA_cohort, "(good guy)"),
               level=5, color="black", underline=0)
xlsx.addLineBreak(sheet, 1) # add one blank line
candidate_good_unimulti_HR0p5[,2] <- formatC(candidate_good_unimulti_HR0p5[,2], format = "e", digits = 2) #as.numeric(as.character
xlsx.addTable(wb, sheet, data = candidate_good_unimulti_HR0p5[, c(1:6)], startCol=2,
              fontColor="darkblue", fontSize=12,
              rowFill=c("white", "white"), row.names = F)
xlsx.addLineBreak(sheet, 2)  # add two blank lines

xlsx.addLineBreak(sheet, 3)
# header of candidate_good_uni_HR0p5
xlsx.addHeader(wb, sheet, value=paste("Table 3B. The candiate genes (uni_CoxHR<0.5) in ", TCGA_cohort, "(good guy)"),
               level=5, color="black", underline=0)
xlsx.addLineBreak(sheet, 1) # add one blank line
candidate_good_uni_HR0p5[,2] <- formatC(candidate_good_uni_HR0p5[,2], format = "e", digits = 2) #as.numeric(as.character
xlsx.addTable(wb, sheet, data = candidate_good_uni_HR0p5[, c(1:6, 9:13)], startCol=2,
              fontColor="darkblue", fontSize=12,
              rowFill=c("white", "white"), row.names = F)
xlsx.addLineBreak(sheet, 2)  # add two blank lines

xlsx.addLineBreak(sheet, 3)
# header of candidate_good_multi_HR0p5
xlsx.addHeader(wb, sheet, value=paste("Table 3C. The candiate genes (multi_CoxHR<0.5) in ", TCGA_cohort, "(good guy)"),
               level=5, color="black", underline=0)
xlsx.addLineBreak(sheet, 1) # add one blank line
candidate_good_multi_HR0p5[,2] <- formatC(candidate_good_multi_HR0p5[,2], format = "e", digits = 2) #as.numeric(as.character
xlsx.addTable(wb, sheet, data = candidate_good_multi_HR0p5[, c(1:6, 9:13)], startCol=2,
              fontColor="darkblue", fontSize=12,
              rowFill=c("white", "white"), row.names = F)
xlsx.addLineBreak(sheet, 2)  # add two blank lines


##..
# save the workbook to an Excel file and write the file to disk.####
setwd(path_ZSWIM2)
xlsx::saveWorkbook(wb, filenamex)
setwd(path_LUAD)
#xlsx.openFile(filenamex) # open file to review


#*** ok and tar until here (refinement ok) ####
# the END of R2Excel ###
#}

# https://david.ncifcrf.gov/conversion.jsp?VFROM=NA DAVID pathway analysis
david_bad <- merge(candidate_bad_uni_HR2p5, candidate_bad_multi_HR2p5, by="gene_id", all=TRUE) #joint by union
# FAM111B (Family With Sequence Similarity 111 Member B) => converting to Entrez Gene ID(374393)
# > david_bad$gene_id
# [1] CYTSB   DKK1    EIF5AL1 FAM111B KIF14   PLCD3   RRM2    TFAP2A 
# [9] TFF1    UCK2    DUSP5
# 
david_good <- merge(candidate_good_uni_HR0p5, candidate_good_multi_HR0p5, by="gene_id", all=TRUE) #joint by union
# > david_good$gene_id
# [1] CRHR2     DBP       LOC284440 MYLIP     NUP210L  
# [6] PDK2      TXNDC11   ZNF709    PXMP4     SLC11A2  
# [11] ZNF682    

# tar
#{
#好用的工具
# ***bash $ TODAY=`date +"%b %d"`;ls -l | grep "$TODAY"
# # list today's files
# https://www.howtogeek.com/248780/how-to-compress-and-extract-files-using-the-tar-command-on-linux/
#   https://www.gnu.org/software/tar/manual/tar.html
# $ info tar # tar -t --list  —remove-file...
# =tar and scp from a list of files .xlsx (candidate genes list); $ 
#   $ tar -czvf ~/marginS_xlsx.tar.gz -T ~/marginS_list.txt   # or  list as many directories
# $ scp  tex@35.201.169.0:~/margin*_xlsx.tar.gz ./
#   } tex@instance-4:$ ~/R/LUAD_Peter_survival$ sudo mv LUAD_survivalAnalysis_marginS*.* ./survivalAnalysis_marginS/
#   tar -xzvf archive.tar.gz -C /tmp # to extract them to /tmp
# tar -xzf archive.tar.gz --overwrite

# $ cd run04_marginS_/
cur_wd <- getwd()
setwd(paste(cur_wd, "run04_marginS_/", sep="/"))
david_bad_list <- paste("LUAD_survivalAnalysis_marginS_", david_bad$gene_id, ".xlsx", sep="")
write.table(david_bad_list, file="david_bad.list", row.names = F, col.names = F, quote = F)
# $ cat  # to show up the file
system("tar -czvf run04_david_bad_xlsx.tar.gz -T david_bad.list ") # tar according this table -T

david_good_list <- paste("LUAD_survivalAnalysis_marginS_", david_good$gene_id, ".xlsx", sep="")
write.table(david_good_list, file="david_good.list", row.names = F, col.names = F, quote = F)
system("tar -czvf run04_david_good_xlsx.tar.gz -T david_good.list") # tar according this table -T

system("ls -al *.gz LUAD_OS_marginS_candidates_Venn.xlsx")
setwd(cur_wd)
# (OK done)
#}

# venn for adeno genes or lung genes
# Peter: lncRNA or ncRNA
#source("https://bioconductor.org/biocLite.R")
#biocLite("GDCRNATools")
install_github("Jialab-UCR/GDCRNATools")

#{from pubmed.mineR text mining
venn_lung_adeno <- list(lung_genes$V1[!duplicated(lung_genes)], adeno_genes$V1[!duplicated(adeno_genes)])
names_lung_adeno <- c("associated_with_lung", "associated_with_adenocarcinoma")
venn(venn_lung_adeno, snames=names_lung_adeno,
     ilabels = T, counts = T, ellipse = FALSE, zcolor = "forestgreen, darkolivegreen1", opacity = 0.6, size = 15, cexil = 3, cexsn = 0.85, borders = TRUE)
#      predefined colors if "style"
title <- paste(c("Text mining (Pubmed) of", "LUAD survival analysis", paste(length(candidates_Bonferroni_pvalue$Gene_id), "candidate genes", sep=" ")), sep = "", collapse = "\n")
text(500,900, labels = title, cex = 0.85) # (0,0) on bottom_left corner

#}
# http://www.biomart.org a collection of database (GO gene ontology)
# install_github("grimbough/biomaRt")
library("biomaRt")
listMarts()
# 'entrezgene', "adenylate kinase" activity => GO:0004017
# https://www.ebi.ac.uk/QuickGO/term/GO:0004017
# # example https://bioconductor.org/packages/release/bioc/vignettes/biomaRt/inst/doc/biomaRt.html
# listDatasets(ensembl)
ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
aks_genes <- getBM(attributes = c('hgnc_symbol'), 
                   filters = 'go', 
                   values = 'GO:0004017', 
                   mart = ensembl)
# #  hgnc_symbol
# 1       RAD50
# 2         AK1
# 3         AK8
# 4         AK2
# 5         AK4
# 6         AK6
# 7         AK7
# 8         AK5
# 
aks_genes %in% unlist(candidates_Bonferroni_pvalue$Gene_id)
# False (none of them within)

#### [mainB process #part B] { margin free ####
# genome-wide scan for margin 0 only cohort ###
# survival_marginFree <- function() {} in TCGA_LUAD_marginFree.R
#source(paste(path_LUAD, "TCGA_LUAD_marginFree.R", sep="")) # survival_marginFree <- function() {} in TCGA_LUAD_marginFree.R
# START: set path, variables and functions...

# then run from here: start_time to end_time

start_time <- Sys.time() # counted in minutes
# file.exists or list.files()
desti_ZSWIM2 <- "ZSWIM2_free_archive.Rda" # appending summary data of error codes or P-value
if (file.exists(desti_ZSWIM2)) {load(file=desti_ZSWIM2)} else
{ZSWIM2 <- data.frame(matrix(data = NA, nrow = LUAD_n, ncol = 2))}
aa <- LUAD_n; bb<- 1

## [2018/04/05] 00:38-08:48, 8 hours, run ZZZ3(20499) to ZNF268(20023)
#main_i loops ##
# #> start_time;end_time for _marginS_
# [1] "2018-04-11 00:19:10 CST"
# [1] "2018-04-17 02:56:45 CST"

# [2018/05/21] for _marginFree_
# stop at  TSPO2 ( 18783 ) (Error in plot.new() : figure margins too large)
aa <- 1784 #GTPBP8, GLT1D1, GLCE, 
# [1] "2018-04-11 00:19:10 CST"
# [1] "2018-04-17 02:56:45 CST"
# 50: In readChar(con, 5L, useBytes = TRUE) :
# cannot open compressed file 'LUAD.mRNA.Exp.C11orf61.Fire.Rda', probable reason 'Permission denied'
# http://catcode.com/teachmod/chmod_2.html  owner group others: read/write/excute (rwx)
# chmod go=r of their permissions: - rw- r--r-- <= - rw- ------
# $ for f in LUAD.mRNA.Exp.*.Fire.Rda; do sudo chmod go=r "$f"; done
# OR (-n, dry run)
# $ rsync -avhn --delete --include '*.pdf' empty_dir/ ~/R/LUAD_Peter_survival

#check file permission:
# $ find . ! -readable
# $ sudo chmod go=rw-

for (main_i in aa:bb) {
  ZSWIM2[main_i, 2] <- survival_marginFree(whole_genome[main_i]) # codes at source("TCGA_LUAD_marginFree.R")
  # gene scan; return() at X2; for loop, we need ZSWIM2 data to be saved
  save(ZSWIM2, file=desti_ZSWIM2)
  print(paste("main_i:", main_i))
}
##}


## ## email for P-value analysisi of "ZSWIM2_archive.Rda (main_i)"####
end_time <- Sys.time()
print(paste("Duration: ", (end_time - start_time)))

# different codes: marginS (option2, (n=256)) and marginFree (option1, (n=245))
# # {
# # surgical margin status: keeping 0 and excluding + margin (as 1; n=11)
# osccCleanNA_freeMargin <- osccCleanNA[osccCleanNA$margin == 0, ] # margin==0
# #osccCleanNA <- osccCleanNA_freeMargin # n=245, LUAD s/p OP with margin free##
# # n=11, margin involved; 11/256 = 5.3% (how about it's survival impact on each individual genes in LUAD?)
# # 
# # option1) marginFree
# # margin free cohort (n=245): #### save as "LUAD_survivalAnalysis_marginS_" .Rda or .xlsx
# osccCleanNA <- osccCleanNA_freeMargin
# # 
# # option2) marginS
# # margin positive and negative cohort (n=256) #### save as "LUAD_survivalAnalysis_marginFree_" .Rda or .xlsx
# #osccCleanNA <- oscc_n256
# #}

# Go back to  [Results] section ####
#.....resume

## Post1 process _marginFree_ ####
## table1 (KM): candidate_sample(KM) ##
# #sort by mpg (ascending) and cyl (descending)
# 楊老師: Bonferroni correction (adjustment) sets the significance cut-off at α/n. of KM P-value; Cox P-value <=0.05 as well.
# https://www.stat.berkeley.edu/~mgoldman/Section0402.pdf
# newdata <- mtcars[order(mpg, -cyl),] 

# _marginFree_, "swap data"
SFree <- "_marginFree_"
load(file="LUAD_OS_marginFree_pvalueKM_candidate_cox.Rda") # as candidate_sample and candidate_cox, n_percent_Bonferroni
attach(candidate_sample)
LUAD_OS_marginS_pvalue_sorted <- candidate_sample[order(p_value, -number),] #(ascending)
detach(candidate_sample)

attach(LUAD_OS_marginS_pvalue_sorted)
# Bonferroni correction (adjustment) sets the significance cut-off at α/n
alpha_LUAD <- 0.05
Bonferroni_cutoff <- alpha_LUAD / (LUAD_n * n_percent_Bonferroni) # as 4.693073e-06
# if no Bonferroni: 
#Bonferroni_cutoff <- alpha_LUAD / 1
# Bonferroni end
LUAD_OS_marginS_pvalue005_sorted <- LUAD_OS_marginS_pvalue_sorted[which(p_value<=alpha_LUAD & !is.na(p_value)), 1:3]
detach(LUAD_OS_marginS_pvalue_sorted)  # keeping drawing
# no correction: n=6260 in _marginS_ and n=6345 in _marginFree_ of LUAD_OS_marginS_pvalue005_sorted


#  plot(OS.km, lty=1, col=c("blue","red"), sub=paste("p-value =", p_OS[j]), main=paste("OS in OSCC(n=", surv_OS$n[1]+surv_OS$n[2],")/", geneName, "cutoff=", i), ylab="Percent Survival", xlab="Days")
#  legend("topright", legend=c(paste("low(",surv_OS$n[1], ")"), paste("high(",surv_OS$n[2], ")")), lty=1:2, col=c("blue","red"))
library(stats)
library(scales)
library(minpack.lm)

attach(LUAD_OS_marginS_pvalue005_sorted)
# reg <- lm(number ~ p_value, data = LUAD_OS_marginS_pvalue005_sorted)
# abline(reg, col="blue")
# or
LUAD_OS_marginS_pvalue005_sorted$z_score <- scale(number, scale = T) # "number" frequency to z-scores ("Z" because the normal distribution is also known as the "Z distribution").
LUAD_OS_marginS_pvalue005_sorted$number_01 <- scales:::rescale(z_score, to = c(0, 1)) # rescaled to range minnew to maxnew (aka. 0 to 1 for binomial glm)

plot(p_value, number_01, type="p", ylab="Z-score", xlab="P-value from KM survival analysis \n (cutoff by Bonferroni correction)", log="x") # log scale x or y
#axis(side=3, at=c(1e-7, 1e-6, 1e-5, 0.01, 0.05)) #1=below, 2=left, 3=above and 4=right
abline(h=0.6, lty=2, col="red")
abline(v=Bonferroni_cutoff, lty=2, col="red") # *** as 4.693073e-06
# run a logistic regression model (categorical 0 vs 1)
g <- glm(number_01 ~ p_value, family=poisson, data=LUAD_OS_marginS_pvalue005_sorted)
# (in this case, generalized linear model with log link)(link = "log"), poisson distribution
curve(predict(g, data.frame(p_value = x), type="response"), add=TRUE, col="blue") # draws a curve based on prediction from logistic regression model

# x# run a non-linear regression: non-linear least squares approach (function nls in R) ; A nice feature of non-linear regression in an applied context is that the estimated parameters have a clear interpretation (Vmax in a Michaelis-Menten model is the maximum rate) which would be harder to get using linear models on transformed data for example.
# # the Levenberg-Marquardt algorithm for nonlinear regression
# # "eyeball" plot to set approximate starting values
# # https://stats.stackexchange.com/questions/160552/why-is-nls-giving-me-singular-gradient-matrix-at-initial-parameter-estimates
# # nls "singular gradient matrix at initial parameter estimates" error, using nlsLM instead
# a_start <- 0.9 #param a is the y value when x=0
# b_start <- 2*log(2)/a_start #b is the decay rate
# m <- nlsLM(number_01 ~ p_value, data=LUAD_OS_marginS_pvalue005_sorted, start=list(a=a_start, b=b_start))
# curve(predict(m, data.frame(p_value = x), type="response", col="green"), add=TRUE, col="blue") # draws a curve based on prediction from logistic regression model
# #lines(p_value, predict(m), lty=2, col="green", lwd=3)

legend("topright", legend=c(paste("LR"), paste("Cutoff")), lty=1:2, col=c("blue","red"), cex=0.7) # box and font size
# figure legend: logistic regression, LR, by Generalized linear model, glm
#detach(LUAD_OS_marginS_pvalue005_sorted)

# then...
# a candidate table, with "number of P-value" under cutoff finding (becoming z-score: bigger, much more significant cutoff "sites")
# => a kind of local minimal or global minimal of curve fitting (Levenberg-Marquardt Optimization)?
# subsetting the candidated genes table, according P-value (Bonferroni_cutoff) and Z-score
# ***after Bonferroni correction => n=33 in _marginS_
# ***after Bonferroni correction => n=33 in _marginFree_
#attach(LUAD_OS_marginS_pvalue005_sorted)
LUAD_OS_marginS_pvalue1e_6_zscore0_6 <- LUAD_OS_marginS_pvalue005_sorted[which(p_value<=Bonferroni_cutoff & z_score>=0.6), 2:5]
LUAD_OS_marginS_pvalue1e_6_zscore0_6[, 2] <- signif(LUAD_OS_marginS_pvalue1e_6_zscore0_6[, 2], 3)
detach(LUAD_OS_marginS_pvalue005_sorted) # n=17
# > colnames(LUAD_OS_marginS_pvalue1e_6_zscore0_6)
# [1] "gene_id"   "p_value"   "z_score"   "number_01"
# 
# # _marginS_ (not here)
# save(LUAD_OS_marginS_pvalue1e_6_zscore0_6, file="LUAD_OS_marginS_pvalue1e_6_zscore0_6.Rda")
# save(LUAD_OS_marginS_pvalue005_sorted, file="LUAD_OS_marginS_pvalue005_sorted.Rda")
# *** OR _marginFree_ (ok)
# 
# save as _marginFree_, "swap data"
LUAD_OS_marginFree_pvalue1e_6_zscore0_6 <- LUAD_OS_marginS_pvalue1e_6_zscore0_6
save(LUAD_OS_marginFree_pvalue1e_6_zscore0_6, file="LUAD_OS_marginFree_pvalue1e_6_zscore0_6.Rda")
LUAD_OS_marginFree_pvalue005_sorted <- LUAD_OS_marginS_pvalue005_sorted
save(LUAD_OS_marginFree_pvalue005_sorted, file="LUAD_OS_marginFree_pvalue005_sorted.Rda")






# the same varible names, "swap data" of _marginS_ and _marginFree_ [2018/05/30]

## Post2 process and add table2 (cox): candidate_cox ####
## "swap data"
# "sig" marking for significant P-value (<=0.05)
# [uni_cox_pvalue, uni_HR, uni_sig]
# [multi_cox_pvalue, multi_HR, multi_sig]
# [exp_pvalue]
# to generate LUAD_OS_marginS_pvalue005KM_sorted_pvalueCox_HR.Rda; n=??

## n=6261 (_marginS_) -> 6345 (_marginFree_)
load(file="LUAD_OS_marginFree_pvalue005_sorted.Rda") # as LUAD_OS_marginFree_pvalue005_sorted
# "swap data" for following code running
LUAD_OS_marginS_pvalue005_sorted <- LUAD_OS_marginFree_pvalue005_sorted 
#x LUAD_OS_marginS_pvalue1e_6_zscore0_6 #n=17
load(file="LUAD_OS_marginFree_pvalueKM_candidate_cox.Rda") # as candidate_cox (a list), since 2018/05/16
## "swap data" for following code running
# loading candidate_sample, candidate_cox, n_percent_Bonferroni


# new variable: KM + Cox, n=6261
LUAD_OS_marginS_pvalue005KM_sorted_pvalueCox_HR <- LUAD_OS_marginS_pvalue005_sorted
#dataframe[,"newName"] <- NA # add more named columns (NA)
LUAD_OS_marginS_pvalue005KM_sorted_pvalueCox_HR[,c("uni_HR", "uni_P_value", "uni_sig","multi_HR", "multi_P_value", "multi_sig")] <- NA
#LUAD_OS_marginS_pvalue005_sorted[,c(colnames(candidate_cox[[1]][8, c(2:4, 6:8)]))] <- NA
#colnames(LUAD_OS_marginS_pvalue005KM_sorted_pvalueCox_HR) <- c(...."uni_HR", "uni_P_value", "uni_sig",
#                                                "multi_HR", "multi_P_value", "multi_sig")

#attach(LUAD_OS_marginS_pvalue005KM_sorted_pvalueCox_HR)
#ipp <- 0
for (ip in 1:nrow(LUAD_OS_marginS_pvalue005KM_sorted_pvalueCox_HR)) {
  pos_gene <- which(whole_genome==LUAD_OS_marginS_pvalue005KM_sorted_pvalueCox_HR$gene_id[ip]) # LUAD_OS_marginS_pvalue005KM_sorted_pvalueCox_HR$gene_id
  # by geneid ?
  if (length(candidate_cox[[pos_gene]])==11) {
    # reference: candidate_cox_ip, which listing as candidate_cox
    LUAD_OS_marginS_pvalue005KM_sorted_pvalueCox_HR[ip, 6:11]  <- candidate_cox[[pos_gene]][8, c(2:4, 6:8)] # taking RNAseg: sig marked and P-values, HRs
    #   ipp <- ipp+1
    print(paste(ip, " out of ", nrow(LUAD_OS_marginS_pvalue005KM_sorted_pvalueCox_HR), " (", whole_genome[pos_gene], ")", sep=""))
    print(paste("..added...", LUAD_OS_marginS_pvalue005KM_sorted_pvalueCox_HR[ip, c(2,6)], sep=";"))
  }
  #?? error on : [2018-05-18 07:26:42] [error] handle_read_frame error: websocketpp.transport:7 (End of File)
  #"exp_pvalue" is a correlation of gene expression vs features (TNM....): 
  # <- candidate_cox[which(gene_id[ip]==whole_genome)][1:7, c(11)] # correlation
  # colnames <-  c("KM_Features", "KM_P_value",  "KM_sig")
  
}
#detach(LUAD_OS_marginS_pvalue005KM_sorted_pvalueCox_HR)
# c( "uni_Features", "uni_HR",       "uni_P_value",  "uni_sig",      "multi_Features", "multi_HR",
#"multi_P_value",  "multi_sig",      "KM_Features", "KM_P_value",  "KM_sig")
# add colname as HRs, "uni_cox"/"multi_cox", sig ... for HRs, P-values, sig of RNAseq(z-score)

##
# save table1 + table2 (ok) [2018/05/18]
# "swap data"
LUAD_OS_marginFree_pvalue005KM_sorted_pvalueCox_HR <- LUAD_OS_marginS_pvalue005KM_sorted_pvalueCox_HR
save(LUAD_OS_marginFree_pvalue005KM_sorted_pvalueCox_HR, file=paste("LUAD_OS", SFree, "pvalue005KM_sorted_pvalueCox_HR.Rda", sep=""))

#save(list(LUAD_OS_marginS_pvalue1e_6_zscore0_6, ???))
# x[i], or might be x[i:j]
# x[i, j]
# x[[i]]; x[[expr]]; it can NOT be x[[i:j]]
# x[[i, j]]
# x$a
# x$"a"




# ** Pickup all significant genes list -> LUAD_OS_marginS_THREE_pvalue005 ####
# [2018/05/30]
# > colnames(LUAD_OS_marginS_pvalue005KM_sorted_pvalueCox_HR)
# [1] "number"        "gene_id"       "p_value"      
# [4] "z_score"       "number_01"     "uni_HR"       
# [7] "uni_P_value"   "uni_sig"       "multi_HR"     
# [10] "multi_P_value" "multi_sig"
LUAD_OS_marginS_THREE_pvalue005 <- subset(LUAD_OS_marginS_pvalue005KM_sorted_pvalueCox_HR, (uni_P_value <= 0.05) & (multi_P_value <= 0.05), 
                                          select=c(number, gene_id, p_value, uni_HR, uni_P_value, multi_HR, multi_P_value))
# n=4132 => n=4464 (_marginFree_)
# "swap data" ###
# # save Bonferroni_cutoff and LUAD_OS_marginFree_THREE_pvalue005.Rda
LUAD_OS_marginFree_THREE_pvalue005 <- LUAD_OS_marginS_THREE_pvalue005
save(LUAD_OS_marginFree_THREE_pvalue005, Bonferroni_cutoff, file=paste("LUAD_OS", SFree, "THREE_pvalue005.Rda", sep=""))


# plot uni_HR, n=4464
attach(LUAD_OS_marginS_THREE_pvalue005)
plot(p_value, uni_HR, type="p", ylab="Cox Uni_HR", xlab="P-value (Bonferroni) from KM survival", log="x")
#axis(side=3, at=c(1e-7, 1e-6, 1e-5, 0.01, 0.05)) #1=below, 2=left, 3=above and 4=right
abline(h=1, lty=2, col="red")
abline(v=Bonferroni_cutoff, lty=2, col="red")
# run a logistic regression model (categorical 0 vs 1)
#g <- glm(uni_HR ~ p_value, family=poisson, data=LUAD_OS_marginS_THREE_pvalue005)
# (in this case, generalized linear model with log link)(link = "log"), poisson distribution
#curve(predict(g, data.frame(p_value = x), type="response"), add=TRUE, col="blue") # draws a curve based on prediction from logistic regression model
#legend("topright", legend=c(paste("LR"), paste("Cutoff")), lty=1:2, col=c("blue","red"), cex=0.7) # box and font size
# figure legend: logistic regression, LR, by Generalized linear model, glm
detach(LUAD_OS_marginS_THREE_pvalue005)


#plot(LUAD_OS_marginS_THREE_pvalue005$p_value,  LUAD_OS_marginS_THREE_pvalue005$uni_HR)
# plot multi_HR, n=4132
attach(LUAD_OS_marginS_THREE_pvalue005)
plot(p_value, multi_HR, type="p", ylab="Cox Multi_HR", xlab="P-value  (Bonferroni) from KM survival", log="x")
#axis(side=3, at=c(1e-7, 1e-6, 1e-5, 0.01, 0.05)) #1=below, 2=left, 3=above and 4=right
abline(h=1, lty=2, col="red")
abline(v=Bonferroni_cutoff, lty=2, col="red")
# run a logistic regression model (categorical 0 vs 1)
#g <- glm(uni_HR ~ p_value, family=poisson, data=LUAD_OS_marginS_THREE_pvalue005)
# (in this case, generalized linear model with log link)(link = "log"), poisson distribution
#curve(predict(g, data.frame(p_value = x), type="response"), add=TRUE, col="blue") # draws a curve based on prediction from logistic regression model
#legend("topright", legend=c(paste("LR"), paste("Cutoff")), lty=1:2, col=c("blue","red"), cex=0.7) # box and font size
# figure legend: logistic regression, LR, by Generalized linear model, glm
detach(LUAD_OS_marginS_THREE_pvalue005)


## Venn_marginFree_: cross matching HR>1 of Uni+Multi, HR<1 of Uni+Multi ####
#library(venn) #https://cran.r-project.org/web/packages/venn/venn.pdf
# venn(x, snames = "", ilabels = FALSE, counts = FALSE, ellipse = FALSE,
#      zcolor = "bw", opacity = 0.3, size = 15, cexil = 0.6, cexsn = 0.85,
#      borders = TRUE, ...)
#library(gplots)
# To get the list of gene present in each Venn compartment we can use the gplots package
#library(gplots) # capture the list of genes from venn

#{ pickup1 from LUAD_OS_marginS_THREE_pvalue005
#* Cox HR (>1 or) >=2.5 (bad guy genes)
LUAD_OS_marginS_uni_CoxHR2p5 <- subset(LUAD_OS_marginS_THREE_pvalue005, (p_value <= Bonferroni_cutoff) & (uni_HR >=2.5), 
                                       select=c(number, gene_id, p_value, uni_HR, uni_P_value, multi_HR, multi_P_value))
# n=10
# {
#   $`uni_Cox HR >=2.5`
#   [1] "TFF1"    "UCK2"    "RRM2"    "EIF5AL1" "DKK1"   
#   [6] "FAM111B" "ECT2"    "KIF20A"  "KIF23"   "DSCC1" 
# }


#..
LUAD_OS_marginS_multi_CoxHR2p5 <- subset(LUAD_OS_marginS_THREE_pvalue005, (p_value <= Bonferroni_cutoff) & (multi_HR >=2.5), 
                                         select=c(number, gene_id, p_value, uni_HR, uni_P_value, multi_HR, multi_P_value))
# n=4
# {
#   $`multi_Cox HR >=2.5`
#   [1] "DUSP5"  "SETD8"  "PLSCR1" "KLHDC5"
# }
#...

# venn1 diagram of HR>=2.5 of Uni & Multi ###
#for list of genes by grouping; 

venn_HR2p5 <- list(LUAD_OS_marginS_uni_CoxHR2p5$gene_id, LUAD_OS_marginS_multi_CoxHR2p5$gene_id)
names_HR2p5 <- c("uni_Cox HR >=2.5", "multi_Cox HR >=2.5")
library(gplots)
tmp <- venn(venn_HR2p5, names=names_HR2p5, show.plot=F) #library(gplots); the group count matrix alone
isect_HR2p5 <- attr(tmp, "intersections")
#isect_HR2p5
detach(package:gplots)

library(venn)
venn(venn_HR2p5, snames=names_HR2p5,
     ilabels = T, counts = T, ellipse = FALSE, zcolor = "red, deeppink", opacity = 0.6, size = 15, cexil = 0.6, cexsn = 0.85, borders = TRUE)
#      predefined colors if "style"
#http://www.stat.columbia.edu/~tzheng/files/Rcolor.pdf
# meta-language 1 0 or -
title <- paste(c("LUAD survival analysis", "(KM P-value <= 0.05)"), collapse = "\n")
#coords <- unlist(getCentroid(getZones(venn_HR2p5, snames="uni_CoxHR>=2p5, multi_CoxHR>=2p5")))
# coords[1], coords[2], 
text(500,900, labels = title, cex = 0.85) # (0,0) on bottom_left corner
# n=6
#{
# # $`uni_Cox HR >=2.5:multi_Cox HR >=2.5`
# [1] "PLCD3"    "CYTSB"    "TFAP2A"   "KIF14"    "C9orf140"
# [6] "RG9MTD2" 
#}

#https://stackoverflow.com/questions/43324180/adding-legend-to-venn-diagram
#legend("top", legend=c("B:multi_Cox HR >=2.5", "A:uni_Cox HR >=2.5"), lty=1:2, col=c("blue","red"), cex=0.7) # box and font size
# figure legend: 



#===
#{ pickup2 from LUAD_OS_marginS_THREE_pvalue005
#* Cox HR <0.4 (good guy genes)

#...
LUAD_OS_marginS_uni_CoxHR0p5 <- subset(LUAD_OS_marginS_THREE_pvalue005, (p_value <= Bonferroni_cutoff) & (uni_HR <0.4), 
                                       select=c(number, gene_id, p_value, uni_HR, uni_P_value, multi_HR, multi_P_value))
# n=5
# {
#  # $`uni_Cox HR <0.4`
# [1] "DBP"     "TXNDC11" "ZNF709"  "PDK2"    "PLEKHB1"
# }

# ...
LUAD_OS_marginS_multi_CoxHR0p5 <- subset(LUAD_OS_marginS_THREE_pvalue005, (p_value <= Bonferroni_cutoff) & (multi_HR <0.4), 
                                         select=c(number, gene_id, p_value, uni_HR, uni_P_value, multi_HR, multi_P_value))
# n=9
# {
#   # $`multi_Cox HR <0.4`
#   # [1] "PXMP4"     "FBP1"      "PLD3"      "ACAD8"     "ZNF14"    
#   # [6] "LIFR"      "LOC151174" "NFATC2"    "MEX3A"    
# }


# venn2 diagram of HR < 0.4 of Uni & Multi ###
#for list of genes by grouping; 
venn_HR0p5 <- list(LUAD_OS_marginS_uni_CoxHR0p5$gene_id, LUAD_OS_marginS_multi_CoxHR0p5$gene_id)
names_HR0p5 <- c("uni_Cox HR < 0.4", "multi_Cox HR < 0.4")
library(gplots)
tmp <- venn(venn_HR0p5, names=names_HR0p5, show.plot=F) #library(gplots); the group count matrix alone
isect_HR0p5 <- attr(tmp, "intersections")
#isect_HR0p5
detach(package:gplots)

library(venn)
venn(venn_HR0p5, snames=names_HR0p5,
     ilabels = T, counts = T, ellipse = FALSE,  zcolor = "forestgreen, darkolivegreen1", opacity = 0.6, size = 15, cexil = 0.6, cexsn = 0.85, borders = TRUE)
#      predefined colors if "style"
# meta-language 1 0 or -, https://cran.r-project.org/web/packages/gplots/vignettes/venn.pdf
title <- paste(c("LUAD survival analysis", "(KM P-value <= 0.05)"), collapse = "\n")
text(500,900, labels = title, cex = 0.85) # (0,0) on bottom_left corner

#{ intersection of them (HR <0.4), n=8
# $`uni_Cox HR <0.4:multi_Cox HR <0.4`
# [1] "CRHR2"        "MYLIP"        "ZNF682"       "SLC11A2"
# [5] "NUP210L"      "LOC284440"    "LOC100130093" "ZKSCAN4"

# #library(VennDiagram)
# VENN.LIST <- list(LUAD_OS_marginS_uni_CoxHR0p5$gene_id, LUAD_OS_marginS_multi_CoxHR0p5$gene_id)
# #venn.plot <- venn.diagram(VENN.LIST , NULL, fill=c("darkmagenta", "darkblue"), alpha=c(0.5,0.5), cex = 2, cat.fontface=4, category.names=c("A", "B"), main="LUAD: Cox HR <0.5")
# # To plot the venn diagram we will use the grid.draw() function to plot the venn diagram
# #grid.draw(venn.plot)
# # We can summarize the contents of each venn compartment, as follows:
# # in 1) ConditionA only, 2) ConditionB only, 3) ConditionA & ConditionB
# lapply(isect_HR0p5, head)
#}


#} venn the end




# Excluding the LUAD cancer driver genes list??? ##
#BioXpress*.csv # data from https://hive.biochemistry.gwu.edu/cgi-bin/prd/bioxpress/servlet.cgi



#*** swap data codes from here ####
## Export r2excel and .Rda ####
# _marginFree_
# sink() for .xlsx export as well :-) https://stackoverflow.com/questions/34038041/how-to-merge-multiple-data-frame-into-one-table-and-export-to-excel
# [2018/05/29] => refinement of LUAD_OS_marginS_candidates_Venn.xlsx: show up Bonferroni_cutoff and 5e-6 (expression style).
# 
load(file="LUAD_OS_marginFree_pvalue1e_6_zscore0_6.Rda") # in LUAD_OS_marginFree_pvalue1e_6_zscore0_6, cutoff by Bonferroni_cutoff
load(file="LUAD_OS_marginFree_pvalueKM_candidate_cox.Rda") # in candidate_cox (a list), candidate_sample, candidate_cox, n_percent_Bonferroni
# "swap data"
LUAD_OS_marginS_pvalue1e_6_zscore0_6 <- LUAD_OS_marginFree_pvalue1e_6_zscore0_6



library("xlsx")
library("r2excel")
# Create an Excel workbook. Both .xls and .xlsx file formats can be used.
filenamex <- paste("LUAD_OS", SFree, "candidates_Venn", ".xlsx", sep = "") # "LUAD_OS_marginFree_candidates_Venn.xlsx"
wb <- createWorkbook(type="xlsx")

# Create a sheet in that workbook
sheet <- xlsx::createSheet(wb, sheetName = paste("Survival_candidates"))
# [add data row by row, start from column 2]
#+++++++++++++++++++++++++++++++
## Add paragraph : Author
library("tis") # by Brian Salzer
# today(), arg must be ti, tis, ts, tif, or tifName
author <- paste("Reported by Tex Li-Hsing Chi. \n",
                "tex@gate.sinica.edu.tw \n", SFree, "\n", Sys.Date(), sep="")
xlsx.addParagraph(wb, sheet, value=author, isItalic=TRUE, colSpan=5, 
                  rowSpan=4, fontColor="darkgray", fontSize=24)
xlsx.addLineBreak(sheet, 3)
# header
xlsx.addHeader(wb, sheet, value=paste("Table 1. The candiate genes expressed in ", TCGA_cohort,
                                      " (ranking by KM P-value, selected by Bonferroni cutoff, ", signif(Bonferroni_cutoff, 3), ") ", "\n", "n= ", nrow(LUAD_OS_marginS_pvalue1e_6_zscore0_6), sep=""),
               level=5, color="black", underline=0)
#xlsx.addHeader(wb, sheet, value=paste("Cutoff at ", round(cutoff1, 3), " (", percent(surv_OS1$n[1]/(surv_OS1$n[1]+surv_OS1$n[2])), ")", sep = ""),
#               level=5, color="red", underline=0) # total n is taken from surv_OS1$n

xlsx.addLineBreak(sheet, 1) # add one blank line

#xlsx.addTable(wb, sheet, data = t(data.frame(c(paste(geneName, "expression"), "", paste(geneName, "expression"), "", "(Optimised)"))), fontSize=12, startCol=4,
#              fontColor="darkblue", row.names = F, col.names = F) #, colSpan=1, rowSpan=1)

# a candidate genes table
# > colnames(LUAD_OS_marginS_pvalue1e_6_zscore0_6)
# [1] "gene_id"   "p_value"   "z_score"   "number_01"
colnames(LUAD_OS_marginS_pvalue1e_6_zscore0_6) <- c("Gene_id",   "P_value",   "z_score_raw",   "Z_score") # Z_score is rescaled as 0-1
# Bonferroni_cutoff; [, c(1,2,4)]
# in scientific notation: formatC of [, c(2)]
attach(LUAD_OS_marginS_pvalue1e_6_zscore0_6)
candidates_Bonferroni_pvalue <- LUAD_OS_marginS_pvalue1e_6_zscore0_6[which(P_value<=Bonferroni_cutoff), c(1,2,4)]
detach(LUAD_OS_marginS_pvalue1e_6_zscore0_6)

attach(candidates_Bonferroni_pvalue) # removal of as.factor (P_value)
candidates_Bonferroni_pvalue$P_value <- formatC(as.numeric(as.character(P_value)), format = "e", digits = 2)
candidates_Bonferroni_pvalue$Z_score <- signif(Z_score, digits=4)
candidates_Bonferroni_pvalue <- candidates_Bonferroni_pvalue[order(P_value, -Z_score), ] #sorting by order(ascending)
detach(candidates_Bonferroni_pvalue)

# content of Table 1
xlsx.addTable(wb, sheet, data = candidates_Bonferroni_pvalue, startCol=2,
              fontColor="darkblue", fontSize=12,
              rowFill=c("white", "white"), row.names = TRUE)

xlsx.addLineBreak(sheet, 5)  # add two blank lines

# Export z-score/P-Value plot by plotFunction() ####
xlsx.addLineBreak(sheet, 1)

# Define function
xlsx.addPlot.candidates<-function( wb, sheet, startRow=NULL,startCol=2,
                                   width=480, height=480,... )
{ # plot of z-score vs p-value digram from "the summary"
  
  png(filename = "plot.png", width = width, height = height,...)
  #{
  # plot(OS.km, lty=1, xscale=12, xmax=60, col=c("blue","red"), sub=paste("P-value =", p_OS1), main=paste("OS in TCGA ", TCGA_cohort, "(n=", surv_OS1$n[1]+surv_OS1$n[2],")/", geneName), ylab="Percent Survival", xlab="Years")
  # legend("topright", legend=c(paste("low(",surv_OS1$n[1], ")"), paste("high(",surv_OS1$n[2], ")")), lty=1:2, col=c("blue","red"))
  load(file="LUAD_OS_marginFree_pvalue005_sorted.Rda") # as LUAD_OS_marginFree_pvalue005_sorted
  # "swap data" for following code running
  LUAD_OS_marginS_pvalue005_sorted <- LUAD_OS_marginFree_pvalue005_sorted
  attach(LUAD_OS_marginS_pvalue005_sorted)
  plot(p_value, number_01, type="p", ylab="Z-score", xlab="P-value from KM survival analysis \n (cutoff by Bonferroni correction)", log="x") # log scale x or y
  #axis(side=3, at=c(1e-7, 1e-6, 1e-5, 0.01, 0.05)) #1=below, 2=left, 3=above and 4=right
  abline(h=0.6, lty=2, col="green")
  abline(v=Bonferroni_cutoff, lty=2, col="red") # *** as 4.693073e-06
  # run a logistic regression model (categorical 0 vs 1)
  g <- glm(number_01 ~ p_value, family=poisson, data=LUAD_OS_marginS_pvalue005_sorted)
  # (in this case, generalized linear model with log link)(link = "log"), poisson distribution
  curve(predict(g, data.frame(p_value = x), type="response"), add=TRUE, col="blue") # draws a curve based on prediction from logistic regression model
  legend("topright", legend=c(paste("LR"), paste("Cutoff Z"), paste("Cutoff B")), lty=1:2, col=c("blue","green","red"), cex=0.9) # box and font size
  # figure legend: logistic regression, LR, by Generalized linear model, glm
  detach(LUAD_OS_marginS_pvalue005_sorted)
  #}
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
} # Define function

# calling
currentRow <- xlsx.addPlot.candidates(wb, sheet) # startRow + a plot
print(paste("The z-score summary plot: ", filenamex, " was exported successfully."))



#..
# Export Venn1 by plotFunction() ####
xlsx.addLineBreak(sheet, 1)

# Define function
xlsx.addPlot.venn<-function( wb, sheet, startRow=NULL,startCol=2,
                             width=480, height=480,... )
{ # plot of venn digram"
  
  png(filename = "plot1.png", width = width, height = height,...)
  #{
  venn(venn_HR2p5, snames=names_HR2p5,
       ilabels = T, counts = T, ellipse = FALSE,  zcolor = "red, deeppink", opacity = 0.6, size = 15, cexil = 3, cexsn = 0.85, borders = TRUE)
  #      predefined colors if "style"; cexil = 0.6 (default text size of counts)
  title <- paste(c("LUAD survival analysis", "KM P-Value <= ", signif(Bonferroni_cutoff, 3)), sep = "", collapse = "\n")
  text(500,900, labels = title, cex = 0.85) # (0,0) on bottom_left corner
  
  #}
  dev.off()
  
  
  #Append plot1 to the sheet
  if(is.null(startRow)){
    rows<- getRows(sheet) #list of row object
    startRow=length(rows)+1
  } 
  # Add the file created previously
  addPicture("plot1.png", sheet=sheet,  startRow = startRow, startColumn = startCol) 
  xlsx.addLineBreak(sheet, round(width/20)+1)
  res<-file.remove("plot1.png")
  
} # end of Define function

# calling
#detach(package:gplots)
library(venn)
xlsx.addPlot.venn(wb, sheet)
print(paste("The Venn diagram (bad guy) and the candidate genes on", filenamex, ", which was exported successfully."))




# *generate (bad) prognostic features of those genes on lists in TCGA LUAD cohort ####
# store KM_sig remrark as a Byte; converted as base10 in list_KM_sigBin
# "swap data"
load(file="LUAD_OS_marginFree_pvalueKM_candidate_cox.Rda") # as candidate_cox (a list), since 2018/05/16; with KM_sig, Remark
# candidate_sample, candidate_cox, n_percent_Bonferroni
load(file="LUAD_OS_marginFree_THREE_pvalue005.Rda") # as as LUAD_OS_marginFree_THREE_pvalue005, Bonferroni_cutoff
LUAD_OS_marginS_THREE_pvalue005 <- LUAD_OS_marginFree_THREE_pvalue005



#.. bad guy gene candidate
# => see "remark=1"; "which" or NOT "which", the order is wrong.
# isect_HR2p5 is a list
#geneid_bad_uni_HR2p5 <- c(isect_HR2p5$`A:B`, isect_HR2p5$`A`) # from venn_HR2p5: A + A:B = group uni_
geneid_bad_uni_HR2p5 <- c(isect_HR2p5[[3]], isect_HR2p5[[1]])
# whole_genome[which(whole_genome %in% geneid_bad_uni_HR2p5)]
# [1] "DKK1"    "DSCC1"   "ECT2"    "EIF5AL1" "FAM111B" "KIF20A"  "KIF23"  
# [8] "RRM2"    "TFF1"    "UCK2" 
candidate_bad_uni_HR2p5 <- cbind(data.frame(gene_id=geneid_bad_uni_HR2p5), LUAD_OS_marginS_THREE_pvalue005[which(LUAD_OS_marginS_THREE_pvalue005$gene_id %in% geneid_bad_uni_HR2p5), 3:7])
#as.binary(candidate_cox[which(whole_genome %in% geneid_bad_uni_HR2p5)]$KM_sig, n=7, logic=T) # store KM_sig remrark as a Byte
# can NOT use "whole_genome" any more
#x list_KM_sig <- candidate_cox[which(whole_genome %in% geneid_bad_uni_HR2p5)] 
list_KM_sig <- candidate_cox[as.numeric(rownames(candidate_bad_uni_HR2p5))] # a list (rownames has whole_genome's right position)
# km <- function(x) {print(list_KM_sig[[x]][8,9])}; lapply(c(1:10), km)
list_KM_sigBin <- data.frame(matrix(NA, nrow = length(geneid_bad_uni_HR2p5), ncol = nrow(candidate_cox[[1]])-1))
colnames(list_KM_sigBin) <- c(rownames(candidate_cox[[1]])[-8])
for (ikm in 1:length(list_KM_sig)) 
{
  #list_KM_sigBin[ikm, 2:nrow(candidate_cox[[1]])] <- 
  list_KM_sigBin[ikm,] <- as.numeric(list_KM_sig[[ikm]]$KM_sig[-8]) # or as .factor
}
candidate_bad_uni_HR2p5 <- cbind(candidate_bad_uni_HR2p5, list_KM_sigBin)
# ex. "RRM2" has good impact on TN and pathology when expression is high

# => see "remark=1"                                                # number, gene_id, p_value, uni_HR, uni_P_value, multi_HR, multi_P_value])
# # isect_HR2p5 is a list
geneid_bad_multi_HR2p5 <- c(isect_HR2p5[[3]], isect_HR2p5[[2]]) # from venn_HR2p5
candidate_bad_multi_HR2p5 <- cbind(data.frame(gene_id=geneid_bad_multi_HR2p5), LUAD_OS_marginS_THREE_pvalue005[(LUAD_OS_marginS_THREE_pvalue005$gene_id %in% geneid_bad_multi_HR2p5), 3:7])
# as.binary
list_KM_sig <- candidate_cox[as.numeric(rownames(candidate_bad_multi_HR2p5))] # a list
list_KM_sigBin <- data.frame(matrix(NA, nrow = length(geneid_bad_multi_HR2p5), ncol = nrow(candidate_cox[[1]])-1))
colnames(list_KM_sigBin) <- c(rownames(candidate_cox[[1]])[-8])
for (ikm in 1:length(list_KM_sig)) 
{
  #list_KM_sigBin[ikm, 2:nrow(candidate_cox[[1]])] <- 
  list_KM_sigBin[ikm,] <- as.numeric(list_KM_sig[[ikm]]$KM_sig[-8]) # or as .factor
}
candidate_bad_multi_HR2p5 <- cbind(candidate_bad_multi_HR2p5, list_KM_sigBin)
# ex. "DUSP5" has good impact on M stage

# => do NOT need to see "remark=1"
geneid_bad_unimulti_HR2p5 <- isect_HR2p5[[3]] #$`A:B`
candidate_bad_unimulti_HR2p5 <- cbind(data.frame(gene_id=geneid_bad_unimulti_HR2p5), LUAD_OS_marginS_THREE_pvalue005[(LUAD_OS_marginS_THREE_pvalue005$gene_id %in% geneid_bad_unimulti_HR2p5), 3:7])


# export tables
xlsx.addLineBreak(sheet, 4) # 
# header of candidate_bad_unimulti_HR2p5
xlsx.addHeader(wb, sheet, value=paste("Table 2A. The ", nrow(candidate_bad_unimulti_HR2p5), " consensus candiate genes (uni_ & multi_CoxHR>2.5) in ", TCGA_cohort, "(bad guy)"),
               level=5, color="black", underline=0)
xlsx.addLineBreak(sheet, 1) # add one blank line
candidate_bad_unimulti_HR2p5[,2] <- formatC(candidate_bad_unimulti_HR2p5[,2], format = "e", digits = 2) #as.numeric(as.character
xlsx.addTable(wb, sheet, data = candidate_bad_unimulti_HR2p5[, c(1:6)], startCol=2,
              fontColor="darkblue", fontSize=12,
              rowFill=c("white", "white"), row.names = F)
xlsx.addLineBreak(sheet, 2)  # add two blank lines

xlsx.addLineBreak(sheet, 3)
# header of candidate_bad_uni_HR2p5
xlsx.addHeader(wb, sheet, value=paste("Table 2B. The candiate genes (uni_CoxHR>2.5) in ", TCGA_cohort, "(bad guy)"),
               level=5, color="black", underline=0)
xlsx.addLineBreak(sheet, 1) # add one blank line
candidate_bad_uni_HR2p5[,2] <- formatC(candidate_bad_uni_HR2p5[,2], format = "e", digits = 2) #as.numeric(as.character
xlsx.addTable(wb, sheet, data = candidate_bad_uni_HR2p5[, c(1:6, 9:13)], startCol=2,
              fontColor="darkblue", fontSize=12,
              rowFill=c("white", "white"), row.names = F)
xlsx.addLineBreak(sheet, 2)  # add two blank lines

xlsx.addLineBreak(sheet, 3)
# header of candidate_bad_multi_HR2p5
xlsx.addHeader(wb, sheet, value=paste("Table 2C. The candiate genes (multi_CoxHR>2.5) in ", TCGA_cohort, "(bad guy)"),
               level=5, color="black", underline=0)
xlsx.addLineBreak(sheet, 1) # add one blank line
candidate_bad_multi_HR2p5[,2] <- formatC(candidate_bad_multi_HR2p5[,2], format = "e", digits = 2) #as.numeric(as.character
xlsx.addTable(wb, sheet, data = candidate_bad_multi_HR2p5[, c(1:6, 9:13)], startCol=2,
              fontColor="darkblue", fontSize=12,
              rowFill=c("white", "white"), row.names = F)
xlsx.addLineBreak(sheet, 2)  # add two blank lines


# Export Venn2 by plotFunction() ####
xlsx.addLineBreak(sheet, 1)

# Define function
xlsx.addPlot.venn<-function( wb, sheet, startRow=NULL,startCol=2,
                             width=480, height=480,... )
{ # plot of venn digram"
  
  
  png(filename = "plot2.png", width = width, height = height,...)
  #{
  venn(venn_HR0p5, snames=names_HR0p5,
       ilabels = T, counts = T, ellipse = FALSE, zcolor = "forestgreen, darkolivegreen1", opacity = 0.6, size = 15, cexil = 3, cexsn = 0.85, borders = TRUE)
  #      predefined colors if "style"
  title <- paste(c("LUAD survival analysis", "KM P-Value <= ", signif(Bonferroni_cutoff, 3)), sep = "", collapse = "\n")
  text(500,900, labels = title, cex = 0.85) # (0,0) on bottom_left corner
  
  #}
  dev.off()
  
  # append plot2
  if(is.null(startRow)){
    rows<- getRows(sheet) #list of row object
    startRow=length(rows)+1
  } 
  addPicture("plot2.png", sheet=sheet,  startRow = startRow , startColumn = startCol) 
  # jump to next venn ?
  xlsx.addLineBreak(sheet, round(width/20)+1)
  res<-file.remove("plot2.png")
} # end of Define function

# calling
#detach(package:gplots)
library(venn)
xlsx.addPlot.venn(wb, sheet)
print(paste("The Venn diagram2 and the candidate genes: ", filenamex, " were exported successfully."))


# *generate (good) prognostic features of those genes on lists in TCGA LUAD cohort ####

# in case, there is only TWO group: isect_HR0p5[[2]] without [[3]]
# > isect_HR0p5
# $`uni_Cox HR < 0.4`[[1]]A
# [1] "DBP"     "C3orf18" "ZNF709"  "MYOZ1"   "PDIK1L"  "DIO1"   
# 
# $`uni_Cox HR < 0.4:multi_Cox HR < 0.4`[[2]]A:B
# [1] "CRHR2"     "NUP210L"   "MYLIP"     "SLC11A2"   "ZNF682"   
# [6] "LOC284440" "PXMP4"     "PDK2" 

# => see "remark=1" 
geneid_good_uni_HR0p5 <- c(isect_HR0p5[[2]], isect_HR0p5[[1]]) #$`A`
candidate_good_uni_HR0p5 <- cbind(data.frame(gene_id=geneid_good_uni_HR0p5), LUAD_OS_marginS_THREE_pvalue005[(LUAD_OS_marginS_THREE_pvalue005$gene_id %in% geneid_good_uni_HR0p5), 3:7])
# number, gene_id, p_value, uni_HR, uni_P_value, multi_HR, multi_P_value])
# the list (order) of gene_id has error ???
list_KM_sig <- candidate_cox[as.numeric(rownames(candidate_good_uni_HR0p5))] # a list
list_KM_sigBin <- data.frame(matrix(NA, nrow = length(geneid_good_uni_HR0p5), ncol = nrow(candidate_cox[[1]])-1))
colnames(list_KM_sigBin) <- c(rownames(candidate_cox[[1]])[-8])
for (ikm in 1:length(list_KM_sig)) 
{
  #list_KM_sigBin[ikm, 2:nrow(candidate_cox[[1]])] <- 
  list_KM_sigBin[ikm,] <- as.numeric(list_KM_sig[[ikm]]$KM_sig[-8]) # or as .factor
}
candidate_good_uni_HR0p5 <- cbind(candidate_good_uni_HR0p5, list_KM_sigBin)
# ex."TXNDC11" higher expression has better prognosis impact as well as smaller T and less N and lower stage

# => see "remark=1" 
geneid_good_multi_HR0p5 <- c(isect_HR0p5[[2]]) #, isect_HR0p5[[2]]) #$`B`
candidate_good_multi_HR0p5 <- cbind(data.frame(gene_id=geneid_good_multi_HR0p5), LUAD_OS_marginS_THREE_pvalue005[(LUAD_OS_marginS_THREE_pvalue005$gene_id %in% geneid_good_multi_HR0p5), 3:7])
# as.binary
list_KM_sig <- candidate_cox[as.numeric(rownames(candidate_good_multi_HR0p5))] # a list
list_KM_sigBin <- data.frame(matrix(NA, nrow = length(geneid_good_multi_HR0p5), ncol = nrow(candidate_cox[[1]])-1))
colnames(list_KM_sigBin) <- c(rownames(candidate_cox[[1]])[-8])
for (ikm in 1:length(list_KM_sig)) 
{
  #list_KM_sigBin[ikm, 2:nrow(candidate_cox[[1]])] <- 
  list_KM_sigBin[ikm,] <- as.numeric(list_KM_sig[[ikm]]$KM_sig[-8]) # or as .factor
}
candidate_good_multi_HR0p5 <- cbind(candidate_good_multi_HR0p5, list_KM_sigBin)
# ex."PLD3" has less N


# => do NOT need to see "remark=1"
geneid_good_unimulti_HR0p5 <- isect_HR0p5[[2]] #$`A:B`
candidate_good_unimulti_HR0p5 <- cbind(data.frame(gene_id=geneid_good_unimulti_HR0p5), LUAD_OS_marginS_THREE_pvalue005[(LUAD_OS_marginS_THREE_pvalue005$gene_id %in% geneid_good_unimulti_HR0p5), 3:7])
# *** using table 2 to interpretate the impact on survival by DEG of this candidate gene list 
# c("CRHR2", "MYLIP", "NUP210L", "ZKSCAN4")
# => see "remark=1" on LUAD_survivalAnalysis_marginS_CRHR2.xlsx => higher CRHR2 expression is associated with less LN metastaese
# => see "remark=1" on LUAD_survivalAnalysis_marginS_CRHR2.xlsx => higher MYLIP expression is associated with smaller tumor size (T)
# => see "remark=1" on LUAD_survivalAnalysis_marginS_CRHR2.xlsx => higher NUP210L expression is associated with less LN metastaese
# => see "remark=1" on LUAD_survivalAnalysis_marginS_CRHR2.xlsx => higher ZKSCAN4 expression is associated with smaller tumor size (T)

# export tables
xlsx.addLineBreak(sheet, 4) # 
# header of candidate_good_unimulti_HR0p5
xlsx.addHeader(wb, sheet, value=paste("Table 3A. The ", nrow(candidate_good_unimulti_HR0p5), " consensus candiate genes (uni_ & multi_CoxHR<0.5) in ", TCGA_cohort, "(good guy)"),
               level=5, color="black", underline=0)
xlsx.addLineBreak(sheet, 1) # add one blank line
candidate_good_unimulti_HR0p5[,2] <- formatC(candidate_good_unimulti_HR0p5[,2], format = "e", digits = 2) #as.numeric(as.character
xlsx.addTable(wb, sheet, data = candidate_good_unimulti_HR0p5[, c(1:6)], startCol=2,
              fontColor="darkblue", fontSize=12,
              rowFill=c("white", "white"), row.names = F)
xlsx.addLineBreak(sheet, 2)  # add two blank lines

xlsx.addLineBreak(sheet, 3)
# header of candidate_good_uni_HR0p5
xlsx.addHeader(wb, sheet, value=paste("Table 3B. The candiate genes (uni_CoxHR<0.5) in ", TCGA_cohort, "(good guy)"),
               level=5, color="black", underline=0)
xlsx.addLineBreak(sheet, 1) # add one blank line
candidate_good_uni_HR0p5[,2] <- formatC(candidate_good_uni_HR0p5[,2], format = "e", digits = 2) #as.numeric(as.character
xlsx.addTable(wb, sheet, data = candidate_good_uni_HR0p5[, c(1:6, 9:13)], startCol=2,
              fontColor="darkblue", fontSize=12,
              rowFill=c("white", "white"), row.names = F)
xlsx.addLineBreak(sheet, 2)  # add two blank lines

xlsx.addLineBreak(sheet, 3)
# header of candidate_good_multi_HR0p5
xlsx.addHeader(wb, sheet, value=paste("Table 3C. The candiate genes (multi_CoxHR<0.5) in ", TCGA_cohort, "(good guy)"),
               level=5, color="black", underline=0)
xlsx.addLineBreak(sheet, 1) # add one blank line
candidate_good_multi_HR0p5[,2] <- formatC(candidate_good_multi_HR0p5[,2], format = "e", digits = 2) #as.numeric(as.character
xlsx.addTable(wb, sheet, data = candidate_good_multi_HR0p5[, c(1:6, 9:13)], startCol=2,
              fontColor="darkblue", fontSize=12,
              rowFill=c("white", "white"), row.names = F)
xlsx.addLineBreak(sheet, 2)  # add two blank lines


##..
# save the workbook to an Excel file and write the file to disk.####
xlsx::saveWorkbook(wb, filenamex)
#xlsx.openFile(filenamex) # open file to review


# the END of R2Excel ###
#}

# https://david.ncifcrf.gov/conversion.jsp?VFROM=NA DAVID pathway analysis
david_bad <- merge(candidate_bad_uni_HR2p5, candidate_bad_multi_HR2p5, by="gene_id", all=TRUE) #joint by union
# ? (Family With Sequence Similarity 111 Member B) => converting to Entrez Gene ID(374393)
# > david_bad$gene_id

# 
david_good <- merge(candidate_good_uni_HR0p5, candidate_good_multi_HR0p5, by="gene_id", all=TRUE) #joint by union






# scp *.xlsx files to tar ####
# failure to apply "sudo" in R ;-)
isect_list <- unlist(c(isect_HR2p5, isect_HR0p5)) # 42 genes
scp_list <- paste("LUAD_survivalAnalysis_marginS_", isect_list, ".xlsx", sep="")
#setwd(path_LUAD)
#if (scp_list %in% dir()) {
# copy file to ./xlsx directory
# dir.create("./xlsx")
# file.copy(from, to, overwrite = recursive, recursive = FALSE,
# copy.mode = TRUE, copy.date = FALSE) showWarnings = TRUE,
file.copy(from=scp_list, to="./xlsx", overwrite = FALSE,
          copy.mode = TRUE, copy.date = FALSE)

#}
write.csv(scp_list, file="scp_list.txt", row.names =F)

# $ sudo xargs -a scp_list.txt -i cp {} ./xlsx
list.files("./xlsx") # ok
#
#
#

## LUAD_survivalAnalysis_marginFree_ vs marginS ####
## analysis of gene list by DAVID
source("https://bioconductor.org/biocLite.R")
biocLite("RDAVIDWebService")
browseVignettes("RDAVIDWebService")

# https://bioconductor.org/packages/release/bioc/html/biomaRt.html
source("https://bioconductor.org/biocLite.R")
biocLite("biomaRt")
## 
# text mining
install_git("cran/pubmed.mineR")
library(pubmed.mineR)

## 
# save tables and survival P-value from  cutoff finder #
# # in OS
# OS <- data.frame(cases_OS, as.numeric(p_OS))
# # removal of duplicated item
# OS <- OS[!duplicated(OS),]
# colnames(OS)[2] <- "p_OS"
# #View(OS[OS$p_OS <= 0.05,1:2])
# # *** appending survival
# debug: Error in if (verbose) cat("Loading objects:\n") : argument is of length zero
# append.Rda(OS[OS$p_OS <= 0.05,1:2], paste("LUAD_survivalAnalysis_marginS_", geneName, ".Rda"))
append.Rda(OS[OS$p_OS <= 0.05,1:2], paste("LUAD_survivalAnalysis_marginFree_", geneName, ".Rda"))

# # in RFS
# RFS <- data.frame(cases_RFS, as.numeric(p_RFS))
# # removal of duplicated item
# RFS <- RFS[!duplicated(RFS),]
# colnames(RFS)[2] <- "p_RFS"
# #View(RFS[RFS$p_RFS <= 0.05,1:2])
# # *** appending survival
# append.Rda(RFS[RFS$p_RFS <= 0.05,1:2], paste("LUAD_survivalAnalysis_marginS_", geneName, ".Rda"))
# append.Rda(RFS[RFS$p_RFS <= 0.05,1:2], paste("LUAD_survivalAnalysis_marginFree_", geneName, ".Rda"))

##
# a new comparison table for impact genes ####
paste("LUAD_survivalAnalysis_marginS_", geneName, ".Rda")
# a new comparison for margin issue
paste("LUAD_survivalAnalysis_marginFree_", geneName, ".Rda")

##
oby <- read.table(header=TRUE, text='
                  Var1 Var2 Freq
                  0      1    36
                  1      1   91
                  0     2    41
                  1       2   88
                  ')
# x oby <- read.table(stdin(), header=TRUE) 
# 
# 
# 
# ----------------------------##
# Spared codes ***
# # 
#ddply(seq(cutoff[1], cutoff[2], length.out = find_repeat),  , failwith(NA, fun(x){...}, quiet = TRUE), .parallel = TRUE)
# error handle; failwith(NA, fun, quiet = TRUE)

# #debug; last task was executed
# dopar <- function(main_i){
#   print(main_i)
#   print(whole_genome[main_i])
# }
# #pforeach(main_i=aa:bb, .cores=2, .errorhandling = c("stop")) ({
# # dopar(main_i)
# # })
# foreach(main_i=aa:bb, .errorhandling = c("stop"), .combine = "c") %dopar% {
#   dopar(main_i)
# }
# #debug
# ## doParallel by Parallelized foreach: pforeach (..., .errorhandling = c("pass", "stop", "remove"))


