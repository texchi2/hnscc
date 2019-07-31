# tobacco_exposure
# => ***** 只要這一個 KM plot 即可：run n=521 cohort, survival analysis without RNAseq, and categorized by tobacco_exposure high/low.
# => table 3, table 4: hazard ratios of features (including tobacco_exposure, without RNAseq) (high/low)
# By Cox HR
## START:  ####
# since [2019/06/06][2019/07/31] tobacco_exposure
TCGA_cohort <- "HNSCC" # cancer type: LUAD or HNSCC, defined
path_cohort <- "~/R/HNSCC_Tex_survival/hnscc_github" # under rstudio-server on GCP instance 4 as well as local Macbook Pro 2009
setwd(path_cohort) # set the working directory under rstudio-server on HNSCC, GCP
load(file="whole_genome.Rda") # the name list of protein coding genome
LUAD_n <- length(whole_genome) # n=20500

source(file=file.path(path_cohort, "TCGA_HNSCC_marginS.R")) # survival_marginS <- function() {} in TCGA_HNSCC_marginS.R
source(file=file.path(path_cohort, "TCGA_HNSCC_marginFree.R")) # survival_marginFree <- function() {} in TCGA_HNSCC_marginFree.R
source(file=file.path(path_cohort, "cutofFinder_func_HNSCC.R")) # cutofFinder_func <- function(geneName) {} in cutofFinder_func.R

library(gmailr) # notify me the script
library(scales)
library(plyr); library(dplyr) #ddply()
library(magrittr)
recipient <- "texchi2@gmail.com"
sender <- "texchi2@gmail.com"

### $setup global variables...(including "margin") ####
library(psych) # for describe()
library(survival)
# # pathologic_T => clinical_T and so on...
coln_osccT <- c("Unique.ID","Gender","ageDx",
                "clinical_T","clinical_N",
                "clinical_M","stage", "margin", "tobacco",
                "OS_IND","OS..months._from.biopsy",
                "RFS_IND", "RFS..months._from.op", "H.score_T", paste("PMM1", "_median", sep=""))


# for tableChi1 (table 2)
contLowHighN <- c("Features",	"Low", "(%)",	"High",	"(%)", "Case no",	"P-value", "Remark", "Matthews")

# Create the Sheet title and subtitle; underline as 0 (default value, no underline), 1 (underline with one line), 2 (underline with two lines)
# Table 3. Univariate/Multivariate Cox proportional hazards regression analyses on OS time
ciUniMti <- c("Features",	"HR",	"CI95%(L)",	"CI95%(H)",	"P-value",	"HR",	"CI95%(L)",	"CI95%(H)",	"P-value")

# with 9 features: in table 3 and table 4
featuresUni <- c("Gender",
                 "Age at diagnosis",
                 #                 "Primary site",
                 "Clinical T Status",
                 "Clinical N Status",
                 "Clinical M Status",                 
                 "Clinical Stage",
                 #"Pathologic T status",
                 #"Pathologic N status",
                 #"Pathologic M status",
                 #"Pathologic Stage",
                 "Surgical Margin status",
                 "Tobacco Exposure", 
                 #                 "Lymphovascular Invasion",
                 #                 "Perineural Invasion",
                 #                 "Extranodal spreading of neck LN",
                 #                 "SCC Histologic Grade",
                 #                 "Radiotherapy",
                 #                 "Chemotherapy",
                 paste("RNAseq(z-score)") # or "IHC score" == "z-score"
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
              "High", # risk of tobacco exposure
              "High") # PMM1 or RNAseq
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
              "Low", # risk of tobacco exposure
              "Low") # PMM1 or RNAseq
# create correlation table 1 with P-value by chisq.test
library(reshape)
library(data.table)


contingencyTCGA <- function(osccCleanNA_conTCGA, geneName) { # no more "run100"; do not need cutoff1 here
  
  #!!!!!# L <- 2 ("Gender"); R <- 9 ("tobacco) ; x 8 ("margin") # in HNSCC
  ## boundary of features column from first to last # 
  L <- which(colnames(osccCleanNA_conTCGA) == "Gender") # "left" feature
  R <- which(colnames(osccCleanNA_conTCGA) == "tobacco") # new "right" feature
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
    # t, table from col of PMM1_median vs col ii (L "Gender" to R "tobacco"); deparse the argument (colnames) by deparse.level 2
    # # {DEBUG 
    # ii<-L-1
    #     ii<-ii+1
    # # DEBUG}
    t<- NULL
    t <- table(osccCleanNA_conTCGA[,osccCleanNA_conTCGAM_pos], osccCleanNA_conTCGA[,ii], useNA = "ifany") #dnn=c(colnames(osccCleanNA_conTCGA[osccCleanNA_conTCGAM_pos])))  
    # useNA = "ifany"; "always" is for "margin" (with n=0 count)
    chiT[ii, 1] <- colnames(osccCleanNA_conTCGA[ii]) # name list of feature variables from L to R
    
    # unbalanced data issue: chisq.test -> fisher.test
    # small sample size (any cell < 5 or 10): fisher.test; https://en.wikipedia.org/wiki/Fisher%27s_exact_test
    # how about Yate's correction of chisq.test?
    # set "correct=TRUE" to turn on Yate's continuity correction of chisq.test
    # https://onlinelibrary.wiley.com/doi/full/10.1002/sim.8294#
    #h0test <- chisq.test(t)
    if (any((chisq.test(t)$expected<5)==TRUE)) { # if warming: Chi-squared approximation may be incorrect.
      h0test <- fisher.test(t) # using Fisher's exact test instead
      flag_fisher <- TRUE
      print("Fisher exact test +1")
    } else {
      h0test <- chisq.test(t) # using Chi squart test is fine
      flag_fisher <- FALSE
      print("chisq.test")
    }
    
    check_p <- h0test$p.value # retrieved in chiT$X2
    if (is.na(check_p)==T) {
      if (flag_fisher==FALSE) {
        check_p <- chisq.test(t[1:2,1:2])$p.value
      } else if (flag_fisher==TRUE) {
        check_p <- fisher.test(t[1:2,1:2])$p.value
      }
    }
    chiT[ii, 2] <- check_p
    
    # debug
    obs <- as.data.frame(chisq.test(t)$observed) # we need this observed table, even using fisher.test    
    # # error when there is only one group
    if (is.na(tryCatch(mdata <- melt(obs, id.vars=c("Var2", "Var1")), error = function(e) return(NA)))) {print(paste("skip at ii= ", ii, sep=""));return(list("skip", chiT, freq))} # back to marginS.R
    # debug for "ZSWIM2"(20482)
    
    # chisq is sum( (o-e)^2/e ); the Na should be removed, You have zero frequencies in 2 counts.
    #warnings(): In chisq.test(t) : Chi-squared approximation may be incorrect
    #=> fisher.test(a) # Fisher exact test is better in small counts. (<5) 
    #table(t); #print(c("ii=",paste(ii)))
    # t
    # 0   5   6 122 123 
    # 2   1   1   1   1 
    #print(ii) # debug
    
    # xtabs(~ PMM1 + variable[ii], data = mydata)    
    # to two-way contingency table of categorical outcome and predictors
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
    
    
    cdata <- as.data.frame(cast(mdata, Var2 ~ Var1)) 
    cdata <- cdata[, c(1:3)] # cdata should has 3 columns (NA column is removed)
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
    cdata <- cbind(cdata, chiT[ii, 1]) # add features name => # now cdata has 4 columns
    #freq[1+(ii-1)*rown,] <- cdata
    # R4> cdata
    #    Var2  0   1 chiT[ii, 1]
    #    1    1 96 216      Gender
    #    2    2 33  81      Gender
    #    3 <NA>  0   0      Gender
    colnames(cdata)[4] <- "Features"
    freq <- rbindlist(list(freq, cdata))
    #plot(ca(as.integer(cdata[-3,])))
    # DEBUG
  }
  freq <- freq[-1,] #removal of 1st row: NA
  name_freq <- colnames(osccCleanNA_conTCGA[L:R])
  # debug
  name_freq <- t(name_freq) # "Gender" "age.at.diagnosis" "T"  "N"  "M"  "stage_2" and "margin" "tobacco"
  
  # array indexing https://cran.r-project.org/doc/manuals/r-release/R-intro.html#Array-indexing
  results <- list (TRUE, chiT, freq) # it x should to be returned/updated the osccCleanNA_conTCGA
  return(results)
} # end of contingencyTCGA function

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
  freq_features <- as.character(freq$Features[seq(1, nrow(freq), 3)]) # pick one per 3 rows
  #c("Gender","age.at.diagnosis", "T", "N", "M", "stage_2","margin", "tobacco"); nrow(freq)==8*3
  # rows need to be selected and reordered
  # ***match freq and osccCleanNA_conBin of colnames
  colnames(osccCleanNA_conBin)[2:(1+nrow(freq)/3)] <- freq_features # [2:9]
  # ***a feature-name mapping (indTbChi) of " freq$Features  <- osccCleanNA" ; 
  # ***it needs to be updated.
  indTbChi <- data.frame(cbind(featuresUni[-length(featuresUni)],
                               freq_features,
                               colnames(osccCleanNA_conBin)[2:(1+nrow(freq)/3)])) # from Gender to tobacco
  colnames(indTbChi) <- c("featuresUni", "freq$Features", "osccCleanNA_conBin")
  # > indTbChi as a data.frame (mapping table)
  # #               featuresUni-1   freq$Features    osccCleanNA_conBin(osccCleanNA)
  # 1                 Gender           Gender       Gender
  # 2       Age at diagnosis age.at.diagnosis        age.at.diagnosis
  # 3    Clinical T status                T         T
  # 4    Clinical N status                N         N
  # 5    Clinical M status                M         M
  # 6       Clinical Stage          stage_2        stage_2
  # 7 Surgical Margin status           margin       margin
  # 8 Tobacco Exposure       tobacco_exposure tobacco_exposure
  
  # rownames(tableChi1) <- featuresUni
  
  
  # reset tableChi1
  tableChi1 <- as.data.frame(setNames(replicate(length(contLowHighN), numeric(0), simplify = F), contLowHighN)) # declare a xlsx output data.frame (dynamic table)
  # colnames(tableChi1) <- contLowHighN
  for (i in 1:(length(featuresUni)-1)) { # without PMM1 on the last row (i in 1:8)
    # generate table 2/tableChi1 by every two rows
    mapi <- as.character(freq$Features) == as.character(indTbChi$osccCleanNA_conBin[i]) # "Gender"...to "Tobacco" in [47:49] of freq
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





## [mainA process #part A] { ####
# genome-wide scan for margin 0 or 1 cohort ##
# TCGA_HNSCC_marginS.R

# # we might compare the difference of P-value between chisq.test and fisher.test
# # ../xlsx and ../xlsx_chisq
# # run01 done on [2019/06/26]


# run KM plot # 即可：run n=521 cohort, survival analysis without RNAseq, and categorized by tobacco_exposure high/low.
marginTag <- "_marginS_" # with margin 0 or 1

# cutoffReturn <- cutofFinder_func(geneName, marginTag) # with return cutoff1 at # of patient
# if (length(cutoffReturn) == 8) {
#   
#   osccCleanNA <- cutoffReturn[[1]] # taken for freq, chiT as well as tableChi1
#   case50_n <- cutoffReturn[[2]] # of patients to cut
#   cutoff1 <- cutoffReturn[[3]] # cutoff1 at RNAseq level
#   surv_OS1 <- cutoffReturn[[4]]
#   p_OS0 <- p_OS1 <- as.numeric(cutoffReturn[[5]]) # optimized P-Value
#   OS <- cutoffReturn[[6]] # cases_OS and p_OS
#   RFS <- cutoffReturn[[7]] # cases_RFS and p_RFS
#   OS_pvalue <- cutoffReturn[[8]]
#   osccCleanNAM_pos <- which(colnames(osccCleanNA) == paste("PMM1", "_median",        sep=""))
#   } else if (length(cutoffReturn) == 1) {
#   return(cutoffReturn) # return with error code
#   }

# to assembly osccCleanNA data.frame
# $ scp tex@35.201.169.0:~/R/HNSCC.clinical.RNAseq.Fire.Rda ./
load(file="~/R/HNSCC.clinical.RNAseq.Fire.Rda") 
# n=521 [2019/07/27] # as clean6_oscc_tobacco (with tobacco_exposure feature created)
oscc <- clean6_oscc_tobacco[, c(1:10, 12, 14:15)]
colnames(oscc) <- c("Unique.ID","Gender","age.at.diagnosis",
                                       "T","N",
                                       "M","stage_2","margin", "tobacco", 
                                       "OS_IND", "OS..months._from.biopsy",
                                       "RFS_IND", "RFS..months._from.op")
#, "H.score_T")
# starting analysis with "oscc" (without RNAseq)
# 
# [2019/06/11][2019/07/30] with tobacco
# > check complete.cases, data cleaning ####
osccClean <- oscc
# commonFeatures <- c(2:6, 8) # common features: gender, age, margin and TNM :-)
commonFeatures <- which(colnames(oscc) %in% c("Gender", "age.at.diagnosis", "T", "N", "M", "margin", "tobacco")) 
#  7 essential features is located at 2 3 4 5 6 8 9;
osccCleanNA <- osccClean[complete.cases(osccClean[, commonFeatures]), ] # copy all features but remove NaN entities in 7 essential features for their "completeness"
# n=415 cases with feature of tobacco_exposure

# *** addNA for counting all NA (e.g. there is 0, no 1, in margin-free cohort) in "M" "stage_2" "margin" or "tobacco, as a level of factor
# the “pathological” case of two different kinds of NAs which are treated differently: exclude = if (useNA == "no") c(NA, NaN)
# "unusual NA comes from addNA() as factor
# table(osccCleanNA$margin)
#   0   1 (NaN) two level factor 
#328  99   (0)
osccCleanNA$margin <- addNA(osccCleanNA$margin, ifany=F) # ifany=="always" add NA as a factor
#   0    1 <NA> => 3 levels of factor in "margin"
# str(osccCleanNA$margin)
# Factor w/ 3 levels "0","1",NA: 1 2 2 1 1 2 1 2 1 1 ...
osccCleanNA$tobacco <- addNA(osccCleanNA$tobacco, ifany=F) # ifany=="always" add NA as a factor
# str(osccCleanNA$tobacco)

# by marginTag:   # surgical margin status
if (marginTag == "_marginS_") {
  
} else if (marginTag == "_marginFree_") {
  osccCleanNA_freeMargin <- osccCleanNA[osccCleanNA$margin == 0, ] # margin==0
  # margin free cohort (n=351): keeping 0(-) and excluding 1(or +) margin ### 
  osccCleanNA <- osccCleanNA_freeMargin
} else if (marginTag == "_marginPlus_") {
  osccCleanNA_PlusMargin <- osccCleanNA[osccCleanNA$margin == 1, ] # margin==1 (close, indeterminate, or positive margins)
  # margin positive cohort (n=109): keeping  1(or +) margin ### 
  osccCleanNA <- osccCleanNA_PlusMargin
}
#
oscc <- oscc0 <- osccCleanNA # syncrhonize them

# checking the completeness of each row: (i.e. rows without NaN)
which(complete.cases(oscc$OS_IND)==F) # if no NaN -> 0

# *** checking column 10 is OS_IND
OS_IND_pos <- which(colnames(oscc) == "OS_IND")
which(complete.cases(oscc[oscc$OS_IND==1, OS_IND_pos])==F) #OS_IND ==1, death event (dead) => result 0, if no NaN

cohort_n <- nrow(oscc) # n=328 or 427 -> 415

p_OS <-  data.frame(matrix(data = NA, nrow = cohort_n, ncol = 1))
cases_OS <-  data.frame(matrix(data = NA, nrow = cohort_n, ncol = 1))
p_RFS <-  data.frame(matrix(data = NA, nrow = cohort_n, ncol = 1))
cases_RFS <-  data.frame(matrix(data = NA, nrow = cohort_n, ncol = 1))




# oscc, n=415
# run survival analysis: OS (modified from cutoff finder, run100) ####
library(survival)
mysurv <- Surv(oscc$OS..months._from.biopsy, oscc$OS_IND==1) #1==death event

# grouping by PMM1_median -> tobacco_exposure "tobacco"
# Test for difference (log-rank test) with P-value
oscc_Tobacco_pos <- which(colnames(oscc) == "tobacco")
if (is.na(tryCatch(surv_OS <- survdiff(mysurv ~ as.vector(unlist(oscc[oscc_Tobacco_pos]), mode="numeric"), data=oscc), error = function(e) return(NA)))) {cat("one group only")}
# OS, error 2 due to ZSWIM2? (one group only) ###
# extract P-value from surv_OS, put in "original" position (unsorted)

#x *** shift here [run100+1]
#x p_OS[exp_geneName_sorted$ix[run100+1], 1] <- 
  p_OS0 <- format(pchisq(surv_OS$chisq, length(surv_OS$n)-1, lower.tail = FALSE), digits=3)
# P-value = "0.0431"
print(paste("KM P-value", p_OS0, "(<=0.05), which is", (p_OS0<=0.05)))
#[coxph] - fits a Cox proportional hazards regression model



# [2019/08/01] 
# [survfit] - Kaplan-Meier curve, or KM plot of OS ####
OS.km <- survfit(mysurv ~ as.vector(unlist(oscc[oscc_Tobacco_pos]), mode="numeric"), data=oscc, type= "kaplan-meier", conf.type = "log-log")
tiff(file=paste("KMplot_OS_tobacco_exposure.tiff", sep = ""))
plot(OS.km, lty=1, xscale=12, xmax=60, col=c("blue","red"), sub=paste("P-Value =", p_OS0), main=paste("OS in TCGA", TCGA_cohort, "(tobacco exposure, n=", surv_OS$n[1]+surv_OS$n[2],")"), ylab="Percent Survival", xlab="Years")
legend("topright", legend=c(paste("low risk (",surv_OS$n[1], ")"), paste("high risk (",surv_OS$n[2], ")")), lty=1:1, col=c("blue","red"))
dev.off()
#  
# [survfit] - Kaplan-Meier curve, or KM plot of OS
# confidence intervals as log hazard or log(-log(survival))
#OS.km <- survfit(mysurv ~ as.vector(unlist(osccCleanNA[, osccCleanNAM_pos]), mode="numeric"), data=osccCleanNA, type= "kaplan-meier", conf.type = "log-log")
# 365.25 days in a year for xscale => 3650 days for 10 years
# 12 months per year for 5 years => 60 months



# 
# table 3, table 4 ####
# hazard ratios of features (including tobacco_exposure, without RNAseq) (high/low)
# By Cox HR [COXPH modelling]
library(rJava) # PASSED
library("xlsx")
library("openssl")
library(r2excel) 
#*** [Multivariate for OS] Table 3 right panle

# selectedFeatures <- colnames(osccCleanNA[commonFeatures])
# colname of oscc = c("Unique.ID","Gender","age.at.diagnosis",
# "margin","T","N",
# "M","stage_2","OS..months._from.biopsy",
# "OS_IND","RFS..months._from.op","RFS_IND","H.score_T")

# ***rename all of colnames (features) as osccT [old coding]
colnames(osccCleanNA) <- coln_osccT # colnames (features) of osccT

osHRmulti <- 0
# 9 + OS/RFS features in HNSCC
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
                 tobacco +
                 #  R.T +
                 #  C.T +
                 #                 presence_of_pathological_nodal_extracapsular_spread +
                 #                 neoplasm_histologic_grade +
                 as.vector(osccCleanNA[, osccCleanNAM_pos], mode="numeric"), 
               data=osccCleanNA) # PMM1_median == IHC score == (H.score_T) in [low, high]
# warning() X matrix deemed to be singular; variable 8:margin
# summary(oscox)

osHRmulti <- cbind(data.frame(summary(oscox)$conf.int)[-2], data.frame(summary(oscox)$coefficients)[5])

# *** removal the "noise" for each new feature specified
skipNA <- rownames(osHRmulti) %in% c(1, "marginNA", "tobaccoNA")
osHRmulti <- osHRmulti[which(!skipNA), ]


# correction of rownames
rownames(osHRmulti) <- featuresUni
#rownames(osHRmulti)[nrow(osHRmulti)] <- paste("PMM1", "_median", sep="")
osHRmulti <- round(osHRmulti, 3) # rounding them by 3 decimal

#View(osHRmulti) # table 3 right panel



#*** [Univariate for OS]  table 3 left panel
# # 9 features in HNSCC, put in coxph one by one
# number of features 13 -> 14 -> 15 (add a feature: margin status, tobacco exposure)
oscox <- 0
## boundary of features column from first to last # 
LL <- which(colnames(osccCleanNA) == "Gender") # "left" feature
RR <- which(colnames(osccCleanNA) == "tobacco") # new "right" feature
pmm1_pos <- which(colnames(osccCleanNA) == paste("PMM1", "_median", sep=""))

features_os <- colnames(osccCleanNA)[c(LL:RR, pmm1_pos)] # features selection ["Gender" to "tobacco", "PMM1"] 9 out of 15
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
## DONE and then export part II ##

## R2Excel export ####
## and save at "xlsx/"
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
filenamex <- paste("xlsx/", TCGA_cohort, "_survivalAnalysis_marginS_", geneName, ".xlsx", sep = "")
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
# Part I:Table 2. The clinicopathological features
#  (tableChi1) from Calling contingencyTCGA and contignecyBin ####
# x using tableChi2 by contigencyBin2 with c("Gender","ageDx", "pathologic_T", "pathologic_N", "pathologic_M", "stage","margin" )
# under best cutoff value (auto choosen which has the smallest P-value)

contiT <- contingencyTCGA(osccCleanNA, geneName) # calling this function (OSCC cohort, PMM1) 2 parameters, no more cutoff1; 
# "margin" at column 8 of oscc
# "tobacco" at column 9 of oscc
#  oscc <- contiT[[1]] # updating PMM1_median
chiT <- contiT[[2]] # extrac it from list by using [[]]; chiT$X2 is the P-value
freq <- contiT[[3]] # well DONE

tableChi1 <- contingencyBin (osccCleanNA, chiT, freq) # calculating the P-value
# ***add rownames for tableChi1
nrow_featuresUni <- length(featuresUni) # aka. 9
nrow_tableChi1 <- as.character(seq(1, 2*(nrow_featuresUni-1))) # aka. 2 *8 個NA (a vector)
nrow_tableChi1[c(T,F)] <- featuresUni[-nrow_featuresUni] # skip last one row RNAseq(z-score)
# > nrow_tableChi1
# [1] "Gender"                 "2"                     
# [3] "Age at diagnosis"       "4"                     
# [5] "Clinical T Status"      "6"                     
# [7] "Clinical N Status"      "8"                     
# [9] "Clinical M Status"      "10"                    
# [11] "Clinical Stage"         "12"                    
# [13] "Surgical Margin status" "14"                    
# [15] "Tobacco Exposure"       "16"   


rownames(tableChi1) <- nrow_tableChi1 # duplicate 'row.names' are not allowed
# 16 rows
# "margin" at column 8 of osccCleanNA

# to save "*" of column "remark" in 300 tableChi1 ###
#?? cut_featuresRemark <- cbind (cut_featuresRemark, tableChi1[7], tableChi1$Remark) # retrieve the P-value and remark
#

# header
xlsx.addHeader(wb, sheet, value=paste("Table 2. The clinicopathological features of ", TCGA_cohort, " cohort and ", geneName,  " expression. (Chi square test)", sep=""),
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

# for xlsx output (*** dynamic tableOS) following new feature
# setNames, This is a convenience function that sets the names on an object and returns the object.
# a.k.a. create a empty data.frame: tableOS1; however, "()" will be assigned as "."
tableOS1 <- as.data.frame(setNames(replicate(length(ciUniMti), numeric(0), simplify = F), ciUniMti)) 
colnames(tableOS1) <- ciUniMti
for (i in 1:nrow(tableOS)) { # eg. 1~9
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
xlsx.addHeader(wb, sheet, value=paste("Table 3. Univariate/Multivariate Cox's proportional hazards regression analyses on OS time of", geneName, "gene expression in ", TCGA_cohort), level=5)
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
xlsx.addHeader(wb, sheet, value=paste("Table 4. Univariate/Multivariate Cox's proportional hazards regression analyses on RFS time of", geneName, "gene expression in ", TCGA_cohort), level=5)
xlsx.addLineBreak(sheet, 1)
xlsx.addTable(wb, sheet, data = t(data.frame(c("Univariate", "\t", "\t", "\t", "Multivariate"))), fontSize=12, startCol=4,
              row.names = F, col.names = F)
#
xlsx.addTable(wb, sheet, data= tableRFS1, startCol=2,
              fontColor="darkblue", fontSize=12,
              rowFill=c("white", "lightblue"), row.names = TRUE
)
xlsx.addLineBreak(sheet, 25) # more space



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
save(list = c("tableChi1", "tableOS1", "tableRFS1", "OS_pvalue", "RFS_pvalue"), file=paste("xlsx/HNSCC_survivalAnalysis_marginS_", geneName, ".Rda", sep=""))
#, error = function(e) return(NA))
#print(paste("Create", paste("LUAD_survivalAnalysis_marginS_", geneName, ".Rda", sep=""), "successfully."))


## 
# Final Return ####
if (nrow(OS_pvalue) > 0)  { # we hit one gene with P-value < 0.05 in KM plot
  return(which.min(OS_pvalue$p_OS))} else {return(0)}

#  ) #system.time end
} # function END (survival_marginS)
##-- -- -- -- --


