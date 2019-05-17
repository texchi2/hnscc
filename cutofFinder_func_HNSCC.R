# # [2019/05/16] Start recoding/debugging, if any; since [2018/10/?], actually
# # on instance-4 for HNSCC
# # Loading .Rda ###
# # Resume:[5a.START]
# geneName <- "TXNDC11"
# # Start the survival analysis for each individual gene
# 
# ## # set path on google drive
# #library(FirebrowseR)
# TCGA_cohort <- "LUAD" # cancer type
# path_LUAD <- "~/R/LUAD_Peter_survival/" # under rstudio-server on GCP
# #path_LUAD <- "~/R/LUAD_Peter_survival/mount" # mount google drive to ubuntu@instances-4 at GCP
# #path_LUAD <- "/Users/apple/Google\ Drive/2018_PhD_dissertation/LUAD_Peter_survival/"
# # path_LUAD <- "/Users/apple/Documents/My\ Tableau\ Repository/Workbooks/"
# setwd(path_LUAD) # set the working directory to the google drive => GCP
# 
# load(file="whole_genome.Rda") # the name list of protein coding genome
# # below 2 lines: it needs to be modified for a function cal
# LUAD_n <- length(whole_genome) #last one 20499: "ZZZ3"
# 
# # source(paste(path_LUAD, "TCGA_LUAD_marginS.R", sep="")) # survival_marginS <- function() {} in TCGA_LUAD_marginS.R
# # source(paste(path_LUAD, "TCGA_LUAD_marginFree.R", sep="")) # survival_marginFree <- function() {} in TCGA_LUAD_marginFree.R
# 
# 
# #install.packages("gmailr") # dependencies ‘curl’, ‘openssl’ are not available for package ‘httr’
# #  'libcurl'
# # make from the source "curl" and its libcurl 
# # $ wget https://github.com/curl/curl/releases/download/curl-7_59_0/curl-7.59.0.tar.gz
# # $ tar -xzvf curl-7.59.0.tar.gz 
# # $ cd curl-7.59.0/
# # $ ./configure # make  # sudo make install
# # then install curl (done)
# # > install.packages("curl", repos="https://cran.r-project.org", type="source") 
# 
# # "openssl"(done), "httr"(done), "git2r"(done)
# # install.packages("git2r", repos="https://cran.r-project.org", type="source")
# #install.packages("devtools")
# #devtools::install_github("hoxo-m/pforeach")
# library(gmailr) # notify me the script progress by email...https://developers.google.com/gmail/api/
# library(scales)
# # library(doParallel) # core n=2 in my CPU
# # library(foreach)
# # #library(pforeach)
# library(plyr); library(dplyr) #ddply()
# library(magrittr)
# recipient <- "texchi2@gmail.com"
# sender <- "texchi2@gmail.com"
# 
# 
# ### $setup global variables...(including "margin") ###
# library("psych") # for describe()
# library(survival)
# 
# # ***colnames of osccT
# coln_osccT <- c("Unique.ID","Gender","ageDx",
#                 "pathologic_T","pathologic_N",
#                 "pathologic_M","stage", "margin",
#                 "OS_IND","OS..months._from.biopsy",
#                 "RFS_IND", "RFS..months._from.op", "H.score_T", paste("PMM1", "_median", sep=""))
# 
# 
# # for tableChi1 (table 2)
# contLowHighN <- c("Features",	"Low", "(%)",	"High",	"(%)", "Case no",	"P-value", "Remark", "Matthews")
# 
# # Create the Sheet title and subtitle; underline as 0 (default value, no underline), 1 (underline with one line), 2 (underline with two lines)
# # Table 3. Univariate/Multivariate Cox proportional hazards regression analyses on OS time
# ciUniMti <- c("Features",	"HR",	"CI95%(L)",	"CI95%(H)",	"P-value",	"HR",	"CI95%(L)",	"CI95%(H)",	"P-value")
# 
# # with 8 features: in table 3 and table 4
# featuresUni <- c("Gender",
#                  "Age at diagnosis",
#                  #                 "Primary site",
#                  #                 "Clinical T Status",
#                  #                 "Clinical N Status",
#                  #                 "Clinical Stage",
#                  "Pathologic T status",
#                  "Pathologic N status",
#                  "Pathologic M status",
#                  "Pathologic Stage",
#                  "Surgical Margin status",
#                  #                 "Lymphovascular Invasion",
#                  #                 "Perineural Invasion",
#                  #                 "Extranodal spreading of neck LN",
#                  #                 "SCC Histologic Grade",
#                  #                 "Radiotherapy",
#                  #                 "Chemotherapy",
#                  paste("RNAseq(z-score)") # "IHC score" == "z-score"
# )
# # z-score calculation = [(value gene X in tumor Y)-(mean gene X in normal)]/(standard deviation X in normal)
# 
# # lower rowname of table 3 and table 4
# colUni_1 <- c("Male",
#               ">65y",
#               #              "higher risk",
#               #              "T3+T4",
#               #              "N1-3",
#               #              "Stage III+IV",
#               "T3+T4",
#               "N1-3",
#               "M1",
#               "Stage III+IV",
#               "Positive", # positive safety margin
#               #              "Yes",# RT
#               #              "Yes",# CT
#               #              "G3+G4",
#               "High") # PMM1
# # upper rowname of table 3 and table 4
# colUni_0 <- c("Female",
#               "<=65y",
#               #              "lower risk",
#               #              "T1+T2",
#               #              "N0",
#               #              "Stage I+II",
#               "T1+T2",
#               "N0",
#               "M0",
#               "Stage I+II",
#               "Negative", # negative safety margin
#               #              "no", # RT
#               #              "no", # CT
#               #              "G1+G2",
#               "Low") # PMM1
# 

# # Define global functions# ###
# # contingency function, Tex's design (a confusion matrix with TP TN FP FN)
# 
# 
# 
# # # #x for DEBUG only {
# # osccCleanNA <- oscc # T,N,M and stage_2
# # # osccCleanNA has n=245 in LUAD (with margin free) or n=256 in LUAD (margin 0 or 1)
# # # geneName
# # #?? cutoff1 <- i <- round(nrow(oscc)/2) #213
# # osccCleanNA_pos <- which(colnames(osccCleanNA) == "H.score_T")
# # exp_geneName <- t(osccCleanNA[, osccCleanNA_pos ]) # z-score of ZZZ3
# # cutoff1 <- quantile(exp_geneName)[3] # at 50%
# 
# # # contingency P-value
# # contiT <- contingencyTCGA(osccCleanNA, geneName, cutoff1) # calling this function (OSCC cohort, PMM1, cutoff);
# # # "margin" at column 8 of oscc
# # #  oscc <- contiT[[1]] # updating PMM1_median
# # chiT <- contiT[[2]] # extrac it from list by using [[]]; chiT$X2 is the P-value
# # freq <- contiT[[3]] # well DONE
# # }
# 
# 
# 
# 
# ## function declare
# # create correlation table 1 with P-value by chisq.test
# library(reshape)
# library(data.table)
# #library(ca) # for Simple correspondence analysis
# contingencyTCGA <- function(osccCleanNA_conTCGA, geneName) { # no more "run100"; do not need cutoff1 here
#   
#   #!!!!!# L <- 2 ("Gender"); R <- 8 ("margin") # in LUAD
#   ## boundary of features column from first to last # 
#   L <- which(colnames(osccCleanNA_conTCGA) == "Gender")
#   R <- which(colnames(osccCleanNA_conTCGA) == "margin")
#   chiT <- data.frame(matrix(data = NA, nrow = R, ncol = 2)) # create a empty data.frame, 2D matrix as R*2
#   #rown = 8
#   freq <- data.frame(matrix(data = NA, nrow = 1, ncol = 4)) #, dimnames = list(c(1:rown+1), c("Var2", "L", "H"))))
#   colnames(freq) <- c("Var2", 0, 1, "Features") # column 4 is for mapi
#   #x freq <- array(NA, dim = c(5, 3, R)) # row, col, and R as 3D array
#   
#   #!!! PMM1 score## at col 13, 14
#   #osccCleanNA_conTCGA_pos <- which(colnames(osccCleanNA_conTCGA) == "H.score_T") # position of PMM1 IHC score
#   osccCleanNA_conTCGAM_pos <- which(colnames(osccCleanNA_conTCGA) == paste("PMM1", "_median", sep=""))
#   #  exp_geneName <- t(osccCleanNA[, osccCleanNA_pos ])
#   
#   # # it should be done at marginS.R
#   # # then [binomial of gene_median] resume the correlation tables
#   # osccCleanNA[osccCleanNAM_pos] <-(osccCleanNA[, osccCleanNA_pos ] >= cutoff1) +0 # binomial after osccCleanNA 
#   # # ***** addNA for counting all NA (e.g. there is 0, no 1, in margin-free cohort) in "M" "stage_2" or "margin"
#   # # the “pathological” case of two different kinds of NAs which are treated differently: exclude = if (useNA == "no") c(NA, NaN)
#   # # "unusual NA comes from addNA() as factor
#   # # [2018/03/14] finallized debug
#   # osccCleanNA$margin <- addNA(osccCleanNA$margin, ifany=F) # ifany=="always" add NA as a factor
#   # #   0    1 <NA> => 3 levels of factor in "margin"
#   # # 245   11    0 
#   # #  osccCleanNA$M <- addNA(osccCleanNA$M, ifany=F) # "always" add NA as a factor
#   # #  osccCleanNA$stage_2 <- addNA(osccCleanNA$stage_2, ifany=F) # "always" add NA as a factor
#   # 
#   
#   # chisq.test(matrix(c(22, 5, 38, 21), ncol = 2), correct=F)$p.value # a example
#   for (ii in L:R){ # generate freq table
#     # Chi-square, Contingent table correlation, binary variables
#     
#     # build a contingency table https://www.rdocumentation.org/packages/base/versions/3.4.3/topics/table
#     # Powered by DataCamp, it might be running online
#     # t, table from col of PMM1_median vs col ii (L "Gender" to R "margin"); deparse the argument (colnames) by deparse.level 2
#     # # {DEBUG 
#     # ii<-L-1
#     #     ii<-ii+1
#     # # DEBUG}
#     t<- NULL
#     t <- table(osccCleanNA_conTCGA[,osccCleanNA_conTCGAM_pos], osccCleanNA_conTCGA[,ii], useNA = "ifany") #dnn=c(colnames(osccCleanNA_conTCGA[osccCleanNA_conTCGAM_pos])))  
#     # useNA = "ifany"; "always" is for "margin" (with n=0 count)
#     chiT[ii,1] <- colnames(osccCleanNA_conTCGA[ii]) # name list of feature variables from L to R
#     check_p <- chisq.test(t)$p.value # retrieved in chiT$X2
#     if (is.na(check_p)==T) {check_p <- chisq.test(t[1:2,1:2])$p.value}
#     chiT[ii,2] <- check_p
#     # chisq is sum( (o-e)^2/e ); the Na should be removed, You have zero frequencies in 2 counts.
#     #warnings(): In chisq.test(t) : Chi-squared approximation may be incorrect
#     #=> fisher.test(a) # Fisher exact test is better in small counts. 
#     #table(t); #print(c("ii=",paste(ii)))
#     # t
#     # 0   5   6 122 123 
#     # 2   1   1   1   1 
#     #print(ii) # debug
#     
#     obs <- as.data.frame(chisq.test(t)$observed)
#     
#     
#     # >table(t)
#     # 50 57 62 76 
#     # 1  1  1  1 
#     # [1] "-0.0765952029885 \n c(127, 129)"
#     #   Var1 Var2 Freq
#     # 1    0    0  127
#     # 2    1    0   129
#     # 3    0    1    0
#     # 4    1    1    0
#     # debug
#     # with error # mdata <- melt(obs, id=c("Var2", "Var1")) # using the colnames of table t
#     # id variables not found in data:" Var2, Var1 [2018/03/27] , only error on "Gene name = ZSWIM2"(20482)
#     
#     #print(paste("Run", run100, geneName, "(", which(whole_genome==geneName)," ): obs", obs, t))
#     # under cutoff finding 100
#     
#     # ***ZSWIM2 skip (?) 
#     # Error when "melt": id variables not found in data: Var2, Var1; (360: mdata <- melt(obs, id=c("Var2", "Var1")) # using the colnames of table t)
#     #if (is.na(tryCatch([expression], error = function(e) return(NA)))) {return(list("skip", chiT, freq))}
#     #if (obs$Freq==c(127,129,"NA","NA")) {print(paste("skip", geneName)); return(list("skip", chiT, freq))} # back to marginS.R
#     # warnings():  the condition has length > 1 and only the first element will be used
#     # [3] "0.238555836961 \n c(170, 75, 9, 2, 0, 0)"
#     #   Var1 Var2 Freq
#     # 1    0    0  170
#     # 2    1    0   75
#     # 3    0    1    9
#     # 4    1    1    2
#     # 5    0 <NA>    0
#     # 6    1 <NA>    0
#     # debug
#     if (is.na(tryCatch(mdata <- melt(obs, id.vars=c("Var2", "Var1")), error = function(e) return(NA)))) {return(list("skip", chiT, freq))} # back to marginS.R
#     # mdata <- melt(obs, id.vars=c("Var2", "Var1")) # error when there is only one group
#     # debug for "ZSWIM2"(20482)
#     
#     cdata <- as.data.frame(cast(mdata, Var2 ~ Var1)) 
#     cdata <- cdata[, c(1:3)] # cdata should has 3 columns
#     # names(dimnames(cdata)) = c("expression", "severity")
#     # NA is necessary in column Var2
#     if (!anyNA(cdata$Var2)) {
#       cdata <- rbind(cdata, c(NA, 0, 0))
#     }
#     if (length(cdata$Var2) < 3) {
#       cdata$ind <- seq_len(nrow(cdata)) # index as 1 2 3...
#       crow <- data.frame("Var2"=2, "0"=0, "1"=0, "ind"=1.5)
#       colnames(crow) <- colnames(cdata)
#       cdata <- rbind(cdata, crow)
#       cdata <- cdata[order(cdata$ind),]
#       cdata <- cdata[,-4] # removal of index
#     }
#     cdata <- cbind(cdata, chiT[ii,1]) # add features name => # cdata has 4 columns
#     # print(cdata)
#     #freq[1+(ii-1)*rown,] <- cdata
#     freq <- rbindlist(list(freq, cdata))
#     #plot(ca(as.integer(cdata[-3,])))
#     # DEBUG
#   }
#   freq <- freq[-1,] #removal of 1st row: NA
#   name_freq <- colnames(osccCleanNA_conTCGA[L:R])
#   name_freq <- t(name_freq) # "Gender" "age.at.diagnosis" "T"  "N"  "M"  "stage_2" and "margin"
#   
#   # array indexing https://cran.r-project.org/doc/manuals/r-release/R-intro.html#Array-indexing
#   results <- list (TRUE, chiT, freq) # it x should to be returned/updated the osccCleanNA_conTCGA
#   return(results)
# } # end of contingencyTCGA function
# #
# 
# 
# 
# 
# # Define function introdcue Matthews Correlation Coefficient,
# # \cite{Sanz-Pamplona2012} The Matthews Correlation Coefficient (MCC) [57] was
# # chosen as measure of classification accuracy [58]. This index combines test
# # sensitivity and specificity. It ranges from ???1 to 1 and its interpretation is
# # similar to the Pearson???s correlation coefficient. In the context of a
# # classification problem it is expected that MCC ranges 
# # from 0 (no prediction ability at all) to +1 (perfect prediction) with negative values near zero
# # possibly occurring in random classifiers due to sample variability. MCC values
# # higher than 0.3 can be considered as indicative of high predictive value as they
# # correspond to more than 65% accuracy in balanced data.
# 
# 
# 
# 
# ##
# library(R.utils) # intToBin()
# library(compositions) # unbinary()
# contingencyBin <- function (osccCleanNA_conBin, chiT, freq) {
#   # to generate tableChi1 (table 2)
#   # enough, geneName is not necessary a parameter (osccCleanNA_conBin <- osccCleanNA)
#   # processing chiT and processing freq; place a "remark" and Matthews
#   # for binary and integer convertion
#   # library(R.utils) # intToBin()
#   # library(compositions) # unbinary()
#   # sigContig <- c(NA, "*") # remarkable when (p<0.05 || scoreContig == c(5,10))
#   # 
#   freq_features <- as.character(freq$Features[seq(1, nrow(freq), 3)]) # every 3 rows, pick one
#   #c("Gender","age.at.diagnosis", "T", "N", "M", "stage_2","margin")
#   # rows need to be selected and reordered
#   # ***match freq and osccCleanNA_conBin of colnames
#   colnames(osccCleanNA_conBin)[2:8] <- freq_features
#   # ***a feature-name mapping (indTbChi) of " freq$Features  <- osccCleanNA" ; it needs to be updated.
#   indTbChi <- data.frame(cbind(featuresUni[-length(featuresUni)],
#                                freq_features,
#                                colnames(osccCleanNA_conBin)[2:8])) # from Gender to margin
#   colnames(indTbChi) <- c("featuresUni", "freq$Features", "osccCleanNA_conBin")
#   # > indTbChi as a data.frame (mapping table)
#   # #               featuresUni-1   freq$Features    osccCleanNA_conBin(osccCleanNA)
#   # 1                 Gender           Gender       Gender
#   # 2       Age at diagnosis age.at.diagnosis        age.at.diagnosis
#   # 3    Pathologic T status                T         T
#   # 4    Pathologic N status                N         N
#   # 5    Pathologic M status                M         M
#   # 6       Pathologic Stage          stage_2        stage_2
#   # 7 Surgical Margin status           margin       margin
#   # #
#   # rownames(tableChi1) <- featuresUni
#   
#   
#   # reset tableChi1
#   tableChi1 <- as.data.frame(setNames(replicate(length(contLowHighN), numeric(0), simplify = F), contLowHighN)) # declare a xlsx output data.frame (dynamic table)
#   # colnames(tableChi1) <- contLowHighN
#   for (i in 1:(length(featuresUni)-1)) { # without PMM1 on the last row (i in 1:7)
#     # generate table 2/tableChi1 by every two rows
#     mapi <- as.character(freq$Features) == as.character(indTbChi$osccCleanNA_conBin[i]) # "Gender"... in [47:49] of freq
#     mapi_pos <- which(mapi == T) #[47:49], every 3
#     var2U <- as.data.frame(freq[mapi_pos[1], ])
#     var2L <- as.data.frame(freq[mapi_pos[2], ])
#     subtotal <- sum(var2U[2:3]) + sum(var2L[2:3])
#     U2 <- var2U[2];  U3 <- var2U[3]
#     L2 <- var2L[2];  L3 <- var2L[3]
#     pContig <- round(chiT[which(chiT$X1 == as.character(indTbChi$osccCleanNA_conBin[i])), 2], 4) # P-value retrieving from chiT$X2
#     
#     abin <- as.numeric(c(U3>U2, L3>L2, L2>U2, L3>U3)) # pattern of Sn=L3/(L2+L3), Sp=U2/(U2+U3)
#     #    abin <- c(1,0,1,0) ; binary 1010 == decimal 10; while binary 0101 == decimal 5.
#     achar <- paste(abin, collapse="") # ex. "0001"
#     #  # a integer score, 10 or 5 and more is significant, which will be marked as "*" in remark
#     scoreContig <- unbinary(achar) # 'structure(list(), *) warning() :-)
#     # Matthews correlation coefficient: 0-1 perfect prediction
#     #[OFF it]:    matthews <- (U2*L3-U3*L2)/sqrt((L3+U3)*(L3+L2)*(U2+U3)*(U2+L2))
#     matthews <- 0
#     
#     # colUni_0 and colUni_1: value of feature => "male"/"female"....
#     row0T2 <- data.frame(t(c(colUni_0[i], as.numeric(c(U2, round(U2/subtotal*100, 1), U3, round(U3/subtotal*100, 1), subtotal, pContig, (as.numeric((!(pContig>0.05) || scoreContig == c(5,10,1,2,4,7,8,11,13,14)))+0), matthews ))))) # remark 0 or 1
#     row1T2 <- data.frame(t(c(colUni_1[i], as.numeric(c(L2, round(L2/subtotal*100, 1), L3, round(L3/subtotal*100, 1), NA, NA, NA, NA)))))
#     #  colnames(row0) <- contLowHighN
#     #  row0[,2:ncol(tableChi1)] <- as.numeric(as.character(row0[,2:ncol(tableChi1)]))
#     tableChi1 <- rbind(tableChi1, row0T2, row1T2)
#   }
#   # tableChi1 <- cbind(tableChi1, chiT[-(1:5), 2]) # can not direct paste the p-value vector
#   colnames(tableChi1) <- contLowHighN
#   
#   return(tableChi1)
# } # end of contingencyBin define ++++++
# 
# 
# 
# 
cutofFinder_func <- function(geneName, marginTag) {
# marginTag == "_marginS_" or == "_marginFree_"
# geneName
# to assembly osccCleanNA data.frame
#
#load(file="TMU_TA51BCDE_T70_clinical_fullFeatures13_dichotomized.Rda") # as oscc (without margin feature)
# or ***
load(file="LUAD_T522_clinical_fullFeatures11_dichotomized.Rda") # as oscc with "margin" features
#, negative n=348 while positive n=18 and 156:NaN)

# load and prepare $$.clinico_mRNA.Fire for survival analysis
#
if (is.na(tryCatch(load(file=paste("LUAD.mRNA.Exp.", geneName, ".Fire.Rda", sep="")), error = function(e) return(NA)))) {return(3)} # there is NO such file
# as LUAD.mRNA.Exp.Fire
# there is NO such file: "LUAD.mRNA.Exp.ZFP91.CNTF.Fire.Rda" %in% dir() => False
#
# "sample_type"          # 59 NT(Solid Tissue Normal) 515 TP(Primary Solid Tumor) or 2 TR(Recurrent Solid Tumor)

LUAD_T_mRNA.Fire <- LUAD.mRNA.Exp.Fire[LUAD.mRNA.Exp.Fire$sample_type %in% c("TP"), c("tcga_participant_barcode","z.score")] # n=517 -2 = 515
# removing duplicated ID %in%  c("TCGA-50-5066","TCGA-50-5946") => TR recurrent sample
# if NULL or more than 30% NA in RNAseq (z.score) => dischard this gene
if (length(LUAD_T_mRNA.Fire) == 0) {return(3)}
if (as.data.frame(colMeans(is.na(LUAD_T_mRNA.Fire)) > 0.3)[2,1] == T) {return(3)}

# inner join by merge
if (is.na(tryCatch(LUAD.clinico_mRNA.Fire <- merge(oscc, LUAD_T_mRNA.Fire, by="tcga_participant_barcode") , error = function(e) return(NA)))) {return(5)} # merge error
#n=515 with sample_type "TP", excluding "NT normal tissue" or "TR"


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
#
osccClean <- oscc
# commonFeatures <- c(4:9) # common features: gender, age, margin and TNM :-)
commonFeatures <- which(colnames(oscc) %in% c("Gender", "age.at.diagnosis", "T", "N", "M", "margin")) #  6 essential features
osccCleanNA <- osccClean[complete.cases(osccClean[, commonFeatures]), ] # copy all features but remove Na or NaN entities in 6 essential features for their "completeness"

# *** addNA for counting all NA (e.g. there is 0, no 1, in margin-free cohort) in "M" "stage_2" or "margin"
# the “pathological” case of two different kinds of NAs which are treated differently: exclude = if (useNA == "no") c(NA, NaN)
# "unusual NA comes from addNA() as factor
osccCleanNA$margin <- addNA(osccCleanNA$margin, ifany=F) # ifany=="always" add NA as a factor
#   0    1 <NA> => 3 levels of factor in "margin"
# 245   11    0


oscc_n256 <- osccCleanNA # n=256; removal of NA cases
# ***_marginS_: margin positive and negative cohort (n=256) ####

osccCleanNA <- oscc_n256

# by marginTag
if (marginTag == "_marginS_") {
  
  } else if (marginTag == "_marginFree_") {
  # surgical margin status: keeping 0 and excluding + margin (as 1; n=11)
  osccCleanNA_freeMargin <- osccCleanNA[osccCleanNA$margin == 0, ] # margin==0
  # margin free cohort (n=245): ### 
  osccCleanNA <- osccCleanNA_freeMargin
  }
#
## 5b. (repeat100) Cutoff finder [osccCleanNA] ###
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
    # => find cutoff1, generate OS P-Value plot and KM plot.
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
    # (for run100) P-value according to KM survival analysis (alone) ####
    for (run100 in seq(cutoff_n[1], cutoff_n[2])){ 
      
      # sorted (by RNAseq) no.78~179 in cases 256 of LUAD
      # use oscc0 for "positioning"; 
      # oscc is the dataset to be analysed here.
    
      # Binominal H.Score_T (RNAseq) by exp_geneName_sorted$x[i]) #RNAseq cutoff
      # exp_geneName_sorted$x[i], cutoff100 is the cutoff value of RNAseq in rank i
      cutoff100 <- exp_geneName_sorted$x[run100+1] #*** shift here
      oscc[osccM_pos] <- (oscc0[oscc0_pos] >= cutoff100) +0 # oscc$PMM1_median <- (oscc0$PMM1 >= i) +0
    #  R4> table(oscc$PMM1_median) # when run100 == 78
    #  0   1 
    #  78 178
    #  *** correction of cases_OS ####
      cases_OS[exp_geneName_sorted$ix[run100+1], 1] <- run100 # group 0 vs 1: 0 has 78 cases ...etc
      cases_RFS[exp_geneName_sorted$ix[run100+1], 1] <- run100 
      #run100 <- i - as.numeric(cutoff_n[1]) +1 # *** incremental (start from 1) by i
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
      # OS, error 2 due to ZSWIM2? (one group only) ###
      # extract P-value from surv_OS, put in "original" position (unsorted)
      # #*** shift here [run100+1]
      p_OS[exp_geneName_sorted$ix[run100+1], 1] <- p_OS0 <- format(pchisq(surv_OS$chisq, length(surv_OS$n)-1, lower.tail = FALSE), digits=3)
    #  *** cases_OS[j] <- surv_OS$n[1] #cutoffs by cases, remaping sorted ###
      #[coxph] - fits a Cox proportional hazards regression model
    
      
      
     #print(paste("run=", run100, "; KM P-value", p_OS0, "(<=0.05), which is", (p_OS0<=0.05))) 
      # [survfit] - Kaplan-Meier curve, or KM plot of OS
       OS.km <- survfit(mysurv ~ as.vector(unlist(oscc[osccM_pos]), mode="numeric"), data=oscc, type= "kaplan-meier", conf.type = "log-log")
    #   jpeg(file=paste("KMplot_OS_", geneName, i, ".jpg", sep = ""))
       plot(OS.km, lty=1, xscale=12, xmax=60, col=c("blue","red"), sub=paste("P-Value =", p_OS0), main=paste("searching OS in TCGA", TCGA_cohort, "(n=", surv_OS$n[1]+surv_OS$n[2],")/", geneName, "cutoff=", cutoff100), ylab="Percent Survival", xlab="Years")
       legend("topright", legend=c(paste("low(",surv_OS$n[1], ")"), paste("high(",surv_OS$n[2], ")")), lty=1:1, col=c("blue","red"))
       # dev.off()
      #  
    # [survfit] - Kaplan-Meier curve, or KM plot of OS
      # confidence intervals as log hazard or log(-log(survival))
      #OS.km <- survfit(mysurv ~ as.vector(unlist(osccCleanNA[, osccCleanNAM_pos]), mode="numeric"), data=osccCleanNA, type= "kaplan-meier", conf.type = "log-log")
      # 365.25 days in a year for xscale => 3650 days for 10 years
      # 12 months per year for 5 years => 60 months
      
      # plot(OS.km, lty=1, xscale=12, xmax=60, col=c("blue","red"), sub=paste("P-value =", p_OS1), main=paste("OS in TCGA", TCGA_cohort, "(n=", surv_OS1$n[1]+surv_OS1$n[2],")/", geneName), ylab="Percent Survival", xlab="Years")
      # legend("topright", legend=c(paste("low(",surv_OS1$n[1], ")"), paste("high(",surv_OS1$n[2], ")")), lty=1:2, col=c("blue","red"))
      #plot(OS.km, lty=1, col=c("blue","red"), sub="p= 0.816", main="OS in osccCleanNA(n=505)/gene level", ylab="Percent Survival", xlab="Days")
      #legend("topright", legend=c('Low', 'High'), lty=1:2, col=c("blue","red"))
      # summary(OS.km, times = seq(0, 3000, 100))
      
      
      # run survival analysis: RFS
      mysurv <- Surv(oscc$RFS..months._from.op, oscc$RFS_IND==1) #1==tumour recurrency
      # Test for difference (log-rank test) with P-value
      # RFS, error 2 = function(e) return(NA)) for "Survdiff.fit error on there is only one group"   # due to ZSWIM2 (one group only) ###
      if (is.na(tryCatch(surv_RFS <- survdiff(mysurv ~ as.vector(unlist(oscc[osccM_pos]), mode="numeric"), data=oscc), error = function(e) return(NA)))) {return(2)}
    
      p_RFS[exp_geneName_sorted$ix[run100], 1] <- format(pchisq(surv_RFS$chisq, length(surv_RFS$n)-1, lower.tail = FALSE), digits=3)
    #  cases_RFS[j] <- surv_RFS$n[1] #cutoffs
      #RFS.km <- survfit(mysurv ~ oscc[osccM_pos], data=oscc, conf.type = "log-log")
      #  jpeg(file=paste("KMplot_RFS_", geneName, i, ".jpg", sep = ""))
      #  plot(RFS.km, lty=1, col=c("blue","red"), sub=paste("p-value =", p_RFS[j]), main=paste("RFS in OSCC(n=", surv_RFS$n[1]+surv_RFS$n[2],")/", geneName, "cutoff=", i), ylab="Percent Survival", xlab="Days")
      #  legend("topright", legend=c(paste("low(",surv_RFS$n[1], ")"), paste("high(",surv_RFS$n[2], ")")), lty=1:2, col=c("blue","red"))
      #  dev.off()
      #  
      #  
    
      
      # library(ggplot2)
      # ggplot(OS, aes(x=rank, y=p_OS)) + geom_point(size=2) +
      #   scale_x_discrete(limit = seq(0,500,50), labels = paste(seq(0,500,50), sep=",")) +
      #   #  scale_y_discrete(breaks = break_y, labels = break_y) + #c("0", "0.05", "0.10", "0.50", "0.99")) +
      #   coord_trans(y = "log10") +
      #   ggtitle(paste("Cumulative P-value plot for OS under", geneName, "expression")) +
      #   xlab("# of patients") + ylab("P-value(log10)") +
      #   geom_hline(yintercept=0.05, linetype="dashed", color = "red", size=1, show.legend = T)
      # #  geom_vline(xintercept=300, color = "green", size=2, show.legend = T)
      # #
      
      #  ++++ Beside KM P-value, contingency P-value is also important ! ++++++
      # Calling for contingency table to find significant clinicopathological features
      # contingency P-value
      # **debug(contingencyTCGA) ###
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
    #x cases_OS[, 1] <- exp_geneName_sorted$ix # unsorted rank
    #x cases_RFS[, 1] <-exp_geneName_sorted$ix
    
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
    
    
    
    # # {5c.Auto cutoff1 => surv_OS1} => select a cutoff=> please choice ONE from OS (in LUAD) or RFS ####
    # auto pickup a cutoff1, one out of cutoff, from OS/RFS tables, which has lowest P-value, to run directly then export
    OS_pvalue <- OS[OS$p_OS <= 0.05, ] # with rank registration
    # [1] cases_OS p_OS     rank     exp  
    cohort_n <- nrow(oscc) # length(oscc) is columne number (variable number)
    #cohort_n <- (surv_OS$n[1]+surv_OS$n[2])
    # case50_n is the cutoff1 by cases number in unsorted.
    if (nrow(OS_pvalue) > 0) {
      p_OS0 <- min(OS_pvalue$p_OS)
      case50_n <- OS_pvalue[which.min(OS_pvalue$p_OS), 1] # there is a hit; original rank; cases_OS
      #  cutoff1 <- exp_geneName[case50_n] # unsorted, original rank###
      cutoff1 <- OS_pvalue[which.min(OS_pvalue$p_OS), 4] # take it's cutoff value (exp of RNAseq)
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
    
    
    
    
    # OS P-value plot ####
    # #{
    # cumulative P value plot/curves
    library(graphics) 
    library(ggplot2) # http://www.cookbook-r.com/Graphs/Legends_(ggplot2)/
    # sorting issue ##
    attach(OS)
    OS_sorted <- OS[order(p_OS), ]
    #asc <- order(exp_geneName, method="radix")
    detach(OS)
    # aes(x=cases_OS)
    
    
    # # set breaks on y axis
    # if (max(OS$p_OS) <= 0.01) {
    #   break_y <- c(0, 0.00001, 0.001, 0.01, 0.05)
    # } else {
    #   break_y <- c(0, 0.05, 0.1, 0.5)  
    # } #c(0, 0.05, 0.10, 0.50, 0.99)
    
    # by auto pickup a cutoff1 point: case50, p_OS0
    g1<- subset(OS, cases_OS==case50_n) #(x=case50_n, y=p_OS0)
    
    if (g1$p_OS<=0.05) {
      ggplot(OS, aes(x=cases_OS, y=p_OS)) + geom_point(size=2) +
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
      ggplot(OS, aes(x=cases_OS, y=p_OS)) + geom_point(size=2) +
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
    
#View(OS[OS$p_OS <= 0.05,1:2])
    #
    #plot(cases_OS, p_OS, type="p", log="y", main=paste("Cumulative P value curves for", geneName), xlab="# of patients", ylab="P-value")
    #axis(side=4, at=c(0.01, 0.05, 0.2, 0.5, 0.99))
    #abline(h=0.05, lty=2, col="red")
    #abline(v=320, col="green")
    #lines(x=c(cases_OS[1], cases_OS[length(cases_OS)]), y=c(0.05, 0.05), col = "red")
    
    ##
    
    # plot(cases_RFS, p_RFS, type="p")
    # library(graphics)
    # # cumulative P value curves
    # library(ggplot2)
    
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

    
    print(paste("Generating P-Value plot of", geneName, "(", which(whole_genome==geneName), ")"))
    
    # [done cutoff finder] ++++++++++++++++++++++++
    # 
    
    
    
    
    
    
    # 5d.[optimized Statistic procedure] [Option-Cmd + E] for running to the end ####
    # The features from column L to R; running to the end of R2Excel export +++++++++++
    # The correlation of gene expression and clinical features  ### contingency tables
    # if (!require(pkg)){ 
    #   install.packages(pkg) 
    # } # Install package automatically if not there
    
    
    

    
    
    
    # Statistics for osccCleanNA by cutoff1
    # KM plot
    ## (6) KM survival curves in a specific cohort ###
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
    # [**binomial of PMM1_median] resume the correlation table 2 ###
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
    
    
    # KM survival curves ####
    #OS1
    # cancer type shouble be defined at TCGA_cohort <- "LUAD" "HNSC"
    mysurv <- Surv(osccCleanNA$OS..months._from.biopsy, osccCleanNA$OS_IND==1) #1==dead
    # Test for difference (log-rank test) between groups (by PMM1_median 0 vs 1)
    tryCatch(surv_OS1 <- survdiff(mysurv ~ as.vector(osccCleanNA[, osccCleanNAM_pos], mode="numeric"), data=osccCleanNA), error = function(e) return(NA)) # PMM1 high or low
    # pchisq gives the distribution function
    p_OS1 <- format(pchisq(surv_OS1$chisq, length(surv_OS1$n)-1, lower.tail = FALSE), digits=3)
    # (p_OS1 == p_OS0) is TRUE
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
    
    # [survfit] - Kaplan-Meier curve: P-Value
    # confidence intervals as log hazard or log(-log(survival))
    OS.km <- survfit(mysurv ~ as.vector(unlist(osccCleanNA[, osccCleanNAM_pos]), mode="numeric"), data=osccCleanNA, type= "kaplan-meier", conf.type = "log-log")
    # 365.25 days in a year for xscale => 3650 days for 10 years
    # 12 months per year for 5 years => 60 months
    
    plot(OS.km, lty=1, xscale=12, xmax=60, col=c("blue","red"), sub=paste("optimized P-Value =", p_OS1), main=paste("OS in TCGA", TCGA_cohort, "(n=", surv_OS1$n[1]+surv_OS1$n[2],")/", geneName), ylab="Percent Survival", xlab="Years")
    legend("topright", legend=c(paste("low(",surv_OS1$n[1], ")"), paste("high(",surv_OS1$n[2], ")")), lty=1:1, col=c("blue","red"))
    #plot(OS.km, lty=1, col=c("blue","red"), sub="p= 0.816", main="OS in osccCleanNA(n=505)/gene level", ylab="Percent Survival", xlab="Days")
    #legend("topright", legend=c('Low', 'High'), lty=1:2, col=c("blue","red"))
    # summary(OS.km, times = seq(0, 3000, 100))
    print(paste("Generating optimized KM plot of", geneName, "(", which(whole_genome==geneName), ")"))

    
    # *** important returns
    # print("return_8")
    return(list(osccCleanNA, case50_n, cutoff1, surv_OS1, p_OS1, OS, RFS, OS_pvalue))

} # end of cutoFinder_func




