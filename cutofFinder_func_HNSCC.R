# # [2019/05/16] Start recoding/debugging, if any; since [2018/10/?], actually
# # on instance-4 for HNSCC
# # Loading .Rda ###
# # Resume:[5a.START]
# geneName <- "TXNDC11"
# # Start the survival analysis for each individual gene
 
cutofFinder_func <- function(geneName, marginTag) {
# marginTag == "_marginS_" or == "_marginFree_" or == "_marginPlus_"
# geneName

# to assembly osccCleanNA data.frame
#load(file="TMU_TA51BCDE_T70_clinical_fullFeatures13_dichotomized.Rda") # as oscc (without margin feature)
# or ***
load(file="~/R/HNSCC.clinical.RNAseq.Fire.Rda") # n=521, as clean6_oscc
#load(file="LUAD_T522_clinical_fullFeatures11_dichotomized.Rda") # as oscc with "margin" features
#, negative n=348 while positive n=18 and 156:NaN)
oscc <- clean6_oscc
# included $$.clinico_mRNA.Fire for survival analysis
# #
# *** (We don't need to load mRNA by each geneName)
# if (is.na(tryCatch(load(file=paste("LUAD.mRNA.Exp.", geneName, ".Fire.Rda", sep="")), error = function(e) return(NA)))) {return(3)} # there is NO such file
# # as LUAD.mRNA.Exp.Fire
# # there is NO such file: "LUAD.mRNA.Exp.ZFP91.CNTF.Fire.Rda" %in% dir() => False
# #
# # "sample_type"          # 59 NT(Solid Tissue Normal) 515 TP(Primary Solid Tumor) or 2 TR(Recurrent Solid Tumor)
# 
# LUAD_T_mRNA.Fire <- LUAD.mRNA.Exp.Fire[LUAD.mRNA.Exp.Fire$sample_type %in% c("TP"), c("tcga_participant_barcode","z.score")] # n=517 -2 = 515
# # removing duplicated ID %in%  c("TCGA-50-5066","TCGA-50-5946") => TR recurrent sample
# # if NULL or more than 30% NA in RNAseq (z.score) => dischard this gene
# if (length(LUAD_T_mRNA.Fire) == 0) {return(3)}
# if (as.data.frame(colMeans(is.na(LUAD_T_mRNA.Fire)) > 0.3)[2,1] == T) {return(3)}
# 
# # inner join by merge
# if (is.na(tryCatch(LUAD.clinico_mRNA.Fire <- merge(oscc, LUAD_T_mRNA.Fire, by="tcga_participant_barcode") , error = function(e) return(NA)))) {return(5)} # merge error
# #n=515 with sample_type "TP", excluding "NT normal tissue" or "TR"

# >generate HNSCC.clinico_mRNA.Fire ####
colnumberRNA <- which(colnames(oscc)==paste("z.score_", geneName, sep=""))
if (identical(colnumberRNA, integer(0))) {return(3)} # there is no 19867 "ZFP91.CNTF" RNAseq data, skip to next gene

HNSCC.clinico_mRNA.Fire <- oscc[, c(1:8, 9, 11, 13, 14, colnumberRNA)] # append RNAseq data column
# rename all of colnames as HNSCC's osccT
# Error in names(x) <- value : 
colnames(HNSCC.clinico_mRNA.Fire) <- c("Unique.ID","Gender","age.at.diagnosis",
                                      "T","N",
                                      "M","stage_2","margin",
                                      "OS_IND", "OS..months._from.biopsy",
                                      "RFS_IND", "RFS..months._from.op", "H.score_T")
# n=521 in HNSCC
oscc <- HNSCC.clinico_mRNA.Fire # starting analysis with "oscc"
# oscc$H.score_T as LUAD.mRNA.Exp.Fire$z.score; expression level: H.score_T as RNAseq z.score
# *** check % RNAseq of a cohort is available in particular gene; e.x. XKRY, 100% is NaN
# Pipes %>%: df %>% map_dbl(mean)
# library(AMR) # for freq()
freq_oscc <- oscc %>% freq("H.score_T", header=F) # https://www.dummies.com/programming/r/how-to-read-the-output-of-str-for-lists-in-r/
#if (is.na(freq_oscc)) {return(5)}
n_freq_oscc <- as.numeric(strsplit(capture.output(str(freq_oscc))[11], " int ")[[1]][2]) # no of NaN
if (n_freq_oscc/nrow(oscc) >= 0.3) {return(5)} # if NaN% > 30%
#table(is.na(oscc$H.score_T))[2] >= nrow(oscc) * 0.7

## a dummy universal variable for binomial (high/low) oscc$geneName_median, all are zero
oscc$PMM1_median <-(oscc$H.score_T >= median(oscc$H.score_T, na.rm=T)) +0 # higher 1 or lower 0
osccM_pos <- which(colnames(oscc) == "PMM1_median") # at column 14
#
# [2019/06/11]
# > check complete.cases, data cleaning ####
osccClean <- oscc
# commonFeatures <- c(2:6, 8) # common features: gender, age, margin and TNM :-)
commonFeatures <- which(colnames(oscc) %in% c("Gender", "age.at.diagnosis", "T", "N", "M", "margin")) #  6 essential features
osccCleanNA <- osccClean[complete.cases(osccClean[, commonFeatures]), ] # copy all features but remove NaN entities in 6 essential features for their "completeness"
# n=427

# *** addNA for counting all NA (e.g. there is 0, no 1, in margin-free cohort) in "M" "stage_2" or "margin", as a level of factor
# the “pathological” case of two different kinds of NAs which are treated differently: exclude = if (useNA == "no") c(NA, NaN)
# "unusual NA comes from addNA() as factor
# table(osccCleanNA$margin)
#   0   1 (NaN) two level factor 
#328  99   (0)
osccCleanNA$margin <- addNA(osccCleanNA$margin, ifany=F) # ifany=="always" add NA as a factor
#   0    1 <NA> => 3 levels of factor in "margin"
# str(osccCleanNA$margin)
# Factor w/ 3 levels "0","1",NA: 1 2 2 1 1 2 1 2 1 1 ...


oscc_n427 <- osccCleanNA # n=427; removal of NA cases; spare data
# ***_marginS_: margin positive and negative cohort (n=427) ####

osccCleanNA <- oscc_n427

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


#
#
## >5b. (repeat100) Cutoff finder [osccCleanNA] ####
oscc <- oscc0 <- osccCleanNA # syncrhonize them

# checking the completeness of each row: (i.e. rows without NaN)
which(complete.cases(oscc$H.score_T)==F) # if no NaN -> 0
# 24: there is one NaN on H.score_T (at TCGA-BA-A6DF)
which(complete.cases(oscc$OS_IND)==F) # if no NaN -> 0

# *** column 9 should be OS_IND
which(complete.cases(oscc[oscc$OS_IND==1, 9])==F) #OS_IND ==1, death event (dead) => result 0, if no NaN



# # (skipped, if HNSCC RFS is copied from OS)
# which(complete.cases(oscc$RFS_IND)==F) # RFS_IND
# # which(complete.cases(oscc$RFS..months._from.op)==F) #n=103; it may be imputed from OS time (oscc$OS..months._from.biopsy)
# osccCleanNA_RFS <- oscc[which(complete.cases(oscc$RFS..months._from.op)==F), ]
# osccCleanNA_RFS$RFS..months._from.op <- osccCleanNA_RFS[osccCleanNA_RFS$RFS_IND==1, 9] # imputed from OS..months._from.biopsy
#

    
    # start 100 round ####
    #     # 100 slicing is NOT right !
    # => find cutoff1, generate OS P-Value plot and KM plot.
    # oscc <- cbind(oscc0, osccCleanNA$H.score_T) # expression(IHC score) of PMM1
    # {debug
    # find_repeat <- 100 # searching 100 slices in the interval
    cohort_n <- nrow(oscc) # n=328 or 427
    # }debug
    oscc0_pos <- which(colnames(oscc0) == "H.score_T") # oscc0$H.score_T as colunm 13
    
    exp_geneName <- t(oscc0[, oscc0_pos]) # horizontal vector of RNAseq, why t?
    #exp_geneName <- oscc0[, oscc0_pos] # vertical vector
    # $x is RNAseq be sorted, $x is its original position
    # sorting RNAseq for 100 slicing ##
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
    # slice by cases
    cutoff_n <- round(quantile(c(1:cohort_n), c(0.30,0.70))) # 30 percentile, 70 percentile
    cutoff_n[2] <- cutoff_n[2] -1 # from 99 to 230 in  328 cases; 129 to 298 in 427 cases
    # *debug, error and stop on "XKRY"(19642)#
    # (ok)Error in quantile.default(exp_geneName, c(0.3, 0.7)) : 
    #   missing values and NaN's not allowed if 'na.rm' is FALSE
    # (OK)Error during wrapup: names() applied to a non-vector
    #debug
    
    
    # (for run100) P-value according to KM survival analysis (alone) ####
    for (run100 in seq(cutoff_n[1], cutoff_n[2])){ 
      #browser()
      print(paste("cutoff run100=", run100, geneName, sep=" "))
      # sorted (by RNAseq) from 129 to 298 in  427 cases of HNSCC
      # use oscc0 for "positioning"; 
      # oscc is the target dataset to be analysed here.
    
      # exp_geneName_sorted$x[i], cutoff100 is the cutoff value of RNAseq in rank i
      cutoff100 <- exp_geneName_sorted$x[run100+1] #*** shift here (column x, ix)
      
      # Binominal H.Score_T (RNAseq) by exp_geneName_sorted$x[i]) #RNAseq cutoff
      oscc[osccM_pos] <- (oscc0[oscc0_pos] >= cutoff100) +0 # oscc$PMM1_median <- (oscc0$PMM1 >= i) +0
    #  R4> table(oscc$PMM1_median) # when run100 == 99
    #  0   1 
    #  99 229

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
      
      
      
      
      
      #  ++++ Beside KM P-value, contingency P-value is also important ! ####
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
    OS <- OS[complete.cases(OS$p_OS), ] # removal of NA, keep n=102 or 170
    
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
    
    
    
    # # {5c.Auto cutoff1 => surv_OS1} => select a cutoff=> please choice ONE from OS or RFS ####
    # auto pickup a cutoff1, one out of cutoff, from OS/RFS tables, which has lowest P-value, to run directly then export
    OS_pvalue <- OS[OS$p_OS <= 0.05, ] # with rank registration
    OS_pvalue <- OS_pvalue[complete.cases(OS_pvalue$exp), ]
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
      case50_n <- round(cohort_n * default50)
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
    
    if (g1$p_OS <= 0.05) {
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
    
    
    
    
    
    
    # >5d.[optimized Statistic procedure] [Option-Cmd + E] for running to the end ####
    # The features from column L to R; running to the end of R2Excel export +++++++++++
    # The correlation of gene expression and clinical features  ### contingency tables
    # if (!require(pkg)){ 
    #   install.packages(pkg) 
    # } # Install package automatically if not there
    
    
    

    
    
    
    # Statistics for osccCleanNA by cutoff1
    # KM plot
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
    # [**binomial of PMM1_median] resume the correlation table 2 ###
    osccCleanNA[osccCleanNAM_pos] <- (osccCleanNA[, osccCleanNA_pos] >= cutoff1) +0 # binomial after osccCleanNA 
    # one_group issue
    #range(osccCleanNA$PMM1_median, na.rm=T) => try the other cutoff1
    library(dplyr) #nth
    ij <- range(osccCleanNA$PMM1_median, na.rm=T)
    OS_pvalue <- OS_pvalue[order(OS_pvalue$p_OS), ] # sorted from small to big
    # #x find cutoff1 smaller 
    # if (ij[1]==ij[2] & ij[2]==0) {
    #   nth_pvalue <- 1
    #   while (OS_pvalue$exp[nth_pvalue]==cutoff1) {
    #     
    #   }
    #   osccCleanNA[osccCleanNAM_pos] <- (osccCleanNA[, osccCleanNA_pos] >= cutoff1) +0
    # }
    # # xfind cutoff1 bigger
      #if (ij[1]==ij[2] & nrow(OS_pvalue)==1
      if (ij[1]==ij[2]) {
        flag_newCutoff1 <- FALSE
        for (nth_pvalue in c(1:nrow(OS_pvalue))) {
          if (OS_pvalue$exp[nth_pvalue]!=cutoff1) {
            cutoff1 <- OS_pvalue$exp[nth_pvalue]
            flag_newCutoff1 <- TRUE # find a new cutoff1
            break
          }
        #if (nth_pvalue > nrow(OS_pvalue)) {return("skip")}
      }
      if (flag_newCutoff1==TRUE) {
        osccCleanNA[osccCleanNAM_pos] <- (osccCleanNA[, osccCleanNA_pos] >= cutoff1) +0 # 266 vs 160 now
      } else if (flag_newCutoff1==FALSE) {return("skip")} # we skip this gene
    }
    
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




