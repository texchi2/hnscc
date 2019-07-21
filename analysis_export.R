# [Results]  ####
#= Analysis of output .Rda of _marginS_ on GCE i-4 or _marginFree_ on GCE i-4Free
# or i-4Plus
# >*** [choice ONE]: _marginFree_ or _marginS_ loading from .Rda
#marginTag <- "_marginS_" #at ./marginS
#marginTag <- "_marginFree_" #at ./marginFree
#marginTag <- "_marginPlus_" #at ./marginPlus
raw <- 
  readline("_margin(S)_, _marginFree(F)_ or _margin(P)lus_ -- process run: ")
select_margin <- function(x) {
  switch(x,
         s = "_marginS_",
         f = "_marginFree_",
         p = "_marginPlus_",
         S = "_marginS_",
         #         F = "_marginFree_",
         P = "_marginPlus_",
         stop("Unknown input")
  )
}
marginTag <- select_margin(raw)
rm(raw)
# [2019/07/03] they are stored at ./xlsx => rename as ./marginS
path_ZSWIM2 <- file.path(path_cohort, gsub("_", "", marginTag)) # e.x. marginS

# ZSWIM2_archive1000_20180408_0042_0933.Rda; 9 hours for 1,000 genes to be scanned
# 3 hours for 1,000 genes to be scanned under GCP Rstudio server
# STRAP to SLC8A1
#which(whole_genome==geneName)
#
#_marginS_
# xmerge ZSWIM2a ZSWIM2b and ZSWIM2c => ZSWIM2abc_archive.Rda [2018/06/19]
# para_seg <- c(20499, 13528, 6833, 1) # ZZZ3-> PLCE1, PNMT-> GARNL3, GAR1-> A1BG
#load(file=file.path(path_cohort, "ZSWIM2_archive.Rda")) # as ZSWIM2
#ZSWIM2 <- ZSWIM2abc #(done, merged)
#save(data=ZSWIM2, file=file.path(path_ZSWIM2, "ZSWIM2_archive.Rda")) # or ZSWIM2abc_archive 
#x OR
# _marginFree_
#load(file="ZSWIM2_free_archive.Rda") # as well as ZSWIM2



#(at the first time) Dissection of ZSWIM2 #
#1) KM plot is not at best cutoff, why??
#2) ZSWIM2 should be saved by appending.
#3)
load(file=file.path(path_cohort, "ZSWIM2_archive.Rda")) # x merged abc
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
sum(ZSWIM2_2$Freq)
#[1] 20500  # scan completely

# reture(0): ok for analysis
text_pie <- paste("Workable genes \n in ", round(ZSWIM2_2$Freq[1]/sum(ZSWIM2_2$Freq)*100), " %") # _marginS_ about 45.9% usable RNAseq data
pie(ZSWIM2_2$Freq, labels=ZSWIM2_2$Var1, main=text_pie)

# deal with return(1), return(2) and return(3) errors, try to solve it
# _marginS_
n_percent_Bonferroni <- ZSWIM2_2$Freq[1]/sum(ZSWIM2_2$Freq) # 0.4593171
error01_sample <- which(ZSWIM2$X3==1) #  32.2% error01: "ZSWIM2 skip from contingencyTCGA()"): one group issue in Melt()
error02_sample <- which(ZSWIM2$X3==2) # 14.2% error02: There has only one group in survdiff.
error03_sample <- which(ZSWIM2$X3==3) # 6.85% error03: There has only one group in survdiff.
# done #






# >A summary table with z-score from all HNSCC_survival*.Rda ####
# _marginFree_ or _marginS_ or _marginPlus_ from .Rda

# get_Rda_pvalue <- function(geneName) {
#   load(file=paste("HNSCC_survivalAnalysis_marginS_", geneName, ".Rda", sep=""))
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
#aa <- ZSWIM2_2$Freq[1]  # n=9416
# however, $ ls HNSCC_sur*.Rda | wc -l
# 16017

#{ declare empty data.frame and list
#*
candidate_sample <- data.frame(matrix(data = NA, nrow = aa, ncol = 3)) # for num, gene ID and it's P-value
colnames(candidate_sample) <- c("number", "gene_id", "p_value") # KM OS P-value only (in HNSCC); 
#"number" position in ZSWIM2
#candidate_sample$number <- data.frame(which(ZSWIM2$X3==0)) # gene "number"
candidate_sample$gene_id <- whole_genome[bb:aa] # retrieving gene name from whole_genome; return() at X2
#

#* a list for all cox survival datas
# http://www.cookbook-r.com/Manipulating_data/Converting_between_data_frames_and_contingency_tables/
# column name might be created by variable: e.x. df.sex[,"per"] <- df.sex$count1/sum(df.sex$count1); # df.sex$per
candidate_cox <- replicate(aa, list()) # empty list(), which doesn't need colnames
#data.table(matrix(data = NA, nrow = aa, ncol = 3)) # for num, gene ID and it's P-value
# colnames(candidate_cox) <- c( "Features", "HR",       "P_value_uni",  "sig",      "Features", "HR",      
#                               "P_value_multi",  "sig",      "Features", "P_value_KM",  "sig" )
# #}


setwd(path_ZSWIM2) # set for a while (for ip loop) at ./xlsx
for (ip in (bb:aa)) {
  geneName <- candidate_sample$gene_id[ip]
  #print(paste("At", path_ZSWIM2, "=> (", ip, ")", geneName), sep="")
  #candidate_sample$p_value <- lapply(unlist(candidate_sample$gene_id[bb:aa]), get_Rda_pvalue) # retrieving P-value by gene name from .Rda; return() at X3
  #[choice ONE]: _marginFree_ or _marginS_ loading from .Rda
  #  _marginS_
  
  load_filename <- file.path(paste("HNSCC_survivalAnalysis", marginTag, geneName, ".Rda", sep=""))
  #OR (automatically defined)
  #_marginFree_
  #load_filename <- paste("HNSCC_survivalAnalysis", SFree, geneName, ".Rda", sep="")
  #
  #
  if (load_filename %in% dir()) {
    
    load(file=load_filename)
    #      if (is.na(tryCatch(load(file=load_filename), error = function(e) return(NA)))) {} # load file with error free :-)
    #  load(file=paste("HNSCC_survivalAnalysis_marginS_", geneName, ".Rda", sep=""))
    # load list = c("tableChi1", "tableOS1", "tableRFS1", "OS_pvalue", "RFS_pvalue")
    # a example: geneName <- "TRIP13"
    print(paste(marginTag, "with how many OS P-values: ", nrow(OS_pvalue)))
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
    
    # *** we don't have RFS in HNSCC cohort :-)
    
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
    print(paste(ip, geneName, ": cox features saved; ncol = ", length(candidate_cox[[ip]]))) # a vector in [] is ok for indexing
    # length == 11 => that is correct !
    #}
  }
} # end of ip for loop
# #[2019/06/26] _marginS_ finished at 19:48 (about 3 hours)


#_marginS_ or _marginFree_ by SFree; saving on ./run04_marginS_, files => 17030
save(candidate_sample, candidate_cox, n_percent_Bonferroni, file=file.path(path_ZSWIM2, paste("HNSCC_OS", marginTag, "pvalueKM_candidate_cox.Rda", sep=""))) #ok; with KM_sig and Remark, and Cox HR
setwd(path_cohort) 
#_marginS_ by marginTag
# saved file="HNSCC_OS_marginS_pvalueKM_candidate_cox.Rda" above
#_marginFree_, _marginPlus_by marginTag
#save(candidate_sample, candidate_cox, n_percent_Bonferroni, file="HNSCC_OS_marginFree_pvalueKM_candidate_cox.Rda") #ok; with KM_sig and Remark


# devtools::install_github(c('jeroenooms/jsonlite', 'rstudio/shiny', 'ramnathv/htmlwidgets', 'timelyportfolio/listviewer'))
# library(listviewer)
# jsonedit( candidate_cox )
# # finished (both) processes




## Post1 process _marginS_ ####
##  table1 (KM): candidate_sample(KM) ##
# #sort by p_value (ascending) and cyl (descending)
# 楊老師: Bonferroni correction (adjustment) sets the significance cut-off at α/n. of KM P-value; Cox P-value <=0.05 as well.
# https://www.stat.berkeley.edu/~mgoldman/Section0402.pdf
# newdata <- mtcars[order(mpg, -cyl),]
# >*** [choice ONE]: _marginFree_ or _marginS_ loading from .Rda
#marginTag <- "_marginS_" #at ./marginS
#marginTag <- "_marginFree_" #at ./marginFree
#marginTag <- "_marginPlus_" #at ./marginPlus 
raw <- 
  readline("_margin(S)_, _marginFree(F)_ or _margin(P)lus_ -- process run: ")
select_margin <- function(x) {
  switch(x,
         s = "_marginS_",
         f = "_marginFree_",
         p = "_marginPlus_",
         S = "_marginS_",
#         F = "_marginFree_",
         P = "_marginPlus_",
         stop("Unknown input")
  )
 }
marginTag <- select_margin(raw)
rm(raw) 
#
load(file=file.path(path_ZSWIM2, paste("HNSCC_OS", marginTag, "pvalueKM_candidate_cox.Rda", sep=""))) # as candidate_sample with candidate_cox and n_percent_Bonferroni 
attach(candidate_sample)
HNSCC_OS_marginS_pvalue_sorted <- candidate_sample[order(p_value, -number),] # sorting by order(ascending)
detach(candidate_sample)

attach(HNSCC_OS_marginS_pvalue_sorted)
# ***no Bonferroni correction (adjustment) sets the significance cut-off at α/n 
# {
alpha_HNSCC <- 0.05
Bonferroni_cutoff <- alpha_HNSCC / (LUAD_n * n_percent_Bonferroni)
# => 4.786521e-06 (LUAD run04); => 5.31011e-06 (2019 run01)
# if no Bonferroni: 
#Bonferroni_cutoff <- alpha_HNSCC / 1
# } Bonferroni end

HNSCC_OS_marginS_pvalue005_sorted <- HNSCC_OS_marginS_pvalue_sorted[which(p_value<=alpha_HNSCC & !is.na(p_value)), 1:3]
detach(HNSCC_OS_marginS_pvalue_sorted)  # keeping drawing
# no correction: n=6601 in _marginS_; n=? in _marginFree_; and n=? in _marginPlus_
# of HNSCC_OS_marginS_pvalue005_sorted


#  plot(OS.km, lty=1, col=c("blue","red"), sub=paste("p-value =", p_OS[j]), main=paste("OS in OSCC(n=", surv_OS$n[1]+surv_OS$n[2],")/", geneName, "cutoff=", i), ylab="Percent Survival", xlab="Days")
#  legend("topright", legend=c(paste("low(",surv_OS$n[1], ")"), paste("high(",surv_OS$n[2], ")")), lty=1:2, col=c("blue","red"))
library(stats)
library(scales)
library(minpack.lm)

attach(HNSCC_OS_marginS_pvalue005_sorted)
# reg <- lm(number ~ p_value, data = HNSCC_OS_marginS_pvalue005_sorted)
# abline(reg, col="blue")
# or
HNSCC_OS_marginS_pvalue005_sorted$z_score <- scale(number, scale = T) # "number" frequency to z-scores ("Z" because the normal distribution is also known as the "Z distribution").
HNSCC_OS_marginS_pvalue005_sorted$number_01 <- scales:::rescale(z_score, to = c(0, 1)) # rescaled to range minnew to maxnew (aka. 0 to 1 for binomial glm)

# number_01 as Z-score
plot(p_value, number_01, type="p", ylab="Z-score", xlab="P-value plot of KM survival analysis", log="x") # log scale x or y
#axis(side=3, at=c(1e-7, 1e-6, 1e-5, 0.01, 0.05)) #1=below, 2=left, 3=above and 4=right
abline(h=0.6, lty=2, col="blue")
abline(v=alpha_HNSCC, lty=2, col="red") # ***as alpha_HNSCC instead of 5.31011e-06

# run a logistic regression model (categorical 0 vs 1)
#g <- glm(number_01 ~ p_value, family=poisson, data=HNSCC_OS_marginS_pvalue005_sorted)
# (in this case, generalized linear model with log link)(link = "log"), poisson distribution
#curve(predict(g, data.frame(p_value = x), type="response"), add=TRUE, col="blue") # draws a curve based on prediction from logistic regression model

# x# run a non-linear regression: non-linear least squares approach (function nls in R) ; A nice feature of non-linear regression in an applied context is that the estimated parameters have a clear interpretation (Vmax in a Michaelis-Menten model is the maximum rate) which would be harder to get using linear models on transformed data for example.
# # the Levenberg-Marquardt algorithm for nonlinear regression
# # "eyeball" plot to set approximate starting values
# # https://stats.stackexchange.com/questions/160552/why-is-nls-giving-me-singular-gradient-matrix-at-initial-parameter-estimates
# # nls "singular gradient matrix at initial parameter estimates" error, using nlsLM instead
# a_start <- 0.9 #param a is the y value when x=0
# b_start <- 2*log(2)/a_start #b is the decay rate
# m <- nlsLM(number_01 ~ p_value, data=HNSCC_OS_marginS_pvalue005_sorted, start=list(a=a_start, b=b_start))
# curve(predict(m, data.frame(p_value = x), type="response", col="green"), add=TRUE, col="blue") # draws a curve based on prediction from logistic regression model
# #lines(p_value, predict(m), lty=2, col="green", lwd=3)

legend("topright", legend=c(paste("Z-score at 0.06"), paste("P-value at ", alpha_HNSCC)), lty=1:2, col=c("blue","red"), cex=0.7) # box and font size
# figure legend: logistic regression, LR, by Generalized linear model, glm
#detach(HNSCC_OS_marginS_pvalue005_sorted)

# then...
# a candidate table, with "number of P-value" under cutoff finding (becoming z-score: bigger, much more significant cutoff "sites")
# => a kind of local minimal or global minimal of curve fitting (Levenberg-Marquardt Optimization)?
# subsetting the candidated genes table, according P-value (Bonferroni_cutoff) and Z-score
# ***after Bonferroni correction => n=26 in _marginS_
# ***after Bonferroni correction => n=33 in _marginFree_
#attach(HNSCC_OS_marginS_pvalue005_sorted)
HNSCC_OS_marginS_pvalue1e_6_zscore0_6 <- HNSCC_OS_marginS_pvalue005_sorted[which(p_value<=alpha_HNSCC & z_score>=0.6), 2:5]
HNSCC_OS_marginS_pvalue1e_6_zscore0_6[, 2] <- signif(HNSCC_OS_marginS_pvalue1e_6_zscore0_6[, 2], 3)
detach(HNSCC_OS_marginS_pvalue005_sorted) # n=17
# > colnames(HNSCC_OS_marginS_pvalue1e_6_zscore0_6)
# [1] "gene_id"   "p_value"   "z_score"   "number_01"
# 
# _marginS_ (ok) as marginTag
save(HNSCC_OS_marginS_pvalue1e_6_zscore0_6, Bonferroni_cutoff, file=file.path(path_ZSWIM2, paste("HNSCC_OS", marginTag, "pvalue1e_6_zscore0_6.Rda", sep=""))) # cutoff by alpha_HNSCC
save(HNSCC_OS_marginS_pvalue005_sorted, file=file.path(path_ZSWIM2, paste("HNSCC_OS", marginTag, "pvalue005_sorted.Rda", sep="")))

# HNSCC_OS_marginFree_pvalue1e_6_zscore0_6 <- HNSCC_OS_marginS_pvalue1e_6_zscore0_6
# save(HNSCC_OS_marginFree_pvalue1e_6_zscore0_6, file="HNSCC_OS_marginFree_pvalue1e_6_zscore0_6.Rda")
# HNSCC_OS_marginFree_pvalue005_sorted <- HNSCC_OS_marginS_pvalue005_sorted
# save(HNSCC_OS_marginFree_pvalue005_sorted, file="HNSCC_OS_marginFree_pvalue005_sorted.Rda")





## Post2 process and add table2 (cox): candidate_cox ####
# "sig" marking for significant P-value (<=0.05)
# [uni_cox_pvalue, uni_HR, uni_sig]
# [multi_cox_pvalue, multi_HR, multi_sig]
# [exp_pvalue]
# to generate HNSCC_OS_marginS_pvalue005KM_sorted_pvalueCox_HR.Rda; n=??

#HNSCC_OS_marginS_pvalue005_sorted # n=6601
load(file=file.path(path_ZSWIM2, paste("HNSCC_OS", marginTag, "pvalue005_sorted.Rda", sep=""))) # as HNSCC_OS_marginS_pvalue005_sorted
load(file=file.path(path_ZSWIM2, paste("HNSCC_OS", marginTag, "pvalueKM_candidate_cox.Rda", sep=""))) # as candidate_cox (a list, Bonferroni_cutoff), since 2018/05/16
# candidate_sample, candidate_cox, n_percent_Bonferroni

# new variable: KM + Cox, n=6601
HNSCC_OS_marginS_pvalue005KM_sorted_pvalueCox_HR <- HNSCC_OS_marginS_pvalue005_sorted
#dataframe[,"newName"] <- NA # add more named columns (NA)
HNSCC_OS_marginS_pvalue005KM_sorted_pvalueCox_HR[,c("uni_HR", "uni_P_value", "uni_sig","multi_HR", "multi_P_value", "multi_sig")] <- NA
#HNSCC_OS_marginS_pvalue005_sorted[,c(colnames(candidate_cox[[1]][8, c(2:4, 6:8)]))] <- NA
#colnames(HNSCC_OS_marginS_pvalue005KM_sorted_pvalueCox_HR) <- c(...."uni_HR", "uni_P_value", "uni_sig",
#                                                "multi_HR", "multi_P_value", "multi_sig")

#attach(HNSCC_OS_marginS_pvalue005KM_sorted_pvalueCox_HR)
#ipp <- 0
for (ip in 1:nrow(HNSCC_OS_marginS_pvalue005KM_sorted_pvalueCox_HR)) {
  pos_gene <- which(whole_genome==HNSCC_OS_marginS_pvalue005KM_sorted_pvalueCox_HR$gene_id[ip]) # HNSCC_OS_marginS_pvalue005KM_sorted_pvalueCox_HR$gene_id
  # by geneid ?
  if (length(candidate_cox[[pos_gene]])==11) {
    # reference: candidate_cox_ip, which listing as candidate_cox
    HNSCC_OS_marginS_pvalue005KM_sorted_pvalueCox_HR[ip, 6:11]  <- candidate_cox[[pos_gene]][8, c(2:4, 6:8)] # taking RNAseg: sig marked and P-values, HRs
    #   ipp <- ipp+1
    print(paste(ip, " out of ", nrow(HNSCC_OS_marginS_pvalue005KM_sorted_pvalueCox_HR), " (", whole_genome[pos_gene], ")", sep=""))
    print(paste("..added...", HNSCC_OS_marginS_pvalue005KM_sorted_pvalueCox_HR[ip, c(2,6)], sep=";"))
  }
  #?? error on : [2018-05-18 07:26:42] [error] handle_read_frame error: websocketpp.transport:7 (End of File)
  #"exp_pvalue" is a correlation of gene expression vs features (TNM....): 
  # <- candidate_cox[which(gene_id[ip]==whole_genome)][1:7, c(11)] # correlation
  # colnames <-  c("KM_Features", "KM_P_value",  "KM_sig")
  
}
#detach(HNSCC_OS_marginS_pvalue005KM_sorted_pvalueCox_HR)
# c( "uni_Features", "uni_HR",       "uni_P_value",  "uni_sig",      "multi_Features", "multi_HR",
#"multi_P_value",  "multi_sig",      "KM_Features", "KM_P_value",  "KM_sig")
# add colname as HRs, "uni_cox"/"multi_cox", sig ... for HRs, P-values, sig of RNAseq(z-score)

##
# save table1 + table2 (ok) [2019/07/01], n= 6601
save(HNSCC_OS_marginS_pvalue005KM_sorted_pvalueCox_HR, file=file.path(path_ZSWIM2, paste("HNSCC_OS", marginTag, "pvalue005KM_sorted_pvalueCox_HR.Rda", sep="")))

#save(list(HNSCC_OS_marginS_pvalue1e_6_zscore0_6, ???))
# x[i], or might be x[i:j]
# x[i, j]
# x[[i]]; x[[expr]]; it can NOT be x[[i:j]]
# x[[i, j]]
# x$a
# x$"a"




# ** Pickup all significant genes list -> HNSCC_OS_marginS_THREE_pvalue005 ####
# _marginS_ on [2019/07/02]
# > colnames(HNSCC_OS_marginS_pvalue005KM_sorted_pvalueCox_HR)
# [1] "number"        "gene_id"       "p_value"      
# [4] "z_score"       "number_01"     "uni_HR"       
# [7] "uni_P_value"   "uni_sig"       "multi_HR"     
# [10] "multi_P_value" "multi_sig"
HNSCC_OS_marginS_THREE_pvalue005 <- subset(HNSCC_OS_marginS_pvalue005KM_sorted_pvalueCox_HR, 
                                           select=c(number, gene_id, p_value, uni_HR, uni_P_value, multi_HR, multi_P_value))
#HNSCC_OS_marginS_THREE_pvalue005 <- subset(HNSCC_OS_marginS_pvalue005KM_sorted_pvalueCox_HR, (uni_P_value <= 0.05) & (multi_P_value <= 0.05), 
#                                           select=c(number, gene_id, p_value, uni_HR, uni_P_value, multi_HR, multi_P_value))
# n=1408, KM P-value only (NOT by Bonferroni_cutoff)
save(HNSCC_OS_marginS_THREE_pvalue005, Bonferroni_cutoff, file=file.path(path_ZSWIM2, paste("HNSCC_OS", marginTag, "THREE_pvalue005.Rda", sep=""))) 
# as HNSCC_OS_marginS_THREE_pvalue005, Bonferroni_cutoff <- 5.31011e-06

##*** RR>reproducible research resume: ####
#load(file=file.path(path_ZSWIM2, paste("HNSCC_OS", marginTag, "THREE_pvalue005.Rda", sep="")))
## <<<

# >pubmed.mineR ####
load(file=file.path(path_ZSWIM2, paste("HNSCC_OS", marginTag, "THREE_pvalue005.Rda", sep="")))
# [GDC data portal]
# https://www.rdocumentation.org/packages/pubmed.mineR/versions/1.0.1
# Excluding Cancer Gene Census (such as TRIM24, ...) from GDC data portal: d/l Cancer Gene Census (575 genes)
# [2019/07/17] downloaded, n=575
# load genes.2019-07-17.json
install.packages("jsonlite")
library(jsonlite)
cc_gene575 <- fromJSON("~/R/genes.2019-07-17.json", simplifyDataFrame=TRUE)
# subsetting by omitting row
cc_Index <- which(HNSCC_OS_marginS_THREE_pvalue005$gene_id %in% cc_gene575$symbol)
View(HNSCC_OS_marginS_THREE_pvalue005$gene_id[cc_Index]) # n=213
HNSCC_OS_marginS_THREE_pvalue005_noCancerGene <- HNSCC_OS_marginS_THREE_pvalue005[-cc_Index, ] # removal
# now n=6388

# ***Is there any updated HNSCC related gene from publication? ####
# go to embase

# [PubMed] search keywoard:"((((marker[Abstract] OR biomarker[Abstract])) AND (HNSCC OR HNSC OR "oral"))) AND (prognostic[Title] OR prognosis[Title])"
install.packages("pubmed.mineR") #  two formats (text and XML) from PubMed
# scp /Users/apple/Downloads/pmc_result_HNSCC_genes1486_MEDLINE.txt tex@35.201.169.0:~/R/
# *** export MEDLINE format (with abstract)
library(pubmed.mineR) # http://ramuigib.blogspot.com
install.packages(c("rentrez", "XML"))
library(rentrez)
library(XML)
# gene_atomization(), official_fn(gene_atomization, abs, filename, terms), 
#abs <- xmlreadabs("~/R/pmc_result_HNSCC_genes1486.xml") #an S4 object of class Abstracts containing journals, abstracts and PMID
#abs <- read.table("~/R/pmc_result_HNSCC_genes1486.xml", sep="|") #an S4 object of class Abstracts containing journals, abstracts and PMID
abs <- readabs("~/R/pmc_result_HNSCC_genes1486_MEDLINE.txt") # Abstracts <- .txt; the slot operator @ 
# # 5.5Mb; objects of class "Abstracts" have 3 slots, Journal, Abstract, and PMID
# readabs(), xmlgene_atomizations(m)
# x pub_hnscc_genes <- xmlgene_atomizations(abs)
#pub_hnscc_genes <- gene_atomization(abs) # n=615; however, "FAU" "T" "GC" "HR" are not gene symbol
## (there is misleading data)
#pub_Index <- which(pub_hnscc_genes[ , 1] %in% whole_genome) # curation1; n=611
#pub_hnscc_genes <- pub_hnscc_genes[pub_Index, ] # n=611
# https://www.rdocumentation.org/packages/pubmed.mineR/versions/1.0.1
# https://omictools.com/pubtator-tool PubTator
#x curation2 <- mapply(Give_Sentences, pub_hnscc_genes[ , 2], abs@Abstract) # or searchabsT
#x curation2 <- mapply(searchabsT, abs, include=pub_hnscc_genes[ , 2]) # or searchabsT
#  getabs(abs, pub_hnscc_genes[1, 2], FALSE)

pubmed_hnscc_PMID <- entrez_search(db="pubmed", 
                                   term="((((((((marker[Abstract] OR biomarker[Abstract])) AND (prognostic[Title] OR prognosis[Title]))))) AND (((Upper Aerodigestive) OR HNSCC OR (oral cancer)))))) NOT ((salivary OR bladder OR breast OR bone OR kidney OR leukoplakia OR dysplasia))", retmax=2000, use_history=TRUE)
# n=1325
# x=> refinement with ((D006258[MeSH Major Topic]) OR D006258[MeSH Terms]) <= Head and Neck Neoplasms
#x pubmed_hnscc_PMID <- entrez_search(db="pubmed", term="((((marker[Abstract] OR biomarker[Abstract])) AND ((D006258[MeSH Major Topic]) OR D006258[MeSH Terms]) AND (prognostic[Title] OR prognosis[Title])", retmax=5000)
# MESH is not universal available :-)

# paper_links <- entrez_link(dbfrom="pubmed", id=25500142, cmd="llinks")
# paper_links$linkouts; linkout_urls(paper_links)
# entrez_fetch(); entrez_summary()
# error: download abstracts (HTTP failure 414, the request is too large. For large requests, try using "web history" as described in the rentrez tutorial)
# "web history" { , use_history=TRUE => web_history object
#upload_ncbi <- entrez_post(db="pubmed", id = pubmed_hnscc_PMID$ids)
#upload_ncbi

pubmed_hnscc_abs <- entrez_fetch(db="pubmed", web_history=pubmed_hnscc_PMID$web_history,
                       rettype="xml", retmax=2000, parsed = T)
#R4> str(pubmed_hnscc_abs)
#Classes 'XMLInternalDocument', 'XMLAbstractDocument' <externalptr> 
# *.txt generated
# /Users/apple/Downloads/dataout.txt => gene without abstract available
# /Users/apple/Downloads/newabs.txt => abstracts
# /Users/apple/Downloads/result.txt => ditto with pubmed_hnscc_PMID_list **
# /Users/apple/Downloads/table.txt => genes list
# }

# x pubmed_hnscc_abs <- entrez_fetch(db = "pubmed", id = pubmed_hnscc_PMID$ids,
#                rettype = "xml", parsed = T)

# >text annotation uisng online PubTator (rentrez::pubtator_function)
#pub_hnscc_PMID <- lapply(abs@PMID, pubtator_function) # Gene, Chemical, Disease and PMID
# ** ditto with result.txt from entrez_fetch
pubmed_hnscc_PMID_list <- lapply(pubmed_hnscc_PMID$ids, pubtator_function) # Gene, Chemical, Disease and PMID
# => n=1325; list of Genes, Diseases, ... and PMID
# e.x. 
# R4> pubmed_hnscc_PMID_list[[555]]$Diseases
#$Diseases
#[1] "cancer>MESH:D009369"                  
#[2] "hematopoiesis>MESH:C536227"           
#[3] "inflammation>MESH:D007249"            
#[4] "tumor>MESH:D009369"                   
#[5] "squamous cell carcinoma>MESH:D002294" 
#[6] "squamous cell carcinomas>MESH:D002294"
#[7] "tumors>MESH:D009369" 
# e.x. pubmed_hnscc_PMID[[555]]$Genes is TATA-binding protein>6908 => Entrez Gene: 6908
#pub_hnscc_table <- pubtator_result_list_to_table(pub_hnscc_PMID) # to exclude HNSCC "disease"
# check genes with HNSCC OR HNSC OR "oral cancer" related
pubmed_hnscc_df <- as.data.frame(matrix(NA, nrow=length(pubmed_hnscc_PMID_list), ncol=3))
colnames(pubmed_hnscc_df) <- c("PMID", "Diseases", "Genes")
for (pubtator_i in 1: length(pubmed_hnscc_PMID_list)) {
  if (pubmed_hnscc_PMID_list[[pubtator_i]] == " No Data ") {next}
  else {pubmed_hnscc_df[pubtator_i, 1] <- pubmed_hnscc_PMID_list[[pubtator_i]]$PMID[1]}

  tryCatch({
    pubmed_hnscc_df[pubtator_i, 2] <- pubmed_hnscc_PMID_list[[pubtator_i]]$Diseases[1]
    pubmed_hnscc_df[pubtator_i, 3] <- pubmed_hnscc_PMID_list[[pubtator_i]]$Genes[1]
    }, error = function(err) {
      # error handler picks up where error was generated
      #print(paste("MY_ERROR:  ", err))
      return()
    }) # END tryCatch
}

pubmed_hnscc_abstractFull <- read.csv(file.path(path_ZSWIM2, "newabs.txt"), header=T) 
# /Users/apple/Downloads/newabs.txt => abstracts
# 29453876: ameloblastoma, Genes of FGF2, Bcl-2
# excluding (salivary OR bladder OR breast OR bone OR kidney OR leukoplakia OR dysplasia) 

#[embase]
# >embase search then exported .csv ####
# https://dev.elsevier.com/documentation/EmbaseAPI.wadl
# https://www.rdocumentation.org/packages/rscopus/versions/0.6.3/topics/scopus_search
#https://cran.r-project.org/web/packages/rscopus/rscopus.pdf
#install.packages("rscopus")
library(rscopus)
library(pubmed.mineR) # http://ramuigib.blogspot.com
#library(RISmed) # https://amunategui.github.io/pubmed-query/; 2017 no more read.ris()
# https://cran.r-project.org/web/packages/RefManageR/RefManageR.pdf

#embase_retrieval(id, identifier = c("lui", "pii", "doi", "embase",
#                                    "pubmed_id", "medline"), http_end = NULL, ...)
# Type of search. See https://dev.elsevier.com/api_docs.html
# [2019/07/20] my elsevier api_key = "9925b63f77af6f11b74b38a3d6a80f41" from web site
# search in embase: (prognosis:ti OR prognostic:ti OR 'tumor marker') AND cancer:ti AND 'head and neck squamous cell carcinoma'
# download records.csv
res_df <- read.csv(file=file.choose(), header=T) # with abstract, PUI as embase
em_hnscc_PMID_list <- lapply(res_df$Medline.PMID, pubtator_function) # list of Gene, Chemical, Disease and PMID
# => n=371 => 134; 503
# e.x. The neuropeptide genes SST, TAC1, HCRT, NPY, and GAL are powerful epigenetic biomarkers in head and neck cancer: a site-specific analysis.
#(PMID:29682090 PMCID:PMC5896056) => The somatostatin (SST)
# parsing and pickup 1st: "PMID", "Diseases", "Genes"
pubmed_hnscc_df <- as.data.frame(matrix(NA, nrow=length(em_hnscc_PMID_list), ncol=3))
colnames(pubmed_hnscc_df) <- c("PMID", "Diseases", "Genes")
for (pubtator_i in 1: length(em_hnscc_PMID_list)) {
  print(pubtator_i)
  if (em_hnscc_PMID_list[[pubtator_i]] == " No Data ") {next}
  else {pubmed_hnscc_df[pubtator_i, 1] <- em_hnscc_PMID_list[[pubtator_i]]$PMID[1]}
  
  tryCatch({
    pubmed_hnscc_df[pubtator_i, 2] <- em_hnscc_PMID_list[[pubtator_i]]$Diseases[1]
    pubmed_hnscc_df[pubtator_i, 3] <- em_hnscc_PMID_list[[pubtator_i]]$Genes[1]
  }, error = function(err) {
    # error handler picks up where error was generated
    #print(paste("MY_ERROR:  ", err))
    return()
  }) # END tryCatch
}
# find published HNSCC genes from Embase
# pubmed_hnscc_df data cleaning: sorting by "Disease", manual curation
pubmed_hnscc_df_sorted <- pubmed_hnscc_df[order(pubmed_hnscc_df$Diseases, na.last=NA), ] # NA removal, n=259
pubmed_hnscc_df_sorted <- pubmed_hnscc_df_sorted[complete.cases(pubmed_hnscc_df_sorted$Genes), ] # n=169 => 64
# n=244
# loading all genes mentioned in each article
pubmed_hnscc_genes <- as.data.frame(matrix(NA, nrow=nrow(pubmed_hnscc_df_sorted)*4, ncol=5))
colnames(pubmed_hnscc_genes) <- c("No", "PMID", "Diseases", "GeneSymbol", "EntrezGene")
No_i <- 1 # a serial number

for (curation_i in 1: nrow(pubmed_hnscc_df_sorted)) {
  print(curation_i)
  No_j <- 1 # index of genes in each article
  # em_hnscc_PMID_list is 1:371 from res_df$Medline.PMID; there is duplicated article (PMID:24565202)
  res_index <- which(res_df$Medline.PMID==pubmed_hnscc_df_sorted$PMID[curation_i]) # always choice res_index[1], first PMID
  for (No_j in 1:length(em_hnscc_PMID_list[[res_index[1]]]$Genes)) {
    pubmed_hnscc_genes[No_i + No_j - 1, ] <- c((No_i + No_j - 1), pubmed_hnscc_df_sorted$PMID[curation_i], 
                                pubmed_hnscc_df_sorted$Diseases[curation_i], 
                                unlist(strsplit(em_hnscc_PMID_list[[res_index[1]]]$Genes[No_j], ">")))
  }
  No_i <- No_i + No_j
# EntrezGene 4313 (MMP2), 4318 (MMP9); 83639 (TEX101) at PMID 23481573
# Cancer Biomark. 2012-2013;12(3):141-8. doi: 10.3233/CBM-130302.
#  [Overexpression of TEX101, a potential novel cancer marker, in head and neck squamous cell carcinoma.]
} # from 1 to 976 entries in [pubmed_hnscc_genes]
# there is NO GLUT4 nor TMSB4X => repeat searching including
#  ('glucose transporter 4' OR 'thymosin beta 4') AND 'hsiao m.':au (yes)
#  'mRNA expression level'
#  if 'Chang W.-M' AND 'hsiao m.':au => 
# resulting pubmed_hnscc_genes, 256 entries
# \Sci Rep. 2017[Parathyroid Hormone-Like Hormone is a
# Poor Prognosis Marker of Head and Neck Cancer and Promotes Cell Growth via RUNX2Regulation.]
pubmed_hnscc_genes <- pubmed_hnscc_genes[complete.cases(pubmed_hnscc_genes$GeneSymbol), ]
# There is n= 785 HNSCC related genes
#  
em_Index <- which(HNSCC_OS_marginS_THREE_pvalue005_noCancerGene$gene_id %in% pubmed_hnscc_genes$GeneSymbol)
# R4> em_Index # n=75
# [1]   26   67  206  284  427  484  530  609  809  971 1020 1067 1137 1178 1193
# [16] 1365 1404 1423 1501 1518 1523 1698 1819 1881 2036 2104 2197 2212 2360 2380
# [31] 2490 2696 2706 2759 2808 2926 2953 2992 3287 3392 3550 3700 3796 3848 3854
# [46] 3938 3985 4076 4095 4125 4213 4246 4352 4360 4604 4667 4671 4677 4731 4742
# [61] 4849 4900 4948 5034 5095 5550 5582 5657 5659 5699 5708 6018 6094 6153 6263
HNSCC_OS_marginS_THREE_pvalue005_noCancerGene <- HNSCC_OS_marginS_THREE_pvalue005_noCancerGene[-em_Index, ] # removal
# 6388 -> n=6313 genes
# Excluding the HNSCC cancer driver genes list: done from Embase
# ?BioXpress*.csv # data from https://hive.biochemistry.gwu.edu/cgi-bin/prd/bioxpress/servlet.cgi
# <done> [2019/07/20] anniversary for discovery of GLUT4 in HNSCC cell line [2014/07/20]


# P-values plot uni_HR, n=6313 ####
HNSCC_OS_marginS_THREE_pvalue005 <- HNSCC_OS_marginS_THREE_pvalue005_noCancerGene
attach(HNSCC_OS_marginS_THREE_pvalue005)
#tiff("Rplot07_cox_uniHR.tiff", units="cm", width=5, height=5, res=300) 
# x# saving as .tiff (by tiff()) => x y axis is too small in .tiff
# adds a loess regression smoother to a scatter plot by smoothScatter
#smoothScatter(p_value, uni_HR, nrpoints = 10000, cex=0.6, xlab="P-value (KM survival)", ylab="Cox's Harzard Ratios (univariate)", bandwidth=c(7, 7), nbin=20)
plot(p_value, uni_HR, type="p", main=paste(TCGA_cohort, "Cox's Harzard Ratios (univariate)"), 
     ylab="Cox's HR", xlab="P-value (KM survival)", log="y", cex=0.1, , cex.main=1.5, cex.lab=1.3, cex.axis=1)
#axis(side=3, at=c(1e-7, 1e-6, 1e-5, 0.01, 0.05)) #1=below, 2=left, 3=above and 4=right
abline(h=2.0, lty=2, col="red") # HR
abline(v=alpha_HNSCC, lty=2, col="blue") #no Bonferroni_cutoff
# run a logistic regression model (categorical 0 vs 1)
#g <- glm(uni_HR ~ p_value, family=poisson, data=HNSCC_OS_marginS_THREE_pvalue005)
# (in this case, generalized linear model with log link)(link = "log"), poisson distribution
#curve(predict(g, data.frame(p_value = x), type="response"), add=TRUE, col="blue") # draws a curve based on prediction from logistic regression model
legend("topright", legend=c(paste("HR at 2.0"), paste("P-value at", alpha_HNSCC)), lty=2:2, col=c("red","blue"), cex=1) # box and font size
# text(500,900, labels=paste(TCGA_cohort, "KM P-value plot"), cex = 0.60) # (0,0) on bottom_left corner
# figure legend: logistic regression, LR, by Generalized linear model, glm
#dev.off()
detach(HNSCC_OS_marginS_THREE_pvalue005)


#plot(HNSCC_OS_marginS_THREE_pvalue005$p_value,  HNSCC_OS_marginS_THREE_pvalue005$uni_HR)
# plot multi_HR, n=1408
attach(HNSCC_OS_marginS_THREE_pvalue005)
plot(p_value, multi_HR, type="p", main=paste(TCGA_cohort, "Cox's Harzard Ratios (multivariate) plot"), 
     ylab="Cox's HR", xlab="P-value (KM survival)", log="y", cex=0.1, , cex.main=1.5, cex.lab=1.3, cex.axis=1)
#axis(side=3, at=c(1e-7, 1e-6, 1e-5, 0.01, 0.05)) #1=below, 2=left, 3=above and 4=right
abline(h=2.0, lty=2, col="red")
abline(v=alpha_HNSCC, lty=2, col="blue")
# run a logistic regression model (categorical 0 vs 1)
#g <- glm(uni_HR ~ p_value, family=poisson, data=HNSCC_OS_marginS_THREE_pvalue005)
# (in this case, generalized linear model with log link)(link = "log"), poisson distribution
#curve(predict(g, data.frame(p_value = x), type="response"), add=TRUE, col="blue") # draws a curve based on prediction from logistic regression model
legend("topright", legend=c(paste("HR at 2.0"), paste("P-value at", alpha_HNSCC)), lty=2:2, col=c("red","blue"), cex=1) # box and font size
# figure legend: logistic regression, LR, by Generalized linear model, glm
detach(HNSCC_OS_marginS_THREE_pvalue005)


## Venn_marginS_: cross matching HR>1 of Uni+Multi, HR<1 of Uni+Multi ####
#library(venn) #https://cran.r-project.org/web/packages/venn/venn.pdf
# venn(x, snames = "", ilabels = FALSE, counts = FALSE, ellipse = FALSE,
#      zcolor = "bw", opacity = 0.3, size = 15, cexil = 0.6, cexsn = 0.85,
#      borders = TRUE, ...)
#library(gplots)
# To get the list of gene present in each Venn compartment we can use the gplots package
#library(gplots) # capture the list of genes from venn





#{ [pickup1] (bad guy genes)#### 
# from HNSCC_OS_marginS_THREE_pvalue005; no Bonferroni_cutoff
# Broader gene candidate (first 100): Cox HR (>1.5 or >=2.5), bad_FC fold change 
# x (uni_P_value <= 0.05) & (multi_P_value <= 0.05)
# Bonferroni_cutoff = 5.31011e-06 is too restricted in this cohort
bad_FC <- 1.8
HNSCC_OS_marginS_uni_CoxHR2p5 <- subset(HNSCC_OS_marginS_THREE_pvalue005, (p_value <= alpha_HNSCC) & (uni_P_value <= 0.05) & (uni_HR >= bad_FC), 
                                        select=c(number, gene_id, p_value, uni_HR, uni_P_value, multi_HR, multi_P_value))
# n=61
#HNSCC_OS_marginS_uni_CoxHR2p5 <- subset(HNSCC_OS_marginS_THREE_pvalue005, (p_value <= Bonferroni_cutoff) & (uni_HR >= 2.5), 
#                                       select=c(number, gene_id, p_value, uni_HR, uni_P_value, multi_HR, multi_P_value))

# multi_HR >=1 or 2.5 # & (uni_P_value <= 0.05) & (multi_P_value <= 0.05) 
HNSCC_OS_marginS_multi_CoxHR2p5 <- subset(HNSCC_OS_marginS_THREE_pvalue005,  (p_value <= alpha_HNSCC) & (multi_P_value <= 0.05) & (multi_HR >= bad_FC), 
                                          select=c(number, gene_id, p_value, uni_HR, uni_P_value, multi_HR, multi_P_value))
# n=71; uni and multi could not be identical gene list :-)


# venn1 diagram of HR>=1 of Uni & Multi ###
#for list of genes by grouping; 
venn_HR2p5 <- list(HNSCC_OS_marginS_uni_CoxHR2p5$gene_id, HNSCC_OS_marginS_multi_CoxHR2p5$gene_id)
# no cutoff by Bonferroni_cutoff
names_HR2p5 <- c(paste("Cox's HR(univariate) >=", bad_FC), paste("Cox's HR(multivariate) >=", bad_FC))
library(gplots)
tmp_bad <- venn(venn_HR2p5, names=names_HR2p5, show.plot=F) #library(gplots); the group count matrix alone
isect_HR2p5 <- attr(tmp_bad, "intersections")
#isect_HR2p5
detach(package:gplots)


# https://www.statmethods.net/advgraphs/parameters.html
library(venn)
# https://cran.r-project.org/web/packages/venn/venn.pdf
tiff("Rplot09_venn_HR1.5.tiff", units="cm", width=5, height=5, res=300) 
# # saving as .tiff (by tiff())
venn(venn_HR2p5, snames=names_HR2p5,
     ilabels = T, counts = T, ellipse = FALSE, zcolor = "red, deeppink", opacity = 0.6, size = 15, cexil = 0.7, cexsn = 0.3, borders = TRUE)
#      predefined colors if "style"
#http://www.stat.columbia.edu/~tzheng/files/Rcolor.pdf
# meta-language 1 0 or -
title <- c(paste(TCGA_cohort, "survival analysis"), paste("(KM P-Value <=", signif(alpha_HNSCC, 3), ")")) #, collapse = "\n")
#coords <- unlist(getCentroid(getZones(venn_HR2p5, snames="uni_CoxHR>=2p5, multi_CoxHR>=2p5")))
# coords[1], coords[2], 
text(500,900, labels = title[1], cex = 0.60) # (0,0) on bottom_left corner
text(500,855, labels = title[2], cex = 0.40) 
# n=47
dev.off() # saving instead of showing

#https://stackoverflow.com/questions/43324180/adding-legend-to-venn-diagram
#legend("top", legend=c("B:multi_Cox HR >=2.5", "A:uni_Cox HR >=2.5"), lty=1:2, col=c("blue","red"), cex=0.7) # box and font size
# figure legend: 



#===
#{ [pickup2]  (good guy genes) ####
# from HNSCC_OS_marginS_THREE_pvalue005
#* Cox HR <0.4 or <0.5 # good_FC <- 0.5

#...
#HNSCC_OS_marginS_uni_CoxHR0p5 <- subset(HNSCC_OS_marginS_THREE_pvalue005, (p_value <= Bonferroni_cutoff) & (uni_P_value <= 0.05) & (multi_P_value <= 0.05) & (uni_HR <0.0), 
#                                       select=c(number, gene_id, p_value, uni_HR, uni_P_value, multi_HR, multi_P_value))
# n=0 while (uni_P_value <= 0.05) & (multi_P_value <= 0.05) & (uni_HR <0.0)
# 
good_FC <- 0.5
HNSCC_OS_marginS_uni_CoxHR0p5 <- subset(HNSCC_OS_marginS_THREE_pvalue005, (p_value <= alpha_HNSCC) & (uni_P_value <= 0.05) & (uni_HR <good_FC), 
                                        select=c(number, gene_id, p_value, uni_HR, uni_P_value, multi_HR, multi_P_value))
print(nrow(HNSCC_OS_marginS_uni_CoxHR0p5))
# n=60
# x R4> HNSCC_OS_marginS_uni_CoxHR0p5$gene_id
# [1] "ZNF557"   "IL19"     "EVPLL"    "ZNF266"   "MYO1H"    "ZNF846"   "MASP1"    "DOT1L"   
# [9] "UTY"      "FAM162B"  "KLRA1"    "FLT3"     "TP53INP1"
# # 

# ...
HNSCC_OS_marginS_multi_CoxHR0p5 <- subset(HNSCC_OS_marginS_THREE_pvalue005, (p_value <= alpha_HNSCC) & (multi_P_value <= 0.05) & (multi_HR <good_FC), 
                                          select=c(number, gene_id, p_value, uni_HR, uni_P_value, multi_HR, multi_P_value))

# n=61



# venn2 diagram of HR < 0.8 of Uni & Multi ###
#for list of genes by grouping; library(gplots)
venn_HR0p5 <- list(HNSCC_OS_marginS_uni_CoxHR0p5$gene_id, HNSCC_OS_marginS_multi_CoxHR0p5$gene_id)
names_HR0p5 <- c(paste("Cox's HR(univariate) <", good_FC), paste("Cox's HR(multivariate) <", good_FC))
library(gplots)
tmp_good <- venn(venn_HR0p5, names=names_HR0p5, show.plot=F) #library(gplots); the group count matrix alone
isect_HR0p5 <- attr(tmp_good, "intersections")
detach(package:gplots)

library(venn)
# venn(venn_HR0p5, snames=names_HR0p5,
#      ilabels = T, counts = T, ellipse = FALSE,  zcolor = "forestgreen, darkolivegreen1", opacity = 0.6, size = 15, cexil = 0.6, cexsn = 0.85, borders = TRUE)
#      predefined colors if "style"
tiff("Rplot09_venn_HR0.5.tiff", units="cm", width=5, height=5, res=300) 
# # saving as .tiff (by tiff())
venn(venn_HR0p5, snames=names_HR0p5,
     ilabels = T, counts = T, ellipse = FALSE, zcolor = "forestgreen, darkolivegreen1", opacity = 0.6, size = 15, cexil = 0.7, cexsn = 0.3, borders = TRUE)
#      predefined colors if "style"
#http://www.stat.columbia.edu/~tzheng/files/Rcolor.pdf
# meta-language 1 0 or -
title <- c(paste(TCGA_cohort, "survival analysis"), paste("(KM P-Value <=", signif(alpha_HNSCC, 3), ")")) #, collapse = "\n")
#coords <- unlist(getCentroid(getZones(venn_HR2p5, snames="uni_CoxHR>=2p5, multi_CoxHR>=2p5")))
# coords[1], coords[2], 
text(500,900, labels = title[1], cex = 0.60) # (0,0) on bottom_left corner
text(500,855, labels = title[2], cex = 0.40) 
# n=47
dev.off() # saving instead of showing

# signif(Bonferroni_cutoff, 3)


#x #library(VennDiagram)
# VENN.LIST <- list(HNSCC_OS_marginS_uni_CoxHR0p5$gene_id, HNSCC_OS_marginS_multi_CoxHR0p5$gene_id)
# #venn.plot <- venn.diagram(VENN.LIST , NULL, fill=c("darkmagenta", "darkblue"), alpha=c(0.5,0.5), cex = 2, cat.fontface=4, category.names=c("A", "B"), main="HNSCC: Cox HR <0.5")
# # To plot the venn diagram we will use the grid.draw() function to plot the venn diagram
# #grid.draw(venn.plot)
# # We can summarize the contents of each venn compartment, as follows:
# # in 1) ConditionA only, 2) ConditionB only, 3) ConditionA & ConditionB
# lapply(isect_HR0p5, head)
#}

# Embase HNSCC genes was removed; saved genes in left, right and intersection
save(isect_HR2p5, isect_HR0p5, 
     file=file.path(path_ZSWIM2, paste(TCGA_cohort, "_OS", marginTag, "candidates_Venn.Rda", sep="")))
#} venn the end






## Export r2excel and .Rda ####
# _marginS_ [2019/07/21]
# sink() for .xlsx export as well :-) https://stackoverflow.com/questions/34038041/how-to-merge-multiple-data-frame-into-one-table-and-export-to-excel
# HNSCC_OS_marginS_candidates_Venn.xlsx: show up Bonferroni_cutoff and 5e-6 (expression style).
# 
load(file=file.path(path_ZSWIM2, "HNSCC_OS_marginS_pvalue1e_6_zscore0_6.Rda")) # in HNSCC_OS_marginS_pvalue1e_6_zscore0_6, cutoff by Bonferroni_cutoff
load(file=file.path(path_ZSWIM2, "HNSCC_OS_marginS_pvalueKM_candidate_cox.Rda")) # in candidate_cox (a list), candidate_sample, candidate_cox, n_percent_Bonferroni
load(file=file.path(path_ZSWIM2, paste(TCGA_cohort, "_OS", marginTag, "candidates_Venn.Rda", sep=""))) # isect_HR2p5, isect_HR0p5

library("xlsx")
library("r2excel")
# Create an Excel workbook. Both .xls and .xlsx file formats can be used.
filenamex <- paste(TCGA_cohort, "_OS", marginTag, "candidates_Venn", ".xlsx", sep = "") 
# "HNSCC_OS_marginS_candidates_Venn.xlsx"
wb <- createWorkbook(type="xlsx")

# Create a sheet in that workbook
sheet <- xlsx::createSheet(wb, sheetName = paste("Survival_candidates"))
# [add data row by row, start from column 2]
#+++++++++++++++++++++++++++++++
## Add paragraph : Author
library("tis") # by Brian Salzer
# today(), arg must be ti, tis, ts, tif, or tifName
select_title <- function(x) {
  switch(x,
         "_marginS_" = "The cohort with surgical margins status (0 or 1)",
         "_marginFree_" = "The cohort with surgical margins free from tumor (0)",
         "_marginPlus_" = "The cohort with surgical margins involving tumor (1)",
         stop("Unknown input")
  )
}
title_candidates_Venn_xlsx <- select_title(marginTag)
author <- paste("Reported by Tex Li-Hsing Chi. \n",
                "tex@gate.sinica.edu.tw \n", title_candidates_Venn_xlsx, "\n", Sys.Date(), sep="")
xlsx.addParagraph(wb, sheet, value=author, isItalic=TRUE, colSpan=5, 
                  rowSpan=4, fontColor="darkgray", fontSize=24)
xlsx.addLineBreak(sheet, 3)
# header
xlsx.addHeader(wb, sheet, value=paste("Table 1. The candiate genes expressed in ", TCGA_cohort,
                                      " (ranking by KM P-value, selected by Z-score > ", signif(0.6, 2), ") ", "\n" , sep=""),
               level=5, color="black", underline=0) # nrow(HNSCC_OS_marginS_pvalue1e_6_zscore0_6)
#xlsx.addHeader(wb, sheet, value=paste("Cutoff at ", round(cutoff1, 3), " (", percent(surv_OS1$n[1]/(surv_OS1$n[1]+surv_OS1$n[2])), ")", sep = ""),
#               level=5, color="red", underline=0) # total n is taken from surv_OS1$n

xlsx.addLineBreak(sheet, 1) # add one blank line

#xlsx.addTable(wb, sheet, data = t(data.frame(c(paste(geneName, "expression"), "", paste(geneName, "expression"), "", "(Optimised)"))), fontSize=12, startCol=4,
#              fontColor="darkblue", row.names = F, col.names = F) #, colSpan=1, rowSpan=1)

# a candidate genes table
# > colnames(HNSCC_OS_marginS_pvalue1e_6_zscore0_6)
# [1] "gene_id"   "p_value"   "z_score"   "number_01"
colnames(HNSCC_OS_marginS_pvalue1e_6_zscore0_6) <- c("Gene_id",   "P_value",   "z_score_raw",   "Z_score") # Z_score is rescaled as 0-1
# # no more Bonferroni Bonferroni_cutoff; [, c(1,2,4)]
attach(HNSCC_OS_marginS_pvalue1e_6_zscore0_6)
candidates_alpha_pvalue <- HNSCC_OS_marginS_pvalue1e_6_zscore0_6[which(P_value<=alpha_HNSCC), c(1,2,4)]
detach(HNSCC_OS_marginS_pvalue1e_6_zscore0_6)

attach(candidates_alpha_pvalue) # removal of as.factor (P_value)
# in scientific notation: formatC of [, c(2)]
candidates_alpha_pvalue$P_value <- formatC(as.numeric(as.character(P_value)), format = "e", digits = 2)
candidates_alpha_pvalue$Z_score <- signif(Z_score, digits=4)
candidates_alpha_pvalue <- candidates_alpha_pvalue[order(P_value, -Z_score), ] #sorting by order(ascending)
detach(candidates_alpha_pvalue)

xlsx.addTable(wb, sheet, data = candidates_alpha_pvalue, startCol=2,
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
  load(file=file.path(path_ZSWIM2, "HNSCC_OS_marginS_pvalue005_sorted.Rda")) # or HNSCC_OS_marginFree_pvalue005_sorted
  # "swap data" for following code running
  
  attach(HNSCC_OS_marginS_pvalue005_sorted) # number_01 is the Z-score rescaled as 0-1
  plot(p_value, number_01, type="p", ylab="Z-score", xlab="P-value (KM survival analysis)", log="x") # log scale x or y
  #axis(side=3, at=c(1e-7, 1e-6, 1e-5, 0.01, 0.05)) #1=below, 2=left, 3=above and 4=right
  abline(h=0.6, lty=2, col="blue")
  abline(v=alpha_HNSCC, lty=2, col="red") #
  # run a logistic regression model (categorical 0 vs 1)
  #g <- glm(number_01 ~ p_value, family=poisson, data=HNSCC_OS_marginS_pvalue005_sorted)
  # (in this case, generalized linear model with log link)(link = "log"), poisson distribution
  #curve(predict(g, data.frame(p_value = x), type="response"), add=TRUE, col="blue") # draws a curve based on prediction from logistic regression model
  legend("topright", legend=c(paste("Z-score at 0.6"), paste("P-value at", alpha_HNSCC)), lty=2:2, col=c("blue","red"), cex=0.9) # box and font size
  # figure legend: logistic regression, LR, by Generalized linear model, glm
  detach(HNSCC_OS_marginS_pvalue005_sorted)
  #}
  dev.off()
  
  #Append plot to the sheet
  if(is.null(startRow)){
    rows<- getRows(sheet) #list of row object
    startRow=length(rows)+1
  } 
  # append to the file created previously
  # or just add 
  # Cox's Harzard Ratios (univariate) plot [Rplot07_cox_uniHR.tiff]
  # Cox's Harzard Ratios (multivariate) plot [Rplot07_cox_multiHR.tiff]
  addPicture("plot.png", sheet=sheet,  startRow = startRow, startColumn = startCol) 
  xlsx.addLineBreak(sheet, round(width/20)+1)
  res<-file.remove("plot.png") # Z-score vs P-value
} # Define function

# calling
currentRow <- xlsx.addPlot.candidates(wb, sheet) # startRow + a plot
#print(paste("The z-score summary plot:", filenamex,"was exported successfully."))



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
  title <- paste(c("HNSCC survival analysis", "KM P-Value <= ", signif(Bonferroni_cutoff, 3)), sep = "", collapse = "\n")
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




# *generate (bad) prognostic features of those genes on lists in TCGA HNSCC cohort ####
# store KM_sig remrark as a Byte; converted as base10 in list_KM_sigBin
load(file=file.path(path_ZSWIM2, "HNSCC_OS_marginS_pvalueKM_candidate_cox.Rda")) # as candidate_cox (a list), since 2019/07/03; with KM_sig, Remark
# candidate_sample, candidate_cox, n_percent_Bonferroni
load(file=file.path(path_ZSWIM2, "HNSCC_OS_marginS_THREE_pvalue005.Rda") )# as as HNSCC_OS_marginS_THREE_pvalue005, Bonferroni_cutoff

#.. bad guy gene candidate
# => see "remark=1"; "which" or NOT "which", the order is wrong.
# isect_HR2p5 is a list, which comes from line 2909
#x geneid_bad_uni_HR2p5 <- c(isect_HR2p5$`A:B`, isect_HR2p5$`A`) # from venn_HR2p5: A + A:B = group uni_
# R4> isect_HR2p5
# $`uni_Cox HR >=1:multi_Cox HR >=1`
# [1] "CAMK2N1" "USP10"   "PGK1"    "SURF4"   "EFNB2"   "EIF2AK1" "STC2"   
# 
# R4> isect_HR2p5
# $`uni_Cox HR >=1:multi_Cox HR >=1`
# [1] "CAMK2N1" "USP10"   "PGK1"    "SURF4"   "EFNB2"   "EIF2AK1" "STC2"   
geneid_bad_uni_HR2p5 <- c(isect_HR2p5[[3]], isect_HR2p5[[1]])
# whole_genome[which(whole_genome %in% geneid_bad_uni_HR2p5)]
#  [1] "CD1A"    "CENPN"   "COX7A2L" "ERCC6"   "HIGD1B"  "HMGN5"   "NEK6"    "NUF2"   
# [9] "PA2G4P4" "PCTP"    "PLOD2"   "RRAS2"   "SFXN1" 
candidate_bad_uni_HR2p5 <- cbind(data.frame(gene_id=geneid_bad_uni_HR2p5), HNSCC_OS_marginS_THREE_pvalue005[which(HNSCC_OS_marginS_THREE_pvalue005$gene_id %in% geneid_bad_uni_HR2p5), 3:7])
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
candidate_bad_multi_HR2p5 <- cbind(data.frame(gene_id=geneid_bad_multi_HR2p5), HNSCC_OS_marginS_THREE_pvalue005[(HNSCC_OS_marginS_THREE_pvalue005$gene_id %in% geneid_bad_multi_HR2p5), 3:7])
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
candidate_bad_unimulti_HR2p5 <- cbind(data.frame(gene_id=geneid_bad_unimulti_HR2p5), HNSCC_OS_marginS_THREE_pvalue005[(HNSCC_OS_marginS_THREE_pvalue005$gene_id %in% geneid_bad_unimulti_HR2p5), 3:7])


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
  title <- paste(c("HNSCC survival analysis", "KM P-Value <= ", signif(Bonferroni_cutoff, 3)), sep = "", collapse = "\n")
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


# *generate (good) prognostic features of those genes on lists in TCGA HNSCC cohort ####
# => see "remark=1" 
geneid_good_uni_HR0p5 <- c(isect_HR0p5[[3]], isect_HR0p5[[1]]) #$`A:B`
candidate_good_uni_HR0p5 <- cbind(data.frame(gene_id=geneid_good_uni_HR0p5), HNSCC_OS_marginS_THREE_pvalue005[(HNSCC_OS_marginS_THREE_pvalue005$gene_id %in% geneid_good_uni_HR0p5), 3:7])
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
candidate_good_multi_HR0p5 <- cbind(data.frame(gene_id=geneid_good_multi_HR0p5), HNSCC_OS_marginS_THREE_pvalue005[(HNSCC_OS_marginS_THREE_pvalue005$gene_id %in% geneid_good_multi_HR0p5), 3:7])
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
candidate_good_unimulti_HR0p5 <- cbind(data.frame(gene_id=geneid_good_unimulti_HR0p5), HNSCC_OS_marginS_THREE_pvalue005[(HNSCC_OS_marginS_THREE_pvalue005$gene_id %in% geneid_good_unimulti_HR0p5), 3:7])
# *** using table 2 to interpretate the impact on survival by DEG of this candidate gene list 
# c("CRHR2", "MYLIP", "NUP210L", "ZKSCAN4")
# => see "remark=1" on HNSCC_survivalAnalysis_marginS_CRHR2.xlsx => higher CRHR2 expression is associated with less LN metastaese
# => see "remark=1" on HNSCC_survivalAnalysis_marginS_CRHR2.xlsx => higher MYLIP expression is associated with smaller tumor size (T)
# => see "remark=1" on HNSCC_survivalAnalysis_marginS_CRHR2.xlsx => higher NUP210L expression is associated with less LN metastaese
# => see "remark=1" on HNSCC_survivalAnalysis_marginS_CRHR2.xlsx => higher ZKSCAN4 expression is associated with smaller tumor size (T)

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
setwd(path_cohort)
#xlsx.openFile(filenamex) # open file to review

# the END of R2Excel ###
#}



# >Analysis finish; tar until here (refinement ok) ####
# R4> options(prompt="R4_plus>")
# tar
#{
#好用的工具
# ***bash $ TODAY=`date +"%b %d"`;ls -l | grep "$TODAY"
# # list today's files
# https://www.howtogeek.com/248780/how-to-compress-and-extract-files-using-the-tar-command-on-linux/
#   https://www.gnu.org/software/tar/manual/tar.html

# =tar and scp from a list of files .xlsx 
# (download
# generated genes list .xlsx and ZSWIM2_free_archive.Rda; 
# analysis -> .tiff, HNSCC_OS_marginS_candidates_Venn.xlsx, HNSCC_OS_marginS_THREE_pvalue005.Rda)
# $ ls HNSCC_survivalAnalysis_marginS_*.* > marginS_list.txt
# $ tar -czvf ~/R/marginPlus_xlsx.tar.gz -T marginPlus_list.txt   # or  list as many directories
# $ info tar # tar -t --list  —remove-file...
# 
# $ scp  tex@35.201.169.0:~/margin*_xlsx.tar.gz ./
# tex@instance-4:$ ~/R/LUAD_Peter_survival$ sudo mv LUAD_survivalAnalysis_marginS*.* ./survivalAnalysis_marginS/
# tar -xvf marginPlus_xlsx.tar HNSCC_survivalAnalysis_marginPlus_*.Rda
# $ tar -xzvf archive.tar.gz -C /tmp # to extract them to /tmp
# $ tar -xzf archive.tar.gz --overwrite
# $ ls -halt | grep -v ".Rda"
