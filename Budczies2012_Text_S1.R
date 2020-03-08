# Cutoff Finder, version 2.0
# 看看範例 cutoff finding\cite{Budczies2012} 是如何寫的: one run one gene (biomarker)
#http://molpath.charite.de/cutoff
# 10th September, 2012
# Jan Budczies, Charite - Universitaetsmedizin Berlin
# jan.budczies@charite.de
#
# Helper functions:
# get.percent()
# calc.auc()
#
# Wrapper function to be called by tomcat:
# get.cutoff()
# 
# Functions for cutoff optimization:
# cutoff.distribution()
# cutoff.outcome()
# cutoff.survival()
#
# Overview plot functions: 
# plot.distribution()
# plot.OR()
# plot.ROC()
# plot.HR()
# plot.time()
#
# Plot functions at cutoff point:
# plot.waterfall()
# plot.kaplanmeier()

library(flexmix)
library(binom)
library(survival)

NUMBER.NA <- 999999.999

get.percent <- function(k, n, method="wilson") {
  res <- as.numeric(binom.confint(k, n, method=method)[, 4:6])
  return(res)
}

calc.auc <- function (S) 
{
  sen <- S[, "sensitivity"]
  spe <- S[, "specificity"]
  if (max(sen) <= 1) {
    sen <- sen * 100
    spe <- spe * 100
  }
  n <- length(sen)
  x <- 1 - spe/100
  y <- sen/100
  index <- order(x)
  x <- x[index]
  y <- y[index]
  x <- c(0, x, 1)
  y <- c(0, y, 1)
  z <- y - x
  a <- sqrt(2) * x + z/sqrt(2)
  b <- z/sqrt(2)
  index <- order(a)
  a <- a[index]
  b <- b[index]
  auc <- try(integrate(function(x) {approx(a, b, x)$y}, 0, sqrt(2), subdivisions = 1000))
  if (length(auc) == 1) auc <- 0
  else auc <- auc$value + sign(auc$value) * 0.5
  return(auc)
}

get.cutoff <- function(type=c("none", "distribution", "outcome_significance", "outcome_euclidean", "outcome_manhattan", "outcome_sensitivity", "outcome_specificity", "survival_significance", "manual"), filename, biomarker=NULL, outcome=NULL, time=NULL, event=NULL, cutoff=NULL, threshold=NULL, plots=c("histogram", "OR", "ROC", "HR", " time", "waterfall", "kaplanmeier"), nmin=10) {
  type.variable <- unlist(strsplit(type, "_"))[1]
  outcome.method <- "none"
  survival.method <- "none"
  if (type.variable == "outcome") outcome.method <- unlist(strsplit(type, "_"))[2]
  if (type.variable == "survival") survival.method <- unlist(strsplit(type, "_"))[2]
  if (!(type == "manual")) cutoff <- NULL 
  msg <- vector()
  pics <- vector()
  n <- length(biomarker)
  msg[length(msg)+1] <- paste(n, "patient data sets loaded.")
  index <- list()
  index$biomarker <- which(biomarker != NUMBER.NA)
  nbiomarker <- length(index$biomarker)
  msg[length(msg)+1] <- paste(nbiomarker, "data sets with biomarker measurements.") 
  index$outcome <- which(outcome %in% 0:1)
  noutcome <- length(index$outcome)
  msg[length(msg)+1] <- paste(noutcome, "data sets with outcome information.")
  index$survival <- intersect(which(time != NUMBER.NA), which(event %in% 0:1))
  nsurvival <- length(index$survival)
  msg[length(msg)+1] <- paste(nsurvival, "data sets with survival information.")
  if (length(biomarker) > 0) {
    if (is.null(names(biomarker)[1])) names(biomarker)[1] <- "biomarker"
    names(biomarker) <- rep(names(biomarker)[1], n)
  }
  if (length(outcome) > 0) {
    if (is.null(names(outcome)[1])) names(outcome)[1] <- "outcome"
    names(outcome) <- rep(names(outcome)[1], length(outcome))
  }
  if (length(time) > 0) {
    if (is.null(names(time)[1])) names(time)[1] <- "survival"
    names(time) <- rep(names(time)[1], length(time))
  }
  biomarker.name <- names(biomarker)[1]
  dir.name <- dirname(filename)
  data.name <- basename(filename)
  type.number <- type
  if (type %in% c("outcome_sensitivity", "outcome_specificity")) type.number <- paste(type, threshold, sep="")
  out.file <- paste(dir.name, paste("CF", data.name, biomarker.name, type.number, sep="_"), sep="/")
  if (nbiomarker < 2*nmin) msg[length(msg)+1] <- paste("<br>ERROR: Number of patients must be greater than ", 2*nmin, ".", sep="")
  else { 
    if (type.variable == "distribution") {
      ind <- index$biomarker
      res <- try(cutoff <- cutoff.distribution(biomarker[ind]))
      status <- attr(res, "class")
      if (!is.null(status) && status == "try-error") msg[length(msg)+1] <- paste("<br>", sub("\n", "", res), sep="")
    }
    OR <- NULL
    if (type.variable == "outcome" || "OR" %in% plots || "ROC" %in% plots) {
      ind <- intersect(index$biomarker, index$outcome)
      if (length(ind) < 2*nmin) msg[length(msg)+1] <- paste("<br>ERROR: Number of patients with outcome information must be greater than ", 2*nmin, ".", sep="")
      else {  
        res <- try(OR <- cutoff.outcome(biomarker[ind], outcome[ind], method=outcome.method, thres.method=threshold))
        status <- attr(res, "class")
        if (!is.null(status) && status == "try-error") msg[length(msg)+1] <- paste("<br>", sub("\n", "", res), sep="")
        if (type.variable == "outcome") {   
          res <- try(cutoff <- OR["optimal", biomarker.name])
          status <- attr(res, "class")
          if (!is.null(status) && status == "try-error") msg[length(msg)+1] <- paste("<br>", sub("\n", "", res), sep="") 
        }  
      }
    }
    HR <- NULL
    if (type.variable == "survival" || "HR" %in% plots || "time" %in% plots) {
      ind <- intersect(index$biomarker, index$survival)
      if (length(ind) < 2*nmin) msg[length(msg)+1] <- paste("<br>ERROR: Number of patients with survival information must be greater than ", 2*nmin, ".", sep="")
      else {
        res <- try(HR <- cutoff.survival(biomarker[ind], time[ind], event[ind], method=survival.method))
        status <- attr(res, "class")
        if (!is.null(status) && status == "try-error") msg[length(msg)+1] <- paste("<br>", sub("\n", "", res), sep="") 
        if (type.variable == "survival") { 
          res <- try(cutoff <- HR["optimal", biomarker.name])
          status <- attr(res, "class")
          if (!is.null(status) && status == "try-error") msg[length(msg)+1] <- paste("<br>", sub("\n", "", res), sep="")   
        }
      }
    }
    if (!is.null(cutoff)) msg[length(msg)+1] <- paste("<br>Optimal cutoff value using method &quot", type, "&quot: ", signif(cutoff, 4), sep="")
    else {
      if ("type" != "none") msg[length(msg)+1] <- paste("<br>WARNING: Could not determine Optimal cutoff value using method &quot", type, "&quot!", sep="")
    }
    if ("histogram" %in% plots) {
      if (type.variable == "distribution") gauss <- TRUE
      else gauss <- FALSE
      res <- try(plot.histogram(marker=biomarker, cutoff=cutoff, gauss=gauss, jpg.file=out.file))
      status <- attr(res, "class")
      if (!is.null(status) && status == "try-error") msg[length(msg)+1] <- paste("<br>", sub("\n", "", res), sep="")
      else pics <- c(pics, res)
    }
    if ("OR" %in% plots) {
      res <- try(plot.OR(OR, marker=biomarker, outcome=outcome, cutoff=cutoff, jpg.file=out.file))
      status <- attr(res, "class")
      if (!is.null(status) && status == "try-error") msg[length(msg)+1] <- paste("<br>", sub("\n", "", res), sep="")
      else pics <- c(pics, res)
    }   
    if ("ROC" %in% plots) {
      res <- try(plot.ROC(OR, marker=biomarker, outcome=outcome, cutoff=cutoff, jpg.file=out.file))
      status <- attr(res, "class")
      if (!is.null(status) && status == "try-error") msg[length(msg)+1] <- paste("<br>", sub("\n", "", res), sep="")
      else pics <- c(pics, res)
    }
    if ("HR" %in% plots) {
      res <- try(plot.HR(HR, marker=biomarker, time=time, event=event, cutoff=cutoff, jpg.file=out.file))
      status <- attr(res, "class")
      if (!is.null(status) && status == "try-error") msg[length(msg)+1] <- paste("<br>", sub("\n", "", res), sep="")
      else pics <- c(pics, res)
    }
    if ("time" %in% plots) {
      res <- try(plot.time(HR, marker=biomarker, time=time, event=event, cutoff=cutoff, jpg.file=out.file))
      status <- attr(res, "class")
      if (!is.null(status) && status == "try-error") msg[length(msg)+1] <- paste("<br>", sub("\n", "", res), sep="")
      else pics <- c(pics, res)
    }
    if (!is.null(cutoff)) {
      ind <- intersect(index$biomarker, index$outcome)
      if ("waterfall" %in% plots && length(ind) >= 2*nmin) {
        res <- try(plot.waterfall(biomarker[ind], outcome[ind], cutoff=cutoff, jpg.file=out.file))
        status <- attr(res, "class")
        if (!is.null(status) && status == "try-error") msg[length(msg)+1] <- paste("<br>", sub("\n", "", res), sep="")   
        else pics <- c(pics, res)
      }
      ind <- intersect(index$biomarker, index$survival)
      if ("kaplanmeier" %in% plots && length(ind) >= 2*nmin) {
        res <- try(plot.kaplanmeier(biomarker[ind], time[ind], event[ind], cutoff=cutoff, jpg.file=out.file))
        status <- attr(res, "class")
        if (!is.null(status) && status == "try-error") msg[length(msg)+1] <- paste("<br>", sub("\n", "", res), sep="")   
        else pics <- c(pics, res)
      }
    }
  }
  result <- c("", sapply(pics, basename))
  for (i in 1:length(msg)) result[1] <- paste(result[1], msg[i], "<br>\n", sep="")
  names(result) <- NULL
  return(result)
}

cutoff.distribution <- function(marker, npoints=1000, nbreaks=100) {
  cat("Optimizing cutoff using method distribution.")
  npatients <- length(marker)
  index <- which(!is.na(marker))
  marker <- marker[index]
  n <- length(marker)
  nmarker <- length(unique(marker))
  if (nmarker < 3) stop("insufficient data")
  median.marker <- median(marker)
  cluster <- 1 + as.numeric(marker > median.marker)
  fit <- flexmix(marker~1, k=2, cluster=cluster)
  param <- parameters(fit)
  tab <- summary(fit)@comptab
  if (nrow(tab) != 2) stop("Could not fit mixture model.")
  min.x <- min(marker)
  max.x <- max(marker)
  d.x <- max.x - min.x  
  x <- seq(min.x, max.x, d.x/npoints)
  npoints <- length(x)
  y1 <- tab["Comp.1", "prior"] * dnorm(x, mean=param[1, "Comp.1"], sd=param[2, "Comp.1"]) 
  y2 <- tab["Comp.2", "prior"] * dnorm(x, mean=param[1, "Comp.2"], sd=param[2, "Comp.2"]) 
  z <- (y1[1:(npoints-1)] - y2[1:(npoints-1)]) * (y1[2:npoints] - y2[2:npoints])
  index <- which(z <= 0)
  if (length(index) == 0) {
    if (y1[1] > y2[1]) cutoff <- min.x 
    else cutoff <- max.x
  }
  else {
    if (length(index) == 1) ind <- index
    else {
      ind <- index[1]
      for (i in 2:length(index)) if (abs(x[ind] - median.marker) > abs(x[index[i]] - median.marker)) ind <- index[i]
    }
    if (ind == 1) ind <- 2
    cutoff <- (x[ind-1] + x[ind]) / 2
  }
  cat("\n")
  return(cutoff)
}

cutoff.outcome <- function(marker, outcome, method="significance", thres.method=NULL, thres.p=0.05, nmin=10, conf.level=95) {
  cat(paste("Optimizing cutoff using method outcome", method, sep="_"))
  if (is.null(names(marker[1]))) marker.name <- "marker"
  else marker.name <- names(marker)[1]
  if (is.null(names(outcome)[1])) outcome.name <- "outcome"
  else outcome.name <- names(outcome)[1]
  npatients <- length(marker)
  index <- intersect(which(!is.na(marker)), which(!is.na(outcome)))
  marker <- marker[index]
  outcome <- outcome[index]
  n <- length(marker)
  nmarker <- length(unique(marker))
  if (nmarker < 3) stop("insufficient data")
  index <- order(marker)
  marker <- marker[index]
  outcome <- outcome[index]
  y <- outcome 
  q <- 1-(1-conf.level/100)/2
  z <- qnorm(q)
  nlow <- nmin:(n-nmin)
  Y <- matrix(nrow=length(nlow), ncol=24)
  colnames(Y) <- c(marker.name, paste(marker.name, "lower", sep="_"), paste(marker.name, "upper", sep="_"), "OR", "OR_lower", "OR_upper", "p_fisher.test", "accuracy", "accuracy_lower", "accuracy_upper", "sensitivity", "sensitivity_lower", "sensitivity_upper", "specificity", "specificity_lower", "specificity_upper", "PPV", "PPV_lower", "PPV_upper", "NPV", "NPV_lower", "NPV_upper", "euclidean", "manhattan")
  rownames(Y) <- nlow
  for (i in nlow) {
    cat(".")
    j <- i-nmin+1
    if (marker[i] != marker[i+1]) {
      x <- c(rep(0, i), rep(1, n-i))
      model <- summary(glm(y ~ x, family="binomial"))
      coef <- model$coefficients["x", ]
      tab <- table(x, y)
      j <- i-nmin+1   
      Y[j, 1] <- (marker[i] + marker[i+1]) / 2
      Y[j, 2] <- marker[i]
      Y[j, 3] <- marker[i+1]
      Y[j, "OR"] <- exp(coef[1])
      Y[j, "OR_lower"] <- exp(coef[1] - z * coef[2])
      Y[j, "OR_upper"] <- exp(coef[1] + z * coef[2])
      Y[j, "p_fisher.test"] <- fisher.test(tab)$p.value
      y.N <- y[1:i]
      y.P <- y[(i + 1):n]
      n.TP <- length(which(y.P == 1))
      n.FP <- length(which(y.P == 0))
      n.TN <- length(which(y.N == 0))
      n.FN <- length(which(y.N == 1))
      for (measure in c("accuracy", "sensitivity", "specificity", "PPV", "NPV")) {
        if (measure == "accuracy") {a <- n.TP + n.FP; b <- n}
        if (measure == "sensitivity") {a <- n.TP; b <- n.TP + n.FN}
        if (measure == "specificity") {a <- n.TN; b <- n.TN + n.FP}
        if (measure == "PPV") {a <- n.TP; b <-n.TP + n.FP}
        if (measure == "NPV") {a <- n.TN; b <-n.TN + n.FN}
        Y[j, paste(measure, c("", "_lower", "_upper"), sep="")] <- 100 * get.percent(a, b)
      }
    }
    else Y[j, "p_fisher.test"] <- 2
  }
  Y[, "euclidean"] <- sqrt(Y[, "sensitivity"]^2 + Y[, "specificity"]^2)
  Y[, "manhattan"] <- Y[, "sensitivity"] + Y[, "specificity"]
  index.optimal <- NULL
  if (method == "significance") index.optimal <- which.min(Y[, "p_fisher.test"])
  if (method == "euclidean") index.optimal <- which.max(Y[, "euclidean"])
  if (method == "manhattan") index.optimal <- which.max(Y[, "manhattan"])
  index.diff <- which(Y[, "p_fisher.test"] < 2)
  if (method == "sensitivity") {
    index <- intersect(which(Y[, "sensitivity"] >= thres.method), index.diff)
    if (length(index) > 0) index.optimal <- index[which.min(Y[index, "sensitivity"])]
  }
  if (method == "specificity") {
    index <- intersect(which(Y[, "specificity"] >= thres.method), index.diff)
    if (length(index) > 0) index.optimal <- index[which.min(Y[index, "specificity"])]
  }  
  if (!is.null(index.optimal)) rownames(Y)[index.optimal] <- "optimal"
  cat("\n")
  return(Y) 
}

cutoff.survival <- function(marker, time, event, method="significance", thres.p=0.05, nmin=10, surv.test="sctest", endtime=max(time), conf.level=95) {
  cat(paste("Optimizing cutoff using method survival", method, sep="_"))
  marker.name <- names(marker)[1]
  survival.name <- names(time)[1]
  npatients <- length(marker)
  index <- intersect(which(!is.na(marker)), intersect(which(!is.na(time)), which(!is.na(event))))
  marker <- marker[index]
  time <- time[index]
  event <- event[index]
  n <- length(marker)
  nevent <- length(which(event == 1))
  nmarker <- length(unique(marker))
  if (nmarker < 3) stop("insufficient data")
  index <- order(marker)
  marker <- marker[index]
  time <- time[index]
  event <- event[index]
  y <- Surv(time, event) 
  q <- 1-(1-conf.level/100)/2
  z <- qnorm(q)
  nlow <- nmin:(n-nmin)
  Y <- matrix(nrow=length(nlow), ncol=11)
  colnames(Y) <- c(marker.name, paste(marker.name, "lower", sep="_"), paste(marker.name, "upper", sep="_"), "HR", "HR_lower", "HR_upper", "low_mean", "low_sd", "high_mean", "high_sd", "p")
  rownames(Y) <- nlow
  for (i in nlow) {
    cat(".")
    j <- i-nmin+1   
    if (marker[i] != marker[i+1]) {
      x <- c(rep(0, i), rep(1, n-i))
      model <- summary(coxph(y ~ x))
      coef <- model$coefficients
      Y[j, 1] <- (marker[i] + marker[i+1]) / 2
      Y[j, 2] <- marker[i]
      Y[j, 3] <- marker[i+1]
      Y[j, "HR"] <- coef[2]
      Y[j, "HR_lower"] <- exp(coef[1] - z * coef[3])
      Y[j, "HR_upper"] <- exp(coef[1] + z * coef[3])
      Y[j, "p"] <- model[[surv.test]]["pvalue"]
      fit <- summary(survfit(y ~ x), rmean=endtime)
      tab <- fit$table
      index.low <- grep(0, rownames(tab))
      Y[j, "low_mean"] <- tab[index.low, "*rmean"]
      Y[j, "low_sd"] <- tab[index.low, "*se(rmean)"]
      index.high <- grep(1, rownames(tab)) 
      Y[j, "high_mean"] <- tab[index.high, "*rmean"]
      Y[j, "high_sd"] <- tab[index.high, "*se(rmean)"]
    }
    else Y[j, "p"] <- 2
  }
  if (method == "significance") {
    index.optimal <- which.min(as.numeric(Y[, "p"]))
    rownames(Y)[index.optimal] <- "optimal"
  }
  index.p <- which(colnames(Y) == "p")
  colnames(Y)[index.p] <- paste("p", surv.test, sep="_")
  cat("\n")
  return(Y)
}

plot.histogram <- function(marker, cutoff=NULL, gauss=TRUE, npoints=1000, nbreaks=100, jpg.width=1000*sqrt(2), jpg.height=1000, jpg.pointsize=30, lwd=3, jpg.file=NULL) {
  cat("Plot histogram\n")
  if (is.null(names(marker)[1])) names(marker)[1] <- "biomarker"
  marker.name <- names(marker)[1]
  npatients <- length(marker)
  index <- which(!is.na(marker))
  marker <- marker[index]
  n <- length(marker)
  nmarker <- length(unique(marker))
  if (nmarker < 3) stop("insufficient data")
  min.x <- min(marker)
  max.x <- max(marker)
  h <- hist(marker, breaks=nbreaks, plot=FALSE)
  if (gauss) {
    median.marker <- median(marker)
    cluster <- 1 + as.numeric(marker > median.marker)
    fit <- flexmix(marker~1, k=2, cluster=cluster)
    param <- parameters(fit)
    tab <- summary(fit)@comptab
    if (nrow(tab) != 2) stop("Could not fit mixture model.")
    d.x <- max.x - min.x  
    x <- seq(min.x, max.x, d.x/npoints)
    npoints <- length(x)
    y1 <- tab["Comp.1", "prior"] * dnorm(x, mean=param[1, "Comp.1"], sd=param[2, "Comp.1"]) 
    y2 <- tab["Comp.2", "prior"] * dnorm(x, mean=param[1, "Comp.2"], sd=param[2, "Comp.2"]) 
    delta <- h$breaks[2] - h$breaks[1]
    y1 <- delta * n * y1
    y2 <- delta * n * y2
    max.y <- max(c(y1, y2, h$counts))
  }
  else max.y <- max(h$counts)
  if (!is.null(cutoff)) {
    main <- paste("cutoff" , "=", signif(cutoff, 4))
    n.high <- length(which(marker > cutoff))
    percent.high <- round(n.high/n * 100, 1)
    n.low <- length(which(marker < cutoff))
    percent.low <- round(n.low/n * 100, 1)
    symbol <- unlist(strsplit(marker.name, " "))[1]
    main <- paste(main, ", ", n.high, " (", percent.high, "%) ", symbol, "+, ", n.low, " (", percent.low, "%) ", symbol, "-", sep="")
  }
  else main <- ""
  out.file <- jpg.file
  if (!is.null(cutoff)) out.file <- paste(out.file, "_cutoff", signif(cutoff, 4), sep="")
  out.file <- paste(out.file, "histogram.jpg", sep="_")
  if (!is.null(jpg.file)) jpeg(out.file, width=jpg.width, height=jpg.height, pointsize=jpg.pointsize, quality=98)
  par(lwd=2)
  hist(marker, breaks=nbreaks, xlim=c(min.x, max.x), ylim=c(0, max.y+1), col="black", lwd=2, xlab=marker.name, main=main, cex.main=1, font.main=1)
  if (gauss) {
    lines(x, y1, col="red", lwd=lwd)
    lines(x, y2, col="red", lwd=lwd)
  }
  if (!is.null(cutoff)) abline(v=cutoff, col="red", lwd=lwd)
  if (!is.null(jpg.file)) dev.off()
  return(out.file)
}

plot.OR <- function(Y, marker=NULL, outcome=NULL, cutoff=NULL, thres.p=0.05, conf.level=95, jpg.width=1000*sqrt(2), jpg.height=1000, jpg.pointsize=30, lwd=3, jpg.file=NULL) {
  cat("Plot OR\n")
  if (is.null(names(marker)[1])) names(marker)[1] <- "biomarker"
  if (is.null(names(outcome)[1])) names(outcome)[1] <- "outcome"
  marker.name <- names(marker)[1]
  outcome.name <- names(outcome)[1] 
  index <- which(as.numeric(Y[, "p_fisher.test"]) <= 1)
  index.sig <- which(as.numeric(Y[, "p_fisher.test"]) < thres.p)
  n.all <- length(index)
  n.sig <- length(index.sig)  
  percent.sig <- 100 * n.sig / n.all
  main <- paste("Significant (p < ", thres.p, ") tests: ", n.sig, " out of ", n.all, " (", round(percent.sig, 1), "%)", sep="")
  Z <- Y[index, ]
  Z[, c("OR", "OR_lower", "OR_upper")] <- log2(Z[, c("OR", "OR_lower", "OR_upper")])
  if ("optimal" %in% rownames(Z)) zopt <- round(Z["optimal", "OR"])
  else zopt <- median(Z[, "OR"])
  zmin <- floor(min(Z[, "OR"])) - 1
  if (zmin < zopt-3) zmin <- zopt-3
  if (zmin > -1) zmin <- -2
  zmax <- ceiling(max(Z[, "OR"])) + 1
  if (zmax > zopt+3) zmax <- zopt+3
  if (zmax < 1) zmax <- 2
  out.file <- jpg.file
  if (!is.null(cutoff)) out.file <- paste(out.file, "_cutoff", signif(cutoff, 4), sep="")
  out.file <- paste(out.file, "OR.jpg", sep="_")
  if (!is.null(jpg.file)) jpeg(out.file, width=jpg.width, height=jpg.height, pointsize=jpg.pointsize, quality=98)
  par(lwd=lwd)
  plot(0, 0, type="n", xlim=range(Z[, 1]), ylim=c(zmin - 0.4, zmax + 0.2), ylab=paste("OR with ", conf.level, "% CI", sep=""), xlab=marker.name, main=main, font.main=1, cex.main=1, yaxt="n")
  at <-  zmin:zmax
  at.pos <- 1:zmax
  at.neg <- abs(zmin):1
  labels <- c(paste(1, 2^at.neg, sep="/"), 1, 2^at.pos)
  axis(2, at=at, labels=labels)
  lines(Z[, marker.name], Z[, "OR"])
  lines(Z[, marker.name], Z[, "OR_lower"], lty=3)
  lines(Z[, marker.name], Z[, "OR_upper"], lty=3)
  points(marker, rep(zmin - 0.4, length(marker)), pch="|")
  if (!is.null(cutoff)) abline(v=cutoff, lty=1)
  abline(h=0, lty=2)
  if (!is.null(jpg.file)) dev.off()
  return(out.file)
}

plot.ROC <- function(Y, marker=NULL, outcome=NULL, cutoff=NULL, jpg.width=1000*sqrt(2), jpg.pointsize=30, lwd=4, jpg.file=NULL) {
  cat("Plot ROC\n")
  if (is.null(names(marker)[1])) names(marker)[1] <- "biomarker"
  if (is.null(names(outcome)[1])) names(outcome)[1] <- "outcome"
  marker.name <- names(marker)[1]
  outcome.name <- names(outcome)[1]
  index <- which(as.numeric(Y[, "p_fisher.test"]) <= 1) 
  Z <- Y[index, ]
  auc <- calc.auc(Z)
  if (auc < 0.985) auc <- round(auc, 2)
  else auc <- 1 - signif(1 - auc, 1)
  sen <- Z[, "sensitivity"]
  spe <- Z[, "specificity"]
  direction <- "positive"
  if (auc < 0) {
    direction <- "negative"
    foo <- sen
    sen <- 100 - spe
    spe <- 100 - foo
    auc <- -auc
  }
  main <- paste(marker.name, "as", direction, "marker", "for", outcome.name)
  legend.auc <- paste("AUC =", auc)
  out.file <- jpg.file 
  if (!is.null(cutoff)) out.file <- paste(out.file, "_cutoff", signif(cutoff, 4), sep="")
  out.file <- paste(out.file, "ROC.jpg", sep="_")
  if (!is.null(jpg.file)) jpeg(out.file, width=jpg.width, height=jpg.width, pointsize=jpg.pointsize*sqrt(2), quality=98)
  par(lwd=lwd)
  plot(c(100, 100 - spe, 0), c(100, sen, 0), type="l", lty=1, lwd=lwd, xlim=c(0, 100), ylim=c(0, 100), xlab="1 - Specificity (%)", ylab="Sensitivity (%)", main=main, cex.main=1, font.main=1)
  abline(0, 1, lty=2, lwd=lwd)
  if (is.null(cutoff)) index.optimal <- NULL  
  else {
    index.lower <- which(Z[, paste(marker.name, "lower", sep="_")] < cutoff)
    index.upper <- which(Z[, paste(marker.name, "upper", sep="_")] > cutoff)
    index.optimal <- intersect(index.lower, index.upper)
  }
  if (length(index.optimal) == 1) {
    legend.cut <- paste("Cutoff = ", signif(cutoff, 4), sep="")
    legend.sen <- paste("Sensitivity = ", round(sen[index.optimal], 1), "%", sep="")
    legend.spe <- paste("Specificity = ", round(spe[index.optimal], 1), "%", sep="")
    points(100 - spe[index.optimal], sen[index.optimal], pch=4, col="red", cex=1.5, lwd=lwd)
    legend("bottomright", legend=c(legend.auc, legend.cut, legend.sen, legend.spe), lty=c(1, -1, -1, -1), pch=c(-1, 4, -1, -1), col=c("black", "red", "white", "white"), inset=0.02, cex=0.9, lwd=lwd)
  }
  else {
    legend("bottomright", legend=legend.auc, inset=0.02, cex=0.9, lwd=lwd)
  }
  if (!is.null(jpg.file)) dev.off()
  return(out.file)
}

plot.HR <- function(Y, marker, time, event, cutoff=NULL, thres.p=0.05, conf.level=95, jpg.width=1000*sqrt(2), jpg.height=1000, jpg.pointsize=30, lwd=3, jpg.file=NULL) {
  cat("Plot HR\n")
  if (is.null(names(marker)[1])) names(marker)[1] <- "biomarker"
  if (is.null(names(time)[1])) names(marker)[1] <- "survival"
  marker.name <- names(marker)[1]
  survival.name <- names(time)[1]
  index.p <- grep("p_", colnames(Y))
  index <- which(as.numeric(Y[, index.p]) <= 1)
  index.sig <- which(as.numeric(Y[, index.p]) < thres.p) 
  n.all <- length(index)
  n.sig <- length(index.sig)  
  percent.sig <- 100 * n.sig / n.all
  main <- paste("Significant (p < ", thres.p, ") tests: ", n.sig, " out of ", n.all, " (", round(percent.sig, 1), "%)", sep="")
  Z <- Y[index, ]
  Z[, c("HR", "HR_lower", "HR_upper")] <- log2(Z[, c("HR", "HR_lower", "HR_upper")])
  if ("optimal" %in% rownames(Z)) zopt <- round(Z["optimal", "HR"])
  else zopt <- median(Z[, "HR"])
  zmin <- floor(min(Z[, "HR"])) - 1
  if (zmin < zopt-3) zmin <- zopt-3
  if (zmin > -1) zmin <- -2
  zmax <- ceiling(max(Z[, "HR"])) + 1
  if (zmax > zopt+3) zmax <- zopt+3
  if (zmax < 1) zmax <- 2
  out.file <- jpg.file 
  if (!is.null(cutoff)) out.file <- paste(out.file, "_cutoff", signif(cutoff, 4), sep="")
  out.file <- paste(out.file, "HR.jpg", sep="_")
  if (!is.null(jpg.file)) jpeg(out.file, width=jpg.width, height=jpg.height, pointsize=jpg.pointsize, quality=98)
  par(lwd=lwd)
  plot(0, 0, type="n", xlim=range(Z[, 1]), ylim=c(zmin - 0.4, zmax + 0.2), ylab=paste("HR with ", conf.level, "% CI", sep=""), xlab=marker.name, main=main, font.main=1, cex.main=1, yaxt="n")
  at <- zmin:zmax
  at.pos <- 1:zmax
  at.neg <- abs(zmin):1
  labels <- c(paste(1, 2^at.neg, sep="/"), 1, 2^at.pos)
  axis(2, at=at, labels=labels)
  lines(Z[, marker.name], Z[, "HR"])
  lines(Z[, marker.name], Z[, "HR_lower"], lty=3)
  lines(Z[, marker.name], Z[, "HR_upper"], lty=3)
  points(marker, rep(zmin - 0.4, length(marker)), pch="|")
  if (!is.null(cutoff)) abline(v=cutoff, lty=1)
  abline(h=0, lty=2)
  if (!is.null(jpg.file)) dev.off()
  return(out.file)
}

plot.waterfall <- function(marker, outcome, cutoff=NULL, conf.level=95, jpg.width=1000*sqrt(2), jpg.height=1000, jpg.pointsize=30, lwd=4, jpg.file=NULL) {
  if (!is.null(cutoff)) {
    cat("Plot waterfall\n")
    marker.name <- names(marker)[1]
    outcome.name <- names(outcome)[1]
    npatients <- length(marker)
    index <- intersect(which(!is.na(marker)), which(!is.na(outcome)))
    marker <- marker[index]
    outcome <- outcome[index]
    n <- length(marker)
    if (n < 3) stop("insufficient data")
    if (is.null(cutoff)) cutoff <- median(marker)
    index <- order(marker)
    marker <- marker[index]
    outcome <- outcome[index]
    q <- 1-(1-conf.level/100)/2
    z <- qnorm(q)
    x <- marker
    y <- outcome
    model <- summary(glm(y ~ x, family="binomial"))
    coef <- model$coefficients["x", ]
    OR <- exp(coef[1])
    OR.lower <- exp(coef[1] - z * coef[2])
    OR.upper <- exp(coef[1] + z * coef[2])
    y <- marker - cutoff
    predicted <- (sign(y) + 1) / 2
    index.ok <- which(y * (1/2 - outcome) < 0)
    color <- rep("red", n)
    color[index.ok] <- "green"
    tab <- matrix(nrow=2, ncol=2)
    tab[1, 1] <- length(intersect(which(predicted == 0), which(outcome == 0)))
    tab[1, 2] <- length(intersect(which(predicted == 0), which(outcome == 1)))
    tab[2, 1] <- length(intersect(which(predicted == 1), which(outcome == 0)))
    tab[2, 2] <- length(intersect(which(predicted == 1), which(outcome == 1)))
    p <- fisher.test(tab)$p.value
    sen <- round(100*get.percent(tab[2, 2], tab[1, 2] + tab[2, 2]), 1)
    spe <- round(100*get.percent(tab[1, 1], tab[1, 1] + tab[2, 1]), 1)
    if (sen[1] + spe[1] < 100)  {
      sen <- 100 - sen
      spe <- 100 - spe
      index.red <- which(color == "red")
      index.green <- which(color == "green")
      color <- rep(NA, length(color))
      color[index.red] <- "green"
      color[index.green] <- "red"
    }
    lab.cut <- paste("Cutoff =", signif(cutoff, 4))
    lab.OR <- paste("OR = ", round(OR, 2), " (", round(OR.lower, 2), "-", round(OR.upper, 2), "), p = ", signif(p, 2), sep="")
    lab.sen <- paste("Sensitivity = ", sen[1], "% (", sen[2], "%-", sen[3], "%)", sep="") 
    lab.spe <- paste("Specificity = ", spe[1], "% (", spe[2], "%-", spe[3], "%)", sep="")
    lab <- c(lab.cut, lab.OR, lab.sen, lab.spe)
    out.file <- paste(jpg.file, "_cutoff", signif(cutoff, 4), sep="")
    out.file <- paste(out.file, "waterfall.jpg", sep="_")
    if (!is.null(jpg.file)) jpeg(out.file, width=jpg.width, height=jpg.height, pointsize=jpg.pointsize, quality=98)
    barplot(as.numeric(y), space=0, col=color, xlab="patients", ylab=marker.name, border=NA, axes=FALSE)
    axis(side=1)
    at <- axTicks(side=2)
    delta <- round(cutoff, 1) - cutoff
    axis(at + delta, label=at + round(cutoff, 1), side=2)
    legend("topleft", legend=paste(c("correct", "wrong"), outcome.name, "classification"), lwd=2, col=c("green", "red"), inset=0.02)
    legend("bottomright", legend=lab, bty="n")
    if (!is.null(jpg.file)) dev.off()
  }
  else out.file <- NULL
  return(out.file)
}

plot.kaplanmeier <- function(marker, time, event, cutoff=NULL, surv.test="sctest", conf.level=95, colors=c("black", "red"), jpg.width=1000*sqrt(2), jpg.height=1000, jpg.pointsize=30, lwd=3, jpg.file=NULL) {
  cat("Plot kaplanmeier\n")
  if (!is.null(cutoff)) {
    marker.name <- names(marker)[1]
    survival.name <- names(time)[1]
    npatients <- length(marker)
    index <- intersect(which(!is.na(marker)), intersect(which(!is.na(time)), which(!is.na(event))))
    marker <- marker[index]
    time <- time[index]
    event <- event[index]
    n <- length(marker)
    nevent <- length(which(event == 1))
    if (n < 3) stop("insufficient data")
    if (is.null(cutoff)) cutoff <- median(marker)
    y <- Surv(time, event)
    x <- rep(NA, length(marker))
    x[which(marker < cutoff)] <- 0
    x[which(marker > cutoff)] <- 1
    fit <- survfit(y ~ x)
    model <- summary(coxph(y ~ x))
    coef <- model$coefficients
    HR <- coef[2]
    q <- 1-(1-conf.level/100)/2
    z <- qnorm(q)
    HR.lower <- exp(coef[1] - z * coef[3])
    HR.upper <- exp(coef[1] + z * coef[3])
    p <- model[[surv.test]]["pvalue"]
    out.file <- paste(jpg.file, "_cutoff", signif(cutoff, 4), sep="")
    out.file <- paste(out.file, "kaplanmeier.jpg", sep="_")
    if (!is.null(jpg.file)) jpeg(out.file, width=jpg.width, height=jpg.height, pointsize=jpg.pointsize, quality=98)
    par(lwd=lwd)
    plot(fit, col=colors, lwd=2, xlab=survival.name, ylab="survival proportion")
    legend("bottomleft", legend=paste(marker.name, c("<", ">"), signif(cutoff, 4)), lty=1, col=colors, inset=0.02)
    legend("bottomright", legend=paste("HR = ", round(HR, 2), " (", round(HR.lower, 2), "-", round(HR.upper, 2), "), p = ", signif(p, 2), sep=""), bty="n")
    if (!is.null(jpg.file)) dev.off()
  }
  else out.file <- NULL
  return(out.file)
}

plot.time <- function(Y, marker, time, event, cutoff=NULL, conf.level=95, jpg.width=1000*sqrt(2), jpg.height=1000, jpg.pointsize=30, lwd=3, jpg.file=NULL) {
  cat("Plot difference in mean survival time\n")
  marker.name <- names(marker)[1]
  survival.name <- names(time)[1]
  index.p <- grep("p_", colnames(Y))
  index <- which(as.numeric(Y[, index.p]) <= 1)
  Z <- Y[index, ]
  x <- Z[,  marker.name]
  y <- Z[, "high_mean"] - Z[, "low_mean"]
  y.sd <- sqrt(Z[, "low_sd"]^2 + Z[, "high_sd"]^2)
  q <- 1-(1-conf.level/100)/2
  z <- qnorm(q)
  y.lower <- y - z*y.sd
  y.upper <- y + z*y.sd
  ymin <- min(y.lower)
  ymax <- max(y.upper)
  dy <- max(y) - min(y)
  if (ymax > max(y) + 2*dy) ymax <- quantile(y.upper, probs=0.9)
  if (ymax > max(y) + 2*dy) ymax <- quantile(y.upper, probs=0.8)
  if (ymin < min(y) - 2*dy) ymin <- quantile(y.lower, probs=0.1)
  if (ymin < min(y) - 2*dy) ymin <- quantile(y.lower, probs=0.2)
  ddy <- ymax - ymin
  ymax <- ymax + 3 * ddy/100
  ymin <- ymin - 6 * ddy/100
  out.file <- jpg.file 
  if (!is.null(cutoff)) out.file <- paste(out.file, "_cutoff", signif(cutoff, 4), sep="")
  out.file <- paste(out.file, "time.jpg", sep="_")
  if (!is.null(jpg.file)) jpeg(out.file, width=jpg.width, height=jpg.height, pointsize=jpg.pointsize, quality=98)
  par(lwd=lwd)
  plot(0, 0, type="n", xlim=range(x), ylim=c(ymin, ymax), ylab=paste("Difference in mean ", survival.name, " with ", conf.level, "% CI", sep=""), xlab=marker.name, main="", font.main=1, cex.main=1)
  lines(x, y)
  lines(x, y.lower, lty=3)
  lines(x, y.upper, lty=3)
  points(marker, rep(ymin, length(marker)), pch="|")
  if (!is.null(cutoff)) abline(v=cutoff, lty=1)
  abline(h=0, lty=2)
  if (!is.null(jpg.file)) dev.off()
  return(out.file)
}
