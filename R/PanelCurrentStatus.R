#' @import stats
NULL

#' Calculate the Conditional Censoring Logistic (CCL) Estimator
#'
#' @description
#' This function computes the scaled coefficients from the conditional censoring logistic 
#' (CCL) estimator for panel current status data. The CCL estimator transforms panel current 
#' status data into a binary outcome analysis, building on existing logistic regression 
#' estimators by incorporating monitoring time information into the working model.
#'
#' @param delta a matrix of event indicators (0/1) where rows correspond to subjects and 
#'        columns correspond to different monitoring times
#' @param ctime a matrix of censoring/monitoring times with the same dimensions as \code{delta}
#' @param predictors a matrix of predictor variables/covariates where rows correspond to subjects
#' @param n.ptb number of bootstrap samples (perturbations) for standard error estimation
#' @param seed an integer specifying the random seed for reproducibility of resampling
#'
#' @export
ccl.fit <- function(delta, ctime, predictors, n.ptb, seed = 1) {
  ccl_fit <- function(delta, ctime, predictors, weights = NULL) {
    fit <- suppressWarnings(
      glm(delta ~ log(ctime) + as.matrix(predictors), 
          weights = weights, family = binomial))
    tmp <- fit$coef[-(1:2)]
    tmp / sqrt(sum(tmp^2))
  }  
  n <- nrow(delta)
  K <- ncol(delta)
  p <- ncol(predictors)
  L <- vector("list", 5)
  for (k in 1:K) {
    L[[k]]$est <- ccl_fit(delta[, k], ctime[, k], predictors)
    L[[k]]$est.ptb <- matrix(NA, p, n.ptb)
    set.seed(seed) 
    W <- matrix(NA, n, n.ptb)
    for (b in 1:n.ptb) {
      W[, b] <- tabulate(sample(1:n, n, replace = TRUE), nbins = n)
      L[[k]]$est.ptb[, b] <- ccl_fit(delta[, k], ctime[, k], predictors, W[, b])
    }
  }
  # optimal combination
  I <- matrix(1, K, 1)
  opt <- rep(NA, p)
  opt.ptb <- matrix(NA, p, n.ptb)
  for (i in 1:p) {
    v <- rep(NA, K)
    M <- matrix(NA, n.ptb, K)
    for (k in 1:K) {
      v[k] <- L[[k]]$est[i]      
      M[, k] <- L[[k]]$est.ptb[i, ]
    }
    S <- cov(M)
    S.inv <- solve(S)
    w <- (S.inv %*% I) / c(t(I) %*% S.inv %*% I)
    w1 <- matrix(w, ncol = 1)
    opt[i] <- v %*% w1
    opt.ptb[i, ] <- M %*% w1 
  }
  est <- opt / sqrt(sum(opt^2))
  est.ptb <- apply(opt.ptb, 2, function(x) x / sqrt(sum(x^2)))
  list(K = K,
       est = est, 
       est.se = apply(est.ptb, 1, sd), 
       est.ptb = est.ptb, 
       W = W)
}

#' Evaluate Prediction Model Performance with ROC Curves
#'
#' @description
#' This function evaluates prediction performance metrics from ROC curve analysis via 
#' kernel smoothing for panel current status data. It calculates area under the curve (AUC) 
#' and true positive rates at specified false positive rates.
#'
#' @param data a data frame in 'long' format containing columns for:
#'        \itemize{
#'          \item{delta: binary event indicator}
#'          \item{ctime: monitoring/examination time}
#'          \item{predictors: covariate values}
#'        }
#' @param fit results from the \code{ccl.fit} function, providing estimated coefficients and bootstrap weights
#' @param t0 pre-specified prediction time point of interest
#' @param h bandwidth for kernel smoothing, controlling the weight assigned to observations 
#'        based on their proximity to \code{t0}
#' 
#' @export
ccl.roc <- function(data, fit, t0, h) {
  sum.I <- function(yy, FUN, Yi, Vi = NULL) {
    if (FUN == "<" | FUN == ">=") { 
      yy <- -yy
      Yi <- -Yi
    }
    pos <- rank(c(yy, Yi), ties.method = 'f')[1:length(yy)] - rank(yy, ties.method = 'f')
    if (substring(FUN, 2, 2) == "=") pos <- length(Yi) - pos
    if (!is.null(Vi)) {
      if (substring(FUN, 2, 2) == "=") tmpind <- order(-Yi) else tmpind <- order(Yi)
      Vi <- apply(as.matrix(Vi)[tmpind, , drop = F], 2, cumsum)
      return(rbind(0, Vi)[pos + 1, ])
    } else return(pos)
  }
  ccl_roc <- function(data, est, t0, h, weights = NULL) {
    if (is.null(weights)) weights <- rep(1, nrow(data))
    k.wts <- dnorm((data$ctime - t0) / h) / h
    si.tilde <- as.matrix(data[, -(1:2)]) %*% est
    tpr.wts <- data$delta * k.wts
    fpr.wts <- (1 - data$delta) * k.wts
    tpr.tilde <- sum.I(si.tilde, "<=", si.tilde, tpr.wts * weights) / 
      sum.I(-100, "<=", si.tilde, tpr.wts * weights)
    fpr.tilde <- sum.I(si.tilde, "<=", si.tilde, fpr.wts * weights) / 
      sum.I(-100, "<=", si.tilde, fpr.wts * weights)
    tpr.prob <- sum.I(-100, "<=", si.tilde, tpr.wts * weights) / sum(k.wts * weights)
    fpr.prob <- sum.I(-100, "<=", si.tilde, fpr.wts * weights) / sum(k.wts * weights)
    ppv.tilde <- (tpr.tilde * tpr.prob) / 
      (tpr.tilde * tpr.prob + fpr.tilde * fpr.prob)
    npv.tilde <- ((1 - fpr.tilde) * fpr.prob) / 
      ((1 - tpr.tilde) * tpr.prob + (1 - fpr.tilde) * fpr.prob)
    tpr.05.tilde <- approx(x = fpr.tilde, y = tpr.tilde, xout = .05)$y
    tpr.10.tilde <- approx(x = fpr.tilde, y = tpr.tilde, xout = .10)$y
    tpr.tilde.sort <- tpr.tilde[order(fpr.tilde, tpr.tilde)]
    fpr.tilde.sort <- fpr.tilde[order(fpr.tilde, tpr.tilde)]
    fpr.tilde.diff <- diff(fpr.tilde.sort)
    auc.tilde.lower <- sum(fpr.tilde.diff * tpr.tilde.sort[1:(nrow(data) - 1)])
    auc.tilde.upper <- sum(fpr.tilde.diff * tpr.tilde.sort[2:nrow(data)])
    c("auc.tilde.lower" = auc.tilde.lower, 
      "auc.tilde.upper" = auc.tilde.upper, 
      "tpr.05.tilde" = tpr.05.tilde, 
      "tpr.10.tilde" = tpr.10.tilde)  
  }  
  n.ptb <- ncol(fit$W)
  roc <- ccl_roc(data, fit$est, t0, h)
  roc.ptb <- matrix(NA, length(roc), n.ptb)
  for (b in 1:n.ptb) {
    weights <- rep(fit$W[, b], fit$K)
    roc.ptb[, b] <- ccl_roc(data, fit$est.ptb[, b], t0, h, weights)
  }
  list(roc = roc, 
       roc.se = apply(roc.ptb, 1, sd),
       roc.ptb = roc.ptb)
}
