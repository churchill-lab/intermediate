#' Mediation Scan
#'
#' For a given QTL haplotype probabilities \code{driver} and target \code{target},
#' the function sequentially tries to add each column \code{m} of \code{mediator} matrix as a covariate
#' and calculates LOD statistic. The low LOD value indicates \code{driver} and
#' \code{target} are conditionally independent given \code{m},
#' i.e. \code{m} is a mediator of causal relationship from \code{driver} to \code{target}.
#'
#'
#' @param target A numeric vector with gene/protein expression
#' @param mediator A matrix, each column is one gene/protein's expression
#' @param driver A matrix, haplotype probabilities at QTL we try to mediate
#' @param annotation A data frame with mediators' annotation, must include columns "chr" and "pos"
#' @param covar A matrix with additive covariates
#' @param method A method to handle missing cases
#' @param verbose If TRUE display information about the progress
#' 
#' @seealso \code{\link{plot.mediation}}, \code{\link{kplot}}
#' 
#' @examples
#' data(Tmem68)
#' # Find and remove Tmem68 from mediators because it is target.
#' m <- match("Tmem68", Tmem68$annotation$symbol)
#' Tmem68$annotation[m,]
#' med_scan <- mediation_scan(target = Tmem68$target,
#'                       mediator = Tmem68$mediator[,-m],
#'                       driver = Tmem68$qtl.geno,
#'                       annotation = Tmem68$annotation[-m,],
#'                       covar = Tmem68$covar,
#'                       method = "double-lod-diff")
#' ggplot2::autoplot(med_scan)
#' 
#' @export

mediation_scan <- function(target, 
                           mediator, 
                           driver, 
                           annotation, 
                           covar=NULL,
                           method=c("double-lod-diff", "ignore", "lod-diff", 
                                    "lod-ratio"), 
                           verbose=TRUE) {

  # calculates log10-Likelihood of linear model y ~ 1 + X
  LL <- function(y, X) {
    -length(y)/2*log10(sum(qr.resid(qr(cbind(X,1)),y)^2))
  }
  
  # Synch sample IDs.
  tmp = synch.samples(pheno = as.matrix(target),
                      probs = driver, expr = mediator, covar = covar)
  target   = tmp$pheno
  driver   = tmp$probs
  mediator = tmp$expr
  covar    = tmp$covar
  rm(tmp)

  # check input
  stopifnot(NROW(target) == NROW(mediator))
  stopifnot(NROW(annotation) == NCOL(mediator))
  stopifnot(NROW(driver) == NROW(target))
  stopifnot(is.null(covar) | NROW(target) == NROW(covar))
  stopifnot(!any(is.na(covar)))
  stopifnot(!any(is.na(driver)))
  stopifnot(all(is.numeric(target)))
  stopifnot(all(is.numeric(mediator)))
  stopifnot(all(is.numeric(driver)))
  stopifnot(all(is.numeric(covar)))
  stopifnot(c("chr", "pos") %in% tolower(names(annotation)))
  method = match.arg(method)

  # data preparation
  mediator <- cbind(mediator) # to ensure 'mediator' is a matrix
  N <- nrow(mediator) # number of individuals
  if (is.null(covar)) covar <- cbind(rep(1, N)) # if no covariates, use just intercept

  if (method == "double-lod-diff") {
    no.na <- !is.na(target)
    LOD0 <- LL(target[no.na], cbind(covar, driver)[no.na,]) - LL(target[no.na], covar[no.na,])
  }

  med_pur <- 
    purrr::transpose(
      list(mediator = as.data.frame(mediator),
           annotation = split(annotation, rownames(annotation))))
  lodfn <- 
    switch(
      method,
      ignore = function(loglik, LOD0) {
        loglik[2] - loglik[1]
      },
      "lod-diff" = function(loglik, LOD0) {
        loglik[4] - loglik[3] - (loglik[2] - loglik[1])
      },
      "double-lod-diff" = function(loglik, LOD0) {
        LOD0 - (loglik[4] - loglik[3] - (loglik[2] - loglik[1]))
      },
      "lod-ratio" = function(loglik, LOD0) {
        (10^loglik[2]-10^loglik[1]) / (10^loglik[4] - 10^loglik[3])
      })
  
  loglik0 <- LL(target[no.na], cbind(covar[no.na,], driver[no.na,])) -
    LL(target[no.na], covar[no.na,])
  
  
  mapfn <- function(x, target, covar, driver, LOD0) {
    no.na <- !is.na(target) & !is.na(x$mediator)
    loglik <- c(LL(target[no.na], cbind(covar[no.na,], x$mediator[no.na])),
                LL(target[no.na], cbind(covar[no.na,], x$mediator[no.na], 
                                        driver[no.na,])),
                LL(target[no.na], covar[no.na,]),
                LL(target[no.na], cbind(covar[no.na,], driver[no.na,])))
    lodfn(loglik, LOD0)
  }
  output <- annotation
  output$lod <- unlist(purrr::map(med_pur, mapfn, target, covar, driver, LOD0))
  attr(output, "targetFit") <- loglik0
  class(output) <- c("mediation_scan", "data.frame")
  return(output)
}
