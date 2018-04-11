#' Mediator LOD
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
#' med_lod <- mediator_lod(mediator = Tmem68$mediator,
#'                       driver = Tmem68$qtl.geno,
#'                       annotation = Tmem68$annotation,
#'                       covar_med = NULL)
#' ggplot2::autoplot(med_lod) +
#'   ggplot2::geom_hline(yintercept = 5, col = "blue")
#' 
#' @export

mediator_lod <- function(mediator, 
                         driver, 
                         annotation, 
                         covar_med=NULL, 
                         verbose=TRUE) {

  # calculates log10-Likelihood of linear model y ~ 1 + X
  LL <- function(y, X) {
    -length(y)/2*log10(sum(qr.resid(qr(cbind(X,1)),y)^2))
  }
  
  # Synch sample IDs.
  tmp = synch.samples(pheno = mediator[,1, drop=FALSE], probs = driver, expr = mediator, covar = covar_med)
  driver   = tmp$probs
  mediator = tmp$expr
  covar_med    = tmp$covar
  rm(tmp)

  # check input
  stopifnot(NROW(annotation) == NCOL(mediator))
  stopifnot(NROW(driver) == NROW(mediator))
  stopifnot(is.null(covar_med) | NROW(mediator) == NROW(covar_med))
  if(!is.null(covar_med)) {
    stopifnot(all(is.numeric(covar_med)))
    stopifnot(!any(is.na(covar_med)))
  }
  stopifnot(!any(is.na(driver)))
  stopifnot(all(is.numeric(mediator)))
  stopifnot(all(is.numeric(driver)))
  stopifnot(c("CHR", "POS") %in% toupper(names(annotation)))

  # data preparation
  mediator <- cbind(mediator) # to ensure 'mediator' is a matrix
  N <- nrow(mediator) # number of individuals
  if (is.null(covar_med)) covar_med <- cbind(rep(1, N)) # if no covariates, use just intercept

  med_pur <- 
    purrr::transpose(
      list(mediator = as.data.frame(mediator),
           annotation = split(annotation, rownames(annotation))))
  lodfn <- function(loglik, LOD0) { loglik[2] - loglik[1] }
  
  mapfn <- function(x, covar_med, driver, LOD0) {
    no.na <- !is.na(x$mediator)
    LL(x$mediator[no.na], cbind(covar_med[no.na,], driver[no.na,])) -
      LL(x$mediator[no.na], covar_med[no.na,])
  }
  output <- annotation
  output$lod <- unlist(purrr::map(med_pur, mapfn, covar_med, driver, LOD0))
  class(output) <- c("mediation_scan", "data.frame")
  return(output)
}
