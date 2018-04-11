#' Fit target relative to driver
#'
#' Fit a model for a target and get detailed results
#' about estimated coefficients and individuals contributions to the LOD score.
#' 
#' @details If \code{kinship} is absent, regression is performed.
#' If \code{kinship} is provided, a linear mixed model is used, with a
#' random effect estimated under the null hypothesis of no driver,
#' and then taken as fixed and known with driver.
#' The default version ignores kinship. See \code{\link[qtl2]{fit1}}
#' for use of \code{kinship}.
#'
#' @param driver A matrix of drivers, individuals x drivers
#' @param target A numeric vector of target values
#' @param kinship Optional kinship matrix.
#' @param addcovar An optional matrix of additive covariates.
#' @param nullcovar An optional matrix of additional additive
#' covariates that are used under the null hypothesis (of no driver)
#' but not under the alternative (with driver).
#' @param intcovar An optional matrix of interactive covariates.
#' @param weights An optional vector of positive weights for the
#' individuals. As with the other inputs, it must have `names`
#' for individual identifiers. Ignored if `kinship` is provided.#' 
#' 
#' @return A list containing
#' * `LR` - The overall likelihood ratio.
#' * `indLR` - Vector of individual contributions to the likelihood ratio.
#' * `df` - Model degrees of freedom.
#' 
#' @export
fitDefault <- function(driver,
                  target,
                  kinship = NULL,
                  addcovar = NULL,
                  intcovar=NULL, weights=NULL,
                  ...) {

  # Construct design matrix X
  if(is.null(driver))
    driver <- matrix(1, length(target), 1)
  
  # Form model matrix from additive covariates.
  if(!is.null(addcovar)) {
    form <- as.formula(paste(" ~ ", paste(colnames(addcovar), collapse = "+")))
    X <- model.matrix(form, data = addcovar)[,-1]
    if(is.null(dim(X))) {
      X <- as.matrix(X)
    }
  } else 
    X <- NULL
  
  # Add interactive covaraties.
  if(!is.null(intcovar)) {
    if(ncol(driver) > 1) {
      int.sub.matrix <- 
        model.matrix(as.formula(paste("~", colnames(intcovar))), intcovar)[,-1]
      driverbyintcovar <- driver[,-1] * int.sub.matrix
      X <- cbind(driver, X, driverbyintcovar)
    } else {
      X <- cbind(driver, X)
    }
  } else {
    if(!is.null(addcovar)){
      X <- cbind(driver, X)
    } else {
      X <- driver
    }
  }
  
  # Calculate log likelihood components
  n <- length(target)
  dX <- ncol(X)
  qrX <- qr(X)
  b <- qr.coef(qrX, target)
  RSS <- target - X %*% b
  RSS <- crossprod(RSS, RSS)

  list(LR = as.vector(- (n/2) - (n/2) * log(2 * pi) - (n/2) * log(RSS/n)),
       indLR = dnorm(target, X %*% b, sqrt(RSS / n), log = TRUE),
       df = dX,
       RSS = RSS)
}

