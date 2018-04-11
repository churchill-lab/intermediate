#' @export
fitQtl2 <- function(driver,
                    target,
                    kinship,
                    addcovar) {
  out <- qtl2::fit1(driver,
             target,
             kinship,
             addcovar)
  
  # Replace lod names with LR
  names(out) <- stringr::str_replace(names(out), "lod", "LR")
  names(out) <- stringr::str_replace(names(out), "_LR", "LR")

  # Rescael to make them likelihoods (or likelihood ratios)
  out$LR <- out$LR * log(10)
  out$indLR <- out$indLR * log(10)
  
  # Add df for later use
  out$df <- ncol(driver) - 1
  
  out
}