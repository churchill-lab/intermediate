#' @export
#'
#' @importFrom purrr simplify_all transpose
#' @importFrom stringr str_split
#' @importFrom dplyr mutate
#'
calcZ <- function(models,
                  S.hat = calcShat(models),
                  ICs = calcICs(models, flavor),
                  flavor = c("B","A"),
                  ...) {
  
  n_ind <- nrow(models$indLR)
  
  LRt <- outer(ICs, ICs, function(x,y) y-x)
  LRt <- LRt[lower.tri(LRt)]
  LRt <- -LRt / (2 * sqrt(n_ind))
  LRt / sqrt(diag(S.hat))
}
#' @export
#'
calcICs <- function(models, flavor = c("B","A")) {

  flavor <- match.arg(flavor)
  n_ind <- nrow(models$indLR)
  -2 * models$LR + models$df * penalty(n_ind, flavor)
}
penalty <- function(n_ind, flavor) {
  switch(flavor, A = 2, B = log(n_ind))
}