#' @export
#'
binomIUCMST <- function(models, flavor = c("B","A")) {

  # Penalize individual log likelihood ratios
  d <- dim(models$indLR)
  if(d[2] <= 1)
    return(data.frame(ref = names(models$LR),
                     alt = "",
                     pv = 1))
  
  flavor <- match.arg(flavor)
  pen <- penalty(d[1], flavor) / (2 * d[1])
  indLR <- as.data.frame(t(t(models$indLR) - models$df * pen))

  # Outer differences of models (columns) of indLR
  tmpfn <- function(d) {
    c(pos = sum(d > 0),
      nz = sum(d != 0))
  }
  LR <- dplyr::bind_cols(
    purrr::map(indLR,
               function(x,y) x - y,
               indLR))

  names(LR) <- paste(rep(names(models$indLR), each = ncol(models$indLR)),
                     rep(names(models$indLR), ncol(models$indLR)),
                     sep = ":")
  # Reduce to unique model comparisons Gj/Gk with j<k.
  LR <- LR[, rep(seq(d[2]), each = d[2]) < rep(seq(d[2]), d[2]), drop = FALSE]

  LR <- t(data.frame(purrr::map(LR, tmpfn),
                      check.names = FALSE))

  ## Organize as left_right, use pos to get pvalue
  pos <- LR[,"pos"]
  if(length(pos) == 1) { # restore colnames
    names(pos) <- rownames(LR)
  }
  nz <- LR[,"nz"]
  names(nz) <- NULL

  LR2 <- dplyr::mutate(
    dplyr::rename(
      dplyr::mutate(
        left_right(pos),
        nz = rep(nz, 2),
        Z = ifelse(Z > 0, Z, nz + Z)),
      pos = Z),
    pv = pbinom(pos - 1, nz, 0.5, lower.tail = FALSE))

  # Compare reference model with all others and get max pvalue.
  dplyr::mutate(
    dplyr::ungroup(
      dplyr::summarize(
        dplyr::group_by(
          dplyr::mutate(LR2,
                        ref = factor(ref, unique(ref))),
          ref),
        alt = alt[which.max(pv)][1],
        pv = max(pv))),
    ref = as.character(ref))
}
