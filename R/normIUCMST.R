#' @export
#' 
normIUCMST <- function(models,
                             Zscores = calcZ(models, ...),
                             ...) {

  if(length(models$LR) == 1)
    return(data.frame(ref = names(models$LR),
                     alt = "",
                     pv = 1))
  
  # Expand to data frame with ref, alt, Z.
  Zscores <- left_right(Zscores)
  
  # Add p-value
  Zscores$pv <- pnorm(Zscores$Z, lower.tail = FALSE)

  # Compare reference model with all others and get max pvalue.
  dplyr::mutate(
    dplyr::ungroup(
      dplyr::summarize(
        dplyr::group_by(
          dplyr::mutate(Zscores, 
                        ref = factor(ref, unique(ref))),
          ref),
        alt = alt[which.max(pv)][1],
        pv = max(pv))),
    ref = as.character(ref))
}
