#' Turn genotype probabilities into labels
#' 
#' @param genos matrix of genotype probabilities at locus
#' @param sdp SNP distribution pattern for plot colors
#' @param allele Driver has alleles if \code{TRUE}, otherwise allele pairs.
#' 
#' @export
pattern_label <- function(genos, allele = TRUE) {
  geno_names <- colnames(genos)

  # Set up allele pair labels and count as color.
  if(allele)
    apply(genos, 1, pattern_allele_pair, geno_names)
  else
    # Pick geno name with highest value; if >1, paste together.
    apply(genos, 1, function(x) paste(geno_names[which.max(x)], collapse = ""))
}
#' @param sdp SNP distribution pattern for plot colors
#' @param label character string from \code{\link{pattern_label}}
#' @param geno_names unique genotype names (alleles or allele pairs)
#' 
#' @export
#' @rdname pattern_label
#' 
pattern_sdp <- function(label, sdp = NULL, geno_names = sort(unique(label))) {
  haplos <- unique(unlist(stringr::str_split(geno_names, "")))
  
  if(is.null(sdp)) {
    label
  } else {
    alt <- haplos[sdp_to_logical(sdp, haplos)]
    as.character(sapply(stringr::str_split(label, ""), 
                 function(x, alt) sum(x %in% alt), alt))
  }
}

pattern_allele_pair <- function(x, geno_names) {
  # Figure allele pair from allele probabilities. 
  # Mostly genos are close to 0, 0.5, 1 and add to 1.
  
  # order by decreasing value
  o <- order(-x)
  # find how many needed to exceed .67
  i <- 1 + sum(cumsum(x[o]) <= .67)
  # Paste together top 2, but first sort to  be in allele order
  paste(rep(geno_names[sort(o[seq_len(i)])], length = 2), collapse = "")
}
sdp_to_logical <- function(sdp, haplos) {
  assertthat::assert_that(!missing(haplos))
  sapply(sdp, function(x, haplos) {
    as.logical(intToBits(x)[seq_along(haplos)])
  }, haplos)
}
