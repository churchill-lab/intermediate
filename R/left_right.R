left_right <- function(Z) {
  out <-
    as.data.frame(
      purrr::simplify_all(
        purrr::transpose(
          stringr::str_split(names(Z), ":"))))
  names(out) <- c("ref", "alt")
  out <- dplyr::mutate(out,
                       ref = as.character(ref),
                       alt = as.character(alt),
                       pair = paste(ref, alt, sep = ":"))
  out$Z <- Z
  out2 <- out[,c(2,1,3,4)]
  out2$Z <- -Z
  names(out2) <- names(out)
  dplyr::bind_rows(out, out2)
}
