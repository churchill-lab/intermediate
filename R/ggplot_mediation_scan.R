#' GG plot of conditioned LOD scores against genomic positions
#'
#' Plot LOD statistics calculated by \code{mediation.scan} against genomic positions
#' using ggplot2.
#'
#' @param x mediation object
#' @param col color of points (default "firebrick4")
#' @param cex character expansion (default 2)
#' @param ylab Y axis label (default "Conditioned LOD")
#' @param ... additional arguments to \code{\link[qtl2ggplot]{ggplot_scan1}}
#' 
#' @seealso \code{\link{identify.mediation}}, \code{\link{kplot}}
#' 
#' @examples
#' data(Tmem68)
#' med <- mediation.scan(target=Tmem68$target,
#'                       mediator=Tmem68$mediator,
#'                       annotation=Tmem68$annotation,
#'                       covar=Tmem68$covar,
#'                       qtl.geno=Tmem68$qtl.geno)
#' ggplot2::autoplot(med)
#' 
#' @export
#' @importFrom qtl2ggplot ggplot_scan1
#' @importFrom ggplot2 geom_hline

ggplot_mediation_scan <- function(x, 
                           col="firebrick4",
                           cex = 1,
                           ylab = "Conditioned LOD",
                           col_target = "blue",
                           ...) {
  if(!is.factor(x$chr))
    x$chr <- factor(x$chr, unique(x$chr))
  x <- dplyr::arrange(x, chr, pos)
  map <- x$pos
  rownames(x) <- names(map) <- x$id
  map <- split(map, x$chr)
  p <- qtl2ggplot::ggplot_scan1(as.matrix(x[,"lod", drop = FALSE]), 
                           map, 
                           lines = FALSE,
                           col = col, cex = cex, 
                           ylab = ylab,
                           ...)
  targetFit <- attr(x, "targetFit")
  if(!is.null(targetFit))
    p <- p + ggplot2::geom_hline(yintercept = targetFit, col = col_target)
  p
}
#' @export
#' @export autoplot.mediation_scan
#' @method autoplot mediation_scan
#' @rdname ggplot_mediation_scan
#'
#' @importFrom ggplot2 autoplot
#'
autoplot.mediation_scan <- function(x, ...) {
  ggplot_mediation_scan(x, ...)
}
