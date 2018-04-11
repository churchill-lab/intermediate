#' @export
#' @rdname mediation_test
autoplot.mediation_test <- function(x, ...)
  ggplot_mediation_test(x, ...)
#' @export
#' @rdname mediation_test
ggplot_mediation_test <- function(x, type = c("pos_lod","pos_pvalue","pvalue_lod","alleles","mediator"),
                               main = params$target,
                               maxPvalue = 0.1, 
                               local_only = FALSE, 
                               significant = TRUE,
                               lod = TRUE,
                               ...) {
  type <- match.arg(type)
  if(is.null(local_only))
    local_only <- FALSE
  if(is.null(significant))
    significant <- TRUE
  
  params <- x$params
  pos_tar <- params$pos
  unmediated <- params$LR
  
  targetFit <- x$targetFit
  x <- x$best
  
  if(lod) {
    unmediated <- unmediated / log(10)
    x$mediation <- x$mediation / log(10)
  }
  
  if(!("symbol" %in% names(x)))
    x <- dplyr::rename(x, symbol = id)
  
  if(tmp <- any(is.na(x$triad))) {
    warning(paste("some triad are NA:", sum(tmp)))
  }
  
  relabel <- levels(x$triad)

  if(type != "pos_lod")
    significant <- TRUE
  if(significant) {
    x <- dplyr::filter(x, 
                       x$pvalue <= maxPvalue)
  } else {
    relabel <- c(relabel, 
                 paste0("n.s. (p>", round(maxPvalue, 2), ")"))
    tmp <- as.character(x$triad)
    tmp[x$pvalue > maxPvalue] <- relabel[5]
    x$triad <- factor(tmp, levels = relabel)
  }
  x <- dplyr::arrange(x, dplyr::desc(triad))
  
  # For expression, use qtl_pos if not missing.
  if(params$data_type == "expression" & !type %in% c("alleles","mediator")) {
    if(local_only)
      x <- dplyr::filter(x, local)
    else {
      x <- dplyr::mutate(x, pos = ifelse(local, pos, qtl_pos))
    }
    
    if(!is.null(x$local) & !is.null(x$qtl_ct)) {
      # Set up plot symbol.
      shapes <- c(17,16,2,1)
      names(shapes) <- c("distal", "local", "distal_info", "local_info")
      x <- dplyr::mutate(x, shape = names(shapes)[1 + local + 2 * (qtl_ct > 1)])
    }
  }

  # Colors
  cols <- c(RColorBrewer::brewer.pal(4, "Dark2"), "#CCCCCC")
  names(cols) <- relabel
  
  if(!(type %in% c("alleles","mediator"))) {
    switch(type,
           pos_pvalue = {
             p <- ggplot2::ggplot(x) +
               ggplot2::aes(x=pos, y=-log10(pvalue)) +
               ggplot2::aes(symbol=symbol, mediation=mediation) +
               ggplot2::facet_grid(~triad, scales = "free_x") +
               ggplot2::xlab("Position (Mbp)") +
               ggplot2::ylab("-log10 of p-value")
             if(!is.null(pos_tar))
               p <- p +
                 ggplot2::geom_vline(xintercept = pos_tar, col = "darkgrey")
           },
           pvalue_lod = {
             p <- ggplot2::ggplot(x) +
               ggplot2::aes(y=mediation, x=-log10(pvalue)) +
               ggplot2::aes(symbol=symbol, position=pos) +
               ggplot2::facet_grid(~triad, scales = "free_x") +
               ggplot2::geom_hline(yintercept = unmediated, col = "darkgrey") +
               ggplot2::xlab("-log10 of p-value") +
               ggplot2::ylab("Mediation LOD")
           },
           pos_lod = {
             p <- ggplot2::ggplot(x) + 
               ggplot2::aes(y=mediation, x=pos) +
               ggplot2::aes(symbol=symbol, pvalue=pvalue) +
               ggplot2::geom_hline(yintercept = unmediated, col = "darkgrey") +
               ggplot2::facet_grid(~triad, scales = "free_x") +
               ggplot2::xlab("Position (Mbp)") +
               ggplot2::ylab("Mediation LOD")
#               ggplot2::scale_color_manual(values = cols)
             if(!is.null(pos_tar))
               p <- p +
                 ggplot2::geom_vline(xintercept = pos_tar, col = "darkgrey")
           })
    if("biotype" %in% names(x)) {
      p <- p + ggplot2::aes(col = biotype)
    }
    if("info" %in% names(x)) {
      p <- p + ggplot2::aes(info = info)
    }
    if(exists("shapes")) {
      p <- p + ggplot2::geom_point(aes(shape = shape), size = 2, alpha = 0.5) +
        ggplot2::aes(chr = chr, qtl_pos = qtl_pos) +
        ggplot2::scale_shape_manual(values = shapes)
    } else {
      p <- p + 
        ggplot2::geom_point(size = 2, alpha = 0.5)
    }
    
    p + 
      ggplot2::ggtitle(main)
  } else {
    # Plot target above mediators
    grid::grid.newpage()
    grid::pushViewport(
      grid::viewport(
        layout = grid::grid.layout(nrow = 2)))
    
    plotfn <- function(x, type, targetCoef, layout.pos.row, ylabel = type) {
      p <- ggplot2::ggplot(allele_prep(x, type)) +
        ggplot2::aes(x=-log10(pvalue), y=value, col = geno, symbol = symbol) +
        ggplot2::facet_wrap(~triad, scales = "free_x") +
        ggplot2::scale_color_manual(values = CCSanger::CCcolors) +
        ggplot2::geom_hline(data = targetCoef, 
                            ggplot2::aes(yintercept = value, col = geno),
                            linetype = "dashed") +
        ggplot2::geom_point(size = 2) +
        ggplot2::ylab(ylabel)
      if(type == "alleles") {
        p <- p + ggplot2::ggtitle(main) +
          ggplot2::theme(axis.title.x = ggplot2::element_blank(),
                         axis.ticks.x = ggplot2::element_blank(),
                         axis.text.x = ggplot2::element_blank())
      } else {
        p <- p + 
          ggplot2::theme(strip.text = ggplot2::element_blank())
      }

      print(p,
            vp = grid::viewport(layout.pos.row = layout.pos.row,
                                layout.pos.col = 1))
    }
    
    targetCoef <- target_prep(targetFit, "alleles")
    switch(type,
           alleles = {
             plotfn(x, "mediator",
                    dplyr::mutate(targetCoef,
                                  value = (value - min(value)) / diff(range(value))),
                    2, "mediator effects")
             plotfn(x, "alleles", targetCoef, 1, "allele effects")
           },
           mediator = {
             plotfn(x, "alleles", targetCoef, 1, "allele effects")
             plotfn(x, "mediator",
                    dplyr::mutate(targetCoef,
                                  value = (value - min(value)) / diff(range(value))),
                    2, "mediator effects")
           })
  }
}
target_prep <- function(targetFit, type, col = CCSanger::CCcolors) {
  targetFit <- as.data.frame(t(targetFit$coef))
  codes <- LETTERS[seq_along(col)]
  m <- switch(type,
              alleles  = match(codes, names(targetFit)),
              mediator = match(paste(codes, "m", sep = "_"), names(targetFit)))
  
  col_names <- names(col)
  targetFit <- targetFit[m]
  names(targetFit) <- col_names
  out <- tidyr::gather(targetFit, geno, value)
  out <- dplyr::mutate(out, geno = factor(geno, col_names))
  switch(type,
         alleles  = dplyr::mutate(out, value = value - mean(value)),
         mediator = dplyr::mutate(out, value = (value - min(value)) / diff(range(value))))
}
allele_prep <- function(x, type, col = CCSanger::CCcolors) {
  codes <- LETTERS[seq_along(col)]
  m <- switch(type,
              alleles  = match(codes, names(x)),
              mediator = match(paste(codes, "m", sep = "_"), names(x)))
  
  col_names <- names(col)
  names(x)[m] <- col_names
  out <- dplyr::group_by(
    tidyr::gather(
      dplyr::select(
        x,
        symbol, triad, pvalue, chr, qtl_pos, dplyr::one_of(col_names)),
      geno, value, -symbol, -triad, -pvalue, -chr, -qtl_pos
    ),
    symbol
  )
  out <- dplyr::mutate(out, geno = factor(geno, col_names))
  out <- switch(type,
                alleles  = dplyr::mutate(out, value = value - mean(value)),
                mediator = dplyr::mutate(out, value = (value - min(value)) / diff(range(value))))
  dplyr::arrange(
    dplyr::ungroup(out),
    pvalue)
}