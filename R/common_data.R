common_data <- function(target, mediator, driver,
                        covar_tar=NULL, covar_med=NULL, kinship=NULL,
                        driver_med=NULL,
                        common = TRUE,
                        minN = 100, minCommon = 0.9) {

  # Make sure all are matrices
  target <- as.matrix(target)
  if(is.null(colnames(target)))
    colnames(target) <- "T"
  if(is.null(rownames(target)))
    rownames(target) <- seq_len(nrow(target))

  driver <- as.matrix(driver)
  if(is.null(rownames(driver)))
    rownames(driver) <- rownames(target)
  if(is.null(colnames(driver)))
    colnames(driver) <- "D"

  mediator <- as.matrix(mediator)
  if(is.null(rownames(mediator)))
    rownames(mediator) <- rownames(target)
  if(is.null(colnames(mediator)))
    colnames(mediator) <- "M"

  if(!is.null(covar_tar)) {
    covar_tar <- as.matrix(covar_tar)
    if(is.null(colnames(covar_tar)))
      colnames(covar_tar) <- paste0("covT", seq_len(ncol(covar_tar)))
  }
  if(!is.null(covar_med)) {
    covar_med <- as.matrix(covar_med)
    if(is.null(colnames(covar_med)))
      colnames(covar_med) <- paste0("covM", seq_len(ncol(covar_med)))
  }
  if(!is.null(driver_med)) {
    driver_med <- as.matrix(driver_med)
    if(is.null(colnames(driver_med)))
      colnames(driver_med) <- paste0("driverM", seq_len(ncol(driver_med)))
  }
  
  # Keep individuals with full records.
  ind2keep <-
    qtl2::get_common_ids(driver, target, covar_tar, covar_med, kinship, driver_med,
                             complete.cases = TRUE)
  
  # Drop mediator columns with too few non-missing data.
  if(enough <- (length(ind2keep) >= minN)) {
    m <- match(ind2keep, rownames(mediator), nomatch = 0)
    ind2keep <- ind2keep[m > 0]
    mediator <- mediator[m,, drop = FALSE]
    
    # This way considers only ind with no missing data.
    # Might want another way if pattern of missing different for some mediators.
    # Then would not do decomp_kinship and need to do common_data for single mediator
    
    # Count as number if not infinite and not missing
    is_num <- function(x) { !is.na(x) & is.finite(x) }
    # Drop mediators with too little data.
    ok_med <- apply(mediator, 2, function(x) sum(is_num(x))) >= minN
    if(enough <- any(ok_med)) {
      mediator <- mediator[, ok_med, drop = FALSE]
      
      # Check for missing across all remaining mediators.
      allMed <- apply(mediator, 1, function(x) any(is_num(x)))
      mediator <- mediator[allMed,, drop = FALSE]
      ind2keep <- ind2keep[allMed]
      if(enough <- (length(ind2keep) >= minN)) {
        common <- common & (sum(ok_med) == 1) | all(is_num(mediator))
      }
    }
  }
  if(!enough) {
    warning(paste0("too few data (", length(ind2keep), ") for analysis"))
    return(NULL)
  }
  driver <- driver[ind2keep,, drop = FALSE]
  target <- target[ind2keep,, drop = FALSE]
  mediator <- mediator[ind2keep,, drop = FALSE]
  if(!is.null(covar_tar))
    covar_tar <- covar_tar[ind2keep,, drop = FALSE]
  if(!is.null(covar_med))
    covar_med <- covar_med[ind2keep,, drop = FALSE]
  if(!is.null(driver_med))
    driver_med <- driver_med[ind2keep,, drop = FALSE]
  if(!is.null(kinship)) {
    kinship <- kinship[ind2keep, ind2keep]
    # Decompose kinship if all in common.
    if(common)
      kinship <- qtl2::decomp_kinship(kinship)
  }
  list(driver = driver,
       target = target,
       mediator = mediator,
       kinship = kinship,
       covar_tar = covar_tar,
       covar_med = covar_med,
       driver_med = driver_med,
       common = common)
}
