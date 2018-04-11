cmst_default <- function(object, driver, target, 
                         kinship, cov_tar, cov_med,
                         driver_med,
                         fitFunction, testFunction,
                         common = TRUE) {
  
  # Force x (= mediator column) to be matrix.
  mediator <- as.matrix(object[[1]])
  rownames(mediator) <- rownames(driver)
  colnames(mediator) <- "mediator"
  if(!is.null(driver_med))
    driver_med <- driver_med[,, object[[2]]$driver]
  
  # Make sure covariates are numeric
  cov_med <- covar_df_mx(cov_med)

  # Fit models
  fits <- med_fits(driver, target, mediator, fitFunction,
                   kinship, cov_tar, cov_med, driver_med,
                   common = common)
  
  combos <- combo_models()
  models <-
    purrr::transpose(
      purrr::map(combos[,1:4],
                 combine_models, fits[c("LR", "indLR", "df")]))
  
  models$LR <- unlist(models$LR)
  models$indLR <- as.data.frame(models$indLR)
  models$df <- unlist(models$df)
  models$coef <- fits$coef
  
  compsLR <- apply(fits$LR * combos[,5:7], 2, sum)

  # CMST on key models. Pick first as best.
  out <- head(dplyr::rename(
    dplyr::filter(
      testFunction(models),
      pv == min(pv)),
    pvalue = pv), n = 1L)
  
  out$ref <- factor(out$ref, c("causal","reactive","independent","correlated"))
  out$alt <- factor(out$alt, c("causal","reactive","independent","correlated"))
  
  # Mediation LR
  out$mediation <- compsLR["mediation"] # t.md_t.m
  
  # Mediator LR
  out$LRmed <- compsLR["mediator"]
  
  # Coefficients
  coef_target <- as.data.frame(t(models$coef$t.md_t.m[seq_len(ncol(driver))]))
  coef_mediator <- as.data.frame(t(models$coef$m.d_m[seq_len(ncol(driver))]))
  names(coef_mediator) <- paste0(names(coef_mediator), "_m")
  
  dplyr::bind_cols(out, coef_target, coef_mediator)
}
combo_models <- function() {
  combos <- 
    matrix(
      0, 5, 7,
      dimnames = list(
        c("t.d_t", "t.md_t.m", "m.d_m", "t.m_t", "m.t_m"),
        c("causal", "reactive", "independent", "correlated",
          "target", "mediator", "mediation")))
  combos[c(3,4), 1] <- 1 # causal: m.d_t.m
  combos[c(1,5), 2] <- 1 # reactive: t.d_m.t
  combos[c(1,3), 3] <- 1 # independent: t.d_m.d
  combos[   2:4, 4] <- 1 # correlated: t.md_m.d
  combos[     1, 5] <- 1 # target contrast: t.d_t
  combos[     3, 6] <- 1 # mediator contrast: m.d_m
  combos[     2, 7] <- 1 # mediation contrast: t.md_t.m
  as.data.frame(combos)
}
combine_models <- function(combos, fits) {
  list(LR = sum(fits$LR * combos),
       indLR = fits$indLR %*% combos,
       df = sum(fits$df * combos),
       coef = fits$coef)
}

cmst_pheno <- function(object, driver, target, 
                       kinship, cov_tar, cov_med,
                       driver_med,
                       fitFunction, testFunction,
                       common = TRUE) {
  
  # Currently, mediation_test uses elements of object[[2]] (columns of annot data frame)
  # to assess TRUE/FALSE on covariate columns. This will likely change.
  
  # Get covariate names appropriate for mediator 
  cov_names <- unlist(object[[2]][colnames(cov_med)])
  if(length(cov_names))
    cov_med <- cov_med[, cov_names, drop = FALSE]
  else
    cov_med <- NULL
  
  cmst_default(object, driver, target, 
               kinship, cov_tar, cov_med,
               driver_med,
               fitFunction, testFunction,
               common)
}  
