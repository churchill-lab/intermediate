med_fits <- function(driver, target, mediator, fitFunction,
                     kinship=NULL, cov_tar=NULL, cov_med=NULL,
                     driver_med = NULL,
                     common = FALSE, ...) {
  
  if(!common) {
    commons <- common_data(target, mediator, driver,
                           cov_tar, cov_med, kinship, driver_med)
    driver <- commons$driver
    target <- commons$target
    mediator <- commons$mediator
    kinship <- commons$kinship
    cov_tar <- commons$cov_tar
    cov_med <- commons$cov_med
    driver_med <- commons$driver_med
  }
  if(is.null(driver_med))
    driver_med <- driver
  
  # Fit mediation models.
  # Transpose list of model fits
  fits <- purrr::transpose(list(
    t.d_t    = fitFunction(driver, target, kinship, cov_tar),
    t.md_t.m = fitFunction(driver, target, kinship, cbind(cov_tar, mediator)),
    m.d_m    = fitFunction(driver_med, mediator, kinship, cov_med),
    t.m_t    = fitFunction(cbind(1, mediator), target, kinship, cov_tar),
    m.t_m    = fitFunction(cbind(1, target), mediator, kinship, cov_med)))
  fits$LR <- unlist(fits$LR)
  fits$indLR <- as.matrix(as.data.frame(fits$indLR))
  fits$df <- unlist(fits$df)
  
  fits
}
