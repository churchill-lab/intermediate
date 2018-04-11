# From qtl2pattern::covar_df_mx
covar_df_mx <- function(addcovar) {
  if(is.null(addcovar))
    return(NULL)
  if(is.data.frame(addcovar)) {
    f <- formula(paste("~", paste(names(addcovar), collapse = "+")))
    addcovar <- model.matrix(f, addcovar)[,-1, drop = FALSE]
  }
  wh_sex(addcovar)
}
# qtl2pattern:::wh_sex
wh_sex <- function(addcovar) {
  # Figure out which column is sex and make sure its name is "sex" 
  m <- match("sexm", tolower(colnames(addcovar)))
  if(is.na(m))
    m <- match("sex", tolower(colnames(addcovar)))
  if(!is.na(m))
    colnames(addcovar)[m] <- "sex"
  
  addcovar
}
