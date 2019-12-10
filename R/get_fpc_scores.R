#' Get FPC scores
#' 
#' @param ds data.frame containing the data
#' @param ycol string denoting the column in ds that contains the longitudinal marker for which FPC scores are desired
#' @param tcol string denoting the column in ds that contains the time
#' @param idcol string denoting the column in ds that contains the identifier of an individual
#' @param options list of options to pass to the FPC fitting procedure
#' @param return_eigenfunction logical to indicate whether the eigenfunctions should be returned as well
#' 
#' @importFrom stringr str_c
#' @seealso fdapace::FPCA
#' @export
get_fpc_scores <- function(ds, ycol, tcol, idcol, options = list(plot = TRUE), return_eigenfunction = FALSE) {
  list_data <- dataframe_to_list(ds, ycol, tcol, idcol)
  fpca_result <- with(list_data, fdapace::FPCA(y_list, t_list, options))
  score_ds <- data.frame(fpca_result$xiEst)
  colnames(score_ds) <- stringr::str_c('xi.', 1:ncol(score_ds))
  score_ds$id <- names(list_data$t_list)
  
  if (return_eigenfunction) {
    phi <- fpca_result$phi
    phi_hat <- approx_phi(fpca_result$phi, fpca_result$workGrid, ds %>% select_(tcol) %>% unlist)
    phi_ds <- data.frame(id = ds %>% select_(idcol), tt = ds %>% select_(tcol), phi = phi_hat)
    list(score_ds = score_ds,
         eigenfunction_ds = phi_ds,
         fpca_result = fpca_result
    )
  } else {
    list(score_ds = score_ds)
  }
}
