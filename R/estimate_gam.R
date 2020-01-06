#' Estimate residual treatment effect after taking into account longitudinal surrogate marker with a GAM' 
#' @param y_t vector of n1 outcome measurements for treatment group
#' @param y_c vector of n0 outcome measurements for control or reference group
#' @param X_t n1 x T matrix of longitudinal surrogate measurements for treatment group
#' @param X_c n0 x T matrix of longitudinal surrogate measurements for control or reference group 
#' @return estimate of residual treatment effect
#' 
#' @examples
#' library(dplyr)
#' library(longsurr)
#' full_data <- 
#' generate_discontinuous_data(n = 50, n_i = 5, delta_s = 0.5, 
#' k = 1, s_y = 0.1, s_x = 0.1)$full_ds
#' 
#' 
#' wide_ds <- full_data %>% 
#' dplyr::select(id, a, tt, x, y) %>%
#' tidyr::spread(tt, x) 
#' 
#' wide_ds_0 <- wide_ds %>% filter(a == 0)
#' wide_ds_1 <- wide_ds %>% filter(a == 1)
#' X_t <- wide_ds_1 %>% dplyr::select(`-1`:`1`) %>% as.matrix
#' y_t <- wide_ds_1 %>% pull(y)
#' X_c <- wide_ds_0 %>% dplyr::select(`-1`:`1`) %>% as.matrix
#' y_c <- wide_ds_0 %>% pull(y)
#' 
#' estimate_gam(y_t = y_t, y_c = y_c, X_t = X_t, X_c = X_c)
#' 
#' @export

estimate_gam <- function(y_t, y_c, X_t, X_c, verbose = FALSE, ...) {
  stopifnot(length(y_t) == nrow(X_t))
  stopifnot(length(y_c) == nrow(X_c))
  stopifnot(ncol(X_t) == ncol(X_c))
  
  fgam_fit <- fgam(y_t ~ af(X_t))
  fgam_yhat = predict(fgam_fit, newdata = list(X_t = X_c), type = 'response')  
  k_fgam_fit <- Rsurrogate::R.s.estimate(sone = predict(fgam_fit, type ='response'),
                             szero = fgam_yhat,
                             yone = y_t,
                             yzero = y_c, var = FALSE, ...)
  if (verbose) {
    unsm_delta_s <- mean(fgam_yhat) - mean(y_c)
    print(glue::glue("Unsmoothed GAM result is {unsm_delta_s}; smoothed GAM results is {k_fgam_fit$delta.s}."))
  }
  k_fgam_fit$delta.s
}
