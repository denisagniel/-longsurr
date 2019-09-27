#' Estimate the surrogate value of a longitudinal marker
#' 
#' @param y_t vector of n1 outcome measurements for treatment group
#' @param y_c vector of n0 outcome measurements for control or reference group
#' @param X_t n1 x T matrix of longitudinal surrogate measurements for treatment group
#' @param X_c n0 x T matrix of longitudinal surrogate measurements for control or reference group 
#' @param method method for dimension-reduction of longitudinal surrogate, either 'gam', 'linear', or 'kernel'
#' @param bootstrap_samples number of bootstrap samples to use for variance estimation. The default is 0, which estimates without providing a variance estimate.
#' 
#' @return a tibble containing estimates, standard errors, and quantile-based confidence intervals for the residual treatment effect \code{Deltahat_s_*} and the proportion of treatment effect explained \code{R_*}
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
#' estimate_surrogate_value(y_t = y_t, y_c = y_c, X_t = X_t, X_c = X_c, method = 'kernel')
#' estimate_surrogate_value(y_t = y_t, y_c = y_c, X_t = X_t, X_c = X_c, method = 'linear', bootstrap_sample = 50)
#' @export

estimate_surrogate_value <- function(y_t, y_c, X_t, X_c, method = c('gam', 'linear', 'kernel'), k = 3, bootstrap_samples = 0, alpha = 0.05) {
  if (method == 'linear') {
    Deltahat_S <- estimate_linear(y_t, y_c, X_t, X_c)
  } else if (method == 'gam') {
    Deltahat_S <- estimate_gam(y_t, y_c, X_t, X_c)
  } else if (method == 'kernel') {
    Deltahat_S <- estimate_kernel(y_t, y_c, X_t, X_c, k)
  }
  Deltahat <- mean(y_t) - mean(y_c)
  if (bootstrap_samples > 0) {
    boot_ests <- purrr::map(1:bootstrap_samples, boot_fn, method, k)
    boot_ests <- dplyr::bind_rows(boot_ests)
    boot_se <- summarise_all(boot_ests, sd, na.rm = TRUE)
    boot_ci_l <- summarise_all(boot_ests, quantile, alpha/2, na.rm = TRUE)
    boot_ci_h <- summarise_all(boot_ests, quantile, 1-alpha/2, na.rm = TRUE)
  } else {
    boot_se <- boot_ci_l <- boot_ci_h <- 
      tibble(Deltahat = NA,
             Deltahat_S = NA,
             R = NA,
             Deltahat_S_se = NA,
             Deltahat_S_ci_l = NA,
             Deltahat_S_ci_h = NA,
             R_se = NA,
             R_ci_l = NA,
             R_ci_h = NA)
  }
  tibble::tibble(
    Deltahat = Deltahat,
    Deltahat_S = Deltahat_S,
    R = 1 - Deltahat_S/Deltahat,
    Deltahat_S_se = boot_se$Deltahat_S,
    Deltahat_S_ci_l = boot_ci_l$Deltahat_S,
    Deltahat_S_ci_h = boot_ci_h$Deltahat_S,
    R_se = boot_se$R,
    R_ci_l = boot_ci_l$R,
    R_ci_h = boot_ci_h$R
  )
}
  
boot_fn <- function(b, method, k) {
    n1 <- length(y_t)
    n0 <- length(y_c)
    
    ind_t <- sample(1:n1, replace = TRUE)
    ind_c <- sample(1:n0, replace = TRUE)
    
    boot_yt <- y_t[ind_t]
    boot_yc <- y_c[ind_c]
    boot_Xt <- X_t[ind_t,]
    boot_Xc <- X_c[ind_c,]
    estimate_surrogate_value(
      boot_yt, 
      boot_yc, 
      boot_Xt, 
      boot_Xc, 
      method, 
      k, 
      bootstrap_samples = 0)
}
