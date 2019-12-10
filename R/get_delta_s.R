#' Estimate residual treatment effect
#' 
#' @param y_t vector of outcomes in the treatment arm.
#' @param y_c vector of outcomes in the control arm.
#' @param X_t matrix of surrogate values in the treatment arm.
#' @param X_c matrix of surrogate values in the control arm.
#' 
#' @return list containing matrices \code{X_t} and \code{X_c}, which are the smoothed surrogate values for the treated and control groups, respectively, for use in downstream analyses 
#' 
#' @examples
#' library(dplyr)
#' library(longsurr)
#' obs_data <- 
#' generate_discontinuous_data(n = 50, n_i = 5, delta_s = 0.5, 
#' k = 1, s_y = 0.1, s_x = 0.1)$obs_ds
#' 
#' head(obs_data)
#' presmooth_X <- presmooth_data(obs_data)
#' 
#' wide_ds <- full_data %>% 
#' dplyr::select(id, a, tt, x, y) %>%
#' tidyr::spread(tt, x) 
#' 
#' y_t <- wide_ds %>%
#' filter(a == 1) %>%
#' pull(y)
#' y_c <- wide_ds %>%
#' filter(a == 0) %>%
#' pull(y)
#' X_t <- presmooth_X$X_t
#' X_c <- presmooth_X$X_c
#' 
#' estimate_surrogate_value(y_t = y_t, y_c = y_c, 
#' X_t = X_t, X_c = X_c, method = 'linear')
#' 
#' @import dplyr
#' @export
#' 


get_delta_s <- function(y_t = NULL, y_c = NULL, X_t = NULL, X_c = NULL) {
  fdX_t <- fdata(X_t)
  fdX_c <- fdata(X_c)

  fgam_fit <- fgam(y_t ~ af(X_t))
  lin_fit <- pfr(y_t ~ lf(X_t))
  mean_fit <- R.s.estimate(sone = rowMeans(X_t),
                           szero = rowMeans(X_c),
                           yone = y_t,
                           yzero = y_c, var = FALSE)
  change_fit <- R.s.estimate(sone = X_t[,ncol(X_t)] - X_t[,1],
                             szero = X_c[,ncol(X_c)] - X_c[,1],
                             yone = y_t,
                             yzero = y_c, var = FALSE)
  

  fgam_yhat = predict(fgam_fit, newdata = list(X_t = X_c), type = 'response')
  lin_yhat = predict(lin_fit, newdata = list(X_t = X_c), type = 'response')
  
  k_lin_fit <- R.s.estimate(sone = predict(lin_fit, type = 'response'),
                            szero = lin_yhat,
                            yone = y_t,
                            yzero = y_c, var = FALSE)
  k_fgam_fit <- R.s.estimate(sone = predict(fgam_fit, type ='response'),
                             szero = fgam_yhat,
                             yone = y_t,
                             yzero = y_c, var = FALSE)
  
  pca2_fit <- fregre.np.cv(fdataobj = fdX_t, y = y_t, 
                           metric = make_semimetric_pca(2))
  pca3_fit <- fregre.np.cv(fdataobj = fdX_t, y = y_t, 
                           metric = make_semimetric_pca(3))
  pca4_fit <- fregre.np.cv(fdataobj = fdX_t, y = y_t, 
                           metric = make_semimetric_pca(4))
  pca10_fit <- fregre.np.cv(fdataobj = fdX_t, y = y_t, 
                           metric = make_semimetric_pca(10))
  
  
  pca2_yhat = predict(pca2_fit, fdX_c)
  pca3_yhat = predict(pca3_fit, fdX_c)
  pca4_yhat = predict(pca4_fit, fdX_c)
  pca10_yhat = predict(pca10_fit, fdX_c)
  
  
  delta_ss <-
    tibble(
           pca2 = mean(pca2_yhat) - mean(y_c),
           pca3 = mean(pca3_yhat) - mean(y_c),
           pca4 = mean(pca4_yhat) - mean(y_c),
           pca10 = mean(pca10_yhat) - mean(y_c),
              fgam = mean(fgam_yhat) - mean(y_c),
           kfgam = k_fgam_fit$delta.s,
              lin = mean(lin_yhat) - mean(y_c),
           klin = k_lin_fit$delta.s,
              mean = mean_fit$delta.s,
              change = change_fit$delta.s)
  delta_ss
}
