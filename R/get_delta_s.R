get_delta_s <- function(y_t = NULL, y_c = NULL, X_t = NULL, X_c = NULL) {
  fdX_t <- fdata(X_t)
  fdX_c <- fdata(X_c)
  
  k_fit <- fregre.np.cv(fdataobj = fdX_t, y = y_t)
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
  
  
  
  k_yhat = predict(k_fit, fdX_c)
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
  
  
  delta_ss <-
    tibble(delta_s_k = mean(k_yhat) - mean(y_c),
              delta_s_fgam = mean(fgam_yhat) - mean(y_c),
           delta_s_kfgam = k_fgam_fit$delta.s,
              delta_s_lin = mean(lin_yhat) - mean(y_c),
           delta_s_klin = k_lin_fit$delta.s,
              delta_s_mean = mean_fit$delta.s,
              delta_s_change = change_fit$delta.s)
  delta_ss
}
