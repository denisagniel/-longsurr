

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
