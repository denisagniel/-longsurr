make_semimetric_pca <- function(k) {
  function(x, y) semimetric.pca(x, y, q = k)
}

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
  
  pca1_fit <- fregre.np.cv(fdataobj = fdX_t, y = y_t, 
                           metric = make_semimetric_pca(1))
  pca2_fit <- fregre.np.cv(fdataobj = fdX_t, y = y_t, 
                           metric = make_semimetric_pca(2))
  pca3_fit <- fregre.np.cv(fdataobj = fdX_t, y = y_t, 
                           metric = make_semimetric_pca(3))
  sspline_fit <- fgam(y_t ~ af(X_t, basistype = 's'))
  
  
  pca1_yhat = predict(pca1_fit, fdX_c)
  pca2_yhat = predict(pca2_fit, fdX_c)
  pca3_yhat = predict(pca3_fit, fdX_c)
  sspline_yhat = predict(sspline_fit, newdata = list(X_t = X_c), type = 'response')
  
  k_sspline_fit <- R.s.estimate(sone = predict(sspline_fit, type ='response'),
                                szero = sspline_yhat,
                                yone = y_t,
                                yzero = y_c, var = FALSE)
  k_pca1_fit <- R.s.estimate(sone = predict(pca1_fit, fdX_t),
                                szero = pca1_yhat,
                                yone = y_t,
                                yzero = y_c, var = FALSE)
  k_pca2_fit <- R.s.estimate(sone = predict(pca2_fit, fdX_t),
                             szero = pca2_yhat,
                             yone = y_t,
                             yzero = y_c, var = FALSE)
  k_pca3_fit <- R.s.estimate(sone = predict(pca3_fit, fdX_t),
                             szero = pca3_yhat,
                             yone = y_t,
                             yzero = y_c, var = FALSE)
  k_k_fit <- R.s.estimate(sone = predict(k_fit, fdX_t),
                             szero = k_yhat,
                             yone = y_t,
                             yzero = y_c, var = FALSE)
  
  
  delta_ss <-
    tibble(k = mean(k_yhat) - mean(y_c),
           pca1 = mean(pca1_yhat) - mean(y_c),
           pca2 = mean(pca2_yhat) - mean(y_c),
           pca3 = mean(pca3_yhat) - mean(y_c),
           kpca1 = k_pca1_fit$delta.s,
           kpca2 = k_pca2_fit$delta.s,
           kpca3 = k_pca3_fit$delta.s,
              fgam = mean(fgam_yhat) - mean(y_c),
           kfgam = k_fgam_fit$delta.s,
           sspline = mean(sspline_yhat) - mean(y_c),
           ksspline = k_sspline_fit$delta.s,
              lin = mean(lin_yhat) - mean(y_c),
           klin = k_lin_fit$delta.s,
              mean = mean_fit$delta.s,
              change = change_fit$delta.s)
  delta_ss
}
