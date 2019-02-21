delta_sum_fn <- function(yt, yc, muhat_ct = NULL, muhat_tc = NULL,
                         delta = NULL,
                         delta_s = NULL,
                         R_s = NULL,
                         se_hat_ct = NULL,
                         se_hat_tc = NULL,
                         delta_var = NULL,
                         delta_s_var = NULL,
                         R_s_var = NULL) {
  # browser()
  n_t <- length(yt)
  n_c <- length(yc)
  muhat_ct <- if_else(is.null(muhat_ct),
                      as.double(NA),
                      muhat_ct)
  muhat_tc <- if_else(is.null(muhat_tc),
                      as.double(NA),
                      muhat_tc)
  delta <- if_else(is.null(delta),
                   mean(yt) - mean(yc),
                   delta)
  delta_s <- if_else(is.null(delta_s),
          mean(yt) - muhat_ct,
          delta_s)
  data.frame(type = c('mu_t',
                      'mu_c',
                      'mu_st',
                      'mu_sc',
                      'delta',
                      'delta_s',
                      'R'),
             est = c(mean(yt),
                     mean(yc),
                     muhat_ct,
                     muhat_tc,
                     delta,
                     delta_s,
                     if_else(is.null(R_s),
                             1 - delta_s/delta,
                             R_s)),
             se = c(sd(yt)/sqrt(n_t),
                    sd(yc)/sqrt(n_c),
                    if_else(is.null(se_hat_ct),
                            as.double(NA),
                            se_hat_ct),
                    if_else(is.null(se_hat_tc),
                            as.double(NA),
                            se_hat_tc),
                    if_else(is.null(delta_var),
                            sqrt(var(yt)/n_t + var(yc)/n_c),
                            delta_var),
                    if_else(is.null(delta_s_var),
                            as.double(NA),
                            delta_s_var),
                    if_else(is.null(R_s_var),
                            as.double(NA),
                            R_s_var)))
}

fit_linear_model <- function(y_t, y_c, X_t, X_c) {
  # browser()
  n_t <- length(y_t)
  n_c <- length(y_c)
  fit_t <- pfr(y_t ~ lf(X_t, k = 30, bs = 'ps'))
  fit_c <- pfr(y_c ~ lf(X_c, k = 30, bs = 'ps'))
# #
# #   predict_t_on_c <- predict(fit_c, newdata = list(X_c = X_t), type = 'response', se.fit = TRUE)
# #   predict_c_on_t <- predict(fit_t, newdata = list(X_t = X_c), type = 'response')
# #
#   yhat_tc <- mean(predict_t_on_c$fit)
  Xbar_t <- colMeans(X_t) %>% matrix(1, ncol(X_t))
  predict_yhat_tc <- predict(fit_c, newdata = list(X_c = Xbar_t), type = 'response', se.fit = TRUE)

  Xbar_c <- colMeans(X_c) %>% matrix(1, ncol(X_c))
  predict_yhat_ct <- predict(fit_t, newdata = list(X_t = Xbar_c), type = 'response', se.fit = TRUE)

    delta_sum_fn(yt = y_t,
                 yc = y_c,
                 muhat_ct = predict_yhat_ct$fit,
                 muhat_tc = predict_yhat_tc$fit,
                 se_hat_ct = predict_yhat_ct$se.fit,
                 se_hat_tc = predict_yhat_tc$se.fit)
}

fit_naive_mean_model <- function(y_t, y_c, X_t, X_c) {

  res <- R.s.estimate(sone = rowMeans(X_t),
                      szero = rowMeans(X_c),
                      yone = y_t,
                      yzero = y_c, var = TRUE)

  with(res, delta_sum_fn(yt = y_t,
               yc = y_c,
               delta = delta,
               delta_s = delta.s,
               R_s = R.s,
               delta_var = delta.var,
               delta_s_var = delta.s.var,
               R_s_var = R.s.var
               )
  )
}

fit_naive_change_model <- function(y_t, y_c, X_t, X_c) {
# browser()
  res <- R.s.estimate(sone = X_t[,ncol(X_t)] - X_t[,1],
                      szero = X_c[,ncol(X_c)] - X_c[,1],
                      yone = y_t,
                      yzero = y_c, var = TRUE)

  with(res, delta_sum_fn(yt = y_t,
                         yc = y_c,
                         delta = delta,
                         delta_s = delta.s,
                         R_s = R.s,
                         delta_var = delta.var,
                         delta_s_var = delta.s.var,
                         R_s_var = R.s.var
  )
  )
}

# fit_fpc_gam <- function(Xi_ds_t, Xi_ds_c, npcs) {

##
### i'm deleting this for now because it doesn't make sense without a common basis
##
#   xifm <- paste0('y ~ ', paste('s(xi.', 1:npcs, ')', collapse = ' + ', sep = ''))
#   # xifm_c <- paste0('y ~ ', paste('s(xi.', 1:npcs, ')', collapse = ' + ', sep = ''))
#
#   fit_t <- gam(as.formula(xifm), data = Xi_ds_t)
#   fit_c <- gam(as.formula(xifm), data = Xi_ds_c)
#
#   # predict_t_on_c <- predict(fit_c, newdata = Xi_ds_t, type = 'response')
#   # predict_c_on_t <- predict(fit_t, newdata = Xi_ds_c, type = 'response')
#   Xbar_t <- colMeans(X_t) %>% matrix(1, ncol(X_t))
#   predict_yhat_tc <- predict(fit_c, newdata = list(X_c = Xbar_t), type = 'response', se.fit = TRUE)
#
#   Xbar_c <- colMeans(X_c) %>% matrix(1, ncol(X_c))
#   predict_yhat_ct <- predict(fit_t, newdata = list(X_t = Xbar_c), type = 'response', se.fit = TRUE)
#
#     data.frame(type = c('mu_t', 'mu_c', 'mu_st', 'mu_sc'),
#                est = c(mean(y_t), mean(y_c), mean(predict_c_on_t), mean(predict_t_on_c)))
#
# }

fit_fgam <- function(y_t, y_c, X_t, X_c) {
  fit_t <- fgam(y_t ~ af(X_t))
  fit_c <- fgam(y_c ~ af(X_c))
  n_t <- length(y_t)
  n_c <- length(y_c)

  Xbar_t <- colMeans(X_t) %>% matrix(1, ncol(X_t))
  predict_yhat_tc <- predict(fit_c, newdata = list(X_c = Xbar_t), type = 'response', se.fit = TRUE)

  Xbar_c <- colMeans(X_c) %>% matrix(1, ncol(X_c))
  predict_yhat_ct <- predict(fit_t, newdata = list(X_t = Xbar_c), type = 'response', se.fit = TRUE)
#
#   data.frame(type = c('mu_t', 'mu_c', 'mu_st', 'mu_sc'),
#              est = c(mean(y_t), mean(y_c), predict_yhat_ct$fit, predict_yhat_tc$fit),
#              se = c(sd(y_t)/sqrt(n_t), sd(y_c)/sqrt(n_c), predict_yhat_ct$se.fit, predict_yhat_tc$se.fit))
  delta_sum_fn(yt = y_t,
               yc = y_c,
               muhat_ct = predict_yhat_ct$fit,
               muhat_tc = predict_yhat_tc$fit,
               se_hat_ct = predict_yhat_ct$se.fit,
               se_hat_tc = predict_yhat_tc$se.fit)

  data.frame(type = c('mu_t', 'mu_c', 'mu_st',
                      'mu_sc', 'delta', 'delta_s', 'R'),
             est = c(mean(y_t), mean(y_c), predict_yhat_ct$fit,
                     predict_yhat_tc$fit, mean(y_t) - mean(y_c),
                     mean(y_t) - predict_yhat_ct$fit,
                     1 - (mean(y_t) - predict_yhat_ct$fit)/
                       (mean(y_t) - mean(y_c))),
             se = c(sd(y_t)/sqrt(n_t), sd(y_c)/sqrt(n_c), predict_yhat_ct$se.fit, predict_yhat_tc$se.fit, sqrt(var(y_t)/n_t + var(y_c)/n_c), NA, NA))

  # predict_t_on_c <- predict(fit_c, newdata = list(X_c = X_t), type = 'response')
  # predict_c_on_t <- predict(fit_t, newdata = list(X_t = X_c), type = 'response')
  #
  # data.frame(type = c('mu_t', 'mu_c', 'mu_st', 'mu_sc'),
  #            est = c(mean(y_t), mean(y_c), mean(predict_c_on_t), mean(predict_t_on_c)))
}

fit_pco_model <- function(y_t, y_c, X_t, X_c) {
  # browser()
  n <- length(y_t)

  dummy <- rep(1, n)
  d_dtw_t <- dist(X_t, method="dtw")
  fit_np_t <- gam(y_t ~ s(dummy, bs="pco", k=4, xt=list(D=d_dtw_t)), method="REML")

  d_dtw_c <- dist(X_c, method="dtw")
  fit_np_c <- gam(y_c ~ s(dummy, bs="pco", k=4, xt=list(D=d_dtw_c)), method="REML")

  dist_t_to_c <- dist(X_t, X_c, method = 'dtw')
  dist_c_to_t <- t(dist_t_to_c)

  pred_data_t_on_c <- pco_predict_preprocess(fit_np_c, newdata=NULL, list(dummy = dist_c_to_t))
  yhat_t_on_c <- predict(fit_np_c, pred_data_t_on_c, se.fit = TRUE)
  pred_data_c_on_t <- pco_predict_preprocess(fit_np_t, newdata=NULL, list(dummy = dist_t_to_c))
  yhat_c_on_t <- predict(fit_np_t, pred_data_c_on_t)

    data.frame(type = c('mu_t', 'mu_c', 'mu_st', 'mu_sc'),
               est = c(mean(y_t), mean(y_c), mean(yhat_c_on_t), mean(yhat_t_on_c)))

}

fit_kernel_model <- function(y_t, y_c, X_t, X_c) {
  # browser()
  fd_Xt <- fdata(X_t)
  fd_Xc <- fdata(X_c)
  fit_t <- fregre.np(fd_Xt, y_t, Ker=AKer.epa)
  fit_c <- fregre.np(fd_Xc, y_c, Ker=AKer.epa)
  n_t <- length(y_t)
  n_c <- length(y_c)

  Xbar_t <- colMeans(X_t) %>% matrix(1, ncol(X_t))
  fd_Xbar_t <- fdata(Xbar_t)
  yhat_t_on_c <- predict.fregre.fd(fit_c, new.fdataobj = fd_Xbar_t, se.fit = TRUE)
  # old_yhat_t_on_c <- predict.fregre.fd(fit_c, new.fdataobj = fd_Xt, se.fit = TRUE)

  Xbar_c <- colMeans(X_c) %>% matrix(1, ncol(X_c))
  fd_Xbar_c <- fdata(Xbar_c)
  yhat_c_on_t <- predict.fregre.fd(fit_t, new.fdataobj = fd_Xbar_c, se.fit = TRUE)
  # old_yhat_c_on_t <- predict.fregre.fd(fit_t, new.fdataobj = fd_Xc, se.fit = TRUE)

  # yhat_t_on_c <- predict.fregre.fd(fit_c, new.fdataobj = fd_Xt)
  # yhat_c_on_t <- predict.fregre.fd(fit_t, new.fdataobj = fd_Xc)
  # data.frame(type = c('mu_t', 'mu_c', 'mu_st', 'mu_sc'),
  #            est = c(mean(y_t), mean(y_c), yhat_c_on_t$fit, yhat_t_on_c$fit),
  #            se = c(sd(y_t)/sqrt(n_t), sd(y_c)/sqrt(n_c), yhat_c_on_t$se.fit, yhat_t_on_c$se.fit))
  data.frame(type = c('mu_t', 'mu_c', 'mu_st',
                      'mu_sc', 'delta', 'delta_s', 'R'),
             est = c(mean(y_t), mean(y_c), yhat_c_on_t$fit,
                     yhat_t_on_c$fit, mean(y_t) - mean(y_c),
                     mean(y_t) - yhat_c_on_t$fit,
                     1 - (mean(y_t) - yhat_c_on_t$fit)/
                       (mean(y_t) - mean(y_c))),
             se = c(sd(y_t)/sqrt(n_t), sd(y_c)/sqrt(n_c), yhat_c_on_t$se.fit, yhat_t_on_c$se.fit, sqrt(var(y_t)/n_t + var(y_c)/n_c), NA, NA))
}

fit_oracle_models <- function(y_t, y_c, X_t, X_c) {
  #------------------------------------------
  ## naive mean method without error or sparsity
  #-------------------------------------------
  oracle_mean_res <-
    fit_naive_mean_model(y_t = y_t, y_c = y_c, X_t = X_t, X_c = X_c) %>%
    mutate(setting = 'oracle_mean')
  #------------------------------------------
  ## naive change method without error or sparsity
  #-------------------------------------------
  oracle_change_res <-
    fit_naive_change_model(y_t = y_t, y_c = y_c, X_t = X_t, X_c = X_c) %>%
    mutate(setting = 'oracle_change')
  #----------------------------------------------
  ## linear performance without error or sparsity
  #----------------------------------------------
  oracle_lin_res <-
    fit_linear_model(y_t = y_t, y_c = y_c, X_t = X_t, X_c = X_c) %>%
    mutate(setting = 'oracle_linear')
  #----------------------------------------------
  ## fgam performance without error or sparsity
  #----------------------------------------------
  oracle_fgam_res <-
    fit_fgam(y_t = y_t, y_c = y_c, X_t = X_t, X_c = X_c) %>%
    mutate(setting = 'oracle_fgam')
  #--------------------------------------
  ## oracle kernel method
  #------------------------------------
  oracle_kernel_res <-
    fit_kernel_model(y_t = y_t, y_c = y_c, X_t = X_t, X_c = X_c) %>%
    mutate(setting = 'oracle_kernel')



  full_res <- full_join(oracle_mean_res, oracle_change_res) %>%
    full_join(oracle_lin_res) %>%
    full_join(oracle_fgam_res) %>%
    full_join(oracle_kernel_res)
  full_res
}

fit_obs_models <- function(y_t, y_c, obs_data,
                           trt_xhat_wide, ctrl_xhat_wide) {
  wide_obs <- obs_data %>%
    group_by(id) %>%
    mutate(time_n = glue('t{rank(tt)}')) %>%
    ungroup %>%
    dplyr::select(id, a, y, time_n, x) %>%
    spread(time_n, x)
  trt_obs_x <- wide_obs %>%
    filter(a == 1) %>%
    arrange(id) %>%
    dplyr::select(t1:t4) %>%
    as.matrix
  ctrl_obs_x <- wide_obs %>%
    filter(a == 0) %>%
    arrange(id) %>%
    dplyr::select(t1:t4) %>%
    as.matrix

  #------------------------------------------
  ## naive mean method on observed data
  #-------------------------------------------
  obs_mean_res <-
    fit_naive_mean_model(y_t = y_t, y_c = y_c, X_t = trt_obs_x, X_c = ctrl_obs_x) %>%
    mutate(setting = 'obs_mean')
  #------------------------------------------
  ## naive change method on observed data
  #-------------------------------------------
  obs_change_res <-
    fit_naive_change_model(y_t = y_t, y_c = y_c, X_t = trt_obs_x, X_c = ctrl_obs_x) %>%
    mutate(setting = 'obs_change')
  #----------------------------------------------
  ## linear on observed data
  #----------------------------------------------
  obs_lin_res <-
    fit_linear_model(y_t = y_t, y_c = y_c, X_t = trt_xhat_wide, X_c = ctrl_xhat_wide) %>%
    mutate(setting = 'obs_linear')

  #----------------------------------------------
  ## fgam performance on observed data
  #----------------------------------------------
  obs_fgam_res <-
    fit_fgam(y_t = y_t, y_c = y_c, X_t = trt_xhat_wide, X_c = ctrl_xhat_wide) %>%
    mutate(setting = 'obs_fgam')

  #--------------------------------------
  ## observed kernel method
  #------------------------------------
  obs_kernel_res <-
    fit_kernel_model(y_t = y_t, y_c = y_c, X_t = trt_xhat_wide, X_c = ctrl_xhat_wide) %>%
    mutate(setting = 'obs_kernel')


  full_res <- full_join(obs_mean_res, obs_change_res) %>%
    full_join(obs_lin_res) %>%
    full_join(obs_fgam_res) %>%
    full_join(obs_kernel_res)
  full_res
}
