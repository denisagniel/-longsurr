lsa <- function(dat, msK = NULL, fpc_fit = NULL) {
  if (is.null(fpc_fit)) {
    fpc <- fpca(dat, 'x', 'tt', 'id', options = list(
      nRegGrid = 150,
      methodSelectK = ifelse(is.null(msK), 'AIC', msK),
      maxK = 60,
      useBinnedCov = FALSE))
  } else fpc <- fpc_fit

  score_dat <- fpc$score_ds %>% merge(dat)
# browser()
  ys0_ds <- score_dat %>% filter(a == 0) %>% select(-x, -tt) %>% unique
  ys1_ds <- score_dat %>% filter(a == 1) %>% select(-x, -tt) %>% unique

  xs_ds <- score_dat %>% mutate(id = as.numeric(id)) %>% inner_join(fpc$phi_ds) %>%
    mutate(x_mu = x - mu)
  # xis <- ys0_ds %>% select(contains('xi')) ## scores among untreated
  # n <- nrow(ys0_ds)
  #
  # ## now we estimate the relationship between the scores and the outcome
  # ## AMONG THE UNTREATED - when we plug the projected scores for the treated
  # ## into this function, we will be estimating the predicted outcome AMONG THE
  # ## UNTREATED if they had longitudinal surrogates like the treated
  #
  # yhat_fns <- lapply(xis, function(xx) {
  #   with(ys0_ds, KernSmooth::locpoly(pnorm(xx), y, bandwidth = 4*sd(pnorm(xx))*n^(-1/5)))
  # })
  # # browser()
  # yhat_fn <- function(new_xi) {
  #   id <- new_xi$id
  #   nx <- new_xi %>% select(-id)
  #   k <- ncol(nx)
  #   yh <- matrix(0, n, k)
  #   for (j in 1:k) {
  #     yh[,j] <- approx(yhat_fns[[j]]$x, yhat_fns[[j]]$y, xout=pnorm(nx[,j]), rule = 2)$y
  #     if (any(yh[,j] == -Inf)) browser()
  #   }
  #   # browser()
  #   tibble(id = id, yh = rowSums(yh))
  # }
  #
  # yhat <- yhat_fn(ys1_ds %>% select(id, contains('xi')))
  #
  # browser()
  n_xi <- ys0_ds %>% select(contains('xi')) %>% ncol
  xi_str <- stringr::str_c('s(xi.', 1:n_xi, collapse = ') + ')
  fm <- stringr::str_c('y ~ ', xi_str, ')')
  gam_fit <- mgcv::gam(as.formula(fm), data = ys1_ds)
  yhat <- predict(gam_fit, newdata = ys0_ds)

  cf_dat <- rbind(ys0_ds, ys1_ds)

  cf_y <- cf_dat$y
  cf_x <- cf_dat %>% select(contains('xi')) %>% as.matrix
  cf_w <- cf_dat$a

  grf_fit <- grf::causal_forest(Y = cf_y, X = cf_x, W = cf_w)
  grf_delta0 <- predict(grf_fit, ys0_ds %>% select(contains('xi')) %>% as.matrix)$predictions

  # n_phi <- min(4, n_xi)
  # phi_str <- stringr::str_c('phi.', 1:n_phi, collapse = ' + ')
  # xfm <- stringr::str_c('x_mu ~ ', phi_str, ' + (-1 + ', phi_str, '| id)')
  # xlm_fit <- lme4::lmer(as.formula(xfm), data = xs_ds)
  # xitilde <- ranef(xlm_fit)[[1]]

  # browser()
  n_b <- dat %>% select(contains('b')) %>% ncol
  b_str <- stringr::str_c('s(b', 1:n_b-1, collapse = ') + ')
  bfm <- stringr::str_c('y ~ ', b_str, ')') %>% as.formula
  oracle_fit_0 <- mgcv::gam(bfm, data = dat %>% filter(a == 1))
  oracle_yh <- predict(oracle_fit_0, newdata = dat %>% filter(a == 0))

  yh_ds <- data.frame(id = dat %>% filter(a == 0) %>% select(id) %>% unlist %>% as.factor, yh_o = oracle_yh) %>% unique %>%
    inner_join(data.frame(id = ys0_ds$id, yh_m = yhat, grf_delta = grf_delta0))

  # browser()
  list(yh_ds = yh_ds,
       # refit_xi = xitilde,
       fpca_fit = fpc)
}


lsa_old <- function(dat) {
  # browser()
  fpc.0 <- fpca(dat %>% filter(a == 0), 'x', 'tt', 'id')
  fpc.1 <- fpca(dat %>% filter(a == 1), 'x', 'tt', 'id')
  ####
  # we are going to use the FPC scores from the untreated (a = 0)
  # when we finally estimate the counterfactual, so we need the scores
  # we estimate in the treated (a = 1) to have the untreated basis.
  # in the following, we take the eigenfunctions from the untreated
  # and we project the treated marker values into the untreated eigenfn
  # space.



  ## treated ids
  all.ids.1 <- dat %>% filter(a == 1) %>% select(id) %>% unlist %>% unique
  sc.l <- list()
  for (idd in all.ids.1) {
    dat.i <- dat %>% filter(id == idd)
    tt.i <- dat.i$tt
    x.i <- dat.i$x
    # find and remove the mean of y.i
    xi.i <- project_scores(x.i,
                           tt.i,
                           fpc.0$fpca_result,
                           fpc.1$fpca_result)
    sc.l[[as.character(idd)]] <- data.frame(xi = t(xi.i))
  }
  sc1_on_0 <- bind_rows(sc.l, .id = 'id')
browser()



  dat_1 <- dat %>% filter(a == 1)
  dat_1 <- dat_1 %>% group_by(id) %>% mutate(yh = predict_fpca(idd = id, tt = tt, score_ds = sc1_on_0, fpca_phi = fpc.0$fpca_result, fpca_mu = fpc.1$fpca_result), yh_old = predict_fpca(idd = id, tt = tt, score_ds = sc1_on_0, fpca_phi = fpc.0$fpca_result, fpca_mu = fpc.0$fpca_result))

  dat_tst <- dat %>% group_by(id) %>% mutate(yh = predict_fpca(idd = id, tt = tt, score_ds = fpc_c$score_ds, fpca_phi = fpc_c$fpca_result))
  new_phi <- with(dat_1, approx_phi(
    phi_mat = fpc_0$fpca_result$phi,
    t_grid = fpc_0$fpca_result$workGrid,
    new_t = tt
  ))
  np <- ncol(new_phi)
  dat_10 <- with(dat_1, data.frame(
    id,
    tt,
    x,
    phi = new_phi
  ))

  phi_str <- stringr::str_c('phi.', 1:7, collapse = ' + ')
  fm <- stringr::str_c('x ~ ', phi_str,
                       ' + (', phi_str, ' | id)')
  lme_fit <- lme4::lmer(as.formula(fm), data = dat_10)

  sc0 <- fpc.0$score_ds
  # ys10_ds <- merge(sc1_on_0, dat) %>% select(-x, -tt) %>% unique
  ys0_ds <- merge(sc0, dat) %>% select(-x, -tt) %>% unique
  xis <- ys0_ds %>% select(contains('xi')) ## scores among untreated
  n <- nrow(ys0_ds)

  ## now we estimate the relationship between the scores and the outcome
  ## AMONG THE UNTREATED - when we plug the projected scores for the treated
  ## into this function, we will be estimating the predicted outcome AMONG THE
  ## UNTREATED if they had longitudinal surrogates like the treated

  yhat_fns <- lapply(xis, function(xx) {
    with(ys0_ds, KernSmooth::locpoly(xx, y, bandwidth = 4*sd(pnorm(xx))*n^(-1/5)))
  })
  # browser()
  yhat_fn <- function(new_xi) {
    id <- new_xi$id
    nx <- new_xi %>% select(-id)
    k <- ncol(nx)
    yh <- matrix(0, n, k)
    for (j in 1:k) {
      yh[,j] <- approx(yhat_fns[[j]]$x, yhat_fns[[j]]$y, xout=nx[,j], rule = 2)$y
    }
    # browser()
    tibble(id = id, yh = rowSums(yh))
  }

  yhat <- yhat_fn(sc1_on_0)

  list(yhat = yhat,
       fpc.0 = fpc.0, fpc.1 = fpc.1)
}


lsa_new <- function(dat, msK = NULL, fpc_fit = NULL) {
  if (is.null(fpc_fit)) {
    fpc <- fpca(dat, 'x', 'tt', 'id', options = list(
      nRegGrid = 150,
      methodSelectK = ifelse(is.null(msK), 'AIC', msK),
      maxK = 60,
      useBinnedCov = FALSE))
  } else fpc <- fpc_fit
  # browser()
  id_unique_dat <- dat %>% select(id, a, y) %>% unique %>%
    mutate(id = as.character(id))
  yh_ds <- fpc$yh_ds %>% inner_join(id_unique_dat)
  yh_x <- yh_ds %>% select(contains('yhat')) %>% as.matrix
  yh_w <- yh_ds$a
  yh_y <- yh_ds$y
  rf_fit <- grf::causal_forest(X = yh_x, W = yh_w, Y = yh_y, seed = pn)
  delta <- grf::estimate_average_effect(rf_fit, target.sample = 'control')
  delta
}
