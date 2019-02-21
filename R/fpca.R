fpca <- function(ds, ycol, tcol, idcol, options = list(plot = TRUE, methodBwCov = 'GCV')) {
  # browser()
  list_data <- dataframe_to_list(ds, ycol, tcol, idcol)
  fpca_result <- with(list_data, fdapace::FPCA(y_list, t_list, options))
  score_ds <- data.frame(fpca_result$xiEst)
  colnames(score_ds) <- stringr::str_c('xi.', 1:ncol(score_ds))
  score_ds$id <- names(list_data$t_list)
  xi <- fpca_result$xiEst
  phi <- fpca_result$phi
  mu <- fpca_result$mu
  yh <- t(t(xi %*% t(phi)) + mu)
  yh_ds <- data.frame(yh)
  colnames(yh_ds) <- stringr::str_c('yhat.', 1:ncol(yh_ds))
  yh_ds$id <- names(list_data$t_list)
  # phi_pl <- plot_phi(fpca_result)
  # mu_pl <- plot_mu(fpca_result)
  phi_hat <- approx_phi(fpca_result$phi, fpca_result$workGrid, ds %>% select_(tcol) %>% unlist)
  muhat <- approx_mu(fpca_result$mu, fpca_result$workGrid, ds %>% select_(tcol) %>% unlist)
  phi_ds <- data.frame(id = ds %>% select_(idcol), tt = ds %>% select_(tcol), mu = muhat, phi = phi_hat)

  Lambda <- diag(fpca_result$lambda)
  ids <- unique(ds$id)
  nn <- length(ids)
  Sigma <- fpca_result$smoothedCov
  diag(Sigma) <- ifelse(diag(Sigma) < 0, 1e-6, diag(Sigma))
  # var_l <- list()
  # for (i in 1:nn) {
  #   id_i <- ids[i]
  #   ds_i <- ds %>% filter(id == id_i)
  #   sigma2 <- fpca_result$sigma2
  #   Sigma_i <- approx_cov(Sigma,
  #                         fpca_result$workGrid,
  #                         ds_i %>% select_(tcol) %>% unlist) +
  #     diag(nrow(ds_i))*sigma2
  #   phi_i <- phi_ds %>%
  #     filter(id == id_i) %>%
  #     select(contains('phi')) %>%
  #     as.matrix
  #   H <- phi_i %*% Lambda
  #   Omega <- Lambda - t(H) %*% solve(Sigma_i) %*% H
  #   # yh_var <- phi_i %*% Omega %*% t(phi_i)
  #   var_l[[i]] <- data.frame(id = as.character(i), t(diag(Omega)))
  #   colnames(var_l[[i]])[-1] <- paste0('xi.var.', 1:ncol(phi_i))
  # }
  # browser()
  # var_ds <- bind_rows(var_l)
  score_ds <- score_ds #%>% inner_join(var_ds)


  list(score_ds = score_ds, phi_ds = phi_ds, yh_ds = yh_ds, fpca_result = fpca_result)
}

approx_phi <- function(phi_mat, t_grid, new_t) {
  apply(phi_mat, 2, function(p) {
    approx(t_grid, p, xout = new_t, rule = 2)$y
  })
}

approx_mu <- function(mu, t_grid, new_t) {
  approx(t_grid, mu, xout = new_t, rule = 2)$y
}

approx_cov <- function(cov_surface, t_grid, new_t) {

  vars <- diag(cov_surface)
  var_approx <- approx(t_grid, vars, xout = new_t, rule = 2)$y

  one.approx <- apply(cov_surface, 2, function(s) {
    approx(t_grid, s, xout = new_t, rule = 2)$y
  })
  cov_approx <- apply(one.approx, 1, function(s) {
    approx(t_grid, s, xout = new_t, rule = 2)$y
  })
  diag(cov_approx) <- var_approx
  cov_approx
}

# project_scores <- function(y, tt, fpca_project, fpc_0) {
#   nwg <- length(fpca_project$workGrid)
#   projecting_phi <- fpca_project$phi
#   mu_0 <- fpc_0$mu
#   projecting_lambda <- fpca_project$lambda
#   sigma_0 <- fpc_0$smoothedCov
#
#
#   phi_i <- approx_phi(projecting_phi,
#                       fpca_project$workGrid,
#                       tt)
#   mu_i <- approx_mu(mu_0,
#                     fpc_0$workGrid,
#                     tt)
#   G_i <- approx_cov(sigma_0,
#                     fpc_0$workGrid,
#                     tt)
#   sigma_i <- G_i + fpc_0$sigma2*diag(length(diag(G_i)))
#   ymu <- y - mu_i
#   projecting_lambda * c(t(phi_i) %*%
#     solve(sigma_i) %*% ymu)
# }

project_scores <- function(y, tt, fpc_project) {
  nwg <- length(fpc_project$workGrid)
  projecting_phi <- fpc_project$phi
  mu_0 <- fpc_project$mu
  projecting_lambda <- fpc_project$lambda
  sigma_0 <- fpc_project$smoothedCov


  phi_i <- approx_phi(projecting_phi,
                      fpc_project$workGrid,
                      tt)
  mu_i <- approx_mu(mu_0,
                    fpc_project$workGrid,
                    tt)
  G_i <- approx_cov(sigma_0,
                    fpc_project$workGrid,
                    tt)
  sigma_i <- G_i + fpc_project$sigma2*diag(length(diag(G_i)))
  ymu <- y - mu_i
  projecting_lambda * c(t(phi_i) %*%
    solve(sigma_i) %*% ymu)
}


dataframe_to_list <- function(ds, ycol, tcol, idcol) {
  id <- dplyr::select_(ds, idcol)
  y_list <- plyr::dlply(ds, .variables = idcol, function(x) {
    unlist(dplyr::select_(x, ycol))
  })
  t_list <- plyr::dlply(ds, .variables = idcol, function(x) {
    unlist(dplyr::select_(x, tcol))
  })
  list(y_list = y_list, t_list = t_list, id = id)
}

# xhat_1 <- function(tt, xi, mu)
