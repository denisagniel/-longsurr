predict_fpca <- function(idd, tt, score_ds, fpca_phi, fpca_mu = NULL) {
  # browser()
  if (is.null(fpca_mu)) fpca_mu <- fpca_phi
  mu_i <- approx_mu(mu = fpca_mu$mu, t_grid = fpca_mu$workGrid, new_t = tt)
  phi_i <- approx_phi(phi_mat = fpca_phi$phi, t_grid = fpca_phi$workGrid, new_t = tt)
  xi_i <- score_ds %>% filter(id == idd) %>% select(-id) %>% unlist
  as.vector(mu_i + phi_i %*% xi_i)
}

