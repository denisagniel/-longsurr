plot_phi <- function(fpca) {
  tt <- fpca$workGrid
  phi <- fpca$phi
  pd <- data.frame(tt = tt,
                   phi = c(phi),
                   phi_n = factor(rep(1:ncol(phi), each = length(tt))))
  qplot(tt, phi, color = phi_n, group = phi_n, geom = 'line', data = pd)
}

plot_mu <- function(fpca) {
  tt <- fpca$workGrid
  mu <- fpca$mu
  pd <- data.frame(tt = tt,
                   mu = mu)
  qplot(tt, mu, geom = 'line', data = pd)
}
