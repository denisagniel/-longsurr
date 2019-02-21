Cx_f<-function(t,s,sig2=1,rho=0.5) { # Matern covariance function of X with nu = 5/2
  d <- abs(outer(t,s,"-"))
  sig2*(1+sqrt(5)*d/rho + 5*d^2/(3*rho^2))*exp(-sqrt(5)*d/rho)
}

matern_cov <- function(t, s, sigma2, p, rho) {
  # browser()
  d <- abs(outer(t,s,"-"))
  last_terms <- lapply(0:p, function(i) {
    factorial(p + i)/factorial(i)/factorial(p-i)*(2*sqrt(2*p+1)*d/rho)^(p-i)
  })
  last_term <- Reduce('+', last_terms)
  sigma2 * (
    exp(-sqrt(2*p + 1)*d/rho)*factorial(p)/factorial(2*p)*last_term
  )
}

estimate_effects <- function(n, w, sx, sy, grid_size, alpha_c = 0, alpha_t = 2, nonlinear = TRUE, fpc_gen = FALSE) {
  # browser()

  grid <- seq(from=0,to=1,length.out = grid_size) ## grid of time points
  # Cx_0 <- Cx_f(grid,grid, sig2 = 0.5) ## covariance of control X
  # Cx_1 <- Cx_f(grid,grid, sig2 = 5) ## covariance of control X
  Cx_0 <- matern_cov(grid, grid, sigma2 = 0.5, p = 0, rho = 2)
  Cx_1 <- matern_cov(grid, grid, sigma2 = 5, p = 2, rho = 0.5)
  # browser()
  phi_0 <- eigen(Cx_0,symmetric = T)$vectors*sqrt(grid_size) ## eigenfunctions
  phi_1 <- eigen(Cx_1,symmetric = T)$vectors*sqrt(grid_size) ## eigenfunctions

  var_eps <- sy ## variance of y
  var_delt_c <- sx ## error variance of x under control
  var_delt_t <- sx ## error variance of x under treatment

  mux_c <- rep(0, grid_size) ## mean of x under control
  mux_t <- rep(w^(1/2), grid_size) ## mean of x under treatment
  beta_c <- sqrt(w) + w*sin(2*pi*grid)
  beta_t <- w^2 + w*sin(2*pi - 2*pi*grid)

  X_c <- mvnfast::rmvn(n,mu = mux_c,sigma = Cx_0)
  X_t <- mvnfast::rmvn(n,mu = mux_t,sigma = Cx_1)
  X_err_c <- X_c + rnorm(n*grid_size,sd = sqrt(var_delt_c))
  X_err_t <- X_t + rnorm(n*grid_size,sd = sqrt(var_delt_t))

  Xi_c <- (X_c-mux_c)%*%phi_0/grid_size
  Xi_t <- (X_t-mux_t)%*%phi_1/grid_size

  eps_c <- rnorm(n,0,sd = sqrt(var_eps))
  eps_t <- rnorm(n,0,sd = sqrt(var_eps))
  if(!fpc_gen) {
    y_c <- c(alpha_c + (X_c + nonlinear*0.1*(X_c - mux_c)^2)%*%beta_c/grid_size + eps_c)
    y_t <- c(alpha_t + (X_t + nonlinear*0.1*(X_t - mux_t)^2) %*% beta_t/grid_size + eps_t)
    y_c_on_t <- c(alpha_t + (X_c + nonlinear*0.1*(X_c - mux_c)^2) %*% beta_t/grid_size)
    y_t_on_c <- c(alpha_c + (X_t + nonlinear*0.1*(X_t - mux_t)^2) %*% beta_c/grid_size)
  } else {
    gamma_c <- c(-1, -2, 0.5, 1, 0.1)
    gamma_t <- c(0, 2, 2, 1, 0.5)
    y_c <- c(alpha_c + Xi_c[,1:5] %*% gamma_c + eps_c)
    y_t <- c(alpha_t + Xi_t[,1:5] %*% gamma_t + eps_t)
    y_c_on_t <- c(alpha_t + Xi_c[,1:5] %*% gamma_t)
    y_t_on_c <- c(alpha_c + Xi_t[,1:5] %*% gamma_c)
  }
  delta <- mean(y_t) - mean(y_c)
  delta_s <- mean(y_c_on_t)
  other_way <- mean(y_t_on_c)
  list(
    delta = delta,
    delta_s = delta_s,
    other_way_delta = other_way
       )
}

generate_data <- function(n, n_i, Cx_0, Cx_1, w, sx, sy, grid_size, alpha_c = 0, alpha_t = 2, nonlinear = TRUE, fpc_gen = FALSE) {
  # browser()

  grid <- seq(from=0,to=1,length.out = grid_size) ## grid of time points
  # Cx_0 <- Cx_f(grid,grid, sig2 = 0.5) ## covariance of control X
  # Cx_1 <- Cx_f(grid,grid, sig2 = 5) ## covariance of control X
  Cx_0 <- matern_cov(grid, grid, sigma2 = 0.5, p = 0, rho = 2)
  Cx_1 <- matern_cov(grid, grid, sigma2 = 5, p = 2, rho = 0.5)
  # browser()
  phi_0 <- eigen(Cx_0,symmetric = T)$vectors*sqrt(grid_size) ## eigenfunctions
  phi_1 <- eigen(Cx_1,symmetric = T)$vectors*sqrt(grid_size) ## eigenfunctions

  var_eps <- sy ## variance of y
  var_delt_c <- sx ## error variance of x under control
  var_delt_t <- sx ## error variance of x under treatment

  mux_c <- rep(0, grid_size) ## mean of x under control
  mux_t <- rep(w^(1/2), grid_size) ## mean of x under treatment
  beta_c <- sqrt(w) + w*sin(2*pi*grid)
  beta_t <- w^2 + w*sin(2*pi - 2*pi*grid)

  X_c <- mvnfast::rmvn(n,mu = mux_c,sigma = Cx_0)
  X_t <- mvnfast::rmvn(n,mu = mux_t,sigma = Cx_1)
  X_err_c <- X_c + rnorm(n*grid_size,sd = sqrt(var_delt_c))
  X_err_t <- X_t + rnorm(n*grid_size,sd = sqrt(var_delt_t))

  Xi_c <- (X_c-mux_c)%*%phi_0/grid_size
  Xi_t <- (X_t-mux_t)%*%phi_1/grid_size

  eps_c <- rnorm(n,0,sd = sqrt(var_eps))
  eps_t <- rnorm(n,0,sd = sqrt(var_eps))
  if(!fpc_gen) {
    y_c <- c(alpha_c + (X_c + nonlinear*0.1*(X_c - mux_c)^2)%*%beta_c/grid_size + eps_c)
    y_t <- c(alpha_t + (X_t + nonlinear*0.1*(X_t - mux_t)^2) %*% beta_t/grid_size + eps_t)
    y_c_on_t <- c(alpha_t + (X_c + nonlinear*0.1*(X_c - mux_c)^2) %*% beta_t/grid_size)
    y_t_on_c <- c(alpha_c + (X_t + nonlinear*0.1*(X_t - mux_t)^2) %*% beta_c/grid_size)
  } else {
    gamma_c <- c(-1, -2, 0.5, 1, 0.1)
    gamma_t <- c(0, 2, 2, 1, 0.5)
    y_c <- c(alpha_c + Xi_c[,1:5] %*% gamma_c + eps_c)
    y_t <- c(alpha_t + Xi_t[,1:5] %*% gamma_t + eps_t)
    y_c_on_t <- c(alpha_t + Xi_c[,1:5] %*% gamma_t)
    y_t_on_c <- c(alpha_c + Xi_t[,1:5] %*% gamma_c)
  }




  full_data <- data.frame(
    id = rep(1:(2*n), each = grid_size),
    tt = rep(grid, 2*n),
    X = c(c(t(X_c)), c(t(X_t))),
    x = c(c(t(X_err_c)), c(t(X_err_t))),
    y = rep(c(y_c, y_t), each = grid_size),
    a = rep(0:1, each = n*grid_size)
  )

  obs_data <- full_data %>%
    group_by(id) %>%
    sample_n(n_i) %>%
    ungroup

  list(y_t = y_t, y_c = y_c, X_t = X_t, X_c = X_c,
         X_err_t = X_err_t,
         X_err_c = X_err_c,
         full_data = full_data,
         obs_data = obs_data,
       y_c_on_t = y_c_on_t,
       y_t_on_c = y_t_on_c,
       Xi_c, Xi_t
         )
}





# generate_mean_shift_data <- function(n, n_i, trt_diff = 0) {
#   tt_0 <- matrix(runif(n*n_i), n, n_i)*10
#   tt_1 <- matrix(runif(n*n_i), n, n_i)*10
#
#   b_00 <- rnorm(n, mean = 5, sd = 4) %>% rep(each = n_i)
#   b_01 <- rnorm(n, mean = 15, sd = 10) %>% rep(each = n_i)
#
#   # browser()
#   x_0 <- b_00 + rnorm(n*n_i, sd = 0.25) %>% matrix(n, n_i)
#   x_1 <- b_01 + rnorm(n*n_i, sd = 0.25) %>% matrix(n, n_i)
#   y_0 <- unique(b_00) + rnorm(n)
#   y_1 <- trt_diff + unique(b_01) + rnorm(n)
#   tibble(
#     id = rep(1:(2*n), each = n_i),
#     tt = c(tt_0, tt_1),
#     a = rep(0:1, each = n*n_i),
#     x = c(x_0, x_1),
#     y = rep(c(y_0, y_1), each = n_i),
#     b0 = c(b_00, b_01)
#   )
# }
#
# generate_spline_data <- function(n, n_i, trt_diff = 0, x_sd = 0.25) {
#   tt_0 <- matrix(runif(n*n_i), n, n_i)*10
#   tt_1 <- matrix(runif(n*n_i), n, n_i)*10
#
#   b_00 <- rnorm(n, mean = 5, sd = 4) %>% rep(each = n_i)
#   b_01 <- rnorm(n, mean = 15, sd = 14) %>% rep(each = n_i)
#   b_10 <- rnorm(n, mean = 0.5, sd = 0.4) %>% rep(each = n_i)
#   b_11 <- rnorm(n, mean = 1.5, sd = 1.4) %>% rep(each = n_i)
#   b_20 <- rnorm(n, mean = 0.25, sd = 0.2) %>% rep(each = n_i)
#   b_21 <- rnorm(n, mean = 0.75, sd = 0.7) %>% rep(each = n_i)
#   # browser()
#   o_basis <- orthogonalsplinebasis::OBasis(
#     knots = orthogonalsplinebasis::expand.knots(seq(0, 10, length = 4))
#     )
#   phi_0 <- evaluate(o_basis, c(tt_0)) ## lazy way to evaluate splines
#   phi_1 <- evaluate(o_basis, c(tt_1))
#
#   # browser()
#   x_0 <- b_00*phi_0[,2] + b_10*phi_0[,3] + b_20*phi_0[,4] + rnorm(n*n_i, sd = x_sd) %>% matrix(n, n_i)
#   x_1 <- b_01*phi_1[,2] + b_11*phi_1[,3] + b_21*phi_1[,4] + rnorm(n*n_i, sd = x_sd) %>% matrix(n, n_i)
#   y_0 <-  unique(b_00 + b_10 + b_20) + rnorm(n)
#   y_1 <- trt_diff + unique(b_01 + b_11 + b_21) + rnorm(n)
#   tibble(
#     id = rep(1:(2*n), each = n_i),
#     tt = c(tt_0, tt_1),
#     a = rep(0:1, each = n*n_i),
#     x = c(x_0, x_1),
#     y = rep(c(y_0, y_1), each = n_i),
#     b0 = c(b_00, b_01),
#     b1 = c(b_10, b_11),
#     b2 = c(b_20, b_21)
#   )
# }
#
# generate_lsa_data <- function(n, n_i, trt_diff = 0) {
#   # browser()
#   tt_0 <- matrix(runif(n*n_i), n, n_i)*10
#   tt_1 <- matrix(runif(n*n_i), n, n_i)*10
#   phi_1 <- function(x) sin(x*2*pi)
#   # phi_2 <- identity
#   # phi_3 <- log
#   phi_2 <- cos
#
#   # phi_0 <- splines::bs(tt_0, df = 5)
#
#  b_00 <- rnorm(n, mean = 10, sd = 4) %>% rep(each = n_i)
#  b_01 <- rnorm(n, mean = 15, sd = 14) %>% rep(each = n_i)
#   b_10 <- rnorm(n, mean = 5, sd = 1) %>% rep(each = n_i)
#   b_11 <- rnorm(n, mean = 7, sd = 5) %>% rep(each = n_i)
#   b_20 <- rnorm(n, mean = 1, sd = 0.5) %>% rep(each = n_i)
#   b_21 <- rnorm(n, mean = 2, sd = 2) %>% rep(each = n_i)
#   b_30 <- rnorm(n, mean = 1, sd = 0.5) %>% rep(each = n_i)
#   b_31 <- rnorm(n, mean = 2, sd = 2) %>% rep(each = n_i)
#   b_40 <- rnorm(n, mean = 1, sd = 0.5) %>% rep(each = n_i)
#   b_41 <- rnorm(n, mean = 2, sd = 2) %>% rep(each = n_i)
#   # browser()
#   x_0 <- tt_0 + sin(tt_0) + b_00 +
#     b_10*phi_1(tt_0) + b_20*phi_2(tt_0) + #b_30*phi_3(tt_0) + b_40*phi_4(tt_0) +
#     rnorm(n*n_i, sd = 0.25) %>% matrix(n, n_i)
#   x_1 <- tt_1 + sin(tt_1) + b_01 +
#     b_11*phi_1(tt_0) + b_21*phi_2(tt_0) + #b_31*phi_3(tt_0) + b_41*phi_4(tt_0) +
#     rnorm(n*n_i, sd = 0.25) %>% matrix(n, n_i)
#   y_0 <- unique(b_00 +
#     2*b_10 + b_20/2 #+ b_30 + b_40
#                 ) + rnorm(n)
#   y_1 <- trt_diff + unique(b_01 +
#     2*b_11 + b_21/2 #+ b_31 + b_41
#                            ) + rnorm(n)
#   tibble(
#     id = rep(1:(2*n), each = n_i),
#     tt = c(tt_0, tt_1),
#     a = rep(0:1, each = n*n_i),
#     x = c(x_0, x_1),
#     y = rep(c(y_0, y_1), each = n_i),
#     b0 = c(b_00, b_01),
#     b1 = c(b_10, b_11),
#     b2 = c(b_20, b_21),
#    b3 = c(b_30, b_31),
#     b4 = c(b_40, b_41)
#   )
# }
#
# # generate_my_data <- function(n = 100) {
# #   n_i <- sample(3:6, size = n)
# #   mu_x <- function(s) s + sin(s)
# #   phi_1 <- function(s) -cos(pi*s/10)/sqrt(5)
# #   phi_2 <- function(s) sin(pi*s/10)/sqrt(5)
# #   lambda1 <- 4
# #   lambda2 <- 1
# #   b_10 <- rnorm(n, 0.5, sd = sqrt(lambda1))
# #   b_11 <- rnorm(n, -0.5, sd = sqrt(lambda1))
# #   b_20 <- rnorm(n, 0.5, sd = sqrt(lambda2))
# #   b_21 <- rnorm(n, -0.5, sd = sqrt(lambda2))
# #   y_1 <- 2*b_11 + b_21/2
# #   y_0 <- 2*b_10 + b_20/2
# # }
