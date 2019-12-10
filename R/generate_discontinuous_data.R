#' Simulate nonsmooth data
#' 
#' @param n sample size
#' @param n_i number of observations per individual
#' @param k tuning parameter which controls the distribution of the parameters that determine the surrogate value
#' @param s_y error standard deviation for the outcome
#' @param s_x error standard deviation for the surrogate
#' @param delta_s residual treatment effect, the parameter of interest
#' @param nt number of possible observation times for the longitudinal surrogate
#' 
#' @import dplyr
#' 
#' @export

generate_discontinuous_data <- function(
  n, n_i, k, s_y = 1, s_x = 1, delta_s, nt = 101) {
  t <- seq(-1, 1, length = nt)
  
  ds_1 <- tibble::tibble(
    id = rep(1:n, each = nt),
    a = 1,
    alpha = rep(runif(n, min = -k, max = k), each = nt),
    beta = rep(runif(n, min = -k, max = k), each = nt),
    gamma = rep(runif(n, min = -2*k, max = 2*k), each = nt),
    omega = rep(runif(n)*2*pi, each = nt),
    epsilon = rep(rnorm(n, sd = s_y), each =nt),
    tt = rep(t, n)
  ) 
  ds_1 <- 
    mutate(ds_1, 
           X = gamma*sin(omega*t) + (alpha + 2*pi)*t + beta,
           x = X + rnorm(n*nt, sd = s_x),
           mu_t = 2*pi*tt,
           r_x = 10*(X > mu_t + 1)*(1-cos(pi*tt)))
  ds_1 <- group_by(ds_1, id) 
  ds_1 <- mutate(ds_1, r = mean(r_x)) 
  ds_1 <- ungroup(ds_1)
  ds_1 <- mutate(ds_1, y = delta_s + r + epsilon)
  
  ds_0 <- tibble(
    id = rep(n + 1:n, each = nt),
    a = 0,
    alpha = rep(runif(n, min = -k/4, max = 3*k/4), each = nt),
    beta = rep(runif(n, min = -k, max = k), each = nt),
    gamma = rep(runif(n, min = -2*k/4, max = 2*3*k/4), each = nt),
    omega = rep(runif(n)*2*pi, each = nt),
    epsilon = rep(rnorm(n, sd = s_y), each =nt),
    tt = rep(t, n)
  ) 
  
  ds_0 <- 
    mutate(ds_0, 
           X = gamma*sin(omega*t) + (alpha + 2*pi)*t + beta,
           x = X + rnorm(n*nt, sd = s_x),
           mu_t = 2*pi*tt,
           r_x = 10*(X > mu_t + 1)*(1-cos(pi*tt)))
  ds_0 <- group_by(ds_0, id) 
  ds_0 <- mutate(ds_0, r = mean(r_x)) 
  ds_0 <- ungroup(ds_0)
  ds_0 <- mutate(ds_0, y = r + epsilon)

  
  ds <- full_join(ds_1, ds_0)
  obs_ds <- group_by(ds, id) 
  obs_ds <- sample_n(obs_ds, n_i)
  obs_ds <- ungroup(obs_ds)
  
  list(full_ds = ds,
       obs_ds = obs_ds)
}
