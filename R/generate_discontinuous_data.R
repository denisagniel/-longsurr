generate_discontinuous_data <- function(n, n_i, k, s_y = 1, s_x = 1, delta_s, nt = 101) {
  t <- seq(-1, 1, length = nt)
  
  ds_1 <- tibble(
    id = rep(1:n, each = nt),
    a = 1,
    alpha = rep(runif(n, min = -k, max = k), each = nt),
    beta = rep(runif(n, min = -k, max = k), each = nt),
    gamma = rep(runif(n, min = -2*k, max = 2*k), each = nt),
    omega = rep(runif(n)*2*pi, each = nt),
    epsilon = rep(rnorm(n, sd = s_y), each =nt),
    tt = rep(t, n)
  ) %>%
    mutate(X = gamma*sin(omega*t) + (alpha + 2*pi)*t + beta,
           x = X + rnorm(n*nt, sd = s_x)) %>%
    mutate(mu_t = 2*pi*tt) %>%
    mutate(r_x = 10*(X > mu_t + 1)*(1-cos(pi*tt))) %>%
    group_by(id) %>%
    mutate(r = mean(r_x)) %>%
    ungroup %>%
    mutate(y = delta_s + r + epsilon)
  
  ds_0 <- tibble(
    id = rep(n + 1:n, each = nt),
    a = 0,
    alpha = rep(runif(n, min = -k/4, max = 3*k/4), each = nt),
    beta = rep(runif(n, min = -k, max = k), each = nt),
    gamma = rep(runif(n, min = -2*k/4, max = 2*3*k/4), each = nt),
    omega = rep(runif(n)*2*pi, each = nt),
    epsilon = rep(rnorm(n, sd = s_y), each =nt),
    tt = rep(t, n)
  ) %>%
    mutate(X = gamma*sin(omega*t) + (alpha + 2*pi)*t + beta,
           x = X + rnorm(n*nt, sd = s_x)) %>%
    mutate(mu_t = 2*pi*tt) %>%
    mutate(r_x = 10*abs(X > mu_t + 1)*(1-cos(pi*tt))) %>%
    group_by(id) %>%
    mutate(r = mean(r_x)) %>%
    ungroup %>%
    mutate(y = r + epsilon)
  
  ds <- ds_1 %>%
    full_join(ds_0)
  
  obs_ds <- ds %>%
    group_by(id) %>%
    sample_n(n_i) %>%
    ungroup
  
  list(full_ds = ds,
       obs_ds = obs_ds)
}
