generate_linear_data <- function(n, n_i, k, s_y, s_x, delta, nt = 101) {
  t <- seq(0, 10, length = nt)
  
  ds_1 <- tibble(
    id = rep(1:n, each = nt),
    a = 1,
    u1 = rep(rnorm(n, sd = 5), each = nt),
    u2 = rep(rnorm(n, mean = 1, sd = 1), each = nt),
    v11 = rep(rnorm(n, sd = 2), each = nt),
    v21 = rep(rnorm(n, sd = 2), each = nt),
    v12 = rep(rnorm(n, sd = 1), each = nt),
    v22 = rep(rnorm(n, sd = 1), each = nt),
    v13 = rep(rnorm(n, sd = 0.5), each = nt),
    v23 = rep(rnorm(n, sd = 0.5), each = nt),
    epsilon = rep(rnorm(n, sd = s_y), each =nt),
    tt = rep(t, n)
  ) %>%
    mutate(X = u1 + u2*t^2*50/333 + 
             v11*sin(pi/5*t) + v21*cos(pi/5*t) +
             v12*sin(2*pi/5*t) + v22*cos(2*pi/5*t) +
             v13*sin(3*pi/5*t) + v23*cos(3*pi/5*t),
           x = X + rnorm(n*nt, sd = s_x)) %>%
    mutate(r_x = X*(t/3)^2) %>%
    group_by(id) %>%
    mutate(r = mean(r_x)) %>%
    ungroup %>%
    mutate(y = delta + r + epsilon)
  
  ds_0 <- tibble(
    id = rep(n + 1:n, each = nt),
    a = 0,
    u1 = rep(rnorm(n, sd = 2), each = nt),
    u2 = rep(rnorm(n, mean = 1, sd = 1), each = nt),
    v11 = rep(rnorm(n, sd = 1), each = nt),
    v21 = rep(rnorm(n, sd = 1), each = nt),
    v12 = rep(rnorm(n, sd = 0.5), each = nt),
    v22 = rep(rnorm(n, sd = 0.5), each = nt),
    v13 = rep(rnorm(n, sd = 1/4), each = nt),
    v23 = rep(rnorm(n, sd = 1/4), each = nt),
    epsilon = rep(rnorm(n, sd = s_y), each =nt),
    tt = rep(t, n)
  ) %>%
    mutate(X = u1 + u2*t + 
             v11*sin(pi/5*t) + v21*cos(pi/5*t) +
             v12*sin(2*pi/5*t) + v22*cos(2*pi/5*t) +
             v13*sin(3*pi/5*t) + v23*cos(3*pi/5*t),
           x = X + rnorm(n*nt, sd = s_x)) %>%
    mutate(r_x = X*(t/3)^2) %>%
    group_by(id) %>%
    mutate(r = mean(r_x)) %>%
    ungroup %>%
    mutate(y = r + epsilon)
  
  ds <- ds_1 %>%
    full_join(ds_0)
  
  obs_ds <- ds %>%
    group_by(id) %>%
    sample_n(n_i)
  
  list(full_ds = ds,
       obs_ds = obs_ds)
}
