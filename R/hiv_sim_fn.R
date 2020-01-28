hiv_sim_fn <- function(s, mean_fn) {
  if (mean_fn == 'kernel') {
    sigma_1 <- k_sigma_1
    sigma_0 <- k_sigma_0
    
    sim_pool <- smoothed_data %>%
      full_join(k_sim_id_data)
    sim_id_data <- k_sim_id_data
    sim_facts <- k_sim_facts
  } else if (mean_fn == 'gam') {
    sigma_1 <- g_sigma_1
    sigma_0 <- g_sigma_0
    
    sim_pool <- smoothed_data %>%
      full_join(g_sim_id_data)
    sim_id_data <- g_sim_id_data
    sim_facts <- g_sim_facts
  } else if (mean_fn == 'linear') {
    sigma_1 <- l_sigma_1
    sigma_0 <- l_sigma_0
    
    sim_pool <- smoothed_data %>%
      full_join(l_sim_id_data)
    sim_id_data <- l_sim_id_data
    sim_facts <- l_sim_facts
  }
  
  sim_x_data <- sim_pool %>%
    dplyr::select(id, tt, xh) %>%
    filter(!is.na(xh))
  
  set.seed(s)
  
  ep_1 <- rnorm(nrow(sim_id_data), sd = sigma_1)
  ep_0 <- rnorm(nrow(sim_id_data), sd = sigma_0)
  rerandomize <- sim_id_data %>%
    sample_frac(replace = TRUE) %>%
    mutate(
      # new_a = sample(sim_id_data$a),
      new_y = if_else(a == 1, y_1 + ep_1, y_0 + ep_0),
      new_id = as.numeric(as.factor(sample(sim_id_data$id))))
  
  rerand_ds <- rerandomize %>%
    inner_join(sim_x_data) %>%
    group_by(new_id) %>%
    sample_n(6) %>%
    ungroup %>%
    transmute(id = new_id,
              a = a,
              y = new_y,
              tt = tt,
              x = xh) %>%
    unique %>%
    arrange(id, tt)
  
  
  c(rr_x1, rr_x0) %<-%
    presmooth_data(obs_data = rerand_ds, 
                   options = 
                     list(plot = FALSE, 
                          # methodBwCov = 'GCV',
                          methodBwMu = 'CV',
                          methodSelectK = 'AIC',
                          useBinnedCov = FALSE,
                          verbose = FALSE,
                          usergrid = FALSE))
  
  
  
  rr_smthds <- as_tibble(rr_x1, rownames = 'id') %>%
    mutate(a = 1) %>%
    full_join(
      as_tibble(rr_x0, rownames = 'id') %>%
        mutate(a = 0)
    ) %>%
    gather(tt, xh, -id, -a) %>%
    mutate(tt = as.numeric(tt),
           id = as.numeric(id)) %>%
    left_join(rerand_ds %>%
                dplyr::select(id, tt, x, a)) %>%
    inner_join(rerand_ds %>%
                 dplyr::select(id, y) %>%
                 unique)
  
  new_grid <- c(3,  7, 11, 15, 19, 23, 27, 31, 35, 39, 43)
  rrX_1 <- rr_x1[,new_grid]
  rrX_0 <- rr_x0[,new_grid]
  
  rrfdX_t <- fdata(rrX_1)
  rrfdX_c <- fdata(rrX_0)
  
  rrtrt_guys <- rerand_ds %>%
    filter(a == 1) %>%
    dplyr::select(id, y) %>%
    unique
  rry_t <- rrtrt_guys %>%
    pull(y)
  rrcontrol_guys <- rerand_ds %>%
    filter(a == 0) %>%
    dplyr::select(id, y) %>%
    unique
  rry_c <- rrcontrol_guys %>%
    pull(y)
  
  out_res <- list(kernel = estimate_surrogate_value(y_t = rry_t, y_c = rry_c, X_t = rrX_1, X_c = rrX_0, method = 'kernel'),
                  gam = estimate_surrogate_value(y_t = rry_t, y_c = rry_c, X_t = rrX_1, X_c = rrX_0, method = 'gam', verbose = TRUE),
                  linear = estimate_surrogate_value(y_t = rry_t, y_c = rry_c, X_t = rrX_1, X_c = rrX_0, method = 'linear', verbose = TRUE)
  ) %>% bind_rows(.id = 'method')
  if (any(is.na(out_res$Deltahat_S))) browser()
  out_res %>%
    mutate(Delta_0 = sim_facts$Delta,
           Delta_S0 = sim_facts$Delta_S,
           R_0 = sim_facts$R)
}

