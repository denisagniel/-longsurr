hiv_sim_fn <- function(s, mean_fn) {
  if (mean_fn == 'kernel') {
    sigma_1 <- k_sigma_1
    sigma_0 <- k_sigma_0
    
    trt_guys <- trt_guys %>%
      mutate(y_1 = predict(kernel_fit_1),
             y_0 = predict(kernel_fit_0, fdX_t))
    control_guys <- control_guys %>%
      mutate(y_0 = predict(kernel_fit_0),
             y_1 = predict(kernel_fit_1, fdX_c))
    
    sim_pool <- smoothed_data %>%
      full_join(trt_guys %>%
                  full_join(control_guys))
    
    sim_id_data <- sim_pool %>%
      dplyr::select(id, a, y_1, y_0) %>%
      unique
    sim_facts <- sim_id_data %>%
      summarise(Delta = mean(y_1[a == 1]) - mean(y_0[a==0]),
                Delta_S = mean(y_1[a==0]) - mean(y_0[a==0]),
                R = 1 - Delta_S/Delta)
  } else if (method == 'gam') {
    sigma_1 <- g_sigma_1
    sigma_0 <- g_sigma_0
    
    trt_guys <- trt_guys %>%
      mutate(y_1 = predict(fgam_fit_1),
             y_0 = predict(fgam_fit_0, newdata = list(X_0 = X_1)))
    control_guys <- control_guys %>%
      mutate(y_0 = predict(fgam_fit_0),
             y_1 = predict(fgam_fit_1, newdata = list(X_1 = X_0)))
    
    sim_pool <- smoothed_data %>%
      full_join(trt_guys %>%
                  full_join(control_guys))
    sim_id_data <- sim_pool %>%
      dplyr::select(id, a, y_1, y_0) %>%
      unique
    sim_facts <- sim_id_data %>%
      summarise(Delta = mean(y_1[a == 1]) - mean(y_0[a==0]),
                Delta_S = mean(y_1[a==0]) - mean(y_0[a==0]),
                R = 1 - Delta_S/Delta)
  } else if (method == 'linear') {
    sigma_1 <- l_sigma_1
    sigma_0 <- l_sigma_0
    
    trt_guys <- trt_guys %>%
      mutate(y_1 = predict(lin_1),
             y_0 = predict(lin_0, newdata = list(X_0 = X_1)))
    control_guys <- control_guys %>%
      mutate(y_0 = predict(lin_0),
             y_1 = predict(lin_1, newdata = list(X_1 = X_0)))
    
    sim_pool <- smoothed_data %>%
      full_join(trt_guys %>%
                  full_join(control_guys))
    sim_id_data <- sim_pool %>%
      dplyr::select(id, a, y_1, y_0) %>%
      unique
    sim_facts <- sim_id_data %>%
      summarise(Delta = mean(y_1[a == 1]) - mean(y_0[a==0]),
                Delta_S = mean(y_1[a==0]) - mean(y_0[a==0]),
                R = 1 - Delta_S/Delta)
  }
  
  sim_x_data <- sim_pool %>%
    dplyr::select(id, tt, xh)
  
  set.seed(s)
  
  ep_1 <- rnorm(nrow(sim_id_data), sd = sigma_1)
  ep_0 <- rnorm(nrow(sim_id_data), sd = sigma_0)
  rerandomize <- sim_id_data %>%
    sample_frac(replace = TRUE) %>%
    mutate(
      # new_a = sample(sim_id_data$a),
      new_y = if_else(a == 1, y_1 + ep_1, y_0 + ep_0),
      new_id = as.numeric(as.factor(sample(sim_id_data$id)))) %>%
    inner_join(analysis_data)
  
  rerand_ds <- rerandomize %>%
    inner_join(sim_x_data) %>%
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
                          verbose = FALSE))
  
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
                dplyr::select(id, tt, x)) %>%
    inner_join(rerand_ds %>%
                 dplyr::select(id, y) %>%
                 unique)
  
  rrX_1 <- rr_x1[,colnames(rr_x1) %in% colnames(rr_x0)]
  rrX_0 <- rr_x0[,colnames(rr_x0) %in% colnames(rr_x1)]
  
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
  out_res %>%
    mutate(Delta_0 = sim_facts$Delta,
           Delta_S0 = sim_facts$Delta_S,
           R_0 = sim_facts$R)
}

