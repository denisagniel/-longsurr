hiv_cv <- function(s, time_list, all_ids, analysis_data, smoothed_data, trt_xhat_wide) {
  training_data <- all_ids %>% sample_frac(0.8) %>%
    inner_join(analysis_data)
  test_data <- analysis_data %>%
    anti_join(training_data)
  train_smth <- smoothed_data %>%
    filter(id %in% training_data$id)
  test_smth <- smoothed_data %>%
    filter(id %in% test_data$id)
  
  tts_out <- map(time_list, function(tts) {
    X_1 <- trt_xhat_wide[rownames(trt_xhat_wide) %in% training_data$id,tts]
    fdX_t <- fdata(X_1)
    
    trt_guys <- training_data %>%
      filter(a == 1) %>%
      dplyr::select(id, y) %>%
      unique %>%
      arrange(id)
    y_t <- trt_guys %>%
      pull(y)
    
    
    k3_fit_1 <- fregre.np.cv(fdataobj = fdX_t, y = y_t,
                             metric = longsurr:::make_semimetric_pca(3))
    k4_fit_1 <- fregre.np.cv(fdataobj = fdX_t, y = y_t,
                             metric = longsurr:::make_semimetric_pca(4))
    k5_fit_1 <- fregre.np.cv(fdataobj = fdX_t, y = y_t,
                             metric = longsurr:::make_semimetric_pca(5))
    
    lin_1 <- pfr(y_t ~ lf(X_1))
    fgam_fit_1 <- fgam(y_t ~ af(X_1))
    
    
    test_X_1 <- trt_xhat_wide[rownames(trt_xhat_wide) %in% test_data$id,tts]
    test_fdX_t <- fdata(test_X_1)
    test_trt_guys <- test_data %>%
      filter(a == 1) %>%
      dplyr::select(id, y) %>%
      unique %>%
      arrange(id)
    test_y_t <- test_trt_guys %>%
      pull(y)
    
    yhat_k3_1 <- predict(k3_fit_1, new.fdataobj = test_fdX_t)
    yhat_k4_1 <- predict(k4_fit_1, new.fdataobj = test_fdX_t)
    yhat_k5_1 <- predict(k5_fit_1, new.fdataobj = test_fdX_t)
    k3_mse1 <- mean((test_y_t - yhat_k3_1)^2)
    k4_mse1 <- mean((test_y_t - yhat_k4_1)^2)
    k5_mse1 <- mean((test_y_t - yhat_k5_1)^2)
    k3_r1 <- cor(test_y_t, yhat_k3_1)
    k4_r1 <- cor(test_y_t, yhat_k4_1)
    k5_r1 <- cor(test_y_t, yhat_k5_1)
    
    
    yhat_g_1 <- predict(fgam_fit_1, newdata = list(X_1 = test_X_1), type = 'response')
    g_mse1 <- mean((test_y_t - yhat_g_1)^2)
    g_r1 <- cor(test_y_t, yhat_g_1)
    
    yhat_l_1 <- predict(lin_1, newdata = list(X_1 = test_X_1), type = 'response')
    l_mse1 <- mean((test_y_t - yhat_l_1)^2)
    l_r1 <- cor(test_y_t, yhat_l_1)
    
    out_l <- tibble(k3_mse = k3_mse1,
                    k4_mse = k4_mse1,
                    k5_mse = k5_mse1,
                    g_mse = g_mse1,
                    l_mse = l_mse1,
                    k3_r = k3_r1,
                    k4_r = k4_r1,
                    k5_r = k5_r1,
                    g_r = g_r1,
                    l_r = l_r1,
                    tt = list(tts))
    out_l
  })
  
  tts_out %>% bind_rows(.id = 'tt_n')
}
