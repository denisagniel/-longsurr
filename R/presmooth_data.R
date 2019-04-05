presmooth_data <- function(obs_data, ...) {
  # browser()
  treatment_arm <- obs_data %>%
    filter(a == 1) %>%
    arrange(id, tt)
  control_arm <- obs_data %>%
    filter(a == 0) %>%
    arrange(id, tt)
  n_trt <- treatment_arm %>%
    count(id) %>%
    nrow
  n_ctrl <- control_arm %>%
    count(id) %>%
    nrow
  
  trt_fpc_fit <- fpca(ds = treatment_arm, ycol = 'x', tcol = 'tt', idcol = 'id', ...)
  ctrl_fpc_fit <- fpca(ds = control_arm, ycol = 'x', tcol = 'tt', idcol = 'id', ...)
# browser()
  trt_xhat <- trt_fpc_fit$yh_ds %>%
    gather(tp, X, -id) %>%
    mutate(id = as.integer(id),
           tt = rep(seq(-1, 1, length = 51), each = n_trt),
           type = 'estimated')

  ctrl_xhat <- ctrl_fpc_fit$yh_ds %>%
    gather(tp, X, -id) %>%
    mutate(id = as.integer(id),
           tt = rep(seq(-1, 1, length = 51), each = n_ctrl),
           type = 'estimated')
  # browser()
  trt_xhat_wide <- trt_xhat %>%
    dplyr::select(-tp, -type) %>%
    spread(tt, X) %>%
    dplyr::select(-id) %>%
    as.matrix
  colnames(trt_xhat_wide) <- colnames(trt_xhat_wide)
  rownames(trt_xhat_wide) <- trt_xhat$id %>% unique

  ctrl_xhat_wide <- ctrl_xhat %>%
    dplyr::select(-tp, -type) %>%
    spread(tt, X) %>%
    dplyr::select(-id) %>%
    as.matrix
  colnames(ctrl_xhat_wide) <- colnames(ctrl_xhat_wide)
  rownames(ctrl_xhat_wide) <- ctrl_xhat$id %>% unique

  list(trt_xhat_wide, ctrl_xhat_wide, trt_fpc_fit$score_ds, ctrl_fpc_fit$score_ds)
}
