#' Pre-smooth sparse longitudinal data
#' 
#' @param obs_data data.frame or tibble containing the observed data, with columns \code{id} identifying the individual measured, \code{tt} identifying the time of the observation, \code{x} the value of the surrogate at time \code{tt}, and \code{a} indicating 1 for treatment arm and 0 for control arm.
#' @param ... additional arguments to pass to the function PCA function
#' 
#' @return list containing matrices \code{X_t} and \code{X_c}, which are the smoothed surrogate values for the treated and control groups, respectively, for use in downstream analyses 
#' 
#' @examples
#' library(dplyr)
#' library(longsurr)
#' obs_data <- 
#' generate_discontinuous_data(n = 50, n_i = 5, delta_s = 0.5, 
#' k = 1, s_y = 0.1, s_x = 0.1)$obs_ds
#' 
#' head(obs_data)
#' presmooth_X <- presmooth_data(obs_data)
#' 
#' wide_ds <- full_data %>% 
#' dplyr::select(id, a, tt, x, y) %>%
#' tidyr::spread(tt, x) 
#' 
#' y_t <- wide_ds %>%
#' filter(a == 1) %>%
#' pull(y)
#' y_c <- wide_ds %>%
#' filter(a == 0) %>%
#' pull(y)
#' X_t <- presmooth_X$X_t
#' X_c <- presmooth_X$X_c
#' 
#' estimate_surrogate_value(y_t = y_t, y_c = y_c, 
#' X_t = X_t, X_c = X_c, method = 'linear')
#' 
#' @import dplyr
#' @export
#' 

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

  list(X_t = trt_xhat_wide, X_c = ctrl_xhat_wide)
}
