library(macpan2)
library(tidyverse)

create_betamac_df <- function(param_df, pad_steps, numrbf, terms = 3){
  SINE_MAT <- sin(2 * pi / wkyear * (1:pad_steps) %o% (1:terms))
  COSINE_MAT <- cos(2 * pi / wkyear * (1:pad_steps) %o% (1:terms))

  b0 <- param_df %>% filter(mat == "b0") %>% pull(current)
  X <- rbf(pad_steps, numrbf)
  c <- param_df %>% filter(mat == "c") %>% pull(current)
  sin_coeffs_mat <- param_df %>%
    filter(mat == "sin_coeffs_mat") %>%
    pull(current) %>%
    matrix(ncol = 3)
  cos_coeffs_mat <- param_df %>%
    filter(mat == "cos_coeffs_mat") %>%
    pull(current) %>%
    matrix(ncol = 3)

  beta_trend <- b0 + X %*% c
  beta <- exp(beta_trend) * (1 + (rowSums(SINE_MAT * (X %*% sin_coeffs_mat)) +
                                  rowSums(COSINE_MAT * (X %*% cos_coeffs_mat))))

  return(data.frame(beta, beta_trend))
}
