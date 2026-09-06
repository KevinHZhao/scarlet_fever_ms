library(macpan2)
library(tidyverse)
library(parallel)
options(
  macpan2_tmb_adfun_args = list(
    inner.control = list(maxit = 100000)
    #, intern = TRUE
  ),
  warn = 1,
  error = function() {
    traceback(2)
    options(error = NULL)
    stop("exiting after script error")
  }
)

source("tmb.R")
wkyear <- 365.25/7
load("SF.RData")

## Read scans of RGWR cases from chatgpt, the date column is the END date of the week,
## to make it consistent I subtract seven days from the year column
## Assuming cases represents number of NEW infecteds each week
RGWR_cases <- read.csv("RGWR_London_scarlet_fever_1901-1954.csv") %>%
  mutate(date = ymd(date) - days(7))

## Test consistency between the deaths from gpt scans and from digitized scans after interpolating missing:
## full_series %>% filter(!is.na(deaths), deaths != interpolated.deaths) %>% View()

full_series <- normalized_scarlet_fever_data %>%
  mutate(birth.trend = approx(x = births$numdate, y = births$birth.trend, xout = numdate)$y,
         pop = approx(x = births$numdate, y = births$pop, xout = numdate)$y,
         date = ymd(date)) %>%
  filter(numdate > 1842.01) %>%
  select(numdate, date, interpolated.deaths, birth.trend, acm_trend, pop) %>%
  left_join(RGWR_cases, by = "date") %>%
  mutate(cases = approx(x = numdate, y = cases, xout = numdate)$y) # missing two weeks of case data in 1939, interpolating it

front_pad <- 479 # Based on visually making early mortality padding look right
end_pad <- 508 # Based on visually making late cases look right
## plots to check these look good:
## plot(c(head(full_series$interpolated.deaths, n = front_pad), full_series$interpolated.deaths), type = "l"); lines(head(full_series$interpolated.deaths, n = front_pad), col = "red")
## plot(c(RGWR_cases$cases, tail(RGWR_cases$cases, n = end_pad)), type = "l"); lines(x = 1:end_pad + length(RGWR_cases$cases), tail(RGWR_cases$cases, n = end_pad), col = "red")
steps <- nrow(full_series) + front_pad + end_pad
sundays <- seq(from = 0, to = steps, by = 1) ## Weekly time step, USING 0:steps CAUSES A BOMB...

first_case_data <- which(!is.na(full_series$cases))[[1]] + front_pad # which time step we will have first case data

CFP_max = c(0.5, 1, 2) * 0.025
CFP_min = c(0.5, 1, 2) * 0.01
CFP_rate = c(0.5, 1, 2) * 0.01
CFP_mid = c(0.9, 1, 1.1) * steps/2
CFP_parms <- expand_grid(CFP_max, CFP_min, CFP_rate, CFP_mid)

numrbf <- 64
terms <- 3
X <- rbf(steps, numrbf)
Xsparse <- macpan2:::sparse_matrix_notation(X, tol = 1e-4)

SINE_MAT = sin(2 * pi / wkyear * (1:steps) %o% (1:terms))
COSINE_MAT = cos(2 * pi / wkyear * (1:steps) %o% (1:terms))

N0 = full_series$pop[1]
# --------
simulator <-
  mp_tmb_update(
    model = spec,
    default = list(
      births = c(rep(full_series$birth.trend[1], front_pad),
                 full_series$birth.trend,
                 rep(full_series$birth.trend[nrow(full_series)], end_pad)
      )
      , pop = c(rep(full_series$pop[1], front_pad),
                full_series$pop,
                rep(full_series$pop[nrow(full_series)], end_pad)
      )
      , N = N0)
    , must_save = c("state", "total_inflow", "beta", "CFP", "infection")
  ) %>%
  mp_tmb_insert(
    phase = "before",
    at = 1L,
    expressions = list(
      sin_weights_mat ~ cbind(group_sums(Xvals * sin_coeffs_mat[Xcols, 0], Xrows, sin_weights_mat[0:(steps - 1), 0]),
                              group_sums(Xvals * sin_coeffs_mat[Xcols, 1], Xrows, sin_weights_mat[0:(steps - 1), 1]),
                              group_sums(Xvals * sin_coeffs_mat[Xcols, 2], Xrows, sin_weights_mat[0:(steps - 1), 2])),
      cos_weights_mat ~ cbind(group_sums(Xvals * cos_coeffs_mat[Xcols, 0], Xrows, cos_weights_mat[0:(steps - 1), 0]),
                              group_sums(Xvals * cos_coeffs_mat[Xcols, 1], Xrows, cos_weights_mat[0:(steps - 1), 1]),
                              group_sums(Xvals * cos_coeffs_mat[Xcols, 2], Xrows, cos_weights_mat[0:(steps - 1), 2])),
      beta_trend ~ b0 + group_sums(Xvals * c[Xcols], Xrows, beta_trend),
      beta_values ~ exp(beta_trend) * (1 + (row_sums(SINE_MAT * sin_weights_mat) + row_sums(COSINE_MAT * cos_weights_mat)))
    ),
    default = list(
      Xvals = Xsparse$values,
      Xrows = Xsparse$row_index,
      Xcols = Xsparse$col_index,
      sin_coeffs_mat = matrix(0, ncol = terms, nrow = numrbf),
      cos_coeffs_mat = matrix(0, ncol = terms, nrow = numrbf),
      sin_weights_mat = matrix(NA_real_, ncol = terms, nrow = steps),
      cos_weights_mat = matrix(NA_real_, ncol = terms, nrow = steps),
      c = rep(0, numrbf),
      b0 = 0,
      beta_trend = rep(NA_real_, steps),
      SINE_MAT = SINE_MAT,
      COSINE_MAT = COSINE_MAT,
      steps = steps,
      first_case_data = first_case_data
    )
  ) %>%
  mp_tmb_insert(
    phase = "after"
    , at = Inf
    , expressions = list(
      weekly ~ rbind_time(D, sundays, 0)
      , weekly ~ block(weekly, 1, 0, n, 1) - block(weekly, 0, 0, n, 1)
      , infects ~ rbind_time(infection, 1:steps, first_case_data + 1)
      # , infects ~ block(infects, first_case_data, 0, steps - first_case_data - 1, 1)
    )
    , integers = nlist(sundays, n = length(sundays) - 1)
  ) %>%
  mp_rk4() %>%
  mp_simulator(time_steps = steps, outputs = c("S", "I", "R", "D", "N", "model_pop", "mu_S", "mu_I", "mu_R", "weekly", "infection", "infects"))

simulator$add$matrices(
  Sp = 1/6.5
  , Ip = 0.0002
  , CFP_max = 0.025
  , CFP_min = 0.01
  , CFP_rate = 0.01
  , CFP_mid = steps/2
  , .mats_to_return = c("Si", "Ii", "beta_trend")
  , log_disp = 1
  , log_disp_cases = 1
)
simulator$insert$expressions(
  S ~ Sp * N,
  I ~ Ip * N,
  R ~ N - S - I,
  .phase = "before"
)
simulator$insert$expressions(
  beta ~ beta_values[time_step(0) - 1],
  # CFP_pointer ~ time_group(CFP_pointer, CFP_changepoints),
  CFP ~ CFP_min + (CFP_max-CFP_min)/(1+exp(CFP_rate*(time_step(0)-CFP_mid))),
  .phase = "during"
)

simulator$add$matrices(
  D_obs = c(full_series$interpolated.deaths[1:front_pad],
            full_series$interpolated.deaths,
            tail(full_series$interpolated.deaths, n = end_pad)),
  cases_obs = na.omit(c(
    full_series$cases,
    tail(RGWR_cases$cases, n = end_pad)
  )),
  log_lik = empty_matrix,
  .mats_to_save = c("weekly", "log_lik", "infection"),
  .mats_to_return = c("D_obs", "cases_obs", "weekly", "log_lik", "infection")
)
## param for std of incidence and penalty params for RBF coeffs (std on a Gaussian random eff)
simulator$add$matrices(
  rbf_beta_sd = 1,
  rbf_fourrier_sd = 1
)

## Objective function
simulator$insert$expressions(
  log_lik ~ sum(
    dnbinom(
      D_obs,  ## observed values
      clamp(weekly), ## Round up to avoid zeros, IS THERE A BETTER WAY? Try clamping
      exp(log_disp)
    )
  ) +
    sum(
      dnbinom(
        cases_obs,
        clamp(infects),
        exp(log_disp_cases)
      )
    ) +
    sum(
      dnorm(c, 0, rbf_beta_sd)
    ) +
    sum(
      dnorm(sin_coeffs_mat, 0, rbf_fourrier_sd)
    ) +
    sum(
      dnorm(cos_coeffs_mat, 0, rbf_fourrier_sd)
    )
  , .at = Inf
  , .phase = "after"
)

## Step 4: specify the objective function (very often
##         this will be minus the sum of the log likelihoods).
simulator$replace$obj_fn(~ -log_lik)

## Step 5: declare (and maybe transform) parameters to be optimized

simulator$add$transformations(
  Logit("Sp"),
  Logit("Ip"),
  Logit("CFP_min"),
  Logit("CFP_max"),
  Log("CFP_rate"),
  Log("rbf_beta_sd"),
  Log("rbf_fourrier_sd"),
  Log("CFP_mid")
  # Logit("CFP_values")
)


simulator_fun <- function(CFP_max, CFP_min, CFP_rate, CFP_mid){
  simulator$replace$params(
    default = c(
      qlogis(1/6.5), ## From AndeMay91
      qlogis(0.0002),
      qlogis(0.01),
      qlogis(0.025),
      1,
      1,
      log(steps/2),
      log(0.01), ## Set initial value of CFP logistic curve rate parameter to 0.01
      0,
      0,
      0,
      rep(0, (2*terms + 1) * numrbf)
    ),
    mat = c(
      "logit_Sp",
      "logit_Ip",
      "logit_CFP_min",
      "logit_CFP_max",
      "log_disp",
      "log_disp_cases",
      "log_CFP_mid",
      "log_CFP_rate",
      "b0",
      "log_rbf_beta_sd",
      "log_rbf_fourrier_sd",
      rep("c", numrbf),
      rep("cos_coeffs_mat", terms * numrbf),
      rep("sin_coeffs_mat", terms * numrbf)
    ),
    col = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            rep(0, numrbf),
            rep(0:(terms-1), each = numrbf),
            rep(0:(terms-1), each = numrbf)),
    row = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0:(numrbf-1),
            rep(0:(numrbf-1), terms),
            rep(0:(numrbf-1), terms))
  )

  simulator$report()
  simulator$print$matrix_dims()

  start_time <- Sys.time()

  out <- simulator$optimize$nlminb(
    control = list(eval.max = 1000000, iter.max = 1000000)
  )

  params <- simulator$current$params_frame()
  results <- simulator$report()
  final <- mp_final(simulator)

  write.csv(params, paste0("output_sensitivity/Params_", CFP_min, "_", CFP_max, "_", CFP_mid, "_", CFP_rate, ".csv"))
  write.csv(results, paste0("output_sensitivity/Results_", CFP_min, "_", CFP_max, "_", CFP_mid, "_", CFP_rate, ".csv"))
  write.csv(final, paste0("output_sensitivity/Final_", CFP_min, "_", CFP_max, "_", CFP_mid, "_", CFP_rate, ".csv"))

  # set.seed(101)
  # nsim <- 1000

  # sdr <- TMB::sdreport(simulator$ad_fun(), getJointPrecision = TRUE)
  # saveRDS(sdr, paste0("sdr_rate_", rate_vals[ID], ".Rds"))

  end_time <- Sys.time()
  print(end_time - start_time)
  return(out)
  #
}

all_res <- mclapply(
  split(CFP_parms, seq_len(nrow(CFP_parms))),
  function(row) {
    simulator_fun(row$CFP_max, row$CFP_min, row$CFP_rate, row$CFP_mid)
  },
  mc.cores = 81
)

saveRDS(all_res, file = "output_sensitivity/all_res.RDS")
