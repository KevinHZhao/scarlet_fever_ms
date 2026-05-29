library(tidyverse)
library(fastbeta)

# Set up ------------------------------------------------------------------

load(file = "../SF.RData")
source("helper_funs/utils.R")
births <- read.csv("../birthrate_1750_1930.csv")
wkyear <- 365.25/7

full_series <- normalized_scarlet_fever_data %>%
  mutate(birth.trend = approx(x = births$numdate, y = births$birth.trend, xout = numdate)$y,
         pop = approx(x = births$numdate, y = births$pop, xout = numdate)$y) %>%
  filter(numdate > 1842.01, numdate < 1930) %>%
  select(numdate, interpolated.deaths, birth.trend, acm_trend, pop)

steps <- nrow(full_series) # Steps in SF series
front_pad <- 479 # SF deaths to repeat at beginning
end_pad <- 5*52 # SF deaths to repeat at end
pad_steps <- steps + front_pad + end_pad # Steps in the final model
numrbf <- 64

CFP_max = c(0.5, 1, 2) * 0.025
CFP_min = c(0.5, 1, 2) * 0.01
CFP_rate = c(0.5, 1, 2) * 0.01
CFP_mid = c(0.9, 1, 1.1) * pad_steps/2
CFP_parms <- expand_grid(CFP_max, CFP_min, CFP_rate, CFP_mid)

# Macpan parameter simplifying --------------------------------------------

read_mp <- function(iCFP_max, iCFP_min, iCFP_rate, iCFP_mid){
  mac_parms <- read.csv(paste0("../output_sensitivity/Params_", iCFP_min, "_", iCFP_max, "_", iCFP_mid, "_", iCFP_rate, ".csv"))
  CFP_min <- mac_parms %>% filter(mat == "logit_CFP_min") %>% pull(current) %>% plogis
  CFP_max <- mac_parms %>% filter(mat == "logit_CFP_max") %>% pull(current) %>% plogis
  CFP_rate <- mac_parms %>% filter(mat == "log_CFP_rate") %>% pull(current) %>% exp
  CFP_mid <- mac_parms %>% filter(mat == "log_CFP_mid") %>% pull(current) %>% exp
  CFP <- CFP_min + (CFP_max-CFP_min)/(1+exp(CFP_rate*((1-8):pad_steps-CFP_mid))) ## CFP curve, with 7 extra time steps

  ## Initial values for S0, I0, R0
  S0_mac <- mac_parms %>% filter(mat=="logit_Sp") %>% pull(current) %>% plogis %>% prod(full_series$pop[1])
  I0_mac <- mac_parms %>% filter(mat=="logit_Ip") %>% pull(current) %>% plogis %>% prod(full_series$pop[1])
  R0_mac <- full_series$pop[1] - S0_mac - I0_mac

  terms <- 3 # Number of Fourier terms

  SINE_MAT <- sin(2 * pi / wkyear * (1:pad_steps) %o% (1:terms))
  COSINE_MAT <- cos(2 * pi / wkyear * (1:pad_steps) %o% (1:terms))

  betamac_df <- create_betamac_df(mac_parms, pad_steps, numrbf)
  ## This df has the fitted betas from macpan for each step (1:pad_steps)

  gamma <- 7/15 # Assume time from infection to death equals latent + infectious period, latent period = 1 day, infectious period = 2 weeks

  mac_final <- read.csv(paste0("../output_sensitivity/Final_", iCFP_min, "_", iCFP_max, "_", iCFP_mid, "_", iCFP_rate, ".csv"))
  mac_results <- read.csv(paste0("../output_sensitivity/Results_", iCFP_min, "_", iCFP_max, "_", iCFP_mid, "_", iCFP_rate, ".csv"))
  mac_muS <- mac_results %>% filter(matrix == "mu_S") %>% pull(value)
  mac_muI <- mac_results %>% filter(matrix == "mu_I") %>% pull(value)
  mac_muR <- mac_results %>% filter(matrix == "mu_R") %>% pull(value)
  mac_outflows <- mac_muS + mac_muI + mac_muR # all-cause mortality used by macpan
  N <- mac_results %>% filter(matrix == "model_pop") %>% pull(value) # Population, equal to full_series$pop

  mac_infection <- mac_results %>% filter(matrix == "infection") %>% pull(value)

  return(betamac_df$beta)
}

full_df <- mclapply(
  split(CFP_parms, seq_len(nrow(CFP_parms))),
  function(row) {
    read_mp(row$CFP_max, row$CFP_min, row$CFP_rate, row$CFP_mid)
  },
  mc.cores = 8
) %>%
  do.call(what = cbind) %>%
  as.data.frame() %>%
  mutate(
    x = c(full_series$numdate[1:front_pad], full_series$numdate, full_series$numdate[nrow(full_series) - end_pad + 1:(end_pad)]),
    week = 1:(pad_steps)
  )

lapply(all_res, FUN = function(x) {x$objective}) |> unlist() |> table()
pdf("compare_fits.pdf")
matplot(full_df %>% select(-x, -week), type = "l", lty = 1, col = 1:ncol(full_df), ylab = "beta", xlab = "week number")
dev.off()
## 5 different plots of beta produced, of which out of the 81 diff param combos
## for which 3 only had one respective param combo, one had two combos, and the
## last had 75 (including the one we used, #41).
