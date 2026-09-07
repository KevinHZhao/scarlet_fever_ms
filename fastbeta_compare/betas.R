library(tidyverse)
library(fastbeta)
library(tikzDevice)

# Set up ------------------------------------------------------------------

load(file = "../CC_code/SF.RData")
source("helper_funs/utils.R")
wkyear <- 365.25/7
RGWR_cases <- read.csv("../CC_code/RGWR_London_scarlet_fever_1901-1954.csv") %>%
  mutate(date = ymd(date) - days(7))

full_series <- normalized_scarlet_fever_data %>%
  mutate(birth.trend = approx(x = births$numdate, y = births$birth.trend, xout = numdate)$y,
         pop = approx(x = births$numdate, y = births$pop, xout = numdate)$y,
         date = ymd(date)) %>%
  filter(numdate > 1842.01) %>%
  select(numdate, date, interpolated.deaths, birth.trend, acm_trend, pop) %>%
  left_join(RGWR_cases, by = "date") %>%
  mutate(cases = approx(x = numdate, y = cases, xout = numdate)$y) # missing two weeks of case data in 1939, interpolating it

steps <- nrow(full_series) # Steps in SF series
front_pad <- 479 # SF deaths to repeat at beginning
end_pad <- 508 # SF deaths to repeat at end
pad_steps <- steps + front_pad + end_pad # Steps in the final model
numrbf <- 64

# Macpan parameter simplifying --------------------------------------------

mac_parms <- read.csv("../CC_code/output/Params.csv")
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

mac_final <- read.csv("../CC_code/output/Final.csv")
mac_results <- read.csv("../CC_code/output/Results.csv")
mac_muS <- mac_results %>% filter(matrix == "mu_S") %>% pull(value)
mac_muI <- mac_results %>% filter(matrix == "mu_I") %>% pull(value)
mac_muR <- mac_results %>% filter(matrix == "mu_R") %>% pull(value)
mac_outflows <- mac_muS + mac_muI + mac_muR # all-cause mortality used by macpan
N <- mac_results %>% filter(matrix == "model_pop") %>% pull(value) # Population, equal to full_series$pop

mac_infection <- mac_results %>% filter(matrix == "infection") %>% pull(value)
# fastbeta setup ----------------------------------------------------------

first_case_data <- which(!is.na(full_series$cases))[[1]] + front_pad # which time step we will have first case data

series_mort <-
  cbind(
    D.obs = full_series$interpolated.deaths[c(1:front_pad, 1:(first_case_data - 1 - front_pad))],
    #Z = mac_infection,
    B = full_series$birth.trend[c(1:front_pad, 1:(first_case_data - 1 - front_pad))],
    mu = (mac_outflows / N)[1:(first_case_data - 1)]
  ) %>%
  ts()

series_case <-
  cbind(
    case.obs = full_series$cases[c((first_case_data - front_pad):nrow(full_series), nrow(full_series) - end_pad + 1:end_pad)],
    #Z = mac_infection,
    B = full_series$birth.trend[c((first_case_data - front_pad):nrow(full_series), nrow(full_series) - end_pad + 1:end_pad)],
    mu = (mac_outflows / N)[first_case_data:pad_steps]
  ) %>%
  ts()

betas_mort <- fastbeta(series_mort,
                  gamma = gamma,
                  init = c(round(S0_mac),
                           round(I0_mac),
                           round(R0_mac)),
                  prob = CFP[1:(first_case_data + 7)],
                  m = 0L,
                  delay = diff(pexp(0L:(8L + 1L), 7/15))
                 )

S0_mac_case <- mac_results %>% filter(time == first_case_data - 1, matrix == "S") %>% pull(value)
I0_mac_case <- mac_results %>% filter(time == first_case_data - 1, matrix == "I") %>% pull(value)
R0_mac_case <- mac_results %>% filter(time == first_case_data - 1, matrix == "R") %>% pull(value)

betas_case <- fastbeta(series_case,
                       gamma = gamma,
                       init = c(round(S0_mac_case),
                                round(I0_mac_case),
                                round(R0_mac_case)),
                       m = 0L
)

betas <- rbind(betas_mort, betas_case) # should be able to just merge the two?
betas[first_case_data-1, "beta"] <- mean(beta[first_case_data-2, "beta"], beta[first_case_data, "beta"]) # interpolate this one missing beta

SI_loess <- stats::loess(
  formula = beta ~ c(1:pad_steps),
  data = data.frame(betas),
  span = 52/pad_steps,
  degree = 2,
  na.action = "na.exclude",
  control = loess.control(surface = "direct")
)

# Comparison of fastbeta and macpan R0's ----------------------------------

padded_term <- # For converting macpan beta to R0...
  c(
    rep(1/gamma, front_pad),
    (full_series$birth.trend/full_series$acm_trend/gamma),
    rep((full_series$birth.trend/full_series$acm_trend/gamma)[nrow(full_series)], end_pad - 1)
  )

padded_pop_term <- # For converting fastbeta beta to R0...
  padded_term * c(
    rep(full_series$pop[1], front_pad),
    full_series$pop,
    rep(full_series$pop[nrow(full_series)], end_pad - 1)
  )

df <- data.frame(x = c(full_series$numdate[1:front_pad], full_series$numdate, full_series$numdate[nrow(full_series) - end_pad + 1:(end_pad - 1)]),
                 week = 1:(pad_steps-1)
                 ) %>%
  mutate(loess_betas = SI_loess$fitted, # fastbeta beta
         loess_R0 = loess_betas * padded_pop_term, # fastbeta R0
         macpan_betas = betamac_df$beta[-nrow(betamac_df)],
         macpan_R0 = macpan_betas * padded_term
  )

## change these values to zoom into specific regions of the plot
lend <- 1
rend <- nrow(df)

pdf("betas.pdf", width = 7, height = 5)
df %>%
  left_join(full_series %>% mutate(x = numdate) %>% select(x, interpolated.deaths), by = "x") %>%
  pivot_longer(
    cols = c(loess_R0, macpan_R0, interpolated.deaths),
    names_to = "Type",
    values_to = "Value") %>%
  mutate(death = Type == "interpolated.deaths",
         Type = factor(Type, levels = c("loess_R0", "macpan_R0", "interpolated.deaths"), labels=c("Fastbeta", "Macpan", "Deaths"))) %>%
  mutate(death = factor(ifelse(death, "RGWR deaths", "Comparison"))) %>%
  mutate(pad = c(rep(TRUE, front_pad*3), rep(FALSE, steps*3), rep(TRUE, (end_pad-1)*3))) %>%
  filter(week <= rend, week >= lend) %>%
  ggplot() +
  geom_line(aes(x=week, y = Value, col = Type)) +
  facet_grid(rows = vars(death), scales = "free_y") +
  scale_colour_manual(values=c("blue", "red", "black")) +
  guides(color=guide_legend("Source")) +
  xlab("Week number") +
  ylab("$\\mathcal{R}_0$") +
  theme_bw() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=15),
        legend.key = element_blank(),
        legend.background = element_blank(),
        legend.position = c(.9,.85)) +
  geom_vline(xintercept = first_case_data)
dev.off()
