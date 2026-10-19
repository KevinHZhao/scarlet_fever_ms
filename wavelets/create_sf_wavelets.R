library(tidyverse)
library(WaveletComp)

load("../ms_data/SF.RData")

## Data from macpan2 fit
final <- read.csv("../ms_data/Final.csv")

full_series <- normalized_scarlet_fever_data %>%
  filter(numdate > year.start, numdate < year.end) %>%
  select(
    numdate,
    deaths,
    cases,
    sqrt.transformed.normalized.deaths,
    sqrt.transformed.normalized.cases,
    pop
  )

steps <- nrow(full_series)
front_pad <- 479
end_pad <- 508

## weekly = fitted weekly deaths by macpan 2
full_series <- full_series %>%
  mutate(
    macpan.deaths = final %>%
      filter(
        matrix %in% c("weekly"),
        between(row, front_pad + 1, front_pad + steps)
      ) %>%
      pull(value),
    sqrt.transformed.macpan.deaths = sqrt(macpan.deaths / pop)
  )

model_case_series <- c(
  head(full_series$cases, n = front_pad),
  full_series$cases,
  tail(full_series$cases, n = end_pad)
) # cases series with padding at the end
first_case_data <- which(!is.na(model_case_series))[[1]] # which time step we will have first case data
no_missing_case_data <- (1:(steps + front_pad))[
  -which(is.na(
    model_case_series
  ))
] -
  first_case_data +
  1

cases_sim <- final %>%
  filter(matrix == "infects") %>%
  pull(value) %>%
  .[no_missing_case_data]
case_series <- full_series %>%
  filter(!is.na(cases)) %>%
  mutate(
    macpan.cases = cases_sim,
    sqrt.transformed.macpan.cases = sqrt(macpan.cases / pop)
  )

SF.sqrt.norm.wavelet <- analyze.wavelet(
  my.data = full_series,
  my.series = "sqrt.transformed.normalized.deaths",
  dt = 1 / 52,
  dj = 1 / 50,
  lowerPeriod = 1 / 2,
  upperPeriod = 10,
  make.pval = TRUE,
  n.sim = 1000,
  date.format = "%Y-%m-%d"
)

saveRDS(SF.sqrt.norm.wavelet, file = "../ms_data/Wavelet-post-1842.rds")

SF.sqrt.norm.macpan.wavelet <- analyze.wavelet(
  my.data = full_series,
  my.series = "sqrt.transformed.macpan.deaths",
  dt = 1 / 52,
  dj = 1 / 50,
  lowerPeriod = 1 / 2,
  upperPeriod = 10,
  make.pval = TRUE,
  n.sim = 1000,
  date.format = "%Y-%m-%d"
)

saveRDS(
  SF.sqrt.norm.macpan.wavelet,
  file = "../ms_data/Wavelet-post-1842-macpan.rds"
)

SF.sqrt.norm.cases.wavelet <- analyze.wavelet(
  my.data = case_series,
  my.series = "sqrt.transformed.normalized.cases",
  dt = 1 / 52,
  dj = 1 / 50,
  lowerPeriod = 1 / 2,
  upperPeriod = 10,
  make.pval = TRUE,
  n.sim = 1000,
  date.format = "%Y-%m-%d"
)

saveRDS(
  SF.sqrt.norm.cases.wavelet,
  file = "../ms_data/Wavelet-post-1901-cases.rds"
)

SF.sqrt.norm.cases.macpan.wavelet <- analyze.wavelet(
  my.data = case_series,
  my.series = "sqrt.transformed.macpan.cases",
  dt = 1 / 52,
  dj = 1 / 50,
  lowerPeriod = 1 / 2,
  upperPeriod = 10,
  make.pval = TRUE,
  n.sim = 1000,
  date.format = "%Y-%m-%d"
)

saveRDS(
  SF.sqrt.norm.cases.macpan.wavelet,
  file = "../ms_data/Wavelet-post-1901-cases-macpan.rds"
)
