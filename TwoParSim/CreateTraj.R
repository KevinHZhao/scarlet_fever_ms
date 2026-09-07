library(macpan2)
library(tidyverse)

params <- read.csv("../CC_code/output/Params.csv")
alphas_raw <- read.csv("../alpha_continuation/alpha_p_data/alphas.csv")$alpha
indices <- read.csv("../alpha_continuation/alpha_p_data/alphas.csv")$index
breaks <- c(1842, read.csv("../breaks/breaks.csv", header = FALSE, comment.char = "#")$V1, 1939)
load("../CC_code/SF.RData")

full_series <- normalized_scarlet_fever_data %>%
  mutate(birth.trend = approx(x = births$numdate, y = births$birth.trend, xout = numdate)$y,
         pop = approx(x = births$numdate, y = births$pop, xout = numdate)$y,
         inner.pop = approx(x = births$numdate, y = births$inner.pop, xout = numdate)$y) %>%
  filter(numdate > 1842.01) %>%
  select(numdate, interpolated.deaths, birth.trend, acm_trend, pop, inner.pop)

steps <- nrow(full_series)
front_pad <- 479
end_pad <- 508
pad_steps <- steps + front_pad + end_pad
numyears <- 64
X <- rbf(pad_steps, numyears)
c <- params %>% filter(mat == "c") %>% pull(current)
b0 <- params %>% filter(mat == "b0") %>% pull(current)

alphas <- rep(NA, length(breaks))
alphas[indices] <- alphas_raw

beta_trend <- exp(b0 + X[front_pad + 1:steps,] %*% c)

firstInds <- match(1:(length(breaks) - 1), findInterval(full_series$numdate, breaks))
lastInds <- c(firstInds - 1, steps)[-1]

saveRDS(firstInds, "firstInds.rds")
saveRDS(lastInds, "lastInds.rds")

ampAdder <- function(df){
  maxSub <- 0.002
  return(seq(from = 0, to = maxSub, length.out = df[2] - df[1] + 1))
}

trajectory <- data.frame(t = full_series$numdate,
                         R0 = full_series$birth.trend/full_series$acm_trend * beta_trend/(7/15 + full_series$acm_trend/full_series$pop)) %>%
  mutate(amp = alphas[findInterval(full_series$numdate, breaks)],
         addedamp = amp -
           apply(data.frame(firstInds, lastInds), MARGIN = 1, ampAdder) %>% unlist)

simTrajectory <- data.frame(t = 1842:1939 + 0.5,
                            R0 = approx(x = trajectory$t, y =  trajectory$R0, xout = 1842:1939 + 0.5, na.rm = FALSE)$y,
                            amp = approx(x = trajectory$t, y =  trajectory$amp, xout = 1842:1939 + 0.5, na.rm = FALSE)$y
                            )

saveRDS(simTrajectory, file = "simTrajectory.rds")

transient1 <- readRDS("Transients/Transient1.rds")
transient2 <- readRDS("Transients/Transient2.rds")

bfd <- tibble()
for (i in 1:10){
  bfd <- bfd %>% bind_rows(readRDS(paste("bfd/bfd_", i, ".rds", sep = "")))
}

usedR0 <- unique(transient1$R0)
usedamp <- unique(transient1$a)

closest_R0 <- function(R0){
  usedR0[which(abs(usedR0 - R0) == min(abs(usedR0 - R0)))]
}

closest_amp <- function(amp){
  if (is.na(amp)) return(NA)
  usedamp[which(abs(usedamp - amp) == min(abs(usedamp - amp)))]
}

trajectory <- trajectory %>%
  mutate(close_R0 = sapply(R0, closest_R0),
         close_amp = sapply(amp, closest_amp)) %>%
  left_join(bfd %>% select(R0, a, S0, I0, period),
            by = join_by(close_R0 == R0, close_amp == a),
            relationship = "many-to-many") %>%
  left_join(transient1 %>%
              select(R0, a, S0, I0, avgtrans, vartrans) %>%
              rename(close_R0 = R0, close_amp = a, avgtrans1 = avgtrans, vartrans1 = vartrans),
            by = c("close_R0", "close_amp", "S0", "I0")) %>%
  left_join(transient2 %>%
              select(R0, a, S0, I0, avgtrans, vartrans) %>%
              rename(close_R0 = R0, close_amp = a, avgtrans2 = avgtrans, vartrans2 = vartrans),
            by = c("close_R0", "close_amp", "S0", "I0"))

saveRDS(trajectory, file = "trajectory.rds")
