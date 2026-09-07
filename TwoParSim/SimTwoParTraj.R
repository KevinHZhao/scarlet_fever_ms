library(tidyverse)
library(deSolve)
library(parallel)
library(doParallel)

simTrajectory <- readRDS("simTrajectory.rds")

times <- seq(0,1000, by = 0.001)

n.cores <- Sys.getenv("SLURM_CPUS_PER_TASK")
registerDoParallel(cores = n.cores)
print(n.cores)
getDoParWorkers()

nS0 <- 100
nI0 <- 100
maxper <- 10

dyn.load(paste("SIRfun", .Platform$dynlib.ext, sep = ""))

for(chckpt in 0:9){
  bfd <-
    foreach(i = seq(1 + chckpt*10, min(10 + chckpt*10, 98), length.out = 10), .combine = "rbind") %:%
    foreach(j = seq(0.8/simTrajectory$R0[i], 1.2/simTrajectory$R0[i], length.out = nS0), .combine = "rbind") %:%
    foreach(k = seq(0.0002/10,  0.0002*10, length.out = nI0), .combine = "rbind") %dopar% {
      if (is.na(simTrajectory$amp[i])) return(data.frame(t = simTrajectory$t[i], S = NA, I = NA, R0 = simTrajectory$R0[i], a = NA, S0 = j, I0 = k, period = NA))
      state <- c(S = j,
                 I = k)
      parameters <- c(mu = 0.02,
                      gam = 24.33,
                      R0 = simTrajectory$R0[i],
                      a = simTrajectory$amp[i])

      rawbfddata <- as.data.frame(ode(y = state, times = times, func = "SIRmap", parms = parameters, dllname="SIRfun", initfunc = "initmod")) %>%
        filter(time %in% 950:1000) %>%
        mutate(R0 = simTrajectory$R0[i], a = simTrajectory$amp[i], S0 = j, I0 = k)
      bfddataper <- rawbfddata %>% tail(n=1)
      period <- length(unique(round(log(rawbfddata$I),4)))
      bfddataper$period <- ifelse(period > maxper, NA, period)
      bfddataper <- bfddataper %>% rename(t = time) %>% mutate(t = simTrajectory$t[i])
      return(bfddataper)
    }
  saveRDS(bfd, file = paste("bfd/sim_trajectory_", chckpt+1, ".rds", sep = ""))
  print(paste("sim_trajectory_", chckpt+1, ".rds created succesfully", sep = ""))
}
