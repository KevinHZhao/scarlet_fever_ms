library(tidyverse)
library(deSolve)
library(parallel)
library(doParallel)

times <- seq(0,1000, by = 0.001)

n.cores <- Sys.getenv("SLURM_CPUS_PER_TASK")
registerDoParallel(cores = n.cores)
print(n.cores)
getDoParWorkers()

nR0 <- 300
namp <- 300
nS0 <- 10
nI0 <- 10
maxper <- 10

dyn.load(paste("SIRfun", .Platform$dynlib.ext, sep = ""))
# REDO WITH I0 between 0.0002/10 to 0.0002*10

for(chckpt in 3:9){
  bfd <-
    foreach(i = seq(1.6 + chckpt*3, 4.5 + chckpt*3, length.out = nR0/10), .combine = "rbind") %:%
    foreach(j = seq(1/300, 1, length.out = namp), .combine = "rbind") %:%
    foreach(k = seq(0.8/i, 1.2/i, length.out = nS0), .combine = "rbind") %:%
    foreach(l = seq(0.0002/10,  0.0002*10, length.out = nI0), .combine = "rbind") %dopar% {
      state <- c(S = k,
                 I = l)
      parameters <- c(mu = 0.02,
                      gam = 24.33,
                      R0 = i,
                      a = j)

      rawbfddata <- as.data.frame(ode(y = state, times = times, func = "SIRmap", parms = parameters, dllname="SIRfun", initfunc = "initmod")) %>%
        filter(time %in% 950:1000) %>%
        mutate(R0 = i, a = j, S0 = k, I0 = l)
      bfddataper <- rawbfddata %>% tail(n=1)
      period <- length(unique(round(log(rawbfddata$I),4)))
      bfddataper$period <- ifelse(period > maxper, NA, period)
      return(bfddataper)
    }
  saveRDS(bfd, file = paste("bfd/bfd_", chckpt+1, ".rds", sep = ""))
  print(paste("bfd_", chckpt+1, ".rds created succesfully", sep = ""))
}
## Took 1m21s doing 30 * 30 * 10 sims and with 0.01 time steps with 10 cores,
## expected to take a little over 3 days with 30 cores and these parms.
