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
ninits <- 100
maxper <- 10

dyn.load(paste("SIRfun", .Platform$dynlib.ext, sep = ""))

for(chckpt in 0:9){
  bfd <-
    foreach(i = seq(11/10 + chckpt*3, 4 + chckpt*3, length.out = nR0/10), .combine = "rbind") %:%
    foreach(j = seq(0, 1, length.out = namp), .combine = "rbind") %:% 
    foreach(k = seq(0.08, 0.12, length.out = ninits), .combine = "rbind") %dopar% {
      state <- c(S = k,
                 I = 0.001)
      parameters <- c(mu = 0.0231,
                      gam = 24.37,
                      R0 = i,
                      a = j)
      
      rawbfddata <- as.data.frame(ode(y = state, times = times, func = "SIRmap", parms = parameters, dllname="SIRfun", initfunc = "initmod")) %>%
        filter(time %in% 950:1000) %>%
        mutate(R0 = i, a = j, S0 = k)
      bfddataper <- rawbfddata %>% tail(n=1)
      period <- length(unique(round(log(rawbfddata$I),4)))
      bfddataper$period <- ifelse(period > maxper, NA, period)
      return(bfddataper)
    }
  saveRDS(bfd, file = paste("bfd_", chckpt+1, ".rds", sep = ""))
  print(paste("bfd_", chckpt+1, ".rds created succesfully", sep = ""))
}
## Took 1m21s doing 30 * 30 * 10 sims and with 0.01 time steps with 10 cores,
## expected to take a little over 3 days with 30 cores and these parms.