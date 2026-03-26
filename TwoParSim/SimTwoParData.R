library(tidyverse)
library(deSolve)
library(parallel)
library(doParallel)

# amps <- c(0.025,0.05,0.095, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)
# ampdata_names <- c("SF_SIR_R0_a025.dat",
#                    "SF_SIR_R0_a05.dat",
#                    "SF_SIR_R0_a095.dat",
#                    "SF_SIR_R0_a2.dat",
#                    "SF_SIR_R0_a3.dat",
#                    "SF_SIR_R0_a4.dat",
#                    "SF_SIR_R0_a5.dat",
#                    "SF_SIR_R0_a6.dat",
#                    "SF_SIR_R0_a7.dat",
#                    "SF_SIR_R0_a8.dat",
#                    "SF_SIR_R0_a9.dat",
#                    "SF_SIR_R0_a10.dat")
# 
# bfdlist <- lapply(ampdata_names,read.table,
#                   col.names=c("time","S","I","R0","log10S","log10I"))
# names(bfdlist) <- amps

# SIR <- function(t, state, parameters) {
#   with(as.list(c(state,parameters)), {
#     beta <- R0 * (gam + mu) * (1+a*cos(2*pi*t))
#     dS <- mu - beta*S*I - mu*S
#     dI = beta*S*I - (gam + mu)*I
#     
#     return(list(c(dS,dI)))
#   })
# }
times <- seq(0,1000, by = 0.001)

n.cores <- detectCores() - 2
my.cluster <- makeCluster(
  n.cores, 
  type = "PSOCK"
)
registerDoParallel(cl = my.cluster)

nR0 <- 10
namp <- 10
ninit <- 10

# paramlist <- expand.grid(R0 = seq(1, 30, length.out = nR0), 
#                          a = seq (0, 1, length.out = namp)) %>%
#   mutate(mu = 0.0231, gam = 24.37) %>%
#   relocate(gam) %>%
#   relocate(mu) %>%
#   t() %>%
#   as.vector() %>%
#   split(sort(rep_len(1:(nR0*namp), nR0*namp*4)))
# 
# clusterEvalQ(my.cluster, {
#   library(tidyverse)
#   library(deSolve)
#   dyn.load("SIRfun.dll")
# })

startT <- Sys.time()
bfdtest <-
  # parLapply(my.cluster, 
  #           X = paramlist, 
  #           fun = ode, 
  #           y = state, 
  #           times = times, 
  #           func = "SIRmap", 
  #           dllname = "SIRfun", 
  #           initfunc = "initmod") %>%
  # filter(time %in% 1950:2000)
  foreach(i = seq(1, 30, length.out = nR0), .combine = "rbind") %:%
    foreach(j = seq(0, 1, length.out = namp), .combine = "rbind") %:%
      foreach(k = seq(0.08, 0.12, length.out = ninit), .packages=c("tidyverse", "deSolve"), .combine = "rbind") %dopar% {
        dyn.load(paste("SIRfun", .Platform$dynlib.ext, sep = ""))
        
        state <- c(S = k,
                   I = 0.001)
        
        parameters <- c(mu = 0.0231,
                        gam = 24.37,
                        R0 = i,
                        a = j)
        return(as.data.frame(ode(y = state, times = times, func = "SIRmap", parms = parameters, dllname="SIRfun", initfunc = "initmod")) %>%
                       filter(time %in% 950:1000) %>%
                       mutate(R0 = i, a = j)
             )
      }
stopCluster(cl = my.cluster)
endT <- Sys.time()
print(endT - startT)
## With times from 0 to 1000, step size 0.001, using foreach loop
## Took roughly 3 minutes for 9 different r0 & a values, vs 47 seconds for 10 diff a values, r0 constant using multiprocessing (10 cores)
## 18.75 minutes for 300 simulations in pure R w/ multiprocessing, 25 seconds for 300 simulations using C dll!?!?!
## 1.13 hours for 30,000 simulations with time from 0 to 2000, step size still 0.001
## 2.134 mins using parLapply for 0 to 2000 with step size 0.001 vs 44.67 seconds using foreach loops, same 0 to 2k steps.
## 1.1646 mins for 1000 sims from 0 to 1000 with step 0.001

head(bfd)

maxper <- 8

last.point.with.period <- function(df, dop=4,
                                   R0lim=5, max.period=maxper) {
  n.cores <- detectCores() - 2
  my.cluster <- makeCluster(
    n.cores, 
    type = "PSOCK"
  )
  registerDoParallel(cl = my.cluster)
  ## list of data frames grouped by amplitude and R0 value
  dfR0alist <- df %>% mutate(index=rep(1:3630000,each=51)) %>% group_split(index)
  dfper <- 
    foreach(df.R0a = dfR0alist, .combine = "rbind", .packages = c("tidyverse")) %dopar% {
      ## compute period of this solution:
      period <- length(unique(round(log(df.R0a$I),dop)))
      df.R0a$period <- ifelse(period > max.period, NA, period)
      return(tail(df.R0a, n = 1))
    }
  stopCluster(cl = my.cluster)
  return(dfper)
}

startT <- Sys.time()
bfd.last <- last.point.with.period(df = bfd)
endT <- Sys.time()
print(endT - startT)
## 17.5 mins for 300*100 simulation data frame by calling the old function on the vector,
## 14.9 seconds for same data frame using parallelization and a few optimizations with 10 cores.

# foreach(i = seq(1, 30, length.out = nR0), .combine = "rbind") %:%
#   foreach(j = seq(0, 1, length.out = namp), .combine = "rbind", .packages = c("tidyverse")) %dopar% {
#     ## data frame with all pts on soln with given R0 and a:
#     df.R0aij <- df %>% filter(R0 == i, a == j)
#     ## compute period of this solution:
#     period <- length(unique(round(log(df.R0aij$I),dop)))
#     if (period > max.period) period <- NA
#     df.last[which(df.last$R0 == i & df.last$a == j),]$period <- period
#     return(df.last[which(df.last$R0 == i & df.last$a == j),])
#   }

# old.last.point.with.period <- function(df, dop=4,
#                                    R0lim=5, max.period=maxper) {
#   ## data frame with only the last pt on each soln:
#   df.last <- df %>% filter(time == max(time))
#   nR0 <- nrow(df.last %>% filter(a == 0)) # number of R0 values
#   namp <- nrow(df.last %>% filter(R0 == 4)) # number of amp values
#   df.last$period <- 0 # add period column
#   for (i in unique(df$R0)) {
#     for(j in unique(df$a)){
#       ## data frame with all pts on soln with given R0 and a:
#       df.R0aij <- df %>% filter(R0 == i, a == j)
#       ## compute period of this solution:
#       period <- length(unique(round(log(df.R0aij$I),dop)))
#       if (period > max.period) period <- NA
#       # }
#       df.last[which(df.last$R0 == i & df.last$a == j),]$period <- period
#     }
#   }
#   return(df.last)
# }

save(bfd, bfd.last, file = "Sims.RData")

# Careful with wrapped ICs
# Ask stephen lee about sharcnet and how he got it sorted out
# ask mark han if there are any question
# Earn compute canada code:
# Compute Canada Identifier (CCI) might be digital alliance identifier
# ijp-293
# Active Roles (-01 since Earn is PI)
# ijp-293-01
# Try using different number of years for converge
# email earn on tuesday 14th asking about meeting
# Try running one of the sims on my own end to see if it matches