library(tidyverse)
library(deSolve)
library(parallel)

set.seed(18301930) # Setting a random seed

bfd <- tibble()
for (i in 1:10){
  bfd <- bfd %>% bind_rows(readRDS(paste("bfd/bfd_", i, ".rds", sep = "")))
}

SIRCsim <- function(run.time, iState, delt = 1e-3, prms){
  return(as.data.frame(ode(y = iState,
                           times = seq(0, run.time, delt),
                           func = "SIRmap",
                           parms = prms,
                           dllname="SIRfun",
                           initfunc = "initmod")
                       ) %>%
           tail(n = 1) %>%
           select(-time) %>%
           unlist(., use.names=FALSE)
         )
}

transP <- function(eq,per,prm,hach=1e-6,rad=1e-5,nsmp=10,trn.tol=1e-4){
  angseq <- seq(0,by=2*pi/nsmp,length=nsmp)

  fun <- function(vec,par=prm){
    return(SIRCsim(run.time=1,iState=vec,delt=1e-3,prms=par))
  }

  trnp <- numeric(0)
  for(i in angseq){
    ini <- eq+c(rad*cos(i),rad*sin(i))
    dFdS<-fun(ini+c(hach,0))-fun(ini)
    dFdI<-fun(ini+c(0,hach))-fun(ini)
    A <- matrix(c(dFdS,dFdI),nrow=2,ncol=2,byrow=FALSE)
    eigs <- eigen(A)$values
    if(Arg(eigs[1]) != 0){
      trnp <- c(trnp,2*pi*per/Arg(eigs[1]))
    }
    if(Re(eigs[1])>0 | Re(eigs[2])>0){
    }
    if(Arg(eigs[1]) == 0){
      trnp <- c(trnp,0)
    }
  }
  # if(mean(trnp)!=0){
  #   if(var(trnp)/mean(trnp)>=trn.tol){
  #     return(NA)
  #   }
  return(c(avgtrans = mean(trnp), vartrans = var(trnp)))
  # }
  # if(mean(trnp)==0){
  #   return(NA)
  # }
}

bfdReader <- function(rowvec){
  eq <- rowvec[2:3]
  per <- rowvec[8]
  prm <- c(mu = 0.02,
           gam = 24.33,
           R0 = rowvec[4],
           a = rowvec[5])
  return(transP(eq, per, prm))
}

uniqueBFD <- bfd %>% distinct(S, I, R0, a, period, .keep_all = TRUE)

numCores <- 12
cl <- makeCluster(numCores)

clusterEvalQ(cl, {
  library(tidyverse)
  library(deSolve)

  dyn.load(paste("SIRfun", .Platform$dynlib.ext, sep = ""))
})

clusterExport(cl, c("uniqueBFD", "SIRCsim", "transP"))

start_time <- Sys.time()
per1 <- uniqueBFD %>% filter(period == 1)
per1 <- per1 %>% mutate(data.frame(t(rbind(parApply(cl, per1, 1, bfdReader)))))
end_time <- Sys.time()

end_time - start_time

start_time <- Sys.time()
per2 <- uniqueBFD %>% filter(period == 2)
per2 <- per2 %>% mutate(data.frame(t(rbind(parApply(cl, per2, 1, bfdReader)))))
end_time <- Sys.time()

end_time - start_time

stopCluster(cl)

full_per1 <- bfd %>%
  filter(period == 1) %>%
  left_join(per1 %>% select(S, I, R0, a, avgtrans, vartrans), by = c("S", "I", "R0", "a"))

full_per2 <- bfd %>%
  filter(period == 2) %>%
  left_join(per2 %>% select(S, I, R0, a, avgtrans, vartrans), by = c("S", "I", "R0", "a"))

# transient_no_na1 <- full_per1 %>% filter(!is.na(transient))
# transient_no_na2 <- full_per2 %>% filter(!is.na(transient))

saveRDS(full_per1, "Transient1.rds")
saveRDS(full_per2, "Transient2.rds")
# write.csv(transient_no_na1, "Transient_no_na1.csv")
# write.csv(transient_no_na2, "Transient_no_na2.csv")
