library(tidyverse)
library(deSolve)
library(parallel)
library(shiny)

# ====Asymptotic sim====

times <- seq(0, 1000, by = 0.001) # Simulate 1000 years of ODE with 0.001 year step size

n.cores <- detectCores() - 2 # use two less than number of cores on this device

nR0 <- 10 # number of diff R_0 to simulate across
namp <- 10 # number of diff alphas to simulate across
possible_gammas <- 365.25/(10:25) # simulate across 10 to 25 day recovery periods

mu = 0.02

# make sure you have SIRfun.so compiled on mac (on windows I think should be SIRfun.dll?)
dyn.load(paste("SIRfun", .Platform$dynlib.ext, sep = ""))
grid <- expand.grid(i = seq(1.6, 31.6, length.out = nR0), j = seq(1/namp, 1, length.out = namp), k = possible_gammas)

startT <- Sys.time()
bfd <-
  mcmapply(
    FUN = function(i, j, k) { # i = R0, j = amp, k = gamma
      state <- c(S = 1/i,
                 I = 0.0002) # start close to the attractor for the SIR model without forcing

      parameters <- c(mu = mu,
                      gam = k,
                      R0 = i,
                      a = j)

      return(as.data.frame(ode(y = state, times = times, func = "SIRmap", parms = parameters, dllname="SIRfun", initfunc = "initmod")) %>%
                     filter(time %in% 950:1000) %>%
                     mutate(R0 = i, a = j, gam = k)
             )
    },
    i = grid$i,
    j = grid$j,
    k = grid$k,
    SIMPLIFY=FALSE,
    mc.cores = n.cores
  )
endT <- Sys.time()
print(paste0("Attractor sim took ", round(as.numeric(endT - startT, units = "secs"), digits = 3), " seconds"))

maxper <- 10

startT <- Sys.time()
bfd.last <- mclapply(
  X = bfd,
  FUN = function(x){
    dop <- 4 # decimals of precision
    period <- length(unique(round(log(x$I), dop)))
    x$period <- ifelse(period > maxper, NA, period)
    return(last(x))
  },
  mc.cores = n.cores
) %>%
  bind_rows()
endT <- Sys.time()
print(paste0("Attractor period calculation from sim results took ", round(as.numeric(endT - startT, units = "secs"), digits = 3), " seconds"))

# ====Transient sim====

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

transP <- function(S,I,per,gam,R0,a,hach=1e-6,rad=1e-5,nsmp=10,trn.tol=1e-4){
  angseq <- seq(0,by=2*pi/nsmp,length=nsmp)

  fun <- function(vec){
    return(SIRCsim(run.time=1,iState=vec,delt=1e-3,prms=c(mu = mu, gam = gam, R0 = R0, a = a)))
  }

  trnp <- numeric(0)
  for(i in angseq){
    ini <- c(S,I)+c(rad*cos(i),rad*sin(i))
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

start_time <- Sys.time()
per1 <- bfd.last %>% filter(period == 1)
per1 <-
  per1 %>%
  cbind(
    mcmapply(
      FUN = transP,
      S = .$S,
      I = .$I,
      per = .$period,
      gam = .$gam,
      R0 = .$R0,
      a = .$a
    ) %>%
    t()
  )
end_time <- Sys.time()

end_time - start_time

start_time <- Sys.time()
per2 <- bfd.last %>% filter(period == 2)
per2 <-
  per2 %>%
  cbind(
    mcmapply(
      FUN = transP,
      S = .$S,
      I = .$I,
      per = .$period,
      gam = .$gam,
      R0 = .$R0,
      a = .$a
    ) %>%
      t()
  )
end_time <- Sys.time()

print(paste0("Transient period calculation from sim results took ", round(as.numeric(end_time - start_time, units = "secs"), digits = 3), " seconds"))

# ====shiny app to display the different plots====

values <- sort(unique(possible_gammas))
ui <- fluidPage(
  radioButtons(
    "status",
    "Select:",
    choices = setNames(values, round(values, 2)),
    inline = TRUE
  ),

  radioButtons(
    inputId = "plot_type",
    label = NULL,
    choices = c("Plot A", "Plot B"),
    selected = "Plot A"
  ),

  plotOutput("plot")
)

server <- function(input, output) {
  output$plot <- renderPlot({

    if (input$plot_type == "Plot A") {

    bfd.last %>%
      mutate(period = factor(period)) %>%
      filter(!is.na(period), gam == input$status) %>%
      ggplot() +
      geom_point(aes(x = R0, y = a, col = period)) +
      ylim(0, 1)

    } else {
      per1 %>%
        filter(vartrans < 1e-2, gam == input$status) %>%
        ggplot() +
        geom_point(aes(x = R0, y = a, col = avgtrans)) +
        ylim(0, 1)
    }
  })
}

shinyApp(ui, server)
