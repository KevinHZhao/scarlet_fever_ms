library(tidyverse)

breaks <- c(1842, read.csv("../breaks/breaks.csv", header = FALSE, comment.char = "#")$V1, 1930)

bfd <- list()
bfd_sin <- list()

for (i in 1:3){
  bfd <- append(
    bfd, list(
      read.table(paste("bfd/bruteforce_", i, ".dat", sep = ""),
                 col.names=c("time","S","I","R0","log10S","log10I"))
    )
  )
  bfd_sin <- append(
    bfd_sin, list(
      read.table(paste("bfd/bruteforce_", i, "_sin.dat", sep = ""),
                 col.names=c("time","S","I","R0","log10S","log10I"))
    )
  )
}

maxper <- 8

last.point.with.period <- function(df, dop=10,
                                   R0lim=4, max.period=maxper) {
  ## data frame with only the last pt on each soln:
  df.last <- subset(df, time==max(time))
  nR0 <- nrow(df.last) # number of R0 values
  df.last$period <- rep(0,nR0) # add period column
  for (i in 1:nR0) {
    ## data frame with all pts on soln with given R0:
    R0i <- df.last[i,"R0"]
    df.R0i <- subset(df,R0==R0i)
    ## compute period of this solution:
    if (R0i < R0lim) {
      period <- 1
    } else {
      period <- length(unique(round(df.R0i[,"I"],dop)))
      if (period > max.period) period <- NA
    }
    df.last$period[i] <- period
  }
  return(df.last)
}

bfd.last <- last.point.with.period(bfd[[1]])
bfd.sin.last <- last.point.with.period(bfd_sin[[1]])

bruteforceplot <- function(df, tsave=maxper, topmar=0, R0lim=1,
                           xlim=c(0,33), ...) {
  df <- subset(df, time>max(time)-tsave)
  dflp <- last.point.with.period(df)
  df$period <- rep(0,nrow(df)) # add period col to df
  for (i in 1:nrow(df)) {
    R0i <- df[i,"R0"]
    df$period[i] <- if (R0i > R0lim)
      dflp[dflp$R0==R0i,"period"] else 1
  }
  with(df,{
    par(mar=c(5,5,topmar,2)) # alter margins
    xlab <- expression(paste(
      "Basic Reproduction Number, ",R[0]))
    ylab <- "Prevalence, I/N"
    plot(R0, log10I, pch=".", yaxt="n", xaxs="i",
         xlim=xlim, las=1, xlab=xlab, ylab=ylab,
         cex.axis=1.3, cex.lab=1.5,
         col=period, cex=period, ...)
    ## add y-axis annotation
    y.ticks <- 0:-8
    y.label <- paste("10", y.ticks, sep="^")
    axis(2, at=y.ticks, las=2, label=parse(text=y.label))
  })
}

doublebfp <- function(df1, df2, i){
  par(mfrow = c(2,1))
  bruteforceplot(df1, ylim=c(-8,-2), topmar=3)
  title(paste("Bifurcation diagram for macpan (top) and sinusoidal (bottom) forcing between", breaks[[i]], "and", breaks[[i+1]]))
  bruteforceplot(df2, ylim=c(-8,-2), topmar=0)
}

pdf("bfd_plots.pdf", width = 14)
for (i in 1:length(bfd)){
  doublebfp(bfd[[i]], bfd_sin[[i]], i)
}
dev.off()

## extract all periods that are not NA:
all.periods <- with(bfd.last, period[!is.na(period)])
unique.periods <- unique(all.periods)
nper <- length(unique.periods)
## Create data frame for list of initial conditions:
ic.set <- bfd.last[1,]
## For each period that occurs, save a final condition:
for (iper in 1:nper) {
  df.iper <- subset(bfd.last,period==unique.periods[iper])
  ## choose the middle row to avoid bifurcation points:
  ic.set[iper,] <- df.iper[round(mean(1:nrow(df.iper))),]
}
## replace original row names with ic number:
row.names(ic.set) <- 1:nrow(ic.set)
## Save these points to be used as initial conditions:
# write.csv(ic.set, "icset_1.csv",
#           row.names=FALSE, quote=FALSE)

## extract all periods that are not NA:
all.periods.sin <- with(bfd.sin.last, period[!is.na(period)])
unique.periods.sin <- unique(all.periods.sin)
nper.sin <- length(unique.periods.sin)
## Create data frame for list of initial conditions:
ic.set.sin <- bfd.sin.last[1,]
## For each period that occurs, save a final condition:
for (iper in 1:nper.sin) {
  df.iper <- subset(bfd.sin.last,period==unique.periods.sin[iper])
  ## choose the middle row to avoid bifurcation points:
  ic.set.sin[iper,] <- df.iper[round(mean(1:nrow(df.iper))),]
}
## replace original row names with ic number:
row.names(ic.set.sin) <- 1:nrow(ic.set.sin)
## Save these points to be used as initial conditions:
# write.csv(ic.set.sin, "icset_1_sin.csv",
#           row.names=FALSE, quote=FALSE)

write(x=NULL,file = "icset.txt")
write_pts <- function(icset, p){
  for (i in 1:nrow(icset)){
    icset <- arrange(icset, period)
    write(
      paste0("set p",i," {init s=",icset[i,5],", init i=",icset[i,6],", R0=", icset[i,4], ", p=", p, ", alpha=0.1, nout=", icset[i,7], "}"),
      file = "icset.txt",
      append = TRUE
    )
  }
}

write_pts(ic.set, 0)
write_pts(ic.set.sin, 1)

## This just considers attractors, try to get attractors and transients like figure 2
## Use papst and earn to estimate alpha

# bfd <- bfd %>% filter(time > max(time) - )
# bfdlp <- last.point.with.period(bfd)
#
# ggplot(bfd)

auto.colnames <- c(
  "ptype", # point type (1=stable, 2=unstable)
  "branch", # branch number
  "p",
  "R0", # first active parameter
  "alpha", # second active parameter
  "period", # (not meaningful for discrete maps)
  "shi", # max S coordinate value
  "ihi", # max I coordinate value
  "slo", # min S coordinate value
  "ilo", # min I coordinate value
  "ev1.re", # real part of first eigenvalue
  "ev1.im", # imaginary part of first eigenvalue
  "ev2.re", # real part of second eigenvalue
  "ev2.im" # imaginary part of second eigenvalue
)

period.set <-
  sort(unique.periods[unique.periods <= maxper])

period.set <- 1:8

colrs = c("black",
           "red",
           "green",
           "blue",
           "cyan",
           "purple",
           "yellow",
           "orange")

macpan.allinfo.auto.data <- lapply(period.set,
                    function(iper)
                      read.table(
                        paste0("first_bif_branches/branch_p",iper,"_allinfo_macpan.dat"),
                        col.names=auto.colnames))

sin.allinfo.auto.data <- lapply(period.set,
                                 function(iper)
                                   read.table(
                                     paste0("first_bif_branches/branch_p",iper,"_allinfo_sin.dat"),
                                     col.names=auto.colnames))

# auto.data.macpan <- lapply(period.set,
#                     function(iper)
#                       read.table(
#                         paste0("first_bif_branches/branch_p",iper,"_macpan.dat"),
#                         col.names=c("R0", "ihi", "ilo", "ptype", "branch", "p")))
#
# auto.data.sin <- lapply(period.set,
#                            function(iper)
#                              read.table(
#                                paste0("first_bif_branches/branch_p",iper,"_sin.dat"),
#                                col.names=c("R0", "ihi", "ilo", "ptype", "branch", "p")))

autoplot <- function(auto.data, period.set, set,
                     xlim=c(0,33), ylim=c(-8,-2), ... ) {
  ## create empty plot with no y-axis annotation
  plot(x=0, type="n", xaxs="i", yaxt="n",
       xlim=xlim, ylim=ylim, cex.axis=1.3, cex.lab=1.5,
       xlab="",#expression(paste("Basic Reproduction Number, ", R[0])),
       ylab=paste0("Prevalence I/N (", set, ")"),
       mar = c(0, 4.1, 4.1, 2.1),
       ...)
  ## add y-axis annotation
  y.ticks <- 0:-8
  y.label <- paste("10", y.ticks, sep="^")
  axis(2, at=y.ticks, las=2, label=parse(text=y.label))
  ## plot the branch associated with each period
  for (iper in period.set) {
    with(auto.data[[iper]],{
      stable.pts <- which(ptype==1)
      unstable.pts <- which(ptype==2)
      points(R0[unstable.pts],ihi[unstable.pts],
             col="grey",pch=".")
      points(R0[stable.pts],ihi[stable.pts],
             col=colrs[iper],pch=19, cex=0.4)
    })
  }
  ## abline(v=17, col="red", lty="dashed")
}

transplot <- function(auto.data.1, auto.data.2, xlim=c(0,33), ylim=c(-8,-2), ... ){
  ## create empty plot with no y-axis annotation
  plot(x=0, type="n", xaxs="i", yaxt="n",
       xlim=xlim, ylim=c(0,30), cex.axis=1.3, cex.lab=1.5,
       xlab=expression(paste("Basic Reproduction Number, ", R[0])),
       ylab="Transient period", ...)
  ## add y-axis annotation
  y.ticks <- 0:30
  axis(2, at=c(5,10,15,20,25,30), las=2, label=c(5,10,15,20,25,30))
  ## plot the branch associated with each period
  for (iper in period.set) {
    with(auto.data.1[[iper]],{
      transients <- ifelse(ev1.im == 0 | ptype == 2, NA, 2*pi*iper/Arg(ev1.re + ev1.im*1i))
      transientpts <- which(!is.na(transients))
      points(R0[transientpts],transients[transientpts],
             col=col2rgb(colrs) %>% t %>% rgb(alpha = 10, maxColorValue = 255) %>% .[iper],pch=".")
    })
    with(auto.data.2[[iper]],{
      transients <- ifelse(ev1.im == 0 | ptype == 2, NA, 2*pi*iper/Arg(ev1.re + ev1.im*1i))
      transientpts <- which(!is.na(transients))
      points(R0[transientpts],transients[transientpts],
             col=col2rgb(colrs) %>% t %>% rgb(alpha = 10, maxColorValue = 255) %>% .[iper],pch=".")
    })
  }
}

pdf("Branch_compare.pdf", width = 8, height = 20)
par(mfrow = c(3,1))
autoplot(macpan.allinfo.auto.data, period.set, "macpan")
autoplot(sin.allinfo.auto.data, period.set, "sinusoidal")
transplot(macpan.allinfo.auto.data, sin.allinfo.auto.data)
dev.off()

## Try looking at phase portrait for when it seems transient period goes to infinity
