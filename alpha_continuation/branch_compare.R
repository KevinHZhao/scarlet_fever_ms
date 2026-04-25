library(tidyverse)
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

period.set <- 1:7

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
                                     tryCatch(
                                     read.table(
                                       paste0("first_bif_branches/branch_p",iper,"_allinfo_macpan.dat"),
                                       col.names=auto.colnames),
                                   error = function(e) {return(NULL)}))

sin.allinfo.auto.data <- lapply(period.set,
                                function(iper)
                                  tryCatch(
                                  read.table(
                                    paste0("first_bif_branches/branch_p",iper,"_allinfo_sin.dat"),
                                    col.names=auto.colnames),
                                error = function(e) {return(NULL)}))

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
    if (!is.null(auto.data[[iper]])){
    with(auto.data[[iper]],{
      stable.pts <- which(ptype==1)
      unstable.pts <- which(ptype==2)
      points(R0[unstable.pts],ihi[unstable.pts],
             col="grey",pch=".")
      points(R0[stable.pts],ihi[stable.pts],
             col=colrs[iper],pch=19, cex=0.4)
    })
    }
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
    if (!is.null(auto.data.1[[iper]])){
    with(auto.data.1[[iper]],{
      transients <- ifelse(ev1.im == 0 | ptype == 2, NA, 2*pi*iper/Arg(ev1.re + ev1.im*1i))
      transientpts <- which(!is.na(transients))
      points(R0[transientpts],transients[transientpts],
             col=col2rgb(colrs) %>% t %>% rgb(alpha = 10, maxColorValue = 255) %>% .[iper],pch=".")
    })
    }
    if (!is.null(auto.data.2[[iper]])){
    with(auto.data.2[[iper]],{
      transients <- ifelse(ev1.im == 0 | ptype == 2, NA, 2*pi*iper/Arg(ev1.re + ev1.im*1i))
      transientpts <- which(!is.na(transients))
      points(R0[transientpts],transients[transientpts],
             col=col2rgb(colrs) %>% t %>% rgb(alpha = 10, maxColorValue = 255) %>% .[iper],pch=".")
    })
    }
  }
}

pdf("branch_compare.pdf", width = 8, height = 20)
par(mfrow = c(3,1))
autoplot(macpan.allinfo.auto.data, period.set, "macpan")
autoplot(sin.allinfo.auto.data, period.set, "sinusoidal")
transplot(macpan.allinfo.auto.data, sin.allinfo.auto.data)
dev.off()

## Try looking at phase portrait for when it seems transient period goes to infinity
