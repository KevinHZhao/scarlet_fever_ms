library(tidyverse); theme_set(theme_bw())

bif_dat <- list()

for (i in 1:8){
  bif_dat <- append(bif_dat, list(read.table(paste0("alpha_p_",i,".dat"))[,1:2]))
}

alphas <- sapply(X = bif_dat, function(df) approx(x = df[,1], y = df[,2], xout = 1)$y)
write.csv(alphas, file = "alphas.csv")

ggplot(bif_dat[[1]]) +
  geom_line(aes(x = V1, y = V2)) +
  xlim(0,1) +
  ylim(0,0.5) +
  xlab("Shape parameter, where 0 = macpan forcing, 1 = sinusoidal forcing") +
  ylab("Amplitude")
