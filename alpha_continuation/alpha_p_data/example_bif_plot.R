library(tidyverse); theme_set(theme_bw())

bif_dat <- list()

indices <- c(1:9, 12, 13)
for (i in indices){
  bif_dat[[i]] <- read.table(paste0("alpha_p_",i,".dat"))[,1:2]
}

alphas <-
  sapply(X = indices, function(i) c(index = i, alpha = approx(x = bif_dat[[i]][,1], y = bif_dat[[i]][,2], xout = 1)$y)) %>%
  t()
write.csv(alphas, file = "alphas.csv")

ggplot(bif_dat[[1]]) +
  geom_line(aes(x = V1, y = V2)) +
  xlim(0,1) +
  ylim(0,0.5) +
  xlab("Shape parameter, where 0 = macpan forcing, 1 = sinusoidal forcing") +
  ylab("Amplitude")
