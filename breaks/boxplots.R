library(tidyverse); theme_set(theme_bw())
library(macpan2)
library(gridExtra)
library(cowplot)
library(latex2exp)

wkyear <- 365.25/7

load("../CC_code/SF.RData")

results <- read.csv("../CC_code/output/Results.csv")
params <- read.csv("../CC_code/output/Params.csv")
final <- read.csv("../CC_code/output/Final.csv")

births <- read.csv("../CC_code/birthrate_1750_1930.csv")
full_series <- normalized_scarlet_fever_data %>%
  mutate(birth.trend = approx(x = births$numdate, y = births$birth.trend, xout = numdate)$y,
         pop = approx(x = births$numdate, y = births$pop, xout = numdate)$y,
         inner.pop = approx(x = births$numdate, y = births$inner.pop, xout = numdate)$y) %>%
  filter(numdate > 1842.01, numdate < 1930) %>%
  select(numdate, interpolated.deaths, birth.trend, acm_trend, pop, inner.pop)

steps <- nrow(full_series)
front_pad <- 479
end_pad <- 5*52
pad_steps <- steps + front_pad + end_pad
numyears <- 64
terms <- 3

X = rbf(pad_steps, numyears)
c <- params %>% filter(mat == "c") %>% pull(current)
S <- results %>% filter(row == "S", time != 0, time <= steps + front_pad) %>% pull(value)
sin_coeffs_mat <- params %>% filter(mat == "sin_coeffs_mat") %>% pull(current) %>% matrix(nrow = numyears)
cos_coeffs_mat <- params %>% filter(mat == "cos_coeffs_mat") %>% pull(current) %>% matrix(nrow = numyears)
b0 <- params %>% filter(mat == "b0") %>% pull(current)

beta_trend <- exp(b0 + X[(front_pad+1):(steps + front_pad),] %*% c)

sin_weights_mat <- X[(front_pad+1):(steps + front_pad),] %*% sin_coeffs_mat
cos_weights_mat <- X[(front_pad+1):(steps + front_pad),] %*% cos_coeffs_mat

breaks <- c(1842, read.csv("breaks.csv", header = FALSE, comment.char = "#")$V1, 1930)

# emd_sin_weights <- apply(sin_weights_mat, MARGIN = 2, FUN = function(x) emd(x, boundary = "wave")$residue)
# emd_cos_weights <- apply(cos_weights_mat, MARGIN = 2, FUN = function(x) emd(x, boundary = "wave")$residue)
# amp_trend = sqrt(emd_sin_weights^2 + emd_cos_weights^2)
# emd_beta <- emd(beta_trend, boundary = "wave")$residue

SINE = sin(2 * pi / wkyear * (1:pad_steps) %o% (1:terms))
COSINE = cos(2 * pi / wkyear * (1:pad_steps) %o% (1:terms))

raw_beta <- beta_trend * (1 + rowSums(
  sin_weights_mat * SINE[(front_pad):(steps + front_pad - 1), ] +
    cos_weights_mat * COSINE[(front_pad):(steps + front_pad - 1), ])
)

betafun <- function(index){
  return(beta_trend[[index]] * (1 +
                                  rowSums(t(sin_weights_mat[index,1:3] * t(sin(2*pi * (1/3650 * (1:(3650)) + 1/wkyear*(front_pad)) %o% (1:3))) +
                                              cos_weights_mat[index,1:3] * t(cos(2*pi * (1/3650 * (1:(3650)) + 1/wkyear*(front_pad)) %o% (1:3))))))
  )
}

# index = which row of the weights to use, from 1 to steps
betaintfun <- function(x, index){
  return(1 + rowSums(matrix(sin_weights_mat[index,1:3], ncol = 3) * sin(2*pi * x * (1:3)) +
                       matrix(cos_weights_mat[index,1:3], ncol = 3) * cos(2*pi * x * (1:3))
  )
  )
}

df <- expand.grid(y = seq(0, 1, length = 3650), x = full_series$numdate[seq(1, steps, by = 52)])

## Create a matrix of betas, with 3650 rows (representing time steps), and
## 88 columns (representing the beta macpan fitted at the beginning of each
## year in our time series).
beta_mat <- sapply(seq(1, steps, by = 52), betafun)

df$z <- c(t(t(sweep(beta_mat, 2, apply(beta_mat, 2, min))) / (apply(beta_mat, 2, max) - apply(beta_mat, 2, min))))
df$rawbeta <- c(beta_mat)

file.create("avg_beta_parms.txt")

pdf("boxplots.pdf", width = 14)
sin_means <- matrix(ncol = (length(breaks)-1), nrow = 3)
cos_means <- matrix(ncol = (length(breaks)-1), nrow = 3)
phases <- c()

monthSeq <- seq.Date(from = as.Date("0000-01-01"),
                     to = as.Date("0001-01-01"),
                     by = "1 month")

nxticks <- 13

for (i in 1:(length(breaks) - 1)) {
  p1 <- ggplot(full_series %>% filter(between(numdate, breaks[[i]], breaks[[i + 1]])) %>% mutate(x = numdate%%1, yr = floor(numdate))) +
    geom_line(aes(x = x, y = interpolated.deaths, group = yr)) +
    geom_line(data = full_series %>% filter(floor(numdate) == breaks[[i]]),
              aes(x = numdate%%1, y = interpolated.deaths),
              colour = "red") +
    geom_line(data = full_series %>% filter(floor(numdate) == breaks[[i + 1]] - 1),
              aes(x = numdate%%1, y = interpolated.deaths),
              colour = "magenta") +
    ggtitle(paste(
      "Observed deaths and fitted transmission rate oscillations between Jan",
      breaks[[i]],
      "(red) and Dec",
      breaks[[i + 1]] - 1,
      "(magenta)."
    )) +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    ) +
    ylab("Observed deaths") +
    scale_x_continuous(breaks = seq(0, 1, length = nxticks), minor_breaks = NULL) +
    coord_cartesian(ylim = c(0,150), xlim = c(-0.1, 1.1), expand = FALSE, clip = 'off')

  current_beta <- df %>% filter(between(x, breaks[[i]], breaks[[i + 1]])) %>% pull(rawbeta)
  indices <- intersect(which(between(full_series$numdate, breaks[[i]], breaks[[i + 1]])), seq(1, steps, by = 52))

  beta_trend_mean <- mean(beta_trend[indices])
  sin_weights_mat_mean <- colMeans(sin_weights_mat[indices, , drop = FALSE])
  cos_weights_mat_mean <- colMeans(cos_weights_mat[indices, , drop = FALSE])

  sin_means[,i] <- sin_weights_mat_mean
  cos_means[,i] <- cos_weights_mat_mean

  abs_macpan_function <- function(t){
    abs(
      sin_weights_mat_mean %*% sapply(X = c(2*pi*t, 4*pi*t, 6*pi*t), FUN = sin) +
        cos_weights_mat_mean %*% sapply(X = c(2*pi*t, 4*pi*t, 6*pi*t), FUN = cos)
      )
  }
  phase = optimize(f = abs_macpan_function, interval = c(0, 1), maximum = TRUE)$maximum
  phases <- c(phases, phase)

  write(paste("Segment ", i, ":\n" ,
              "b0=", beta_trend_mean, ";\n",
              "s1=", sin_weights_mat_mean[[1]], ";\ns2=", sin_weights_mat_mean[[2]], ";\ns3=", sin_weights_mat_mean[[3]],
              ";\nc1=", cos_weights_mat_mean[[1]], ";\nc2=", cos_weights_mat_mean[[2]], ";\nc3=", cos_weights_mat_mean[[3]],
              ";\nphase=", phase, ";",
              sep = ""),
        file = "avg_beta_parms.txt",
        sep="\n",
        append=TRUE
  )

  avg_beta_df <- tibble(y = unique(df$y), average_beta =
                          beta_trend_mean +
                          colSums(sin_weights_mat_mean * t(sin(2*pi * (unique(df$y) + 1/wkyear*(front_pad)) %o% (1:3))) +
                                    cos_weights_mat_mean * t(cos(2*pi * (unique(df$y) + 1/wkyear*(front_pad)) %o% (1:3)))
                          )
  ) %>%
    mutate(z = (average_beta - min(average_beta)) / (max(average_beta - min(average_beta))))


  beta_dist <- function(x, index) {
    return(
      (sapply(x, betaintfun, index = index) -
         (1 + colSums(sin_weights_mat_mean * t(sin(2*pi * (x + 1/wkyear*(front_pad)) %o% (1:3))) +
                        cos_weights_mat_mean * t(cos(2*pi * (x + 1/wkyear*(front_pad)) %o% (1:3)))
         )
         )
      )^2
    )
  }

  distances <- sqrt(
    unlist(
      sapply(
        X = indices,
        function(y) integrate(function(x) beta_dist(x, y), lower = 0, upper = 1)
      )[1,]
    )
  )

  p2 <- ggplot(df %>% filter(between(x, breaks[[i]], breaks[[i + 1]]))) +
    geom_line(data = avg_beta_df,
              aes(x = y, y = z),
              colour = "grey",
              linewidth = 3) +
    geom_line(aes(x = y, y = z, group = x)) +
    geom_line(data = df %>% filter(x == min(x[x > breaks[[i]]])),
              aes(x = y, y = z),
              colour = "red") +
    geom_line(data = df %>% filter(x == max(x[x < breaks[[i + 1]]])),
              aes(x = y, y = z),
              colour = "magenta") +
    # scale_fill_manual(values = c("red", rep("black", times = sum(between(unique(df$x), breaks[[i]], breaks[[i + 1]])) - 2), "blue")) +
    ylab("Normalized transmission rate") +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    ) +
    coord_cartesian(ylim=c(0,1), clip="off") +
    annotate("text",
             x = df$y[seq(1, 3650, length = nxticks)],
             y = 1.1,
             size = 3,
             label = paste("SD =", df %>%
                             filter(between(x, breaks[[i]], breaks[[i + 1]]),
                                    y %in% y[seq(1, 3650, length = nxticks)]) %>%
                             select(-rawbeta) %>%
                             pivot_wider(names_from = y, values_from = z) %>%
                             select(-x) %>%
                             sapply(sd) %>%
                             as.vector() %>%
                             round(digits = 3)
             )
    ) +
    annotate("text",
             x = -0.05,
             y = 0.9,
             size = 3,
             label = TeX(paste("$\\beta_H/\\gamma = $", round(max(current_beta)/(7/15), 3)))
    ) +
    annotate("text",
             x = -0.05,
             y = 0.8,
             size = 3,
             label = TeX(paste("$\\beta_L/\\gamma = $", round(min(current_beta)/(7/15), 3)))
    ) +
    annotate("text",
             x = -0.05,
             y = 0.7,
             size = 3,
             label = TeX(paste("$\\langle\\beta\\rangle/\\gamma = $", round(mean(current_beta)/(7/15), 3)))
    ) +
    annotate("text",
             x = -0.05,
             y = 0.6,
             size = 3,
             label = TeX(paste("$\\alpha = $", round((max(current_beta) - min(current_beta))/(2*mean(current_beta)), 3)))
    ) +
    annotate("text",
             x = 1.05,
             y = 0.9,
             size = 3,
             label = paste("[", signif(min(distances), 3), ",", signif(max(distances), 3), "]", sep = "")
    ) +
    annotate("text",
             x = 1.05,
             y = 0.8,
             size = 3,
             label = paste("Average dist:", signif(mean(distances), 2))
    ) +
    annotate("text",
             x = 1.05,
             y = 0.7,
             size = 3,
             label = paste("SD:", signif(sd(distances), 2))
    ) +
    scale_x_continuous(breaks = seq(0, 1, length = nxticks),
                       limits = c(-0.1,1.1), expand = c(0, 0),
                       minor_breaks = NULL)

  p3 <-
    ggplot(df %>% filter(between(x, breaks[[i]], breaks[[i + 1]]), y %in% y[seq(1, 3650, length = nxticks)])) +
    geom_boxplot(outlier.colour = "red", aes(y = z, group = y, x = y)) +
    xlab("Time in year") +
    ylab("Normalized transmission rate") +
    scale_x_continuous(breaks = seq(0, 1, length = nxticks),
                       labels = format(monthSeq, format = "%b"),
                       limits = c(-0.1,1.1),
                       expand = c(0, 0),
                       minor_breaks = NULL)

  grid::grid.newpage()
  grid::grid.draw(rbind(ggplotGrob(p1), ggplotGrob(p2), ggplotGrob(p3)))
}

dev.off()

saveRDS(object = list(sin_means = sin_means, cos_means = cos_means, phases = phases), file = "avg_beta_parms.Rds")
