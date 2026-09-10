library(tidyverse)

simTrajectory <- readRDS("simTrajectory.rds")
traj_trans1 <- readRDS("Transients/Transient1_trajectory.rds")

yrly_traj <- tibble()
for (i in 1:10){
  yrly_traj <- yrly_traj %>%
    bind_rows(readRDS(paste("bfd/sim_trajectory_", i, ".rds", sep = "")))
}
yrly_traj <- yrly_traj %>%
  left_join(simTrajectory, by = c(t = "t", R0 = "R0", a = "amp")) %>%
  relocate(t) %>%
  left_join(traj_trans1 %>%
              distinct(S, I, R0, a, avgtrans, vartrans) %>%
              select(R0, a, S, I, avgtrans, vartrans) %>%
              rename(avgtrans1 = avgtrans, vartrans1 = vartrans),
            by = c("R0", "a", "S", "I"))

saveRDS(yrly_traj, file = "yearly_trajectory.rds")
