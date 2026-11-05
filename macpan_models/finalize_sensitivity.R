library(tidyverse)

max_n <- 0

for (n in 1:10) {
  # Adjust range as needed
  dir_name <- file.path(paste0("output_sensitivity_", n))
  if (dir.exists(dir_name)) {
    max_n <- n # Update max_n if the directory exists
  }
}

all_res <- readRDS(paste0("output_sensitivity_", max_n, "/all_res.RDS"))
objs <- sapply(all_res, function(res) round(res$objective, digits = 6))
best_ind <- which(objs == sort(unique(objs))[4])[[1]] # some of the fits actually gives an unrealistic fit with R0 that grows very fast after 1920, so go for lowest min without these issues
CFP_parms <- unlist(readRDS(paste0(
  "output_sensitivity_",
  max_n,
  "/CFP_parms.RDS"
))[best_ind, ])

# it seems like the 32211.25 and 31974.42 minimum is singular convergence (fake min)
# objs[[49]]
# objs[[22]]

table(objs)

file.copy(
  paste0(
    "output_sensitivity_",
    max_n,
    "/Final_",
    CFP_parms[["CFP_min"]],
    "_",
    CFP_parms[["CFP_max"]],
    "_",
    CFP_parms[["CFP_mid"]],
    "_",
    CFP_parms[["CFP_rate"]],
    "_",
    CFP_parms[["case_reporting_rate"]],
    ".csv"
  ),
  "../ms_data/Final.csv",
  overwrite = TRUE
)
file.copy(
  paste0(
    "output_sensitivity_",
    max_n,
    "/Results_",
    CFP_parms[["CFP_min"]],
    "_",
    CFP_parms[["CFP_max"]],
    "_",
    CFP_parms[["CFP_mid"]],
    "_",
    CFP_parms[["CFP_rate"]],
    "_",
    CFP_parms[["case_reporting_rate"]],
    ".csv"
  ),
  "../ms_data/Results.csv",
  overwrite = TRUE
)
file.copy(
  paste0(
    "output_sensitivity_",
    max_n,
    "/Params_",
    CFP_parms[["CFP_min"]],
    "_",
    CFP_parms[["CFP_max"]],
    "_",
    CFP_parms[["CFP_mid"]],
    "_",
    CFP_parms[["CFP_rate"]],
    "_",
    CFP_parms[["case_reporting_rate"]],
    ".csv"
  ),
  "../ms_data/Params.csv",
  overwrite = TRUE
)
