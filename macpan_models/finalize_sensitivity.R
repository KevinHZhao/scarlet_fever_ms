library(tidyverse)

max_n <- 0

for (n in 1:10) {  # Adjust range as needed
    dir_name <- file.path(paste0("output_sensitivity_", n))
    if (dir.exists(dir_name)) {
        max_n <- n  # Update max_n if the directory exists
    }
}

all_res <- readRDS(paste0("output_sensitivity_", max_n, "/all_res.RDS"))
best_ind <- which.min(sapply(all_res, function(res) res$objective))
CFP_parms <- unlist(readRDS(paste0("output_sensitivity_", max_n, "/CFP_parms.RDS"))[best_ind, ])

objectives <- sapply(all_res, \(res) round(res$objective, 5)) # it seems like the 32211.25 and 31974.42 minimum is singular convergence (fake min)
# objectives[[49]]
# objectives[[22]]

sapply(all_res, \(res) round(res$objective, 5)) %>%
    table()

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
      ".csv"
   ), 
   "../ms_data/Params.csv", 
   overwrite = TRUE
)
