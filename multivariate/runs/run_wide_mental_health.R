# Multivariate wide/path G-theory: mental_health
source("multivariate/src/wide_models.R")
OUTPUT_DIR <- "multivariate/outputs/wide"
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

dat_bundle <- mvow_prepare_data(
  task = "mental_health",
  data_path = "data/mh_labeling_final.csv"
)

results <- mvow_run_all_models(dat_bundle, extra_tries = 5L)

save(results, file = file.path(OUTPUT_DIR, "mental_health_results.RData"))
cat("\nSaved to", file.path(OUTPUT_DIR, "mental_health_results.RData"), "\n")
