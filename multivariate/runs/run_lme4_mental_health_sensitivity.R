# Multivariate lme4 G-theory: mental_health_sensitivity
source("multivariate/src/lme4_models.R")
OUTPUT_DIR <- "multivariate/outputs/lme4"
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

results <- mv_lme4_run_all_models(
  task = "mental_health",
  data_path = "data/mh_labeling_final.csv",
  trait_subset = MG_MH_SENSITIVITY_TRAITS
)

save(results, file = file.path(OUTPUT_DIR, "mental_health_sensitivity_results.RData"))
cat("\nSaved to", file.path(OUTPUT_DIR, "mental_health_sensitivity_results.RData"), "\n")
