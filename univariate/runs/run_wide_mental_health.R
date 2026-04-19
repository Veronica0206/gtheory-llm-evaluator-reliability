# Wide/path G-theory: mental_health
source("univariate/src/wide_models.R")
OUTPUT_DIR <- "univariate/outputs/wide"
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

dat <- ow_prepare_data(task = "mental_health", data_path = "data/mh_labeling_final.csv")
cat(sprintf("Dataset: mental_health | Items: %d | Rows: %d\n", nlevels(dat$item_id), nrow(dat)))

results <- ow_run_all_models(dat, extra_tries = 0L)

save(results, file = file.path(OUTPUT_DIR, "mental_health_results.RData"))
cat("\nSaved to", file.path(OUTPUT_DIR, "mental_health_results.RData"), "\n")
