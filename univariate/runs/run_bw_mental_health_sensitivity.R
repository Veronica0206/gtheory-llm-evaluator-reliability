# BW G-theory (nested models): mental_health_sensitivity
MAX_ITEMS <- NULL  # full dataset; set to 5L for fast testing

source("univariate/src/bw_models.R")
OUTPUT_DIR <- "univariate/outputs/bw"
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

dat <- bw_prepare_data(task = "mental_health_3level", data_path = "data/mh_labeling_final.csv")
cat(sprintf("Dataset: mental_health_sensitivity | Items: %d | Rows: %d\n", nlevels(dat$item_id), nrow(dat)))

results <- bw_run_all_models(dat, extra_tries = 5L, max_items = MAX_ITEMS)

save(results, file = file.path(OUTPUT_DIR, "mental_health_sensitivity_results.RData"))
cat("\nSaved to", file.path(OUTPUT_DIR, "mental_health_sensitivity_results.RData"), "\n")
