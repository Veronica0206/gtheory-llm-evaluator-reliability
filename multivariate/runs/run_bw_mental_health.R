# Multivariate BW G-theory (nested models): mental_health (6 traits)
MAX_ITEMS <- 5L  # full dataset is slow; set NULL for all items

source("multivariate/src/bw_models.R")
OUTPUT_DIR <- "multivariate/outputs/bw"
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

dat_bundle <- mvbw_prepare_data(task = "mental_health",
                                data_path = "data/mh_labeling_final.csv")
cat(sprintf("Dataset: mental_health | Items: %d | Rows: %d | Traits: %d (%s)\n",
            nlevels(dat_bundle$data$item_id), nrow(dat_bundle$data),
            dat_bundle$p, paste(dat_bundle$traits, collapse = ", ")))

results <- mvbw_run_all_models(dat_bundle, extra_tries = 5L,
                               max_items = MAX_ITEMS)

save(results, file = file.path(OUTPUT_DIR, "mental_health_results.RData"))
cat("\nSaved to", file.path(OUTPUT_DIR, "mental_health_results.RData"), "\n")
