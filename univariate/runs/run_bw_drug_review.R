# BW G-theory (nested models): drug_review
MAX_ITEMS <- NULL  # full dataset; set to 5L for fast testing

source("univariate/src/bw_models.R")
OUTPUT_DIR <- "univariate/outputs/bw"
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

dat <- bw_prepare_data(task = "drug_review", data_path = "data/drug_labeling_final.csv")
cat(sprintf("Dataset: drug_review | Items: %d | Rows: %d\n", nlevels(dat$item_id), nrow(dat)))

results <- bw_run_all_models(dat, extra_tries = 5L, max_items = MAX_ITEMS)

save(results, file = file.path(OUTPUT_DIR, "drug_review_results.RData"))
cat("\nSaved to", file.path(OUTPUT_DIR, "drug_review_results.RData"), "\n")
