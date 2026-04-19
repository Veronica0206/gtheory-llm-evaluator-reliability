# Multivariate BW G-theory (nested models): drug_review
MAX_ITEMS <- 5L  # full dataset is slow; set NULL for all items

source("multivariate/src/bw_models.R")
OUTPUT_DIR <- "multivariate/outputs/bw"
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

dat_bundle <- mvbw_prepare_data(task = "drug_review",
                                data_path = "data/drug_labeling_final.csv")
cat(sprintf("Dataset: drug_review | Items: %d | Rows: %d | Traits: %d (%s)\n",
            nlevels(dat_bundle$data$item_id), nrow(dat_bundle$data),
            dat_bundle$p, paste(dat_bundle$traits, collapse = ", ")))

results <- mvbw_run_all_models(dat_bundle, extra_tries = 5L,
                               max_items = MAX_ITEMS)

save(results, file = file.path(OUTPUT_DIR, "drug_review_results.RData"))
cat("\nSaved to", file.path(OUTPUT_DIR, "drug_review_results.RData"), "\n")
