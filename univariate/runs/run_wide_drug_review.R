# Wide/path G-theory: drug_review
source("univariate/src/wide_models.R")
OUTPUT_DIR <- "univariate/outputs/wide"
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

dat <- ow_prepare_data(task = "drug_review", data_path = "data/drug_labeling_final.csv")
cat(sprintf("Dataset: drug_review | Items: %d | Rows: %d\n", nlevels(dat$item_id), nrow(dat)))

results <- ow_run_all_models(dat, extra_tries = 0L)

save(results, file = file.path(OUTPUT_DIR, "drug_review_results.RData"))
cat("\nSaved to", file.path(OUTPUT_DIR, "drug_review_results.RData"), "\n")
