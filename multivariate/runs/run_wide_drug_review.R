# Multivariate wide/path G-theory: drug_review
source("multivariate/src/wide_models.R")
OUTPUT_DIR <- "multivariate/outputs/wide"
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

dat_bundle <- mvow_prepare_data(
  task = "drug_review",
  data_path = "data/drug_labeling_final.csv"
)

results <- mvow_run_all_models(dat_bundle, extra_tries = 5L)

save(results, file = file.path(OUTPUT_DIR, "drug_review_results.RData"))
cat("\nSaved to", file.path(OUTPUT_DIR, "drug_review_results.RData"), "\n")
