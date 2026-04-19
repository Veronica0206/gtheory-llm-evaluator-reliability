# Multivariate lme4 G-theory: drug_review
source("multivariate/src/lme4_models.R")
OUTPUT_DIR <- "multivariate/outputs/lme4"
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

results <- mv_lme4_run_all_models(
  task = "drug_review",
  data_path = "data/drug_labeling_final.csv"
)

save(results, file = file.path(OUTPUT_DIR, "drug_review_results.RData"))
cat("\nSaved to", file.path(OUTPUT_DIR, "drug_review_results.RData"), "\n")
