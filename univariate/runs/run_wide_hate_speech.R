# Wide/path G-theory: hate_speech
source("univariate/src/wide_models.R")
OUTPUT_DIR <- "univariate/outputs/wide"
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

dat <- ow_prepare_data(task = "hate_speech", data_path = "data/hate_labeling_final.csv")
cat(sprintf("Dataset: hate_speech | Items: %d | Rows: %d\n", nlevels(dat$item_id), nrow(dat)))

results <- ow_run_all_models(dat, extra_tries = 0L)

save(results, file = file.path(OUTPUT_DIR, "hate_speech_results.RData"))
cat("\nSaved to", file.path(OUTPUT_DIR, "hate_speech_results.RData"), "\n")
