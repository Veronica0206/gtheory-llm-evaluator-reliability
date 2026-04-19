# lme4 G-theory: hate_speech
source("univariate/src/lme4_models.R")
OUTPUT_DIR <- "univariate/outputs/lme4"
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

models <- c("full", "nested_seed", "nested_seed_temperature")
results <- list()
for (m in models) {
  cat(sprintf("\n=== Fitting lme4 [%s] on hate_speech ===\n", m))
  results[[m]] <- umc_fit_lme4_model("hate_speech", m)
  cat(sprintf("  Status: singular=%s | -2LL: %.4f | AIC: %.4f | BIC: %.4f\n",
              results[[m]]$summary$singular,
              results[[m]]$summary$minus2LL,
              results[[m]]$summary$AIC,
              results[[m]]$summary$BIC))
}

save(results, file = file.path(OUTPUT_DIR, "hate_speech_results.RData"))
cat("\nSaved to", file.path(OUTPUT_DIR, "hate_speech_results.RData"), "\n")
