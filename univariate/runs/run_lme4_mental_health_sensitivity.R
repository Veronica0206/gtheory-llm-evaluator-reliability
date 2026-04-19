# lme4 G-theory: mental_health_sensitivity
source("univariate/src/lme4_models.R")
OUTPUT_DIR <- "univariate/outputs/lme4"
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

models <- c("full", "nested_seed", "nested_seed_temperature")
results <- list()
for (m in models) {
  cat(sprintf("\n=== Fitting lme4 [%s] on mental_health_sensitivity ===\n", m))
  results[[m]] <- umc_fit_lme4_model("mental_health_sensitivity", m)
  cat(sprintf("  Status: singular=%s | -2LL: %.4f | AIC: %.4f | BIC: %.4f\n",
              results[[m]]$summary$singular,
              results[[m]]$summary$minus2LL,
              results[[m]]$summary$AIC,
              results[[m]]$summary$BIC))
}

save(results, file = file.path(OUTPUT_DIR, "mental_health_sensitivity_results.RData"))
cat("\nSaved to", file.path(OUTPUT_DIR, "mental_health_sensitivity_results.RData"), "\n")
