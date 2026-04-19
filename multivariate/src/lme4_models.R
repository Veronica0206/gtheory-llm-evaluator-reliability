# ---------------------------------------------------------------------------
# lme4_models.R (MULTIVARIATE)
#
# Multivariate G-theory via lme4 with double-bar (||) syntax.
# Suppresses cross-trait correlations within each random-effect group.
# Fits all 3 model types: full, nested_seed, nested_seed_temperature.
#
# Expected working directory: 01. Research/
# ---------------------------------------------------------------------------

source("multivariate/src/shared/multivariate_workflow.R")

suppressPackageStartupMessages(library(lme4))


mv_lme4_fit_model <- function(task, data_path, model_name,
                               trait_subset = NULL, use_REML = FALSE) {
  dat <- mg_prepare_multivariate_data(
    task = task, data_path = data_path, trait_subset = trait_subset
  )
  dat <- mg_add_trait_indicators(dat)
  trait_cols <- mg_trait_columns(dat)

  model_groups <- mg_multivariate_model_groups(include_full = TRUE)
  if (!model_name %in% names(model_groups)) {
    stop("Unknown model: ", model_name,
         ". Available: ", paste(names(model_groups), collapse = ", "))
  }

  groups <- model_groups[[model_name]]
  formula_obj <- mg_make_multivariate_formula(dat, groups = groups)
  cat(sprintf("Fitting multivariate lme4 [%s] — %d traits, %d groups\n",
              model_name, length(trait_cols), length(groups)))

  control <- lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000))
  fit <- lmer(formula_obj, data = dat, REML = use_REML, control = control)

  singular <- lme4::isSingular(fit, tol = 1e-4)
  minus2LL <- -2 * as.numeric(logLik(fit))

  cat(sprintf("  Singular: %s | -2LL: %.4f | AIC: %.4f | BIC: %.4f\n",
              singular, minus2LL, AIC(fit), BIC(fit)))

  list(
    model_name = model_name,
    fit        = fit,
    singular   = singular,
    minus2LL   = minus2LL,
    AIC        = AIC(fit),
    BIC        = BIC(fit)
  )
}


mv_lme4_run_all_models <- function(task, data_path,
                                    trait_subset = NULL, use_REML = FALSE) {
  models <- c("full", "nested_seed", "nested_seed_temperature")
  results <- list()

  for (m in models) {
    cat(sprintf("\n=== Multivariate lme4 [%s] ===\n", m))
    results[[m]] <- mv_lme4_fit_model(
      task = task, data_path = data_path,
      model_name = m, trait_subset = trait_subset, use_REML = use_REML
    )
  }

  results
}
