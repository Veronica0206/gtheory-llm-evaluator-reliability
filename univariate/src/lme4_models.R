source("univariate/compare_common_methods.R")  # provides umc_prepare_dataset, umc_model_specs, etc.


umc_fit_lme4_model <- function(dataset_name, model_name, use_REML = FALSE) {
  ds <- umc_prepare_dataset(dataset_name)
  dat <- ds$data
  specs <- umc_model_specs()

  if (!model_name %in% names(specs)) {
    stop("Unknown model: ", model_name)
  }

  fit <- lmer(
    mg_univariate_model_formulas()[[model_name]],
    data = dat,
    REML = use_REML,
    control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000))
  )

  comp_tbl <- mg_univariate_component_table(fit)
  metrics_tbl <- umc_metrics_from_component_table(
    comp_tbl,
    dat,
    modeled_facets = specs[[model_name]]$modeled_facets
  )

  summary_tbl <- tibble(
    dataset = ds$label,
    task = ds$task,
    model = model_name,
    method = "lme4",
    estimation = if (use_REML) "REML" else "ML",
    status = 0L,
    singular = lme4::isSingular(fit, tol = 1e-4),
    minus2LL = -2 * as.numeric(logLik(fit)),
    AIC = AIC(fit),
    BIC = BIC(fit)
  ) %>%
    bind_cols(metrics_tbl)

  components_out <- umc_standardize_components(comp_tbl, "lme4", ds$label) %>%
    mutate(model = model_name, .before = method)

  list(
    dataset = ds$label,
    task = ds$task,
    model = model_name,
    method = "lme4",
    fit = fit,
    summary = summary_tbl,
    components = components_out
  )
}


args <- commandArgs(trailingOnly = TRUE)
if (sys.nframe() == 0L && length(args) >= 2) {
  dataset_name <- args[[1]]
  model_name <- args[[2]]
  use_REML <- if (length(args) >= 3) as.logical(args[[3]]) else FALSE
  res <- umc_fit_lme4_model(dataset_name = dataset_name, model_name = model_name, use_REML = use_REML)
  prefix <- paste0(res$dataset, "_", res$model, "_lme4")
  umc_write_outputs(res, prefix)
  print(res$summary, width = Inf)
  print(res$components, n = Inf, width = Inf)
}
