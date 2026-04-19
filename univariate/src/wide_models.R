suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tibble))

source("univariate/src/shared/univariate_workflow.R")
source("univariate/src/fit_gtheory_openmx_path.R")


ow_log <- function(...) {
  cat(...)
  flush.console()
}


ow_timed <- function(label, expr) {
  t0 <- proc.time()[["elapsed"]]
  value <- force(expr)
  elapsed <- proc.time()[["elapsed"]] - t0
  ow_log(sprintf("%s completed in %.2fs\n", label, elapsed))
  value
}


ow_prepare_data <- function(task, data_path) {
  mg_prepare_univariate_data(task = task, data_path = data_path)
}


ow_modeled_facets <- function(facet_spec) {
  names(facet_spec)[vapply(facet_spec, function(x) identical(x$role, "facet"), logical(1))]
}


ow_make_facet_spec <- function(modeled_facets, nested_within = NULL) {
  all_facets <- c("evaluator", "prompt_type", "temperature", "seed")
  modeled_facets <- intersect(modeled_facets, all_facets)
  nested_facets <- setdiff(all_facets, modeled_facets)

  spec <- list()
  for (f in modeled_facets) {
    spec[[f]] <- list(role = "facet")
  }

  for (f in nested_facets) {
    within_f <- nested_within[[f]] %||% modeled_facets
    spec[[f]] <- list(role = "nested", within = within_f)
  }

  spec[all_facets]
}


ow_active_component_names <- function(modeled_facets) {
  comp_maps <- build_component_maps(
    crossed_facets = modeled_facets,
    item_interaction_orders = 1L,
    facet_interaction_orders = 2L,
    n_rep_groups = 1L
  )
  c(names(comp_maps$item), names(comp_maps$facet), "epsilon")
}


ow_full_start_values <- function(dat) {
  full_spec <- ow_model_specs()[["full"]]$facet_spec
  full_prep <- prepare_openmx_data_design(
    dat = as.data.frame(dat),
    facet_spec = full_spec,
    item_col = "item_id",
    score_col = "severity",
    n_rep_groups = 1L
  )
  estimate_raw_start_values(full_prep)
}


ow_clean_start_value <- function(x, floor = 1e-6) {
  val <- as.numeric(x[[1]])
  if (!is.finite(val) || is.na(val)) {
    val <- 0
  }
  max(val, floor)
}


ow_start_values_for_model <- function(dat, model_type, full_starts = NULL) {
  specs <- ow_model_specs()
  spec <- specs[[model_type]]

  if (is.null(full_starts)) {
    full_starts <- ow_full_start_values(dat)
  }

  if (identical(model_type, "full")) {
    return(full_starts)
  }

  reduced_prep <- prepare_openmx_data_design(
    dat = as.data.frame(dat),
    facet_spec = spec$facet_spec,
    item_col = "item_id",
    score_col = "severity",
    n_rep_groups = 1L
  )
  reduced_starts <- estimate_raw_start_values(reduced_prep)
  starts <- lapply(reduced_starts, ow_clean_start_value)
  starts[["epsilon"]] <- as.numeric(starts[["epsilon"]]) * 8

  if (identical(model_type, "nested_seed")) {
    temp_floor <- max(as.numeric(starts[["epsilon"]]) * 5e-3, 5e-4)
    for (nm in c("tau_temperature", "temperature",
                 "evaluator_temperature", "prompt_type_temperature")) {
      if (nm %in% names(starts)) {
        starts[[nm]] <- max(as.numeric(starts[[nm]]), temp_floor)
      }
    }
  }

  starts
}


ow_model_specs <- function() {
  list(
    full = list(
      facet_spec = ow_make_facet_spec(c("evaluator", "prompt_type", "temperature", "seed"))
    ),
    nested_seed = list(
      facet_spec = ow_make_facet_spec(c("evaluator", "prompt_type", "temperature"))
    ),
    nested_seed_temperature = list(
      facet_spec = ow_make_facet_spec(c("evaluator", "prompt_type"))
    )
  )
}


ow_component_table <- function(path_result) {
  tibble(
    component = names(path_result$vc_list),
    variance = as.numeric(unlist(path_result$vc_list))
  ) %>%
    mutate(prop_total = variance / sum(variance)) %>%
    arrange(desc(variance))
}


ow_fit_model <- function(dat,
                         model_type = c("full", "nested_seed", "nested_seed_temperature"),
                         extra_tries = 5L,
                         full_starts = NULL) {
  model_type <- match.arg(model_type)
  specs <- ow_model_specs()
  spec <- specs[[model_type]]
  starts <- ow_start_values_for_model(dat, model_type, full_starts = full_starts)

  ow_log("\n", strrep("=", 50), "\n", sep = "")
  ow_log(sprintf("Fitting wide/path OpenMx model [%s]\n", model_type))
  if (identical(model_type, "full")) {
    ow_log("  Start values: full-model raw-data effect decomposition\n")
  } else {
    ow_log("  Start values: reduced-data effect decomposition\n")
  }

  fit <- ow_timed(
    sprintf("  Wide/path optimization [%s]", model_type),
    fit_gtheory_openmx_path_design(
      dat = as.data.frame(dat),
      facet_spec = spec$facet_spec,
      score_col = "severity",
      n_rep_groups = 1L,
      model_name = paste0("gtheory_", model_type, "_wide"),
      extra_tries = extra_tries,
      start_values = starts,
      use_lme4_starts = FALSE
    )
  )

  comp_tbl <- ow_component_table(fit)
  mf <- fit$model_fit

  ow_log(sprintf(
    "  Status: %d | -2LL: %.4f | AIC: %.4f | BIC: %.4f\n",
    mf$status, mf$minus2LL, mf$AIC, mf$BIC
  ))

  list(
    model_type = model_type,
    fit = fit$fit,
    status = mf$status,
    minus2LL = mf$minus2LL,
    AIC = mf$AIC,
    BIC = mf$BIC,
    components = comp_tbl,
    path_result = fit
  )
}


ow_run_all_models <- function(dat, extra_tries = 5L) {
  full_starts <- ow_timed(
    "  Wide start-value setup [full decomposition]",
    ow_full_start_values(dat)
  )

  results <- list()

  for (mt in c("full", "nested_seed", "nested_seed_temperature")) {
    res <- ow_fit_model(dat, model_type = mt, extra_tries = extra_tries,
                        full_starts = full_starts)
    results[[mt]] <- res

    ow_log("\nVariance components:\n")
    print(res$components, n = Inf, width = Inf)
  }

  ow_log("\n\n", strrep("=", 50), "\n", sep = "")
  ow_log("MODEL COMPARISON\n")
  ow_log(strrep("=", 50), "\n\n", sep = "")

  summary_tbl <- bind_rows(lapply(names(results), function(mt) {
    r <- results[[mt]]
    tibble(
      model = mt,
      status = r$status,
      minus2LL = r$minus2LL,
      AIC = r$AIC,
      BIC = r$BIC,
      n_components = nrow(r$components)
    )
  }))
  print(summary_tbl, width = Inf)

  results
}
