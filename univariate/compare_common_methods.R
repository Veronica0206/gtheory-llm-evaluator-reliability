suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(lme4))

source("univariate/src/shared/univariate_workflow.R")


UMC_OUTPUT_DIR <- "univariate/outputs"


umc_all_facets <- function() {
  c("evaluator", "prompt_type", "temperature", "seed")
}


umc_full_spec <- function() {
  list(
    evaluator = list(role = "facet"),
    prompt_type = list(role = "facet"),
    temperature = list(role = "facet"),
    seed = list(role = "facet")
  )
}


umc_make_facet_spec <- function(modeled_facets, nested_within = NULL) {
  all_facets <- umc_all_facets()
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


umc_model_specs <- function() {
  list(
    full = list(
      modeled_facets = c("evaluator", "prompt_type", "temperature", "seed"),
      facet_spec = umc_make_facet_spec(c("evaluator", "prompt_type", "temperature", "seed"))
    ),
    nested_seed = list(
      modeled_facets = c("evaluator", "prompt_type", "temperature"),
      facet_spec = umc_make_facet_spec(c("evaluator", "prompt_type", "temperature"))
    ),
    nested_seed_temperature = list(
      modeled_facets = c("evaluator", "prompt_type"),
      facet_spec = umc_make_facet_spec(c("evaluator", "prompt_type"))
    ),
    nested_prompt_temperature_seed = list(
      modeled_facets = c("evaluator"),
      facet_spec = umc_make_facet_spec(c("evaluator"))
    ),
    prompt_only = list(
      modeled_facets = c("prompt_type"),
      facet_spec = umc_make_facet_spec(c("prompt_type"))
    ),
    item_only = list(
      modeled_facets = character(0),
      facet_spec = NULL
    )
  )
}


umc_model_names <- function(include_item_only = TRUE) {
  nms <- names(umc_model_specs())
  if (!include_item_only) {
    nms <- setdiff(nms, "item_only")
  }
  nms
}


umc_dataset_specs <- function() {
  list(
    hate_speech = list(
      task = "hate_speech",
      data_path = "data/hate_labeling_final.csv",
      label = "hate_speech"
    ),
    mental_health = list(
      task = "mental_health",
      data_path = "data/mh_labeling_final.csv",
      label = "mental_health"
    ),
    mental_health_sensitivity = list(
      task = "mental_health_3level",
      data_path = "data/mh_labeling_final.csv",
      label = "mental_health_sensitivity"
    ),
    drug_review = list(
      task = "drug_review",
      data_path = "data/drug_labeling_final.csv",
      label = "drug_review"
    )
  )
}


umc_prepare_dataset <- function(dataset_name) {
  specs <- umc_dataset_specs()
  if (!dataset_name %in% names(specs)) {
    stop(
      "Unknown dataset '", dataset_name, "'. ",
      "Choose from: ", paste(names(specs), collapse = ", "), "."
    )
  }

  spec <- specs[[dataset_name]]
  dat <- mg_prepare_univariate_data(
    task = spec$task,
    data_path = spec$data_path
  )

  list(
    dataset = dataset_name,
    task = spec$task,
    label = spec$label,
    data_path = spec$data_path,
    data = dat
  )
}


umc_key_components <- function() {
  c(
    "item_id",
    "evaluator",
    "prompt_type",
    "temperature",
    "seed",
    "item_id:evaluator",
    "item_id:prompt_type",
    "item_id:temperature",
    "item_id:seed",
    "evaluator:prompt_type",
    "evaluator:temperature",
    "evaluator:seed",
    "prompt_type:temperature",
    "prompt_type:seed",
    "temperature:seed",
    "Residual"
  )
}


umc_item_component_name <- function(facet) {
  switch(
    facet,
    evaluator = "item_id:evaluator",
    prompt_type = "item_id:prompt_type",
    temperature = "item_id:temperature",
    seed = "item_id:seed",
    stop("Unknown facet: ", facet)
  )
}


umc_facet_component_name <- function(facets) {
  facets <- sort(facets)
  key <- paste(facets, collapse = "|")
  switch(
    key,
    "evaluator|prompt_type" = "evaluator:prompt_type",
    "evaluator|temperature" = "evaluator:temperature",
    "evaluator|seed" = "evaluator:seed",
    "prompt_type|temperature" = "prompt_type:temperature",
    "prompt_type|seed" = "prompt_type:seed",
    "seed|temperature" = "temperature:seed",
    stop("Unknown facet pair: ", key)
  )
}


umc_design_counts <- function(dat) {
  c(
    evaluator = nlevels(dat$evaluator),
    prompt_type = nlevels(dat$prompt_type),
    temperature = nlevels(dat$temperature),
    seed = nlevels(dat$seed)
  )
}


umc_component_divisor <- function(component, dat) {
  counts <- umc_design_counts(dat)
  parts <- setdiff(strsplit(component, ":", fixed = TRUE)[[1]], "item_id")
  parts <- intersect(parts, names(counts))

  if (length(parts) == 0) {
    return(1)
  }

  prod(counts[parts])
}


umc_get_component_variance <- function(comp_tbl, component_name) {
  val <- comp_tbl$variance[comp_tbl$component == component_name]
  if (length(val) == 0) 0 else as.numeric(val[1])
}


umc_metrics_from_component_table <- function(comp_tbl, dat, modeled_facets = umc_all_facets()) {
  sigma_tau <- umc_get_component_variance(comp_tbl, "item_id")
  counts <- umc_design_counts(dat)

  sigma_delta <- comp_tbl %>%
    filter(grepl("^item_id:", component)) %>%
    summarise(val = sum(variance / vapply(component, umc_component_divisor, numeric(1), dat = dat))) %>%
    pull(val)
  sigma_delta <- sigma_delta + umc_get_component_variance(comp_tbl, "Residual") / prod(counts)

  sigma_abs_extra <- comp_tbl %>%
    filter(component != "item_id", component != "Residual", !grepl("^item_id:", component)) %>%
    summarise(val = sum(variance / vapply(component, umc_component_divisor, numeric(1), dat = dat))) %>%
    pull(val)

  tibble(
    erho2 = sigma_tau / (sigma_tau + sigma_delta),
    phi = sigma_tau / (sigma_tau + sigma_delta + sigma_abs_extra),
    sigma_tau = sigma_tau,
    sigma_delta = sigma_delta,
    sigma_abs_extra = sigma_abs_extra
  )
}


umc_standardize_components <- function(tbl, method_name, dataset_name) {
  tibble(component = umc_key_components()) %>%
    left_join(tbl %>% select(component, variance), by = "component") %>%
    mutate(
      dataset = dataset_name,
      method = method_name,
      variance = ifelse(is.na(variance), 0, variance)
    )
}


umc_write_outputs <- function(result, prefix) {
  dir.create(UMC_OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
  save(result, file = file.path(UMC_OUTPUT_DIR, paste0(prefix, ".RData")))
  write.csv(result$summary, file.path(UMC_OUTPUT_DIR, paste0(prefix, "_summary.csv")), row.names = FALSE)
  write.csv(result$components, file.path(UMC_OUTPUT_DIR, paste0(prefix, "_components.csv")), row.names = FALSE)
  invisible(result)
}
