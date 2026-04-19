# ---------------------------------------------------------------------------
# bw_models.R
#
# OpenMx Between-Within (RAM/path) specification for nested G-theory models.
# Each variance component is a between-level submodel with a latent variable,
# connected to the observation level via primaryKey/joinKey paths.
#
# Two models (nested facets only — full crossed model uses wide_models.R):
#   1. nested_seed             — seed nested in e×p×t   (12 components)
#   2. nested_seed_temperature — seed & temp nested     (9 components)
#
# Expected working directory: 01. Research/
# ---------------------------------------------------------------------------

suppressPackageStartupMessages(library(OpenMx))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(lme4))

source("univariate/src/shared/univariate_workflow.R")
source("univariate/src/shared/openmx_threads.R")

mxOption(model = NULL, key = "Default optimizer", "CSOLNP", reset = FALSE)


bw_resolve_num_threads <- function(num_threads = NULL) {
  omx_resolve_num_threads(num_threads)
}


bw_configure_openmx <- function(num_threads = NULL, parallel_diagnostics = FALSE) {
  omx_configure_threads(
    num_threads = num_threads,
    parallel_diagnostics = parallel_diagnostics,
    log_fn = bw_log
  )
}


# ---------------------------------------------------------------------------
# Data preparation: add composite key columns needed by all models
# ---------------------------------------------------------------------------

bw_prepare_data <- function(task, data_path) {
  dat <- mg_prepare_univariate_data(task = task, data_path = data_path)
  dat %>%
    mutate(
      item_eval      = factor(paste(item_id, evaluator, sep = ":")),
      item_prompt    = factor(paste(item_id, prompt_type, sep = ":")),
      item_temp      = factor(paste(item_id, temperature, sep = ":")),
      item_seed      = factor(paste(item_id, seed, sep = ":")),
      eval_prompt    = factor(paste(evaluator, prompt_type, sep = ":")),
      eval_temp      = factor(paste(evaluator, temperature, sep = ":")),
      eval_seed      = factor(paste(evaluator, seed, sep = ":")),
      prompt_temp    = factor(paste(prompt_type, temperature, sep = ":")),
      prompt_seed    = factor(paste(prompt_type, seed, sep = ":")),
      temp_seed      = factor(paste(temperature, seed, sep = ":")),
      ept_cell       = factor(paste(evaluator, prompt_type, temperature, sep = ":")),
      epts_cell      = factor(paste(evaluator, prompt_type, temperature, seed, sep = ":"))
    )
}


bw_subset_items <- function(dat, max_items = NULL) {
  if (is.null(max_items)) {
    return(dat)
  }

  all_items <- levels(dat$item_id)
  if (length(all_items) <= max_items) {
    return(dat)
  }

  keep_items <- all_items[seq_len(max_items)]
  dat %>% filter(item_id %in% keep_items) %>% droplevels()
}


bw_log <- function(...) {
  cat(...)
  flush.console()
}


bw_timed <- function(label, expr) {
  t0 <- proc.time()[["elapsed"]]
  value <- force(expr)
  elapsed <- proc.time()[["elapsed"]] - t0
  bw_log(sprintf("%s completed in %.2fs\n", label, elapsed))
  value
}


bw_design_dims <- function() {
  c("item_id", "evaluator", "prompt_type", "temperature", "seed")
}


bw_apply_mode_matrix <- function(arr, mat, dim_idx) {
  perm <- c(dim_idx, setdiff(seq_along(dim(arr)), dim_idx))
  arr_perm <- aperm(arr, perm)
  arr_mat <- matrix(arr_perm, nrow = dim(arr_perm)[1])
  res_mat <- mat %*% arr_mat
  res_arr <- array(res_mat, dim = c(nrow(mat), dim(arr_perm)[-1]))
  aperm(res_arr, order(perm))
}


bw_project_balanced_component <- function(arr, subset_dims, all_dims) {
  out <- arr
  dim_sizes <- setNames(dim(arr), all_dims)

  for (d in all_dims) {
    n_d <- dim_sizes[[d]]
    P_d <- matrix(1 / n_d, n_d, n_d)
    H_d <- diag(n_d) - P_d
    projector <- if (d %in% subset_dims) H_d else P_d
    out <- bw_apply_mode_matrix(out, projector, match(d, all_dims))
  }

  out
}


bw_full_component_dims <- function() {
  list(
    item_id = "item_id",
    evaluator = "evaluator",
    prompt_type = "prompt_type",
    temperature = "temperature",
    seed = "seed",
    "item_id:evaluator" = c("item_id", "evaluator"),
    "item_id:prompt_type" = c("item_id", "prompt_type"),
    "item_id:temperature" = c("item_id", "temperature"),
    "item_id:seed" = c("item_id", "seed"),
    "evaluator:prompt_type" = c("evaluator", "prompt_type"),
    "evaluator:temperature" = c("evaluator", "temperature"),
    "evaluator:seed" = c("evaluator", "seed"),
    "prompt_type:temperature" = c("prompt_type", "temperature"),
    "prompt_type:seed" = c("prompt_type", "seed"),
    "temperature:seed" = c("temperature", "seed")
  )
}


bw_calibrate_start_values <- function(starts,
                                      total_variance,
                                      reserve_frac = 0.05,
                                      component_floor = 1e-5,
                                      residual_floor = 1e-4) {
  total_variance <- as.numeric(total_variance[[1]])
  if (!is.finite(total_variance) || is.na(total_variance) || total_variance <= 0) {
    total_variance <- 1e-3
  }

  component_names <- setdiff(names(starts), "Residual")
  component_vals <- vapply(component_names, function(nm) {
    val <- as.numeric(starts[[nm]][[1]])
    if (!is.finite(val) || is.na(val)) {
      val <- 0
    }
    max(val, component_floor)
  }, numeric(1))

  max_component_total <- max(total_variance * (1 - reserve_frac), component_floor)
  if (length(component_vals) > 0 && sum(component_vals) > max_component_total) {
    component_vals <- component_vals * (max_component_total / sum(component_vals))
  }

  calibrated <- as.list(component_vals)
  residual_val <- max(
    total_variance - sum(component_vals),
    total_variance * reserve_frac,
    residual_floor
  )
  calibrated[["Residual"]] <- residual_val
  calibrated
}


bw_generic_start_values <- function(dat) {
  total_var <- var(as.numeric(dat$severity))
  if (!is.finite(total_var) || is.na(total_var) || total_var <= 0) {
    total_var <- 1
  }

  component_names <- c(names(bw_full_component_dims()), "Residual")
  raw_starts <- setNames(as.list(rep(total_var / length(component_names), length(component_names))),
                         component_names)
  bw_calibrate_start_values(raw_starts, total_variance = total_var)
}


bw_moment_start_values <- function(dat) {
  dims <- bw_design_dims()
  arr <- with(dat, tapply(severity, dat[dims], mean))

  if (is.null(arr) || anyNA(arr)) {
    bw_log("  Moment-based starts unavailable (missing design cells), using generic starts\n")
    return(bw_generic_start_values(dat))
  }

  centered_arr <- arr - mean(arr, na.rm = TRUE)
  total_variance <- mean(centered_arr^2, na.rm = TRUE)
  component_dims <- bw_full_component_dims()

  starts <- list()
  modeled_sum <- array(0, dim = dim(centered_arr), dimnames = dimnames(centered_arr))
  for (nm in names(component_dims)) {
    effect_arr <- bw_project_balanced_component(centered_arr, component_dims[[nm]], dims)
    starts[[nm]] <- mean(effect_arr^2, na.rm = TRUE)
    modeled_sum <- modeled_sum + effect_arr
  }

  resid_arr <- centered_arr - modeled_sum
  starts[["Residual"]] <- mean(resid_arr^2, na.rm = TRUE)
  bw_calibrate_start_values(starts, total_variance = total_variance)
}


bw_compute_start_values <- function(dat, method = c("moments", "lme4", "generic")) {
  method <- match.arg(method)

  if (method == "moments") {
    bw_log("  Computing moment-based start values ... ")
    starts <- bw_moment_start_values(dat)
    bw_log("done\n")
    return(starts)
  }

  if (method == "lme4") {
    return(bw_lme4_start_values(dat))
  }

  bw_log("  Using generic calibrated starts\n")
  bw_generic_start_values(dat)
}


# ---------------------------------------------------------------------------
# Model component definitions
#
# Each component is a list with:
#   name      — OpenMx submodel name
#   key       — column in data used as primaryKey / joinKey
#   lme4_name — equivalent random effect grouping term
# ---------------------------------------------------------------------------

# Model 1: Nested seed (12 components)
# Crossed: evaluator, prompt, temperature.
# Nested: seed within evaluator/prompt/temperature.
bw_nested_seed_components <- function() {
  list(
    list(name = "item",        key = "item_id",       lme4_name = "item_id"),
    list(name = "eval",        key = "evaluator",     lme4_name = "evaluator"),
    list(name = "prompt",      key = "prompt_type",    lme4_name = "prompt_type"),
    list(name = "temp",        key = "temperature",    lme4_name = "temperature"),
    list(name = "item_eval",   key = "item_eval",      lme4_name = "item_id:evaluator"),
    list(name = "item_prompt", key = "item_prompt",     lme4_name = "item_id:prompt_type"),
    list(name = "item_temp",   key = "item_temp",       lme4_name = "item_id:temperature"),
    list(name = "eval_prompt", key = "eval_prompt",     lme4_name = "evaluator:prompt_type"),
    list(name = "eval_temp",   key = "eval_temp",       lme4_name = "evaluator:temperature"),
    list(name = "prompt_temp", key = "prompt_temp",     lme4_name = "prompt_type:temperature"),
    list(name = "seed_in_ept", key = "epts_cell",       lme4_name = "seed:evaluator:prompt_type:temperature")
  )
}

# Model 3: Nested seed + temperature (9 components)
# Crossed: evaluator, prompt.
# Nested: temperature within evaluator/prompt; seed within evaluator/prompt/temperature.
bw_nested_seed_temp_components <- function() {
  list(
    list(name = "item",        key = "item_id",        lme4_name = "item_id"),
    list(name = "eval",        key = "evaluator",      lme4_name = "evaluator"),
    list(name = "prompt",      key = "prompt_type",     lme4_name = "prompt_type"),
    list(name = "item_eval",   key = "item_eval",       lme4_name = "item_id:evaluator"),
    list(name = "item_prompt", key = "item_prompt",      lme4_name = "item_id:prompt_type"),
    list(name = "eval_prompt", key = "eval_prompt",      lme4_name = "evaluator:prompt_type"),
    list(name = "temp_in_ep",  key = "ept_cell",        lme4_name = "temperature:evaluator:prompt_type"),
    list(name = "seed_in_ept", key = "epts_cell",        lme4_name = "seed:evaluator:prompt_type:temperature")
  )
}


# ---------------------------------------------------------------------------
# Compute starting values from lme4 full model
# ---------------------------------------------------------------------------

bw_lme4_start_values <- function(dat) {
  formula_full <- as.formula(paste(
    "severity ~ 1 +",
    "(1|item_id) + (1|evaluator) + (1|prompt_type) +",
    "(1|temperature) + (1|seed) +",
    "(1|item_id:evaluator) + (1|item_id:prompt_type) +",
    "(1|item_id:temperature) + (1|item_id:seed) +",
    "(1|evaluator:prompt_type) + (1|evaluator:temperature) +",
    "(1|evaluator:seed) + (1|prompt_type:temperature) +",
    "(1|prompt_type:seed) + (1|temperature:seed)"
  ))

  cat("  Computing lme4 start values ... ")
  fit <- tryCatch(
    lmer(formula_full, data = dat, REML = FALSE,
         control = lmerControl(optimizer = "bobyqa",
                               optCtrl = list(maxfun = 100000))),
    error = function(e) NULL
  )

  if (is.null(fit)) {
    cat("lme4 failed, using generic starts\n")
    total_var <- var(as.numeric(dat$severity))
    generic <- total_var / 16
    nms <- c("item_id", "evaluator", "prompt_type", "temperature", "seed",
             "item_id:evaluator", "item_id:prompt_type",
             "item_id:temperature", "item_id:seed",
             "evaluator:prompt_type", "evaluator:temperature",
             "evaluator:seed", "prompt_type:temperature",
             "prompt_type:seed", "temperature:seed", "Residual")
    return(setNames(rep(generic, length(nms)), nms))
  }

  vc <- as.data.frame(VarCorr(fit))
  sv <- setNames(vc$vcov, vc$grp)
  cat("done\n")
  sv
}

# Map lme4 full-model variances to starting values for each model type.
# For nested terms, sum only the facet-level crossed components that collapse
# into the nested term. Item-by-facet terms are dropped in the nested models,
# not re-expressed as nested facet variance, so they are excluded here.
bw_map_starts <- function(lme4_sv, model_type) {
  resid_var <- max(lme4_sv[["Residual"]] %||% 0.1, 1e-6)
  eps <- resid_var * 1e-4
  g <- function(name) max(lme4_sv[[name]] %||% 0, eps)

  if (model_type == "full") {
    list(
      item        = g("item_id"),
      eval        = g("evaluator"),
      prompt      = g("prompt_type"),
      temp        = g("temperature"),
      seed        = g("seed"),
      item_eval   = g("item_id:evaluator"),
      item_prompt = g("item_id:prompt_type"),
      item_temp   = g("item_id:temperature"),
      item_seed   = g("item_id:seed"),
      eval_prompt = g("evaluator:prompt_type"),
      eval_temp   = g("evaluator:temperature"),
      eval_seed   = g("evaluator:seed"),
      prompt_temp = g("prompt_type:temperature"),
      prompt_seed = g("prompt_type:seed"),
      temp_seed   = g("temperature:seed"),
      residual    = g("Residual")
    )
  } else if (model_type == "nested_seed") {
    list(
      item        = g("item_id"),
      eval        = g("evaluator"),
      prompt      = g("prompt_type"),
      temp        = g("temperature"),
      item_eval   = g("item_id:evaluator"),
      item_prompt = g("item_id:prompt_type"),
      item_temp   = g("item_id:temperature"),
      eval_prompt = g("evaluator:prompt_type"),
      eval_temp   = g("evaluator:temperature"),
      prompt_temp = g("prompt_type:temperature"),
      # seed nested in ept: absorbs seed + facet-level seed interactions only
      seed_in_ept = g("seed") + g("evaluator:seed") + g("prompt_type:seed") +
                    g("temperature:seed"),
      residual    = g("Residual")
    )
  } else if (model_type == "nested_seed_temperature") {
    list(
      item        = g("item_id"),
      eval        = g("evaluator"),
      prompt      = g("prompt_type"),
      item_eval   = g("item_id:evaluator"),
      item_prompt = g("item_id:prompt_type"),
      eval_prompt = g("evaluator:prompt_type"),
      # temp nested in ep: absorbs temp + facet-level temp interactions only
      temp_in_ep  = g("temperature") + g("evaluator:temperature") +
                    g("prompt_type:temperature"),
      # seed nested in ept: absorbs seed + facet-level seed interactions only
      seed_in_ept = g("seed") + g("evaluator:seed") + g("prompt_type:seed") +
                    g("temperature:seed"),
      residual    = g("Residual")
    )
  } else {
    stop("Unknown model_type: ", model_type)
  }
}


# ---------------------------------------------------------------------------
# Build a single between-level submodel (one variance component)
# ---------------------------------------------------------------------------

bw_make_submodel <- function(comp_name, key_col, dat, start_var) {
  # Use levels from the observation data to ensure factor consistency
  key_vals <- levels(dat[[key_col]])
  key_data <- data.frame(key = factor(key_vals, levels = key_vals))
  names(key_data) <- key_col

  lat_name <- paste0("lat_", comp_name)
  var_label <- paste0("sigma2_", comp_name)

  mxModel(
    name = comp_name,
    type = "RAM",
    latentVars = lat_name,
    mxData(key_data, type = "raw", primaryKey = key_col),
    mxPath(from = lat_name, arrows = 2,
           free = TRUE, values = start_var, lbound = 0,
           labels = var_label)
  )
}


# ---------------------------------------------------------------------------
# Build the observation-level model with all submodels
# ---------------------------------------------------------------------------

bw_build_model <- function(dat, components, starts, model_name = "gtheory_bw") {
  key_cols <- vapply(components, function(comp) comp$key, character(1))
  obs_data <- dat %>%
    mutate(severity = as.numeric(severity)) %>%
    mutate(across(all_of(key_cols), factor)) %>%
    as.data.frame()

  submodels <- lapply(components, function(comp) {
    sv <- starts[[comp$name]] %||% 0.1
    bw_make_submodel(comp$name, comp$key, obs_data, start_var = sv)
  })

  join_paths <- lapply(components, function(comp) {
    lat_name <- paste0(comp$name, ".lat_", comp$name)
    mxPath(from = lat_name, to = "severity",
           free = FALSE, values = 1, joinKey = comp$key)
  })

  mxModel(
    name = model_name, type = "RAM",
    manifestVars = "severity",
    submodels,
    mxData(obs_data, type = "raw"),
    mxPath(from = "one", to = "severity", free = TRUE,
           values = mean(obs_data$severity)),
    mxPath(from = "severity", arrows = 2,
           free = TRUE, values = starts[["residual"]] %||% 0.1,
           lbound = 1e-6, labels = "sigma2_residual"),
    join_paths
  )
}


# ---------------------------------------------------------------------------
# Extract variance components from a fitted model
# ---------------------------------------------------------------------------

bw_extract_components <- function(fit, components) {
  params <- omxGetParameters(fit)

  vc <- list()
  for (comp in components) {
    label <- paste0("sigma2_", comp$name)
    vc[[comp$lme4_name]] <- unname(params[label])
  }
  vc[["Residual"]] <- unname(params["sigma2_residual"])

  tibble(
    component = names(vc),
    variance  = as.numeric(unlist(vc))
  ) %>%
    mutate(prop_total = variance / sum(variance)) %>%
    arrange(desc(variance))
}


# ---------------------------------------------------------------------------
# Fit a single model
# ---------------------------------------------------------------------------

bw_fit_model <- function(dat, model_type = "nested_seed",
                         extra_tries = 20L, max_items = NULL,
                         start_values = NULL) {

  dat <- bw_subset_items(dat, max_items = max_items)

  components <- switch(model_type,
    nested_seed             = bw_nested_seed_components(),
    nested_seed_temperature = bw_nested_seed_temp_components(),
    stop("Unknown model_type: ", model_type)
  )

  # Compute start values if not provided
  if (is.null(start_values)) {
    start_values <- bw_compute_start_values(dat, method = "moments")
  }
  starts <- bw_map_starts(start_values, model_type)

  model_name <- paste0("gtheory_bw_", model_type)
  mx_model <- bw_build_model(dat, components, starts, model_name)

  bw_log(sprintf("Fitting OpenMx BW [%s] — %d components, %d items, %d rows\n",
                 model_type, length(components) + 1L, nlevels(dat$item_id), nrow(dat)))

  if (!is.null(extra_tries)) {
    fit <- bw_timed(
      sprintf("  OpenMx BW optimization [%s]", model_type),
      mxTryHard(mx_model, extraTries = extra_tries,
                OKstatuscodes = 0L,
                jitterDistrib = "runif", loc = 1, scale = 0.25)
    )
  } else {
    fit <- bw_timed(
      sprintf("  OpenMx BW optimization [%s]", model_type),
      mxRun(mx_model)
    )
  }

  status <- fit$output$status$code
  minus2LL <- fit$output$Minus2LogLikelihood
  bw_log(sprintf("  Status: %d | -2LL: %.4f\n", status, minus2LL))

  comp_tbl <- bw_extract_components(fit, components)

  list(
    model_type = model_type,
    fit        = fit,
    status     = status,
    minus2LL   = minus2LL,
    AIC        = AIC(fit),
    BIC        = BIC(fit),
    components = comp_tbl
  )
}


# ---------------------------------------------------------------------------
# Run both nested models on one dataset
# ---------------------------------------------------------------------------

bw_run_all_models <- function(dat, extra_tries = 20L, max_items = NULL,
                              start_method = "moments") {

  # Always compute start values on the full dataset
  start_values <- bw_timed(
    sprintf("Start-value computation [%s]", start_method),
    bw_compute_start_values(dat, method = start_method)
  )

  results <- list()
  for (mt in c("nested_seed", "nested_seed_temperature")) {
    bw_log("\n", strrep("=", 50), "\n", sep = "")
    res <- bw_fit_model(dat, model_type = mt,
                        extra_tries = extra_tries, max_items = max_items,
                        start_values = start_values)
    results[[mt]] <- res

    bw_log("\nVariance components:\n")
    print(res$components, n = Inf, width = Inf)
  }

  # Summary
  bw_log("\n\n", strrep("=", 50), "\n", sep = "")
  bw_log("MODEL COMPARISON\n")
  bw_log(strrep("=", 50), "\n\n", sep = "")
  summary_tbl <- bind_rows(lapply(names(results), function(mt) {
    r <- results[[mt]]
    tibble(model = mt, status = r$status, minus2LL = r$minus2LL,
           AIC = r$AIC, BIC = r$BIC,
           n_components = nrow(r$components))
  }))
  print(summary_tbl, width = Inf)

  results
}
