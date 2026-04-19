suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(lme4))


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

MG_DEFAULT_MAX_MODELS  <- 8L
MG_DEFAULT_MAX_PROMPTS <- 6L

MG_MH_SEVERITY_MAP <- c(
  NORMAL = 1,
  STRESS = 2,
  ANXIETY = 3,
  DEPRESSION = 4,
  BIPOLAR = 5,
  PERSONALITY_DISORDER = 6,
  SUICIDAL = 7
)

# 3-level collapsed scale (Option B: severity-ordered tertile split)
#   1 = normal  : NORMAL, STRESS          (original levels 1-2; 29+3  = 32 items)
#   2 = mild    : ANXIETY, DEPRESSION     (original levels 3-4; 17+21 = 38 items)
#   3 = severe  : BIPOLAR, PD, SUICIDAL   (original levels 5-7;  4+2+24 = 30 items)
# Motivation: merges sparse categories (BIPOLAR n=4, PD n=2) with adjacent
# conditions while preserving the ordinal severity gradient.
MG_MH_SEVERITY_MAP_3L <- c(
  NORMAL               = 1L,
  STRESS               = 1L,
  ANXIETY              = 2L,
  DEPRESSION           = 2L,
  BIPOLAR              = 3L,
  PERSONALITY_DISORDER = 3L,
  SUICIDAL             = 3L
)


mg_fit_with_warnings <- function(expr) {
  warnings <- character(0)
  value <- withCallingHandlers(
    tryCatch(expr, error = function(e) e),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )

  list(value = value, warnings = unique(warnings))
}


mg_write_csv <- function(tbl, output_prefix, suffix) {
  if (is.null(output_prefix) || is.null(tbl)) {
    return(invisible(NULL))
  }

  out_path <- paste0(output_prefix, "_", suffix, ".csv")
  dir.create(dirname(out_path), showWarnings = FALSE, recursive = TRUE)
  write.csv(tbl, out_path, row.names = FALSE)

  if (!file.exists(out_path)) {
    warning("Failed to write output file: ", out_path)
  }
  invisible(out_path)
}


MG_UNIVARIATE_TASKS <- c(
  "hate_speech",         # expects pre-existing `severity` column
  "drug_review",         # expects pre-existing `severity` column
  "mental_health",       # derives `severity` from `label` via MG_MH_SEVERITY_MAP
  "mental_health_3level" # derives `severity` from `label` via MG_MH_SEVERITY_MAP_3L
)

mg_prepare_univariate_data <- function(task, data_path = NULL, data = NULL) {
  if (!task %in% MG_UNIVARIATE_TASKS) {
    stop(
      "Unsupported univariate task '", task, "'. ",
      "Must be one of: ", paste(MG_UNIVARIATE_TASKS, collapse = ", "), "."
    )
  }

  if (is.null(data)) {
    if (is.null(data_path)) {
      stop("Provide either `data` or `data_path`.")
    }
    data <- read.csv(data_path, stringsAsFactors = FALSE)
  }

  dat <- as_tibble(data)

  if (task == "mental_health") {
    dat$severity <- unname(MG_MH_SEVERITY_MAP[dat$label])
    if (any(is.na(dat$severity))) {
      stop("Mental-health severity mapping produced NA values.")
    }
  } else if (task == "mental_health_3level") {
    dat$severity <- unname(MG_MH_SEVERITY_MAP_3L[dat$label])
    if (any(is.na(dat$severity))) {
      stop("Mental-health 3-level severity mapping produced NA values. ",
           "Unknown labels: ",
           paste(unique(dat$label[is.na(unname(MG_MH_SEVERITY_MAP_3L[dat$label]))]),
                 collapse = ", "))
    }
  } else {
    # hate_speech and drug_review: severity column must already exist in the data
    if (!"severity" %in% names(dat)) {
      stop(
        "Task '", task, "' requires a pre-existing `severity` column in the data, ",
        "but none was found. Available columns: ",
        paste(names(dat), collapse = ", "), "."
      )
    }
  }

  dat %>%
    mutate(
      severity = as.numeric(severity),
      item_id = factor(item_id),
      evaluator = factor(evaluator),
      prompt_type = factor(prompt_type),
      temperature = factor(temperature),
      seed = factor(seed)
    )
}


# ---------------------------------------------------------------------------
# Shared model group specifications (used by both univariate & multivariate)
# ---------------------------------------------------------------------------

#' Returns the 6 nested model structures as named lists of random-effect
#' grouping terms.  Both univariate and multivariate workflows call this so
#' that the two model sets stay in sync.
mg_model_group_specs <- function() {
  list(
    full = c(
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
      "temperature:seed"
    ),
    nested_seed = c(
      "item_id",
      "evaluator",
      "prompt_type",
      "temperature",
      "item_id:evaluator",
      "item_id:prompt_type",
      "item_id:temperature",
      "evaluator:prompt_type",
      "evaluator:temperature",
      "prompt_type:temperature",
      "evaluator:prompt_type:temperature:seed"
    ),
    nested_seed_temperature = c(
      "item_id",
      "evaluator",
      "prompt_type",
      "item_id:evaluator",
      "item_id:prompt_type",
      "evaluator:prompt_type",
      "evaluator:prompt_type:temperature",
      "evaluator:prompt_type:temperature:seed"
    ),
    nested_prompt_temperature_seed = c(
      "item_id",
      "evaluator",
      "item_id:evaluator"
    ),
    prompt_only = c(
      "item_id",
      "prompt_type",
      "item_id:prompt_type"
    ),
    item_only = c("item_id")
  )
}


mg_make_univariate_formula <- function(groups, response = "severity") {
  rhs <- paste(sprintf("(1 | %s)", groups), collapse = " + ")
  as.formula(paste(response, "~ 1 +", rhs))
}


mg_univariate_model_formulas <- function() {
  list(
    full = as.formula(
      paste(
        "severity ~ 1 +",
        "(1 | item_id) +",
        "(1 | evaluator) +",
        "(1 | prompt_type) +",
        "(1 | temperature) +",
        "(1 | seed) +",
        "(1 | item_id:evaluator) +",
        "(1 | item_id:prompt_type) +",
        "(1 | item_id:temperature) +",
        "(1 | item_id:seed) +",
        "(1 | evaluator:prompt_type) +",
        "(1 | evaluator:temperature) +",
        "(1 | evaluator:seed) +",
        "(1 | prompt_type:temperature) +",
        "(1 | prompt_type:seed) +",
        "(1 | temperature:seed)"
      )
    ),
    # Nest seed within all crossed facets (evaluator × prompt × temperature).
    # Seed has no main effect, no 2-way interactions; its variance is captured
    # entirely by the 4-way cell evaluator:prompt_type:temperature:seed.
    # All other facets remain fully crossed.
    nested_seed = as.formula(
      paste(
        "severity ~ 1 +",
        "(1 | item_id) +",
        "(1 | evaluator) +",
        "(1 | prompt_type) +",
        "(1 | temperature) +",
        "(1 | item_id:evaluator) +",
        "(1 | item_id:prompt_type) +",
        "(1 | item_id:temperature) +",
        "(1 | evaluator:prompt_type) +",
        "(1 | evaluator:temperature) +",
        "(1 | prompt_type:temperature) +",
        "(1 | evaluator:prompt_type:temperature:seed)"
      )
    ),
    # Nest both seed and temperature within evaluator × prompt.
    # Temperature has no main effect or separate interactions with item;
    # its variance is captured by evaluator:prompt_type:temperature.
    # Seed is further nested within that: evaluator:prompt_type:temperature:seed.
    nested_seed_temperature = as.formula(
      paste(
        "severity ~ 1 +",
        "(1 | item_id) +",
        "(1 | evaluator) +",
        "(1 | prompt_type) +",
        "(1 | item_id:evaluator) +",
        "(1 | item_id:prompt_type) +",
        "(1 | evaluator:prompt_type) +",
        "(1 | evaluator:prompt_type:temperature) +",
        "(1 | evaluator:prompt_type:temperature:seed)"
      )
    ),
    # Nest prompt, temperature, and seed within evaluator.
    # Only evaluator is a crossed facet; prompt is nested within evaluator,
    # temperature within evaluator:prompt, seed within evaluator:prompt:temperature.
    nested_prompt_temperature_seed = as.formula(
      paste(
        "severity ~ 1 +",
        "(1 | item_id) +",
        "(1 | evaluator) +",
        "(1 | item_id:evaluator) +",
        "(1 | evaluator:prompt_type) +",
        "(1 | evaluator:prompt_type:temperature) +",
        "(1 | evaluator:prompt_type:temperature:seed)"
      )
    ),
    prompt_only = as.formula(
      paste(
        "severity ~ 1 +",
        "(1 | item_id) +",
        "(1 | prompt_type) +",
        "(1 | item_id:prompt_type)"
      )
    ),
    item_only = as.formula("severity ~ 1 + (1 | item_id)")
  )
}


# Shared model-row collector (used by both univariate & multivariate fitting)
mg_collect_model_row <- function(model_name, fit_obj, warnings) {
  if (inherits(fit_obj, "error")) {
    return(tibble(
      model = model_name,
      converged = FALSE,
      singular = NA,
      n_par = NA_integer_,
      minus2LL = NA_real_,
      AIC = NA_real_,
      BIC = NA_real_,
      warning_text = paste(c(warnings, fit_obj$message), collapse = " | ")
    ))
  }

  ll <- logLik(fit_obj)
  tibble(
    model = model_name,
    converged = TRUE,
    singular = lme4::isSingular(fit_obj, tol = 1e-4),
    n_par = attr(ll, "df"),
    minus2LL = -2 * as.numeric(ll),
    AIC = AIC(fit_obj),
    BIC = BIC(fit_obj),
    warning_text = paste(warnings, collapse = " | ")
  )
}


mg_fit_univariate_models <- function(dat, model_formulas = mg_univariate_model_formulas()) {
  control <- lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000))

  fits <- list()
  rows <- list()
  idx <- 1

  for (model_name in names(model_formulas)) {
    cat(sprintf("Fitting univariate %s...\n", model_name))
    fit_res <- mg_fit_with_warnings(
      lmer(model_formulas[[model_name]], data = dat, REML = FALSE, control = control)
    )
    fits[[model_name]] <- fit_res$value
    rows[[idx]] <- mg_collect_model_row(model_name, fit_res$value, fit_res$warnings)
    idx <- idx + 1
  }

  list(
    fits = fits,
    summary = bind_rows(rows) %>% arrange(BIC)
  )
}


mg_select_best_model <- function(summary_tbl) {
  usable <- summary_tbl %>%
    filter(converged, !is.na(BIC)) %>%
    arrange(BIC)

  if (nrow(usable) == 0) {
    stop("No converged models were available for selection.")
  }

  best <- usable[1, , drop = FALSE]

  if (isTRUE(best$singular)) {
    message(
      "Note: best model '", best$model, "' is singular (one or more ",
      "variance components estimated at zero). Estimates are valid but ",
      "confidence intervals may be unreliable."
    )
  }

  best
}


mg_map_univariate_facet_group <- function(component) {
  dplyr::case_when(
    component == "item_id" ~ "item",
    component %in% c("evaluator", "item_id:evaluator") ~ "evaluator",
    component %in% c("prompt_type", "item_id:prompt_type") ~ "prompt",
    component %in% c("temperature", "item_id:temperature") ~ "temperature",
    component %in% c("seed", "item_id:seed") ~ "seed",
    component == "evaluator:prompt_type" ~ "evaluator_prompt",
    component == "evaluator:temperature" ~ "evaluator_temperature",
    component == "evaluator:seed" ~ "evaluator_seed",
    component == "prompt_type:temperature" ~ "prompt_temperature",
    component == "prompt_type:seed" ~ "prompt_seed",
    component == "temperature:seed" ~ "temperature_seed",
    component == "Residual" ~ "residual",
    TRUE ~ "other"
  )
}


mg_map_univariate_component_type <- function(component) {
  dplyr::case_when(
    component == "item_id" ~ "object_of_measurement",
    grepl("^item_id:", component) ~ "item_interaction",
    component == "Residual" ~ "residual",
    TRUE ~ "facet_component"
  )
}


mg_univariate_component_table <- function(fit_obj) {
  vc <- as.data.frame(VarCorr(fit_obj))

  comp_tbl <- bind_rows(
    tibble(
      component = vc$grp[vc$grp != "Residual"],
      variance = as.numeric(vc$vcov[vc$grp != "Residual"])
    ),
    tibble(component = "Residual", variance = sigma(fit_obj)^2)
  ) %>%
    group_by(component) %>%
    summarise(variance = sum(variance), .groups = "drop") %>%
    mutate(
      facet_group = mg_map_univariate_facet_group(component),
      component_type = mg_map_univariate_component_type(component),
      prop_total = variance / sum(variance)
    ) %>%
    arrange(desc(variance))

  comp_tbl
}


mg_univariate_facet_table <- function(component_tbl) {
  component_tbl %>%
    group_by(facet_group) %>%
    summarise(
      variance = sum(variance),
      prop_total = variance / sum(component_tbl$variance),
      .groups = "drop"
    ) %>%
    arrange(desc(variance))
}


mg_get_component_variance <- function(component_tbl, component_name) {
  val <- component_tbl$variance[component_tbl$component == component_name]
  if (length(val) == 0) 0 else as.numeric(val[1])
}


mg_univariate_design_counts <- function(dat) {
  c(
    evaluator = nlevels(dat$evaluator),
    prompt_type = nlevels(dat$prompt_type),
    temperature = nlevels(dat$temperature),
    seed = nlevels(dat$seed)
  )
}


mg_component_divisor <- function(component, dat) {
  counts <- mg_univariate_design_counts(dat)
  parts <- setdiff(strsplit(component, ":", fixed = TRUE)[[1]], "item_id")
  parts <- intersect(parts, names(counts))

  if (length(parts) == 0) {
    return(1)
  }

  prod(counts[parts])
}


mg_predict_univariate_g <- function(
    fit_obj,
    dat,
    n_models,
    n_prompts,
    n_temperatures = nlevels(dat$temperature),
    n_seeds = nlevels(dat$seed)) {
  comp_tbl <- mg_univariate_component_table(fit_obj)

  sigma_tau <- mg_get_component_variance(comp_tbl, "item_id")
  sigma_delta <- comp_tbl %>%
    filter(component_type == "item_interaction") %>%
    summarise(val = sum(variance / vapply(component, mg_component_divisor, numeric(1), dat = dat))) %>%
    pull(val)
  sigma_delta <- sigma_delta + mg_get_component_variance(comp_tbl, "Residual") / (n_models * n_prompts * n_temperatures * n_seeds)

  sigma_tau / (sigma_tau + sigma_delta)
}


# Φ (dependability coefficient) for the observed design.
# Extends Eρ² by adding absolute error from non-item facet main effects and
# cross-facet interactions that don't involve item_id.  These contribute to
# error for absolute (criterion-referenced) decisions but cancel for relative
# (norm-referenced) decisions captured by Eρ².
mg_predict_univariate_phi <- function(
    fit_obj,
    dat,
    n_models       = nlevels(dat$evaluator),
    n_prompts      = nlevels(dat$prompt_type),
    n_temperatures = nlevels(dat$temperature),
    n_seeds        = nlevels(dat$seed)) {
  comp_tbl <- mg_univariate_component_table(fit_obj)

  sigma_tau <- mg_get_component_variance(comp_tbl, "item_id")

  # Relative error (shared with Eρ² formula)
  sigma_delta <- comp_tbl %>%
    filter(component_type == "item_interaction") %>%
    summarise(val = sum(variance / vapply(component, mg_component_divisor, numeric(1), dat = dat))) %>%
    pull(val)
  sigma_delta <- sigma_delta + mg_get_component_variance(comp_tbl, "Residual") /
    (n_models * n_prompts * n_temperatures * n_seeds)

  # Additional absolute error: non-item main effects + cross-facet interactions
  sigma_abs_extra <- comp_tbl %>%
    filter(component_type == "facet_component") %>%
    summarise(val = sum(variance / vapply(component, mg_component_divisor, numeric(1), dat = dat))) %>%
    pull(val)

  denom <- sigma_tau + sigma_delta + sigma_abs_extra
  if (denom == 0) NA_real_ else sigma_tau / denom
}


# ---------------------------------------------------------------------------
# Bootstrap confidence intervals for the observed-design G-coefficient and Φ
# ---------------------------------------------------------------------------

#' Item-level bootstrap (percentile method) for the observed-design Eρ².
#'
#' Resamples items with replacement B times, re-fits the best-model formula
#' via lmer on each resample (preserving the full crossed structure within
#' each selected item), recomputes Eρ² for the observed design, and returns
#' percentile CIs.
#'
#' @param fit_obj  Fitted lmer object (best model on full data).
#' @param dat      Prepared data frame (output of mg_prepare_univariate_data).
#' @param n_models,n_prompts,n_temperatures,n_seeds  Design counts for G.
#'   Defaults to the observed levels in dat.
#' @param B        Number of bootstrap resamples (default 2000).
#' @param seed     RNG seed for reproducibility (default 42).
#' @param conf     Confidence level (default 0.95).
#'
#' @return A one-row tibble with observed_g, ci_lower, ci_upper,
#'   observed_phi, phi_ci_lower, phi_ci_upper, conf_level,
#'   n_boot, n_used, n_error, n_singular, n_nonconverged, ci_status, seed.
mg_bootstrap_g_ci <- function(
    fit_obj,
    dat,
    n_models     = nlevels(dat$evaluator),
    n_prompts    = nlevels(dat$prompt_type),
    n_temperatures = nlevels(dat$temperature),
    n_seeds      = nlevels(dat$seed),
    B            = 2000L,
    seed         = 42L,
    conf         = 0.95,
    min_valid    = 100L) {

  set.seed(seed)

  items   <- levels(dat$item_id)
  n_items <- length(items)
  formula_b <- formula(fit_obj)
  control   <- lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000))

  # Point estimates from the original fit
  observed_g   <- mg_predict_univariate_g(
    fit_obj, dat, n_models, n_prompts, n_temperatures, n_seeds
  )
  observed_phi <- mg_predict_univariate_phi(
    fit_obj, dat, n_models, n_prompts, n_temperatures, n_seeds
  )

  boot_g    <- numeric(B)
  boot_phi  <- numeric(B)
  n_error   <- 0L
  n_singular <- 0L
  n_nonconverged <- 0L
  used_mask <- rep(FALSE, B)

  classify_convergence_warning <- function(warnings) {
    if (length(warnings) == 0) {
      return(FALSE)
    }

    any(grepl(
      "failed to converge|degenerate Hessian|unable to evaluate scaled gradient|false convergence|convergence code",
      warnings,
      ignore.case = TRUE
    ))
  }

  for (b in seq_len(B)) {
    if (b %% 200 == 0) {
      cat(sprintf("  Bootstrap: %d / %d\n", b, B))
    }

    # Resample items with replacement; duplicate items get new unique ids
    # so that lmer treats each copy as a separate (but structurally identical)
    # item in the random-effects grouping.
    selected <- sample(items, n_items, replace = TRUE)
    boot_dat  <- bind_rows(lapply(seq_along(selected), function(k) {
      rows           <- dat[dat$item_id == selected[k], , drop = FALSE]
      rows$item_id   <- factor(paste0("i", k))
      rows
    }))

    fit_res <- mg_fit_with_warnings(
      suppressMessages(
        lmer(formula_b, data = boot_dat, REML = FALSE, control = control)
      )
    )

    fit_b <- fit_res$value
    if (inherits(fit_b, "error")) {
      boot_g[b]  <- NA_real_
      n_error    <- n_error + 1L
      next
    }

    if (classify_convergence_warning(fit_res$warnings)) {
      boot_g[b] <- NA_real_
      n_nonconverged <- n_nonconverged + 1L
      next
    }

    comp_b  <- mg_univariate_component_table(fit_b)
    if (mg_get_component_variance(comp_b, "item_id") < 1e-6) {
      boot_g[b] <- NA_real_
      n_singular <- n_singular + 1L
      next
    }

    s_tau   <- mg_get_component_variance(comp_b, "item_id")
    s_delta <- comp_b %>%
      filter(component_type == "item_interaction") %>%
      summarise(val = sum(variance / vapply(component, mg_component_divisor, numeric(1), dat = dat))) %>%
      pull(val)
    s_delta <- s_delta + mg_get_component_variance(comp_b, "Residual") /
      (n_models * n_prompts * n_temperatures * n_seeds)

    boot_g[b] <- if ((s_tau + s_delta) == 0) NA_real_ else
      s_tau / (s_tau + s_delta)

    # Φ: add absolute error from non-item main effects + cross-facet interactions
    s_abs_extra <- comp_b %>%
      filter(component_type == "facet_component") %>%
      summarise(val = sum(variance / vapply(component, mg_component_divisor, numeric(1), dat = dat))) %>%
      pull(val)

    phi_denom   <- s_tau + s_delta + s_abs_extra
    boot_phi[b] <- if (phi_denom == 0) NA_real_ else s_tau / phi_denom

    used_mask[b] <- !is.na(boot_g[b])
  }

  n_used <- sum(used_mask)
  ci_status <- "ok"

  if (n_error > 0L) {
    message(sprintf(
      "Bootstrap: %d / %d resamples errored and were excluded.",
      n_error, B
    ))
  }
  if (n_singular > 0L) {
    message(sprintf(
      "Bootstrap: %d / %d resamples were singular and were excluded.",
      n_singular, B
    ))
  }
  if (n_nonconverged > 0L) {
    message(sprintf(
      "Bootstrap: %d / %d resamples had convergence warnings and were excluded.",
      n_nonconverged, B
    ))
  }
  if (n_used < B * 0.8) {
    warning(sprintf(
      "Bootstrap: only %d / %d resamples were usable (<80%%). CIs may be unreliable.",
      n_used, B
    ))
    ci_status <- "low_success_rate"
  }

  probs <- c((1 - conf) / 2, 1 - (1 - conf) / 2)
  if (n_used < min_valid) {
    warning(sprintf(
      "Bootstrap: only %d usable resamples (< min_valid = %d). CI not reported.",
      n_used, min_valid
    ))
    ci_g   <- c(NA_real_, NA_real_)
    ci_phi <- c(NA_real_, NA_real_)
    ci_status <- "insufficient_valid_resamples"
  } else {
    ci_g   <- quantile(boot_g[used_mask],   probs = probs, na.rm = TRUE)
    ci_phi <- quantile(boot_phi[used_mask], probs = probs, na.rm = TRUE)
  }

  tibble(
    observed_g     = observed_g,
    ci_lower       = unname(ci_g[1]),
    ci_upper       = unname(ci_g[2]),
    observed_phi   = observed_phi,
    phi_ci_lower   = unname(ci_phi[1]),
    phi_ci_upper   = unname(ci_phi[2]),
    conf_level     = conf,
    n_boot         = as.integer(B),
    n_used         = as.integer(n_used),
    n_error        = as.integer(n_error),
    n_singular     = as.integer(n_singular),
    n_nonconverged = as.integer(n_nonconverged),
    ci_status      = ci_status,
    seed           = as.integer(seed)
  )
}


# D-study prediction grid: varies n_models and n_prompts while holding
# n_temperatures and n_seeds fixed at observed counts.  Temperature and seed
# are treated as replication facets whose levels are cheap to add and should
# be retained at their observed design size.
mg_univariate_prediction_grid <- function(
    fit_obj,
    dat,
    model_range = seq_len(max(MG_DEFAULT_MAX_MODELS, nlevels(dat$evaluator))),
    prompt_range = seq_len(max(MG_DEFAULT_MAX_PROMPTS, nlevels(dat$prompt_type))),
    target_g = c(0.80, 0.85, 0.90, 0.95)) {
  n_temp <- nlevels(dat$temperature)
  n_seed <- nlevels(dat$seed)

  grid <- expand.grid(
    n_models = model_range,
    n_prompts = prompt_range,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  ) %>%
    as_tibble() %>%
    mutate(
      n_temperatures = n_temp,
      n_seeds = n_seed,
      predicted_g = vapply(
        seq_len(n()),
        function(i) {
          mg_predict_univariate_g(
            fit_obj = fit_obj,
            dat = dat,
            n_models = n_models[[i]],
            n_prompts = n_prompts[[i]],
            n_temperatures = n_temperatures[[i]],
            n_seeds = n_seeds[[i]]
          )
        },
        numeric(1)
      ),
      total_cells = n_models * n_prompts * n_temperatures * n_seeds
    )

  targets <- bind_rows(lapply(target_g, function(g_cut) {
    eligible <- grid %>%
      filter(predicted_g >= g_cut) %>%
      arrange(total_cells, n_models, n_prompts)

    if (nrow(eligible) == 0) {
      return(tibble(
        target_g = g_cut,
        n_models = NA_integer_,
        n_prompts = NA_integer_,
        n_temperatures = n_temp,
        n_seeds = n_seed,
        predicted_g = NA_real_,
        total_cells = NA_integer_
      ))
    }

    eligible[1, ] %>%
      mutate(target_g = g_cut) %>%
      select(target_g, everything())
  }))

  list(grid = grid, targets = targets)
}


mg_univariate_prediction_points <- function(
    fit_obj,
    dat,
    target_tbl = NULL,
    target_g = c(0.80, 0.85, 0.90, 0.95)) {
  baseline_tbl <- tibble(
    scenario = "baseline_1x1x1x1",
    target_g = NA_real_,
    n_models = 1L,
    n_prompts = 1L,
    n_temperatures = 1L,
    n_seeds = 1L,
    predicted_g = mg_predict_univariate_g(
      fit_obj = fit_obj,
      dat = dat,
      n_models = 1L,
      n_prompts = 1L,
      n_temperatures = 1L,
      n_seeds = 1L
    )
  )

  observed_tbl <- tibble(
    scenario = "observed_design",
    target_g = NA_real_,
    n_models = nlevels(dat$evaluator),
    n_prompts = nlevels(dat$prompt_type),
    n_temperatures = nlevels(dat$temperature),
    n_seeds = nlevels(dat$seed),
    predicted_g = mg_predict_univariate_g(
      fit_obj = fit_obj,
      dat = dat,
      n_models = nlevels(dat$evaluator),
      n_prompts = nlevels(dat$prompt_type),
      n_temperatures = nlevels(dat$temperature),
      n_seeds = nlevels(dat$seed)
    )
  )

  if (is.null(target_tbl)) {
    target_tbl <- mg_univariate_prediction_grid(
      fit_obj = fit_obj,
      dat = dat,
      target_g = target_g
    )$targets
  }

  target_points <- target_tbl %>%
    transmute(
      scenario = paste0("minimum_design_for_G_", format(target_g, nsmall = 2)),
      target_g,
      n_models,
      n_prompts,
      n_temperatures,
      n_seeds,
      predicted_g
    )

  bind_rows(baseline_tbl, observed_tbl, target_points)
}


mg_run_univariate_manuscript <- function(
    task,
    data_path = NULL,
    data = NULL,
    output_prefix = NULL,
    model_formulas = mg_univariate_model_formulas(),
    target_g = c(0.80, 0.85, 0.90, 0.95),
    model_range = NULL,
    prompt_range = NULL,
    n_boot = 2000L,
    boot_seed = 42L) {
  dat <- mg_prepare_univariate_data(task = task, data_path = data_path, data = data)
  fit_bundle <- mg_fit_univariate_models(dat = dat, model_formulas = model_formulas)
  best_row <- mg_select_best_model(fit_bundle$summary)
  best_fit <- fit_bundle$fits[[best_row$model]]

  component_tbl <- mg_univariate_component_table(best_fit)
  facet_tbl <- mg_univariate_facet_table(component_tbl)
  prediction <- mg_univariate_prediction_grid(
    fit_obj = best_fit,
    dat = dat,
    model_range = if (is.null(model_range)) seq_len(max(MG_DEFAULT_MAX_MODELS, nlevels(dat$evaluator))) else model_range,
    prompt_range = if (is.null(prompt_range)) seq_len(max(MG_DEFAULT_MAX_PROMPTS, nlevels(dat$prompt_type))) else prompt_range,
    target_g = target_g
  )

  prediction_points <- mg_univariate_prediction_points(
    fit_obj = best_fit,
    dat = dat,
    target_tbl = prediction$targets,
    target_g = target_g
  )

  observed_g <- mg_predict_univariate_g(
    fit_obj = best_fit,
    dat = dat,
    n_models = nlevels(dat$evaluator),
    n_prompts = nlevels(dat$prompt_type)
  )

  observed_phi <- mg_predict_univariate_phi(
    fit_obj = best_fit,
    dat = dat,
    n_models = nlevels(dat$evaluator),
    n_prompts = nlevels(dat$prompt_type)
  )

  # Bootstrap CI for the observed-design G-coefficient and Φ
  boot_ci <- NULL
  if (n_boot > 0L) {
    cat(sprintf("Running %d bootstrap resamples for %s...\n", n_boot, task))
    boot_ci <- mg_bootstrap_g_ci(
      fit_obj       = best_fit,
      dat           = dat,
      n_models      = nlevels(dat$evaluator),
      n_prompts     = nlevels(dat$prompt_type),
      n_temperatures = nlevels(dat$temperature),
      n_seeds       = nlevels(dat$seed),
      B             = n_boot,
      seed          = boot_seed
    )
    cat(sprintf(
      "  %s: Eρ² = %.3f (%.0f%% CI: %.3f, %.3f)  Φ = %.3f (%.0f%% CI: %.3f, %.3f)\n",
      task,
      boot_ci$observed_g,
      boot_ci$conf_level * 100,
      boot_ci$ci_lower,
      boot_ci$ci_upper,
      boot_ci$observed_phi,
      boot_ci$conf_level * 100,
      boot_ci$phi_ci_lower,
      boot_ci$phi_ci_upper
    ))
  }

  primary_source <- facet_tbl %>%
    filter(!facet_group %in% c("item")) %>%
    arrange(desc(variance)) %>%
    slice(1)

  primary_source_name <- if (nrow(primary_source) == 0) NA_character_ else primary_source$facet_group[[1]]
  primary_source_prop <- if (nrow(primary_source) == 0) NA_real_ else primary_source$prop_total[[1]]

  # Identify variance components estimated at zero (boundary estimates)
  boundary_comps <- component_tbl$component[component_tbl$variance == 0]
  boundary_str <- if (length(boundary_comps) == 0) NA_character_ else paste(boundary_comps, collapse = "; ")

  overview_tbl <- tibble(
    task = task,
    analysis = "univariate",
    best_model = best_row$model,
    best_model_BIC = best_row$BIC,
    best_model_singular = best_row$singular,
    boundary_components = boundary_str,
    observed_items = nlevels(dat$item_id),
    observed_models = nlevels(dat$evaluator),
    observed_prompts = nlevels(dat$prompt_type),
    observed_temperatures = nlevels(dat$temperature),
    observed_seeds = nlevels(dat$seed),
    observed_design_g = observed_g,
    g_ci_lower = if (!is.null(boot_ci)) boot_ci$ci_lower else NA_real_,
    g_ci_upper = if (!is.null(boot_ci)) boot_ci$ci_upper else NA_real_,
    g_ci_conf  = if (!is.null(boot_ci)) boot_ci$conf_level else NA_real_,
    g_ci_status = if (!is.null(boot_ci)) boot_ci$ci_status else NA_character_,
    observed_design_phi = observed_phi,
    phi_ci_lower = if (!is.null(boot_ci)) boot_ci$phi_ci_lower else NA_real_,
    phi_ci_upper = if (!is.null(boot_ci)) boot_ci$phi_ci_upper else NA_real_,
    bootstrap_n_used = if (!is.null(boot_ci)) boot_ci$n_used else NA_integer_,
    bootstrap_n_error = if (!is.null(boot_ci)) boot_ci$n_error else NA_integer_,
    bootstrap_n_singular = if (!is.null(boot_ci)) boot_ci$n_singular else NA_integer_,
    bootstrap_n_nonconverged = if (!is.null(boot_ci)) boot_ci$n_nonconverged else NA_integer_,
    primary_source = primary_source_name,
    primary_source_prop_total = primary_source_prop
  )

  mg_write_csv(fit_bundle$summary, output_prefix, "model_bic")
  mg_write_csv(component_tbl, output_prefix, "best_model_components")
  mg_write_csv(facet_tbl, output_prefix, "best_model_facets")
  mg_write_csv(prediction$grid, output_prefix, "prediction_grid")
  mg_write_csv(prediction$targets, output_prefix, "prediction_targets")
  mg_write_csv(prediction_points, output_prefix, "prediction_points")
  mg_write_csv(overview_tbl, output_prefix, "overview")
  mg_write_csv(boot_ci, output_prefix, "bootstrap_ci")

  list(
    data = dat,
    fits = fit_bundle$fits,
    model_table = fit_bundle$summary,
    best_model = best_row,
    best_fit = best_fit,
    component_table = component_tbl,
    facet_table = facet_tbl,
    prediction_grid = prediction$grid,
    prediction_targets = prediction$targets,
    prediction_points = prediction_points,
    bootstrap_ci = boot_ci,
    overview = overview_tbl
  )
}
