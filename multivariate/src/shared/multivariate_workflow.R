source("univariate/src/shared/univariate_workflow.R")

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(lme4))


MG_MH_TRAITS <- c(
  depression_present = "depression",
  anxiety_present = "anxiety",
  suicidal_present = "suicidal",
  stress_present = "stress",
  bipolar_present = "bipolar",
  personality_disorder_present = "personality_disorder"
)

MG_DRUG_TRAITS <- c(
  ai_efficacy = "efficacy",
  ai_safety = "safety",
  ai_burden = "burden",
  ai_cost = "cost"
)

MG_DRUG_LEVELS <- c("NEGATIVE", "NEUTRAL", "POSITIVE")
MG_DRUG_SCORE_MAP <- setNames(seq_along(MG_DRUG_LEVELS), MG_DRUG_LEVELS)

# This matches the most defensible 3-condition subset available in the data
# when the manuscript sensitivity analysis is limited to three comorbidities.
MG_MH_SENSITIVITY_TRAITS <- c("depression", "anxiety", "suicidal")


mg_coerce_binary_flag <- function(x) {
  if (is.logical(x) || is.numeric(x)) {
    return(as.integer(x))
  }

  x_chr <- as.character(x)
  as.integer(ifelse(
    x_chr %in% c("True", "TRUE", "1"), 1L,
    ifelse(x_chr %in% c("False", "FALSE", "0"), 0L, NA_integer_)
  ))
}


mg_stack_mental_health_data <- function(dat, trait_subset = NULL) {
  trait_map <- MG_MH_TRAITS

  if (!is.null(trait_subset)) {
    trait_map <- trait_map[trait_map %in% trait_subset]
  }

  missing <- setdiff(names(trait_map), names(dat))
  if (length(missing) > 0) {
    stop("Missing mental-health trait columns: ", paste(missing, collapse = ", "))
  }

  dat %>%
    mutate(across(all_of(names(trait_map)), mg_coerce_binary_flag)) %>%
    pivot_longer(
      cols = all_of(names(trait_map)),
      names_to = "trait_col",
      values_to = "y_num"
    ) %>%
    mutate(
      task = "mental_health",
      trait = factor(unname(trait_map[trait_col]), levels = unname(trait_map)),
      outcome_type = "binary",
      y_num = as.numeric(y_num)
    ) %>%
    select(
      task, item_id, evaluator, prompt_type, temperature, seed,
      trait, outcome_type, y_num, everything()
    )
}


mg_stack_drug_review_data <- function(dat, trait_subset = NULL) {
  trait_map <- MG_DRUG_TRAITS

  if (!is.null(trait_subset)) {
    trait_map <- trait_map[trait_map %in% trait_subset]
  }

  missing <- setdiff(names(trait_map), names(dat))
  if (length(missing) > 0) {
    stop("Missing drug-review trait columns: ", paste(missing, collapse = ", "))
  }

  dat %>%
    mutate(across(all_of(names(trait_map)), as.character)) %>%
    pivot_longer(
      cols = all_of(names(trait_map)),
      names_to = "trait_col",
      values_to = "y_label"
    ) %>%
    mutate(
      task = "drug_review",
      trait = factor(unname(trait_map[trait_col]), levels = unname(trait_map)),
      outcome_type = "ordinal_3",
      y_label = factor(y_label, levels = MG_DRUG_LEVELS, ordered = TRUE),
      y_num = as.numeric(unname(MG_DRUG_SCORE_MAP[as.character(y_label)]))
    ) %>%
    select(
      task, item_id, evaluator, prompt_type, temperature, seed,
      trait, outcome_type, y_num, everything()
    )
}


mg_prepare_multivariate_data <- function(task, data_path = NULL, data = NULL, trait_subset = NULL) {
  if (is.null(data)) {
    if (is.null(data_path)) {
      stop("Provide either `data` or `data_path`.")
    }
    data <- read.csv(data_path, stringsAsFactors = FALSE)
  }

  long_dat <- switch(
    task,
    mental_health = mg_stack_mental_health_data(data, trait_subset = trait_subset),
    drug_review = mg_stack_drug_review_data(data, trait_subset = trait_subset),
    stop("Unsupported multivariate task: ", task)
  )

  long_dat %>%
    mutate(
      item_id = factor(item_id),
      evaluator = factor(evaluator),
      prompt_type = factor(prompt_type),
      temperature = factor(temperature),
      seed = factor(seed),
      trait = factor(trait)
    ) %>%
    arrange(item_id, evaluator, prompt_type, temperature, seed, trait)
}


mg_add_trait_indicators <- function(dat) {
  out <- dat
  for (tr in levels(dat$trait)) {
    col_name <- paste0("trait_", tr)
    out[[col_name]] <- as.numeric(dat$trait == tr)
  }
  out
}


mg_trait_columns <- function(dat) {
  paste0("trait_", levels(dat$trait))
}


# Reuse the shared model group specs defined in univariate_workflow.R.
# By default, the `full` model (all 15 random-effect groups) is excluded
# for multivariate fits because the number of variance parameters scales
# as n_traits x n_groups, making it computationally infeasible for 4+
# traits.  Pass include_full = TRUE to override.
mg_multivariate_model_groups <- function(include_full = FALSE) {
  specs <- mg_model_group_specs()
  if (!include_full) {
    specs[["full"]] <- NULL
  }
  specs
}


# Uses || (double-bar) to suppress cross-trait correlations within each
# random-effect grouping factor.  This is deliberate: G-theory decomposes
# variance components independently per facet, and estimating full
# trait x trait correlation matrices is infeasible with many traits.
mg_make_multivariate_formula <- function(dat, groups, response = "y_num") {
  trait_cols <- mg_trait_columns(dat)
  fixed_part <- paste0("0 + ", paste(trait_cols, collapse = " + "))
  rand_terms <- paste(sprintf("(0 + %s || %s)", paste(trait_cols, collapse = " + "), groups), collapse = " + ")
  as.formula(paste(response, "~", fixed_part, "+", rand_terms))
}


# Row collector reused from univariate_workflow.R: mg_collect_model_row()


mg_fit_multivariate_models <- function(dat, model_groups = mg_multivariate_model_groups()) {
  control <- lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000))
  dat <- mg_add_trait_indicators(dat)

  fits <- list()
  rows <- list()
  idx <- 1

  for (model_name in names(model_groups)) {
    formula_obj <- mg_make_multivariate_formula(dat, groups = model_groups[[model_name]])
    cat(sprintf("Fitting multivariate %s...\n", model_name))
    fit_res <- mg_fit_with_warnings(
      lmer(formula_obj, data = dat, REML = FALSE, control = control)
    )
    fits[[model_name]] <- fit_res$value
    rows[[idx]] <- mg_collect_model_row(model_name, fit_res$value, fit_res$warnings)
    idx <- idx + 1
  }

  list(
    data = dat,
    fits = fits,
    summary = bind_rows(rows) %>% arrange(BIC)
  )
}


mg_map_multivariate_facet_group <- function(component) {
  mg_map_univariate_facet_group(component)
}


mg_multivariate_component_table <- function(fit_obj, traits) {
  vc <- as.data.frame(VarCorr(fit_obj))

  canonicalize_component <- function(x) {
    x <- sub("\\.[0-9]+$", "", x)
    gsub("\\.", ":", x)
  }

  rand_tbl <- vc %>%
    filter(grp != "Residual", !is.na(var1), is.na(var2)) %>%
    transmute(
      trait = sub("^trait_", "", var1),
      component = canonicalize_component(grp),
      variance = as.numeric(vcov)
    )

  # lmer estimates a single pooled residual variance across all traits.
  # We assign it to every trait so that per-trait G-coefficients use a
  # consistent denominator.  This is a modelling assumption, not a bug.
  resid_tbl <- tibble(
    trait = traits,
    component = "Residual",
    variance = rep(sigma(fit_obj)^2, length(traits))
  )

  bind_rows(rand_tbl, resid_tbl) %>%
    mutate(
      facet_group = mg_map_multivariate_facet_group(component),
      component_type = case_when(
        component == "item_id" ~ "object_of_measurement",
        grepl("^item_id:", component) ~ "item_interaction",
        component == "Residual" ~ "residual",
        TRUE ~ "facet_component"
      )
    ) %>%
    group_by(trait) %>%
    mutate(prop_total = variance / sum(variance)) %>%
    ungroup() %>%
    arrange(trait, desc(variance))
}


mg_multivariate_facet_table <- function(component_tbl) {
  component_tbl %>%
    group_by(trait, facet_group) %>%
    summarise(
      variance = sum(variance),
      prop_total = variance / sum(component_tbl$variance[component_tbl$trait == first(trait)]),
      .groups = "drop"
    ) %>%
    arrange(trait, desc(variance))
}


mg_multivariate_mean_facet_table <- function(facet_tbl) {
  facet_tbl %>%
    group_by(facet_group) %>%
    summarise(
      mean_variance = mean(variance),
      mean_prop_total = mean(prop_total),
      .groups = "drop"
    ) %>%
    arrange(desc(mean_variance))
}


mg_get_trait_component_variance <- function(component_tbl, trait_name, component_name) {
  val <- component_tbl$variance[component_tbl$trait == trait_name & component_tbl$component == component_name]
  if (length(val) == 0) 0 else as.numeric(val[1])
}


mg_predict_multivariate_g <- function(
    fit_obj,
    dat,
    n_models,
    n_prompts,
    n_temperatures = nlevels(dat$temperature),
    n_seeds = nlevels(dat$seed)) {
  traits <- levels(dat$trait)
  component_tbl <- mg_multivariate_component_table(fit_obj, traits = traits)

  bind_rows(lapply(traits, function(tr) {
    sigma_tau <- mg_get_trait_component_variance(component_tbl, tr, "item_id")
    sigma_delta <-
      mg_get_trait_component_variance(component_tbl, tr, "item_id:evaluator") / n_models +
      mg_get_trait_component_variance(component_tbl, tr, "item_id:prompt_type") / n_prompts +
      mg_get_trait_component_variance(component_tbl, tr, "item_id:temperature") / n_temperatures +
      mg_get_trait_component_variance(component_tbl, tr, "item_id:seed") / n_seeds +
      mg_get_trait_component_variance(component_tbl, tr, "Residual") /
      (n_models * n_prompts * n_temperatures * n_seeds)

    tibble(
      trait = tr,
      predicted_g = sigma_tau / (sigma_tau + sigma_delta)
    )
  }))
}


mg_multivariate_prediction_grid <- function(
    fit_obj,
    dat,
    model_range = seq_len(max(MG_DEFAULT_MAX_MODELS, nlevels(dat$evaluator))),
    prompt_range = seq_len(max(MG_DEFAULT_MAX_PROMPTS, nlevels(dat$prompt_type))),
    target_g = c(0.80, 0.85, 0.90, 0.95)) {
  n_temp <- nlevels(dat$temperature)
  n_seed <- nlevels(dat$seed)

  grid_rows <- lapply(model_range, function(nm) {
    lapply(prompt_range, function(np) {
      trait_g <- mg_predict_multivariate_g(
        fit_obj = fit_obj,
        dat = dat,
        n_models = nm,
        n_prompts = np,
        n_temperatures = n_temp,
        n_seeds = n_seed
      )

      tibble(
        n_models = nm,
        n_prompts = np,
        n_temperatures = n_temp,
        n_seeds = n_seed,
        mean_g = mean(trait_g$predicted_g, na.rm = TRUE),
        min_g = min(trait_g$predicted_g, na.rm = TRUE),
        max_g = max(trait_g$predicted_g, na.rm = TRUE),
        total_cells = nm * np * n_temp * n_seed
      )
    }) %>%
      bind_rows()
  }) %>%
    bind_rows()

  targets <- bind_rows(lapply(target_g, function(g_cut) {
    bind_rows(lapply(c("mean_g", "min_g"), function(metric_name) {
      eligible <- grid_rows %>%
        filter(.data[[metric_name]] >= g_cut) %>%
        arrange(total_cells, n_models, n_prompts)

      if (nrow(eligible) == 0) {
        return(tibble(
          criterion = metric_name,
          target_g = g_cut,
          n_models = NA_integer_,
          n_prompts = NA_integer_,
          n_temperatures = n_temp,
          n_seeds = n_seed,
          achieved_g = NA_real_,
          total_cells = NA_integer_
        ))
      }

      eligible[1, ] %>%
        transmute(
          criterion = metric_name,
          target_g = g_cut,
          n_models,
          n_prompts,
          n_temperatures,
          n_seeds,
          achieved_g = .data[[metric_name]],
          total_cells
        )
    }))
  }))

  list(grid = grid_rows, targets = targets)
}


mg_multivariate_prediction_points <- function(
    fit_obj,
    dat,
    target_tbl = NULL,
    target_g = c(0.80, 0.85, 0.90, 0.95)) {
  baseline_trait_g <- mg_predict_multivariate_g(
    fit_obj = fit_obj,
    dat = dat,
    n_models = 1L,
    n_prompts = 1L,
    n_temperatures = 1L,
    n_seeds = 1L
  )

  observed_trait_g <- mg_predict_multivariate_g(
    fit_obj = fit_obj,
    dat = dat,
    n_models = nlevels(dat$evaluator),
    n_prompts = nlevels(dat$prompt_type),
    n_temperatures = nlevels(dat$temperature),
    n_seeds = nlevels(dat$seed)
  )

  baseline_tbl <- tibble(
    scenario = "baseline_1x1x1x1",
    criterion = "observed",
    target_g = NA_real_,
    n_models = 1L,
    n_prompts = 1L,
    n_temperatures = 1L,
    n_seeds = 1L,
    achieved_g = mean(baseline_trait_g$predicted_g, na.rm = TRUE),
    min_trait_g = min(baseline_trait_g$predicted_g, na.rm = TRUE)
  )

  observed_tbl <- tibble(
    scenario = "observed_design",
    criterion = "observed",
    target_g = NA_real_,
    n_models = nlevels(dat$evaluator),
    n_prompts = nlevels(dat$prompt_type),
    n_temperatures = nlevels(dat$temperature),
    n_seeds = nlevels(dat$seed),
    achieved_g = mean(observed_trait_g$predicted_g, na.rm = TRUE),
    min_trait_g = min(observed_trait_g$predicted_g, na.rm = TRUE)
  )

  if (is.null(target_tbl)) {
    target_tbl <- mg_multivariate_prediction_grid(
      fit_obj = fit_obj,
      dat = dat,
      target_g = target_g
    )$targets
  }

  target_points <- target_tbl %>%
    transmute(
      scenario = paste0("minimum_design_for_G_", format(target_g, nsmall = 2)),
      criterion,
      target_g,
      n_models,
      n_prompts,
      n_temperatures,
      n_seeds,
      achieved_g,
      min_trait_g = achieved_g
    )

  bind_rows(baseline_tbl, observed_tbl, target_points)
}


# ---------------------------------------------------------------------------
# Composite G-coefficient functions
# ---------------------------------------------------------------------------

# Per-trait Φ (dependability) coefficients at the observed design.
# Analogous to mg_predict_univariate_phi() but operating on the long-format
# multivariate component table.
mg_predict_multivariate_phi <- function(
    fit_obj,
    dat,
    n_models       = nlevels(dat$evaluator),
    n_prompts      = nlevels(dat$prompt_type),
    n_temperatures = nlevels(dat$temperature),
    n_seeds        = nlevels(dat$seed)) {
  traits        <- levels(dat$trait)
  component_tbl <- mg_multivariate_component_table(fit_obj, traits = traits)

  bind_rows(lapply(traits, function(tr) {
    sigma_tau <- mg_get_trait_component_variance(component_tbl, tr, "item_id")

    sigma_delta <-
      mg_get_trait_component_variance(component_tbl, tr, "item_id:evaluator")   / n_models +
      mg_get_trait_component_variance(component_tbl, tr, "item_id:prompt_type") / n_prompts +
      mg_get_trait_component_variance(component_tbl, tr, "item_id:temperature") / n_temperatures +
      mg_get_trait_component_variance(component_tbl, tr, "item_id:seed")        / n_seeds +
      mg_get_trait_component_variance(component_tbl, tr, "Residual") /
        (n_models * n_prompts * n_temperatures * n_seeds)

    sigma_abs_extra <-
      mg_get_trait_component_variance(component_tbl, tr, "evaluator")               / n_models +
      mg_get_trait_component_variance(component_tbl, tr, "prompt_type")             / n_prompts +
      mg_get_trait_component_variance(component_tbl, tr, "temperature")             / n_temperatures +
      mg_get_trait_component_variance(component_tbl, tr, "seed")                    / n_seeds +
      mg_get_trait_component_variance(component_tbl, tr, "evaluator:prompt_type")   / (n_models * n_prompts) +
      mg_get_trait_component_variance(component_tbl, tr, "evaluator:temperature")   / (n_models * n_temperatures) +
      mg_get_trait_component_variance(component_tbl, tr, "evaluator:seed")          / (n_models * n_seeds) +
      mg_get_trait_component_variance(component_tbl, tr, "prompt_type:temperature") / (n_prompts * n_temperatures) +
      mg_get_trait_component_variance(component_tbl, tr, "prompt_type:seed")        / (n_prompts * n_seeds) +
      mg_get_trait_component_variance(component_tbl, tr, "temperature:seed")        / (n_temperatures * n_seeds)

    phi_denom <- sigma_tau + sigma_delta + sigma_abs_extra
    tibble(
      trait       = tr,
      predicted_g = sigma_tau / (sigma_tau + sigma_delta),
      predicted_phi = if (phi_denom == 0) NA_real_ else sigma_tau / phi_denom
    )
  }))
}


# Compute per-trait σ²_τ, σ²_δ, and σ²_abs_extra vectors from a component
# table (already extracted from a fitted lmer object).
.mg_per_trait_variances <- function(component_tbl, traits,
                                     n_models, n_prompts, n_temperatures, n_seeds) {
  sigma_tau_vec <- vapply(traits, function(tr)
    mg_get_trait_component_variance(component_tbl, tr, "item_id"), numeric(1))

  sigma_delta_vec <- vapply(traits, function(tr)
    mg_get_trait_component_variance(component_tbl, tr, "item_id:evaluator")   / n_models +
    mg_get_trait_component_variance(component_tbl, tr, "item_id:prompt_type") / n_prompts +
    mg_get_trait_component_variance(component_tbl, tr, "item_id:temperature") / n_temperatures +
    mg_get_trait_component_variance(component_tbl, tr, "item_id:seed")        / n_seeds +
    mg_get_trait_component_variance(component_tbl, tr, "Residual") /
      (n_models * n_prompts * n_temperatures * n_seeds),
    numeric(1))

  sigma_abs_extra_vec <- vapply(traits, function(tr)
    mg_get_trait_component_variance(component_tbl, tr, "evaluator")               / n_models +
    mg_get_trait_component_variance(component_tbl, tr, "prompt_type")             / n_prompts +
    mg_get_trait_component_variance(component_tbl, tr, "temperature")             / n_temperatures +
    mg_get_trait_component_variance(component_tbl, tr, "seed")                    / n_seeds +
    mg_get_trait_component_variance(component_tbl, tr, "evaluator:prompt_type")   / (n_models * n_prompts) +
    mg_get_trait_component_variance(component_tbl, tr, "evaluator:temperature")   / (n_models * n_temperatures) +
    mg_get_trait_component_variance(component_tbl, tr, "evaluator:seed")          / (n_models * n_seeds) +
    mg_get_trait_component_variance(component_tbl, tr, "prompt_type:temperature") / (n_prompts * n_temperatures) +
    mg_get_trait_component_variance(component_tbl, tr, "prompt_type:seed")        / (n_prompts * n_seeds) +
    mg_get_trait_component_variance(component_tbl, tr, "temperature:seed")        / (n_temperatures * n_seeds),
    numeric(1))

  list(tau = sigma_tau_vec, delta = sigma_delta_vec, abs_extra = sigma_abs_extra_vec)
}


#' Compute normalized composite weights for a multivariate G-theory analysis.
#'
#' Uses the diagonal approximation (cross-trait covariances = 0), consistent
#' with the || model specification. Two schemes are returned:
#'   equal   — w_k = 1/K (uniform)
#'   optimal — w_k ∝ σ²_τ[k] / σ²_δ[k] (maximizes Eρ²_c; normalized to sum = 1)
#'
#' @param component_tbl  Output of mg_multivariate_component_table().
#' @param traits         Character vector of trait names (ordered).
#' @param n_models,...   Design counts for computing σ²_δ per trait.
#' @return Named list with elements "equal" and "optimal", each a named
#'   numeric vector of weights summing to 1.
mg_composite_weights <- function(component_tbl, traits,
                                  n_models, n_prompts, n_temperatures, n_seeds) {
  K    <- length(traits)
  pvec <- .mg_per_trait_variances(component_tbl, traits,
                                   n_models, n_prompts, n_temperatures, n_seeds)

  # Equal weights
  w_equal <- setNames(rep(1 / K, K), traits)

  # Optimal weights: proportional to signal-to-noise ratio σ²_τ / σ²_δ.
  # For diagonal Σ_τ and Σ_δ this is the leading generalized eigenvector of
  # the pair (Σ_τ, Σ_δ), normalized to sum to 1 rather than unit L2 norm so
  # the composite score remains interpretable as a weighted average.
  snr       <- ifelse(pvec$delta > 0, pvec$tau / pvec$delta, 0)
  snr_safe  <- pmax(snr, 0)
  w_optimal <- if (sum(snr_safe) > 0) {
    setNames(snr_safe / sum(snr_safe), traits)
  } else {
    w_equal  # fall back to equal weights if all SNR = 0
  }

  list(equal = w_equal, optimal = w_optimal)
}


#' Composite Eρ²_c and Φ_c given weight vector and per-trait components.
#'
#' Diagonal approximation:
#'   Eρ²_c = (Σ_k w_k² σ²_τ[k]) / (Σ_k w_k² σ²_τ[k] + Σ_k w_k² σ²_δ[k])
#'   Φ_c   = (Σ_k w_k² σ²_τ[k]) / (Σ_k w_k² (σ²_τ[k] + σ²_δ[k] + σ²_abs[k]))
#'
#' @param w              Named numeric vector of composite weights (sum = 1).
#' @param component_tbl  Output of mg_multivariate_component_table().
#' @param traits         Character vector of trait names.
#' @param n_models,...   Design counts.
#' @return Named list: erho2_c, phi_c.
mg_composite_g_from_weights <- function(w, component_tbl, traits,
                                         n_models, n_prompts, n_temperatures, n_seeds) {
  pvec <- .mg_per_trait_variances(component_tbl, traits,
                                   n_models, n_prompts, n_temperatures, n_seeds)
  w2   <- w^2

  num       <- sum(w2 * pvec$tau)
  denom_rel <- sum(w2 * pvec$delta)
  denom_abs <- sum(w2 * pvec$abs_extra)

  erho2_c <- if ((num + denom_rel) == 0) NA_real_ else num / (num + denom_rel)
  phi_c   <- if ((num + denom_rel + denom_abs) == 0) NA_real_ else
               num / (num + denom_rel + denom_abs)

  list(erho2_c = erho2_c, phi_c = phi_c)
}


#' Bootstrap CI for composite Eρ²_c and Φ_c (and optionally the difference
#' Eρ²_c_optimal - Eρ²_holistic from a univariate "holistic" analysis).
#'
#' Both equal-weight and optimal-weight composites are bootstrapped in the
#' same loop.  When \code{fit_obj_uni} / \code{dat_uni} are supplied the
#' univariate holistic model is re-fitted on the SAME item resample in each
#' iteration (joint bootstrap), producing paired difference statistics.
#'
#' @param fit_obj_multi  Fitted lmer object (best multivariate model, with
#'   trait indicator columns in \code{dat_multi}).
#' @param dat_multi      Long-format multivariate data (output of
#'   mg_prepare_multivariate_data; WITHOUT trait indicator columns — they are
#'   added inside this function).
#' @param fit_obj_uni,dat_uni  Optional univariate counterparts for joint
#'   bootstrap difference CI.
#' @param n_models,...  Design counts (default = observed levels).
#' @param B,seed,conf,min_valid  Bootstrap parameters.
#' @return A tibble with one row per (weighting, metric) combination containing
#'   point estimates and percentile CIs; plus difference rows if univariate
#'   arguments are provided.
mg_bootstrap_composite_ci <- function(
    fit_obj_multi,
    dat_multi,
    fit_obj_uni    = NULL,
    dat_uni        = NULL,
    n_models       = nlevels(dat_multi$evaluator),
    n_prompts      = nlevels(dat_multi$prompt_type),
    n_temperatures = nlevels(dat_multi$temperature),
    n_seeds        = nlevels(dat_multi$seed),
    B              = 2000L,
    seed           = 42L,
    conf           = 0.95,
    min_valid      = 100L) {

  set.seed(seed)

  traits         <- levels(dat_multi$trait)
  items          <- levels(dat_multi$item_id)
  n_items        <- length(items)
  formula_multi  <- formula(fit_obj_multi)
  formula_uni    <- if (!is.null(fit_obj_uni)) formula(fit_obj_uni) else NULL
  do_diff        <- !is.null(fit_obj_uni) && !is.null(dat_uni)
  control        <- lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000))

  # Observed composite G (point estimates)
  obs_comp_tbl  <- mg_multivariate_component_table(fit_obj_multi, traits = traits)
  w_obs         <- mg_composite_weights(obs_comp_tbl, traits, n_models, n_prompts,
                                        n_temperatures, n_seeds)
  obs_eq        <- mg_composite_g_from_weights(w_obs$equal,   obs_comp_tbl, traits,
                                               n_models, n_prompts, n_temperatures, n_seeds)
  obs_opt       <- mg_composite_g_from_weights(w_obs$optimal, obs_comp_tbl, traits,
                                               n_models, n_prompts, n_temperatures, n_seeds)
  obs_uni_g     <- if (do_diff)
    mg_predict_univariate_g(fit_obj_uni, dat_uni, n_models, n_prompts, n_temperatures, n_seeds)
    else NA_real_

  # Bootstrap storage
  boot_erho2_eq   <- numeric(B);  boot_phi_eq   <- numeric(B)
  boot_erho2_opt  <- numeric(B);  boot_phi_opt  <- numeric(B)
  boot_diff_eq    <- numeric(B);  boot_diff_opt <- numeric(B)

  n_error_multi  <- 0L;  n_error_uni    <- 0L
  n_singular     <- 0L;  n_nonconverged <- 0L
  used_mask      <- rep(FALSE, B)   # TRUE when multivariate half succeeded
  diff_used_mask <- rep(FALSE, B)   # TRUE when BOTH multivariate AND univariate succeeded

  classify_convergence_warning <- function(warnings) {
    if (length(warnings) == 0) return(FALSE)
    any(grepl(
      "failed to converge|degenerate Hessian|unable to evaluate scaled gradient|false convergence|convergence code",
      warnings, ignore.case = TRUE
    ))
  }

  for (b in seq_len(B)) {
    if (b %% 200 == 0) cat(sprintf("  Composite bootstrap: %d / %d\n", b, B))

    selected <- sample(items, n_items, replace = TRUE)

    # --- Multivariate resample ---
    boot_multi <- bind_rows(lapply(seq_along(selected), function(k) {
      rows           <- dat_multi[dat_multi$item_id == selected[k], , drop = FALSE]
      rows$item_id   <- paste0("i", k)
      rows
    }))
    boot_multi$item_id <- factor(boot_multi$item_id)
    boot_multi <- mg_add_trait_indicators(boot_multi)

    fit_res_m <- mg_fit_with_warnings(suppressMessages(
      lmer(formula_multi, data = boot_multi, REML = FALSE, control = control)
    ))

    if (inherits(fit_res_m$value, "error")) {
      n_error_multi <- n_error_multi + 1L
      next
    }
    if (classify_convergence_warning(fit_res_m$warnings)) {
      n_nonconverged <- n_nonconverged + 1L
      next
    }

    comp_b    <- mg_multivariate_component_table(fit_res_m$value, traits = traits)
    tau_b     <- vapply(traits, function(tr)
      mg_get_trait_component_variance(comp_b, tr, "item_id"), numeric(1))

    if (all(tau_b < 1e-6)) {
      n_singular <- n_singular + 1L
      next
    }

    # Compute composite G for this resample
    w_b      <- mg_composite_weights(comp_b, traits, n_models, n_prompts,
                                     n_temperatures, n_seeds)
    g_eq_b   <- mg_composite_g_from_weights(w_b$equal,   comp_b, traits,
                                             n_models, n_prompts, n_temperatures, n_seeds)
    g_opt_b  <- mg_composite_g_from_weights(w_b$optimal, comp_b, traits,
                                             n_models, n_prompts, n_temperatures, n_seeds)

    boot_erho2_eq[b]  <- g_eq_b$erho2_c
    boot_phi_eq[b]    <- g_eq_b$phi_c
    boot_erho2_opt[b] <- g_opt_b$erho2_c
    boot_phi_opt[b]   <- g_opt_b$phi_c

    # --- Univariate resample (joint bootstrap for difference CI) ---
    if (do_diff) {
      boot_uni <- bind_rows(lapply(seq_along(selected), function(k) {
        rows           <- dat_uni[dat_uni$item_id == selected[k], , drop = FALSE]
        rows$item_id   <- paste0("i", k)
        rows
      }))
      boot_uni$item_id <- factor(boot_uni$item_id)

      fit_res_u <- mg_fit_with_warnings(suppressMessages(
        lmer(formula_uni, data = boot_uni, REML = FALSE, control = control)
      ))

      if (!inherits(fit_res_u$value, "error")) {
        comp_u  <- mg_univariate_component_table(fit_res_u$value)
        g_uni_b <- if (mg_get_component_variance(comp_u, "item_id") >= 1e-6)
          mg_predict_univariate_g(fit_res_u$value, boot_uni, n_models, n_prompts,
                                  n_temperatures, n_seeds)
          else NA_real_
      } else {
        g_uni_b <- NA_real_
        n_error_uni <- n_error_uni + 1L
      }

      # A paired draw is only valid when BOTH sides succeeded.
      # diff_used_mask is tracked separately from used_mask so the composite
      # CIs (multivariate-only) are not penalised by univariate failures, and
      # so the difference CI diagnostics reflect the true paired count.
      both_ok <- !is.na(boot_erho2_eq[b]) && !is.na(g_uni_b)
      boot_diff_eq[b]  <- if (both_ok) boot_erho2_eq[b]  - g_uni_b else NA_real_
      boot_diff_opt[b] <- if (both_ok) boot_erho2_opt[b] - g_uni_b else NA_real_
      diff_used_mask[b] <- both_ok
    }

    used_mask[b] <- !is.na(boot_erho2_eq[b])
  }

  n_used      <- sum(used_mask)
  n_diff_used <- if (do_diff) sum(diff_used_mask) else NA_integer_

  if (n_used < B * 0.8) {
    warning(sprintf(
      "Composite bootstrap: only %d / %d resamples usable (<80%%).",
      n_used, B
    ))
  }
  if (do_diff && n_diff_used < B * 0.8) {
    warning(sprintf(
      "Composite bootstrap (difference): only %d / %d paired resamples usable (<80%%).",
      n_diff_used, B
    ))
  }

  probs  <- c((1 - conf) / 2, 1 - (1 - conf) / 2)
  mk_ci  <- function(v, mask, n) {
    if (n < min_valid) c(NA_real_, NA_real_) else
      unname(quantile(v[mask], probs = probs, na.rm = TRUE))
  }

  ci_eq_g    <- mk_ci(boot_erho2_eq,  used_mask,      n_used)
  ci_eq_phi  <- mk_ci(boot_phi_eq,    used_mask,      n_used)
  ci_opt_g   <- mk_ci(boot_erho2_opt, used_mask,      n_used)
  ci_opt_phi <- mk_ci(boot_phi_opt,   used_mask,      n_used)
  ci_diff_eq  <- if (do_diff) mk_ci(boot_diff_eq,  diff_used_mask, n_diff_used) else c(NA_real_, NA_real_)
  ci_diff_opt <- if (do_diff) mk_ci(boot_diff_opt, diff_used_mask, n_diff_used) else c(NA_real_, NA_real_)

  # Composite-only rows use n_used; difference rows use n_diff_used
  rows <- list(
    tibble(weighting = "equal",   metric = "erho2_c",
           observed = obs_eq$erho2_c,  ci_lower = ci_eq_g[1],   ci_upper = ci_eq_g[2],
           n_pairs_used = as.integer(n_used)),
    tibble(weighting = "equal",   metric = "phi_c",
           observed = obs_eq$phi_c,    ci_lower = ci_eq_phi[1],  ci_upper = ci_eq_phi[2],
           n_pairs_used = as.integer(n_used)),
    tibble(weighting = "optimal", metric = "erho2_c",
           observed = obs_opt$erho2_c, ci_lower = ci_opt_g[1],  ci_upper = ci_opt_g[2],
           n_pairs_used = as.integer(n_used)),
    tibble(weighting = "optimal", metric = "phi_c",
           observed = obs_opt$phi_c,   ci_lower = ci_opt_phi[1], ci_upper = ci_opt_phi[2],
           n_pairs_used = as.integer(n_used))
  )

  if (do_diff) {
    rows <- c(rows, list(
      tibble(weighting = "equal",   metric = "diff_erho2_c_vs_holistic",
             observed = obs_eq$erho2_c  - obs_uni_g,
             ci_lower = ci_diff_eq[1],  ci_upper = ci_diff_eq[2],
             n_pairs_used = as.integer(n_diff_used)),
      tibble(weighting = "optimal", metric = "diff_erho2_c_vs_holistic",
             observed = obs_opt$erho2_c - obs_uni_g,
             ci_lower = ci_diff_opt[1], ci_upper = ci_diff_opt[2],
             n_pairs_used = as.integer(n_diff_used))
    ))
  }

  bind_rows(rows) %>%
    mutate(conf_level     = conf,
           n_error_multi  = as.integer(n_error_multi),
           n_error_uni    = as.integer(n_error_uni),
           n_singular     = as.integer(n_singular),
           n_nonconverged = as.integer(n_nonconverged),
           seed           = as.integer(seed))
}


mg_run_multivariate_manuscript <- function(
    task,
    data_path     = NULL,
    data          = NULL,
    trait_subset  = NULL,
    output_prefix = NULL,
    model_groups  = mg_multivariate_model_groups(),
    target_g      = c(0.80, 0.85, 0.90, 0.95),
    model_range   = NULL,
    prompt_range  = NULL,
    ref_fit_uni   = NULL,    # univariate best-fit lmer for joint bootstrap difference CI
    ref_dat_uni   = NULL,    # univariate prepared data for joint bootstrap
    n_boot        = 2000L,
    boot_seed     = 42L) {
  dat <- mg_prepare_multivariate_data(
    task = task,
    data_path = data_path,
    data = data,
    trait_subset = trait_subset
  )

  fit_bundle <- mg_fit_multivariate_models(dat = dat, model_groups = model_groups)
  best_row <- mg_select_best_model(fit_bundle$summary)
  best_fit <- fit_bundle$fits[[best_row$model]]

  component_tbl <- mg_multivariate_component_table(best_fit, traits = levels(dat$trait))
  facet_tbl <- mg_multivariate_facet_table(component_tbl)
  mean_facet_tbl <- mg_multivariate_mean_facet_table(facet_tbl)
  prediction <- mg_multivariate_prediction_grid(
    fit_obj = best_fit,
    dat = dat,
    model_range = if (is.null(model_range)) seq_len(max(MG_DEFAULT_MAX_MODELS, nlevels(dat$evaluator))) else model_range,
    prompt_range = if (is.null(prompt_range)) seq_len(max(MG_DEFAULT_MAX_PROMPTS, nlevels(dat$prompt_type))) else prompt_range,
    target_g = target_g
  )

  prediction_points <- mg_multivariate_prediction_points(
    fit_obj = best_fit,
    dat = dat,
    target_tbl = prediction$targets,
    target_g = target_g
  )

  observed_trait_g <- mg_predict_multivariate_g(
    fit_obj = best_fit,
    dat = dat,
    n_models = nlevels(dat$evaluator),
    n_prompts = nlevels(dat$prompt_type)
  )

  observed_trait_phi <- mg_predict_multivariate_phi(
    fit_obj = best_fit,
    dat = dat,
    n_models = nlevels(dat$evaluator),
    n_prompts = nlevels(dat$prompt_type)
  )

  # Composite G and Φ (equal and optimal weights)
  n_m <- nlevels(dat$evaluator); n_p <- nlevels(dat$prompt_type)
  n_t <- nlevels(dat$temperature); n_s <- nlevels(dat$seed)

  comp_w   <- mg_composite_weights(component_tbl, levels(dat$trait), n_m, n_p, n_t, n_s)
  comp_eq  <- mg_composite_g_from_weights(comp_w$equal,   component_tbl, levels(dat$trait), n_m, n_p, n_t, n_s)
  comp_opt <- mg_composite_g_from_weights(comp_w$optimal, component_tbl, levels(dat$trait), n_m, n_p, n_t, n_s)

  composite_g_tbl <- tibble(
    weighting  = c("equal",         "equal",    "optimal",        "optimal"),
    metric     = c("erho2_c",       "phi_c",    "erho2_c",        "phi_c"),
    observed   = c(comp_eq$erho2_c, comp_eq$phi_c, comp_opt$erho2_c, comp_opt$phi_c)
  )

  cat(sprintf(
    "  Composite Eρ²_c: equal = %.3f, optimal = %.3f\n",
    comp_eq$erho2_c, comp_opt$erho2_c
  ))

  # Composite bootstrap CI (and joint difference CI if univariate ref is provided)
  composite_boot_ci <- NULL
  if (n_boot > 0L) {
    cat(sprintf("Running %d composite bootstrap resamples for %s...\n", n_boot, task))
    composite_boot_ci <- mg_bootstrap_composite_ci(
      fit_obj_multi  = best_fit,
      dat_multi      = dat,
      fit_obj_uni    = ref_fit_uni,
      dat_uni        = ref_dat_uni,
      n_models       = n_m,
      n_prompts      = n_p,
      n_temperatures = n_t,
      n_seeds        = n_s,
      B              = n_boot,
      seed           = boot_seed
    )
  }

  primary_source <- mean_facet_tbl %>%
    filter(!facet_group %in% c("item")) %>%
    arrange(desc(mean_variance)) %>%
    slice(1)

  primary_source_name <- if (nrow(primary_source) == 0) NA_character_ else primary_source$facet_group[[1]]
  primary_source_prop <- if (nrow(primary_source) == 0) NA_real_ else primary_source$mean_prop_total[[1]]

  # Identify boundary components (zero variance) across any trait
  boundary_comps <- component_tbl %>%
    filter(variance == 0) %>%
    distinct(component) %>%
    pull(component)
  boundary_str <- if (length(boundary_comps) == 0) NA_character_ else paste(boundary_comps, collapse = "; ")

  overview_tbl <- tibble(
    task = task,
    analysis = "multivariate",
    best_model = best_row$model,
    best_model_BIC = best_row$BIC,
    best_model_singular = best_row$singular,
    boundary_components = boundary_str,
    observed_items = nlevels(dat$item_id),
    observed_models = nlevels(dat$evaluator),
    observed_prompts = nlevels(dat$prompt_type),
    observed_temperatures = nlevels(dat$temperature),
    observed_seeds = nlevels(dat$seed),
    n_traits = nlevels(dat$trait),
    observed_mean_g   = mean(observed_trait_g$predicted_g,   na.rm = TRUE),
    observed_min_g    = min(observed_trait_g$predicted_g,    na.rm = TRUE),
    observed_mean_phi = mean(observed_trait_phi$predicted_phi, na.rm = TRUE),
    composite_erho2_equal   = comp_eq$erho2_c,
    composite_phi_equal     = comp_eq$phi_c,
    composite_erho2_optimal = comp_opt$erho2_c,
    composite_phi_optimal   = comp_opt$phi_c,
    primary_source = primary_source_name,
    primary_source_mean_prop = primary_source_prop
  )

  mg_write_csv(fit_bundle$summary,     output_prefix, "model_bic")
  mg_write_csv(component_tbl,          output_prefix, "best_model_components")
  mg_write_csv(facet_tbl,              output_prefix, "best_model_facets_by_trait")
  mg_write_csv(mean_facet_tbl,         output_prefix, "best_model_facets_mean")
  mg_write_csv(observed_trait_g,       output_prefix, "observed_design_g_by_trait")
  mg_write_csv(observed_trait_phi,     output_prefix, "observed_design_phi_by_trait")
  mg_write_csv(composite_g_tbl,        output_prefix, "composite_g")
  mg_write_csv(composite_boot_ci,      output_prefix, "composite_bootstrap_ci")
  mg_write_csv(prediction$grid,        output_prefix, "prediction_grid")
  mg_write_csv(prediction$targets,     output_prefix, "prediction_targets")
  mg_write_csv(prediction_points,      output_prefix, "prediction_points")
  mg_write_csv(overview_tbl,           output_prefix, "overview")

  list(
    data = dat,
    fits = fit_bundle$fits,
    model_table = fit_bundle$summary,
    best_model = best_row,
    best_fit = best_fit,
    component_table = component_tbl,
    facet_table = facet_tbl,
    mean_facet_table = mean_facet_tbl,
    observed_trait_g = observed_trait_g,
    observed_trait_phi = observed_trait_phi,
    composite_g = composite_g_tbl,
    composite_bootstrap_ci = composite_boot_ci,
    prediction_grid = prediction$grid,
    prediction_targets = prediction$targets,
    prediction_points = prediction_points,
    overview = overview_tbl
  )
}
