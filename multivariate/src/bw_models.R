# ---------------------------------------------------------------------------
# bw_models.R (MULTIVARIATE)
#
# OpenMx Between-Within (RAM/path) specification for multivariate nested
# G-theory models. Each variance component is a between-level submodel
# with p latent variables and a p×p covariance matrix, connected to
# p manifest variables at the observation level via primaryKey/joinKey.
#
# Two models (nested facets only — full crossed uses wide_models.R):
#   1. nested_seed             — seed nested in e×p×t
#   2. nested_seed_temperature — seed & temp nested in e×p
#
# Expected working directory: 01. Research/
# ---------------------------------------------------------------------------

suppressPackageStartupMessages(library(OpenMx))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(lme4))

source("multivariate/src/shared/multivariate_workflow.R")

mxOption(model = NULL, key = "Default optimizer", "CSOLNP", reset = FALSE)


mvbw_log <- function(...) {
  cat(...)
  flush.console()
}

mvbw_timed <- function(label, expr) {
  t0 <- proc.time()[["elapsed"]]
  value <- force(expr)
  elapsed <- proc.time()[["elapsed"]] - t0
  mvbw_log(sprintf("%s completed in %.2fs\n", label, elapsed))
  value
}


# ---------------------------------------------------------------------------
# Data preparation: wide format with p trait columns + composite keys
# ---------------------------------------------------------------------------

mvbw_prepare_data <- function(task, data_path, trait_subset = NULL) {
  raw <- read.csv(data_path, stringsAsFactors = FALSE) %>% as_tibble()

  if (task == "mental_health") {
    trait_map <- MG_MH_TRAITS
    if (!is.null(trait_subset)) trait_map <- trait_map[trait_map %in% trait_subset]
    for (col in names(trait_map)) {
      raw[[trait_map[[col]]]] <- as.numeric(mg_coerce_binary_flag(raw[[col]]))
    }
    traits <- unname(trait_map)
  } else if (task == "drug_review") {
    trait_map <- MG_DRUG_TRAITS
    if (!is.null(trait_subset)) trait_map <- trait_map[trait_map %in% trait_subset]
    for (col in names(trait_map)) {
      raw[[trait_map[[col]]]] <- as.numeric(MG_DRUG_SCORE_MAP[as.character(raw[[col]])])
    }
    traits <- unname(trait_map)
  } else {
    stop("Unsupported multivariate task: ", task)
  }

  dat <- raw %>%
    mutate(
      item_id     = factor(item_id),
      evaluator   = factor(evaluator),
      prompt_type = factor(prompt_type),
      temperature = factor(temperature),
      seed        = factor(seed),
      item_eval   = factor(paste(item_id, evaluator, sep = ":")),
      item_prompt = factor(paste(item_id, prompt_type, sep = ":")),
      item_temp   = factor(paste(item_id, temperature, sep = ":")),
      item_seed   = factor(paste(item_id, seed, sep = ":")),
      eval_prompt = factor(paste(evaluator, prompt_type, sep = ":")),
      eval_temp   = factor(paste(evaluator, temperature, sep = ":")),
      eval_seed   = factor(paste(evaluator, seed, sep = ":")),
      prompt_temp = factor(paste(prompt_type, temperature, sep = ":")),
      prompt_seed = factor(paste(prompt_type, seed, sep = ":")),
      temp_seed   = factor(paste(temperature, seed, sep = ":")),
      ept_cell    = factor(paste(evaluator, prompt_type, temperature, sep = ":")),
      epts_cell   = factor(paste(evaluator, prompt_type, temperature, seed, sep = ":"))
    )

  list(data = dat, traits = traits, p = length(traits))
}


mvbw_subset_items <- function(dat, max_items = NULL) {
  if (is.null(max_items)) return(dat)
  all_items <- levels(dat$item_id)
  if (length(all_items) <= max_items) return(dat)
  keep_items <- all_items[seq_len(max_items)]
  dat %>% filter(item_id %in% keep_items) %>% droplevels()
}


# ---------------------------------------------------------------------------
# Nested component definitions (same keys as univariate)
# ---------------------------------------------------------------------------

mvbw_nested_seed_components <- function() {
  list(
    list(name = "item",        key = "item_id"),
    list(name = "eval",        key = "evaluator"),
    list(name = "prompt",      key = "prompt_type"),
    list(name = "temp",        key = "temperature"),
    list(name = "item_eval",   key = "item_eval"),
    list(name = "item_prompt", key = "item_prompt"),
    list(name = "item_temp",   key = "item_temp"),
    list(name = "eval_prompt", key = "eval_prompt"),
    list(name = "eval_temp",   key = "eval_temp"),
    list(name = "prompt_temp", key = "prompt_temp"),
    list(name = "seed_in_ept", key = "epts_cell")
  )
}

mvbw_nested_seed_temp_components <- function() {
  list(
    list(name = "item",        key = "item_id"),
    list(name = "eval",        key = "evaluator"),
    list(name = "prompt",      key = "prompt_type"),
    list(name = "item_eval",   key = "item_eval"),
    list(name = "item_prompt", key = "item_prompt"),
    list(name = "eval_prompt", key = "eval_prompt"),
    list(name = "temp_in_ep",  key = "ept_cell"),
    list(name = "seed_in_ept", key = "epts_cell")
  )
}


# ---------------------------------------------------------------------------
# Start values: p independent univariate lme4 models → diagonal p×p matrices
# ---------------------------------------------------------------------------

mvbw_lme4_start_values <- function(dat, traits) {
  formula_str <- paste(
    "y ~ 1 +",
    "(1|item_id) + (1|evaluator) + (1|prompt_type) +",
    "(1|temperature) + (1|seed) +",
    "(1|item_id:evaluator) + (1|item_id:prompt_type) +",
    "(1|item_id:temperature) + (1|item_id:seed) +",
    "(1|evaluator:prompt_type) + (1|evaluator:temperature) +",
    "(1|evaluator:seed) + (1|prompt_type:temperature) +",
    "(1|prompt_type:seed) + (1|temperature:seed)"
  )

  grp_names <- c("item_id", "evaluator", "prompt_type", "temperature", "seed",
                  "item_id:evaluator", "item_id:prompt_type",
                  "item_id:temperature", "item_id:seed",
                  "evaluator:prompt_type", "evaluator:temperature",
                  "evaluator:seed", "prompt_type:temperature",
                  "prompt_type:seed", "temperature:seed", "Residual")

  p <- length(traits)
  sv <- setNames(
    lapply(grp_names, function(g) setNames(rep(NA_real_, p), traits)),
    grp_names
  )

  for (tr in traits) {
    mvbw_log(sprintf("  lme4 start values [%s] ... ", tr))
    tmp_dat <- dat
    tmp_dat$y <- as.numeric(tmp_dat[[tr]])

    fit <- tryCatch(
      lmer(as.formula(formula_str), data = tmp_dat, REML = FALSE,
           control = lmerControl(optimizer = "bobyqa",
                                 optCtrl = list(maxfun = 100000))),
      error = function(e) NULL
    )

    if (is.null(fit)) {
      mvbw_log("failed, using generic\n")
      generic <- var(tmp_dat$y, na.rm = TRUE) / 16
      for (g in grp_names) sv[[g]][[tr]] <- generic
    } else {
      mvbw_log("done\n")
      vc <- as.data.frame(VarCorr(fit))
      vc_map <- setNames(vc$vcov, vc$grp)
      for (g in grp_names) sv[[g]][[tr]] <- vc_map[[g]] %||% 0
    }
  }

  sv
}


mvbw_map_starts <- function(lme4_sv, model_type, traits) {
  p <- length(traits)
  resid_vars <- pmax(lme4_sv[["Residual"]], 1e-6)
  eps <- resid_vars * 1e-4

  diag_mat <- function(grp_name) {
    v <- pmax(lme4_sv[[grp_name]] %||% rep(0, p), eps)
    m <- matrix(0, p, p, dimnames = list(traits, traits))
    diag(m) <- v
    m
  }

  sum_diag_mats <- function(grp_names) {
    Reduce(`+`, lapply(grp_names, diag_mat))
  }

  if (model_type == "nested_seed") {
    list(
      item        = diag_mat("item_id"),
      eval        = diag_mat("evaluator"),
      prompt      = diag_mat("prompt_type"),
      temp        = diag_mat("temperature"),
      item_eval   = diag_mat("item_id:evaluator"),
      item_prompt = diag_mat("item_id:prompt_type"),
      item_temp   = diag_mat("item_id:temperature"),
      eval_prompt = diag_mat("evaluator:prompt_type"),
      eval_temp   = diag_mat("evaluator:temperature"),
      prompt_temp = diag_mat("prompt_type:temperature"),
      seed_in_ept = sum_diag_mats(c("seed", "evaluator:seed",
                                     "prompt_type:seed", "temperature:seed")),
      residual    = diag_mat("Residual")
    )
  } else if (model_type == "nested_seed_temperature") {
    list(
      item        = diag_mat("item_id"),
      eval        = diag_mat("evaluator"),
      prompt      = diag_mat("prompt_type"),
      item_eval   = diag_mat("item_id:evaluator"),
      item_prompt = diag_mat("item_id:prompt_type"),
      eval_prompt = diag_mat("evaluator:prompt_type"),
      temp_in_ep  = sum_diag_mats(c("temperature", "evaluator:temperature",
                                     "prompt_type:temperature")),
      seed_in_ept = sum_diag_mats(c("seed", "evaluator:seed",
                                     "prompt_type:seed", "temperature:seed")),
      residual    = diag_mat("Residual")
    )
  } else {
    stop("Unknown model_type: ", model_type)
  }
}


# ---------------------------------------------------------------------------
# Build a single between-level submodel (p latents, p×p covariance)
# ---------------------------------------------------------------------------

mvbw_make_submodel <- function(comp_name, key_col, dat, traits, start_cov) {
  p <- length(traits)
  key_vals <- levels(dat[[key_col]])
  key_data <- data.frame(key = factor(key_vals, levels = key_vals))
  names(key_data) <- key_col

  lat_names <- paste0("lat_", comp_name, "_", traits)

  # Variance paths
  var_paths <- lapply(seq_len(p), function(i) {
    mxPath(from = lat_names[i], arrows = 2,
           free = TRUE, values = start_cov[i, i], lbound = 0,
           labels = paste0("s_", comp_name, "_", traits[i], "_", traits[i]))
  })

  # Covariance paths
  cov_paths <- list()
  if (p > 1) {
    idx <- 1
    for (i in 1:(p - 1)) {
      for (j in (i + 1):p) {
        cov_paths[[idx]] <- mxPath(
          from = lat_names[i], to = lat_names[j], arrows = 2,
          free = TRUE, values = start_cov[i, j],
          labels = paste0("s_", comp_name, "_", traits[i], "_", traits[j])
        )
        idx <- idx + 1
      }
    }
  }

  mxModel(
    name = comp_name, type = "RAM",
    latentVars = lat_names,
    mxData(key_data, type = "raw", primaryKey = key_col),
    var_paths, cov_paths
  )
}


# ---------------------------------------------------------------------------
# Build the observation-level model with all submodels
# ---------------------------------------------------------------------------

mvbw_build_model <- function(dat, components, starts, traits,
                             model_name = "mgtheory_bw") {
  p <- length(traits)
  key_cols <- vapply(components, function(comp) comp$key, character(1))

  obs_data <- dat %>%
    mutate(across(all_of(traits), as.numeric)) %>%
    mutate(across(all_of(key_cols), factor)) %>%
    as.data.frame()

  # Between-level submodels
  submodels <- lapply(components, function(comp) {
    sv <- starts[[comp$name]]
    if (is.null(sv)) sv <- diag(0.1, p)
    mvbw_make_submodel(comp$name, comp$key, obs_data, traits, start_cov = sv)
  })

  # joinKey paths: each submodel latent → corresponding manifest trait
  join_paths <- unlist(lapply(components, function(comp) {
    lapply(seq_len(p), function(k) {
      lat_name <- paste0(comp$name, ".lat_", comp$name, "_", traits[k])
      mxPath(from = lat_name, to = traits[k],
             free = FALSE, values = 1, joinKey = comp$key)
    })
  }), recursive = FALSE)

  # Means for each trait
  mean_paths <- lapply(seq_len(p), function(k) {
    mxPath(from = "one", to = traits[k], free = TRUE,
           values = mean(obs_data[[traits[k]]], na.rm = TRUE))
  })

  # Residual p×p covariance at observation level
  resid_start <- starts[["residual"]]
  if (is.null(resid_start)) resid_start <- diag(0.1, p)

  resid_var_paths <- lapply(seq_len(p), function(i) {
    mxPath(from = traits[i], arrows = 2,
           free = TRUE, values = resid_start[i, i], lbound = 1e-6,
           labels = paste0("s_resid_", traits[i], "_", traits[i]))
  })

  resid_cov_paths <- list()
  if (p > 1) {
    idx <- 1
    for (i in 1:(p - 1)) {
      for (j in (i + 1):p) {
        resid_cov_paths[[idx]] <- mxPath(
          from = traits[i], to = traits[j], arrows = 2,
          free = TRUE, values = resid_start[i, j],
          labels = paste0("s_resid_", traits[i], "_", traits[j])
        )
        idx <- idx + 1
      }
    }
  }

  mxModel(
    name = model_name, type = "RAM",
    manifestVars = traits,
    submodels,
    mxData(obs_data, type = "raw"),
    mean_paths,
    resid_var_paths,
    resid_cov_paths,
    join_paths
  )
}


# ---------------------------------------------------------------------------
# Extract p×p covariance matrices from fitted model
# ---------------------------------------------------------------------------

mvbw_extract_components <- function(fit, components, traits) {
  params <- omxGetParameters(fit)
  p <- length(traits)

  rows <- list()
  idx <- 1

  extract_one <- function(comp_name) {
    for (i in seq_len(p)) {
      for (j in i:p) {
        label <- paste0("s_", comp_name, "_", traits[i], "_", traits[j])
        rows[[idx]] <<- tibble(
          component = comp_name, trait_i = traits[i], trait_j = traits[j],
          type = if (i == j) "variance" else "covariance",
          value = as.numeric(unname(params[label])),
          label = label
        )
        idx <<- idx + 1
      }
    }
  }

  for (comp in components) extract_one(comp$name)
  extract_one("resid")
  bind_rows(rows)
}


mvbw_variance_summary <- function(comp_tbl, traits) {
  comp_tbl %>%
    filter(type == "variance") %>%
    group_by(trait_i) %>%
    mutate(total_var = sum(value), prop_total = value / total_var) %>%
    ungroup() %>%
    select(trait = trait_i, component, variance = value, prop_total, total_var) %>%
    arrange(trait, desc(variance))
}


# ---------------------------------------------------------------------------
# Fit a single multivariate nested model
# ---------------------------------------------------------------------------

mvbw_fit_model <- function(dat_bundle, model_type = "nested_seed",
                           extra_tries = 5L, max_items = NULL,
                           lme4_sv = NULL) {
  dat    <- dat_bundle$data
  traits <- dat_bundle$traits
  p      <- dat_bundle$p

  dat <- mvbw_subset_items(dat, max_items = max_items)

  components <- switch(model_type,
    nested_seed             = mvbw_nested_seed_components(),
    nested_seed_temperature = mvbw_nested_seed_temp_components(),
    stop("Unknown model_type: ", model_type)
  )

  if (is.null(lme4_sv)) {
    lme4_sv <- mvbw_lme4_start_values(dat, traits)
  }
  starts <- mvbw_map_starts(lme4_sv, model_type, traits)

  n_cov_params <- (length(components) + 1L) * p * (p + 1) / 2
  model_name <- paste0("mgtheory_bw_", model_type)
  mx_model <- mvbw_build_model(dat, components, starts, traits, model_name)

  mvbw_log(sprintf(
    "Fitting MV BW [%s] — %d components, %d traits, %.0f cov params, %d items, %d rows\n",
    model_type, length(components) + 1L, p, n_cov_params,
    nlevels(dat$item_id), nrow(dat)
  ))

  if (!is.null(extra_tries)) {
    fit <- mvbw_timed(
      sprintf("  MV BW optimization [%s]", model_type),
      mxTryHard(mx_model, extraTries = extra_tries,
                OKstatuscodes = c(0L, 6L),
                jitterDistrib = "runif", loc = 1, scale = 0.25)
    )
  } else {
    fit <- mvbw_timed(
      sprintf("  MV BW optimization [%s]", model_type),
      mxRun(mx_model)
    )
  }

  status <- fit$output$status$code
  minus2LL <- fit$output$Minus2LogLikelihood
  mvbw_log(sprintf("  Status: %d | -2LL: %.4f\n", status, minus2LL))

  comp_tbl <- mvbw_extract_components(fit, components, traits)
  var_summary <- mvbw_variance_summary(comp_tbl, traits)

  list(
    model_type  = model_type,
    fit         = fit,
    status      = status,
    minus2LL    = minus2LL,
    AIC         = AIC(fit),
    BIC         = BIC(fit),
    components  = comp_tbl,
    var_summary = var_summary,
    traits      = traits
  )
}


# ---------------------------------------------------------------------------
# Run both nested models
# ---------------------------------------------------------------------------

mvbw_run_all_models <- function(dat_bundle, extra_tries = 5L,
                                max_items = NULL) {
  dat    <- dat_bundle$data
  traits <- dat_bundle$traits

  # Compute lme4 start values on full dataset
  lme4_sv <- mvbw_timed(
    "  lme4 start values (all traits)",
    mvbw_lme4_start_values(dat, traits)
  )

  results <- list()
  for (mt in c("nested_seed", "nested_seed_temperature")) {
    mvbw_log("\n", strrep("=", 60), "\n", sep = "")
    res <- mvbw_fit_model(dat_bundle, model_type = mt,
                          extra_tries = extra_tries,
                          max_items = max_items,
                          lme4_sv = lme4_sv)
    results[[mt]] <- res

    mvbw_log("\nPer-trait variance decomposition:\n")
    print(res$var_summary, n = 40, width = Inf)
  }

  # Summary
  mvbw_log("\n\n", strrep("=", 60), "\n", sep = "")
  mvbw_log("MODEL COMPARISON\n")
  mvbw_log(strrep("=", 60), "\n\n", sep = "")
  summary_tbl <- bind_rows(lapply(names(results), function(mt) {
    r <- results[[mt]]
    tibble(model = mt, status = r$status, minus2LL = r$minus2LL,
           AIC = r$AIC, BIC = r$BIC, n_traits = length(r$traits))
  }))
  print(summary_tbl, width = Inf)

  results
}
