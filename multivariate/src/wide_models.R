# ---------------------------------------------------------------------------
# wide_models.R (MULTIVARIATE)
#
# Multivariate G-theory via OpenMx path spec — full crossed model only.
# Extension of univariate wide_models.R to p traits.
# Each variance component has a p×p covariance matrix Σ_k and a C×C kernel K_k.
#
# Wide/path engine cannot model nested facets separately.
# For multivariate nested models, use bw_models.R.
#
# Expected working directory: 01. Research/
# ---------------------------------------------------------------------------

suppressPackageStartupMessages(library(OpenMx))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(tidyr))

source("multivariate/src/shared/multivariate_workflow.R")
source("univariate/src/fit_gtheory_openmx_path.R")

mxOption(model = NULL, key = "Default optimizer", "CSOLNP", reset = FALSE)


mvow_log <- function(...) {
  cat(...)
  flush.console()
}

mvow_timed <- function(label, expr) {
  t0 <- proc.time()[["elapsed"]]
  value <- force(expr)
  elapsed <- proc.time()[["elapsed"]] - t0
  mvow_log(sprintf("%s completed in %.2fs\n", label, elapsed))
  value
}


mvow_manifest_names <- function(traits, base_conds) {
  unlist(lapply(traits, function(tr) paste0(tr, "_", base_conds)), use.names = FALSE)
}


# ---------------------------------------------------------------------------
# Data preparation: wide item × (trait × condition)
# ---------------------------------------------------------------------------

mvow_prepare_data <- function(task, data_path, trait_subset = NULL) {
  raw <- read.csv(data_path, stringsAsFactors = FALSE)

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

  p <- length(traits)
  crossed_facets <- c("evaluator", "prompt_type", "temperature", "seed")

  # Get condition structure from first trait
  raw[["_score_"]] <- raw[[traits[1]]]
  invisible(capture.output(
    base_prep <- prepare_openmx_data(
      dat = raw, item_col = "item_id", score_col = "_score_",
      crossed_facets = crossed_facets, rep_facets = character(0),
      n_rep_groups = 1
    )
  ))
  base_conds <- base_prep$cond_cols
  C <- length(base_conds)

  # Build multi-trait wide data
  combined <- base_prep$wide_data[, "item", drop = FALSE]
  all_cond_cols <- character(0)

  for (tr in traits) {
    raw[["_score_"]] <- raw[[tr]]
    invisible(capture.output(
      tr_prep <- prepare_openmx_data(
        dat = raw, item_col = "item_id", score_col = "_score_",
        crossed_facets = crossed_facets, rep_facets = character(0),
        n_rep_groups = 1
      )
    ))
    tr_wide <- tr_prep$wide_data[, base_conds, drop = FALSE]
    tr_names <- paste0(tr, "_", base_conds)
    names(tr_wide) <- tr_names
    all_cond_cols <- c(all_cond_cols, tr_names)
    combined <- cbind(combined, tr_wide)
  }

  mvow_log(sprintf("=== Multivariate Data Prepared ===\n"))
  mvow_log(sprintf("Task: %s | Traits: %d (%s)\n", task, p, paste(traits, collapse = ", ")))
  mvow_log(sprintf("Items: %d | Conditions: %d | Manifest vars: %d\n",
                    nrow(combined), C, length(all_cond_cols)))

  list(
    wide_data      = combined,
    all_cond_cols  = all_cond_cols,
    base_cond_cols = base_conds,
    traits         = traits,
    p              = p,
    C              = C,
    base_prep      = base_prep,
    n_p            = base_prep$n_p,
    n_R            = base_prep$n_R
  )
}


# ---------------------------------------------------------------------------
# Build multivariate RAM group model (one variance component → p×p cov)
# ---------------------------------------------------------------------------

mvow_build_ram_group <- function(group_name, raw_data, traits, base_cond_cols,
                                 kernels, comp_maps, sv,
                                 include_facet_part = FALSE,
                                 n_items = 1L, n_R = 1L, tol = 1e-8) {
  p <- length(traits)
  C <- length(base_cond_cols)
  manifest_vars <- mvow_manifest_names(traits, base_cond_cols)
  latent_vars <- character(0)
  path_parts <- list()

  add_mv_component <- function(component_name, kernel_name, scale = 1) {
    basis <- path_kernel_basis(scale * kernels[[kernel_name]], tol = tol)
    if (is.null(basis)) return(invisible(NULL))

    n_basis <- ncol(basis)
    sv_mat <- sv[[component_name]]
    if (is.null(sv_mat)) sv_mat <- diag(0.01, p)

    for (j in seq_len(n_basis)) {
      lv_names <- paste0("lv_", component_name, "_", j, "_", traits)
      latent_vars <<- c(latent_vars, lv_names)

      # Zero means
      path_parts[[length(path_parts) + 1]] <<- mxPath(
        from = "one", to = lv_names, arrows = 1, free = FALSE, values = rep(0, p)
      )

      # Variance paths
      for (ti in seq_len(p)) {
        path_parts[[length(path_parts) + 1]] <<- mxPath(
          from = lv_names[ti], arrows = 2,
          free = TRUE, values = sv_mat[ti, ti], lbound = 0,
          labels = paste0("s_", component_name, "_", traits[ti], "_", traits[ti])
        )
      }

      # Covariance paths
      if (p > 1) {
        for (ti in 1:(p - 1)) {
          for (tj in (ti + 1):p) {
            path_parts[[length(path_parts) + 1]] <<- mxPath(
              from = lv_names[ti], to = lv_names[tj], arrows = 2,
              free = TRUE, values = sv_mat[ti, tj],
              labels = paste0("s_", component_name, "_", traits[ti], "_", traits[tj])
            )
          }
        }
      }

      # Loading paths
      for (ti in seq_len(p)) {
        trait_manifests <- paste0(traits[ti], "_", base_cond_cols)
        path_parts[[length(path_parts) + 1]] <<- mxPath(
          from = lv_names[ti], to = trait_manifests,
          arrows = 1, free = FALSE, values = basis[, j]
        )
      }
    }
  }

  # Item-side components
  for (nm in names(comp_maps$item)) {
    add_mv_component(nm, unname(comp_maps$item[[nm]]))
  }

  # Facet-side components (grandmean only)
  if (include_facet_part) {
    for (nm in names(comp_maps$facet)) {
      add_mv_component(nm, unname(comp_maps$facet[[nm]]), scale = n_items)
    }
  }

  # Residual: Σ_epsilon ⊗ I_C
  sv_eps <- sv[["epsilon"]]
  if (is.null(sv_eps)) sv_eps <- diag(0.01, p)

  for (c_idx in seq_len(C)) {
    cond <- base_cond_cols[c_idx]
    cond_manifests <- paste0(traits, "_", cond)

    for (ti in seq_len(p)) {
      path_parts[[length(path_parts) + 1]] <- mxPath(
        from = cond_manifests[ti], arrows = 2,
        free = TRUE, values = sv_eps[ti, ti] / n_R, lbound = 1e-6,
        labels = paste0("s_epsilon_", traits[ti], "_", traits[ti])
      )
    }
    if (p > 1) {
      for (ti in 1:(p - 1)) {
        for (tj in (ti + 1):p) {
          path_parts[[length(path_parts) + 1]] <- mxPath(
            from = cond_manifests[ti], to = cond_manifests[tj], arrows = 2,
            free = TRUE, values = sv_eps[ti, tj] / n_R,
            labels = paste0("s_epsilon_", traits[ti], "_", traits[tj])
          )
        }
      }
    }
  }

  # Manifest means
  if (include_facet_part) {
    for (ti in seq_len(p)) {
      trait_manifests <- paste0(traits[ti], "_", base_cond_cols)
      trait_mean <- mean(as.numeric(unlist(raw_data[, trait_manifests])), na.rm = TRUE)
      path_parts[[length(path_parts) + 1]] <- mxPath(
        from = "one", to = trait_manifests, arrows = 1, free = TRUE,
        values = rep(trait_mean, C),
        labels = rep(paste0("mean_", traits[ti]), C)
      )
    }
  } else {
    path_parts[[length(path_parts) + 1]] <- mxPath(
      from = "one", to = manifest_vars, arrows = 1, free = FALSE,
      values = rep(0, length(manifest_vars))
    )
  }

  mxModel(
    group_name, type = "RAM",
    manifestVars = manifest_vars, latentVars = latent_vars,
    mxData(observed = as.data.frame(raw_data), type = "raw"),
    path_parts
  )
}


# ---------------------------------------------------------------------------
# Moment-based start values: per-trait diagonal
# ---------------------------------------------------------------------------

mvow_moment_start_values <- function(dat_bundle) {
  traits <- dat_bundle$traits
  p <- dat_bundle$p
  wide <- dat_bundle$wide_data
  base_conds <- dat_bundle$base_cond_cols
  prep <- dat_bundle$base_prep
  crossed <- prep$crossed_facets

  comp_maps <- build_component_maps(
    crossed_facets = crossed,
    item_interaction_orders = 1L,
    facet_interaction_orders = 2L,
    n_rep_groups = 1
  )
  all_comp_names <- c(names(comp_maps$item), names(comp_maps$facet), "epsilon")

  # Per-trait: tau from between-item variance, rest split equally
  trait_starts <- list()
  for (tr in traits) {
    tr_cols <- paste0(tr, "_", base_conds)
    tr_mat <- as.matrix(wide[, tr_cols, drop = FALSE])
    grand_mean <- mean(tr_mat)
    total_var <- mean((tr_mat - grand_mean)^2)
    tau <- var(rowMeans(tr_mat))
    other_var <- max(total_var - tau, 1e-6)
    per_other <- other_var / max(length(all_comp_names) - 1, 1)

    sv <- list()
    for (nm in all_comp_names) {
      sv[[nm]] <- if (nm == "tau") tau else per_other
    }
    trait_starts[[tr]] <- sv
  }

  # Build p×p start matrices (diagonal)
  start_mats <- list()
  for (nm in all_comp_names) {
    m <- matrix(0, p, p, dimnames = list(traits, traits))
    for (ti in seq_len(p)) {
      m[ti, ti] <- max(trait_starts[[traits[ti]]][[nm]], 1e-6)
    }
    start_mats[[nm]] <- m
  }
  start_mats
}


# ---------------------------------------------------------------------------
# Fit full multivariate model
# ---------------------------------------------------------------------------

mvow_fit_full <- function(dat_bundle, extra_tries = 5L) {
  traits <- dat_bundle$traits
  p <- dat_bundle$p
  C <- dat_bundle$C
  base_conds <- dat_bundle$base_cond_cols
  prep <- dat_bundle$base_prep
  n_items <- dat_bundle$n_p
  n_R <- dat_bundle$n_R

  crossed <- prep$crossed_facets
  facet_info <- prep$facet_info
  facet_sizes <- sapply(crossed, function(f) facet_info[[f]]$n)
  short_map <- facet_short_map(crossed)
  names(facet_sizes) <- unname(short_map[crossed])

  comp_maps <- build_component_maps(
    crossed_facets = crossed,
    item_interaction_orders = 1L,
    facet_interaction_orders = 2L,
    n_rep_groups = 1
  )
  kernels <- build_structure_matrices(facet_sizes, n_g = 1)
  sv <- mvow_moment_start_values(dat_bundle)

  # Contrast/grandmean split
  dat_mat <- as.matrix(dat_bundle$wide_data[, dat_bundle$all_cond_cols, drop = FALSE])
  mean_basis <- rep(1 / sqrt(n_items), n_items)
  full_basis <- qr.Q(qr(cbind(mean_basis, diag(n_items))), complete = TRUE)
  contrast_basis <- full_basis[, -1, drop = FALSE]
  mean_data <- crossprod(mean_basis, dat_mat)
  contrast_data <- t(contrast_basis) %*% dat_mat

  n_params <- (length(comp_maps$item) + length(comp_maps$facet) + 1) * p * (p + 1) / 2
  mvow_log(sprintf("=== Multivariate Full Model ===\n"))
  mvow_log(sprintf("Traits: %d | Conditions: %d | Manifest: %d | Params: %.0f\n",
                    p, C, p * C, n_params + p))

  contrast_model <- mvow_build_ram_group(
    "contrast", contrast_data, traits, base_conds,
    kernels, comp_maps, sv,
    include_facet_part = FALSE, n_items = n_items, n_R = n_R
  )
  mean_model <- mvow_build_ram_group(
    "grandmean", mean_data, traits, base_conds,
    kernels, comp_maps, sv,
    include_facet_part = TRUE, n_items = n_items, n_R = n_R
  )

  model <- mxModel("mv_gtheory_full", contrast_model, mean_model,
                    mxFitFunctionMultigroup(c("contrast", "grandmean")))

  mvow_log("Fitting...\n")
  fit <- mvow_timed("  Multivariate full optimization",
    if (!is.null(extra_tries)) {
      mxTryHard(model, extraTries = extra_tries,
                OKstatuscodes = c(0L, 6L),
                jitterDistrib = "runif", loc = 1, scale = 0.25)
    } else {
      mxRun(model)
    }
  )

  status <- fit$output$status$code
  minus2LL <- fit$output$Minus2LogLikelihood
  params <- omxGetParameters(fit)
  n_par <- length(params)
  mvow_log(sprintf("  Status: %d | -2LL: %.4f\n", status, minus2LL))

  # Extract p×p covariance matrices
  all_comp_names <- c(names(comp_maps$item), names(comp_maps$facet), "epsilon")
  vc_mats <- list()
  for (nm in all_comp_names) {
    m <- matrix(0, p, p, dimnames = list(traits, traits))
    for (ti in seq_len(p)) {
      m[ti, ti] <- unname(params[paste0("s_", nm, "_", traits[ti], "_", traits[ti])])
    }
    if (p > 1) {
      for (ti in 1:(p - 1)) {
        for (tj in (ti + 1):p) {
          m[ti, tj] <- unname(params[paste0("s_", nm, "_", traits[ti], "_", traits[tj])])
          m[tj, ti] <- m[ti, tj]
        }
      }
    }
    if (nm == "epsilon") m <- m * n_R
    vc_mats[[nm]] <- m
  }

  # Per-trait variance summary
  diag_rows <- list()
  for (nm in all_comp_names) {
    for (ti in seq_len(p)) {
      diag_rows <- c(diag_rows, list(tibble(
        trait = traits[ti], component = nm, variance = vc_mats[[nm]][ti, ti]
      )))
    }
  }
  var_summary <- bind_rows(diag_rows) %>%
    group_by(trait) %>%
    mutate(total_var = sum(variance), prop_total = variance / total_var) %>%
    ungroup() %>%
    arrange(trait, desc(variance))

  list(
    model_type  = "full",
    fit         = fit,
    status      = status,
    minus2LL    = minus2LL,
    AIC         = minus2LL + 2 * n_par,
    BIC         = minus2LL + n_par * log(n_items),
    vc_mats     = vc_mats,
    var_summary = var_summary,
    traits      = traits,
    params      = params
  )
}


# ---------------------------------------------------------------------------
# Run (full model only — nested requires BW)
# ---------------------------------------------------------------------------

mvow_run_all_models <- function(dat_bundle, extra_tries = 5L) {
  res <- mvow_fit_full(dat_bundle, extra_tries = extra_tries)

  mvow_log("\nPer-trait variance decomposition:\n")
  print(res$var_summary, n = 60, width = Inf)

  list(full = res)
}
