# ============================================================
# G-Theory via OpenMx RAM/Path Specification
#
# This file provides a path-based OpenMx version of the manuscript
# variance-component models. It reuses the same data preparation and
# kernel construction helpers from fit_gtheory_openmx.R, but converts
# each covariance component kernel into a RAM latent block using an
# eigendecomposition:
#
#   sigma^2 * K = L %*% diag(sigma^2) %*% t(L)
#
# where the columns of L are fixed loading patterns derived from the
# positive-eigenvalue basis of K. This keeps the model in mxPath()
# syntax while preserving the same implied covariance structure as the
# matrix-algebra version.
# ============================================================

if (!exists("prepare_openmx_data") ||
    !exists("build_component_maps") ||
    !exists("build_structure_matrices") ||
    !exists("sanitize_openmx_name")) {
  source("univariate/src/fit_gtheory_openmx.R")
}

suppressPackageStartupMessages(library(OpenMx))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tibble))


path_kernel_basis <- function(kernel, tol = 1e-8) {
  kernel_sym <- (kernel + t(kernel)) / 2
  eig <- eigen(kernel_sym, symmetric = TRUE)
  keep <- which(eig$values > tol)

  if (length(keep) == 0) {
    return(NULL)
  }

  basis <- eig$vectors[, keep, drop = FALSE] %*%
    diag(sqrt(eig$values[keep]), nrow = length(keep))
  basis
}


add_manifest_means_paths <- function(manifest_vars, free, value = 0, label = "grand_mean") {
  if (free) {
    mxPath(
      from = "one",
      to = manifest_vars,
      arrows = 1,
      free = TRUE,
      values = rep(value, length(manifest_vars)),
      labels = rep(label, length(manifest_vars))
    )
  } else {
    mxPath(
      from = "one",
      to = manifest_vars,
      arrows = 1,
      free = FALSE,
      values = rep(0, length(manifest_vars))
    )
  }
}


build_ram_group_model <- function(group_name,
                                  raw_data,
                                  cond_cols,
                                  kernels,
                                  comp_maps,
                                  sv,
                                  include_facet_part = FALSE,
                                  cov_type = "standard",
                                  n_items = 1L,
                                  n_R = 1L,
                                  tol = 1e-8) {
  manifest_vars <- cond_cols
  latent_vars <- character(0)
  path_parts <- list()

  add_component_paths <- function(component_name, kernel_name, scale = 1, lbound = 0) {
    kernel_basis <- path_kernel_basis(scale * kernels[[kernel_name]], tol = tol)
    if (is.null(kernel_basis)) {
      return(invisible(NULL))
    }

    n_basis <- ncol(kernel_basis)
    lv_names <- paste0("lv_", component_name, "_", seq_len(n_basis))
    latent_vars <<- c(latent_vars, lv_names)

    path_parts[[length(path_parts) + 1]] <<- mxPath(
      from = "one",
      to = lv_names,
      arrows = 1,
      free = FALSE,
      values = rep(0, n_basis)
    )

    path_parts[[length(path_parts) + 1]] <<- mxPath(
      from = lv_names,
      arrows = 2,
      free = TRUE,
      values = rep(if (!is.null(sv[[component_name]])) sv[[component_name]] else 0.01, n_basis),
      lbound = rep(lbound, n_basis),
      labels = rep(paste0("s2_", component_name), n_basis)
    )

    for (j in seq_len(n_basis)) {
      path_parts[[length(path_parts) + 1]] <<- mxPath(
        from = lv_names[[j]],
        to = manifest_vars,
        arrows = 1,
        free = FALSE,
        values = kernel_basis[, j]
      )
    }
  }

  for (component_name in names(comp_maps$item)) {
    add_component_paths(
      component_name = component_name,
      kernel_name = unname(comp_maps$item[[component_name]])
    )
  }

  if (cov_type == "clustered") {
    add_component_paths(
      component_name = "clustered",
      kernel_name = "K_clustered"
    )
  }

  if (include_facet_part) {
    for (component_name in names(comp_maps$facet)) {
      add_component_paths(
        component_name = component_name,
        kernel_name = unname(comp_maps$facet[[component_name]]),
        scale = n_items
      )
    }
  }

  path_parts[[length(path_parts) + 1]] <- mxPath(
    from = manifest_vars,
    arrows = 2,
    free = TRUE,
    values = rep(if (!is.null(sv[["epsilon"]])) sv[["epsilon"]] / n_R else 0.01, length(manifest_vars)),
    lbound = rep(1e-6, length(manifest_vars)),
    labels = rep("s2_epsilon", length(manifest_vars))
  )

  path_parts[[length(path_parts) + 1]] <- add_manifest_means_paths(
    manifest_vars = manifest_vars,
    free = include_facet_part,
    value = as.numeric(mean(raw_data)),
    label = "grand_mean"
  )

  mxModel(
    group_name,
    type = "RAM",
    manifestVars = manifest_vars,
    latentVars = latent_vars,
    mxData(observed = as.data.frame(raw_data), type = "raw"),
    path_parts
  )
}


fit_gtheory_openmx_path <- function(
    prep,
    model_name = "GTheoryPath",
    cov_type = "standard",
    eval_families = NULL,
    item_interaction_orders = 1L,
    facet_interaction_orders = 2L,
    extra_tries = 50L,
    use_matrix_warm_start = TRUE,
    original_data = NULL,
    item_col = "item_id",
    score_col = "severity",
    use_lme4_starts = FALSE,
    tol = 1e-8
) {
  safe_model_name <- sanitize_openmx_name(model_name)

  facet_info <- prep$facet_info
  crossed <- prep$crossed_facets
  cond_cols <- prep$cond_cols
  n_c <- length(cond_cols)
  n_R <- prep$n_R
  n_rg <- prep$n_rep_groups
  n_items <- prep$n_p

  facet_sizes <- sapply(crossed, function(f) facet_info[[f]]$n)
  short_map <- facet_short_map(crossed)
  names(facet_sizes) <- unname(short_map[crossed])

  cat(sprintf("=== OpenMx G-Theory (Path): %s ===\n", model_name))
  cat(sprintf("Model: %s | Facets: %s | Manifest: %d\n",
              cov_type, paste(crossed, collapse = " × "), n_c))
  if (!identical(model_name, safe_model_name)) {
    cat(sprintf("OpenMx model name sanitized to: %s\n", safe_model_name))
  }

  comp_maps <- build_component_maps(
    crossed_facets = crossed,
    item_interaction_orders = item_interaction_orders,
    facet_interaction_orders = facet_interaction_orders,
    n_rep_groups = n_rg
  )

  kernels <- build_structure_matrices(facet_sizes, n_rg, eval_families = eval_families)

  dat_mat <- as.matrix(prep$wide_data[, cond_cols, drop = FALSE])
  if (nrow(dat_mat) < 2) {
    stop("OpenMx fitting requires at least two items.")
  }

  mean_basis <- rep(1 / sqrt(n_items), n_items)
  full_basis <- qr.Q(qr(cbind(mean_basis, diag(n_items))), complete = TRUE)
  contrast_basis <- full_basis[, -1, drop = FALSE]
  mean_data <- crossprod(mean_basis, dat_mat)
  contrast_data <- t(contrast_basis) %*% dat_mat

  if (!is.null(prep$start_values)) {
    sv <- prep$start_values
    cat("Start values: oracle\n")
  } else {
    sv <- NULL

    if (isTRUE(use_lme4_starts) && !is.null(original_data)) {
      sv <- tryCatch(
        estimate_lme4_start_values(
          prep = prep,
          original_data = original_data,
          item_col = item_col,
          score_col = score_col,
          cov_type = cov_type,
          eval_families = eval_families,
          item_interaction_orders = item_interaction_orders,
          facet_interaction_orders = facet_interaction_orders,
          use_REML = FALSE
        ),
        error = function(e) {
          cat(sprintf("lme4 start estimation failed: %s\n", e$message))
          NULL
        }
      )
    }

    if (!is.null(sv)) {
      cat("Start values: lme4 warm start\n")
    } else {
      sv <- tryCatch(
        estimate_raw_start_values(
          prep,
          cov_type = cov_type,
          eval_families = eval_families,
          item_interaction_orders = item_interaction_orders,
          facet_interaction_orders = facet_interaction_orders
        ),
        error = function(e) {
          cat(sprintf("Raw-data start estimation failed: %s\n", e$message))
          NULL
        }
      )

      if (!is.null(sv)) {
        cat("Start values: raw-data effect decomposition\n")
      } else {
        sv <- estimate_simple_start_values(
          prep,
          cov_type = cov_type,
          item_interaction_orders = item_interaction_orders,
          facet_interaction_orders = facet_interaction_orders
        )
        cat("Start values: method-of-moments (approximate fallback)\n")
      }
    }
  }

  run_path_fit <- function(model, extra_tries, scale = 0.25, wtgcsv = c("best", "prev", "initial")) {
    mxTryHard(
      model,
      extraTries = extra_tries,
      OKstatuscodes = c(0L, 6L),
      silent = TRUE,
      scale = scale,
      wtgcsv = wtgcsv,
      exhaustive = FALSE
    )
  }

  contrast_model <- build_ram_group_model(
    group_name = "contrast",
    raw_data = contrast_data,
    cond_cols = cond_cols,
    kernels = kernels,
    comp_maps = comp_maps,
    sv = sv,
    include_facet_part = FALSE,
    cov_type = cov_type,
    n_items = n_items,
    n_R = n_R,
    tol = tol
  )

  mean_model <- build_ram_group_model(
    group_name = "grandmean",
    raw_data = mean_data,
    cond_cols = cond_cols,
    kernels = kernels,
    comp_maps = comp_maps,
    sv = sv,
    include_facet_part = TRUE,
    cov_type = cov_type,
    n_items = n_items,
    n_R = n_R,
    tol = tol
  )

  model <- mxModel(
    safe_model_name,
    contrast_model,
    mean_model,
    mxFitFunctionMultigroup(c("contrast", "grandmean"))
  )

  cat("Fitting...\n")
  fit <- run_path_fit(model, extra_tries = extra_tries, scale = 0.20)

  status <- fit$output$status$code
  if (status == 6) {
    cat("Path model returned status 6; rerunning from previous solution with tighter jitter.\n")
    fit <- run_path_fit(fit, extra_tries = max(5L, extra_tries), scale = 0.05, wtgcsv = c("prev", "best"))
    status <- fit$output$status$code
  }

  if (status != 0 && use_matrix_warm_start && is.null(prep$start_values)) {
    cat("Path model still not status 0; fitting matrix OpenMx model for warm starts.\n")
    matrix_fit <- fit_gtheory_openmx(
      prep = prep,
      model_name = paste0(model_name, "_matrixWarm"),
      cov_type = cov_type,
      eval_families = eval_families,
      item_interaction_orders = item_interaction_orders,
      facet_interaction_orders = facet_interaction_orders,
      extra_tries = max(3L, min(extra_tries, 10L))
    )

    prep_warm <- prep
    prep_warm$start_values <- matrix_fit$vc_list

    cat("Refitting path model from matrix-based start values.\n")
    return(fit_gtheory_openmx_path(
      prep = prep_warm,
      model_name = model_name,
      cov_type = cov_type,
      eval_families = eval_families,
      item_interaction_orders = item_interaction_orders,
      facet_interaction_orders = facet_interaction_orders,
      extra_tries = max(5L, extra_tries),
      use_matrix_warm_start = FALSE,
      original_data = original_data,
      item_col = item_col,
      score_col = score_col,
      use_lme4_starts = FALSE,
      tol = tol
    ))
  }

  if (status != 0) warning(sprintf("Optimizer status: %d", status))
  cat(sprintf("Optimizer status: %d\n", status))

  params <- omxGetParameters(fit)

  vc_list <- list()
  all_component_names <- c(names(comp_maps$item), names(comp_maps$facet), "epsilon")
  if (cov_type == "clustered") {
    all_component_names <- c(all_component_names, "clustered")
  }

  for (component_name in all_component_names) {
    param_name <- paste0("s2_", component_name)
    value <- if (param_name %in% names(params)) as.numeric(params[[param_name]]) else 0
    if (component_name == "epsilon") {
      value <- value * n_R
    }
    vc_list[[component_name]] <- value
  }

  s2_tau <- vc_list[["tau"]]
  sigma2_delta <- 0
  item_error_components <- setdiff(names(comp_maps$item), "tau")
  for (component_name in item_error_components) {
    kernel_name <- unname(comp_maps$item[[component_name]])
    kernel_shorts <- strsplit(sub("^K_", "", kernel_name), "_")[[1]]
    denom <- prod(unname(facet_sizes[kernel_shorts]))
    sigma2_delta <- sigma2_delta + vc_list[[component_name]] / denom
  }
  if (cov_type == "clustered") {
    n_fam <- length(unique(eval_families))
    sigma2_delta <- sigma2_delta + vc_list[["clustered"]] / n_fam
  }
  sigma2_delta <- sigma2_delta + vc_list[["epsilon"]] / (prod(unname(facet_sizes)) * n_R)

  g_coef <- s2_tau / (s2_tau + sigma2_delta)

  minus2LL <- fit$output$Minus2LogLikelihood
  n_par <- length(params)
  aic <- minus2LL + 2 * n_par
  bic <- minus2LL + n_par * log(prep$n_p)

  obs_cov <- cov(dat_mat)
  exp_cov <- mxGetExpected(fit$contrast, "covariance")
  rownames(exp_cov) <- cond_cols
  colnames(exp_cov) <- cond_cols
  resid_cov <- obs_cov - exp_cov
  rmsr <- sqrt(mean(resid_cov[lower.tri(resid_cov, diag = TRUE)]^2))

  cat(sprintf("\nG-coefficient: %.4f\n", g_coef))
  cat(sprintf("RMSR: %.6f | -2LL: %.2f | AIC: %.2f | nPar: %d | Status: %d\n",
              rmsr, minus2LL, aic, n_par, status))

  list(
    fit = fit,
    params = params,
    vc_list = vc_list,
    s2_epsilon = vc_list[["epsilon"]],
    g_coef = g_coef,
    rmsr = rmsr,
    model_name = safe_model_name,
    cov_type = cov_type,
    specification = "RAM_path",
    obs_cov = obs_cov,
    exp_cov = exp_cov,
    model_fit = list(
      minus2LL = minus2LL,
      AIC = aic,
      BIC = bic,
      n_par = n_par,
      status = status
    ),
    design = list(
      crossed = crossed,
      facet_sizes = facet_sizes,
      nested = prep$rep_facets,
      nested_within = prep$nested_within %||% setNames(vector("list", length(prep$rep_facets)), prep$rep_facets),
      facet_roles = prep$facet_roles %||% setNames(rep("facet", length(crossed)), crossed),
      n_p = prep$n_p,
      n_R = n_R,
      n_per_group = prep$n_per_group,
      n_rep_groups = n_rg
    )
  )
}


fit_gtheory_openmx_path_design <- function(
    dat,
    facet_spec,
    item_col = "item_id",
    score_col = "severity",
    n_rep_groups = 1L,
    model_name = "GTheoryPath",
    cov_type = "standard",
    eval_families = NULL,
    item_interaction_orders = 1L,
    facet_interaction_orders = 2L,
    extra_tries = 50L,
    start_values = NULL,
    use_matrix_warm_start = TRUE,
    use_lme4_starts = FALSE,
    tol = 1e-8
) {
  prep <- prepare_openmx_data_design(
    dat = dat,
    facet_spec = facet_spec,
    item_col = item_col,
    score_col = score_col,
    n_rep_groups = n_rep_groups
  )
  prep$start_values <- start_values

  fit_gtheory_openmx_path(
    prep = prep,
    model_name = model_name,
    cov_type = cov_type,
    eval_families = eval_families,
    item_interaction_orders = item_interaction_orders,
    facet_interaction_orders = facet_interaction_orders,
    extra_tries = extra_tries,
    use_matrix_warm_start = use_matrix_warm_start,
    original_data = dat,
    item_col = item_col,
    score_col = score_col,
    use_lme4_starts = use_lme4_starts,
    tol = tol
  )
}
