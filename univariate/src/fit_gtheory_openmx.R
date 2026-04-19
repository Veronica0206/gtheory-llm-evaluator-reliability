# ============================================================
# G-Theory via OpenMx: Matched Variance-Component Models
#
# Fits the same variance decomposition used in the lme4 reference:
#   - item universe score (tau)
#   - facet main effects
#   - item x facet interactions
#   - facet x facet interactions
#   - residual
#   - optional clustered/shared bias term
#
# For balanced designs, the item-by-condition matrix is orthogonally
# decomposed into:
#   - n-1 item-contrast rows with covariance A
#   - 1 grand-mean row with covariance A + nB
#
# This lets OpenMx estimate both item-side and facet-side components
# without constructing a giant observation-level covariance matrix.
# ============================================================

library(OpenMx)
library(dplyr)
library(tidyr)
library(tibble)

mxOption(model = NULL, key = "Default optimizer", value = "CSOLNP", reset = FALSE)


# ============================================================
# Small utilities
# ============================================================

enumerate_subsets <- function(x, min_size = 1, max_size = length(x)) {
  if (length(x) == 0 || min_size > max_size) {
    return(list())
  }

  subsets <- list()
  idx <- 1
  for (size in seq.int(min_size, max_size)) {
    current <- combn(x, size, simplify = FALSE)
    for (subset in current) {
      subsets[[idx]] <- subset
      idx <- idx + 1
    }
  }

  subsets
}


`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0) y else x
}


facet_short_map <- function(crossed_facets) {
  facet_shorts <- substr(crossed_facets, 1, 1)
  if (anyDuplicated(facet_shorts)) {
    stop("OpenMx facet shorthand requires unique first-letter facet abbreviations.")
  }
  setNames(facet_shorts, crossed_facets)
}


canonical_component_name <- function(facets, item_side = FALSE) {
  if (length(facets) == 0) {
    return("tau")
  }
  if (item_side) {
    return(paste(c("tau", facets), collapse = "_"))
  }
  paste(facets, collapse = "_")
}


kernel_name_for_subset <- function(facets, short_map) {
  if (length(facets) == 0) {
    return("K_tau")
  }
  paste0("K_", paste(unname(short_map[facets]), collapse = "_"))
}


sanitize_item_orders <- function(item_interaction_orders, n_crossed, n_rep_groups) {
  max_estimable <- if (n_rep_groups > 1) n_crossed else max(n_crossed - 1, 1)
  if (is.null(item_interaction_orders)) {
    return(seq_len(max_estimable))
  }

  orders <- sort(unique(as.integer(item_interaction_orders)))
  orders <- orders[orders >= 1 & orders <= max_estimable]
  orders
}


sanitize_facet_orders <- function(facet_interaction_orders, n_crossed) {
  if (is.null(facet_interaction_orders)) {
    return(if (n_crossed >= 2) 2:n_crossed else integer(0))
  }

  orders <- sort(unique(as.integer(facet_interaction_orders)))
  orders <- orders[orders >= 2 & orders <= n_crossed]
  orders
}


build_component_maps <- function(crossed_facets,
                                 item_interaction_orders = 1L,
                                 facet_interaction_orders = 2L,
                                 n_rep_groups = 1) {
  n_crossed <- length(crossed_facets)
  short_map <- facet_short_map(crossed_facets)

  item_orders <- sanitize_item_orders(item_interaction_orders, n_crossed, n_rep_groups)
  facet_orders <- sanitize_facet_orders(facet_interaction_orders, n_crossed)

  item_components <- c(tau = kernel_name_for_subset(character(0), short_map))
  for (order in item_orders) {
    for (subset in combn(crossed_facets, order, simplify = FALSE)) {
      item_components[[canonical_component_name(subset, item_side = TRUE)]] <- kernel_name_for_subset(subset, short_map)
    }
  }

  facet_components <- c()
  for (facet in crossed_facets) {
    facet_components[[facet]] <- kernel_name_for_subset(facet, short_map)
  }
  for (order in facet_orders) {
    for (subset in combn(crossed_facets, order, simplify = FALSE)) {
      facet_components[[canonical_component_name(subset, item_side = FALSE)]] <- kernel_name_for_subset(subset, short_map)
    }
  }

  list(
    item = item_components,
    facet = facet_components,
    short_map = short_map,
    item_orders = item_orders,
    facet_orders = facet_orders
  )
}


build_true_start_values <- function(prep, variance_components, clustered = NULL,
                                    item_interaction_orders = 1L,
                                    facet_interaction_orders = 2L) {
  comp_maps <- build_component_maps(
    crossed_facets = prep$crossed_facets,
    item_interaction_orders = item_interaction_orders,
    facet_interaction_orders = facet_interaction_orders,
    n_rep_groups = prep$n_rep_groups
  )

  starts <- list()
  active_names <- c(names(comp_maps$item), names(comp_maps$facet), "epsilon")
  for (nm in active_names) {
    starts[[nm]] <- if (nm %in% names(variance_components)) variance_components[[nm]] else 0
  }

  if (!is.null(clustered)) {
    starts[["clustered"]] <- if (!is.null(clustered$s2)) clustered$s2 else 0
  }

  starts
}


sanitize_openmx_name <- function(name) {
  safe <- gsub("[^A-Za-z0-9_]", "_", name)
  safe <- gsub("_+", "_", safe)
  safe <- gsub("^_+|_+$", "", safe)

  if (!nzchar(safe)) {
    safe <- "OpenMxModel"
  }
  if (!grepl("^[A-Za-z_]", safe)) {
    safe <- paste0("M_", safe)
  }

  safe
}


matrix_interaction_order <- function(mat_name) {
  if (mat_name %in% c("K_tau", "K_identity", "K_clustered")) {
    return(NA_integer_)
  }
  length(strsplit(sub("^K_", "", mat_name), "_")[[1]])
}


kernel_to_facets <- function(kernel_name, short_to_full) {
  if (kernel_name == "K_tau") {
    return(character(0))
  }
  shorts <- strsplit(sub("^K_", "", kernel_name), "_")[[1]]
  unname(short_to_full[shorts])
}


apply_mode_matrix <- function(arr, mat, dim_idx) {
  perm <- c(dim_idx, setdiff(seq_along(dim(arr)), dim_idx))
  arr_perm <- aperm(arr, perm)
  arr_mat <- matrix(arr_perm, nrow = dim(arr_perm)[1])
  res_mat <- mat %*% arr_mat
  res_arr <- array(res_mat, dim = c(nrow(mat), dim(arr_perm)[-1]))
  aperm(res_arr, order(perm))
}


project_balanced_component <- function(arr, subset_dims, all_dims) {
  out <- arr
  dim_sizes <- setNames(dim(arr), all_dims)

  for (d in all_dims) {
    n_d <- dim_sizes[[d]]
    P_d <- matrix(1 / n_d, n_d, n_d)
    H_d <- diag(n_d) - P_d
    projector <- if (d %in% subset_dims) H_d else P_d
    out <- apply_mode_matrix(out, projector, match(d, all_dims))
  }

  out
}


build_cell_mean_array <- function(prep) {
  n_items <- prep$n_p
  crossed <- prep$crossed_facets
  facet_sizes <- vapply(crossed, function(f) prep$facet_info[[f]]$n, numeric(1))
  arr_dims <- c(item = n_items, setNames(facet_sizes, crossed))
  if (prep$n_rep_groups > 1) {
    arr_dims <- c(arr_dims, rep_group = prep$n_rep_groups)
  }

  arr <- array(NA_real_, dim = arr_dims)
  dat_mat <- as.matrix(prep$wide_data[, prep$cond_cols, drop = FALSE])

  # Validate: balanced design requires no missing cells
  n_missing <- sum(is.na(dat_mat))
  if (n_missing > 0) {
    warning(sprintf(
      "build_cell_mean_array: %d missing values in wide data (%d items x %d conditions). Estimates may be biased.",
      n_missing, nrow(dat_mat), ncol(dat_mat)
    ))
  }

  for (j in seq_along(prep$cond_cols)) {
    map_row <- prep$cond_map[j, ]
    fixed_idx <- integer(0)

    for (f in crossed) {
      fixed_idx <- c(
        fixed_idx,
        match(map_row[[paste0(f, "_label")]], prep$facet_info[[f]]$labels)
      )
    }
    if (prep$n_rep_groups > 1) {
      fixed_idx <- c(
        fixed_idx,
        match(map_row[["rep_group"]], paste0("g", seq_len(prep$n_rep_groups)))
      )
    }

    subs <- c(list(arr, seq_len(n_items)), as.list(fixed_idx))
    subs$value <- dat_mat[, j]
    arr <- do.call("[<-", subs)
  }

  arr
}


project_clustered_effect <- function(arr, crossed_facets, eval_families) {
  if (is.null(eval_families) || !("evaluator" %in% crossed_facets)) {
    return(NULL)
  }

  arr_dims <- names(dim(arr))
  eval_dim_idx <- match("evaluator", arr_dims)
  if (is.na(eval_dim_idx)) {
    return(NULL)
  }

  fam_levels <- sort(unique(eval_families))
  n_fam <- length(fam_levels)
  fam_dims <- dim(arr)
  fam_dims[eval_dim_idx] <- n_fam
  fam_arr <- array(0, dim = fam_dims)
  fam_dim_names <- arr_dims
  fam_dim_names[eval_dim_idx] <- "family"
  names(dim(fam_arr)) <- fam_dim_names

  for (f_idx in seq_along(fam_levels)) {
    eval_idx <- which(eval_families == fam_levels[[f_idx]])
    fam_slice <- apply(arr, setdiff(seq_along(dim(arr)), eval_dim_idx), function(x) mean(x[eval_idx]))

    subs <- vector("list", length(dim(arr)) + 1)
    subs[[1]] <- fam_arr
    for (d in seq_along(dim(arr))) {
      subs[[d + 1]] <- if (d == eval_dim_idx) f_idx else TRUE
    }
    subs$value <- fam_slice
    fam_arr <- do.call("[<-", subs)
  }

  fam_effect <- project_balanced_component(
    fam_arr,
    subset_dims = c("item", "family"),
    all_dims = names(dim(fam_arr))
  )

  expanded <- array(0, dim = dim(arr))
  names(dim(expanded)) <- arr_dims
  for (e_idx in seq_along(eval_families)) {
    fam_idx <- match(eval_families[[e_idx]], fam_levels)
    subs_src <- vector("list", length(dim(arr)) + 1)
    subs_src[[1]] <- fam_effect
    subs_dst <- vector("list", length(dim(arr)) + 1)
    subs_dst[[1]] <- expanded
    for (d in seq_along(dim(arr))) {
      if (d == eval_dim_idx) {
        subs_src[[d + 1]] <- fam_idx
        subs_dst[[d + 1]] <- e_idx
      } else {
        subs_src[[d + 1]] <- TRUE
        subs_dst[[d + 1]] <- TRUE
      }
    }
    slice_val <- do.call("[", c(subs_src, drop = FALSE))
    subs_dst$value <- slice_val
    expanded <- do.call("[<-", subs_dst)
  }

  expanded
}


calibrate_start_values <- function(starts,
                                   total_variance,
                                   n_R = 1L,
                                   reserve_frac = 0.05,
                                   tau_floor = 1e-4,
                                   other_floor = 1e-5,
                                   epsilon_floor = 1e-4) {
  if (length(starts) == 0) {
    return(starts)
  }

  total_variance <- as.numeric(total_variance[[1]])
  if (!is.finite(total_variance) || is.na(total_variance) || total_variance <= 0) {
    total_variance <- 1e-3
  }

  n_R <- max(as.numeric(n_R[[1]]), 1)
  nm_all <- names(starts)

  component_names <- setdiff(nm_all, "epsilon")
  component_vals <- vapply(component_names, function(nm) {
    floor_val <- if (identical(nm, "tau")) tau_floor else other_floor
    val <- as.numeric(starts[[nm]][[1]])
    if (!is.finite(val) || is.na(val)) {
      val <- 0
    }
    max(val, floor_val)
  }, numeric(1))

  max_component_total <- max(total_variance * (1 - reserve_frac), tau_floor)
  if (length(component_vals) > 0 && sum(component_vals) > max_component_total) {
    component_vals <- component_vals * (max_component_total / sum(component_vals))
  }

  for (nm in names(component_vals)) {
    starts[[nm]] <- component_vals[[nm]]
  }

  remaining_var <- max(total_variance - sum(component_vals), total_variance * reserve_frac, epsilon_floor / n_R)
  starts[["epsilon"]] <- max(remaining_var * n_R, epsilon_floor)

  starts
}


estimate_raw_start_values <- function(prep,
                                      cov_type = "standard",
                                      eval_families = NULL,
                                      item_interaction_orders = 1L,
                                      facet_interaction_orders = 2L) {
  mean_arr <- build_cell_mean_array(prep)
  arr_dims <- names(dim(mean_arr))
  crossed <- prep$crossed_facets
  comp_maps <- build_component_maps(
    crossed_facets = crossed,
    item_interaction_orders = item_interaction_orders,
    facet_interaction_orders = facet_interaction_orders,
    n_rep_groups = prep$n_rep_groups
  )
  short_to_full <- setNames(names(comp_maps$short_map), comp_maps$short_map)

  starts <- list()
  modeled_sum <- array(0, dim = dim(mean_arr))
  names(dim(modeled_sum)) <- arr_dims
  centered_arr <- mean_arr - mean(mean_arr, na.rm = TRUE)
  total_variance <- mean(centered_arr^2, na.rm = TRUE)

  for (component_name in names(comp_maps$item)) {
    kernel_name <- unname(comp_maps$item[[component_name]])
    subset_dims <- c("item", kernel_to_facets(kernel_name, short_to_full))
    effect_arr <- project_balanced_component(mean_arr, subset_dims, arr_dims)
    starts[[component_name]] <- mean(effect_arr^2, na.rm = TRUE)
    modeled_sum <- modeled_sum + effect_arr
  }

  for (component_name in names(comp_maps$facet)) {
    kernel_name <- unname(comp_maps$facet[[component_name]])
    subset_dims <- kernel_to_facets(kernel_name, short_to_full)
    effect_arr <- project_balanced_component(mean_arr, subset_dims, arr_dims)
    starts[[component_name]] <- mean(effect_arr^2, na.rm = TRUE)
    modeled_sum <- modeled_sum + effect_arr
  }

  if (cov_type == "clustered") {
    cluster_arr <- project_clustered_effect(mean_arr, crossed, eval_families)
    if (!is.null(cluster_arr)) {
      starts[["clustered"]] <- mean(cluster_arr^2, na.rm = TRUE)
    }
  }

  resid_arr <- mean_arr - modeled_sum
  starts[["epsilon"]] <- mean(resid_arr^2, na.rm = TRUE) * prep$n_R

  calibrate_start_values(
    starts = starts,
    total_variance = total_variance,
    n_R = prep$n_R
  )
}


estimate_simple_start_values <- function(prep,
                                         cov_type = "standard",
                                         item_interaction_orders = 1L,
                                         facet_interaction_orders = 2L) {
  dat_mat <- as.matrix(prep$wide_data[, prep$cond_cols, drop = FALSE])
  centered_cov <- cov(dat_mat)
  comp_maps <- build_component_maps(
    crossed_facets = prep$crossed_facets,
    item_interaction_orders = item_interaction_orders,
    facet_interaction_orders = facet_interaction_orders,
    n_rep_groups = prep$n_rep_groups
  )

  sv <- list()
  diag_mean <- mean(diag(centered_cov))
  offdiag_mean <- if (ncol(dat_mat) > 1) mean(centered_cov[lower.tri(centered_cov)]) else diag_mean
  tau_start <- max(offdiag_mean, 1e-4)
  other_start <- max((diag_mean - tau_start) / max(length(comp_maps$item) + length(comp_maps$facet), 1), 1e-4)

  sv[["tau"]] <- tau_start
  for (nm in setdiff(names(comp_maps$item), "tau")) {
    sv[[nm]] <- other_start
  }
  for (nm in names(comp_maps$facet)) {
    sv[[nm]] <- other_start
  }
  sv[["epsilon"]] <- max(other_start * prep$n_R, 1e-4)
  if (cov_type == "clustered") {
    sv[["clustered"]] <- other_start
  }

  calibrate_start_values(
    starts = sv,
    total_variance = diag_mean,
    n_R = prep$n_R
  )
}


estimate_lme4_start_values <- function(prep,
                                       original_data = NULL,
                                       item_col = "item_id",
                                       score_col = "severity",
                                       cov_type = "standard",
                                       eval_families = NULL,
                                       item_interaction_orders = 1L,
                                       facet_interaction_orders = 2L,
                                       use_REML = FALSE) {
  if (is.null(original_data)) {
    stop("`original_data` is required for lme4-based start values.")
  }

  if (!exists("fit_gtheory_lme4")) {
    stop("`fit_gtheory_lme4()` is not available.")
  }

  clustered <- NULL
  if (identical(cov_type, "clustered") && !is.null(eval_families)) {
    family_col <- ".openmx_eval_family"
    original_data[[family_col]] <- eval_families[match(original_data[["evaluator"]], prep$facet_info[["evaluator"]]$levels)]
    clustered <- list(family_col = family_col)
  }

  fit_lme4 <- fit_gtheory_lme4(
    dat = original_data,
    item_col = item_col,
    score_col = score_col,
    crossed_facets = prep$crossed_facets,
    clustered = clustered,
    item_interaction_orders = item_interaction_orders,
    facet_interaction_orders = facet_interaction_orders,
    use_REML = use_REML
  )

  vc <- fit_lme4$vc
  starts <- list(
    tau = vc[["s2_tau"]] %||% 0,
    epsilon = vc[["s2_epsilon"]] %||% 0
  )

  for (facet in prep$crossed_facets) {
    starts[[facet]] <- vc[[paste0("s2_", facet)]] %||% 0
    starts[[paste0("tau_", facet)]] <- vc[[paste0("s2_tau_", facet)]] %||% 0
  }

  facet_orders <- sanitize_facet_orders(facet_interaction_orders, length(prep$crossed_facets))
  if (length(facet_orders) > 0) {
    for (subset in enumerate_subsets(prep$crossed_facets, min(facet_orders), max(facet_orders))) {
      if (length(subset) %in% facet_orders) {
        starts[[paste(subset, collapse = "_")]] <- vc[[paste0("s2_", paste(subset, collapse = "_"))]] %||% 0
      }
    }
  }

  item_orders <- sanitize_item_orders(item_interaction_orders, length(prep$crossed_facets), prep$n_rep_groups)
  if (length(item_orders) > 0) {
    for (subset in enumerate_subsets(prep$crossed_facets, min(item_orders), max(item_orders))) {
      if (length(subset) %in% item_orders) {
        starts[[paste(c("tau", subset), collapse = "_")]] <- vc[[paste0("s2_tau_", paste(subset, collapse = "_"))]] %||% 0
      }
    }
  }

  if (identical(cov_type, "clustered")) {
    starts[["clustered"]] <- vc[["s2_clustered"]] %||% 0
  }

  starts
}


# ============================================================
# Step 0: Flexible design specification helpers
# ============================================================

normalize_openmx_design <- function(facet_spec, n_rep_groups = 1L) {
  if (is.null(names(facet_spec)) || any(names(facet_spec) == "")) {
    stop("`facet_spec` must be a named vector or named list.")
  }

  get_role <- function(x) {
    if (is.list(x)) {
      role <- x$role %||% x$type
    } else {
      role <- x
    }
    as.character(role[[1]])
  }

  get_within <- function(x) {
    if (is.list(x) && !is.null(x$within)) as.character(x$within) else NULL
  }

  roles <- vapply(facet_spec, get_role, character(1))
  valid_roles <- c("facet", "nested")
  if (!all(roles %in% valid_roles)) {
    bad <- unique(roles[!roles %in% valid_roles])
    stop("Unsupported role(s) in `facet_spec`: ", paste(bad, collapse = ", "),
         ". Valid roles are: ", paste(valid_roles, collapse = ", "), ".")
  }

  crossed_facets <- names(roles)[roles == "facet"]
  rep_facets <- names(roles)[roles == "nested"]
  nested_within <- setNames(lapply(facet_spec, get_within), names(facet_spec))

  if (length(crossed_facets) == 0) {
    stop("At least one facet must have role = 'facet'.")
  }

  for (f in rep_facets) {
    within_f <- nested_within[[f]]
    if (is.null(within_f)) {
      nested_within[[f]] <- crossed_facets
      next
    }

    bad_within <- setdiff(within_f, crossed_facets)
    if (length(bad_within) > 0) {
      stop(
        "Nested facet `", f, "` declares `within = ",
        paste(within_f, collapse = ", "),
        "` but the following variables are not modeled facets: ",
        paste(bad_within, collapse = ", "), "."
      )
    }
  }

  unsupported_nested <- names(which(vapply(rep_facets, function(f) {
    within_f <- nested_within[[f]]
    !is.null(within_f) && !setequal(sort(within_f), sort(crossed_facets))
  }, logical(1))))

  if (length(unsupported_nested) > 0) {
    warning(
      "Custom `within` nesting is not yet modeled separately in the OpenMx kernel code. ",
      "The following nested facets will be treated as replications within all modeled facets: ",
      paste(unsupported_nested, collapse = ", "), "."
    )
  }

  list(
    crossed_facets = crossed_facets,
    rep_facets = rep_facets,
    n_rep_groups = as.integer(n_rep_groups),
    facet_roles = roles,
    nested_within = nested_within
  )
}


summarize_openmx_design <- function(design) {
  nested_labels <- vapply(names(design$nested_within), function(f) {
    within_f <- design$nested_within[[f]]
    if (is.null(within_f)) {
      sprintf("%s [nested]", f)
    } else {
      sprintf("%s [nested within %s]", f, paste(within_f, collapse = " x "))
    }
  }, character(1))

  list(
    facets = design$crossed_facets,
    nested = design$rep_facets,
    nested_labels = nested_labels
  )
}


prepare_openmx_data_design <- function(
    dat,
    facet_spec,
    item_col = "item_id",
    score_col = "severity",
    n_rep_groups = 1L
) {
  design <- normalize_openmx_design(facet_spec, n_rep_groups = n_rep_groups)

  prep <- prepare_openmx_data(
    dat = dat,
    item_col = item_col,
    score_col = score_col,
    crossed_facets = design$crossed_facets,
    rep_facets = design$rep_facets,
    n_rep_groups = design$n_rep_groups
  )

  prep$facet_roles <- design$facet_roles
  prep$nested_within <- design$nested_within
  prep$facet_spec <- facet_spec
  prep$design_summary <- summarize_openmx_design(design)
  prep
}


# ============================================================
# Step 1: Prepare data
# ============================================================

prepare_openmx_data <- function(
    dat,
    item_col       = "item_id",
    score_col      = "severity",
    crossed_facets = c("evaluator", "prompt_type"),
    rep_facets     = c("temperature", "seed"),
    n_rep_groups   = 2
) {

  # Validate columns
  needed <- c(score_col, item_col, crossed_facets, rep_facets)
  missing <- setdiff(needed, names(dat))
  if (length(missing) > 0) stop("Missing columns: ", paste(missing, collapse = ", "))

  df <- dat
  names(df)[names(df) == score_col] <- "score"
  names(df)[names(df) == item_col] <- "item"

  # Assign short labels to each crossed facet
  facet_info <- list()
  for (f in crossed_facets) {
    levels_f <- sort(unique(df[[f]]))
    n_f <- length(levels_f)
    short <- paste0(substr(f, 1, 1), seq_len(n_f))
    label_map <- setNames(short, as.character(levels_f))
    df[[paste0(f, "_label")]] <- unname(label_map[as.character(df[[f]])])
    facet_info[[f]] <- list(levels = levels_f, n = n_f, labels = short, map = label_map)
  }

  # Handle replications
  if (length(rep_facets) > 0 && n_rep_groups > 1) {
    cell_counts <- df %>%
      group_by(across(all_of(c("item", crossed_facets)))) %>%
      summarise(n = n(), .groups = "drop")

    unique_counts <- unique(cell_counts$n)
    if (length(unique_counts) != 1) {
      stop("Replication counts vary across item x crossed-facet cells; equal replication is required.")
    }

    n_R_total <- unique_counts[[1]]
    if (n_R_total %% n_rep_groups != 0) {
      stop(sprintf(
        "Total replications per cell (%d) must divide evenly by n_rep_groups (%d).",
        n_R_total, n_rep_groups
      ))
    }

    n_per_group <- n_R_total / n_rep_groups

    # Assign rep_groups within each (item × crossed facets) cell
    df <- df %>%
      group_by(across(all_of(c("item", crossed_facets)))) %>%
      arrange(across(all_of(rep_facets)), .by_group = TRUE) %>%
      mutate(rep_idx = row_number(),
             rep_group = paste0("g", ((rep_idx - 1) %% n_rep_groups) + 1)) %>%
      ungroup()

    # Aggregate to item-level means within (crossed facets × rep_group)
    group_cols <- c("item", paste0(crossed_facets, "_label"), "rep_group")
    item_means <- df %>%
      group_by(across(all_of(group_cols))) %>%
      summarise(mean_score = mean(score, na.rm = TRUE), .groups = "drop")

    # Build condition label
    label_cols <- c(paste0(crossed_facets, "_label"), "rep_group")
    item_means$condition <- apply(item_means[, label_cols], 1, paste, collapse = "_")

    # Pivot wide
    wide <- item_means %>%
      select(item, condition, mean_score) %>%
      pivot_wider(names_from = condition, values_from = mean_score)

  } else {
    # No replication facets or n_rep_groups = 1: each cell is a single observation
    n_rep_groups <- 1
    n_R_total <- 1
    n_per_group <- 1

    # Build condition from crossed facet labels
    label_cols <- paste0(crossed_facets, "_label")
    df$condition <- apply(df[, label_cols, drop = FALSE], 1, paste, collapse = "_")

    # For no rep_facets: each (item, condition) should be unique
    # For rep_facets with n_rep_groups=1: average over rep_facets
    if (length(rep_facets) > 0) {
      n_R_total <- df %>%
        group_by(across(all_of(c("item", crossed_facets)))) %>%
        summarise(n = n(), .groups = "drop") %>%
        pull(n) %>%
        mean()
      n_per_group <- n_R_total

      wide <- df %>%
        group_by(item, condition) %>%
        summarise(mean_score = mean(score, na.rm = TRUE), .groups = "drop") %>%
        pivot_wider(names_from = condition, values_from = mean_score)
    } else {
      wide <- df %>%
        select(item, condition, score) %>%
        rename(mean_score = score) %>%
        pivot_wider(names_from = condition, values_from = mean_score)
    }
  }

  cond_cols <- sort(setdiff(names(wide), "item"))
  wide <- wide %>% select(item, all_of(cond_cols))
  n_p <- nrow(wide)

  # Build condition map
  cond_map <- tibble(condition = cond_cols)
  for (f in crossed_facets) {
    short_prefix <- substr(f, 1, 1)
    cond_map[[paste0(f, "_label")]] <- sub(
      paste0(".*?(", short_prefix, "[0-9]+).*"), "\\1", cond_cols
    )
  }
  if (n_rep_groups > 1) {
    cond_map$rep_group <- sub(".*_(g[0-9]+)$", "\\1", cond_cols)
  }

  cat("=== OpenMx G-Theory: Data Preparation ===\n")
  cat(sprintf("Items:          %d\n", n_p))
  for (f in crossed_facets) {
    fi <- facet_info[[f]]
    cat(sprintf("%-15s %d levels\n", paste0(f, ":"), fi$n))
  }
  cat(sprintf("Rep facets:     %s\n", if (length(rep_facets) > 0) paste(rep_facets, collapse = ", ") else "none"))
  cat(sprintf("Total reps/cell: %g\n", n_R_total))
  cat(sprintf("Rep groups:     %d (n_per_group = %g)\n", n_rep_groups, n_per_group))
  cat(sprintf("Manifest vars:  %d\n", length(cond_cols)))
  cat(sprintf("Data:           %d items x %d manifest\n\n", n_p, length(cond_cols)))

  list(
    wide_data     = as.data.frame(wide),
    cond_cols     = cond_cols,
    cond_map      = cond_map,
    facet_info    = facet_info,
    crossed_facets = crossed_facets,
    rep_facets    = rep_facets,
    n_p           = n_p,
    n_R           = n_R_total,
    n_rep_groups  = n_rep_groups,
    n_per_group   = n_per_group,
    start_values  = NULL
  )
}


# ============================================================
# Step 2: Build structure matrices (general)
# ============================================================

build_structure_matrices <- function(facet_sizes, n_g = 1, eval_families = NULL) {
  # facet_sizes: named vector of crossed-facet sizes in crossed-facet order
  # n_g: number of rep_groups (1 if no grouped replications)
  # Returns condition-space kernels used by both item-side and facet-side terms

  fnames <- names(facet_sizes)
  K <- length(fnames)

  # Pre-build I and J matrices for each facet + rep_group
  I_mats <- lapply(facet_sizes, diag)
  J_mats <- lapply(facet_sizes, function(n) matrix(1, n, n))
  Ig <- diag(n_g)
  Jg <- matrix(1, n_g, n_g)

  # Enumerate all 2^K subsets of facets
  mats <- list()

  for (bits in 0:(2^K - 1)) {
    # Decode which facets are "in" the subset (use I) vs "out" (use J)
    in_subset <- c()
    kron_list <- list()
    for (k in seq_len(K)) {
      if (bitwAnd(bits, 2^(k-1)) > 0) {
        kron_list[[k]] <- I_mats[[k]]
        in_subset <- c(in_subset, fnames[k])
      } else {
        kron_list[[k]] <- J_mats[[k]]
      }
    }

    # Build name
    if (length(in_subset) == 0) {
      mat_name <- "K_tau"
    } else {
      mat_name <- paste0("K_", paste(in_subset, collapse = "_"))
    }

    # Build Kronecker product
    result <- kron_list[[1]]
    if (K >= 2) {
      for (k in 2:K) {
        result <- kronecker(result, kron_list[[k]])
      }
    }
    if (n_g > 1) {
      result <- kronecker(result, Jg)  # cross-group covariance
    }

    mats[[mat_name]] <- result
  }

  # Residual: always the identity matrix (diagonal)
  n_total <- prod(facet_sizes) * max(n_g, 1)
  mats[["K_identity"]] <- diag(n_total)

  # Clustered evaluator structure
  if (!is.null(eval_families)) {
    n_e <- facet_sizes[1]  # assumes evaluator is the first facet
    B_cluster <- matrix(0, n_e, n_e)
    for (j in seq_len(n_e)) {
      for (jp in seq_len(n_e)) {
        if (eval_families[j] == eval_families[jp]) {
          B_cluster[j, jp] <- 1
        }
      }
    }
    # Kronecker with J for all other facets and rep_groups
    result <- B_cluster
    if (K >= 2) {
      for (k in 2:K) {
        result <- kronecker(result, J_mats[[k]])
      }
    }
    if (n_g > 1) {
      result <- kronecker(result, Jg)
    }
    mats[["K_clustered"]] <- result
  }

  mats
}


# ============================================================
# Step 3: Fit model
# ============================================================

fit_gtheory_openmx <- function(
    prep,
    model_name     = "GTheory",
    cov_type       = "standard",   # "standard", "clustered", "unstructured"
    eval_families  = NULL,
    item_interaction_orders = 1L,
    facet_interaction_orders = 2L,
    extra_tries = 50L,
    original_data = NULL,
    item_col = "item_id",
    score_col = "severity",
    use_lme4_starts = FALSE
) {
  safe_model_name <- sanitize_openmx_name(model_name)

  facet_info <- prep$facet_info
  crossed    <- prep$crossed_facets
  n_pg       <- prep$n_per_group
  n_rg       <- prep$n_rep_groups
  n_R        <- prep$n_R
  cond_cols  <- prep$cond_cols
  n_c        <- length(cond_cols)
  K          <- length(crossed)
  n_items     <- prep$n_p

  # Facet sizes in order
  facet_sizes <- sapply(crossed, function(f) facet_info[[f]]$n)
  short_map <- facet_short_map(crossed)
  names(facet_sizes) <- unname(short_map[crossed])

  cat(sprintf("=== OpenMx G-Theory: %s ===\n", model_name))
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

  # Build condition kernels
  kernels <- build_structure_matrices(facet_sizes, n_rg, eval_families = eval_families)

  # Data matrix: rows are items, columns are condition means
  dat_mat <- as.matrix(prep$wide_data[, cond_cols, drop = FALSE])
  if (nrow(dat_mat) < 2) {
    stop("OpenMx fitting requires at least two items.")
  }

  # Orthonormal transform across item rows:
  # one grand-mean row with covariance A + n*B, and n-1 contrast rows with covariance A.
  mean_basis <- rep(1 / sqrt(n_items), n_items)
  full_basis <- qr.Q(qr(cbind(mean_basis, diag(n_items))), complete = TRUE)
  contrast_basis <- full_basis[, -1, drop = FALSE]
  mean_data <- crossprod(mean_basis, dat_mat)
  contrast_data <- t(contrast_basis) %*% dat_mat

  centered_cov <- cov(dat_mat)

  # ---- Start values (method of moments or oracle) ----
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

  build_group_model <- function(group_name, raw_data, include_facet_part = FALSE) {
    model_parts <- list()

    add_kernel <- function(kernel_name) {
      model_parts[[length(model_parts) + 1]] <<- mxMatrix(
        "Full", n_c, n_c, free = FALSE, values = kernels[[kernel_name]], name = kernel_name
      )
    }

    add_parameter <- function(component_name, start_value, lbound = 0) {
      model_parts[[length(model_parts) + 1]] <<- mxMatrix(
        "Full", 1, 1, free = TRUE, values = start_value, lbound = lbound,
        labels = paste0("s2_", component_name),
        name = paste0("par_", component_name)
      )
    }

    active_kernel_names <- unique(c(
      unname(comp_maps$item),
      if (include_facet_part) unname(comp_maps$facet) else NULL,
      "K_identity",
      if (cov_type == "clustered") "K_clustered" else NULL
    ))
    for (kernel_name in active_kernel_names) {
      add_kernel(kernel_name)
    }

    for (component_name in names(comp_maps$item)) {
      add_parameter(component_name, if (!is.null(sv[[component_name]])) sv[[component_name]] else 0.01)
    }
    if (include_facet_part) {
      for (component_name in names(comp_maps$facet)) {
        add_parameter(component_name, if (!is.null(sv[[component_name]])) sv[[component_name]] else 0.01)
      }
    }
    add_parameter("epsilon", if (!is.null(sv[["epsilon"]])) sv[["epsilon"]] / n_R else 0.01, lbound = 1e-6)
    if (cov_type == "clustered") {
      add_parameter("clustered", if (!is.null(sv[["clustered"]])) sv[["clustered"]] else 0.01)
    }

    item_terms <- paste0("par_", names(comp_maps$item), " %x% ", unname(comp_maps$item))
    if (cov_type == "clustered") {
      item_terms <- c(item_terms, "par_clustered %x% K_clustered")
    }
    item_terms <- c(item_terms, "par_epsilon %x% K_identity")
    item_algebra <- eval(parse(text = sprintf(
      'mxAlgebra(%s, name = "A_cov")',
      paste(item_terms, collapse = " + ")
    )))
    model_parts[[length(model_parts) + 1]] <- item_algebra

    if (include_facet_part) {
      facet_terms <- paste0("par_", names(comp_maps$facet), " %x% ", unname(comp_maps$facet))
      facet_algebra <- eval(parse(text = sprintf(
        'mxAlgebra(%s, name = "B_cov")',
        paste(facet_terms, collapse = " + ")
      )))
      model_parts[[length(model_parts) + 1]] <- facet_algebra
      model_parts[[length(model_parts) + 1]] <- mxMatrix(
        "Full", 1, 1, free = FALSE, values = n_items, name = "n_items_scalar"
      )
      total_algebra <- mxAlgebra(A_cov + n_items_scalar * B_cov, name = "expectedCov")
      model_parts[[length(model_parts) + 1]] <- total_algebra
      mean_mat <- mxMatrix(
        "Full", 1, n_c, free = TRUE,
        values = as.numeric(mean(raw_data)),
        labels = rep("grand_mean", n_c),
        name = "expectedMean",
        dimnames = list(NULL, cond_cols)
      )
    } else {
      total_algebra <- mxAlgebra(A_cov, name = "expectedCov")
      model_parts[[length(model_parts) + 1]] <- total_algebra
      mean_mat <- mxMatrix(
        "Zero", 1, n_c, name = "expectedMean",
        dimnames = list(NULL, cond_cols)
      )
    }
    model_parts[[length(model_parts) + 1]] <- mean_mat

    mxModel(
      group_name,
      model_parts,
      mxExpectationNormal(
        covariance = "expectedCov",
        means = "expectedMean",
        dimnames = cond_cols
      ),
      mxFitFunctionML(),
      mxData(as.data.frame(raw_data), type = "raw")
    )
  }

  contrast_model <- build_group_model("contrast", contrast_data, include_facet_part = FALSE)
  mean_model <- build_group_model("grandmean", mean_data, include_facet_part = TRUE)

  model <- mxModel(
    safe_model_name,
    contrast_model,
    mean_model,
    mxFitFunctionMultigroup(c("contrast", "grandmean"))
  )

  # ---- Fit ----
  cat("Fitting...\n")
  fit <- mxTryHard(model, extraTries = extra_tries, silent = TRUE)

  status <- fit$output$status$code
  if (status != 0) warning(sprintf("Optimizer status: %d", status))
  cat(sprintf("Optimizer status: %d\n", status))

  # ---- Extract parameters ----
  params <- omxGetParameters(fit)

  # Build variance component table on the original variance scale
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

  # ---- G-coefficient ----
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

  # ---- Model fit diagnostics ----
  minus2LL <- fit$output$Minus2LogLikelihood
  n_par <- length(params)
  aic <- minus2LL + 2 * n_par
  bic <- minus2LL + n_par * log(prep$n_p)

  obs_cov <- cov(dat_mat)
  exp_cov <- mxEval(contrast.expectedCov, fit)
  rownames(exp_cov) <- cond_cols
  colnames(exp_cov) <- cond_cols
  resid_cov <- obs_cov - exp_cov
  rmsr <- sqrt(mean(resid_cov[lower.tri(resid_cov, diag = TRUE)]^2))

  # ---- Print ----
  cat(sprintf("\nG-coefficient: %.4f\n", g_coef))
  cat(sprintf("RMSR: %.6f | -2LL: %.2f | AIC: %.2f | nPar: %d | Status: %d\n",
              rmsr, minus2LL, aic, n_par, status))

  # ---- Return ----
  list(
    fit        = fit,
    params     = params,
    vc_list    = vc_list,
    s2_epsilon = vc_list[["epsilon"]],
    g_coef     = g_coef,
    rmsr       = rmsr,
    model_name = safe_model_name,
    cov_type   = cov_type,
    obs_cov    = obs_cov,
    exp_cov    = exp_cov,
    model_fit  = list(
      minus2LL = minus2LL, AIC = aic, BIC = bic,
      n_par = n_par, status = status
    ),
    design = list(
      crossed = crossed, facet_sizes = facet_sizes,
      nested = prep$rep_facets,
      nested_within = prep$nested_within %||% setNames(vector("list", length(prep$rep_facets)), prep$rep_facets),
      facet_roles = prep$facet_roles %||% setNames(rep("facet", length(crossed)), crossed),
      n_p = prep$n_p, n_R = n_R, n_per_group = n_pg, n_rep_groups = n_rg
    )
  )
}


fit_gtheory_openmx_design <- function(
    dat,
    facet_spec,
    item_col = "item_id",
    score_col = "severity",
    n_rep_groups = 1L,
    model_name = "GTheory",
    cov_type = "standard",
    eval_families = NULL,
    item_interaction_orders = 1L,
    facet_interaction_orders = 2L,
    extra_tries = 50L,
    start_values = NULL,
    use_lme4_starts = FALSE
) {
  prep <- prepare_openmx_data_design(
    dat = dat,
    facet_spec = facet_spec,
    item_col = item_col,
    score_col = score_col,
    n_rep_groups = n_rep_groups
  )
  prep$start_values <- start_values

  fit_gtheory_openmx(
    prep = prep,
    model_name = model_name,
    cov_type = cov_type,
    eval_families = eval_families,
    item_interaction_orders = item_interaction_orders,
    facet_interaction_orders = facet_interaction_orders,
    extra_tries = extra_tries,
    original_data = dat,
    item_col = item_col,
    score_col = score_col,
    use_lme4_starts = use_lme4_starts
  )
}


# ============================================================
# Step 4: lme4 helper — build and fit formula dynamically
# ============================================================

fit_gtheory_lme4 <- function(
    dat,
    item_col       = "item_id",
    score_col      = "severity",
    crossed_facets = c("evaluator", "prompt_type"),
    clustered      = NULL,  # list(facet = "evaluator", family_col = "eval_family")
    item_interaction_orders = 1L,
    facet_interaction_orders = 2L,
    use_REML       = TRUE
) {
  suppressPackageStartupMessages(library(lme4))
  suppressPackageStartupMessages(library(dplyr))

  df <- dat
  # Convert to factors
  df[[item_col]] <- factor(df[[item_col]])
  for (f in crossed_facets) {
    df[[f]] <- factor(df[[f]])
  }

  # Build formula: severity ~ 1 + (1|item) + (1|facet1) + ... + (1|item:facet1) + ... + (1|facet1:facet2) + ...
  re_terms <- c()
  n_crossed <- length(crossed_facets)

  # Main effects
  re_terms <- c(re_terms, paste0("(1 | ", item_col, ")"))
  for (f in crossed_facets) {
    re_terms <- c(re_terms, paste0("(1 | ", f, ")"))
  }

  # Item × facet-subset interactions
  cell_counts <- df %>%
    count(across(all_of(c(item_col, crossed_facets))), name = "n")
  has_within_cell_replication <- any(cell_counts$n > 1)
  max_estimable_item_order <- if (has_within_cell_replication) n_crossed else max(n_crossed - 1, 1)
  item_orders <- sort(unique(as.integer(item_interaction_orders)))
  item_orders <- item_orders[item_orders >= 1 & item_orders <= max_estimable_item_order]

  if (length(item_orders) > 0) {
    for (subset in enumerate_subsets(crossed_facets, min(item_orders), max(item_orders))) {
      if (length(subset) %in% item_orders) {
        re_terms <- c(re_terms, paste0("(1 | ", paste(c(item_col, subset), collapse = ":"), ")"))
      }
    }
  }

  # Facet × facet interactions
  facet_orders <- sort(unique(as.integer(facet_interaction_orders)))
  facet_orders <- facet_orders[facet_orders >= 2 & facet_orders <= n_crossed]

  if (length(facet_orders) > 0) {
    for (subset in enumerate_subsets(crossed_facets, min(facet_orders), max(facet_orders))) {
      if (length(subset) %in% facet_orders) {
        re_terms <- c(re_terms, paste0("(1 | ", paste(subset, collapse = ":"), ")"))
      }
    }
  }

  # Clustered: item × eval_family
  if (!is.null(clustered)) {
    fam_col <- clustered$family_col
    df[[fam_col]] <- factor(df[[fam_col]])
    re_terms <- c(re_terms, paste0("(1 | ", item_col, ":", fam_col, ")"))
  }

  formula_str <- paste(score_col, "~ 1 +", paste(re_terms, collapse = " + "))
  formula_obj <- as.formula(formula_str)

  cat(sprintf("lme4 formula: %s\n", formula_str))

  fit <- lmer(
    formula_obj, data = df,
    REML = use_REML,
    control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000))
  )

  # Extract variance components
  vc <- as.data.frame(VarCorr(fit))
  get_vc <- function(grp) {
    val <- vc$vcov[vc$grp == grp]
    if (length(val) == 0) 0 else as.numeric(val)
  }

  vc_results <- list()
  vc_results$s2_tau <- get_vc(item_col)
  for (f in crossed_facets) {
    vc_results[[paste0("s2_", f)]] <- get_vc(f)
  }

  if (length(item_orders) > 0) {
    for (subset in enumerate_subsets(crossed_facets, min(item_orders), max(item_orders))) {
      if (length(subset) %in% item_orders) {
        key <- paste0("s2_tau_", paste(subset, collapse = "_"))
        grp <- paste(c(item_col, subset), collapse = ":")
        vc_results[[key]] <- get_vc(grp)
      }
    }
  }

  if (length(facet_orders) > 0) {
    for (subset in enumerate_subsets(crossed_facets, min(facet_orders), max(facet_orders))) {
      if (length(subset) %in% facet_orders) {
        key <- paste0("s2_", paste(subset, collapse = "_"))
        grp <- paste(subset, collapse = ":")
        vc_results[[key]] <- get_vc(grp)
      }
    }
  }
  vc_results$s2_epsilon <- sigma(fit)^2

  if (!is.null(clustered)) {
    vc_results$s2_clustered <- get_vc(paste0(item_col, ":", clustered$family_col))
  }

  list(fit = fit, vc = vc_results)
}


# ============================================================
# Step 5: Model comparison
# ============================================================

compare_openmx_models <- function(results_list) {
  cat("=== Model Comparison ===\n")
  cat(sprintf("%-25s %12s %12s %12s %6s %8s\n",
              "Model", "-2LL", "AIC", "BIC", "nPar", "RMSR"))
  cat(strrep("-", 80), "\n")

  for (res in results_list) {
    mf <- res$model_fit
    cat(sprintf("%-25s %12.2f %12.2f %12.2f %6d %8.4f\n",
                res$cov_type, mf$minus2LL, mf$AIC, mf$BIC,
                mf$n_par, res$rmsr))
  }

  if (length(results_list) >= 2) {
    cat("\nLikelihood Ratio Tests (each vs first):\n")
    base <- results_list[[1]]
    for (i in 2:length(results_list)) {
      ext <- results_list[[i]]
      chi2 <- base$model_fit$minus2LL - ext$model_fit$minus2LL
      df_diff <- ext$model_fit$n_par - base$model_fit$n_par
      if (df_diff > 0 && chi2 > 0) {
        p_val <- pchisq(chi2, df = df_diff, lower.tail = FALSE)
        cat(sprintf("  %s vs %s: chi2=%.2f, df=%d, p=%.4f %s\n",
                    ext$cov_type, base$cov_type,
                    chi2, df_diff, p_val,
                    if (p_val < 0.05) "***" else "ns"))
      }
    }
  }
}
