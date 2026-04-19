# ============================================================
# G-Theory Utilities: Shared functions for G-coefficient
# computation and configuration defaults.
#
# Sourced by sim_layer*.R, validate_single_run.R, and
# any script needing lme4 G-coefficient formulas.
# ============================================================


#' Compute the relative-error G-coefficient from lme4 variance components.
#'
#' Works for any combination of crossed facets, with optional clustered bias.
#' Facet sizes are extracted from the facet_spec or supplied directly.
#'
#' @param vc Named list of variance components from fit_gtheory_lme4()$vc.
#' @param facet_spec Named list with elements like list(n = 4, type = "crossed").
#'   Only crossed facets contribute item-interaction terms to the denominator.
#' @param n_rep_per_cell Total number of replications per (item x crossed-facet) cell.
#'   Defaults to 1 if no replication facets are present.
#' @param n_families Number of evaluator families (for clustered models). NULL if none.
#' @return Numeric G-coefficient.
compute_lme4_g <- function(vc, facet_spec, n_rep_per_cell = 1, n_families = NULL) {
  s2_tau <- vc$s2_tau
  if (is.null(s2_tau) || is.na(s2_tau)) return(NA_real_)

  crossed_names <- names(facet_spec)[vapply(facet_spec, function(x) x$type == "crossed", logical(1))]
  crossed_sizes <- vapply(crossed_names, function(f) facet_spec[[f]]$n, numeric(1))

  # Item x single-facet interaction terms: s2_tau_<facet> / n_<facet>
  sigma2_delta <- 0
  for (f in crossed_names) {
    vc_name <- paste0("s2_tau_", f)
    if (vc_name %in% names(vc) && !is.na(vc[[vc_name]])) {
      sigma2_delta <- sigma2_delta + vc[[vc_name]] / crossed_sizes[[f]]
    }
  }

  # Clustered contribution
  if (!is.null(n_families) && "s2_clustered" %in% names(vc) && !is.na(vc$s2_clustered)) {
    sigma2_delta <- sigma2_delta + vc$s2_clustered / n_families
  }

  # Residual contribution: s2_epsilon / (product of crossed sizes * n_rep_per_cell)
  n_total <- prod(crossed_sizes) * n_rep_per_cell
  sigma2_delta <- sigma2_delta + vc$s2_epsilon / n_total

  s2_tau / (s2_tau + sigma2_delta)
}


#' Default simulation configuration with validation.
#'
#' Provides default values for sim_relative_bias.R configuration,
#' avoiding fragile exists() checks.
#'
#' @param overrides Named list of values to override defaults.
#' @return Named list of configuration values.
sim_bias_config <- function(overrides = list()) {
  defaults <- list(
    BIAS_N_REPS = 100,
    BIAS_N_ITEMS = 100,
    BIAS_USE_TRUE_STARTS = TRUE,
    BIAS_L3_CLUSTERED_S2 = 0.05,
    BIAS_RUN_LAYERS = c("L1", "L2", "L3")
  )

  for (nm in names(overrides)) {
    if (nm %in% names(defaults)) {
      defaults[[nm]] <- overrides[[nm]]
    } else {
      warning(sprintf("Unknown config key '%s' ignored.", nm))
    }
  }

  defaults
}
