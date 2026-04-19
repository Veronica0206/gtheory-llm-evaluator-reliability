omx_resolve_num_threads <- function(num_threads = NULL) {
  if (!is.null(num_threads)) {
    resolved <- suppressWarnings(as.integer(num_threads[[1]]))
    if (is.finite(resolved) && !is.na(resolved) && resolved >= 1L) {
      return(resolved)
    }
  }

  env_candidates <- c(
    Sys.getenv("OPENMX_NUM_THREADS", unset = ""),
    Sys.getenv("OMP_NUM_THREADS", unset = "")
  )
  for (candidate in env_candidates) {
    resolved <- suppressWarnings(as.integer(candidate))
    if (is.finite(resolved) && !is.na(resolved) && resolved >= 1L) {
      return(resolved)
    }
  }

  resolved <- suppressWarnings(parallel::detectCores(logical = FALSE))
  if (!is.finite(resolved) || is.na(resolved) || resolved < 1L) {
    resolved <- suppressWarnings(parallel::detectCores())
  }
  if (!is.finite(resolved) || is.na(resolved) || resolved < 1L) {
    resolved <- 1L
  }

  as.integer(resolved)
}


omx_configure_threads <- function(num_threads = NULL,
                                  parallel_diagnostics = FALSE,
                                  log_fn = NULL) {
  resolved_threads <- omx_resolve_num_threads(num_threads)

  Sys.setenv(
    OPENMX_NUM_THREADS = resolved_threads,
    OMP_NUM_THREADS = resolved_threads
  )
  OpenMx::mxOption(model = NULL, key = "Number of Threads", value = resolved_threads)
  OpenMx::mxOption(
    model = NULL,
    key = "Parallel diagnostics",
    value = if (isTRUE(parallel_diagnostics)) "Yes" else "No"
  )

  if (is.function(log_fn)) {
    log_fn(sprintf("Configured OpenMx threads: %d\n", resolved_threads))
  }

  invisible(resolved_threads)
}
