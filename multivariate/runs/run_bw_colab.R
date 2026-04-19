parse_flag <- function(args, name, default = NULL) {
  prefix <- paste0("--", name, "=")
  match <- grep(paste0("^", prefix), args, value = TRUE)
  if (length(match) == 0) {
    return(default)
  }
  sub(prefix, "", match[[1]])
}


parse_nullable_int <- function(x, default = NULL) {
  if (is.null(x) || identical(x, "") || identical(toupper(x), "NULL")) {
    return(default)
  }
  as.integer(x)
}


parse_flag_bool <- function(args, name, default = FALSE) {
  raw <- parse_flag(args, name, if (isTRUE(default)) "TRUE" else "FALSE")
  identical(toupper(raw), "TRUE")
}


default_data_path <- function(task) {
  switch(
    task,
    mental_health = "data/mh_labeling_final.csv",
    mental_health_sensitivity = "data/mh_labeling_final.csv",
    drug_review = "data/drug_labeling_final.csv",
    stop("Unsupported multivariate task: ", task)
  )
}


args <- commandArgs(trailingOnly = TRUE)

task <- parse_flag(args, "task", "mental_health")
data_path <- parse_flag(args, "data-path", default_data_path(task))
max_items <- parse_nullable_int(parse_flag(args, "max-items", "NULL"))
extra_tries <- parse_nullable_int(parse_flag(args, "extra-tries", "5"), 5L)
output_dir <- parse_flag(args, "output-dir", "multivariate/outputs/bw")
output_name <- parse_flag(args, "output-name", task)
num_threads <- parse_nullable_int(parse_flag(args, "num-threads", "NULL"))
parallel_diagnostics <- parse_flag_bool(args, "parallel-diagnostics", FALSE)

if (!file.exists("multivariate/src/bw_models.R")) {
  stop(
    "Expected to run from the project root (the directory that contains ",
    "'multivariate/src/bw_models.R')."
  )
}

source("univariate/src/shared/openmx_threads.R")
source("multivariate/src/bw_models.R")

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

resolved_threads <- omx_configure_threads(
  num_threads = num_threads,
  parallel_diagnostics = parallel_diagnostics,
  log_fn = mvbw_log
)

cat("Multivariate BW Colab runner\n")
cat("  task        :", task, "\n")
cat("  data_path   :", data_path, "\n")
cat("  max_items   :", if (is.null(max_items)) "NULL" else max_items, "\n")
cat("  extra_tries :", extra_tries, "\n")
cat("  output_dir  :", output_dir, "\n")
cat("  output_name :", output_name, "\n")
cat("  num_threads :", resolved_threads, "\n")
cat("  parallel_diagnostics :", parallel_diagnostics, "\n\n")

if (identical(task, "mental_health_sensitivity")) {
  dat_bundle <- mvbw_prepare_data(
    task = "mental_health",
    data_path = data_path,
    trait_subset = MG_MH_SENSITIVITY_TRAITS
  )
} else {
  dat_bundle <- mvbw_prepare_data(task = task, data_path = data_path)
}

cat(sprintf("Dataset prepared | Items: %d | Rows: %d | Traits: %d (%s)\n",
            nlevels(dat_bundle$data$item_id), nrow(dat_bundle$data),
            dat_bundle$p, paste(dat_bundle$traits, collapse = ", ")))

results <- mvbw_run_all_models(
  dat_bundle,
  extra_tries = extra_tries,
  max_items = max_items
)

out_file <- file.path(output_dir, paste0(output_name, "_results.RData"))
save(results, file = out_file)

cat("\nSaved results to", out_file, "\n")
