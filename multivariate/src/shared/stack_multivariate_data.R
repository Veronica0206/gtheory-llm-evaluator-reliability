source("univariate/openmx_matrix/fit_gtheory_openmx.R")

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(tibble))


# ============================================================
# Multivariate outcome stacking helpers
#
# Mental health:
#   - six binary comorbidity flags
#   - returns long data with y_bin / y_num in {0,1}
#
# Drug review:
#   - four ordinal aspect ratings
#   - returns long data with ordered labels and y_num in {1,2,3}
# ============================================================

MH_TRAITS <- c(
  depression_present = "depression",
  anxiety_present = "anxiety",
  suicidal_present = "suicidal",
  stress_present = "stress",
  bipolar_present = "bipolar",
  personality_disorder_present = "personality_disorder"
)

DRUG_TRAITS <- c(
  ai_efficacy = "efficacy",
  ai_safety = "safety",
  ai_burden = "burden",
  ai_cost = "cost"
)

DRUG_LEVELS <- c("NEGATIVE", "NEUTRAL", "POSITIVE")
DRUG_SCORE_MAP <- setNames(seq_along(DRUG_LEVELS), DRUG_LEVELS)


coerce_tf_to_binary <- function(x) {
  if (is.logical(x)) {
    return(as.integer(x))
  }

  if (is.numeric(x)) {
    return(as.integer(x))
  }

  x_chr <- as.character(x)
  out <- ifelse(x_chr %in% c("True", "TRUE", "1"), 1L,
                ifelse(x_chr %in% c("False", "FALSE", "0"), 0L, NA_integer_))
  as.integer(out)
}


stack_mental_health_multivariate <- function(dat) {
  needed <- names(MH_TRAITS)
  missing <- setdiff(needed, names(dat))
  if (length(missing) > 0) {
    stop("Missing mental-health trait columns: ", paste(missing, collapse = ", "))
  }

  dat %>%
    mutate(across(all_of(needed), coerce_tf_to_binary)) %>%
    pivot_longer(
      cols = all_of(needed),
      names_to = "trait_col",
      values_to = "y_bin"
    ) %>%
    mutate(
      task = "mental_health",
      trait = factor(unname(MH_TRAITS[trait_col]), levels = unname(MH_TRAITS)),
      y_bin = as.integer(y_bin),
      y_num = as.numeric(y_bin),
      outcome_type = "binary"
    ) %>%
    select(
      task, item_id, evaluator, prompt_type, temperature, seed,
      provider, model_string, trait, outcome_type, y_bin, y_num, everything()
    )
}


stack_drug_review_multivariate <- function(dat) {
  needed <- names(DRUG_TRAITS)
  missing <- setdiff(needed, names(dat))
  if (length(missing) > 0) {
    stop("Missing drug-review trait columns: ", paste(missing, collapse = ", "))
  }

  dat %>%
    mutate(across(all_of(needed), as.character)) %>%
    pivot_longer(
      cols = all_of(needed),
      names_to = "trait_col",
      values_to = "y_label"
    ) %>%
    mutate(
      task = "drug_review",
      trait = factor(unname(DRUG_TRAITS[trait_col]), levels = unname(DRUG_TRAITS)),
      y_label = factor(y_label, levels = DRUG_LEVELS, ordered = TRUE),
      y_num = as.numeric(unname(DRUG_SCORE_MAP[as.character(y_label)])),
      outcome_type = "ordinal_3"
    ) %>%
    select(
      task, item_id, evaluator, prompt_type, temperature, seed,
      provider, model_string, trait, outcome_type, y_label, y_num, everything()
    )
}


read_and_stack_multivariate_data <- function(task = c("mental_health", "drug_review")) {
  task <- match.arg(task)

  if (task == "mental_health") {
    dat <- read.csv("data/mh_labeling_final.csv", stringsAsFactors = FALSE)
    return(stack_mental_health_multivariate(dat))
  }

  dat <- read.csv("data/drug_labeling_final.csv", stringsAsFactors = FALSE)
  stack_drug_review_multivariate(dat)
}


if (sys.nframe() == 0) {
  mh_long <- read_and_stack_multivariate_data("mental_health")
  drug_long <- read_and_stack_multivariate_data("drug_review")

  dir.create("outputs/data_model_fits", showWarnings = FALSE, recursive = TRUE)
  write.csv(mh_long, "outputs/data_model_fits/mh_multivariate_long.csv", row.names = FALSE)
  write.csv(drug_long, "outputs/data_model_fits/drug_multivariate_long.csv", row.names = FALSE)

  cat("Saved: outputs/data_model_fits/mh_multivariate_long.csv\n")
  cat("Saved: outputs/data_model_fits/drug_multivariate_long.csv\n")
}
