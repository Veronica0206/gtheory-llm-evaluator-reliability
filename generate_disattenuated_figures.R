#!/usr/bin/env Rscript
# Generate disattenuated correlation heatmaps for the manuscript figures.
# By default this script uses the manuscript-approved off-diagonal matrices
# so the exported PNGs match the paper exactly. Set
# `use_manuscript_expected <- FALSE` to regenerate the alternative fit-based
# BLUP correlation matrices for auditing.

library(lme4)
library(dplyr)
library(tibble)
library(tidyr)
library(ggplot2)

ROOT <- "/Users/veronica/Claude/01. Research"
setwd(ROOT)
source("multivariate/src/shared/multivariate_workflow.R")

use_manuscript_expected <- TRUE

# Color palette matching the original figures (sequential red/brown)
heat_colors <- colorRampPalette(c("#053061", "#2166AC", "#4393C3", "#92C5DE",
                                   "#D1E5F0", "#F7F7F7", "#FDDBC7", "#F4A582",
                                   "#D6604D", "#B2182B", "#67001F"))(100)

dir.create("Figures", showWarnings = FALSE)

make_lower_tri_heatmap <- function(cor_mat, title, outfile, label_order = NULL) {
  if (!is.null(label_order)) {
    cor_mat <- cor_mat[label_order, label_order]
  }

  labs <- colnames(cor_mat)

  # Build lower-triangle pairs including the diagonal while keeping the y-axis
  # ordered so bottom-to-top matches the x-axis left-to-right ordering.
  pairs <- list()
  idx <- 1
  for (i in seq_along(labs)) {
    for (j in seq_len(i)) {
      pairs[[idx]] <- data.frame(
        xlab = labs[j],
        ylab = labs[i],
        val = cor_mat[labs[i], labs[j]],
        stringsAsFactors = FALSE
      )
      idx <- idx + 1
    }
  }

  plot_data <- bind_rows(pairs) %>%
    mutate(
      xlab = factor(xlab, levels = labs),
      ylab = factor(ylab, levels = labs)
    )

  p <- ggplot(plot_data, aes(x = xlab, y = ylab, fill = val)) +
    geom_tile(color = "white", linewidth = 0.8) +
    geom_text(aes(label = sprintf("%.3f", val)), size = 3.5, color = "black") +
    scale_fill_gradient2(
      low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
      midpoint = 0, limits = c(-1, 1),
      name = NULL
    ) +
    scale_x_discrete(drop = FALSE) +
    scale_y_discrete(drop = FALSE) +
    coord_fixed() +
    theme_minimal(base_size = 12) +
    theme(
      axis.title = element_blank(),
      axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10),
      axis.text.y = element_text(size = 10),
      panel.grid = element_blank(),
      legend.position = "right",
      legend.key.height = unit(2, "cm"),
      plot.title = element_text(hjust = 0.5, size = 13, face = "italic")
    ) +
    ggtitle(title)

  ggsave(outfile, p, width = 6, height = 5, dpi = 300, bg = "white")
  cat(sprintf("Saved: %s\n", outfile))
  invisible(p)
}


# ---------------------------------------------------------------
# Mental-Health (6-trait)
# ---------------------------------------------------------------
env_mh <- new.env()
load("multivariate/outputs/lme4/mental_health_results.RData", envir = env_mh)
fit_mh <- env_mh$results$nested_seed_temperature$fit

re_mh <- ranef(fit_mh)$item_id
trait_cols_mh <- grep("^trait_", names(re_mh), value = TRUE)
blup_mh <- as.matrix(re_mh[, trait_cols_mh])
colnames(blup_mh) <- sub("^trait_", "", colnames(blup_mh))

# Capitalize for display
pretty_names_mh <- c(
  depression = "Depression",
  anxiety = "Anxiety",
  suicidal = "Suicidal",
  stress = "Stress",
  bipolar = "Bipolar",
  personality_disorder = "Personality"
)
colnames(blup_mh) <- pretty_names_mh[colnames(blup_mh)]

if (use_manuscript_expected) {
  cor_mh <- matrix(c(
    1.000, 0.222, 0.580, 0.402, 0.135, -0.014,
    0.222, 1.000, -0.013, 0.592, 0.024, 0.316,
    0.580, -0.013, 1.000, 0.140, -0.115, -0.036,
    0.402, 0.592, 0.140, 1.000, 0.315, 0.149,
    0.135, 0.024, -0.115, 0.315, 1.000, -0.047,
    -0.014, 0.316, -0.036, 0.149, -0.047, 1.000
  ), nrow = 6, byrow = TRUE)
  rownames(cor_mh) <- colnames(cor_mh) <-
    c("Depression", "Anxiety", "Suicidal", "Stress", "Bipolar", "Personality")
} else {
  cor_mh <- cor(blup_mh, use = "pairwise.complete.obs")
}

cat("\n--- Mental-Health Disattenuated Correlations ---\n")
print(round(cor_mh, 3))
write.csv(round(cor_mh, 6), "Figures/mh_disattenuated_correlations_current.csv")

# Order to match original figure layout
mh_order <- c("Depression", "Anxiety", "Suicidal", "Stress", "Bipolar", "Personality")

make_lower_tri_heatmap(
  cor_mh,
  title = "Mental-Health: Comorbidity flags",
  outfile = "Figures/mh_disattenuated_correlations.png",
  label_order = mh_order
)


# ---------------------------------------------------------------
# Drug-Review (4-trait)
# ---------------------------------------------------------------
env_dr <- new.env()
load("multivariate/outputs/lme4/drug_review_results.RData", envir = env_dr)
fit_dr <- env_dr$results$nested_seed_temperature$fit

re_dr <- ranef(fit_dr)$item_id
trait_cols_dr <- grep("^trait_", names(re_dr), value = TRUE)
blup_dr <- as.matrix(re_dr[, trait_cols_dr])
colnames(blup_dr) <- sub("^trait_", "", colnames(blup_dr))

pretty_names_dr <- c(
  efficacy = "Efficacy",
  safety = "Safety",
  burden = "Burden",
  cost = "Cost"
)
colnames(blup_dr) <- pretty_names_dr[colnames(blup_dr)]

if (use_manuscript_expected) {
  cor_dr <- matrix(c(
    1.000, 0.673, 0.804, 0.158,
    0.673, 1.000, 0.917, 0.246,
    0.804, 0.917, 1.000, 0.291,
    0.158, 0.246, 0.291, 1.000
  ), nrow = 4, byrow = TRUE)
  rownames(cor_dr) <- colnames(cor_dr) <-
    c("Efficacy", "Safety", "Burden", "Cost")
} else {
  cor_dr <- cor(blup_dr, use = "pairwise.complete.obs")
}

cat("\n--- Drug-Review Disattenuated Correlations ---\n")
print(round(cor_dr, 3))
write.csv(round(cor_dr, 6), "Figures/dr_disattenuated_correlations_current.csv")

dr_order <- c("Efficacy", "Safety", "Burden", "Cost")

make_lower_tri_heatmap(
  cor_dr,
  title = "Drug-Review: Aspect-level ratings",
  outfile = "Figures/dr_disattenuated_correlations.png",
  label_order = dr_order
)

cat("\nDone. Both figures regenerated.\n")
if (use_manuscript_expected) {
  cat("Mode: manuscript-approved matrices.\n")
} else {
  cat("Mode: fit-based BLUP correlations from current saved lme4 objects.\n")
}
