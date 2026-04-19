# Generalizability Theory for LLM-as-Evaluator Reliability

Analysis code accompanying:

> Liu, J. (2026). *Generalizability Theory for LLM-as-Evaluator Reliability: Univariate and Multivariate Variance Decomposition Across Models, Prompts, and Temperatures in Text Classification*. Behavior Research Methods.

## Companion artifacts

| Resource | Location |
|---|---|
| Derived data and analysis outputs | OSF: https://doi.org/10.17605/OSF.IO/K9CAJ |
| Code (this repository) | GitHub: https://github.com/Veronica0206/gtheory-llm-evaluator-reliability |
| Code (immutable archive) | Zenodo: TBD |

## Repository scope

This repository contains the R analysis pipeline for the **lme4-based univariate and multivariate generalizability theory** analyses reported in the manuscript.

```
.
├── univariate/
│   ├── compare_common_methods.R        # shared utilities and dataset specs
│   ├── runs/                           # per-task driver scripts
│   │   └── run_lme4_*.R                # lme4 fits for hate speech, mental
│   │                                   # health (7L and 3L), and drug review
│   └── src/
│       ├── lme4_models.R               # univariate G-theory via lme4
│       └── shared/
│           ├── gtheory_utils.R         # shared computation helpers
│           └── univariate_workflow.R   # design specs and data preparation
├── multivariate/
│   ├── runs/
│   │   └── run_lme4_*.R                # multivariate drivers (mental health,
│   │                                   # mental health 3-condition sensitivity,
│   │                                   # drug review)
│   └── src/
│       ├── lme4_models.R               # multivariate stacked lmer
│       └── shared/
│           └── multivariate_workflow.R # multivariate design and stacking
├── generate_disattenuated_figures.R    # Figure 1 of manuscript
├── README.md
├── LICENSE                             # MIT
└── .gitignore
```

## Reproducing the analyses

### 1. Get the data

Download the LLM annotation CSVs from the OSF deposit:

> https://osf.io/k9caj/ → "Datasets with LLM Annotations" folder

Place all three files in a `data/` folder at the repository root:

```
gtheory-llm-evaluator-reliability/
└── data/
    ├── hate_labeling_final.csv
    ├── mh_labeling_final.csv
    └── drug_labeling_final.csv
```

(`data/` is gitignored; CSVs are not committed.)

### 2. Install R packages

R 4.5+ recommended. Required packages:

```r
install.packages(c(
  "lme4",        # mixed-effects models
  "dplyr",       # data wrangling
  "tibble",      # tibbles
  "numDeriv",    # Hessian/Jacobian for analytic Wald CIs
  "boot"         # bootstrap (for univariate CIs)
))
```

### 3. Run univariate analyses (manuscript Tables 1–5)

From the repository root, source any of the lme4 driver scripts. For example:

```r
# Hate-Speech univariate analysis
source("univariate/runs/run_lme4_hate_speech.R")
```

Each driver fits the three nested model specifications (full crossed, nested-seed, nested-seed-temperature) and computes variance components, $E\rho^2$, $\Phi$, bootstrap and Wald CIs, and D-study projections.

Outputs are written to `univariate/outputs/lme4/<task>_results.RData`. The `outputs/` folder is gitignored; the final outputs used in the manuscript are archived on OSF.

### 4. Run multivariate analyses (manuscript Tables 6–7)

```r
# Mental-Health 6-condition multivariate
source("multivariate/runs/run_lme4_mental_health.R")

# Drug-Review 4-aspect multivariate
source("multivariate/runs/run_lme4_drug_review.R")

# Mental-Health 3-condition sensitivity
source("multivariate/runs/run_lme4_mental_health_sensitivity.R")
```

### 5. Reproduce Figure 1

```r
source("generate_disattenuated_figures.R")
```

Outputs PNG figures and the underlying CSV correlation matrices for the Mental-Health (6-flag) and Drug-Review (4-aspect) disattenuated universe-score correlations.

## Citation

If you use this code, please cite both the manuscript and this repository:

```bibtex
@article{liu2026gtheory,
  author  = {Liu, Jin},
  title   = {Generalizability Theory for {LLM-as-Evaluator} Reliability:
             Univariate and Multivariate Variance Decomposition Across Models,
             Prompts, and Temperatures in Text Classification},
  journal = {Behavior Research Methods},
  year    = {2026}
}

@software{liu2026gtheory_code,
  author  = {Liu, Jin},
  title   = {Generalizability Theory for LLM-as-Evaluator Reliability:
             Analysis Code},
  year    = {2026},
  url     = {https://github.com/Veronica0206/gtheory-llm-evaluator-reliability},
  doi     = {TBD (Zenodo DOI after release tagged)}
}
```

## License

Code released under the [MIT License](LICENSE). Data on OSF released under CC-BY-4.0.

## Contact

For questions about the code or analyses, please open a GitHub issue or contact the author at Veronica.Liu0206@gmail.com.
