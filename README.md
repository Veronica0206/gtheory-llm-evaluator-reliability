# Generalizability Theory for LLM-as-Evaluator Reliability

Code and results accompanying the revised manuscript by **Jin Liu and Robert A. Perera**:

> *Generalizability Theory for LLM-as-Evaluator Reliability: Univariate and Multivariate Variance Decomposition Across Models, Prompts, and Temperatures in Text Classification.*

The study examines annotation reliability across evaluators, prompts, temperatures and seeds in Hate-Speech, Mental-Health and Drug-Review tasks. The revised primary analysis uses OpenMx for Gaussian ML/REML estimation, including jointly estimated cross-trait covariance components.

## Start here

**Download this repository and open [index.html](index.html) in a browser.** GitHub’s file view displays the HTML source. The reports contain the analysis code and saved outputs; no standalone Rmd, R, Python or notebook source files are distributed in this revision.

| Section | Contents |
|---|---|
| [Labelling](labelling/index.html) | Pointers to the supplementary prompts and archived annotations; no labelling notebooks |
| [Primary analysis: OpenMx](primary_analysis/index.html) | Five reports covering seven outcome sets, model comparisons, variance/covariance estimates, reliability and D-study intervals |
| [Supportive analysis: lme4 and MoM](supportive_analysis/index.html) | Five completed lme4 point-estimate reports and univariate method-of-moments diagnostics |
| [Other analysis](other_analysis/index.html) | Saved agreement, entropy and calibration analyses; archived facet-severity intervals and cross-trait diagnostics |

## Paired analysis reports

| Analysis | Primary OpenMx | Supportive lme4 |
|---|---|---|
| Hate-Speech, univariate | [HTML](primary_analysis/uni_hate_speech.html) | [HTML](supportive_analysis/lme4/uni_hate_speech.html) |
| Mental-Health, seven- and three-level univariate | [HTML](primary_analysis/uni_mental_health.html) | [HTML](supportive_analysis/lme4/uni_mental_health.html) |
| Drug-Review, univariate | [HTML](primary_analysis/uni_drug_review.html) | [HTML](supportive_analysis/lme4/uni_drug_review.html) |
| Mental-Health, multivariate | [HTML](primary_analysis/mv_mental_health.html) | [HTML](supportive_analysis/lme4/mv_mental_health.html) |
| Drug-Review, multivariate | [HTML](primary_analysis/mv_drug_review.html) | [HTML](supportive_analysis/lme4/mv_drug_review.html) |

The OpenMx reports are the revised primary results. All seven outcome sets reached **1,000 accepted bootstrap fits and 100 completed jackknife fits**. Attempt counts, unavailable estimates and sparse-tail diagnostics remain visible. The final report refresh reused the completed fits. The two multivariate reports now also provide the trait-specific variance percentages in manuscript Table 5, their calculation, and full-precision CSV downloads.

The supplementary lme4 reports were rerun on 10 September 2026. They compare all **11 specifications** by ML and select the within-Family 2 BIC winner for REML point estimates; all seven outcome sets selected `f2_complete_drift`. They report components, reliability and D-study projections without bootstrap, jackknife or confidence intervals.

The multivariate reports use one stacked model with trait-specific diagonal source variances and a shared residual variance, following the covariance restrictions of the original submission. They cover Mental-Health’s six binary flags and separate three OR groups, and Drug-Review’s four aspects. The Mental-Health groups are stress; anxiety or depression; and bipolar disorder, personality disorder or suicidal ideation. Groups can overlap. Cross-trait covariances are fixed at zero, so comparison with the primary full-covariance OpenMx analysis is a covariance-restriction sensitivity.

The manuscript’s retained facet-severity intervals and cross-trait diagnostics are published separately in [Other analysis](other_analysis/index.html). These focused historical reports preserve their original code and saved outputs; they are distinct from the new point-only lme4 fits.

## Data and reproduction

- [OSF data deposit](https://doi.org/10.17605/OSF.IO/K9CAJ).
- [Data checksums, source extraction, environments and rerun instructions](reproducing.html).
- Each OpenMx or lme4 G-theory HTML includes its complete authoring source through **Code → Download Rmd**, as well as code, numerical outputs and session information. The source can be recovered locally from the HTML; no separate Rmd upload is needed.
- The five lme4 HTMLs additionally contain **106 downloadable output files**, with exact saved CSV and RDS bytes. Use each report’s “Download the complete point-estimate outputs” section.
- The agreement reports preserve the saved notebook outputs and identify historical wording/provenance limitations. Labelling notebooks are omitted; the Supplementary Materials contain the prompts.

The annotation CSVs on OSF match the primary report fingerprints. Supportive and ancillary reruns also need three sampled-item CSVs, which are not in the current public deposit; their filenames and checksums are listed in the reproduction guide.

Input data and optional fit caches are not duplicated here. This GitHub update does not update the OSF deposit. The repository includes the completed lme4 rerun; packaging preserves those finished results without additional G-theory fitting, bootstrap resampling or annotation. Focused descriptive source checks verify retained outputs. The two multivariate OpenMx reports have received a reporting-only refresh; existing model estimates and bootstrap intervals are preserved.

## Earlier public version

The previous public R implementation is preserved at the [pre-openmx-r1 tag](https://github.com/Veronica0206/gtheory-llm-evaluator-reliability/tree/pre-openmx-r1), based on commit `97dfc8d57999b0f16fc9a5c579d6f4eee037d2ba`. It is distinct from the current supplementary lme4 HTML reports. Existing Git history is retained.

## Citation and license

Please cite the manuscript and identify the repository commit used in your work. The associated data deposit DOI is [10.17605/OSF.IO/K9CAJ](https://doi.org/10.17605/OSF.IO/K9CAJ).

Code is released under the [MIT License](LICENSE). Source datasets retain their respective licensing and attribution requirements; consult the OSF deposit and original corpus sources. Questions and reproducibility reports are welcome through GitHub issues.
