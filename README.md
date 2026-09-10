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
| [Supportive analysis: lme4 and MoM](supportive_analysis/index.html) | Five earlier R1 lme4 reports and an index of their method-of-moments diagnostics |
| [Other analysis](other_analysis/index.html) | Saved agreement, entropy and calibration analyses; links to facet-severity and descriptive diagnostics |

## Paired analysis reports

| Analysis | Primary OpenMx | Supportive lme4 |
|---|---|---|
| Hate-Speech, univariate | [HTML](primary_analysis/uni_hate_speech.html) | [HTML](supportive_analysis/lme4/uni_hate_speech.html) |
| Mental-Health, seven- and three-level univariate | [HTML](primary_analysis/uni_mental_health.html) | [HTML](supportive_analysis/lme4/uni_mental_health.html) |
| Drug-Review, univariate | [HTML](primary_analysis/uni_drug_review.html) | [HTML](supportive_analysis/lme4/uni_drug_review.html) |
| Mental-Health, multivariate | [HTML](primary_analysis/mv_mental_health.html) | [HTML](supportive_analysis/lme4/mv_mental_health.html) |
| Drug-Review, multivariate | [HTML](primary_analysis/mv_drug_review.html) | [HTML](supportive_analysis/lme4/mv_drug_review.html) |

The OpenMx reports are the revised primary results. All seven outcome sets reached **1,000 accepted bootstrap fits and 100 completed jackknife fits**. Attempt counts, unavailable estimates and sparse-tail diagnostics remain visible. The final report refresh reused the completed fits.

The earlier multivariate lme4 reports fit traits separately and use diagonal composites. Their Mental-Health three-condition subset (depression, anxiety and suicidal ideation) differs from OpenMx’s separate three OR-group model. They provide supportive reference analyses and retained diagnostics, not a comparison differing only in software.

## Data and reproduction

- [OSF data deposit](https://doi.org/10.17605/OSF.IO/K9CAJ).
- [Data checksums, source extraction, environments and rerun instructions](reproducing.html).
- Each R HTML includes its complete authoring source through **Code → Download Rmd**, as well as code, numerical outputs and session information. The source can be recovered locally from the HTML; no separate Rmd upload is needed.
- The agreement reports preserve the saved notebook outputs and identify historical wording/provenance limitations. Labelling notebooks are omitted; the Supplementary Materials contain the prompts.

The annotation CSVs on OSF match the primary report fingerprints. Supportive and ancillary reruns also need three sampled-item CSVs, which are not in the current public deposit; their filenames and checksums are listed in the reproduction guide.

Input data and optional fit caches are not duplicated here. This GitHub update does not update the OSF deposit. No new annotations, model fits or resamples were generated for this publication package.

## Earlier public version

The previous public R implementation is preserved at the [pre-openmx-r1 tag](https://github.com/Veronica0206/gtheory-llm-evaluator-reliability/tree/pre-openmx-r1), based on commit `97dfc8d57999b0f16fc9a5c579d6f4eee037d2ba`. It is distinct from the supportive R1 HTML reports. Existing Git history is retained.

## Citation and license

Please cite the manuscript and identify the repository commit used in your work. The associated data deposit DOI is [10.17605/OSF.IO/K9CAJ](https://doi.org/10.17605/OSF.IO/K9CAJ).

Code is released under the [MIT License](LICENSE). Source datasets retain their respective licensing and attribution requirements; consult the OSF deposit and original corpus sources. Questions and reproducibility reports are welcome through GitHub issues.
