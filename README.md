# 2020-cov-tracing

Code for modelling combined isolation, tracing and physical distancing measures. 

_Note: this is working repository, so non-archived code and data are likely to change over time_

### Quick start guide

First, set local path in R to GitHub directory, e.g.:
`
setwd("~/Documents/GitHub/2020-cov-tracing/")
`
Data loading and main model run script are in `scripts/contact_model.r`. This calls the following R file:

> `R/model_functions.R` - Function to sample individual contact networks, calculate secondary infections and effect of control measures.

### Archived versions

• Code to accompany [medRxiv V1 pre-print](https://www.medrxiv.org/content/10.1101/2020.04.23.20077024v1) is in `V1_code`, with file paths as above. An independent reproducibility CODECHECK has been performed on this version: http://doi.org/10.5281/zenodo.3767060

• Code to accompany peer-reviewed [paper published](https://www.thelancet.com/journals/laninf/article/PIIS1473-3099(20)30457-6/fulltext) in _Lancet Infectious Diseases_ is in `V2_code`, with file paths as above.

### Citation

[Kucharski AJ, Klepac P, Conlan AJK et al. Effectiveness of isolation, testing, contact tracing and physical distancing on reducing transmission of SARS-CoV-2 in different settings. Lancet Infectious Diseases, 2020](https://www.thelancet.com/journals/laninf/article/PIIS1473-3099(20)30457-6/fulltext)
