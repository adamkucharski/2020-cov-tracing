# 2020-cov-tracing

Code for modelling combined isolation, tracing and physical distancing measures. _Note: this is working repository, so code and data are likely to change over time_

### Quick start guide

First, set local path in R to GitHub directory, e.g.:
`
setwd("~/Documents/GitHub/2020-cov-tracing/V2_code")
`

Data loading and main model run script are in `scripts/contact_model.r`. This calls the following R file:

> `R/model_functions.R` - Function to sample individual contact networks, calculate secondary infections and effect of control measures.
