# 2020-cov-tracing

Code for modelling combined isolation, tracing and physical distancing measures.

### Guide to files for `stoch_model`

Data loading and model run script is in `scripts/contact_model.r`. Note: need to change local path within this file to run. Calls the following R file:

> `R/model_functions.R` - Function to sample individual contact networks, calculate secondary infections and effect of control measures.
