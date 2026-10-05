# Run ONCE before submitting the array job, so that the tasks do not all try to compile at the same time.
#   Rscript scripts/compile_model.R
library(cmdstanr); library(here)
cmdstan_model(here("stan", "single_area.stan"), include_paths = here("stan"), force_recompile = TRUE)
cat("compiled\n")
