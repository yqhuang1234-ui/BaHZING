# Run BaHZING_Model with package-default hyperparameters (n.chains = 3,
# n.adapt, n.iter.burnin, n.iter.sample, exposure_standardization,
# counterfactual_profiles, ROPE_range all left at their function defaults),
# using the same exposures/covariates as tests/testthat/test-BaHZING_Model.R,
# varying only `parallel` and a fixed, explicitly-assigned `seed`.
#
# config.yml records the ACTUAL VALUE of every BaHZING_Model argument used
# for this run - not the word "default" - by introspecting
# formals(BaHZING_Model) at run time. This matters because "default" is not
# a stable fact: if a future commit changes a default value, an old run's
# config.yml that just said "left at function defaults" would become
# ambiguous about what actually ran. Reading the value at run time means
# each run's config.yml stays an accurate historical record regardless of
# what later commits change the defaults to - while the script itself still
# exercises "whatever the current defaults are" each time it's run, which is
# the point of this harness.
#
# Each call writes to its own runs/<date>/<git_sha>_<variant>_seed<seed>/
# folder, tagged with the git commit that produced it, so results stay
# traceable and comparable across code versions (e.g. serial vs. parallel
# today, either vs. a refactored BaHZING_Model later). See compare_runs.R to
# diff two run folders. Not part of the package build (dev/ is
# .Rbuildignore'd) or CI; for the fast, deterministic, CI-safe equivalence
# check, see tests/testthat/test-BaHZING_Model_parallel.R.
#
# Usage:
#   Rscript dev/parallel_validation/run_variant.R <variant_name> <parallel:TRUE|FALSE> [seed]
#
# Example:
#   Rscript dev/parallel_validation/run_variant.R serial_default   FALSE 20260728
#   Rscript dev/parallel_validation/run_variant.R parallel_default TRUE  20260728

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript run_variant.R <variant_name> <parallel:TRUE|FALSE> [seed]")
}
variant <- args[[1]]
is_parallel <- as.logical(args[[2]])
seed <- if (length(args) >= 3) as.integer(args[[3]]) else 20260728L

pkgload::load_all(".", quiet = TRUE)
library(yaml)

git_sha <- tryCatch(
  trimws(system("git rev-parse --short HEAD", intern = TRUE)),
  error = function(e) "unknown"
)

# The actual value BaHZING_Model's own formal default currently evaluates
# to, for any argument this script doesn't explicitly override. Every
# default here (NULL, a number, FALSE, or a simple c(...) call) is a
# self-contained literal with no external dependency, so eval() is safe.
model_default <- function(arg_name) eval(formals(BaHZING_Model)[[arg_name]])

# BaHZING_Model resolves n.cores = NULL to an actual core count internally
# (parallel::detectCores(logical=FALSE), clamped to n.chains) but doesn't
# return that value - mirrored here (same 3-line calculation) purely so
# config.yml can record what core count this run actually used, not just
# "auto-detected".
resolve_n_cores <- function(n_chains) {
  detected_cores <- parallel::detectCores(logical = FALSE)
  if (is.na(detected_cores) || detected_cores < 1) detected_cores <- 1
  max(1, min(n_chains, detected_cores))
}

n_chains <- 3

config <- list(
  variant = variant,
  parallel = is_parallel,
  seed = seed,
  git_sha = git_sha,
  x = c("soft_drinks_dietnum", "diet_soft_drinks_dietnum"),
  covar = c("consent_age"),
  n.chains = n_chains,
  n.adapt = model_default("n.adapt"),
  n.iter.burnin = model_default("n.iter.burnin"),
  n.iter.sample = model_default("n.iter.sample"),
  exposure_standardization = model_default("exposure_standardization"),
  counterfactual_profiles = model_default("counterfactual_profiles"),
  q = model_default("q"),
  return_all_estimates = model_default("return_all_estimates"),
  ROPE_range = model_default("ROPE_range"),
  n.cores = if (is_parallel) resolve_n_cores(n_chains) else model_default("n.cores"),
  verbose = FALSE,
  note = paste(
    "x and covar match tests/testthat/test-BaHZING_Model.R.",
    "n.chains explicitly set to 3 (same as function default) so parallel",
    "actually runs chains concurrently. Every other value above is",
    "BaHZING_Model's own current default, read via formals() at the time",
    "this run happened (not hardcoded here), so it stays accurate even if",
    "a later commit changes what the default is."
  )
)

run_dir <- file.path(
  "dev/parallel_validation/runs",
  format(Sys.Date(), "%Y-%m-%d"),
  paste0(git_sha, "_", variant, "_seed", seed)
)
dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
write_yaml(config, file.path(run_dir, "config.yml"))

data("iHMP_Reduced")
formatted_data <- Format_BaHZING(iHMP_Reduced)

model_args <- list(
  formatted_data = formatted_data,
  x = config$x,
  covar = config$covar,
  n.chains = config$n.chains,
  seed = config$seed,
  parallel = config$parallel,
  verbose = FALSE
)

message("Running variant '", variant, "' (parallel = ", is_parallel, ", seed = ", seed, ")...")
t <- system.time({
  results <- do.call(BaHZING_Model, model_args)
})

saveRDS(results, file.path(run_dir, "results.rds"))
write.csv(
  data.frame(variant = variant, parallel = is_parallel, elapsed_sec = t[["elapsed"]]),
  file.path(run_dir, "timing.csv"),
  row.names = FALSE
)

message("Done. Run folder: ", run_dir, " (elapsed: ", round(t[["elapsed"]], 1), "s)")
