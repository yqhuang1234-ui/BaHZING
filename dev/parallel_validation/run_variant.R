# Run BaHZING_Model with package-default hyperparameters (n.chains = 3,
# n.adapt, n.iter.burnin, n.iter.sample, exposure_standardization,
# counterfactual_profiles, ROPE_range all left at their function defaults),
# using the same exposures/covariates as tests/testthat/test-BaHZING_Model.R,
# varying only `parallel` and a fixed, explicitly-assigned `seed`.
#
# Each call writes to its own runs/<git_sha>_<variant>_seed<seed>/ folder,
# tagged with the git commit that produced it, so results stay traceable and
# comparable across code versions (e.g. serial vs. parallel today, either vs.
# a refactored BaHZING_Model later). See compare_runs.R to diff two run
# folders. Not part of the package build (dev/ is .Rbuildignore'd) or CI; for
# the fast, deterministic, CI-safe equivalence check, see
# tests/testthat/test-BaHZING_Model_parallel.R.
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

config <- list(
  variant = variant,
  parallel = is_parallel,
  seed = seed,
  git_sha = git_sha,
  x = c("soft_drinks_dietnum", "diet_soft_drinks_dietnum"),
  covar = c("consent_age"),
  n.chains = 3,
  note = paste(
    "x and covar match tests/testthat/test-BaHZING_Model.R.",
    "n.chains explicitly set to 3 (same as function default) so parallel",
    "actually runs chains concurrently. n.adapt, n.iter.burnin,",
    "n.iter.sample, exposure_standardization, counterfactual_profiles,",
    "ROPE_range left at BaHZING_Model function defaults (not the test",
    "file's tiny overrides)."
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
