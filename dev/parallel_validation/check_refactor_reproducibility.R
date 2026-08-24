# One-time reproducibility check: does the current code - after this
# session's full BaHZING_Model() refactor (taxa_levels generalization, the
# JAGS-generator restructuring, the narrowest-level naming-consistency
# change) - produce identical output to the pre-refactor baseline
# (git_sha f00040d), under the exact same config and seed?
#
# Runs run_variant.R fresh (current code, current git_sha, the same full
# config as the baseline run: parallel=TRUE, seed=20260728, x/covar
# matching tests/testthat/test-BaHZING_Model.R), saved under
# dev/parallel_validation/runs/<today>/<current_git_sha>_post_session_refactor_seed20260728/
# (run_variant.R's own existing folder-naming convention - nothing new
# invented here), then diffs it against the baseline folder via the
# existing compare_runs.R.
#
# Usage:
#   Rscript dev/parallel_validation/check_refactor_reproducibility.R
#
# Takes ~10 minutes (one real BaHZING_Model() run, n.chains=3,
# n.iter.sample=5000, parallel=TRUE - the same cost as the original
# baseline run recorded in its timing.csv).

# Every path below is relative to the repository root (same convention
# run_variant.R itself requires) - fail fast with a clear message rather
# than a confusing "not found" if run from somewhere else, e.g. from
# inside dev/parallel_validation/ itself.
if (!file.exists("dev/parallel_validation/run_variant.R")) {
  stop(
    "This script must be run from the repository root, e.g.:\n",
    "  cd /Users/huangyanqi/Documents/bahzing_usc_github/BaHZING\n",
    "  Rscript dev/parallel_validation/check_refactor_reproducibility.R\n",
    "Current working directory: ", getwd()
  )
}

baseline_dir <- "dev/parallel_validation/runs/2026-07-28/f00040d_parallel_test_covar_seed20260728"
if (!dir.exists(baseline_dir)) {
  stop("Baseline run folder not found: ", baseline_dir)
}

variant <- "post_session_refactor"
seed <- 20260728

message("Running current code with the baseline's exact config (parallel = TRUE, seed = ", seed, ")...")
status <- system2("Rscript", c("dev/parallel_validation/run_variant.R", variant, "TRUE", seed))
if (status != 0) {
  stop("run_variant.R failed (exit status ", status, ")")
}

git_sha <- trimws(system("git rev-parse --short HEAD", intern = TRUE))
today <- format(Sys.Date(), "%Y-%m-%d")
new_dir <- file.path("dev/parallel_validation/runs", today,
                     paste0(git_sha, "_", variant, "_seed", seed))
if (!dir.exists(new_dir)) {
  stop("Expected new run folder not found: ", new_dir)
}

message("\nComparing against baseline (", basename(baseline_dir), ")...")
system2("Rscript", c("dev/parallel_validation/compare_runs.R", baseline_dir, new_dir))
