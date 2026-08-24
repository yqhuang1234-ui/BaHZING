# Compare the results.rds from any two dev/parallel_validation/runs/* folders
# produced by run_variant.R (e.g. serial vs. parallel today, or either against
# a run produced by a later refactor of BaHZING_Model).
#
# Usage:
#   Rscript dev/parallel_validation/compare_runs.R <run_dir_a> <run_dir_b>

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: Rscript compare_runs.R <run_dir_a> <run_dir_b>")
}
dir_a <- args[[1]]
dir_b <- args[[2]]

results_a <- readRDS(file.path(dir_a, "results.rds"))
results_b <- readRDS(file.path(dir_b, "results.rds"))

config_a <- yaml::read_yaml(file.path(dir_a, "config.yml"))
config_b <- yaml::read_yaml(file.path(dir_b, "config.yml"))

timing_a <- read.csv(file.path(dir_a, "timing.csv"))
timing_b <- read.csv(file.path(dir_b, "timing.csv"))

message("Comparing:")
message("  A: ", dir_a, " (git ", config_a$git_sha, ", variant ", config_a$variant, ")")
message("  B: ", dir_b, " (git ", config_b$git_sha, ", variant ", config_b$variant, ")")

if (!identical(dim(results_a), dim(results_b))) {
  message("Dimension mismatch: A = ", paste(dim(results_a), collapse = "x"),
          ", B = ", paste(dim(results_b), collapse = "x"))
} else {
  numeric_cols <- c("estimate", "bci_lcl", "bci_ucl", "pdir", "prope", "pmap")
  max_abs_diff <- sapply(numeric_cols, function(col) {
    max(abs(results_a[[col]] - results_b[[col]]))
  })

  summary_df <- data.frame(
    metric = c(numeric_cols, "identical", "elapsed_sec_A", "elapsed_sec_B"),
    value = c(
      unname(max_abs_diff),
      identical(results_a[, numeric_cols], results_b[, numeric_cols]),
      timing_a$elapsed_sec,
      timing_b$elapsed_sec
    )
  )
  print(summary_df)
}
