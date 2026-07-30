# Normalizes JAGS text for structural comparison: strips "#..."-to-end-of-line
# comments (JAGS ignores these entirely - they carry no statistical meaning,
# and the original hand-written model text is inconsistent in its comment
# wording/capitalization from level to level, e.g. "#zero inflation" vs
# "#Zero inflation" vs "#zero infl" vs "#zero component") and strips all
# whitespace entirely, not just collapsing it - the original text is also
# inconsistent about e.g. spaces after commas (genus's block omits the space
# genus/family/order/class all otherwise use) - so two model strings that
# define the same JAGS model compare equal regardless of any cosmetic
# formatting differences.
normalize_jags <- function(text) {
  text <- gsub("#[^\n]*", "", text)
  text <- gsub("\\s+", "", text)
  text
}

read_fixture <- function(name) {
  paste(readLines(test_path("fixtures", name)), collapse = "\n")
}

test_that(".bahzing_build_model_text reproduces the legacy with-covariates model", {
  generated <- .bahzing_build_model_text(default_taxa_levels, has_covar = TRUE)
  legacy <- read_fixture("legacy_model_with_covar.txt")
  expect_equal(normalize_jags(generated), normalize_jags(legacy))
})

test_that(".bahzing_build_model_text reproduces the legacy without-covariates model", {
  generated <- .bahzing_build_model_text(default_taxa_levels, has_covar = FALSE)
  legacy <- read_fixture("legacy_model_without_covar.txt")
  expect_equal(normalize_jags(generated), normalize_jags(legacy))
})

test_that(".bahzing_build_model_text places the terminal (no-parent) block at the broadest level for a custom hierarchy", {
  custom_levels <- c("Phylum", "Family", "Genus", "ASV")
  generated <- .bahzing_build_model_text(custom_levels, has_covar = FALSE)

  # Broadest level (Phylum) should have no parent-derived mean - a flat
  # dnorm(0, ...) prior, and no mu node keyed by its own index (p.r).
  # ("mu.phylum[f.r,...]" legitimately appears in Family's block below,
  # since Family borrows from Phylum - only Phylum's own index (p.r) should
  # never appear as a mu subscript.)
  expect_true(grepl("phylum\\.beta\\[p\\.r,p\\] ~ dnorm\\(0, phylum\\.tau\\[p\\.r\\]\\)", generated))
  expect_false(grepl("mu\\.[a-z]+\\[p\\.r", generated))

  # Narrowest level (ASV, this hierarchy's data-likelihood role) should have
  # the actual data likelihood and borrow its prior mean from its parent,
  # Genus - not from a hardcoded "genus"/"Genus.R"/"GenusData".
  expect_true(grepl("Y\\[i,r\\] ~ dnegbin", generated))
  expect_true(grepl("asv\\.beta\\[r,1:P\\]", generated))
  expect_true(grepl("mu\\.asv\\[r,p\\] <- inprod\\(genus\\.beta\\[1:Genus\\.R,p\\], GenusData\\[r,1:Genus\\.R\\]\\)", generated))

  # Mid-hierarchy level (Family) should borrow from its parent (Phylum)
  expect_true(grepl("mu\\.phylum\\[f\\.r,p\\] <- inprod\\(phylum\\.beta\\[1:Phylum\\.R,p\\], PhylumData\\[f\\.r,1:Phylum\\.R\\]\\)", generated))
})
