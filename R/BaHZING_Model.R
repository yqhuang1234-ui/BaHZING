#' BaHZING_Model Function
#' This function implements the BaHZING model for microbiome data analysis.
#' @import R2jags
#' @import rjags
#' @import pscl
#' @import dplyr
#' @import tidyr
#' @import phyloseq
#' @import stringr
#' @importFrom utils globalVariables
#' @importFrom stats quantile update
#' @importFrom bayestestR p_direction p_rope p_map
#' @importFrom glue glue
#' @param formatted_data An object containing formatted microbiome data.
#' @param x A vector of column names of the exposures.
#' @param covar An optional vector of the column names of covariates.
#' @param n.chains An optional integer specifying the number of parallel chains
#' for the model in jags.model function. Default is 3.
#' @param n.adapt An optional integer specifying the number of iterations for
#' adaptation in jags.model function. Default is 100.
#' @param n.iter.burnin An optional integer specifying number of iterations in
#' update function. Default is 1,000.
#' @param n.iter.sample An optional integer specifying the number of iterations
#' in coda.samples function. Default is 5,000.
#' @param exposure_standardization Method for standardizing the exposures.
#' Should be one of "standard_normal" (the default), "quantile", or "none". If
#' "none", exposures are not standardized before analysis, and counterfactual
#' profiles must be specified by the user.
#' @param counterfactual_profiles A 2xP matrix or a vector with length of 2; P
#' is the number of exposures in x. If a 2xP matrix is provided,
#' the effect estimates for the mixture are interpreted as the estimated change
#' in the outcome when changing each exposure p in 1:P is changed from
#' `counterfactual_profiles[1,p]` to `counterfactual_profiles[2,p]`. If a vector of
#' length 2 is provided, the effect estimates for the mixture are interpreted as
#' the estimated change in the outcome when changing each exposure from
#' `counterfactual_profiles[1]` to `counterfactual_profiles[2]`. If
#' exposure_standardization = "standard_normal", then the default is c(-0.5, 0.5),
#' and the effect estimate is calculated based on increasing all exposures in
#' the mixture by one standard deviation. If exposure_standardization = "quantile",
#' then the default is c(0,1), and the effect estimate is calculated based on
#' increasing all exposures in the mixture by one quantile (where the number of
#' quantiles is based on the parameter q).
#' @param q An integer specifying the number of quantiles. Only required if
#' exposure_standardization = "quantile". If exposure_standardization =
#' "quantile" and q is not specified, then a default of q = 4 is used.
#' @param verbose If TRUE (default), function returns information a data quality
#' check.
#' @param return_all_estimates If FALSE (default), results do not include
#' the dispersion and omega estimates from the BaHZING model.
#' @param seed Optional seed used to generate per-chain JAGS RNG initial values.
#' @param parallel Logical; if TRUE, chains can be run in parallel with
#' parallel::mclapply on supported platforms.
#' @param n.cores Number of cores to use when parallel is TRUE. If NULL,
#' cores are detected automatically.
#' @param ROPE_range Region of practical equivalence (ROPE) for calculating
#' p_rope. Default is c(-0.1, 0.1).
#' @return A data frame containing results of the Bayesian analysis, with the
#' following columns:
#' - taxa_full: Full Taxa information, including all levels of the taxonomy.
#' Taxanomic levels are split by two underscores ('__').
#' - taxa_name: Taxa name, which is the last level of the taxonomy.
#' - domain: domain of the taxa.
#' - exposure: Exposure name (either one of  the individual exposures, or the
#' mixture).
#' - component: Zero inflated model estimate or the Count model estimate.
#' - estimate: Point estimate of the posterior distributions.
#' - bci_lcl: 95% Bayesian Credible Interval Lower Limit. Calculated as the
#' equal tailed interval of posterior distributions using the quantiles method.
#' - bci_ucl: 95% Bayesian Credible Interval Upper Limit. Calculated as the
#' equal tailed interval of posterior distributions using the quantiles method.
#' - p_direction: The Probability of Direction, calculated with `bayestestR`. A
#' higher value suggests a higher probability that the estimate is strictly
#' positive or negative. In other words, the closer the value to 1, the higher
#' the probability that the estimate is non-zero. Values can not be less than
#' 50%. From `bayestestR`: also known as the Maximum Probability of Effect
#' (MPE). This can be interpreted as the probability that a parameter (described
#' by its posterior distribution) is strictly positive or negative (whichever
#' is the most probable). Although differently expressed, this index is fairly
#' similar (i.e., is strongly correlated) to the frequentist p-value.
#' - p_rope: The probability that the estimate is not within the Region of
#' practical equivalence (ROPE), calculated with `bayestestR`. The proportion
#' of the whole posterior distribution that doesn't lie within the `ROPE_range`.
#' - p_map: Bayesian equivalent of the p-value, calculated with `bayestestR`.
#'  From `bayestestR`:  p_map is related to the odds that a parameter (described
#'  by its posterior distribution) has against the null hypothesis (h0) using
#'   Mills' (2014, 2017) Objective Bayesian Hypothesis Testing framework. It
#'   corresponds to the density value at the null (e.g., 0) divided by the
#'   density at the Maximum A Posteriori (MAP).
#' @export
#' @name BaHZING_Model

# Declare global variables
globalVariables(c("LibrarySize", "X2.5.", "X97.5.", "Mean",
                  "Exposure.Index", "taxa_index", "taxa_full",
                  "component", "estimate", "bci_lcl", "bci_ucl",
                  "domain", "taxa_name", "pdir","prope","pmap",
                  "name"))

# JAGS model text generators -------------------------------------------------
# BaHZING's hierarchical shrinkage model has one recurring pattern: every
# taxonomic level's exposure effect is drawn from a normal distribution
# centered on its parent level's corresponding effect, borrowed via a binary
# taxonomy indicator matrix (inprod against e.g. GenusData). The narrowest
# level in taxa_levels (species, by default) is the one exception with a
# parent to borrow from AND its own data likelihood; the broadest level
# (phylum, by default) is the other exception, with no parent to borrow from
# at all. These generators build the JAGS text for an arbitrary-length
# taxa_levels hierarchy from these three roles (narrowest/data-likelihood,
# mid-hierarchy, broadest/terminal), replacing what used to be two ~180-line
# hand-maintained, near-identical literal JAGS strings (one per
# with/without-covariates variant) that assumed exactly 6 fixed level names.
#
# All generators use glue::glue(.open="<<", .close=">>", .trim=FALSE): custom
# delimiters because JAGS's own "{ }" block syntax would otherwise collide
# with glue's default interpolation markers, and .trim=FALSE because glue's
# default (.trim=TRUE) strips leading whitespace and the trailing newline -
# harmless for a block that ends the whole template, but it silently glues
# adjacent fragments onto the same line wherever one block's output is
# spliced into another.

#' Build the precision-prior and g-estimation lines for one taxonomic level
#'
#' Generates the JAGS lines for one level's precision priors (count-model and
#' zero-inflation dispersion) and its g-estimation block (the mixture
#' contrast computed from the low/high counterfactual `profiles`). Shared by
#' every level via [.bahzing_level_block()] and [.bahzing_species_block()].
#'
#' @param level Character. The level's variable-name prefix, e.g. `"genus"`.
#'   Used for the g-estimation nodes (`<level>.psi`, `<level>.eta.*`), which
#'   are always named after the level itself.
#' @param idx Character. The JAGS loop index variable for this level, e.g.
#'   `"g.r"`.
#' @param tau_prefix Character. Prefix used for the precision nodes
#'   (`<tau_prefix>tau`/`<tau_prefix>sigma`). Defaults to `paste0(level,".")`
#'   (e.g. genus's precision nodes are `genus.tau`/`genus.sigma`) - overridden
#'   to `""` by [.bahzing_species_block()], since the narrowest/species-role
#'   level's precision nodes are unprefixed (`tau`/`sigma`) in BaHZING's
#'   original model, unlike every other level.
#' @return A length-1 character string ending in a trailing newline.
#' @keywords internal
#' @noRd
.bahzing_precision_gestimation_block <- function(level, idx, tau_prefix = paste0(level, ".")) {
  glue::glue(
"      # prior on precision
      <<tau_prefix>>tau[<<idx>>] <- 1/(<<tau_prefix>>sigma[<<idx>>]*<<tau_prefix>>sigma[<<idx>>])
      <<tau_prefix>>sigma[<<idx>>] ~ dunif(0,3)
      <<tau_prefix>>tau.zero[<<idx>>] <- 1/(<<tau_prefix>>sigma.zero[<<idx>>]*<<tau_prefix>>sigma.zero[<<idx>>])
      <<tau_prefix>>sigma.zero[<<idx>>] ~ dunif(0,3)

      # g-estimation
      <<level>>.eta.low[<<idx>>] <- inprod(<<level>>.beta[<<idx>>,1:P], profiles[1,1:P])
      <<level>>.eta.high[<<idx>>] <- inprod(<<level>>.beta[<<idx>>,1:P], profiles[2,1:P])
      <<level>>.psi[<<idx>>] <- <<level>>.eta.high[<<idx>>]-<<level>>.eta.low[<<idx>>]
      # zero-inflation
      <<level>>.eta.low.zero[<<idx>>] <- inprod(<<level>>.beta.zero[<<idx>>,1:P], profiles[1,1:P])
      <<level>>.eta.high.zero[<<idx>>] <- inprod(<<level>>.beta.zero[<<idx>>,1:P], profiles[2,1:P])
      <<level>>.psi.zero[<<idx>>] <- <<level>>.eta.high.zero[<<idx>>]-<<level>>.eta.low.zero[<<idx>>]
", .open = "<<", .close = ">>", .trim = FALSE)
}

#' Build the exposure-effect prior for one taxonomic level
#'
#' Generates `for(p in 1:P) { <level>.beta[idx,p] ~ dnorm(...) }`: each
#' level's exposure effect is drawn from a normal distribution centered on a
#' value borrowed from the parent level, computed via `inprod()` over the
#' parent's beta and the taxonomy indicator matrix - both the count-model and
#' zero-inflation components. If `parent` is `NULL` (only the broadest/
#' terminal level, phylum by default), there's nothing to borrow from, so the
#' prior falls back to a fixed `dnorm(0, tau)` instead.
#'
#' Shared by every level including the narrowest/species-role level, via
#' [.bahzing_level_block()] and [.bahzing_species_block()] respectively - see
#' `mu_owner` below for why that level needs one extra parameter to reuse
#' this.
#'
#' @param level Character. This level's variable-name prefix (e.g.
#'   `"genus"`), used for `<level>.beta`/`<level>.tau`.
#' @param idx Character. This level's JAGS loop index variable (e.g.
#'   `"g.r"`).
#' @param parent Character or `NULL`. The parent level's variable-name prefix
#'   to borrow the prior mean from (e.g. `"family"` for genus). `NULL` only
#'   for the broadest/terminal level, which has no parent.
#' @param parent_R Character or `NULL`. The JAGS scalar holding the parent
#'   level's taxon count (e.g. `"Family.R"`), used to size the `inprod()`.
#'   Required whenever `parent` is supplied.
#' @param parent_data Character or `NULL`. The taxonomy indicator matrix
#'   linking this level's taxa to the parent's (e.g. `"FamilyData"`).
#'   Required whenever `parent` is supplied.
#' @param mu_owner Character. Whose name labels the intermediate "borrowed
#'   mean" node (`mu.<mu_owner>`). Defaults to `parent`, matching the
#'   original model's convention for mid-hierarchy levels (e.g. genus's node
#'   is `mu.family`, named after *its* parent). The narrowest/species-role
#'   level is the one exception in the original model - its node is
#'   `mu.species` (named after *itself*, not `mu.genus`) - so
#'   [.bahzing_species_block()] overrides this explicitly to its own level
#'   name. This is a naming inconsistency inherited from the original model,
#'   preserved rather than "fixed", since `mu.*` is never a monitored/
#'   extracted node - only its label differs, not the value.
#' @param comment Character or `NULL`. If supplied, prepended as a standalone
#'   comment line above the `for(p in 1:P)` loop.
#' @return A length-1 character string ending in a trailing newline.
#' @keywords internal
#' @noRd
.bahzing_exposure_prior_block <- function(level, idx, parent = NULL, parent_R = NULL,
                                            parent_data = NULL, mu_owner = parent,
                                            comment = NULL, tau_prefix = paste0(level, ".")) {
  # plain paste0, not glue - glue's default .trim strips leading whitespace
  # and the trailing newline, which would collapse this onto the same line
  # as whatever follows it.
  comment_line <- if (!is.null(comment)) paste0("      # ", comment, "\n") else ""
  body <- if (!is.null(parent)) {
    glue::glue(
"        <<level>>.beta[<<idx>>,p] ~ dnorm(mu.<<mu_owner>>[<<idx>>,p], <<tau_prefix>>tau[<<idx>>])
        mu.<<mu_owner>>[<<idx>>,p] <- inprod(<<parent>>.beta[1:<<parent_R>>,p], <<parent_data>>[<<idx>>,1:<<parent_R>>])
        #Zero inflation component
        <<level>>.beta.zero[<<idx>>,p] ~ dnorm(mu.<<mu_owner>>.zero[<<idx>>,p], <<tau_prefix>>tau.zero[<<idx>>])
        mu.<<mu_owner>>.zero[<<idx>>,p] <- inprod(<<parent>>.beta.zero[1:<<parent_R>>,p], <<parent_data>>[<<idx>>,1:<<parent_R>>])
", .open = "<<", .close = ">>", .trim = FALSE)
  } else {
    glue::glue(
"        <<level>>.beta[<<idx>>,p] ~ dnorm(0, <<tau_prefix>>tau[<<idx>>])
        #Zero inflation component
        <<level>>.beta.zero[<<idx>>,p] ~ dnorm(0, <<tau_prefix>>tau.zero[<<idx>>])
", .open = "<<", .close = ">>", .trim = FALSE)
  }
  glue::glue(
"<<comment_line>>      for(p in 1:P) {
<<body>>      }
", .open = "<<", .close = ">>", .trim = FALSE)
}

#' Build one mid-hierarchy or broadest/terminal level's full JAGS block
#'
#' Assembles one level's complete `for(idx in 1:level_R) { ... }` block: a
#' comment header, the exposure-effect prior
#' ([.bahzing_exposure_prior_block()]), and the precision-prior/g-estimation
#' lines ([.bahzing_precision_gestimation_block()]). Handles both mid-
#' hierarchy levels (`parent` supplied - borrows its prior mean from the
#' parent level) and the broadest/terminal level (`parent = NULL` - the top
#' of the hierarchy, with nothing to borrow from).
#'
#' Not used for the narrowest/species-role level - that level has its own
#' likelihood/dispersion/intercept/covariate structure with no equivalent at
#' any other level, so it gets its own function,
#' [.bahzing_species_block()], which reuses the same two shared sub-block
#' generators this function uses.
#'
#' @param level Character. This level's variable-name prefix, e.g. `"genus"`.
#' @param level_label Character. Human-readable label for this level's
#'   comment header, e.g. `"Genus"` (produces `# Genus level`).
#' @param idx Character. This level's JAGS loop index variable, e.g. `"g.r"`.
#' @param level_R Character. The JAGS scalar holding this level's taxon
#'   count, e.g. `"Genus.R"` - bounds the `for(idx in 1:level_R)` loop.
#' @param parent Character or `NULL`. The parent level's variable-name prefix
#'   (e.g. `"family"` for genus). `NULL` only for the broadest/terminal level.
#' @param parent_R Character or `NULL`. The parent level's taxon-count
#'   scalar (e.g. `"Family.R"`). Required whenever `parent` is supplied.
#' @param parent_data Character or `NULL`. The taxonomy indicator matrix
#'   linking this level to the parent (e.g. `"FamilyData"`). Required
#'   whenever `parent` is supplied.
#' @return A length-1 character string: this level's complete JAGS block.
#' @keywords internal
#' @noRd
.bahzing_level_block <- function(level, level_label, idx, level_R,
                                  parent = NULL, parent_R = NULL, parent_data = NULL) {
  exposure_prior  <- .bahzing_exposure_prior_block(level, idx, parent, parent_R, parent_data)
  precision_gestim <- .bahzing_precision_gestimation_block(level, idx)

  glue::glue(
"    # <<level_label>> level
    for(<<idx>> in 1:<<level_R>>) {
<<exposure_prior>><<precision_gestim>>    }

", .open = "<<", .close = ">>", .trim = FALSE)
}

#' Build the narrowest/species-role level's JAGS block - the data-likelihood level
#'
#' The narrowest level in `taxa_levels` is tied directly to the observed
#' data, so unlike [.bahzing_level_block()] (mid-hierarchy/terminal levels,
#' which are pure prior/hierarchy blocks), this one also includes: the
#' zero-inflated negative binomial likelihood (`Y ~ dnegbin(...)`,
#' `zero ~ dbern(...)`), the dispersion prior (`disp`), the intercept priors
#' (`alpha`/`alpha.zero`), and - when `has_covar` is `TRUE` - the covariate
#' terms and their priors (`delta`/`delta.zero`). This is what collapses the
#' original model's separate with/without-covariates literal text into one
#' function.
#'
#' Reuses the same two shared sub-block generators every other level uses:
#' [.bahzing_exposure_prior_block()] (with `mu_owner = level` and a
#' `"prior on exposure effects"` comment, matching the original model's
#' species-specific naming/comment conventions) and
#' [.bahzing_precision_gestimation_block()] (with `tau_prefix = ""`, since
#' this level's precision nodes are unprefixed in the original model, unlike
#' every other level - see that function's docs).
#'
#' @param has_covar Logical. Whether covariates are included in the model.
#'   When `TRUE`, adds the `delta`/`delta.zero` covariate terms to the
#'   likelihood and their `dnorm(0, 1.0E-02)` priors; when `FALSE`, omits
#'   them entirely (matching the original model's with/without-covariates
#'   variants). Confirmed by diffing the two original literal model strings:
#'   covariate terms only ever appear in this block - every other level's
#'   text is identical between the with/without-covariates variants.
#' @param level Character. This level's variable-name prefix (`"species"` by
#'   default).
#' @param parent Character. The parent level's variable-name prefix to
#'   borrow the prior mean from (`"genus"` by default).
#' @param parent_R Character. The parent level's taxon-count scalar
#'   (`"Genus.R"` by default).
#' @param parent_data Character. The taxonomy indicator matrix linking this
#'   level to the parent (`"GenusData"` by default).
#' @return A length-1 character string: the complete data-likelihood-level
#'   JAGS block.
#' @keywords internal
#' @noRd
.bahzing_species_block <- function(has_covar, level = "species", parent = "genus",
                                    parent_R = "Genus.R", parent_data = "GenusData") {
  covar_lambda <- if (has_covar) " + inprod(delta[r, 1:Q], W[i,1:Q])" else ""
  covar_pi     <- if (has_covar) " + inprod(delta.zero[r, 1:Q], W[i,1:Q])" else ""
  covar_prior  <- if (has_covar)
"
      # prior on covariate effects
      for(q in 1:Q) {
        delta[r,q] ~ dnorm(0, 1.0E-02)
        delta.zero[r,q] ~ dnorm(0, 1.0E-02)
      }
" else ""

  exposure_prior <- .bahzing_exposure_prior_block(level, "r", parent = parent,
                                                    parent_R = parent_R, parent_data = parent_data,
                                                    mu_owner = level, tau_prefix = "",
                                                    comment = "prior on exposure effects")
  precision_gestim <- .bahzing_precision_gestimation_block(level, "r", tau_prefix = "")

  glue::glue(
"    for(r in 1:R) {
      for(i in 1:N) {
        Y[i,r] ~ dnegbin(mu[i,r], disp[r])
        mu[i,r] <- disp[r]/(disp[r]+(1-zero[i,r])*lambda[i,r]) - 0.000001*zero[i,r]
        log(lambda[i,r]) <- alpha[r] + inprod(<<level>>.beta[r,1:P], X.q[i,1:P])<<covar_lambda>> + log(L[i,1])

        # zero-inflation
        zero[i,r] ~ dbern(pi[i,r])
        logit(pi[i,r]) <- alpha.zero[r] + inprod(<<level>>.beta.zero[r,1:P], X.q[i,1:P])<<covar_pi>> + log(L[i,1])
      }
      # prior on dispersion parameter
      disp[r] ~ dunif(0,50)

      # prior on intercept
      alpha[r] ~ dnorm(0, 1.0E-02)
      alpha.zero[r] ~ dnorm(0, 1.0E-02)
<<covar_prior>>
<<exposure_prior>>
<<precision_gestim>>    }

", .open = "<<", .close = ">>", .trim = FALSE)
}

#' Assemble the complete JAGS model text for BaHZING_Model()
#'
#' The top-level entry point for the templating system: builds every
#' taxonomic level's block (the narrowest level via
#' [.bahzing_species_block()]; every other level via
#' [.bahzing_level_block()], each passed its own parent level's info, with
#' the broadest level in `taxa_levels` passed `parent = NULL`) and
#' concatenates them narrowest-to-broadest into one complete `model { ... }`
#' string, matching the original hand-written text's declaration order,
#' ready to pass to `jags.model(file = textConnection(...))`.
#'
#' @param taxa_levels Character vector naming the taxonomic hierarchy,
#'   ordered broadest to narrowest (e.g. `default_taxa_levels`). Must have
#'   already been validated with `validate_taxa_levels()`.
#' @param has_covar Logical. Whether the model includes covariates - passed
#'   straight through to [.bahzing_species_block()], the only level that
#'   varies its text based on this (see that function for why).
#' @return A length-1 character string: the complete JAGS model text, from
#'   `model {` through the final closing `}`.
#' @keywords internal
#' @noRd
.bahzing_build_model_text <- function(taxa_levels, has_covar) {
  n <- length(taxa_levels)
  blocks <- vector("list", n)
  for (i in seq_len(n)) {
    level <- tolower(taxa_level_name(taxa_levels, i))
    parent <- if (i > 1) taxa_level_name(taxa_levels, i - 1) else NULL
    if (i == n) {
      # Narrowest level: the data-likelihood ("species") role. Always has a
      # parent (n >= 2 is enforced by validate_taxa_levels()).
      blocks[[i]] <- .bahzing_species_block(has_covar, level,
                       tolower(parent), paste0(parent, ".R"), paste0(parent, "Data"))
    } else {
      idx <- paste0(substr(level, 1, 1), ".r")
      level_label <- taxa_level_name(taxa_levels, i)
      blocks[[i]] <- .bahzing_level_block(level, level_label, idx, paste0(level_label, ".R"),
                       if (!is.null(parent)) tolower(parent) else NULL,
                       if (!is.null(parent)) paste0(parent, ".R") else NULL,
                       if (!is.null(parent)) paste0(parent, "Data") else NULL)
    }
  }
  # Emit narrowest -> broadest (reverse of taxa_levels' broadest-first
  # storage order) to match the original text's declaration order. JAGS
  # registers stochastic nodes in declaration order, which can affect the
  # sampler's RNG draw sequence for a given seed even when the model is
  # mathematically equivalent - preserving order keeps a seeded before/after
  # diff exactly reproducible rather than merely statistically similar.
  paste0("model {\n", paste(rev(blocks), collapse = ""), "  }")
}

BaHZING_Model <- function(formatted_data,
                          x,
                          covar = NULL,
                          exposure_standardization = NULL,
                          n.chains = 3,
                          n.adapt = 100,
                          n.iter.burnin = 1000,
                          n.iter.sample = 5000,
                          counterfactual_profiles = NULL,
                          q = NULL,
                          verbose = TRUE,
                          return_all_estimates = FALSE,
                          seed = NULL,
                          parallel = FALSE,
                          n.cores = NULL,
                          ROPE_range = c(-0.1, 0.1)) {

  # 1. Check input data ----
  # Extract metadata file from formatted data
  exposure_covar_dat <- data.frame(formatted_data$Table)

  # Create covariate dataframe
  if (!is.null(covar)){
    W <- data.frame(exposure_covar_dat[covar])
    Q <- ncol(W)
  }

  # Create exposure dataframe
  if(!all(x %in% colnames(exposure_covar_dat))) {
    stop("Not all exposured are found in the formatted data")
  }
  X <- exposure_covar_dat[x]
  P <- ncol(X)

  # Set exposure_standardization if missing
  if(is.null(exposure_standardization)) {
    exposure_standardization = "standard_normal"
  }
  # Check exposure_standardization is within the bounds
  if(!(exposure_standardization %in% c("standard_normal", "quantile", "none"))){
    stop("exposure_standardization must be either standard_normal, quantile, or none")
  }

  # Give warning if exposure_standardization == "standard_normal" and q is specified
  if(exposure_standardization == "standard_normal" & !is.null(q)){
    message("Note: q is not required when exposure_standardization is standard_normal. q will be ignored.")
  }

  # If standardization is "none", then check to make sure that counterfactual_profiles is specified
  if(exposure_standardization  == "none" & is.null(counterfactual_profiles)){
    stop("counterfactual_profiles must be speficied if exposure_standardization is none.")
  }

  # 2. Counterfactual profiles -----------------------------------
  ## 2.a Set counterfactual_profiles if missing --------
  if(is.null(counterfactual_profiles)) {
    if(exposure_standardization == "standard_normal") {
      counterfactual_profiles <- c(-0.5,0.5)
    } else {
      counterfactual_profiles <- c(0,1)
    }
  }

  ## 2.b Check counterfactual profiles structure ----
  if (is.matrix(counterfactual_profiles)) { # Checks if matrix
    if (nrow(counterfactual_profiles) != 2) {
      stop("counterfactual_profiles must have 2 rows when provided as a matrix.")
    }
    if (ncol(counterfactual_profiles) != P) {
      stop(paste0("When provided as a matrix, the number of columns in counterfactual_profiles must be equal to the number of exposures in the model."))
    }
    if (!is.numeric(counterfactual_profiles)) {
      stop("counterfactual_profiles must be numeric.")
    }
  } else { # Checks if is vector
    if (is.vector(counterfactual_profiles)) {
      if (!is.numeric(counterfactual_profiles)) {
        stop("counterfactual_profiles must be numeric.")
      }
      if (length(counterfactual_profiles) != 2) {
        stop(paste0("counterfactual_profiles must have 2 elements when provided as a vector."))
      }
    } else {
      stop("counterfactual_profiles must be either a numeric 2xP matrix or a numeric vector with length P.")
    }
  }

  ## 2.c Set profiles --------------------------------
  if(is.matrix(counterfactual_profiles)) {
    profiles = counterfactual_profiles
  } else {
    profiles <- rbind(rep(counterfactual_profiles[1], P),
                      rep(counterfactual_profiles[2], P))
  }

  chain_inits <- if (!is.null(seed)) {
    lapply(seq_len(n.chains), function(i) {
      list(.RNG.name = "base::Wichmann-Hill",
           .RNG.seed = seed + i)
    })
  } else {
    NULL
  }

  if (isTRUE(parallel) && is.null(n.cores)) {
    detected_cores <- parallel::detectCores(logical = FALSE)
    if (is.na(detected_cores) || detected_cores < 1) {
      detected_cores <- 1
    }
    n.cores <- max(1, min(n.chains, detected_cores))
  }

  # Give warning if using qualtiles but counterfactuals < 0
  if(exposure_standardization=="quantiles" &
     any(counterfactual_profiles<0 |
         counterfactual_profiles > q+1)){
    warning("Note: Quantiles are used, but counterfactual_profiles includes values outside of range. Estimates will be calculated based on the specified counterfactual_profiles, but results may not be interpretable.")
  }


  # 3. Scale exposures --------------------------------
  # Set default q if not provided
  if(is.null(q) & exposure_standardization == "quantile") {
    q = 4
  }

  # If using quantiles, quantize X
  if(exposure_standardization=="quantile") {
    probs <- seq(0, 1, length.out = q + 1)
    X.q <- apply(X, 2, function(v) {
      cut(v, breaks = c(-Inf, quantile(v, probs = probs, include.lowest = FALSE)), labels = FALSE)
    })
  }

  #If not quantized and not standardized, scale X
  if(exposure_standardization=="standard_normal") {
    X.q <- scale(X)
  }

  #If none, no scaling
  if(exposure_standardization == "none") {
    X.q <- X
  }

  # 4. Format microbiome matricies ----
  taxa_levels <- formatted_data$taxa_levels
  if (is.null(taxa_levels)) {
    stop("formatted_data has no $taxa_levels; regenerate it with the current Format_BaHZING().")
  }
  validate_taxa_levels(taxa_levels)
  narrowest_level <- taxa_level_name(taxa_levels, length(taxa_levels))

  # Which columns of exposure_covar_dat are taxon-count (outcome) data,
  # read directly from Format_BaHZING()'s own authoritative record rather
  # than inferred via string-matching (e.g. assuming every such column
  # contains "k__", which silently breaks if Kingdom/Domain isn't part of
  # the input taxonomy at all - Format_BaHZING() already computes this
  # exact list when it builds these columns, so there's no need to guess).
  taxon_columns <- formatted_data$taxon_columns
  if (is.null(taxon_columns) || length(taxon_columns) == 0) {
    stop("formatted_data has no $taxon_columns; regenerate it with the current Format_BaHZING().")
  }

  #Create outcome dataframe
  Y <- exposure_covar_dat[, taxon_columns]
  N <- nrow(Y)
  R <- ncol(Y)

  # One data.frame + taxon-count per level except the narrowest (species-role)
  # level, e.g. GenusData/Genus.R, FamilyData/Family.R, ... - read from the
  # binary incidence matrices Format_BaHZING() built, keyed generically via
  # hierarchy_matrix_name() so this never has to hardcode a level name.
  level_data <- list()
  for (i in seq_len(length(taxa_levels) - 1)) {
    lvl <- taxa_level_name(taxa_levels, i)
    mat_name <- hierarchy_matrix_name(taxa_levels, i)
    if (is.null(formatted_data[[mat_name]])) {
      stop(sprintf("formatted_data is missing '%s'; regenerate it with Format_BaHZING().", mat_name))
    }
    df <- as.data.frame(t(formatted_data[[mat_name]]))
    level_data[[lvl]] <- list(data = df, R = ncol(df))
  }

  ## Create Library Size Offset
  L <- exposure_covar_dat[, taxon_columns]
  L <- L %>%
    mutate(LibrarySize=rowSums(across(everything())))
  L <- L %>%
    select(LibrarySize)

  # 5. Return "Sanity" Messages ----
  if(verbose == TRUE){
    message("#### Checking input data ####")
    message("Exposure and Covariate Data:")
    message(paste0("- Total sample size: ", N))
    message(paste0("- Number of exposures: ", P))

    message("Microbiome Data:")
    # Narrowest (of the non-species-role levels) to broadest, e.g.
    # genus, family, order, class, phylum for the default hierarchy -
    # rev() of level_data's broadest-first insertion order above.
    for (lvl in rev(names(level_data))) {
      message(paste0("- Number of unique ", tolower(lvl), " in data: ", level_data[[lvl]]$R))
    }

    message("#### Running BaHZING with the following parameters #### ")
    if (exposure_standardization == "standard_normal"){
      message("Exposure standardization: Standard Normal")
    }
    if (exposure_standardization == "none"){
      message("Exposure standardization: None")
    }
    if (exposure_standardization == "quantile"){
      message(paste0("Exposure standardization: Quantiles, with q = ", q))
    }
  }

  # 6. Run Model ----
  BHRM.microbiome <- .bahzing_build_model_text(taxa_levels, has_covar = !is.null(covar))

  ### Run JAGs Estimation ----
  # set up for JAGs based on taxonomy
  jdata <- list(N=N, Y=Y, R=R, X.q=X.q, P=P, profiles=profiles, L=L)
  for (lvl in names(level_data)) {
    jdata[[paste0(lvl, "Data")]] <- level_data[[lvl]]$data
    jdata[[paste0(lvl, ".R")]] <- level_data[[lvl]]$R
  }
  if (!is.null(covar)) {
    jdata$Q <- Q
    jdata$W <- W
  }

  # Narrowest -> broadest prefix order, matching the original hardcoded var.s
  # (e.g. species, genus, family, order, class, phylum for the default
  # hierarchy): narrowest_level first, then level_data's names (broadest ->
  # narrowest insertion order, built in section 4) reversed.
  level_prefixes <- tolower(c(narrowest_level, rev(names(level_data))))
  var.s <- c(paste0(level_prefixes, ".beta"), paste0(level_prefixes, ".beta.zero"),
             paste0(level_prefixes, ".psi"), paste0(level_prefixes, ".psi.zero"),
             "disp")

  if (isTRUE(parallel) && n.cores > 1 && n.chains > 1) {
    model.fit <- parallel::mclapply(seq_len(n.chains), function(chain_index) {
      chain_model <- jags.model(file=textConnection(BHRM.microbiome),
                                data=jdata,
                                n.chains=1,
                                n.adapt=n.adapt,
                                quiet=F,
                                inits=chain_inits[[chain_index]])
      update(chain_model, n.iter=n.iter.burnin, progress.bar="text")
      coda.samples(model=chain_model,
                   variable.names=var.s,
                   n.iter=n.iter.sample,
                   thin=1,
                   progress.bar="text")
    }, mc.cores = n.cores)
    model.fit <- do.call(coda::mcmc.list, lapply(model.fit, function(x) x[[1]]))
  } else {
    model.fit <- jags.model(file=textConnection(BHRM.microbiome),
                            data=jdata,
                            n.chains=n.chains,
                            n.adapt=n.adapt,
                            quiet=F,
                            inits=chain_inits)
    update(model.fit, n.iter=n.iter.burnin, progress.bar="text")
    model.fit <- coda.samples(model=model.fit,
                              variable.names=var.s,
                              n.iter=n.iter.sample,
                              thin=1,
                              progress.bar="text")
  }

  # 7. summarize results -------------------------------------------------------
  ## Calculate Mean, SD, and quantiles ----
  r <- summary(model.fit)
  results <- data.frame(round(r$statistics[,1:2],3), round(r$quantiles[,c(1,5)],3))

  ## Calculate HPD intervals (not included in current version) ------
  # x1 <- HPDinterval(model.fit[[1]], prob = 0.95)  %>% as.data.frame()

  ## Calculate p-values  ------
  post_dist <-  as.data.frame(model.fit[[1]])[, grep("beta|zero|psi", colnames(model.fit[[1]]))]

  ### p_direction ----
  pdir <- apply(post_dist, 2, function(x){
    p_direction(x = x, threshold = 0.05) %>%
      as.numeric()})

  ### p_rope ----
  prope <- apply(post_dist, 2, function(x){
    p_rope(x = x, rope = ROPE_range)$p_ROPE})

  ### p_map ----
  pmap <- apply(post_dist, 2, function(x){
    p_map(x = x) %>%
      as.numeric()})

  # Get dataframe of p-values
  p_value_df <- data.frame(name = names(pdir),
                           pdir = pdir,
                           prope = prope,
                           pmap = pmap)

  # # "Significant" example
  # bayestestR::describe_posterior(post_dist$`family.beta[1,1]`)
  # (pdirection <- bayestestR::p_direction(post_dist$`family.beta[1,1]`, as_p = TRUE)) %>% as.numeric()
  # (prope <- bayestestR::p_rope(post_dist$`family.beta[1,1]`, rope = c(-0.1, 0.1), as_p = TRUE) %>% as.numeric())
  # (pmap <- bayestestR::p_map(post_dist$`family.beta[1,1]`)) %>% as.numeric()
  # (ps <- bayestestR::p_significance(post_dist$`family.beta[1,1]`, rope = c(-0.1, 0.1), as_p = TRUE))
  #
  # # "Non-significant" example
  # bayestestR::describe_posterior(post_dist$`species.beta[189,2]`)
  # bayestestR::p_direction(post_dist$`species.beta[189,2]`, as_p = TRUE)
  # bayestestR::p_rope(post_dist$`species.beta[189,2]`, rope = c(-0.1, 0.1))

  ## Create "component" variable ----
  results <- results %>%
    mutate(
      component=case_when(
        grepl("zero",rownames(.)) ~ "Zero-inflation model coefficients",
        grepl("beta",rownames(.)) ~ "Count model coefficients",
        grepl("psi",rownames(.))  ~ "Count model coefficients",
        grepl("disp",rownames(.)) ~ "Dispersion",
        grepl("omega",rownames(.)) ~ "Omega",
        TRUE ~ "Other"))

  # Calculate significance based on Bayesian Interval, rename variables (removed- jg 02_13_25 in place of p-values)
  results <- results %>%
    # mutate(bci_inc_zero = case_when(
    #   grepl("disp",rownames(.)) ~ NA_character_,
    #   grepl("omega",rownames(.)) ~ NA_character_,
    #   (X2.5.<0 & X97.5.<0) | (X2.5.>0 & X97.5.>0) ~ "BCI ",
    #   TRUE ~ "N.S.")) %>%
    rename(estimate = Mean,
           bci_lcl = X2.5.,
           bci_ucl = X97.5.)

  # Combine results df with p-value df
  results$name = rownames(results)
  results2 <- left_join(results, p_value_df, by = "name")
  rownames(results2) <- results2$name
  results2 <- results2 %>% select(-name)

  # Calculate odds ratios-- removed --
  # results2 <- results2 %>%
  #   mutate(OR=exp(Mean),
  #          OR.ll=exp(X2.5.),
  #          OR.ul=exp(X97.5.))

  #Format output names
  # One taxon-name vector per level, keyed by level name (proper case, e.g.
  # "Genus") - narrowest level uses Y's colnames (the actual outcome/species
  # columns), every other level uses its level_data data.frame's colnames
  # (built in section 4).
  level_labels <- lapply(taxa_levels, function(lvl) {
    if (identical(lvl, narrowest_level)) colnames(Y) else colnames(level_data[[lvl]]$data)
  })
  names(level_labels) <- taxa_levels
  exposure <- colnames(X)
  results2$taxa_index <- str_remove(rownames(results2),"..*\\[")
  results2$taxa_index <- str_remove(results2$taxa_index,",.*$")
  results2$taxa_index <- str_remove(results2$taxa_index,"]")
  results2$taxa_index <- as.numeric(results2$taxa_index)
  results2$Exposure.Index <- str_remove(rownames(results2),"..*\\,")
  results2 <- results2 %>%
    mutate(Exposure.Index=ifelse(
      grepl("disp",Exposure.Index)|
        grepl("omega",Exposure.Index)|
        grepl("psi",Exposure.Index),
      NA,Exposure.Index))
  results2$Exposure.Index <- str_remove(results2$Exposure.Index,"]")
  results2$Exposure.Index <- as.numeric(results2$Exposure.Index)

  # Resolve which taxonomic level a JAGS parameter row belongs to, from its
  # rowname's "<level>."-prefixed parameter name (e.g. "genus.beta[12,1]" ->
  # "Genus") - anchored to the start of the parameter name so a level whose
  # name is a substring of another's can't silently mis-attribute (`stop()`s
  # instead). "disp"/"omega" rows have no level prefix at all - they're
  # dispersion/legacy parameters computed per narrowest-level taxon, so they
  # resolve to the narrowest level, matching the original model's behavior.
  resolve_level_for_row <- function(rowname) {
    if (grepl("disp", rowname) || grepl("omega", rowname)) {
      return(narrowest_level)
    }
    lower_levels <- tolower(taxa_levels)
    hits <- taxa_levels[vapply(lower_levels, function(p) {
      grepl(paste0("^", p, "\\."), rowname)
    }, logical(1))]
    if (length(hits) == 0) {
      stop(sprintf("Unrecognized taxonomic parameter '%s' (taxa_levels: %s).",
                   rowname, paste(taxa_levels, collapse = ", ")))
    }
    if (length(hits) > 1) {
      stop(sprintf("Ambiguous taxonomic parameter '%s' matches multiple levels: %s.",
                   rowname, paste(hits, collapse = ", ")))
    }
    hits
  }

  row_levels <- vapply(rownames(results2), resolve_level_for_row, character(1), USE.NAMES = FALSE)
  results2$domain <- row_levels
  results2$taxa_full <- mapply(function(lvl, idx) level_labels[[lvl]][idx],
                               row_levels, results2$taxa_index)
  results2$exposure <- paste0(exposure[results2$Exposure.Index])

  # Modify Exposure variable
  results2 <- results2 %>%
    mutate(exposure=case_when(
      grepl("psi",rownames(.)) ~ "Mixture",
      grepl("disp",rownames(.)) ~ "Dispersion",
      grepl("omega",rownames(.)) ~ "Omega",
      TRUE ~ exposure))

  # Get Taxa name from full taxonomic string
  results2 <- results2 %>%
    mutate(taxa_name = sub(".*__", "", taxa_full))

  # Remove "disp" and "omega" estimates
  if(!return_all_estimates){
    results2 <- results2 %>%
      filter(!grepl("disp",rownames(results2)),
             !grepl("omega",rownames(results2)))
  }

  # Remove rownames
  rownames(results2) <- NULL

  # Select final variables
  results2 <- results2 %>%
    select(taxa_full, taxa_name, domain, exposure,component,
           estimate,bci_lcl,bci_ucl,pdir,prope,pmap)

  return(results2)
}
