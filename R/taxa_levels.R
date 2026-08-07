#' Default taxonomic hierarchy used by BaHZING
#'
#' Ordered broadest to narrowest. `Format_BaHZING()` and `BaHZING_Model()`
#' both default to this when no custom `taxa_levels` is supplied, reproducing
#' BaHZING's original fixed 6-level hierarchy exactly. `"Kingdom"` is
#' deliberately not part of this vector: it is a separate, always-present
#' rank handled as a hand-written special case wherever taxonomy prefixes are
#' assigned, since every taxon has a kingdom but hierarchy depth below it can
#' vary by dataset.
#'
#' @keywords internal
#' @noRd
default_taxa_levels <- c("Phylum", "Class", "Order", "Family", "Genus", "Species")

#' Validate a taxa_levels vector
#'
#' Checks structural invariants that both `Format_BaHZING()` and
#' `BaHZING_Model()` rely on: at least two levels (hierarchical shrinkage
#' needs a parent-child pair to borrow strength across), no missing/duplicate
#' names, names usable as R identifiers (they get pasted directly into
#' generated JAGS variable names and data.frame column names), and no
#' first-letter collision either among the levels themselves or with the
#' fixed `"Kingdom"` rank. First letters double as both the taxonomy prefix
#' codes (`p__`, `c__`, ...) and the generated JAGS variable-name prefixes, so
#' a collision would silently merge two levels' data instead of erroring.
#'
#' @param taxa_levels Character vector, broadest to narrowest.
#' @return Invisibly `TRUE` if valid; otherwise `stop()`s with a description
#'   of the violated invariant.
#' @keywords internal
#' @noRd
validate_taxa_levels <- function(taxa_levels) {
  if (!is.character(taxa_levels) || anyNA(taxa_levels)) {
    stop("taxa_levels must be a character vector with no NAs.")
  }
  if (length(taxa_levels) < 2) {
    stop("taxa_levels must have at least 2 levels for hierarchical shrinkage.")
  }
  if (any(duplicated(taxa_levels))) {
    stop("taxa_levels must not contain duplicate level names.")
  }
  if (any(taxa_levels != make.names(taxa_levels))) {
    stop("taxa_levels entries must be valid R identifiers (used to build data/variable names).")
  }
  all_names <- c("Kingdom", taxa_levels)
  first_letters <- tolower(substr(all_names, 1, 1))
  if (any(duplicated(first_letters))) {
    stop("taxa_levels (plus the fixed 'Kingdom' rank) must have distinct first letters: ",
         "first letters drive both the taxonomy prefix codes (p__, c__, ...) and the ",
         "generated JAGS variable prefixes, and colliding first letters would silently ",
         "merge two levels.")
  }
  invisible(TRUE)
}

#' Look up a taxonomic level name by position, with a bounds check
#'
#' @param taxa_levels Character vector, broadest to narrowest.
#' @param i Integer position to look up.
#' @return The level name at position `i`.
#' @keywords internal
#' @noRd
taxa_level_name <- function(taxa_levels, i) {
  if (length(i) != 1 || is.na(i) || i < 1 || i > length(taxa_levels)) {
    stop(sprintf("Invalid taxonomic level index %s: taxa_levels has %d levels (%s).",
                 i, length(taxa_levels), paste(taxa_levels, collapse = " > ")))
  }
  taxa_levels[[i]]
}

#' Name of the binary incidence matrix linking adjacent levels i and i+1
#'
#' Single owner of the `"<narrower>.<broader>.Matrix"` naming convention
#' (e.g. `"Genus.Family.Matrix"`), so `Format_BaHZING()` (which writes these
#' matrices) and `BaHZING_Model()` (which reads them back) can never drift
#' out of sync on the name.
#'
#' @param taxa_levels Character vector, broadest to narrowest.
#' @param i Integer position of the broader level in the pair; the narrower
#'   level is `i + 1`.
#' @return A length-1 character string, the matrix's key name.
#' @keywords internal
#' @noRd
hierarchy_matrix_name <- function(taxa_levels, i) {
  broader <- taxa_level_name(taxa_levels, i)
  narrower <- taxa_level_name(taxa_levels, i + 1)
  paste0(narrower, ".", broader, ".Matrix")
}
