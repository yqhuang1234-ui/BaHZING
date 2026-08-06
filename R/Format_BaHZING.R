#' Format_BaHZING Function
#'
#' This function takes a phyloseq object and performs formatting operations on it,
#' including modifying the taxonomic table, uniting taxonomic levels, and creating
#' matrices based on taxonomic information.
#' @details
#' The Format_BaHZING function is the core function of the Format_BaHZING package.
#' It takes a phyloseq object as input and performs various formatting operations
#' to prepare the data for analysis. The function modifies the taxonomic table to add
#' taxonomic prefixes (e.g., "d__" for Kingdom), unites taxonomic levels, and creates
#' matrices based on taxonomic information. The formatted data is then returned as a list
#' containing different data frames for further analysis.
#'
#' The package relies on the `phyloseq`, `dplyr`, and `stringr` packages for data manipulation.
#'
#' The main function `Format_BaHZING` is exported and can be accessed by other packages or scripts
#' that depend on the functionalities provided by this package.
#'
#'
#' @param phyloseq.object A phyloseq object.
#' @param taxa_levels Character vector naming the taxonomic hierarchy to use,
#' ordered broadest to narrowest (e.g. `c("Phylum","Class","Order","Family","Genus","Species")`,
#' the default). These must be column names in `phyloseq.object`'s `tax_table`
#' (aside from the narrowest level, which is created automatically if absent).
#' Override this for a dataset with a different taxonomic depth or naming
#' convention (e.g. a 4-level 16S scheme) - the hierarchical matrices and
#' downstream model in `BaHZING_Model()` adapt to whatever is supplied here.
#' The default reproduces BaHZING's original fixed 6-level hierarchy exactly.
#' @return A list with the following elements:
#'   - `Table`: Formatted microbiome data as a data frame.
#'   - `taxa_levels`: the taxonomic hierarchy used (echoes the `taxa_levels` argument).
#'   - One binary incidence matrix per adjacent pair of levels, named
#'     `"<narrower>.<broader>.Matrix"` (e.g. `Species.Genus.Matrix`,
#'     `Genus.Family.Matrix`, ... for the default hierarchy).
#' @details The column names 'Kingdom' (or 'Domain') plus whatever is passed
#' in `taxa_levels` should be present in the tax_table of the phyloseq object.
#'
#' @import phyloseq
#' @import dplyr
#' @import stringr
#' @importFrom utils globalVariables
#' @importFrom phyloseq sample_data otu_table tax_table
#' @importFrom dplyr full_join select mutate %>%
#' @export
#' @name Format_BaHZING

# Declare global variables
utils::globalVariables(c("Domain"))

Format_BaHZING <- function(phyloseq.object, taxa_levels = default_taxa_levels) {
  validate_taxa_levels(taxa_levels)

  Kingdom <- NULL    # Variable to store Kingdom taxonomic level

  # Check if taxa are stored as rows in the phyloseq object
  if (phyloseq::taxa_are_rows(phyloseq.object)) {
    # If taxa are stored as rows, transpose the phyloseq object
    # to convert the taxa into columns, making further processing easier
    phyloseq.object <- t(phyloseq.object)
  }

  # Downstream (colnames(table) <- ASV.names) renames otu.table's columns
  # purely by POSITION, assuming otu.table column i and taxa.table row i
  # describe the same taxon - true because both are extracted from the same
  # phyloseq.object, and phyloseq keeps every component of one object
  # aligned to that object's own canonical taxon order. That alignment is
  # never itself re-verified though, so check it explicitly here, on the
  # phyloseq accessors directly (BEFORE data.frame() below, whose default
  # check.names=TRUE sanitizes column names - turning ";"/" " into "." -
  # but leaves row names untouched, which would make an otherwise-correctly
  # -aligned pair look mismatched if compared after coercion): a malformed
  # or unusually-constructed phyloseq object would otherwise silently
  # mislabel every taxon instead of erroring.
  stopifnot(identical(colnames(phyloseq::otu_table(phyloseq.object)),
                      rownames(phyloseq::tax_table(phyloseq.object))))

  # Extract metadata, OTU table, and taxonomic table from the phyloseq object
  meta.data <- data.frame(phyloseq::sample_data(phyloseq.object))  # Data frame containing metadata information
  otu.table <- data.frame(phyloseq::otu_table(phyloseq.object))    # Data frame containing OTU (Operational Taxonomic Unit) table
  taxa.table <- data.frame(phyloseq::tax_table(phyloseq.object))  # Data frame containing taxonomic table

  # Check if the taxonomic table contains more than one taxonomic level
  if (length(colnames(taxa.table)) < 2) {
    # If the table has less than two taxonomic levels, raise an error and stop the function
    stop("Need > 1 taxonomic level")
  }

  # Add taxonomic prefixes to the Kingdom/Domain level (hand-written special
  # case: every taxon has exactly one kingdom, unlike the rest of the
  # hierarchy below it, which can vary in depth/naming by dataset).
  if ("Kingdom" %in% colnames(taxa.table)) {
    # If 'Kingdom' column is present in the taxonomic table, add 'k__' prefix.
    # Strip any existing "<letter>__" prefix first (e.g. "d__" left over from
    # a SILVA/QIIME2-style source that wasn't also named "Domain") - without
    # this, a value already carrying a different prefix would get "k__"
    # stacked on top instead of replaced (e.g. "k__d__Bacteria").
    taxa.table <- taxa.table %>%
      mutate(Kingdom=sub("^[A-Za-z]__", "", Kingdom)) %>%
      mutate(Kingdom=ifelse(grepl("k__",Kingdom),Kingdom,paste0("k__", Kingdom)))
  }

  if ("Domain" %in% colnames(taxa.table)) {
    # If 'Kingdom' column is present in the taxonomic table, add 'k__' prefix
    taxa.table <- taxa.table %>%
      rename(Kingdom=Domain)%>%
      mutate(Kingdom=case_when(grepl("d__",Kingdom) ~ str_replace(Kingdom,"d__","k__"),
                               grepl("k__",Kingdom) ~ Kingdom,
                               !grepl("d__",Kingdom) & !grepl("k__",Kingdom) ~ paste0("k__", Kingdom)))
  }

  #If missing kingdom level information, stop and produce error message
  if ("k__NA" %in% taxa.table$Kingdom) {
    stop("Missing Kindom-level information for an ASV/OTU. Remove unidentified bacteria and rerun.")
  }

  narrowest_level <- taxa_levels[[length(taxa_levels)]]

  # Add taxonomic prefixes to every level in taxa_levels, e.g. "p__" for
  # Phylum, "g__" for Genus - derived from each level's own first letter, so
  # this generalizes to custom level names for free. Any existing
  # "<letter>__" prefix is stripped first - values coming from a real
  # classifier pipeline (QIIME2/GTDB-Tk/etc.) commonly already carry the
  # standard p__/c__/o__/f__/g__/s__ convention baked in, which only happens
  # to match what gets computed here when taxa_levels uses the standard
  # level names. With custom names (e.g. "Alpha" -> "a__"), not stripping
  # first would stack the new prefix on top of the old one instead of
  # replacing it (e.g. "a__p__Firmicutes_A" instead of "a__Firmicutes_A").
  for (lvl in taxa_levels) {
    if (lvl %in% colnames(taxa.table)) {
      prefix <- paste0(tolower(substr(lvl, 1, 1)), "__")
      col <- taxa.table[[lvl]]
      col <- sub("^[A-Za-z]__", "", col)
      col <- ifelse(grepl(prefix, col), col, paste0(prefix, col))
      col <- ifelse(grepl(paste0(prefix, "NA"), col), NA, col)
      taxa.table[[lvl]] <- col
    }
  }

  #If narrowest level not present, create it.
  if (!(narrowest_level %in% colnames(taxa.table))) {
    prefix <- paste0(tolower(substr(narrowest_level, 1, 1)), "__")
    taxa.table[[narrowest_level]] <- paste0("unclassified", seq_len(nrow(taxa.table)))
    taxa.table[[narrowest_level]] <- paste0(prefix, taxa.table[[narrowest_level]])
    taxa.table[[narrowest_level]] <- ifelse(grepl(paste0(prefix, "NA"), taxa.table[[narrowest_level]]),
                                             NA, taxa.table[[narrowest_level]])
  }

  #Fill in any NAs with unclassified, and remember how many were filled per
  #level (needed below to keep the narrowest level's duplicate-name
  #numbering from colliding with its "unclassifiedN" numbering).
  fill_counts <- setNames(integer(length(taxa_levels)), taxa_levels)
  for (lvl in taxa_levels) {
    if (lvl %in% colnames(taxa.table)) {
      prefix <- paste0(tolower(substr(lvl, 1, 1)), "__")
      fill_idx <- sort(unique(which(!grepl(prefix, taxa.table[[lvl]]))))
      fill_counts[lvl] <- length(fill_idx)
      if (fill_counts[lvl] > 0) {
        taxa.table[[lvl]][fill_idx] <- paste0(prefix, "unclassified", seq_len(fill_counts[lvl]))
      }
    }
  }

  #Real (classified) names can collide at the narrowest level's resolution -
  #broader levels already got unique "unclassifiedN" labels above, but two
  #distinct narrowest-level taxa can otherwise share the same classified name.
  if (narrowest_level %in% colnames(taxa.table)) {
    dup_idx <- which(duplicated(taxa.table[[narrowest_level]]))
    if (length(dup_idx) > 0) {
      low <- fill_counts[[narrowest_level]] + 1
      high <- fill_counts[[narrowest_level]] + length(dup_idx)
      duplicate.names <- paste(taxa.table[[narrowest_level]][dup_idx], "_", low:high)
      duplicate.names <- str_replace_all(duplicate.names, " ", "")
      taxa.table[[narrowest_level]][dup_idx] <- duplicate.names
    }
  }

  #Create a name vector for all taxa levels in taxa table
  taxa.names <- colnames(taxa.table)
  # Concatenate every level's lineage (Kingdom + all levels up to and
  # including itself) into that level's own column, e.g. Genus becomes
  # "k__Bacteria_p__..._g__Blautia". Must proceed narrowest -> broadest: each
  # step overwrites its own target column in place, and only reads pristine
  # data because narrower levels are finalized before broader ranges (which
  # still include them) are built - the range shrinks by exactly the column
  # just written, every step.
  kingdom_col <- if ("Kingdom" %in% colnames(taxa.table)) "Kingdom" else character(0)
  n_levels <- length(taxa_levels)
  for (i in n_levels:1) {
    lvl <- taxa_levels[[i]]
    lineage_cols <- c(kingdom_col, taxa_levels[seq_len(i)])
    taxa.table[[lvl]] <- apply(taxa.table[lineage_cols], 1, paste, collapse = "_")
  }

  ASV.names <- taxa.table[[narrowest_level]]

  # Rename the columns of the otu.table with the values from ASV.names
  table <- otu.table
  colnames(table) <- ASV.names

  # Add an 'id' column to the table and set its values to row names of the table
  table$id <- rownames(table)
  rownames(table) <- NULL

  # Add an 'id' column to the meta.data and set its values to row names of the meta.data
  meta.data$id <- rownames(meta.data)
  rownames(meta.data) <- NULL

  # Perform a full join between meta.data and table using the 'id' column as the key
  table <- dplyr::full_join(meta.data, table, by="id")

  # Set the row names of the table to be the 'id' column values
  rownames(table) <- table$id

  # Remove the 'id' column from the table
  table <- table %>%
    select(-id)

  # Build a binary incidence matrix for every adjacent pair of levels (e.g.
  # Genus x Species, Family x Genus, ...), keyed by hierarchy_matrix_name()
  # so BaHZING_Model() can look each one up generically by the same
  # convention this uses to write it.
  hierarchy_matrices <- list()
  for (i in seq_len(length(taxa_levels) - 1)) {
    broader <- taxa_levels[[i]]
    narrower <- taxa_levels[[i + 1]]
    if (narrower %in% colnames(taxa.table) && broader %in% colnames(taxa.table)) {
      unique_broader <- unique(taxa.table[[broader]])
      unique_narrower <- unique(taxa.table[[narrower]])

      m <- matrix(0, nrow = length(unique_broader), ncol = length(unique_narrower),
                  dimnames = list(unique_broader, unique_narrower))
      m[cbind(match(taxa.table[[broader]], unique_broader),
              match(taxa.table[[narrower]], unique_narrower))] <- 1

      hierarchy_matrices[[hierarchy_matrix_name(taxa_levels, i)]] <- m
    }
  }

  # Create a list named 'Object' to store the formatted table, the
  # taxa_levels used to build it, every hierarchy matrix, and the exact set
  # of Table's columns that are taxon-count (outcome) data - so
  # BaHZING_Model() can select them directly instead of having to guess via
  # string-matching against an implementation detail (e.g. assuming every
  # outcome column contains "k__", which silently breaks if Kingdom/Domain
  # isn't part of the input taxonomy at all).
  #
  # make.names(..., unique = TRUE), not the raw ASV.names: BaHZING_Model()
  # (and any other caller) reconstructs Table via data.frame(formatted_data$Table),
  # whose default check.names=TRUE silently sanitizes non-syntactic column
  # names (e.g. a literal space in a species name becomes "."). table's own
  # colnames were set directly via colnames(table) <- ASV.names, which
  # bypasses that sanitization - so the raw ASV.names would otherwise not
  # match what the columns are actually named by the time a consumer reads
  # them back out. Applying the same sanitization here keeps this list
  # consistent with reality instead of with an intermediate representation
  # nothing downstream actually uses.
  Object <- list()
  Object[["Table"]] <- list(table)
  Object[["taxa_levels"]] <- taxa_levels
  Object[["taxon_columns"]] <- make.names(ASV.names, unique = TRUE)
  for (nm in names(hierarchy_matrices)) {
    Object[[nm]] <- hierarchy_matrices[[nm]]
  }

  # Return the 'Object' list containing the matrices representing relationships between taxonomic levels
  return(Object)
}
