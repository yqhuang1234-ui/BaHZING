test_that("test Format_BaHZING", {
  # Load data
  data("iHMP_Reduced")
  # Format microbiome data
  formatted_data <- Format_BaHZING(iHMP_Reduced)
  # Test format data
  testthat::expect_equal(object = length(formatted_data), expected = 7)
  testthat::expect_true(is.list(formatted_data))
  testthat::expect_true("Table" %in% names(formatted_data))
  testthat::expect_equal(formatted_data$taxa_levels, default_taxa_levels)
})



test_that("Error thrown when taxonomic levels are less than 2", {

  data("iHMP_Reduced")

  PS <- iHMP_Reduced
  # Modify the taxonomic table to have only one level
  tax_table(PS) <- tax_table(tax_table(PS)[,1])

  # Expect that the function stops with an error message
  testthat::expect_error(Format_BaHZING(PS), "Need > 1 taxonomic level")
  })


test_that("Format_BaHZING generalizes to a custom, differently-sized taxa_levels", {
  data("iHMP_Reduced")

  custom_levels <- c("Phylum", "Family", "Genus", "Species")
  formatted_data <- Format_BaHZING(iHMP_Reduced, taxa_levels = custom_levels)

  # 1 Table + 1 taxa_levels + 3 adjacent-pair matrices for a 4-level hierarchy
  testthat::expect_equal(length(formatted_data), 5)
  testthat::expect_equal(formatted_data$taxa_levels, custom_levels)
  testthat::expect_true(all(c("Family.Phylum.Matrix", "Genus.Family.Matrix",
                              "Species.Genus.Matrix") %in% names(formatted_data)))
  # No Class/Order levels in this mapping, so no matrix should reference them
  testthat::expect_false(any(grepl("Class|Order", names(formatted_data))))
})

test_that("If species level not present, create species column", {

  data("iHMP_Reduced")

  PS <- iHMP_Reduced
  # phyloseq object without a 'Species' column
  tax_table(PS) <- tax_table(PS)[ , !colnames(tax_table(PS)) %in% "Species"]

  formatted_data <- Format_BaHZING(PS)

  # Extract the taxonomic table from the result
  taxa_table_result <- formatted_data$Table[[1]]

  # Check if 'Species' column is created
  testthat::expect_true(grepl("s__unclassified",
                              colnames(taxa_table_result)[ncol(taxa_table_result)]))
  testthat::expect_true(grepl("s__unclassified",
                              colnames(taxa_table_result)[ncol(taxa_table_result)-10]))
  testthat::expect_true(grepl("s__unclassified",
                              colnames(taxa_table_result)[ncol(taxa_table_result)-20]))
})

