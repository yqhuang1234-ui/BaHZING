test_that("default_taxa_levels is valid and matches BaHZING's original 6-level scheme", {
  expect_equal(default_taxa_levels, c("Phylum", "Class", "Order", "Family", "Genus", "Species"))
  expect_true(validate_taxa_levels(default_taxa_levels))
})

test_that("validate_taxa_levels rejects structural violations", {
  expect_error(validate_taxa_levels(c("Only.One")), "at least 2 levels")
  expect_error(validate_taxa_levels(c("Genus", NA)), "no NAs")
  expect_error(validate_taxa_levels(c(1, 2)), "character vector")
  expect_error(validate_taxa_levels(c("Genus", "Genus")), "duplicate")
  expect_error(validate_taxa_levels(c("Genus", "Not Valid")), "valid R identifiers")
  # "Kingdom" (fixed) and "Klass" collide on first letter "k"
  expect_error(validate_taxa_levels(c("Klass", "Genus")), "distinct first letters")
  # Two custom levels colliding with each other
  expect_error(validate_taxa_levels(c("Genus", "Group")), "distinct first letters")
})

test_that("validate_taxa_levels accepts a custom differently-sized hierarchy", {
  expect_true(validate_taxa_levels(c("Phylum", "Family", "Genus", "ASV")))
})

test_that("taxa_level_name looks up by position and bounds-checks", {
  expect_equal(taxa_level_name(default_taxa_levels, 1), "Phylum")
  expect_equal(taxa_level_name(default_taxa_levels, 6), "Species")
  expect_error(taxa_level_name(default_taxa_levels, 0), "Invalid taxonomic level index")
  expect_error(taxa_level_name(default_taxa_levels, 7), "Invalid taxonomic level index")
  expect_error(taxa_level_name(default_taxa_levels, NA), "Invalid taxonomic level index")
})

test_that("hierarchy_matrix_name matches the legacy naming convention", {
  expect_equal(hierarchy_matrix_name(default_taxa_levels, 1), "Class.Phylum.Matrix")
  expect_equal(hierarchy_matrix_name(default_taxa_levels, 2), "Order.Class.Matrix")
  expect_equal(hierarchy_matrix_name(default_taxa_levels, 3), "Family.Order.Matrix")
  expect_equal(hierarchy_matrix_name(default_taxa_levels, 4), "Genus.Family.Matrix")
  expect_equal(hierarchy_matrix_name(default_taxa_levels, 5), "Species.Genus.Matrix")
})
