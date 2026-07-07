test_that("biomass_tree returns correct structure", {
  result <- biomass_tree(
    species = c("lodgepole pine", "white spruce", "trembling aspen"),
    dbh     = c(15, 22, 18)
  )
  expect_length(result, 3)
  expect_true(all(result > 0))
  expect_true(all(is.finite(result)))
})

test_that("biomass_tree DBH-only increases with DBH", {
  dbh <- c(10, 20, 30, 40)
  b   <- biomass_tree("lodgepole pine", dbh)
  expect_true(all(diff(b) > 0))
})

test_that("biomass_tree height-included > DBH-only for most species", {
  # Adding height generally changes the estimate
  b_dbh <- biomass_tree("douglas-fir", 30)
  b_ht  <- biomass_tree("douglas-fir", 30, height = 28)
  expect_false(isTRUE(all.equal(b_dbh, b_ht)))
})

test_that("biomass_tree Lambert2005 vs Ung2008 differ for lodgepole pine", {
  b_L <- biomass_tree("lodgepole pine", 25, paper_source = "Lambert2005")
  b_U <- biomass_tree("lodgepole pine", 25, paper_source = "Ung2008")
  expect_false(isTRUE(all.equal(b_L, b_U)))
})

test_that("biomass_tree is case-insensitive for species", {
  b1 <- biomass_tree("Lodgepole Pine", 20)
  b2 <- biomass_tree("lodgepole pine", 20)
  expect_equal(b1, b2)
})

test_that("biomass_tree returns NA with warning for unknown species (DBH-only)", {
  expect_warning(
    result <- biomass_tree("unicorn tree", 20),
    "No coefficients"
  )
  expect_true(is.na(result))
})

test_that("biomass_tree uses generic with warning for unknown species (height)", {
  expect_warning(
    result <- biomass_tree("unicorn tree", 20, height = 15),
    "not recognised"
  )
  expect_true(is.finite(result) && result > 0)
})

test_that("biomass_tree errors when lengths differ", {
  # scalar species is fine (recycled), but two species vs four dbh must error
  expect_error(
    biomass_tree(c("lodgepole pine", "white spruce"), c(20, 25, 30, 35)),
    "same length"
  )
})

test_that("biomass_components returns correct structure", {
  result <- biomass_components(
    species = c("douglas-fir", "western hemlock"),
    dbh     = c(30, 25),
    height  = c(28, 22)
  )
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 2)
  expect_named(result, c("wood", "bark", "branches", "foliage", "total"))
  expect_equal(result$total, result$wood + result$bark +
                             result$branches + result$foliage)
})

test_that("biomass_components total matches biomass_tree", {
  species <- c("lodgepole pine", "trembling aspen", "western redcedar")
  dbh     <- c(15, 20, 35)
  height  <- c(12, 16, 30)

  totals_tree <- biomass_tree(species, dbh, height)
  totals_comp <- biomass_components(species, dbh, height)$total

  expect_equal(totals_tree, totals_comp)
})

test_that("biomass_coefs dataset is accessible and well-formed", {
  expect_s3_class(biomass_coefs, "data.frame")
  expect_true(all(c("species", "component", "paper_source",
                    "height_included", "a1", "a2", "a3") %in%
                    names(biomass_coefs)))
  expect_true(all(biomass_coefs$component %in%
                    c("wood", "bark", "branches", "foliage")))
  expect_true(all(biomass_coefs$paper_source %in%
                    c("Lambert2005", "Ung2008")))
  # a1 and a2 always positive
  expect_true(all(biomass_coefs$a1 > 0))
  expect_true(all(biomass_coefs$a2 > 0))
})
