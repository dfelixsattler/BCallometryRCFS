test_that("ht_naslund returns correct scalar value", {
  expect_equal(ht_naslund(20, a = 1.5, b = 0.05),
               1.3 + 20^2 / (1.5 + 0.05 * 20)^2)
})

test_that("ht_curtis returns correct scalar value", {
  expect_equal(ht_curtis(20, a = 35, b = 0.9),
               1.3 + 35 * (20 / 21)^0.9)
})

test_that("ht_logistic returns correct scalar value", {
  expect_equal(ht_logistic(20, a = 40, b = 5, c = 0.1),
               1.3 + 40 / (1 + 5 * exp(-0.1 * 20)))
})

test_that("ht_korf returns correct scalar value", {
  expect_equal(ht_korf(20, a = 40, b = 3, c = 0.5),
               1.3 + 40 * exp(-3 * 20^(-0.5)))
})

test_that("ht_weibull returns correct scalar value", {
  expect_equal(ht_weibull(20, a = 40, b = 0.05, c = 1.2),
               1.3 + 40 * (1 - exp(-0.05 * 20^1.2)))
})

test_that("ht_richards returns correct scalar value", {
  expect_equal(ht_richards(20, a = 40, b = 0.05, c = 2),
               1.3 + 40 * (1 - exp(-0.05 * 20))^2)
})

test_that("prediction functions return vectors of correct length", {
  dbh <- c(10, 20, 30, 40)
  expect_length(ht_naslund(dbh, 1.5, 0.05),    4)
  expect_length(ht_curtis(dbh, 35, 0.9),        4)
  expect_length(ht_logistic(dbh, 40, 5, 0.1),   4)
  expect_length(ht_korf(dbh, 40, 3, 0.5),       4)
  expect_length(ht_weibull(dbh, 40, 0.05, 1.2), 4)
  expect_length(ht_richards(dbh, 40, 0.05, 2),  4)
})

test_that("ht_from_dbh dispatches to correct equation", {
  dbh <- c(10, 20, 30)
  expect_equal(ht_from_dbh(dbh, "naslund",  a = 1.5,  b = 0.05),
               ht_naslund(dbh, 1.5, 0.05))
  expect_equal(ht_from_dbh(dbh, "curtis",   a = 35,   b = 0.9),
               ht_curtis(dbh, 35, 0.9))
  expect_equal(ht_from_dbh(dbh, "logistic", a = 40,   b = 5,    c = 0.1),
               ht_logistic(dbh, 40, 5, 0.1))
  expect_equal(ht_from_dbh(dbh, "korf",     a = 40,   b = 3,    c = 0.5),
               ht_korf(dbh, 40, 3, 0.5))
  expect_equal(ht_from_dbh(dbh, "weibull",  a = 40,   b = 0.05, c = 1.2),
               ht_weibull(dbh, 40, 0.05, 1.2))
  expect_equal(ht_from_dbh(dbh, "richards", a = 40,   b = 0.05, c = 2),
               ht_richards(dbh, 40, 0.05, 2))
})

test_that("ht_from_dbh is case-insensitive for model name", {
  dbh <- 20
  expect_equal(ht_from_dbh(dbh, "NASLUND", a = 1.5, b = 0.05),
               ht_naslund(dbh, 1.5, 0.05))
  expect_equal(ht_from_dbh(dbh, "Korf", a = 40, b = 3, c = 0.5),
               ht_korf(dbh, 40, 3, 0.5))
})

test_that("ht_from_dbh errors on unknown model", {
  expect_error(ht_from_dbh(10, "bogus", a = 1, b = 1),
               "Unknown model")
})

test_that("hd_start_values returns named list for all models", {
  set.seed(1)
  dbh <- runif(50, 5, 60)
  ht  <- ht_korf(dbh, 38, 2.5, 0.4) + rnorm(50, 0, 1)

  for (m in c("naslund", "curtis", "logistic", "korf", "weibull", "richards")) {
    sv <- hd_start_values(dbh, ht, model = m)
    expect_type(sv, "list")
    expect_true(length(sv) >= 2)
    expect_true(all(names(sv) %in% c("a", "b", "c")))
  }
})

test_that("fit_hd_model returns correct structure (fixed effects)", {
  set.seed(42)
  dbh    <- runif(200, 5, 60)
  height <- ht_naslund(dbh, a = 1.2, b = 0.04) + rnorm(200, 0, 1.5)
  dat    <- data.frame(DBH = dbh, HEIGHT = height)

  result <- fit_hd_model(dat, model = "naslund")

  expect_named(result, c("model", "fit", "coefficients", "r2_marg"))
  expect_equal(result$model, "naslund")
  expect_s3_class(result$fit, "nls")
  expect_true(result$r2_marg > 0.8)
  expect_named(result$coefficients, c("a", "b"))
})

test_that("fit_hd_model works with custom column names", {
  set.seed(7)
  dbh    <- runif(100, 5, 50)
  height <- ht_naslund(dbh, a = 1.5, b = 0.12) + rnorm(100, 0, 1)
  dat    <- data.frame(diam = dbh, ht = height)

  result <- fit_hd_model(dat, model = "naslund",
                         dbh_col = "diam", height_col = "ht")
  expect_equal(result$model, "naslund")
  expect_named(result$coefficients, c("a", "b"))
  expect_true(result$r2_marg > 0.8)
})

test_that("fit_hd_model errors on missing columns", {
  dat <- data.frame(D = 1:20, H = seq(5, 40, length.out = 20))
  expect_error(fit_hd_model(dat), "Column\\(s\\) not found")
})

test_that("fit_hd_model errors when too few valid rows", {
  dat <- data.frame(DBH = c(NA, NA, 5), HEIGHT = c(10, NA, 2.0))
  expect_error(fit_hd_model(dat), "Fewer than 10 valid")
})

test_that("fit_hd_model errors when nested_col supplied without group_col", {
  set.seed(1)
  dat <- data.frame(DBH = runif(50, 5, 40),
                    HEIGHT = runif(50, 5, 30),
                    TREE = paste0("T", 1:50))
  expect_error(fit_hd_model(dat, nested_col = "TREE"),
               "requires 'group_col'")
})

test_that("fit_hd_model single-level ME returns random_effects_site", {
  skip_if_not_installed("nlme")
  set.seed(42)
  sites  <- rep(paste0("S", 1:5), each = 20)
  dbh    <- runif(100, 5, 50)
  height <- ht_naslund(dbh, a = 1.2, b = 0.04) + rnorm(100, 0, 1.5)
  dat    <- data.frame(DBH = dbh, HEIGHT = height, SITE = sites)
  res    <- fit_hd_model(dat, model = "naslund", group_col = "SITE")
  expect_true("random_effects_site" %in% names(res))
  expect_false("random_effects_tree" %in% names(res))
  expect_equal(nrow(res$random_effects_site), 5L)
  expect_true(all(c("SITE_IDENTIFIER", "a") %in% names(res$random_effects_site)))
})

test_that("fit_hd_model two-level ME returns both random effect tables", {
  skip_if_not_installed("nlme")
  set.seed(7)
  sites  <- rep(paste0("S", 1:5), each = 40)
  trees  <- paste0(sites, "_T", rep(1:10, 20))   # 10 trees per site, 4 obs each
  dbh    <- runif(200, 5, 50)
  height <- ht_naslund(dbh, a = 1.2, b = 0.04) + rnorm(200, 0, 1.5)
  dat    <- data.frame(DBH = dbh, HEIGHT = height, SITE = sites, TREE = trees)
  res    <- fit_hd_model(dat, model = "naslund",
                         group_col = "SITE", nested_col = "TREE")
  expect_true("random_effects_site" %in% names(res))
  expect_true("random_effects_tree" %in% names(res))
  expect_equal(nrow(res$random_effects_site), 5L)
  expect_equal(nrow(res$random_effects_tree), 50L)   # 5 sites * 10 trees
  expect_true(all(c("unitreeid", "a") %in% names(res$random_effects_tree)))
})

# ---- ht_predict() -----------------------------------------------------------

.make_trees <- function() {
  data.frame(
    SPECIES = c("Pl", "Sx", "Hw", "Fd", "Bl", "At"),
    SP0     = c("PL", "S",  "H",  "F",  "B",  "AT"),
    DBH     = c(15,   22,   18,   30,   12,    8),
    stringsAsFactors = FALSE
  )
}

.make_coefs_sp <- function() {
  data.frame(
    SPECIES = c("Pl", "Sx", "Hw"),
    model   = "naslund",
    a = c(1.2, 1.0, 1.5), b = c(0.04, 0.05, 0.03),
    stringsAsFactors = FALSE
  )
}

.make_coefs_sp0 <- function() {
  data.frame(
    SP0   = c("PL", "S",  "H",  "F",  "B",  "AT"),
    model = "naslund",
    a = c(1.2,  1.0,  1.5,  1.1,  1.3,  0.9),
    b = c(0.04, 0.05, 0.03, 0.045, 0.035, 0.06),
    stringsAsFactors = FALSE
  )
}

test_that("ht_predict exact match: all trees matched at first level", {
  trees    <- .make_trees()[1:3, ]      # Pl, Sx, Hw — all in coefs_sp
  coefs_sp <- .make_coefs_sp()
  result   <- ht_predict(trees, models = list(species = coefs_sp))
  expect_true(all(!is.na(result$HT_PREDICTED)))
  expect_true(all(result$HT_MODEL_SOURCE == "species"))
  expect_true(all(result$HT_PREDICTED > 1.3))
})

test_that("ht_predict fallback: unmatched trees cascade to next level", {
  trees    <- .make_trees()             # 6 species; Pl/Sx/Hw in sp, F/B/AT in sp0
  coefs_sp <- .make_coefs_sp()
  coefs_sp0 <- .make_coefs_sp0()
  result   <- ht_predict(trees,
                          models = list(species = coefs_sp, sp0 = coefs_sp0))
  expect_true(all(!is.na(result$HT_PREDICTED)))
  expect_equal(result$HT_MODEL_SOURCE[result$SPECIES %in% c("Pl","Sx","Hw")],
               rep("species", 3))
  expect_equal(result$HT_MODEL_SOURCE[result$SPECIES %in% c("Fd","Bl","At")],
               rep("sp0", 3))
})

test_that("ht_predict source_col records exact level name used", {
  trees <- .make_trees()[4:5, ]         # Fd (SP0=F), Bl (SP0=B)
  coefs_sp0 <- .make_coefs_sp0()
  result <- ht_predict(trees, models = list(sp0_fallback = coefs_sp0))
  expect_true(all(result$HT_MODEL_SOURCE == "sp0_fallback"))
})

test_that("ht_predict warns and returns NA for unmatched trees", {
  trees    <- .make_trees()[4:5, ]      # Fd, Bl — not in species-level table
  coefs_sp <- .make_coefs_sp()          # only Pl/Sx/Hw
  expect_warning(
    result <- ht_predict(trees, models = list(species = coefs_sp)),
    regexp = "could not be matched"
  )
  expect_true(all(is.na(result$HT_PREDICTED)))
  expect_true(all(is.na(result$HT_MODEL_SOURCE)))
})

test_that("ht_predict supports custom column names", {
  trees <- data.frame(d = c(15, 22), SP0 = c("PL", "S"))
  coefs <- data.frame(SP0 = c("PL", "S"), model = "naslund",
                      a = c(1.2, 1.0), b = c(0.04, 0.05))
  result <- ht_predict(trees, models = list(sp0 = coefs),
                       dbh_col = "d", ht_col = "ht", source_col = "src")
  expect_true("ht"  %in% names(result))
  expect_true("src" %in% names(result))
  expect_true(all(!is.na(result$ht)))
})

test_that("ht_predict works with three-parameter (korf) model", {
  trees <- data.frame(SP0 = c("PL", "S"), DBH = c(15, 22))
  coefs <- data.frame(SP0 = c("PL", "S"), model = "korf",
                      a = c(38, 40), b = c(2.5, 2.0), c = c(0.4, 0.5))
  result <- ht_predict(trees, models = list(sp0 = coefs))
  expect_true(all(!is.na(result$HT_PREDICTED)))
  expect_true(all(result$HT_PREDICTED > 1.3))
})

test_that("ht_predict stops with informative error on non-named models", {
  trees <- .make_trees()
  expect_error(ht_predict(trees, models = list(.make_coefs_sp())),
               "named list")
})

test_that("ht_predict warns and skips level with missing group column in data", {
  trees <- data.frame(DBH = 15, SP0 = "PL")   # no SPECIES col; SP0 "PL" is in coefs_sp0
  coefs_sp  <- .make_coefs_sp()                # requires SPECIES — should be skipped
  coefs_sp0 <- .make_coefs_sp0()               # requires SP0 — should match
  expect_warning(
    result <- ht_predict(trees,
                          models = list(species = coefs_sp, sp0 = coefs_sp0)),
    regexp = "not found in data"
  )
  expect_false(is.na(result$HT_PREDICTED))
  expect_equal(result$HT_MODEL_SOURCE, "sp0")
})

# ---- fit_hd_models_by_group() -----------------------------------------------

.make_group_dat <- function(n = 300, seed = 1) {
  set.seed(seed)
  sp0  <- sample(c("F", "PL", "H"), n, replace = TRUE)
  spt  <- ifelse(sp0 == "H", "D", "C")
  dbh  <- runif(n, 5, 50)
  ht   <- ht_naslund(dbh, a = 1.2, b = 0.04) + rnorm(n, 0, 2)
  site <- sample(paste0("S", 1:10), n, replace = TRUE)
  data.frame(SPECIES_SP0 = sp0, SPECIES_SP_TYPE = spt,
             SITE_IDENTIFIER = site, unitreeid = paste0(site, "_T", seq_len(n)),
             DBH = dbh, HEIGHT = ht, stringsAsFactors = FALSE)
}

test_that("fit_hd_models_by_group returns fits list with correct structure", {
  dat <- .make_group_dat()
  res <- fit_hd_models_by_group(dat)
  expect_type(res, "list")
  expect_named(res, c("fits", "summary"))
  # All SP0 groups present
  sp0_in_data  <- sort(unique(dat$SPECIES_SP0))
  sp0_in_fits  <- names(res$fits)[vapply(res$fits, function(m) m$level == "sp0", logical(1))]
  expect_setequal(sp0_in_fits, sp0_in_data)
})

test_that("fit_hd_models_by_group fit elements have expected fields", {
  dat <- .make_group_dat()
  res <- fit_hd_models_by_group(dat)
  for (m in res$fits) {
    expect_true(all(c("fit", "method", "level", "group", "n") %in% names(m)))
    expect_true(m$level %in% c("sp0", "sp_type"))
  }
})

test_that("fit_hd_models_by_group fits succeed when data are adequate", {
  dat <- .make_group_dat()
  res <- fit_hd_models_by_group(dat)
  sp0_fits <- res$fits[vapply(res$fits, function(m) m$level == "sp0", logical(1))]
  expect_true(all(vapply(sp0_fits, function(m) m$method != "none", logical(1))))
})

test_that("fit_hd_models_by_group SP_TYPE fallback triggers for sparse SP0", {
  dat  <- .make_group_dat()
  # Keep only 5 obs of SP0 "F" so it falls below min_n
  f_rows <- which(dat$SPECIES_SP0 == "F")
  dat2   <- dat[-f_rows[6:length(f_rows)], ]  # keep only first 5 F rows
  res    <- fit_hd_models_by_group(dat2, min_n = 10L)
  # "F" should have method "none"; its SP_TYPE "C" should appear as fallback
  f_fit  <- res$fits[["F"]]
  expect_equal(f_fit$method, "none")
  expect_true("C" %in% names(res$fits))
  expect_equal(res$fits[["C"]]$level, "sp_type")
})

test_that("fit_hd_models_by_group summary has correct columns and rows", {
  dat <- .make_group_dat()
  res <- fit_hd_models_by_group(dat)
  expect_true(all(c("Level", "Group", "Method", "n", "Marg_R2") %in%
                  names(res$summary)))
  expect_equal(nrow(res$summary), length(res$fits))
})

test_that("fit_hd_models_by_group errors on missing required column", {
  dat <- .make_group_dat()
  dat$SPECIES_SP_TYPE <- NULL
  expect_error(fit_hd_models_by_group(dat), "Column\\(s\\) not found")
})
