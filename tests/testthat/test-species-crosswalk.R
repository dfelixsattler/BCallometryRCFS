# Tests for species crosswalk functions and lookup table consistency.
#
# Four categories:
#   A. Lookup table integrity  — every row in bc_species_codes and plant_codes
#      resolves without NA through bc_species_to_sp0() and sp0_to_plant_code()
#   B. bc_psp_to_biomass_name — the three SP0-ambiguous overrides (SE, SS, BL)
#      and a set of common PSP codes that must resolve correctly
#   C. BEC-dependent disambiguation via species_correction() feeding into the
#      full crosswalk chain
#   D. Direct known-pair unit tests for each exported crosswalk function,
#      including edge-case NA / unknown-code behaviour

# ---- A. Lookup table integrity -----------------------------------------------

test_that("every species_code in bc_species_codes resolves to a non-NA SP0", {
  codes  <- bc_species_codes$species_code
  result <- suppressWarnings(bc_species_to_sp0(codes))
  failed <- codes[is.na(result)]
  expect_equal(
    length(failed), 0L,
    info = paste("Codes with no SP0:", paste(failed, collapse = ", "))
  )
})

test_that("every SP0 in plant_codes (taper species) resolves via bc_species_to_plant_code", {
  # Use the detailed codes that anchor each SP0 group as representatives
  detail_codes <- c("FD", "BL", "CW", "DR", "EP", "HW", "LW", "MB",
                    "PL", "PW", "PY", "SW", "AT")
  result <- suppressWarnings(bc_species_to_osm(detail_codes))
  failed <- detail_codes[is.na(result)]
  expect_equal(
    length(failed), 0L,
    info = paste("Codes with no plant_code:", paste(failed, collapse = ", "))
  )
})

test_that("bc_species_codes SP0 column contains only values present in plant_codes or NA", {
  # All non-NA SP0 codes from bc_species_codes must appear somewhere in
  # plant_codes (as the taper-equation species anchor).
  sp0_in_pc    <- unique(plant_codes$sp0[!is.na(plant_codes$sp0)])
  sp0_in_codes <- unique(bc_species_codes$sp0)
  orphaned     <- setdiff(sp0_in_codes, sp0_in_pc)
  expect_equal(
    length(orphaned), 0L,
    info = paste("SP0 codes in bc_species_codes with no plant_codes entry:",
                 paste(orphaned, collapse = ", "))
  )
})

# ---- B. bc_psp_to_biomass_name — override cases and common PSP codes ---------

test_that("bc_psp_to_biomass_name resolves the three SP0-ambiguous overrides", {
  # SE: shares SP0 "S" with white spruce (primary) — must give engelmann spruce
  expect_equal(bc_species_to_biomass_name("SE", "ESSF"), "engelmann spruce")
  # SS: shares SP0 "S" with white spruce (primary) — must give sitka spruce
  expect_equal(bc_species_to_biomass_name("SS", "CWH"),  "sitka spruce")
  # BL: shares SP0 "B" with pacific silver fir (primary) — must give alpine fir
  expect_equal(bc_species_to_biomass_name("BL", "SBS"),  "alpine fir")
})

test_that("bc_psp_to_biomass_name resolves common BC PSP codes correctly", {
  species  <- c("PLI",  "SW",           "FDI",         "HW",
                "AT",             "DR",       "CW",             "EP")
  bec      <- c("SBS",  "SBS",          "IDF",         "CWH",
                "SBS",            "CWH",      "CWH",            "SBS")
  expected <- c("lodgepole pine", "white spruce", "douglas-fir",
                "western hemlock", "trembling aspen", "red alder",
                "western redcedar", "white birch")
  expect_equal(bc_species_to_biomass_name(species, bec), expected)
})

test_that("bc_psp_to_biomass_name vectorises correctly over bec_zone", {
  # Scalar bec_zone should recycle
  result <- bc_species_to_biomass_name(c("PLI", "SW", "AT"), bec_zone = "SBS")
  expect_equal(result, c("lodgepole pine", "white spruce", "trembling aspen"))
})

test_that("bc_psp_to_biomass_name returns NA with warning for unknown codes", {
  expect_warning(
    result <- bc_species_to_biomass_name("ZZZ", "SBS"),
    regexp = "no biomass name"
  )
  expect_true(is.na(result))
})

# ---- C. BEC-dependent disambiguation via species_correction -----------------

test_that("generic 'S' resolves to the correct spruce by BEC zone", {
  # Interior non-ESSF (SBS) → SW → white spruce
  expect_equal(bc_species_to_biomass_name("S", "SBS"),  "white spruce")
  # ESSF → SE → engelmann spruce
  expect_equal(bc_species_to_biomass_name("S", "ESSF"), "engelmann spruce")
  # Coastal CWH → SS → sitka spruce
  expect_equal(bc_species_to_biomass_name("S", "CWH"),  "sitka spruce")
})

test_that("generic 'B' resolves to correct fir by BEC zone", {
  # Interior (SBS) → BL → alpine fir
  expect_equal(bc_species_to_biomass_name("B", "SBS"),  "alpine fir")
  # Coastal (CWH, non-mm subzone) → BA → pacific silver fir
  expect_equal(bc_species_to_biomass_name("B", "CWH"),  "pacific silver fir")
})

test_that("generic 'H' resolves to correct hemlock by BEC zone", {
  # Coastal or SBS → HW → western hemlock
  expect_equal(bc_species_to_biomass_name("H", "CWH"),  "western hemlock")
  expect_equal(bc_species_to_biomass_name("H", "SBS"),  "western hemlock")
  # MH / ESSF → HM → western hemlock (HM has no separate biomass eq.; falls back)
  expect_equal(bc_species_to_biomass_name("H", "MH"),   "western hemlock")
})

test_that("species_correction + bc_species_to_sp0 round-trip for PLI and FDI", {
  # PLI (interior lodgepole): correction → PL; SP0 → PL
  corr <- species_correction("PLI", "SBS")
  expect_equal(corr, "PL")
  expect_equal(bc_species_to_sp0(corr), "PL")

  # FDI (interior Douglas-fir): correction → FD; SP0 → F
  corr2 <- species_correction("FDI", "IDF")
  expect_equal(corr2, "FD")
  expect_equal(bc_species_to_sp0(corr2), "F")
})

# ---- D. Direct known-pair unit tests for each exported function --------------

# -- species_correction -------------------------------------------------------

test_that("species_correction applies unconditional remaps correctly", {
  # The documented example (from ?species_correction)
  species_raw <- c("FDI", "FDC", "HXM", "B",   "S",   "SXE", "PLI", "L")
  bec         <- c("CWH", "IDF", "MH",  "SBS", "CWH", "ESSF","SBS", "ICH")
  expected    <- c("FD",  "FD",  "HM",  "BL",  "SS",  "SE",  "PL",  "LW")
  expect_equal(species_correction(species_raw, bec), expected)
})

test_that("species_correction handles scalar bec_zone recycling", {
  # All trees in SBS: PLI → PL, S → SW, B → BL
  expect_equal(
    species_correction(c("PLI", "S", "B"), bec_zone = "SBS"),
    c("PL", "SW", "BL")
  )
})

test_that("species_correction passes through already-correct codes unchanged", {
  # Codes that need no correction stay as-is
  expect_equal(species_correction(c("HW", "SW", "BL", "AT"), "SBS"),
               c("HW", "SW", "BL", "AT"))
})

# -- bc_species_to_sp0 --------------------------------------------------------

test_that("bc_species_to_sp0 maps known codes to correct SP0", {
  # From the function's own @examples
  result <- bc_species_to_sp0(c("FDI", "FDC", "SW", "SE", "SS", "HW", "BL", "PLI"))
  expect_equal(result, c("F", "F", "S", "S", "S", "H", "B", "PL"))
})

test_that("bc_species_to_sp0 returns NA with warning for unknown code", {
  expect_warning(
    result <- bc_species_to_sp0("ZZZ"),
    regexp = "unrecognised"
  )
  expect_true(is.na(result))
})

test_that("bc_species_to_sp0 vectorises and preserves length", {
  codes  <- c("AT", "CW", "DR", "EP", "L", "MB", "PA", "PY", "Y")
  result <- bc_species_to_sp0(codes)
  expect_length(result, length(codes))
  expect_equal(result, c("AT", "C", "D", "E", "L", "MB", "PA", "PY", "Y"))
})

# -- bc_species_to_plant_code -------------------------------------------------

test_that("bc_species_to_osm maps known detailed codes to correct PlantCodes", {
  result <- bc_species_to_osm(c("FDI", "FDC", "HW", "BL", "PLI", "SW"))
  expect_equal(result, c("PSME", "PSME", "TSHE", "ABAM", "PICO", "PIGL"))
})

test_that("bc_species_to_osm returns NA with warning for unknown code", {
  expect_warning(
    result <- bc_species_to_osm("ZZZ"),
    regexp = "no OSM PlantCode"
  )
  expect_true(is.na(result))
})

# -- bc_species_to_plant_code -------------------------------------------------
