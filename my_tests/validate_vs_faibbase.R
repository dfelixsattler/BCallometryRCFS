# validate_vs_faibbase.R
#
# Numerical comparison of BCallometryR volume functions against the
# FAIBBase reference implementation.
#
# Requires FAIBBase source at C:/FAIBBase and its dependencies
# (data.table, fpCompare -- both already installed in C:/Rlibs).
#
# Run from RStudio after opening the BCallometryR project.

library(devtools)

# Load BCallometryR (re-loads any in-session changes)
load_all("C:/BCallometryR", quiet = TRUE)

# Load FAIBBase from source -- this brings in treeVolume() / treeProfile()
# Note: FAIBBase has many spatial dependencies; load_all skips non-essential ones.
load_all("C:/FAIBBase", quiet = TRUE)

# ---------------------------------------------------------------------------
# Test grid: varied species, BEC zones, DBH, height, with and without btop
# ---------------------------------------------------------------------------
trees <- expand.grid(
  bec_zone = c("CWH", "SBS", "IDF", "BWBS", "CDF", "ICH", "ESSF"),
  species  = c("H", "F", "S", "PL", "B", "C", "D"),
  dbh      = c(10, 20, 35, 55, 80),
  stringsAsFactors = FALSE
)

# Use a simple H-D relationship to assign realistic heights
# H = 1.3 + 3 * DBH^0.6, with +-15% random variation
set.seed(42)
trees$height <- round(1.3 + 3 * trees$dbh^0.6 *
                        runif(nrow(trees), 0.85, 1.15), 1)
trees$height <- pmax(trees$height, 5)    # floor at 5 m

# Add some broken-top trees (~15 % of the dataset)
trees$btop <- NA_real_
btop_idx <- sample(nrow(trees), round(nrow(trees) * 0.15))
trees$btop[btop_idx] <- round(trees$height[btop_idx] *
                                runif(length(btop_idx), 0.4, 0.85), 1)

cat("Test grid: ", nrow(trees), " trees\n", sep = "")
cat("Species   : ", paste(unique(trees$species), collapse = ", "), "\n", sep = "")
cat("BEC zones : ", paste(unique(trees$bec_zone), collapse = ", "), "\n", sep = "")
cat("DBH range : ", min(trees$dbh), "--", max(trees$dbh), " cm\n", sep = "")
cat("Height range: ", min(trees$height), "--", max(trees$height), " m\n\n", sep = "")

# ---------------------------------------------------------------------------
# Helper: run FAIBBase treeVolume() over the grid
# ---------------------------------------------------------------------------
run_faibbase <- function(volume_name, stump_h = 0.3, utop = 10) {
  # FAIBBase known bug: when no part of a tree is merchantable (all DIB < utop),
  # treeVolume("MER") returns numeric(0) and silently drops that tree from the
  # result vector.  Work-around: run single-tree calls for MER and substitute 0
  # for numeric(0) returns (0 is the correct merchantable volume).
  if (volume_name == "MER") {
    vapply(seq_len(nrow(trees)), function(i) {
      result <- treeVolume(
        taperEquationForm = "KBEC",
        FIZorBEC    = trees$bec_zone[i],
        species     = trees$species[i],
        DBH         = trees$dbh[i],
        height      = trees$height[i],
        BTOPHeight  = trees$btop[i],
        volumeName  = "MER",
        stumpHeight = stump_h,
        UTOPDIB     = utop
      )
      if (length(result) == 0L) 0 else result
    }, numeric(1))
  } else {
    treeVolume(
      taperEquationForm = "KBEC",
      FIZorBEC          = trees$bec_zone,
      species           = trees$species,
      DBH               = trees$dbh,
      height            = trees$height,
      BTOPHeight        = trees$btop,
      volumeName        = volume_name,
      stumpHeight       = stump_h,
      UTOPDIB           = utop
    )
  }
}

# ---------------------------------------------------------------------------
# Helper: run BCallometryR tree_volume() over the grid
# ---------------------------------------------------------------------------
run_bcallometry <- function(volume_type, stump_h = 0.3, utop = 10) {
  tree_volume(
    bec_zone     = trees$bec_zone,
    species      = trees$species,
    dbh          = trees$dbh,
    height       = trees$height,
    volume_type  = volume_type,
    stump_height = stump_h,
    utop_dib     = utop,
    btop_height  = trees$btop
  )
}

# ---------------------------------------------------------------------------
# Compare
# ---------------------------------------------------------------------------
compare <- function(label, ref, our) {
  stopifnot(length(ref) == length(our))
  diff    <- our - ref
  rel_err <- ifelse(abs(ref) > 1e-9, abs(diff) / ref * 100, NA_real_)
  cat(sprintf("%-10s  n=%d  max|diff|=%.2e m3  max rel err=%.4f %%  match(1e-4)=%s\n",
              label, sum(!is.na(diff)),
              max(abs(diff), na.rm = TRUE),
              max(rel_err, na.rm = TRUE),
              isTRUE(all.equal(ref, our, tolerance = 1e-4))))
}

cat("=== Comparison: BCallometryR vs FAIBBase (default stump=0.3, utop=10) ===\n\n")
cat("Note: FAIBBase bug -- treeVolume('MER') silently drops trees with no\n")
cat("      merchantable portion (all DIB < utop).  run_faibbase() patches\n")
cat("      this by running single-tree calls and substituting 0.\n\n")

wsv_ref <- run_faibbase("WSV");   wsv_bc  <- run_bcallometry("WSV")
mer_ref <- run_faibbase("MER");   mer_bc  <- run_bcallometry("MER")
stp_ref <- run_faibbase("STUMP"); stp_bc  <- run_bcallometry("STUMP")

compare("WSV",   wsv_ref, wsv_bc)
compare("MER",   mer_ref, mer_bc)
compare("STUMP", stp_ref, stp_bc)

# ---------------------------------------------------------------------------
# NMR sanity check: WSV ≈ MER + NMR + STUMP  (internal consistency)
# ---------------------------------------------------------------------------
nmr_bc  <- run_bcallometry("NMR")
residual <- wsv_bc - (mer_bc + nmr_bc + stp_bc)
cat(sprintf("\nNMR identity (WSV - MER - NMR - STUMP):  max|resid|=%.2e m3\n",
            max(abs(residual), na.rm = TRUE)))

# ---------------------------------------------------------------------------
# Custom merchantable limits
# ---------------------------------------------------------------------------
cat("\n=== Custom limits: stump=0.5 m, utop_dib=7 cm ===\n\n")
wsv_ref2 <- run_faibbase("WSV",   0.5, 7); wsv_bc2 <- run_bcallometry("WSV",   0.5, 7)
mer_ref2 <- run_faibbase("MER",   0.5, 7); mer_bc2 <- run_bcallometry("MER",   0.5, 7)
stp_ref2 <- run_faibbase("STUMP", 0.5, 7); stp_bc2 <- run_bcallometry("STUMP", 0.5, 7)
compare("WSV",   wsv_ref2, wsv_bc2)
compare("MER",   mer_ref2, mer_bc2)
compare("STUMP", stp_ref2, stp_bc2)

# ---------------------------------------------------------------------------
# Spot-check the largest differences (if any)
# ---------------------------------------------------------------------------
spot_check <- function(label, ref, our) {
  big <- which(abs(our - ref) > 1e-4)
  if (length(big) > 0) {
    cat(sprintf("\n--- Trees with %s |diff| > 1e-4 ---\n", label))
    print(cbind(trees[big, c("bec_zone","species","dbh","height","btop")],
                faibbase    = ref[big],
                bcallometry = our[big],
                diff        = our[big] - ref[big]))
  } else {
    cat(sprintf("\nAll %s differences < 1e-4 m3.\n", label))
  }
}
spot_check("WSV", wsv_ref, wsv_bc)
spot_check("MER", mer_ref, mer_bc)
