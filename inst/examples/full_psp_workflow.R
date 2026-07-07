# =============================================================================
# BCallometryR — Complete PSP Workflow
#
# This script demonstrates the end-to-end workflow for a BC PSP dataset that
# has raw inventory species codes, DBH for all trees, and heights measured on
# only a subset of trees.
#
# Steps covered:
#   1. Load / simulate cleaned PSP data (your starting point in practice)
#   2. Species code crosswalk
#   3. Fit H-D models per species (mixed-effects with random plot effect;
#      automatic fallback to fixed-effects if nlme is unavailable or fails)
#   4. Predict missing heights
#   5. Biomass estimation (total and component)
#   6. Volume estimation (WSV and gross merchantable)
#   7. Plot-level summaries
# =============================================================================

library(BCallometryR)

# =============================================================================
# 1. Cleaned PSP data
#    In practice: replace this block with  trees <- read.csv("your_psp.csv")
#    or similar.  The data frame must contain at minimum:
#      PLOT_ID   — plot identifier
#      BEC_ZONE  — BEC zone code (used for species correction and taper eq.)
#      SPECIES   — raw BC inventory code, e.g. "PLI", "SW", "FDI", "HW", "AT"
#      DBH       — diameter at breast height (cm); required for all trees
#      HEIGHT    — total height (m); NA for unmeasured trees (~30-40% typical)
# =============================================================================

set.seed(42)

n_plots  <- 15
bec_pool <- c("SBS", "SBS", "SBS", "SBS", "SBS",   # 5 SBS
              "IDF", "IDF", "IDF", "IDF", "IDF",   # 5 IDF
              "CWH", "CWH", "CWH", "CWH", "CWH")   # 5 CWH

sp_by_bec <- list(
  SBS = c("PLI", "SW",  "AT"),   # lodgepole pine / white spruce / trembling aspen
  IDF = c("FDI", "PLI"),         # interior Douglas-fir / lodgepole pine
  CWH = c("HW",  "SS")           # western hemlock / Sitka spruce
)

# Naslund H-D parameters keyed by SP0 code (for height simulation only)
sp0_hd <- list(
  PL = list(a = 1.85, b = 0.130),
  S  = list(a = 1.70, b = 0.115),   # covers both SW (SBS) and SS (CWH)
  AT = list(a = 1.60, b = 0.120),
  F  = list(a = 1.55, b = 0.110),
  H  = list(a = 1.45, b = 0.100)
)

plot_meta <- data.frame(
  PLOT_ID  = paste0("PSP-", formatC(seq_len(n_plots), width = 3, flag = "0")),
  BEC_ZONE = bec_pool,
  a_shift  = rnorm(n_plots, mean = 0, sd = 0.25)   # site-productivity variation
)

trees <- do.call(rbind, lapply(seq_len(n_plots), function(i) {
  bec     <- plot_meta$BEC_ZONE[i]
  sps     <- sp_by_bec[[bec]]
  n_trees <- sample(25:40, 1)
  sp      <- sample(sps, n_trees, replace = TRUE)
  dbh     <- round(runif(n_trees, min = 5, max = 65), 1)

  # Build SP0 for height simulation lookup only
  sp_corr_i <- species_correction(sp, bec_zone = bec)
  sp0_i     <- bc_species_to_sp0(sp_corr_i)

  ht_true <- mapply(function(d, s0) {
    p  <- sp0_hd[[s0]]
    ht <- ht_naslund(d, a = p$a + plot_meta$a_shift[i], b = p$b)
    round(pmax(ht + rnorm(1, 0, 1.2), 1.5), 1)
  }, dbh, sp0_i)

  # ~35% of trees have no measured height (typical PSP scenario)
  unm          <- sample(n_trees, size = round(n_trees * 0.35))
  ht_obs       <- ht_true
  ht_obs[unm]  <- NA

  data.frame(
    PLOT_ID  = plot_meta$PLOT_ID[i],
    BEC_ZONE = bec,
    TREE_ID  = paste0(plot_meta$PLOT_ID[i], "-",
                      formatC(seq_len(n_trees), width = 2, flag = "0")),
    SPECIES  = sp,        # raw PSP codes — this is all you'd have in practice
    DBH      = dbh,
    HEIGHT   = ht_obs     # NA = not measured
  )
}))

cat(sprintf(
  "PSP dataset: %d plots | %d trees | %d heights measured (%.0f%%) | %d missing\n\n",
  n_plots, nrow(trees),
  sum(!is.na(trees$HEIGHT)), 100 * mean(!is.na(trees$HEIGHT)),
  sum(is.na(trees$HEIGHT))
))


# =============================================================================
# 2. Species code crosswalk
#    species_correction()    — BEC-zone-aware standardisation
#                              "PLI" → "PL", "FDI" → "FD"
#                              generic "S" → "SW"/"SE"/"SS" by zone
#                              generic "B" → "BL"/"BA" by zone
#    bc_species_to_sp0()     — corrected code → SP0 (for volume taper eq.)
#    bc_psp_to_biomass_name()— PSP code + BEC → common name (for biomass eq.)
# =============================================================================

trees$SPECIES_CORR <- species_correction(trees$SPECIES, trees$BEC_ZONE)
trees$SPECIES_SP0  <- bc_species_to_sp0(trees$SPECIES_CORR)
trees$SPECIES_NAME <- bc_species_to_biomass_name(trees$SPECIES, trees$BEC_ZONE)

cat("--- Species crosswalk summary ---\n")
xwalk <- unique(trees[, c("SPECIES", "BEC_ZONE", "SPECIES_CORR",
                           "SPECIES_SP0", "SPECIES_NAME")])
xwalk <- xwalk[order(xwalk$BEC_ZONE, xwalk$SPECIES), ]
print(xwalk, row.names = FALSE)
cat("\n")


# =============================================================================
# 3. Fit H-D models
#    One model per SP0 species group, fitted on measured trees only.
#    Strategy:
#      (a) Try nlme mixed-effects with a random plot effect on parameter 'a'.
#          This accounts for plot-level site productivity variation and gives
#          more reliable fixed-effect estimates when data span many plots.
#      (b) If nlme is unavailable or convergence fails, fall back to pooled
#          fixed-effects nls across all plots.
# =============================================================================

measured <- trees[!is.na(trees$HEIGHT) & !is.na(trees$SPECIES_SP0), ]
sp0_list <- sort(unique(trees$SPECIES_SP0[!is.na(trees$SPECIES_SP0)]))

cat("--- H-D model fitting ---\n")
cat(sprintf("  %-4s  %-15s  %5s  %6s\n", "SP0", "Method", "n", "Marg.R2"))
cat(sprintf("  %-4s  %-15s  %5s  %6s\n", "----", "---------------", "-----", "------"))

hd_fits <- lapply(sp0_list, function(sp0) {
  sub <- measured[measured$SPECIES_SP0 == sp0, ]
  if (nrow(sub) < 10) {
    cat(sprintf("  %-4s  %-15s  %5d  %6s\n", sp0, "too few trees", nrow(sub), "—"))
    return(list(sp0 = sp0, fit = NULL, method = "none"))
  }

  # (a) Mixed-effects (requires nlme; random effect on asymptote 'a' by plot)
  n_plots_sp <- length(unique(sub$PLOT_ID))
  fit <- NULL
  method <- "none"

  if (n_plots_sp >= 5 && requireNamespace("nlme", quietly = TRUE)) {
    fit <- tryCatch(
      fit_hd_model(sub, model = "naslund", group_col = "PLOT_ID"),
      error = function(e) NULL
    )
    if (!is.null(fit)) method <- "mixed-effects"
  }

  # (b) Fallback: fixed effects
  if (is.null(fit)) {
    fit <- tryCatch(
      fit_hd_model(sub, model = "naslund"),
      error = function(e) NULL
    )
    if (!is.null(fit)) method <- "fixed-effects"
  }

  if (!is.null(fit)) {
    cat(sprintf("  %-4s  %-15s  %5d  %6.4f\n",
                sp0, method, nrow(sub), fit$r2_marginal))
  } else {
    cat(sprintf("  %-4s  %-15s  %5d  %6s\n", sp0, "FAILED", nrow(sub), "—"))
  }

  list(sp0 = sp0, fit = fit, method = method)
})
names(hd_fits) <- sp0_list
cat("\n")


# =============================================================================
# 4. Predict missing heights
#    Fixed-effect coefficients are used for population-average prediction.
#    (For remeasurement plots that were in the training data, plot-level BLUPs
#    from the nlme object would reduce bias: predict(fit$fit, level = 1).)
# =============================================================================

trees$HEIGHT_FILLED <- trees$HEIGHT

for (m in hd_fits) {
  if (is.null(m$fit)) next
  cf      <- m$fit$coefficients
  missing <- is.na(trees$HEIGHT) & !is.na(trees$SPECIES_SP0) &
             trees$SPECIES_SP0 == m$sp0
  if (!any(missing)) next
  trees$HEIGHT_FILLED[missing] <- round(
    ht_from_dbh(trees$DBH[missing], model = "naslund",
                a = cf["a"], b = cf["b"]),
    1
  )
}

n_filled  <- sum(!is.na(trees$HEIGHT_FILLED) & is.na(trees$HEIGHT))
n_still_na <- sum(is.na(trees$HEIGHT_FILLED))
cat(sprintf("Heights filled: %d  |  Still missing (no model): %d\n\n",
            n_filled, n_still_na))


# =============================================================================
# 5. Biomass estimation
#    biomass_tree()       — total aboveground biomass (kg); uses height when
#                           available, falls back to DBH-only if still NA
#    biomass_components() — wood / bark / branches / foliage breakdown
# =============================================================================

# Total aboveground biomass — use filled height where available
trees$BIOMASS_KG <- biomass_tree(
  species = trees$SPECIES_NAME,
  dbh     = trees$DBH,
  height  = trees$HEIGHT_FILLED   # NA triggers DBH-only equation automatically
)

# Component breakdown (also uses filled height)
comp <- biomass_components(
  species = trees$SPECIES_NAME,
  dbh     = trees$DBH,
  height  = trees$HEIGHT_FILLED
)
trees$BIO_WOOD_KG   <- comp$wood
trees$BIO_BARK_KG   <- comp$bark
trees$BIO_BRANCH_KG <- comp$branches
trees$BIO_FOLIAGE_KG<- comp$foliage

cat("--- Biomass summary by SP0 ---\n")
cat(sprintf("  %-4s  %5s  %8s  %8s\n", "SP0", "n", "mean(kg)", "SD(kg)"))
cat(sprintf("  %-4s  %5s  %8s  %8s\n", "----", "-----", "--------", "--------"))
for (sp0 in sp0_list) {
  sub <- trees[!is.na(trees$SPECIES_SP0) & trees$SPECIES_SP0 == sp0 &
               !is.na(trees$BIOMASS_KG), ]
  cat(sprintf("  %-4s  %5d  %8.1f  %8.1f\n",
              sp0, nrow(sub),
              mean(sub$BIOMASS_KG), sd(sub$BIOMASS_KG)))
}
cat("\n")


# =============================================================================
# 6. Volume estimation
#    tree_volume() implements the Kozak (2002) KBEC variable-exponent taper
#    equation.  Height is required — trees without a filled height are skipped.
# =============================================================================

has_ht <- !is.na(trees$HEIGHT_FILLED) & !is.na(trees$SPECIES_SP0)

trees$WSV_M3 <- NA_real_
trees$MER_M3 <- NA_real_

trees$WSV_M3[has_ht] <- tree_volume(
  bec_zone = trees$BEC_ZONE[has_ht],
  species  = trees$SPECIES_SP0[has_ht],
  dbh      = trees$DBH[has_ht],
  height   = trees$HEIGHT_FILLED[has_ht]
)

trees$MER_M3[has_ht] <- tree_volume(
  bec_zone    = trees$BEC_ZONE[has_ht],
  species     = trees$SPECIES_SP0[has_ht],
  dbh         = trees$DBH[has_ht],
  height      = trees$HEIGHT_FILLED[has_ht],
  volume_type = "MER"
)

cat("--- Volume summary by SP0 ---\n")
cat(sprintf("  %-4s  %5s  %9s  %9s\n", "SP0", "n", "WSV(m3)", "MER(m3)"))
cat(sprintf("  %-4s  %5s  %9s  %9s\n", "----", "-----", "---------", "---------"))
for (sp0 in sp0_list) {
  sub <- trees[!is.na(trees$SPECIES_SP0) & trees$SPECIES_SP0 == sp0 &
               !is.na(trees$WSV_M3), ]
  cat(sprintf("  %-4s  %5d  %9.4f  %9.4f\n",
              sp0, nrow(sub),
              mean(sub$WSV_M3), mean(sub$MER_M3)))
}
cat("\n")


# =============================================================================
# 7. Plot-level summaries
#    Expansion factor converts per-tree values to per-hectare.
#    Adjust PLOT_AREA_HA to match your plot design.
# =============================================================================

PLOT_AREA_HA <- 0.04        # 400 m² (20 m × 20 m) fixed-area subplot
EXP          <- 1 / PLOT_AREA_HA

plot_summary <- do.call(rbind, lapply(unique(trees$PLOT_ID), function(pid) {
  sub <- trees[trees$PLOT_ID == pid, ]
  data.frame(
    PLOT_ID       = pid,
    BEC_ZONE      = sub$BEC_ZONE[1],
    N_TREES       = nrow(sub),
    BIOMASS_MG_HA = round(sum(sub$BIOMASS_KG,  na.rm = TRUE) * EXP / 1000, 2),
    WSV_M3_HA     = round(sum(sub$WSV_M3,       na.rm = TRUE) * EXP,       1),
    MER_M3_HA     = round(sum(sub$MER_M3,       na.rm = TRUE) * EXP,       1)
  )
}))

plot_summary <- plot_summary[order(plot_summary$BEC_ZONE, plot_summary$PLOT_ID), ]

cat("--- Plot-level summary ---\n")
print(plot_summary, row.names = FALSE)

cat(sprintf(
  "\nMean by BEC zone (Mg/ha biomass | m³/ha WSV | m³/ha MER):\n"
))
for (bz in sort(unique(plot_summary$BEC_ZONE))) {
  ps <- plot_summary[plot_summary$BEC_ZONE == bz, ]
  cat(sprintf("  %-5s  %6.2f Mg/ha  |  %6.1f m³/ha WSV  |  %6.1f m³/ha MER\n",
              bz,
              mean(ps$BIOMASS_MG_HA),
              mean(ps$WSV_M3_HA),
              mean(ps$MER_M3_HA)))
}

cat("\nDone.\n")
