# =============================================================================
# BCallometryR — Tree Biomass Estimation with PSP Data
#
# This script demonstrates the biomass workflow using simulated PSP data:
#   1. Simulate a multi-plot, multi-species PSP dataset
#   2. DBH-only biomass (most common for large PSP archives)
#   3. DBH + height biomass (for trees with measured heights)
#   4. Component breakdown: wood, bark, branches, foliage
#   5. Compare Lambert2005 vs Ung2008 equations
#   6. Plot-level aggregation to Mg/ha
# =============================================================================

library(BCallometryR)

set.seed(2025)

# -----------------------------------------------------------------------------
# 1. Simulate PSP data
# -----------------------------------------------------------------------------
n_plots      <- 12
species_list <- c("lodgepole pine", "white spruce",
                  "trembling aspen", "western hemlock")

sp_params <- list(
  "lodgepole pine"  = list(a = 1.85, b = 0.130),
  "white spruce"    = list(a = 1.70, b = 0.115),
  "trembling aspen" = list(a = 1.60, b = 0.120),
  "western hemlock" = list(a = 1.45, b = 0.100)
)

plot_meta <- data.frame(
  PLOT_ID  = paste0("PSP-", formatC(1:n_plots, width = 3, flag = "0")),
  BEC_ZONE = sample(c("SBS", "ICH", "CWH", "IDF"), n_plots, replace = TRUE),
  a_shift  = rnorm(n_plots, mean = 0, sd = 0.20)
)

sim_trees <- do.call(rbind, lapply(seq_len(n_plots), function(i) {
  n_trees <- sample(25:40, 1)
  sp      <- sample(species_list, n_trees, replace = TRUE,
                    prob = c(0.40, 0.30, 0.20, 0.10))
  dbh     <- round(runif(n_trees, min = 5, max = 60), 1)

  height <- mapply(function(d, s) {
    p  <- sp_params[[s]]
    ht <- ht_naslund(d, a = p$a + plot_meta$a_shift[i], b = p$b)
    round(pmax(ht + rnorm(1, 0, sd = 1.0), 1.5), 1)
  }, dbh, sp)

  unmeasured          <- sample(n_trees, size = round(n_trees * 0.70))
  height_obs          <- height
  height_obs[unmeasured] <- NA

  data.frame(
    PLOT_ID   = plot_meta$PLOT_ID[i],
    BEC_ZONE  = plot_meta$BEC_ZONE[i],
    TREE_ID   = paste0(plot_meta$PLOT_ID[i], "-",
                       formatC(seq_len(n_trees), width = 2, flag = "0")),
    SPECIES   = sp,
    DBH       = dbh,
    HEIGHT    = height_obs
  )
}))

cat("Simulated PSP dataset\n")
cat("  Plots    :", n_plots, "\n")
cat("  Trees    :", nrow(sim_trees), "\n")
cat("  Measured :", sum(!is.na(sim_trees$HEIGHT)),
    sprintf("(%.0f%%)\n", 100 * mean(!is.na(sim_trees$HEIGHT))))

# -----------------------------------------------------------------------------
# 2. DBH-only biomass
# -----------------------------------------------------------------------------
# Most common scenario for historic PSP archives where heights were not
# systematically collected.

sim_trees$BIOMASS_DBH <- biomass_tree(
  species = sim_trees$SPECIES,
  dbh     = sim_trees$DBH
)

cat("\n--- DBH-only biomass summary by species ---\n")
for (sp in species_list) {
  sub <- sim_trees[sim_trees$SPECIES == sp, ]
  cat(sprintf("  %-20s  n=%3d  mean DBH=%5.1f cm  mean biomass=%6.1f kg\n",
              sp, nrow(sub), mean(sub$DBH),
              mean(sub$BIOMASS_DBH, na.rm = TRUE)))
}

# -----------------------------------------------------------------------------
# 3. DBH + height biomass (for trees with measured heights)
# -----------------------------------------------------------------------------
has_ht <- !is.na(sim_trees$HEIGHT)

bio_ht <- rep(NA_real_, nrow(sim_trees))
bio_ht[has_ht] <- biomass_tree(
  species = sim_trees$SPECIES[has_ht],
  dbh     = sim_trees$DBH[has_ht],
  height  = sim_trees$HEIGHT[has_ht]
)
sim_trees$BIOMASS_HT <- bio_ht

both      <- sim_trees[has_ht, ]
diff_pct  <- 100 * (both$BIOMASS_HT - both$BIOMASS_DBH) / both$BIOMASS_DBH

cat(sprintf("\n--- DBH-only vs DBH+height (%d trees with both) ---\n", nrow(both)))
cat(sprintf("  Mean difference: %.1f%%  (SD: %.1f%%)\n",
            mean(diff_pct, na.rm = TRUE),
            sd(diff_pct,   na.rm = TRUE)))

# -----------------------------------------------------------------------------
# 4. Component breakdown (wood, bark, branches, foliage)
# -----------------------------------------------------------------------------
cat("\n--- Component biomass: first 5 trees ---\n")
comp <- biomass_components(
  species = sim_trees$SPECIES[1:5],
  dbh     = sim_trees$DBH[1:5]
)
print(cbind(
  species = sim_trees$SPECIES[1:5],
  dbh     = sim_trees$DBH[1:5],
  round(comp, 2)
), row.names = FALSE)

# Wood fraction
comp_all   <- biomass_components(species = sim_trees$SPECIES,
                                 dbh     = sim_trees$DBH)
wood_frac  <- comp_all$wood / comp_all$total
cat(sprintf("\nMean wood fraction across all trees: %.1f%%\n",
            100 * mean(wood_frac, na.rm = TRUE)))

# -----------------------------------------------------------------------------
# 5. Lambert 2005 vs Ung 2008
# -----------------------------------------------------------------------------
cat("\n--- Lambert2005 vs Ung2008 at DBH = 20, 40, 60 cm ---\n")
dbh_check <- c(20, 40, 60)
for (sp in c("lodgepole pine", "white spruce", "trembling aspen")) {
  bl <- biomass_tree(sp, dbh_check, paper_source = "Lambert2005")
  bu <- biomass_tree(sp, dbh_check, paper_source = "Ung2008")
  cat(sprintf("  %-20s  L2005: %s  U2008: %s  diff%%: %s\n",
              sp,
              paste(round(bl, 1), collapse = "/"),
              paste(round(bu, 1), collapse = "/"),
              paste(round(100 * (bu - bl) / bl, 1), collapse = "/")))
}

# -----------------------------------------------------------------------------
# 6. Plot-level aggregation (Mg/ha)
# -----------------------------------------------------------------------------
PLOT_AREA_HA <- 0.04           # 20 m × 20 m subplot
EXP_FACTOR   <- 1 / PLOT_AREA_HA

cat("\n--- Plot-level biomass (Mg/ha) ---\n")
plot_bio <- do.call(rbind, lapply(unique(sim_trees$PLOT_ID), function(pid) {
  sub   <- sim_trees[sim_trees$PLOT_ID == pid, ]
  total <- sum(sub$BIOMASS_DBH, na.rm = TRUE) * EXP_FACTOR / 1000  # kg → Mg
  data.frame(
    PLOT_ID       = pid,
    BEC_ZONE      = sub$BEC_ZONE[1],
    N_TREES       = nrow(sub),
    BIOMASS_MG_HA = round(total, 1)
  )
}))
print(plot_bio[order(plot_bio$BIOMASS_MG_HA, decreasing = TRUE), ],
      row.names = FALSE)
cat(sprintf("\nMean: %.1f Mg/ha  |  SD: %.1f  |  Range: %.1f–%.1f\n",
            mean(plot_bio$BIOMASS_MG_HA),
            sd(plot_bio$BIOMASS_MG_HA),
            min(plot_bio$BIOMASS_MG_HA),
            max(plot_bio$BIOMASS_MG_HA)))

# -----------------------------------------------------------------------------
# References
# -----------------------------------------------------------------------------
cat("\n--- Equation sources ---\n")
print(biomass_citations())
