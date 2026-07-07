# =============================================================================
# BCallometryR — Tree Volume Estimation with PSP Data
#
# This script demonstrates the Kozak (2002) KBEC volume workflow:
#   1. Simulate a multi-plot PSP dataset with BEC zone and species codes
#   2. Compute whole-stem volume (WSV)
#   3. Compare WSV, merchantable (MER), and stump volume
#   4. Handle broken-top trees
#   5. Plot-level growing stock (m³/ha)
#   6. Bridging workflow: predict missing heights, then compute volume
# =============================================================================

library(BCallometryR)

set.seed(2025)

# -----------------------------------------------------------------------------
# 1. Simulate PSP data
# -----------------------------------------------------------------------------
n_plots   <- 10
bec_zones <- c("SBS", "CWH", "IDF")
sp_by_bec <- list(SBS = c("PL","S"), CWH = c("H","S"), IDF = c("D","PL"))

sp_hd <- list(
  PL = list(a = 1.85, b = 0.130),
  S  = list(a = 1.65, b = 0.105),
  H  = list(a = 1.45, b = 0.100),
  D  = list(a = 1.55, b = 0.110)
)

plot_meta <- data.frame(
  PLOT_ID  = paste0("PSP-", formatC(1:n_plots, width = 3, flag = "0")),
  BEC_ZONE = rep(bec_zones, length.out = n_plots),
  a_shift  = rnorm(n_plots, mean = 0, sd = 0.20)
)

sim_trees <- do.call(rbind, lapply(seq_len(n_plots), function(i) {
  bec     <- plot_meta$BEC_ZONE[i]
  sps     <- sp_by_bec[[bec]]
  n_trees <- sample(20:35, 1)
  sp      <- sample(sps, n_trees, replace = TRUE)
  dbh     <- round(runif(n_trees, min = 7, max = 65), 1)

  height <- mapply(function(d, s) {
    p  <- sp_hd[[s]]
    ht <- ht_naslund(d, a = p$a + plot_meta$a_shift[i], b = p$b)
    round(pmax(ht + rnorm(1, 0, sd = 1.0), 2.0), 1)
  }, dbh, sp)

  data.frame(
    PLOT_ID  = plot_meta$PLOT_ID[i],
    BEC_ZONE = bec,
    TREE_ID  = paste0(plot_meta$PLOT_ID[i], "-",
                      formatC(seq_len(n_trees), width = 2, flag = "0")),
    SPECIES  = sp,
    DBH      = dbh,
    HEIGHT   = height
  )
}))

cat("Simulated dataset:", nrow(sim_trees), "trees in",
    n_plots, "plots —", paste(bec_zones, collapse = ", "), "\n\n")

# -----------------------------------------------------------------------------
# 2. Whole-stem volume (WSV)
# -----------------------------------------------------------------------------
sim_trees$WSV <- tree_volume(
  bec_zone = sim_trees$BEC_ZONE,
  species  = sim_trees$SPECIES,
  dbh      = sim_trees$DBH,
  height   = sim_trees$HEIGHT
)

cat("--- WSV summary by species ---\n")
for (sp in c("D", "H", "PL", "S")) {
  sub <- sim_trees[sim_trees$SPECIES == sp & !is.na(sim_trees$WSV), ]
  if (nrow(sub) == 0) next
  cat(sprintf("  %-4s  n=%3d  DBH: %4.1f–%4.1f cm  WSV: %.3f–%.3f m³\n",
              sp, nrow(sub),
              min(sub$DBH), max(sub$DBH),
              min(sub$WSV), max(sub$WSV)))
}

# -----------------------------------------------------------------------------
# 3. All three volume types
# -----------------------------------------------------------------------------
cat("\n--- WSV vs MER vs STUMP for a sample of trees ---\n")

sim_trees$MER   <- tree_volume(sim_trees$BEC_ZONE, sim_trees$SPECIES,
                               sim_trees$DBH, sim_trees$HEIGHT,
                               volume_type = "MER")
sim_trees$STUMP <- tree_volume(sim_trees$BEC_ZONE, sim_trees$SPECIES,
                               sim_trees$DBH, sim_trees$HEIGHT,
                               volume_type = "STUMP")

sample_idx <- order(sim_trees$DBH)[round(seq(1, nrow(sim_trees), length.out = 8))]
print(round(sim_trees[sample_idx, c("SPECIES","DBH","HEIGHT","WSV","MER","STUMP")], 4),
      row.names = FALSE)

# MER as a fraction of WSV
valid     <- !is.na(sim_trees$WSV) & sim_trees$WSV > 0
mer_ratio <- sim_trees$MER[valid] / sim_trees$WSV[valid]
cat(sprintf("\nMER/WSV ratio  mean: %.3f  range: %.3f–%.3f\n",
            mean(mer_ratio, na.rm = TRUE),
            min(mer_ratio, na.rm = TRUE),
            max(mer_ratio, na.rm = TRUE)))

# -----------------------------------------------------------------------------
# 4. Broken-top trees
# -----------------------------------------------------------------------------
cat("\n--- Broken-top example (H, 35 cm DBH, 30 m height, CWH) ---\n")
v_full <- tree_volume("CWH", "H", 35, 30)
btops  <- seq(10, 28, by = 2)
# tree_volume iterates over trees; use sapply to vary btop_height for one tree
v_btop <- sapply(btops, function(bt)
  tree_volume("CWH", "H", 35, 30, btop_height = bt))

cat(sprintf("  Full WSV: %.4f m³\n", v_full))
print(data.frame(btop_height = btops,
                 WSV         = round(v_btop, 4),
                 pct_loss    = round(100 * (v_full - v_btop) / v_full, 1)),
      row.names = FALSE)

# -----------------------------------------------------------------------------
# 5. Plot-level growing stock (m³/ha)
# -----------------------------------------------------------------------------
PLOT_AREA_HA <- 0.04
EXP_FACTOR   <- 1 / PLOT_AREA_HA

cat("\n--- Plot-level growing stock ---\n")
plot_vol <- do.call(rbind, lapply(unique(sim_trees$PLOT_ID), function(pid) {
  sub <- sim_trees[sim_trees$PLOT_ID == pid & !is.na(sim_trees$WSV), ]
  data.frame(
    PLOT_ID   = pid,
    BEC_ZONE  = sub$BEC_ZONE[1],
    N_TREES   = nrow(sub),
    WSV_M3_HA = round(sum(sub$WSV,   na.rm = TRUE) * EXP_FACTOR, 1),
    MER_M3_HA = round(sum(sub$MER,   na.rm = TRUE) * EXP_FACTOR, 1)
  )
}))
print(plot_vol[order(plot_vol$WSV_M3_HA, decreasing = TRUE), ],
      row.names = FALSE)
cat(sprintf("\nMean WSV: %.1f  |  Mean MER: %.1f  m³/ha\n",
            mean(plot_vol$WSV_M3_HA), mean(plot_vol$MER_M3_HA)))

# -----------------------------------------------------------------------------
# 6. Bridging workflow: predict missing heights → volume
# -----------------------------------------------------------------------------
cat("\n--- Bridging workflow: H-D prediction → volume ---\n")

# Introduce 60% missing heights
set.seed(99)
sim_trees$HEIGHT_OBS <- sim_trees$HEIGHT
sim_trees$HEIGHT_OBS[sample(nrow(sim_trees),
                            round(nrow(sim_trees) * 0.60))] <- NA
cat(sprintf("Heights available: %d / %d\n",
            sum(!is.na(sim_trees$HEIGHT_OBS)), nrow(sim_trees)))

# Fit per-species H-D models
fits <- lapply(c("D", "H", "PL", "S"), function(sp) {
  sub <- sim_trees[!is.na(sim_trees$HEIGHT_OBS) & sim_trees$SPECIES == sp, ]
  if (nrow(sub) < 5) return(NULL)
  fit <- tryCatch(fit_hd_model(sub, model = "naslund"), error = function(e) NULL)
  list(sp = sp, fit = fit)
})

# Predict missing heights
sim_trees$HEIGHT_FILLED <- sim_trees$HEIGHT_OBS
for (f in fits) {
  if (is.null(f$fit)) next
  cf      <- f$fit$coefficients
  missing <- is.na(sim_trees$HEIGHT_OBS) & sim_trees$SPECIES == f$sp
  sim_trees$HEIGHT_FILLED[missing] <- round(
    ht_from_dbh(sim_trees$DBH[missing],
                model = "naslund", a = cf["a"], b = cf["b"]),
    1
  )
}
cat(sprintf("After H-D fill: %d / %d heights available\n",
            sum(!is.na(sim_trees$HEIGHT_FILLED)), nrow(sim_trees)))

# Compute volume on the filled heights
complete <- !is.na(sim_trees$HEIGHT_FILLED)
sim_trees$WSV_FILLED <- NA_real_
sim_trees$WSV_FILLED[complete] <- tree_volume(
  bec_zone = sim_trees$BEC_ZONE[complete],
  species  = sim_trees$SPECIES[complete],
  dbh      = sim_trees$DBH[complete],
  height   = sim_trees$HEIGHT_FILLED[complete]
)

cat(sprintf("Volume computed for %d trees  (%.0f%% coverage)\n",
            sum(!is.na(sim_trees$WSV_FILLED)),
            100 * mean(!is.na(sim_trees$WSV_FILLED))))

# -----------------------------------------------------------------------------
# References
# -----------------------------------------------------------------------------
cat("\n--- Equation sources ---\n")
print(volume_citations())
