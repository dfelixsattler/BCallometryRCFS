# =============================================================================
# BCallometryR — Height-Diameter modelling with simulated PSP data
#
# This script simulates a small Permanent Sample Plot (PSP) dataset with
# multiple plots, species, and measurement cycles, then demonstrates:
#   1. Exploratory data summary
#   2. Fitting a fixed-effects H-D model (single species)
#   3. Comparing model forms for best fit
#   4. Fitting a mixed-effects model (random plot effect)
#   5. Predicting missing heights — the most common operational use
# =============================================================================

library(BCallometryR)

set.seed(2024)

# -----------------------------------------------------------------------------
# 1. Simulate PSP data
# -----------------------------------------------------------------------------
# Design: 15 plots, 20-35 trees per plot, 3 species (Fd, Pl, Hw),
#         "true" H-D relationship varies by species and has a plot-level
#         random effect on the asymptote (as one would expect in real PSPs).

n_plots   <- 15
species_params <- list(
  Fd = list(a = 1.60, b = 0.110),   # Naslund coefficients by species
  Pl = list(a = 1.80, b = 0.130),
  Hw = list(a = 1.40, b = 0.095)
)

# Random plot-level shifts on parameter 'a' (site productivity variation)
plot_meta <- data.frame(
  PLOT_ID  = paste0("PSP-", formatC(1:n_plots, width = 3, flag = "0")),
  BEC_ZONE = sample(c("ICH", "SBS", "ESSF", "CWH"), n_plots, replace = TRUE),
  a_shift  = rnorm(n_plots, mean = 0, sd = 0.25)   # site-level random effect
)

sim_trees <- do.call(rbind, lapply(seq_len(n_plots), function(i) {
  n_trees <- sample(20:35, 1)
  sp      <- sample(names(species_params), n_trees, replace = TRUE,
                    prob = c(0.5, 0.3, 0.2))
  dbh     <- round(runif(n_trees, min = 5, max = 65), 1)

  height <- mapply(function(d, s) {
    p     <- species_params[[s]]
    a_eff <- p$a + plot_meta$a_shift[i]          # plot-level variation
    ht    <- ht_naslund(d, a = a_eff, b = p$b)
    ht    + rnorm(1, 0, sd = 1.2)                # measurement noise
  }, dbh, sp)

  # Introduce ~20% missing heights (unmeasured trees — very common in PSPs)
  missing_idx          <- sample(n_trees, size = round(n_trees * 0.20))
  height[missing_idx]  <- NA

  data.frame(
    PLOT_ID  = plot_meta$PLOT_ID[i],
    BEC_ZONE = plot_meta$BEC_ZONE[i],
    TREE_ID  = paste0(plot_meta$PLOT_ID[i], "-", formatC(seq_len(n_trees), width = 2, flag = "0")),
    SPECIES  = sp,
    DBH      = dbh,
    HEIGHT   = round(pmax(height, 1.5), 1)       # floor at 1.5 m
  )
}))

cat("Simulated PSP dataset\n")
cat("  Plots    :", n_plots, "\n")
cat("  Trees    :", nrow(sim_trees), "\n")
cat("  Measured :", sum(!is.na(sim_trees$HEIGHT)), "\n")
cat("  Missing  :", sum( is.na(sim_trees$HEIGHT)), "\n\n")
print(table(sim_trees$SPECIES, sim_trees$BEC_ZONE))

# -----------------------------------------------------------------------------
# 2. Exploratory summary — observed H-D by species
# -----------------------------------------------------------------------------
measured <- sim_trees[!is.na(sim_trees$HEIGHT), ]

cat("\n--- Observed H-D summary (measured trees) ---\n")
for (sp in names(species_params)) {
  sub <- measured[measured$SPECIES == sp, ]
  cat(sprintf("  %-4s  n=%3d  DBH: %4.1f–%4.1f cm   HT: %4.1f–%4.1f m\n",
              sp, nrow(sub),
              min(sub$DBH), max(sub$DBH),
              min(sub$HEIGHT), max(sub$HEIGHT)))
}

# -----------------------------------------------------------------------------
# 3. Fixed-effects model — single species (Fd)
# -----------------------------------------------------------------------------
fd_data <- measured[measured$SPECIES == "Fd", ]

cat("\n--- Fixed-effects model: Fd, Naslund ---\n")
fd_naslund <- fit_hd_model(fd_data, model = "naslund")
cat("  Coefficients :", round(fd_naslund$coefficients, 4), "\n")
cat("  Marginal R²  :", round(fd_naslund$r2_marginal, 4), "\n")

# -----------------------------------------------------------------------------
# 4. Compare model forms for Fd
# -----------------------------------------------------------------------------
cat("\n--- Model comparison: Fd (fixed effects) ---\n")
cat(sprintf("  %-10s  %s\n", "Model", "Marginal R²"))
cat(sprintf("  %-10s  %s\n", "-----", "-----------"))

model_forms <- c("naslund", "curtis", "korf", "weibull", "richards")
comparison  <- lapply(model_forms, function(m) {
  fit <- tryCatch(
    fit_hd_model(fd_data, model = m),
    error = function(e) NULL
  )
  if (is.null(fit)) return(data.frame(model = m, r2 = NA))
  data.frame(model = m, r2 = fit$r2_marginal)
})
comparison <- do.call(rbind, comparison)
comparison <- comparison[order(-comparison$r2, na.last = TRUE), ]

for (i in seq_len(nrow(comparison))) {
  cat(sprintf("  %-10s  %.4f\n", comparison$model[i],
              ifelse(is.na(comparison$r2[i]), 0, comparison$r2[i])))
}

best_model <- comparison$model[1]
cat(sprintf("\n  Best form: %s\n", best_model))

# -----------------------------------------------------------------------------
# 5. Mixed-effects model — all species separately, random plot effect
# -----------------------------------------------------------------------------
if (requireNamespace("nlme", quietly = TRUE)) {
  cat("\n--- Mixed-effects model: Fd, Naslund + random plot effect ---\n")
  fd_me <- tryCatch(
    fit_hd_model(fd_data, model = "naslund", group_col = "PLOT_ID"),
    error = function(e) { cat("  nlme failed:", conditionMessage(e), "\n"); NULL }
  )
  if (!is.null(fd_me)) {
    cat("  Fixed coefficients:", round(fd_me$coefficients, 4), "\n")
    cat("  Marginal R²       :", round(fd_me$r2_marginal, 4), "\n")
  }
} else {
  cat("\n  (Install 'nlme' to demonstrate mixed-effects fitting)\n")
}

# -----------------------------------------------------------------------------
# 6. Predict missing heights — the key operational use case
# -----------------------------------------------------------------------------
cat("\n--- Predicting missing heights ---\n")

# Fit one model per species pooling all measured trees.
# Pooling across measurement cycles is correct — consistent with FAIBCompiler,
# which develops regional H-D coefficient tables from the full archive rather
# than refitting per cruise cycle.  The random plot effect (if used) accounts
# for repeated visits to the same plots.
species_fits <- lapply(names(species_params), function(sp) {
  sp_data <- measured[measured$SPECIES == sp, ]
  n_plots <- length(unique(sp_data$PLOT_ID))
  fit     <- NULL
  method  <- "none"

  # Mixed-effects preferred when enough plots available
  if (n_plots >= 5 && requireNamespace("nlme", quietly = TRUE)) {
    fit <- tryCatch(
      fit_hd_model(sp_data, model = "naslund", group_col = "PLOT_ID"),
      error = function(e) NULL
    )
    if (!is.null(fit)) method <- "mixed-effects"
  }
  # Fixed-effects fallback
  if (is.null(fit)) {
    fit <- tryCatch(
      fit_hd_model(sp_data, model = "naslund"),
      error = function(e) NULL
    )
    if (!is.null(fit)) method <- "fixed-effects"
  }
  list(species = sp, fit = fit, method = method)
})

# Apply predictions — record source for audit trail
sim_trees$HEIGHT_PRED   <- sim_trees$HEIGHT
sim_trees$HEIGHT_METHOD <- ifelse(!is.na(sim_trees$HEIGHT), "measured", NA_character_)

for (sf in species_fits) {
  if (is.null(sf$fit)) next
  cf      <- sf$fit$coefficients
  missing <- is.na(sim_trees$HEIGHT) & sim_trees$SPECIES == sf$species
  if (any(missing)) {
    sim_trees$HEIGHT_PRED[missing] <- round(
      ht_from_dbh(sim_trees$DBH[missing],
                  model = "naslund",
                  a     = cf["a"],
                  b     = cf["b"]),
      1
    )
    sim_trees$HEIGHT_METHOD[missing] <- paste0("naslund_", sf$method)
  }
}

n_filled <- sum(is.na(sim_trees$HEIGHT) & !is.na(sim_trees$HEIGHT_PRED))
cat(sprintf("  Heights predicted for %d previously unmeasured trees\n", n_filled))
cat(sprintf("  Trees still missing (no species model): %d\n",
            sum(is.na(sim_trees$HEIGHT_PRED))))

# Audit: count by source
cat("\n  Height source breakdown:\n")
tbl <- table(sim_trees$HEIGHT_METHOD, useNA = "ifany")
for (nm in names(tbl)) {
  cat(sprintf("    %-30s  %d\n", nm, tbl[[nm]]))
}

# Quick sanity check on predictions
predicted_only <- sim_trees[is.na(sim_trees$HEIGHT) & !is.na(sim_trees$HEIGHT_PRED), ]
cat("\n  Sample predictions (first 10 unmeasured trees):\n")
cat(sprintf("  %-14s  %-7s  %6s  %10s\n",
            "TREE_ID", "SPECIES", "DBH", "PRED_HT"))
cat(sprintf("  %-14s  %-7s  %6s  %10s\n",
            "--------------", "-------", "------", "----------"))
for (i in seq_len(min(10, nrow(predicted_only)))) {
  r <- predicted_only[i, ]
  cat(sprintf("  %-14s  %-7s  %6.1f  %10.1f\n",
              r$TREE_ID, r$SPECIES, r$DBH, r$HEIGHT_PRED))
}

cat("\nDone.\n")
