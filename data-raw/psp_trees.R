# =============================================================================
# BCallometryR — psp_trees dataset
#
# Simulated permanent-sample-plot (PSP) tree list for use in package
# vignettes and examples.  Mimics the structure of a real BC PSP database:
# raw inventory species codes, BEC zone, DBH for all trees, and measured
# heights for roughly 25% of trees (the rest are NA, as in real PSPs where
# only a subsample of trees is height-measured each visit).
#
# All five species have taper equations in taper_coefs_kbec and biomass
# equations in biomass_coefs, so the dataset supports every vignette workflow
# (H-D modelling, biomass estimation, volume estimation).
#
# Columns
#   SITE_IDENTIFIER  character  Plot identifier, e.g. "PSP-001"
#   BEC_ZONE         character  BEC zone of the plot
#   unitreeid        character  Unique tree identifier, e.g. "PSP-001-01"
#   SPECIES          character  Raw BC inventory species code (PSP/VRI)
#   DBH              numeric    Diameter at breast height (cm)
#   HEIGHT           numeric    Measured total height (m); NA for unmeasured,
#                               broken-top, or dead trees
#   HEIGHT_TRUE      numeric    True simulated total height (m); for validation
#   BTOP             logical    TRUE if the tree has a broken top
#   LV_D             character  Tree status: "L" = live, "D" = dead (standing)
# =============================================================================

library(BCallometryR)   # for ht_naslund
library(usethis)

set.seed(42)

# ---------------------------------------------------------------------------
# Simulation parameters
# ---------------------------------------------------------------------------

# Naslund H-D parameters keyed by PSP species code
sp_params <- list(
  PLI = list(a = 1.85, b = 0.130),   # lodgepole pine
  SW  = list(a = 1.65, b = 0.105),   # white spruce
  HW  = list(a = 1.45, b = 0.100),   # western hemlock
  FDI = list(a = 1.60, b = 0.110),   # interior Douglas-fir
  AT  = list(a = 1.55, b = 0.115)    # trembling aspen
)

# Species that commonly occur in each BEC zone
sp_by_bec <- list(
  SBS = c("PLI", "SW",  "AT"),
  ICH = c("PLI", "HW",  "FDI"),
  CWH = c("HW",  "FDI", "SW"),
  IDF = c("FDI", "PLI", "AT")
)

n_plots <- 20   # 5 plots per BEC zone

plot_meta <- data.frame(
  SITE_IDENTIFIER = paste0("PSP-", formatC(seq_len(n_plots), width = 3, flag = "0")),
  BEC_ZONE        = rep(names(sp_by_bec), each = 5),
  a_shift         = rnorm(n_plots, mean = 0, sd = 0.22),
  stringsAsFactors = FALSE
)

# ---------------------------------------------------------------------------
# Generate tree records
# ---------------------------------------------------------------------------

psp_trees <- do.call(rbind, lapply(seq_len(n_plots), function(i) {
  bec     <- plot_meta$BEC_ZONE[i]
  sps     <- sp_by_bec[[bec]]
  n_trees <- sample(20:35, 1)

  sp  <- sample(sps, n_trees, replace = TRUE)
  dbh <- round(runif(n_trees, min = 5, max = 65), 1)

  # True heights from Naslund model with plot random effect + noise
  ht_true <- mapply(function(d, s) {
    p   <- sp_params[[s]]
    eff <- p$a + plot_meta$a_shift[i]
    ht  <- ht_naslund(d, a = eff, b = p$b) + rnorm(1, 0, 1.0)
    round(pmax(ht, 1.5), 1)
  }, dbh, sp)

  # ~25% of trees height-measured; rest are NA (typical PSP protocol)
  unmeasured          <- sample(n_trees, size = round(n_trees * 0.75))
  ht_obs              <- ht_true
  ht_obs[unmeasured]  <- NA

  data.frame(
    SITE_IDENTIFIER = plot_meta$SITE_IDENTIFIER[i],
    BEC_ZONE        = bec,
    unitreeid       = paste0(plot_meta$SITE_IDENTIFIER[i], "-",
                             formatC(seq_len(n_trees), width = 2, flag = "0")),
    SPECIES         = sp,
    DBH             = dbh,
    HEIGHT          = ht_obs,
    HEIGHT_TRUE     = ht_true,
    BTOP            = FALSE,
    stringsAsFactors = FALSE
  )
}))

# Reset row names
rownames(psp_trees) <- NULL

# ---------------------------------------------------------------------------
# Add broken-top trees (~5% of trees, distributed across plots)
# Broken-top trees cannot have total height measured or imputed;
# HEIGHT is set to NA regardless of measurement status.
# ---------------------------------------------------------------------------

set.seed(99)
btop_idx <- sample(nrow(psp_trees), size = round(nrow(psp_trees) * 0.05))
psp_trees$BTOP[btop_idx]   <- TRUE
psp_trees$HEIGHT[btop_idx] <- NA   # not measurable as total height

# ---------------------------------------------------------------------------
# Add live/dead status (~12% dead, standing dead snags)
# Dead trees: DBH still measured, but HEIGHT is NA (not measured after death).
# Dead trees are excluded from H-D model fitting (matching FAIBCompiler's
# TREE_EXTANT_CODE == "L" filter) but retained for biomass and volume.
# Dead trees do not overlap with broken-top trees (kept separate for clarity).
# ---------------------------------------------------------------------------

set.seed(123)
psp_trees$LV_D <- "L"
live_non_btop  <- which(!psp_trees$BTOP)
dead_idx       <- sample(live_non_btop, size = round(length(live_non_btop) * 0.12))
psp_trees$LV_D[dead_idx]    <- "D"
psp_trees$HEIGHT[dead_idx]  <- NA   # dead trees not height-measured

# ---------------------------------------------------------------------------
# Save
# ---------------------------------------------------------------------------

usethis::use_data(psp_trees, overwrite = TRUE)
