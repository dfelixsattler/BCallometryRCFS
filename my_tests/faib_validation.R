# =============================================================================
# BCallometryR vs FAIBCompiler validation — multi-site coastal PSP dataset
# Trees: faib_tree_detail_subtest.csv  |  BEC info: faib_plot_header.csv
# =============================================================================
devtools::load_all("C:/BCallometryR", quiet = TRUE)
suppressPackageStartupMessages(library(data.table))

# -----------------------------------------------------------------------------
# 0. Load and join
# -----------------------------------------------------------------------------
trees <- fread("my_tests/faib_tree_detail_subtest.csv")
plots <- fread("my_tests/faib_plot_header.csv")

trees <- merge(trees,
               plots[, .(SITE_IDENTIFIER, BEC_ZONE, BEC_SBZ)],
               by = "SITE_IDENTIFIER", all.x = TRUE)

trees[, DBH     := as.numeric(DBH)]
trees[, HEIGHT  := as.numeric(HEIGHT)]
trees[, VOL_WSV := as.numeric(VOL_WSV)]
trees[, VOL_MER := as.numeric(VOL_MER)]
trees[, BTOP    := BROKEN_TOP_IND == "Y"]
trees[, MEAS    := HEIGHT_SOURCE == "Field measured"]

cat(sprintf(
  "Loaded %d trees | %d sites | measured: %d | est-DBH: %d | btop: %d | no-HT-source: %d\n\n",
  nrow(trees),
  length(unique(trees$SITE_IDENTIFIER)),
  sum(trees$MEAS, na.rm = TRUE),
  sum(trees$HEIGHT_SOURCE == "Estimated based on DBH", na.rm = TRUE),
  sum(trees$BTOP, na.rm = TRUE),
  sum(trees$HEIGHT_SOURCE == "", na.rm = TRUE)))

# =============================================================================
# 1. SPECIES CROSSWALK
# =============================================================================
cat("=== 1. Species crosswalk ===\n")
trees[, SPECIES_CORR    := species_correction(SPECIES, BEC_ZONE)]
trees[, SPECIES_SP0     := bc_species_to_sp0(SPECIES_CORR)]
trees[, SPECIES_SP_TYPE := bc_species_to_sp_type(SPECIES_CORR)]
trees[, SPECIES_NAME    := bc_species_to_biomass_name(SPECIES, BEC_ZONE)]

xwalk <- unique(trees[, .(SPECIES, SPECIES_CORR, SPECIES_SP0,
                           SPECIES_SP_TYPE, SPECIES_NAME)])[order(SPECIES)]
print(xwalk, row.names = FALSE)

unmapped <- trees[is.na(SPECIES_SP0), .N, SPECIES][order(-N)]
if (nrow(unmapped) > 0) {
  cat("\nSpecies with no SP0 mapping (excluded from volume/biomass):\n")
  print(unmapped, row.names = FALSE)
}
cat("\n")

# =============================================================================
# 2. VOLUME — using FAIB's own HEIGHT values (tests taper equation only)
#    Restricted to field-measured trees with known SP0 and BEC zone.
# =============================================================================
cat("=== 2. Volume (FAIB heights, field-measured trees — tests taper equation) ===\n")

vol_test <- trees[MEAS == TRUE & !is.na(SPECIES_SP0) & !is.na(BEC_ZONE) &
                  !is.na(DBH) & !is.na(HEIGHT) & DBH > 1 & HEIGHT > 1.4]

vol_test[, VOL_WSV_BC := tree_volume(BEC_ZONE, SPECIES_SP0, DBH, HEIGHT,
                                      volume_type = "WSV")]
vol_test[, VOL_MER_BC := tree_volume(BEC_ZONE, SPECIES_SP0, DBH, HEIGHT,
                                      volume_type = "MER")]

vol_test[, WSV_diff := VOL_WSV_BC - VOL_WSV]
vol_test[, MER_diff := VOL_MER_BC - VOL_MER]
vol_test[, WSV_pct  := WSV_diff / VOL_WSV * 100]
vol_test[, MER_pct  := MER_diff / VOL_MER * 100]

cat(sprintf("Trees compared: %d\n", nrow(vol_test)))
cat(sprintf("WSV — mean abs diff: %.6f m3  |  RMSE: %.6f m3  |  mean abs %%: %.3f%%\n",
            mean(abs(vol_test$WSV_diff), na.rm = TRUE),
            sqrt(mean(vol_test$WSV_diff^2, na.rm = TRUE)),
            mean(abs(vol_test$WSV_pct), na.rm = TRUE)))
cat(sprintf("MER — mean abs diff: %.6f m3  |  RMSE: %.6f m3  |  mean abs %%: %.3f%%\n\n",
            mean(abs(vol_test$MER_diff), na.rm = TRUE),
            sqrt(mean(vol_test$MER_diff^2, na.rm = TRUE)),
            mean(abs(vol_test$MER_pct[is.finite(vol_test$MER_pct)]), na.rm = TRUE)))

wsv_sp <- vol_test[, .(n       = .N,
                        WSV_mae = round(mean(abs(WSV_diff), na.rm = TRUE), 6),
                        WSV_pct = round(mean(abs(WSV_pct),  na.rm = TRUE), 3)),
                   by = .(SPECIES, SPECIES_SP0)][order(SPECIES)]
cat("WSV by species:\n")
print(wsv_sp, row.names = FALSE)
cat("\n")

# =============================================================================
# 3. HEIGHT IMPUTATION
#    Blank FAIB-estimated heights, then re-impute with BCallometryR H-D models
#    fit on all field-measured trees across all sites in this dataset.
# =============================================================================
cat("=== 3. Height imputation (BCallometryR vs FAIB) ===\n")
cat("Blanking FAIB-estimated heights and re-imputing.\n")
cat("BCallometryR trains on field-measured trees from this dataset;\n")
cat("FAIB uses its full provincial PSP network - differences expected.\n\n")

trees_blank <- copy(trees)
trees_blank[HEIGHT_SOURCE == "Estimated based on DBH", HEIGHT := NA]

measured_df <- as.data.frame(
  trees_blank[MEAS == TRUE & !is.na(HEIGHT) & !BTOP &
              !is.na(SPECIES_SP0) & !is.na(DBH) & DBH > 0 & HEIGHT > 1.3])

cat(sprintf("Training: %d measured trees across %d sites\n",
            nrow(measured_df), length(unique(measured_df$SITE_IDENTIFIER))))

hd_result <- fit_hd_models_by_group(measured_df)
cat("\nH-D model summary:\n")
print(hd_result$summary, row.names = FALSE)
cat("\n")

trees_imp <- as.data.table(ht_impute(as.data.frame(trees_blank),
                                     hd_result, impute_btop = TRUE))

e <- trees$HEIGHT_SOURCE == "Estimated based on DBH" &
     !is.na(trees$HEIGHT) & !is.na(trees_imp$HT_PROJ)

ht_diff <- trees_imp$HT_PROJ[e] - trees$HEIGHT[e]
cat(sprintf("Trees compared: %d\n", sum(e, na.rm = TRUE)))
cat(sprintf("Mean diff (m): %+.2f  |  Mean abs diff: %.2f  |  RMSE: %.2f m\n\n",
            mean(ht_diff, na.rm = TRUE),
            mean(abs(ht_diff), na.rm = TRUE),
            sqrt(mean(ht_diff^2, na.rm = TRUE))))

ht_sp <- data.table(SPECIES     = trees$SPECIES[e],
                    SPECIES_SP0 = trees$SPECIES_SP0[e],
                    HT_FAIB     = trees$HEIGHT[e],
                    diff        = ht_diff
)[, .(n       = .N,
      HT_mean = round(mean(HT_FAIB, na.rm = TRUE), 1),
      bias    = round(mean(diff,     na.rm = TRUE), 2),
      mae     = round(mean(abs(diff),na.rm = TRUE), 2),
      rmse    = round(sqrt(mean(diff^2, na.rm = TRUE)), 2)),
  by = .(SPECIES, SPECIES_SP0)][order(SPECIES)]
cat("Height comparison by species:\n")
print(ht_sp, row.names = FALSE)
cat("\n")

# =============================================================================
# 4. VOLUME — full BCallometryR pipeline (BC-imputed heights for est. trees)
# =============================================================================
cat("=== 4. Volume — full BC pipeline (BC-imputed heights) ===\n")

trees_imp[, VOL_WSV_BC_full := tree_volume(BEC_ZONE, SPECIES_SP0, DBH,
                                             HT_PROJ, volume_type = "WSV")]

ev <- e & !is.na(trees_imp$VOL_WSV_BC_full) &
      !is.na(trees$VOL_WSV) & trees$VOL_WSV > 0

wsv_diff <- trees_imp$VOL_WSV_BC_full[ev] - trees$VOL_WSV[ev]
wsv_pct  <- wsv_diff / trees$VOL_WSV[ev] * 100

cat(sprintf("Trees compared: %d\n", sum(ev, na.rm = TRUE)))
cat(sprintf("Mean abs diff: %.5f m3  |  RMSE: %.5f m3  |  Mean abs %%: %.1f%%\n\n",
            mean(abs(wsv_diff), na.rm = TRUE),
            sqrt(mean(wsv_diff^2, na.rm = TRUE)),
            mean(abs(wsv_pct), na.rm = TRUE)))

vol_sp <- data.table(SPECIES     = trees$SPECIES[ev],
                     SPECIES_SP0 = trees$SPECIES_SP0[ev],
                     WSV_diff    = wsv_diff,
                     WSV_pct     = wsv_pct
)[, .(n       = .N,
      mae_m3  = round(mean(abs(WSV_diff), na.rm = TRUE), 5),
      rmse_m3 = round(sqrt(mean(WSV_diff^2, na.rm = TRUE)), 5),
      mae_pct = round(mean(abs(WSV_pct), na.rm = TRUE), 1)),
  by = .(SPECIES, SPECIES_SP0)][order(SPECIES)]
cat("Volume (BC pipeline) by species:\n")
print(vol_sp, row.names = FALSE)
cat("\n")

# =============================================================================
# 5. BIOMASS — BCallometryR only (no FAIB reference in this file)
# =============================================================================
cat("=== 5. Biomass (BCallometryR, height-included) ===\n")
has_sp <- !is.na(trees_imp$SPECIES_NAME)
trees_imp[, BIOMASS_KG := NA_real_]
trees_imp[has_sp, BIOMASS_KG := biomass_tree(SPECIES_NAME, DBH,
                                               height = HT_PROJ)]

bio_sp <- trees_imp[has_sp & !is.na(BIOMASS_KG),
                    .(n       = .N,
                      mean_kg = round(mean(BIOMASS_KG), 1),
                      sd_kg   = round(sd(BIOMASS_KG),   1)),
                    by = .(SPECIES, SPECIES_SP0)][order(SPECIES)]
cat("Biomass by species (kg/tree):\n")
print(bio_sp, row.names = FALSE)
