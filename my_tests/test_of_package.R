setwd("C:/BCallometryR")

# build the r package
library(devtools)

remove.packages("BCallometryR")   # removes from the first writable lib
devtools::install("C:/BCallometryR", build_vignettes = TRUE)

# when the tarball has been build, run this in place of the above
# remove.packages("BCallometryR")
# devtools::install_github("dfelixsattler/BCallometryR")
# install.packages("C:/BCallometryR/BCallometryR_0.1.0.tar.gz", repos = NULL, type = "source")

library(BCallometryR)
?BCallometryR

###########################################################
# Full Workflow
###########################################################
trees <- psp_trees

# examine the data
cat(sprintf(
  "%d plots | %d trees | %d heights measured (%.0f%%) | %d missing\n",
  length(unique(trees$SITE_IDENTIFIER)), nrow(trees),
  sum(!is.na(trees$HEIGHT)), 100 * mean(!is.na(trees$HEIGHT)),
  sum(is.na(trees$HEIGHT))
))
#> 20 plots | 588 trees | 142 heights measured (24%) | 446 missing

# species cross walk
trees$SPECIES_CORR    <- species_correction(trees$SPECIES, trees$BEC_ZONE)
trees$SPECIES_SP0     <- bc_species_to_sp0(trees$SPECIES_CORR)
trees$SPECIES_SP_TYPE <- bc_species_to_sp_type(trees$SPECIES_CORR)  # "C" or "D"

trees$SPECIES_NAME    <- bc_species_to_biomass_name(trees$SPECIES, trees$BEC_ZONE)

# Show unique crosswalk mappings
xwalk <- unique(trees[, c("SPECIES", "BEC_ZONE", "SPECIES_CORR",
                          "SPECIES_SP0", "SPECIES_SP_TYPE", "SPECIES_NAME")])
print(xwalk[order(xwalk$BEC_ZONE, xwalk$SPECIES), ], row.names = FALSE)

# fit the HD model by group
measured  <- trees[!is.na(trees$HEIGHT) & !trees$BTOP, ]
hd_result <- fit_hd_models_by_group(measured)
print(hd_result$summary, row.names = FALSE)

# impute heights for trees with missing
trees <- ht_impute(trees, hd_result)

cat(sprintf("Heights filled: %d  |  Still missing: %d\n",
            sum(!is.na(trees$HT_PROJ) & is.na(trees$HEIGHT)),
            sum(is.na(trees$HT_PROJ))))

# no for biomass
trees$BIOMASS_KG <- biomass_tree(
  species = trees$SPECIES_NAME,
  dbh     = trees$DBH,
  height  = trees$HT_PROJ
)

comp <- biomass_components(
  species = trees$SPECIES_NAME,
  dbh     = trees$DBH,
  height  = trees$HT_PROJ
)

trees$BIO_WOOD_KG    <- comp$wood
trees$BIO_BARK_KG    <- comp$bark
trees$BIO_BRANCH_KG  <- comp$branches
trees$BIO_FOLIAGE_KG <- comp$foliage


bio_summary <- do.call(rbind, lapply(sp0_list, function(sp0) {
  sub <- trees[!is.na(trees$SPECIES_SP0) & trees$SPECIES_SP0 == sp0 &
                 !is.na(trees$BIOMASS_KG), ]
  data.frame(SP0      = sp0,
             n        = nrow(sub),
             mean_kg  = round(mean(sub$BIOMASS_KG), 1),
             sd_kg    = round(sd(sub$BIOMASS_KG),   1))
}))
