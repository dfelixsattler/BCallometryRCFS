# BCallometryRCFS

An R package for tree-level allometric calculations from BC permanent sample
plot (PSP) and non-PSP inventory data: species code crosswalk, height–diameter
modelling, height imputation, biomass estimation, and stem volume.

## Installation

```r
# Install from GitHub
remotes::install_github("your-org/BCallometryRCFS")
```

## Data pipeline

The diagram below shows the end-to-end workflow from raw field data to
estimated heights, volumes, and biomass.

```mermaid
flowchart TD
    A[("Field data\nDBH · HEIGHT · SPECIES · BEC_ZONE · LV_D")]

    A --> SC

    subgraph SC["1 · Species code crosswalk"]
        direction TB
        SC1["species_correction(species, bec_zone)\ntranslate & BEC-disambiguate raw codes"]
        SC2["bc_species_to_sp0()\nSP0 group code  ·  e.g. F, H, S, PL"]
        SC3["bc_species_to_sp_type()\nSP_TYPE  ·  C (conifer) or D (deciduous)"]
        SC4["bc_species_to_biomass_name()\ncommon name  ·  e.g. douglas-fir"]
        SC1 --> SC2
        SC1 --> SC3
        SC1 --> SC4
    end

    SC2 --> HD

    subgraph HD["2 · H-D modelling  (alive, measured trees only)"]
        direction TB
        HD1["fit_hd_models_by_group(measured)\nNaslund / Curtis · mixed-effects or fixed-effects\nSP0 level → SP_TYPE fallback"]
        HD2["ht_impute(all_trees, hd_result, impute_btop = TRUE)\nHT_PROJ for every tree\nBLUP-calibrated per plot"]
        HD1 --> HD2
    end

    SC4 --> BIO
    HD2 --> VOL
    HD2 --> BIO

    subgraph VOL["3 · Volume"]
        VOL1["tree_volume(bec_zone, sp0, dbh, ht_proj)\nWSV · MER · NMR · STUMP  (m³)\nKozak KBEC or KFIZ3 taper equation"]
    end

    subgraph BIO["4 · Biomass"]
        BIO1["biomass_tree(species_name, dbh, height)\nTotal aboveground dry biomass (kg)\nLambert 2005 / Ung 2008 national equations"]
        BIO2["biomass_components()\nwood · bark · branches · foliage (kg)"]
    end
```

> **Notes**
> - Dead trees (`LV_D == "D"`) and broken-top trees (`BTOP == TRUE`) are
>   excluded from H-D model training but flow through to volume and biomass.
> - `impute_btop = TRUE` estimates total height for broken-top trees from DBH,
>   matching the BC MoF compilation routine.
> - Use `taper_eq = "KFIZ3"` with FIZ zone codes for pre-BEC inventory data.

## Quick start

```r
library(BCallometryRCFS)

trees <- psp_trees

# 1. Species crosswalk
trees$SPECIES_CORR    <- species_correction(trees$SPECIES, trees$BEC_ZONE)
trees$SPECIES_SP0     <- bc_species_to_sp0(trees$SPECIES_CORR)
trees$SPECIES_SP_TYPE <- bc_species_to_sp_type(trees$SPECIES_CORR)
trees$SPECIES_NAME    <- bc_species_to_biomass_name(trees$SPECIES, trees$BEC_ZONE)

# 2. Fit H-D models and impute heights
measured  <- trees[!is.na(trees$HEIGHT) & trees$LV_D == "L" & !trees$BTOP, ]
hd_result <- fit_hd_models_by_group(measured)
trees     <- ht_impute(trees, hd_result, impute_btop = TRUE)

# 3. Volume
trees$WSV_M3 <- tree_volume(trees$BEC_ZONE, trees$SPECIES_SP0,
                             trees$DBH, trees$HT_PROJ)

# 4. Biomass
trees$BIOMASS_KG <- biomass_tree(trees$SPECIES_NAME, trees$DBH,
                                  height = trees$HT_PROJ)
```

## Vignettes

| Vignette | Topic |
|---|---|
| `vignette("full_psp_workflow")` | End-to-end PSP workflow |
| `vignette("hd_psp_workflow")` | H-D modelling in depth |
| `vignette("biomass_psp_workflow")` | Biomass equations |
| `vignette("volume_psp_workflow")` | Kozak taper volume |

## Citation

```r
citation("BCallometryRCFS")
```
