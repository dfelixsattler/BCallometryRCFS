# =============================================================================
# BCallometryR — plant_codes dataset
#
# Builds a cross-reference table mapping OSM PlantCodes (USDA-style codes
# used in the OSM.Allometry .NET library) to:
#   - BC SP0 codes  (used in BCallometryR::taper_coefs_kbec / tree_volume)
#   - common names  (used in BCallometryR::biomass_coefs / biomass_tree)
#
# Sources:
#   OSM PlantCodes enum (OSM.CommonModels, C. Hennigar, NRCan)
#   USDA PLANTS database: https://plants.usda.gov
#   BC Ministry of Forests SP0 species codes (Timber Supply Review conventions)
#   Lambert et al. (2005) and Ung et al. (2008) species lists
# =============================================================================

library(usethis)

# Helper to build one row
.row <- function(plant_code, latin_name, common_name, sp0) {
  data.frame(plant_code  = plant_code,
             latin_name  = latin_name,
             common_name = common_name,
             sp0         = sp0,
             stringsAsFactors = FALSE)
}
NA_ch <- NA_character_

plant_codes <- rbind(

  # ---- Species in taper_coefs_kbec (BC SP0 codes) -------------------------
  # AC — Black cottonwood (primary BC species behind SP0 "AC")
  .row("POBAT",  "Populus trichocarpa",             "black cottonwood",   "AC"),
  # AC — Balsam poplar (interior BC; also compiles as "AC" in volume system)
  .row("POBA2",  "Populus balsamifera",             "balsam poplar",      "AC"),

  # AT — Trembling aspen
  .row("POTR5",  "Populus tremuloides",             "trembling aspen",    "AT"),

  # B — True firs (Abies spp.); all BC Abies species compile as SP0 "B"
  # Amabilis fir is the primary coastal species (BA → "B")
  .row("ABAM",   "Abies amabilis",                  "pacific silver fir", "B"),
  # Subalpine fir is the primary interior species (BL → "B")
  .row("ABLA",   "Abies lasiocarpa",                "subalpine fir",      "B"),
  # Balsam fir (BB → "B"; mainly eastern Canada but included in BC compiler)
  .row("ABBA",   "Abies balsamea",                  "balsam fir",         "B"),

  # C — Western redcedar
  .row("THPL",   "Thuja plicata",                   "western redcedar",   "C"),

  # D — Red alder (SP0 "D" = alder group in BC volume compilation)
  .row("ALRU2",  "Alnus rubra",                     "red alder",          "D"),

  # E — Paper birch (SP0 "E" = birch group)
  .row("BEPA",   "Betula papyrifera",               "white birch",        "E"),

  # F — Douglas-fir (SP0 "F" = Douglas-fir group; FD, FDI, FDC all → "F")
  .row("PSME",   "Pseudotsuga menziesii",           "douglas-fir",        "F"),

  # H — Western hemlock
  .row("TSHE",   "Tsuga heterophylla",              "western hemlock",    "H"),

  # L — Western larch (no Lambert/Ung biomass equations for this species)
  .row("LAOC",   "Larix occidentalis",              NA_ch,                "L"),

  # MB — Bigleaf maple (no Lambert/Ung biomass equations)
  .row("ACMA3",  "Acer macrophyllum",               NA_ch,                "MB"),

  # PA — Whitebark pine (no Lambert/Ung biomass equations)
  .row("PIAL",   "Pinus albicaulis",                NA_ch,                "PA"),

  # PL — Lodgepole pine
  .row("PICO",   "Pinus contorta",                  "lodgepole pine",     "PL"),

  # PW — Western white pine (no Lambert/Ung biomass equations)
  .row("PIMO3",  "Pinus monticola",                 NA_ch,                "PW"),

  # PY — Ponderosa pine (no Lambert/Ung biomass equations)
  .row("PIPO",   "Pinus ponderosa",                 NA_ch,                "PY"),

  # S — Spruce composite: three species all compile as "S" in BC
  # White spruce is the primary interior mapping
  .row("PIGL",   "Picea glauca",                    "white spruce",       "S"),
  .row("PISI",   "Picea sitchensis",                "sitka spruce",       "S"),
  .row("PIEN",   "Picea engelmannii",               "engelmann spruce",   "S"),

  # Y — Alaska yellow-cedar (no Lambert/Ung biomass equations)
  .row("CHNO",   "Chamaecyparis nootkatensis",      NA_ch,                "Y"),

  # ---- Species in biomass_coefs only (no BC taper equation) ---------------

  .row("ABBI3",  "Abies bifolia",                   "alpine fir",         NA_ch),
  .row("TIAM",   "Tilia americana",                 "basswood",           NA_ch),
  .row("FRNI",   "Fraxinus nigra",                  "black ash",          NA_ch),
  .row("PRSE2",  "Prunus serotina",                 "black cherry",       NA_ch),
  .row("PIMA",   "Picea mariana",                   "black spruce",       NA_ch),
  .row("TSCA",   "Tsuga canadensis",                "eastern hemlock",    NA_ch),
  .row("JUVI",   "Juniperus virginiana",            "eastern redcedar",   NA_ch),
  .row("THOC2",  "Thuja occidentalis",              "eastern white-cedar",NA_ch),
  .row("PIST",   "Pinus strobus",                   "eastern white pine", NA_ch),
  .row("BEPO",   "Betula populifolia",              "grey birch",         NA_ch),
  .row("CARYA",  "Carya spp.",                      "hickory",            NA_ch),
  .row("OSVI",   "Ostrya virginiana",               "hop-hornbeam",       NA_ch),
  .row("PIBA2",  "Pinus banksiana",                 "jack pine",          NA_ch),
  .row("POGR4",  "Populus grandidentata",           "largetooth aspen",   NA_ch),
  .row("FRPE",   "Fraxinus pennsylvanica",          "red ash",            NA_ch),
  .row("ACRU",   "Acer rubrum",                     "red maple",          NA_ch),
  .row("QURU",   "Quercus rubra",                   "red oak",            NA_ch),
  .row("PIRE",   "Pinus resinosa",                  "red pine",           NA_ch),
  .row("PIRU",   "Picea rubens",                    "red spruce",         NA_ch),
  .row("ACSA2",  "Acer saccharinum",                "silver maple",       NA_ch),
  .row("ACSA3",  "Acer saccharum",                  "sugar maple",        NA_ch),
  .row("LALA",   "Larix laricina",                  "tamarack larch",     NA_ch),
  .row("ULAM",   "Ulmus americana",                 "white elm",          NA_ch),
  .row("QUAL",   "Quercus alba",                    "white oak",          NA_ch),
  .row("BEAL2",  "Betula alleghaniensis",           "yellow birch",       NA_ch),

  # Generic / unknown categories
  .row("T_B",    "Unknown hardwood",                "hardwood",           NA_ch),
  .row("T_C",    "Unknown softwood",                "softwood",           NA_ch),

  # Biomass-only (DBH + height equations only)
  .row("FAGR",   "Fagus grandifolia",               "beech",              NA_ch)
)

rownames(plant_codes) <- NULL

cat(sprintf("plant_codes: %d rows, %d species with sp0, %d species with common_name\n",
            nrow(plant_codes),
            sum(!is.na(plant_codes$sp0)),
            sum(!is.na(plant_codes$common_name))))

usethis::use_data(plant_codes, overwrite = TRUE)
