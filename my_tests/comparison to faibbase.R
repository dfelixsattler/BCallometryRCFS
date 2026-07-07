setwd("C:/BCallometryR")
library(devtools)
load_all("C:/BCallometryR")
load_all("C:/FAIBBase")

source("C:/BCallometryR/my_tests/validate_vs_faibbase.R")

install.packages(c("bcmaps", "raster"), lib = "C:/Rlibs",
                 repos = "https://cran.rstudio.com")

devtools::load_all("C:/FAIBBase", quiet = TRUE)
devtools::load_all("C:/BCallometryR", quiet = TRUE)

# intact tree
treeProfile("KBEC", "CWH", "H", height=27.4, DBH=30.7)$volume_summary

# broken top
treeProfile("KBEC", "CWH", "H", height=27.4, DBH=30.7, BTOPHeight=15)$volume_summary
