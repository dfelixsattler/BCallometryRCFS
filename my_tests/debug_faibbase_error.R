library(data.table)
library(fpCompare)
devtools::load_all("C:/FAIBBase", quiet = TRUE)

cat("=== Testing failing tree: BWBS H DBH=10 ht=14.4 ===\n")

withCallingHandlers(
  treeProfile(taperEquationForm = "KBEC", FIZorBEC = "BWBS", species = "H",
              height = 14.4, DBH = 10, stumpHeight = 0.3, UTOPDIB = 10,
              BTOPHeight = NA),
  error = function(e) {
    cat("Error message:", conditionMessage(e), "\n")
    calls <- sys.calls()
    cat("Call stack:\n")
    for (i in seq_along(calls)) {
      cat(sprintf("  [%d] %s\n", i, deparse(calls[[i]])[1]))
    }
  }
)

cat("\n=== Step-by-step trace ===\n")

# Replicate the internals manually
taperEquationForm <- "KBEC"
FIZorBEC <- "BWBS"
species <- "H"
height <- 14.4
DBH <- 10
stumpHeight <- 0.3
breastHeight <- 1.3
UTOPDIB <- 10
BTOPHeight <- NA

cat("1. DIB_stump:\n")
dib_stump <- DIB_ICalculator(taperEquationForm, FIZorBEC = FIZorBEC, species = species,
                              height_I = stumpHeight, heightTotal = height, DBH = DBH,
                              volMultiplier = 1)
cat("  DIB_stump =", dib_stump, "\n")

cat("2. DIB_BH:\n")
dib_bh <- DIB_ICalculator(taperEquationForm, FIZorBEC = FIZorBEC, species = species,
                           height_I = breastHeight, heightTotal = height, DBH = DBH,
                           volMultiplier = 1)
cat("  DIB_BH =", dib_bh, "  (length =", length(dib_bh), ")\n")

cat("3. Profile heights:\n")
ht_seq <- seq(stumpHeight, height, by = 0.1)
cat("  n heights =", length(ht_seq), "  range:", min(ht_seq), "-", max(ht_seq), "\n")

cat("4. Vector DIB_ICalculator call:\n")
dib_vec <- DIB_ICalculator(taperEquationForm, FIZorBEC = FIZorBEC, species = species,
                            height_I = ht_seq, heightTotal = height, DBH = DBH,
                            volMultiplier = 1)
cat("  length(dib_vec) =", length(dib_vec), "  (vs n heights =", length(ht_seq), ")\n")
cat("  any NA:", any(is.na(dib_vec)), "  any Inf:", any(is.infinite(dib_vec)), "\n")

cat("5. BH floor check:\n")
below_bh <- ht_seq < breastHeight
cat("  rows below BH:", sum(below_bh), "\n")
cat("  DIBs below BH < DIB_BH:", sum(dib_vec[below_bh] < dib_bh), "\n")
