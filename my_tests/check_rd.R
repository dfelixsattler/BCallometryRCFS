library(BCallometryR)

# Test a few different httpd paths RStudio might request
paths <- c(
  "/library/BCallometryR/html/BCallometryR-package.html",
  "/library/BCallometryR/html/00Index.html",
  "/library/BCallometryR/doc/index.html"
)

for (p in paths) {
  cat("\n--- Path:", p, "---\n")
  tryCatch({
    res <- tools:::httpd(p, NULL, NULL)
    cat("Response names:", paste(names(res), collapse=", "), "\n")
    if (is.character(res$payload))
      cat("Payload (first 200):", substr(res$payload, 1, 200), "\n")
    else
      cat("Payload type:", class(res$payload), "\n")
  }, error = function(e) cat("ERROR:", conditionMessage(e), "\n"))
}
