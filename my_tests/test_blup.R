library(BCallometryR)
set.seed(42)
n_plots <- 12; n_per <- 25
plot_ids <- rep(paste0("P", 1:n_plots), each = n_per)
dbh    <- runif(n_plots * n_per, 5, 60)
a_true <- rep(rnorm(n_plots, 1.5, 0.3), each = n_per)
height <- ht_naslund(dbh, a = a_true, b = 0.05) + rnorm(length(dbh), 0, 1)
train  <- data.frame(PLOT_ID = plot_ids, DBH = dbh, HEIGHT = height)
fit    <- fit_hd_model(train, model = "naslund", group_col = "PLOT_ID")
cat("Fixed-effect a:", round(fit$coefficients["a"], 4), "\n")

# New plot: true a = 1.9 (well above average)
calib <- data.frame(
  DBH    = c(15, 28, 42),
  HEIGHT = ht_naslund(c(15, 28, 42), a = 1.9, b = 0.05) + rnorm(3, 0, 0.3)
)
res <- calibrate_hd_blup(fit, calib)
cat("Calibrated a   :", round(res$coefficients["a"], 4), "\n")
cat("u_hat          :", round(res$u_hat, 4), "(expected ~+0.4)\n")
cat("n_calib        :", res$n_calib, "\n")

# Show improvement vs population average
new_dbh <- 50
h_pop   <- ht_from_dbh(new_dbh, "naslund", a = fit$coefficients["a"],
                        b = fit$coefficients["b"])
h_cal   <- ht_from_dbh(new_dbh, "naslund", a = res$coefficients["a"],
                        b = res$coefficients["b"])
h_true  <- ht_naslund(new_dbh, a = 1.9, b = 0.05)
cat(sprintf("\nDBH=50cm:  pop-avg=%.1fm  calibrated=%.1fm  true=%.1fm\n",
            h_pop, h_cal, h_true))
