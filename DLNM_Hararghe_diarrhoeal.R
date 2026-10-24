# ==============================================================================
#  Temperature, rainfall and diarrhoeal mortality (all ages and under five)
#  Hararghe HDSS verbal autopsy data, eastern Ethiopia, weekly 2015-2023
#  Distributed lag non-linear models (DLNM), quasi-Poisson time-series regression
#
#  Outcomes:   deaths_total (all ages), deaths_under5 (under five)
#  Subgroups:  under-five male and female
#  Exposures:  one temperature variable (Tmin, Tavg or Tmax; chosen by Q-AIC)
#              + weekly total rainfall
#
# ==============================================================================

# ---- 0. Setup ----------------------------------------------------------------
suppressPackageStartupMessages({
  library(dlnm)
  library(splines)
})

setwd("D:/APHRC/Hararghe HDSS-Ethiopia")
data_file <- "DiarrClimate_Harar_Weekly.csv"
out_dir   <- "Ethiopia_Hararghe"
dir.create(out_dir, showWarnings = FALSE)

# ---- Loading Data -----------------------------------------------------------------
dat <- read.csv(data_file, stringsAsFactors = FALSE)
dat$Week <- as.Date(dat$Week)
dat$time <- seq_len(nrow(dat))

# ---- Model settings -------------------------------------------------------
maxlag    <- 3                         # lag 0-3 weeks
lag_knots <- 2                         # one internal knot on the lag scale
nyears    <- as.numeric(diff(range(dat$Week))) / 365.25
time_df   <- round(4 * nyears)         # 4 df per year = 36

temp_labels <- c(Tmin = "Minimum temperature", Tavg = "Mean temperature",
                 Tmax = "Maximum temperature")
temp_axis   <- c(Tmin = "Weekly mean minimum temperature (\u00b0C)",
                 Tavg = "Weekly mean temperature (\u00b0C)",
                 Tmax = "Weekly mean maximum temperature (\u00b0C)")
outcome_labels <- c(deaths_total = "All ages", deaths_under5 = "Under five",
                    under5_male = "Under-five male", under5_female = "Under-five female")

# Cross-basis function
make_cb <- function(x, term, lagknots = lag_knots) {
  cb <- crossbasis(x, lag = maxlag,
                   argvar = list(fun = "ns", knots = equalknots(x, nk = 2)),
                   arglag = list(fun = "ns", knots = lagknots))
  attr(cb, "term") <- term
  cb
}


cpx <- function(cb, model, ...) {
  b <- coef(model)
  i <- grep(paste0("^", attr(cb, "term"), "v[0-9]+\\.l[0-9]+$"), names(b))
  stopifnot(length(i) == ncol(cb), !anyNA(b[i]))
  crosspred(cb, coef = b[i], vcov = vcov(model)[i, i, drop = FALSE],
            model.link = "log", ...)
}

# Main model: outcome ~ cb(temperature) + cb(rainfall) + ns(time)
fit_model <- function(outcome, cbT, cbR, data = dat, tdf = time_df) {
  f <- as.formula(paste(outcome, "~ cbT + cbR + ns(time, df = tdf)"))
  glm(f, family = quasipoisson(), data = data)
}

#--------------------------------------------------------------------------
# models selection using Q-AIC and  Q-BIC 
#------------------------------------------------------------------------------
qaic <- function(m) {
  ll  <- sum(dpois(m$y, m$fitted.values, log = TRUE))
  phi <- summary(m)$dispersion
  k   <- summary(m)$df[3]
  c(QAIC = -2 * ll + 2 * phi * k,
    QBIC = -2 * ll + log(length(m$y)) * phi * k)
}

fmt_rr <- function(f, l, u) sprintf("%.2f (%.2f\u2013%.2f)", f, l, u)

cb_rain <- make_cb(dat$Rain, "cbR")
cb_temp_all <- lapply(c(Tmin = "Tmin", Tavg = "Tavg", Tmax = "Tmax"),
                      function(v) make_cb(dat[[v]], "cbT"))


sel <- do.call(rbind, lapply(c("deaths_total", "deaths_under5"), function(y) {
  do.call(rbind, lapply(names(cb_temp_all), function(v) {
    m <- fit_model(y, cb_temp_all[[v]], cb_rain)
    data.frame(Outcome = outcome_labels[y],
               Model = paste(v, "+ Rainfall"), Temp = v,
               QAIC = round(qaic(m)["QAIC"], 1), QBIC = round(qaic(m)["QBIC"], 1))
        
  }))
}))
sel$Selected <- ave(sel$QAIC, sel$Outcome, FUN = function(z) z == min(z))
sel$Selected <- ifelse(sel$Selected == 1, "Yes", "")
print(sel)
write.csv(sel[, names(sel) != "Temp"], file.path(out_dir, "TableS1_model_selection.csv"), row.names = FALSE)

best_temp <- setNames(sel$Temp[sel$Selected == "Yes"], c("deaths_total", "deaths_under5"))
cat("Selected temperature variable:\n"); print(best_temp)


#--------------------------------------------------------------------
# All ages and children under 5
#----------------------------------------------------------------------
find_mmt <- function(cb, model, x, lo = .01, hi = .99) {
  grid <- seq(quantile(x, lo), quantile(x, hi), by = 0.01)
  cp <- cpx(cb, model, at = grid, cen = median(x))
  grid[which.min(cp$allRRfit)]
}

pred_grid <- function(cb, model, x, cen, by)
  cpx(cb, model, from = min(x), to = max(x), by = by, cen = cen, bylag = 1)

rr_table <- function(cb, model, x, cen, label, unit, probs) {
  at <- as.numeric(quantile(x, probs))
  cp <- cpx(cb, model, at = at, cen = cen, bylag = 1)
  lag_txt <- sapply(seq_len(ncol(cp$matRRfit)), function(j)
    fmt_rr(cp$matRRfit[, j], cp$matRRlow[, j], cp$matRRhigh[, j]))
  if (is.null(dim(lag_txt))) lag_txt <- matrix(lag_txt, nrow = 1)
  colnames(lag_txt) <- paste0("Lag ", 0:maxlag)
  data.frame(Exposure = label, Percentile = paste0("P", probs * 100),
             Value = sprintf(paste0("%.1f", unit), at), lag_txt,
             check.names = FALSE, row.names = NULL)
}

# Percentiles at which RRs are reported. Different sets for the two exposures because
# rainfall is strongly right-skewed: its lower percentiles (P10 = 0.02 mm, P25 = 0.65 mm)
# are effectively equal to the 0 mm reference and would give uninformative RRs of ~1.
probs_t <- c(.10, .25, .75, .95, .99)        # temperature (reference: MMT)
probs_r <- c(.50, .75, .90, .95, .99)        # rainfall (reference: 0 mm)
cen_rain <- 0


run_outcome <- function(y) {
  tv  <- best_temp[[y]]
  x   <- dat[[tv]]
  cbT <- cb_temp_all[[tv]]
  
  # Fit the primary model
  mod <- fit_model(y, cbT, cb_rain)
  
  # Estimate the minimum mortality temperature (MMT)
  mmt <- find_mmt(cbT, mod, x)
  
  # Primary lag-specific relative risks
  tab <- rbind(
    rr_table(
      cbT, mod, x, mmt,
      temp_labels[[tv]], " °C", probs_t
    ),
    rr_table(
      cb_rain, mod, dat$Rain, cen_rain,
      "Rainfall", " mm", probs_r
    )
  )
  
  # Residual diagnostics
  res <- residuals(mod, type = "deviance")
  
  full_res <- rep(NA, nrow(dat))
  full_res[as.numeric(names(res))] <- res
  
  pa <- pacf(
    full_res,
    lag.max = 8,
    plot = FALSE,
    na.action = na.pass
  )
  
  # Exposure values at the 95th percentile
  p95_t <- as.numeric(quantile(x, 0.95, na.rm = TRUE))
  p95_r <- as.numeric(quantile(dat$Rain, 0.95, na.rm = TRUE))
  
  # Return primary analysis results
  list(
    outcome = y,
    label = outcome_labels[[y]],
    tv = tv,
    x = x,
    cbT = cbT,
    model = mod,
    mmt = mmt,
    mmt_pct = mean(x <= mmt, na.rm = TRUE),
    
    cp_T = pred_grid(cbT, mod, x, mmt, 0.05),
    cp_R = pred_grid(cb_rain, mod, dat$Rain, cen_rain, 0.25),
    
    lr_T = cpx(cbT, mod, at = p95_t, cen = mmt, bylag = 0.1),
    lr_R = cpx(
      cb_rain, mod,
      at = p95_r,
      cen = cen_rain,
      bylag = 0.1
    ),
    
    p95_t = p95_t,
    p95_r = p95_r,
    tab = tab,
    n = length(mod$y),
    dispersion = summary(mod)$dispersion,
    pacf_max = max(abs(pa$acf), na.rm = TRUE),
    pacf_bound = 1.96 / sqrt(length(res)),
    knots_T = attr(cbT, "argvar")$knots
  )
}

# Run primary analyses
res_list <- list()

for (y in c("deaths_total", "deaths_under5")) {
  res_list[[y]] <- run_outcome(y)
  
  cat(
    "\n====", outcome_labels[[y]], "====\n",
    "MMT:", round(res_list[[y]]$mmt, 2),
    "(percentile",
    round(100 * res_list[[y]]$mmt_pct), ")\n"
  )
  
  print(res_list[[y]]$tab, right = FALSE)
  cat("Modelled weeks:", res_list[[y]]$n, "\n")
}

A <- res_list$deaths_total
U <- res_list$deaths_under5

# Export primary lag-specific RR results
tab_rr <- rbind(
  cbind(Outcome = A$label, A$tab),
  cbind(Outcome = U$label, U$tab)
)

write.csv(
  tab_rr,
  file.path(out_dir, "Table 3 and 4_RR_lag_specific.csv"),
  row.names = FALSE
)

#------------------------------------------------------------------------
# Under-five sex-stratified analysis
#---------------------------------------------------------------------
tvU <- U$tv; xU <- U$x; cbTU <- U$cbT
mod_m <- fit_model("under5_male",   cbTU, cb_rain)
mod_f <- fit_model("under5_female", cbTU, cb_rain)

sex_tab <- do.call(rbind, lapply(list(Male = mod_m, Female = mod_f), function(m)
  cbind(Sex = NA, rbind(
    rr_table(cbTU, m, xU, U$mmt, temp_labels[[tvU]], " \u00b0C", probs_t),
    rr_table(cb_rain, m, dat$Rain, cen_rain, "Rainfall", " mm", probs_r)))))
sex_tab$Sex <- rep(c("Male", "Female"), each = nrow(sex_tab) / 2); rownames(sex_tab) <- NULL

print(sex_tab, right = FALSE)
write.csv(sex_tab, file.path(out_dir, "Table5_under5_by_sex.csv"), row.names = FALSE)

cp_T_m <- pred_grid(cbTU, mod_m, xU, U$mmt, 0.05)
cp_T_f <- pred_grid(cbTU, mod_f, xU, U$mmt, 0.05)
cp_R_m <- pred_grid(cb_rain, mod_m, dat$Rain, cen_rain, 0.25)
cp_R_f <- pred_grid(cb_rain, mod_f, dat$Rain, cen_rain, 0.25)

#-----------------------------------------------------------------------
# FIGURES
#-------------------------------------------------------------------------
col_t <- "#1b7837"; fill_t <- adjustcolor(col_t, 0.22)
col_r <- "#b35806"; fill_r <- adjustcolor(col_r, 0.22)
lty_m <- 1; lty_f <- 5   # male solid, female dashed

save_both <- function(name, w, h, expr) {
  tiff(file.path(out_dir, paste0(name, ".tiff")), width = w, height = h, units = "in",
       res = 600, compression = "lzw")
  eval.parent(substitute(expr)); dev.off()
}
band <- function(x, lo, hi, fill) polygon(c(x, rev(x)), c(lo, rev(hi)), col = fill, border = NA)

ycap <- 6   # upper limit of plotted y axis; CIs above it are cut
lag_plot <- function(cp, l, col, fill, xlab, main, ylim = c(0, 3)) {
  j <- l + 1
  plot(cp$predvar, cp$matRRfit[, j], type = "n", ylim = ylim, xlab = xlab,
       ylab = "Relative risk", main = main, font.main = 1)
  band(cp$predvar, cp$matRRlow[, j], pmin(cp$matRRhigh[, j], ylim[2] * 2), fill)
  lines(cp$predvar, cp$matRRfit[, j], col = col, lwd = 2)
  abline(h = 1, lty = 2, col = "grey40")
}
tax <- function(r) temp_axis[[r$tv]]
tshort <- function(r) { z <- sub("Weekly mean ", "", sub(" \\(.*", "", tax(r))); paste0(toupper(substr(z, 1, 1)), substring(z, 2)) }


# Figures 2 and 3: lag-specific curves 
lag_fig <- function(r, name) save_both(name, 14, 7, {
  par(mfrow = c(2, 4), mar = c(4.5, 4.5, 3, 1), las = 1)
  for (l in 0:3) lag_plot(r$cp_T, l, col_t, fill_t, tax(r), paste0(tshort(r), ": lag ", l))
  for (l in 0:3) lag_plot(r$cp_R, l, col_r, fill_r, "Weekly total rainfall (mm)", paste0("Rainfall: lag ", l))
})
lag_fig(A, "Figure2_lag_specific_all_ages")
lag_fig(U, "Figure3_lag_specific_under5")


#----------------------------------------------------------------------------------
# Figure 4 under-five lag-specific curves by sex 
#--------------------------------------------------------------------------------
save_both("Figure4_under5_sex_lag_specific", 14, 12, {
  par(mfrow = c(4, 4), mar = c(4.5, 4.5, 3, 1), las = 1)
  for (l in 0:3) lag_plot(cp_T_m, l, col_t, fill_t, tax(U), paste0("Male, ", tolower(tshort(U)), ": lag ", l), c(0, 5))
  for (l in 0:3) lag_plot(cp_R_m, l, col_r, fill_r, "Rainfall (mm)", paste0("Male, rainfall: lag ", l), c(0, 5))
  for (l in 0:3) lag_plot(cp_T_f, l, col_t, fill_t, tax(U), paste0("Female, ", tolower(tshort(U)), ": lag ", l), c(0, 5))
  for (l in 0:3) lag_plot(cp_R_f, l, col_r, fill_r, "Rainfall (mm)", paste0("Female, rainfall: lag ", l), c(0, 5))
})


