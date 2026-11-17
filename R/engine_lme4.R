## The lme4 engine: maximum likelihood via lme4::lmer(). This file holds all
## code that touches merMod objects; everything else in the package works on
## the canonical tables it returns.

.engine_check.mpa_engine_lme4 <- function(engine, ...) {
  if (!requireNamespace("lme4", quietly = TRUE)) stop("Package `lme4` is required.")
  if (engine$df == "satterthwaite" && !requireNamespace("lmerTest", quietly = TRUE))
    stop("`df = \"satterthwaite\"` requires the `lmerTest` package.")
  invisible(TRUE)
}

.engine_packages.mpa_engine_lme4 <- function(engine, ...) c("lme4", if (engine$df == "satterthwaite") "lmerTest")

.engine_label.mpa_engine_lme4 <- function(engine, native, ...) paste0(if (engine$REML) "REML" else "ML", " fit (lme4)")

.engine_fit.mpa_engine_lme4 <- function(engine, formula, data, start = NULL, ...) {
  ## the control is kept as an unevaluated call so that update() on the fit works anywhere
  ctrl_call <- if (is.null(engine$control))
    bquote(lme4::lmerControl(optimizer = .(engine$optimizer), calc.derivs = .(engine$check_derivs)))
  else engine$control
  ## the formula's environment carries the data, so that update() on the merMod
  ## (which lme4 evaluates in environment(formula)) works from anywhere
  dat <- data
  fenv <- new.env(parent = globalenv()); fenv$dat <- dat
  environment(formula) <- fenv
  cl <- bquote(lme4::lmer(.(formula), data = dat, REML = .(engine$REML), control = .(ctrl_call)))
  if (!is.null(start)) cl$start <- start
  for (nm in names(engine$extra)) cl[[nm]] <- engine$extra[[nm]]
  eval(cl)
}

## Starting values for the variance parameters from a previous fit. lme4 stores
## theta blockwise, one block per random term (k(k+1)/2 parameters for k
## correlated effects; `||` terms are expanded into one-effect blocks). Blocks
## that exist identically in the previous fit (same grouping factor, same
## effects) take its estimates; new blocks keep lme4's defaults (1 on the
## diagonal, 0 off the diagonal). So the previous step's fit can serve as a
## start even when a step adds a random slope.
.engine_start.mpa_engine_lme4 <- function(engine, formula, data, old_fit, ...) {
  if (!inherits(old_fit$fit, "merMod")) {
    message("`start_from`: the previous fit was not estimated with lme4; starting values not used.")
    return(NULL)
  }
  probe <- lme4::lFormula(formula, data = data, REML = engine$REML)
  reTrms <- probe$reTrms
  new_cn <- reTrms$cnms; old_cn <- old_fit$fit@cnms
  theta_new <- reTrms$theta
  theta_old <- lme4::getME(old_fit$fit, "theta")
  seg <- function(cn) { k <- lengths(cn); split(seq_len(sum(k * (k + 1) / 2)), rep(seq_along(k), k * (k + 1) / 2)) }
  seg_new <- seg(new_cn); seg_old <- seg(old_cn)
  if (length(theta_new) != sum(lengths(seg_new))) return(NULL)      # unexpected layout: let lme4 start itself
  key <- function(cn) paste(names(cn), vapply(cn, paste, character(1), collapse = "+"))
  k_new <- key(new_cn); k_old <- key(old_cn)
  used <- logical(length(k_old)); matched <- 0L
  for (b in seq_along(k_new)) {
    j <- which(k_old == k_new[b] & !used)[1]
    if (!is.na(j)) { theta_new[seg_new[[b]]] <- theta_old[seg_old[[j]]]; used[j] <- TRUE; matched <- matched + 1L }
  }
  if (matched == 0L) {
    message("`start_from`: no random term of the previous fit matches; starting values not used.")
    return(NULL)
  }
  if (matched < length(k_new))
    message("`start_from`: random-effects structure differs; starting values taken over for ", matched,
            " of ", length(k_new), " random terms.")
  list(theta = unname(theta_new))
}

.engine_diagnostics.mpa_engine_lme4 <- function(engine, native, ...) {
  msgs <- native@optinfo$conv$lme4$messages
  msgs <- msgs[!grepl("singular", msgs)]   # singularity is reported separately
  if (is.null(msgs)) msgs <- character(0)
  list(converged = length(msgs) == 0, singular = lme4::isSingular(native), messages = msgs,
       details = list(optimizer = engine$optimizer, check_derivs = engine$check_derivs))
}

.engine_fit_stats.mpa_engine_lme4 <- function(engine, native, ...) {
  ll <- stats::logLik(native)
  list(logLik = as.numeric(ll), df = attr(ll, "df"), AIC = stats::AIC(native), BIC = stats::BIC(native))
}

.engine_tables.mpa_engine_lme4 <- function(engine, native, level, ...) {
  list(fixed = .lme4_fixed_table(native, engine$df, level),
       variance = .lme4_variance_table(native),
       correlations = .lme4_correlation_table(native),
       n_units = vapply(native@flist, function(g) length(levels(g)), integer(1)),
       labels = list(estimate = if (engine$REML) "REML estimates" else "ML estimates",
                     test = if (engine$df == "z") "Wald z-tests" else "Satterthwaite df",
                     interval = paste0(round(100 * level), "% Wald CI")))
}

.lme4_fixed_table <- function(fit, df = "z", level = 0.95) {
  if (df == "satterthwaite") {
    cf <- as.data.frame(stats::coef(summary(lmerTest::as_lmerModLmerTest(fit))))
    names(cf) <- c("estimate", "se", "df", "t", "p")
  } else {
    cf <- as.data.frame(stats::coef(summary(fit)))
    names(cf) <- c("estimate", "se", "z")
    cf$p <- 2 * stats::pnorm(-abs(cf$z))
  }
  q <- stats::qnorm(1 - (1 - level) / 2)
  cf$ci_lower <- cf$estimate - q * cf$se
  cf$ci_upper <- cf$estimate + q * cf$se
  cf <- cbind(term = rownames(cf), cf)
  rownames(cf) <- NULL
  cf
}

## lme4 names the groups of "||" terms perceiver, perceiver.1, ... and turns
## "perceiver:target" into "perceiver.target"; map them back to the grouping factors
.lme4_unit_names <- function(grp, fit) {
  known <- names(fit@flist)
  g <- gsub("\\.\\d+$", "", grp)
  fix <- !(g %in% c(known, "Residual"))
  g[fix] <- vapply(g[fix], function(z) {
    cand <- known[gsub(":", ".", known, fixed = TRUE) == z]
    if (length(cand) == 1) cand else z
  }, character(1))
  g
}

.lme4_variance_table <- function(fit) {
  vc <- as.data.frame(lme4::VarCorr(fit))
  vc <- vc[is.na(vc$var2), c("grp", "var1", "vcov", "sdcor")]
  vc$var1[is.na(vc$var1)] <- ifelse(vc$grp[is.na(vc$var1)] == "Residual", "", "(Intercept)")
  names(vc) <- c("unit", "effect", "variance", "sd")
  vc$unit <- .lme4_unit_names(vc$unit, fit)
  fe <- names(lme4::fixef(fit))
  vc <- vc[order(match(vc$unit, c(names(fit@flist), "Residual")), vc$effect != "(Intercept)", match(vc$effect, fe)), ]
  rownames(vc) <- NULL
  vc
}

.lme4_correlation_table <- function(fit) {
  vc <- as.data.frame(lme4::VarCorr(fit))
  vc <- vc[!is.na(vc$var2), c("grp", "var1", "var2", "sdcor")]
  if (!nrow(vc)) return(NULL)
  names(vc) <- c("unit", "effect1", "effect2", "correlation")
  vc$unit <- .lme4_unit_names(vc$unit, fit)
  rownames(vc) <- NULL
  vc
}

## L %*% beta with Wald intervals; L has the fixed effects as columns (by name)
.engine_lincomb.mpa_engine_lme4 <- function(engine, native, L, level, ...) {
  beta <- lme4::fixef(native)
  L <- L[, names(beta), drop = FALSE]
  V <- as.matrix(stats::vcov(native))
  est <- drop(L %*% beta)
  se <- sqrt(rowSums((L %*% V) * L))
  q <- stats::qnorm(1 - (1 - level) / 2)
  data.frame(estimate = est, se = se, lower = est - q * se, upper = est + q * se)
}

.engine_print_diagnostics.mpa_engine_lme4 <- function(engine, x, ...) {
  d <- x$diagnostics
  other <- if (identical(engine$optimizer, "bobyqa")) "nloptwrap" else "bobyqa"
  if (!d$converged) cat("WARNING: convergence problems: ", paste(d$messages, collapse = "; "),
                        "\n         (a gradient only slightly above lme4's tolerance is usually harmless; ",
                        "refit with the other optimizer, optimizer = \"", other, "\", to check)\n", sep = "")
  else if (!d$singular) cat("Converged without warnings.\n")
  if (isFALSE(engine$check_derivs)) cat("NOTE: fitted with check_derivs = FALSE - lme4's gradient-based convergence checks were skipped.\n")
  if (d$singular) cat("WARNING: singular fit (some variance components estimated as zero or correlations as +/-1).\n",
                      "         mpa_prune(<this fit>) removes the zero-variance terms; correlated = FALSE drops the correlations.\n", sep = "")
  s <- x$fit_stats
  cat("logLik = ", round(s$logLik, 1), ", AIC = ", round(s$AIC, 1), ", BIC = ", round(s$BIC, 1), "\n", sep = "")
  invisible(NULL)
}
