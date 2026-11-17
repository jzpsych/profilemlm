## Tests of the engine layer: engine objects, the canonical fit object, and the
## agreement of the lme4 engine with profilemlm 0.3.0 (before the engine layer
## existed). The reference values in fixtures/reference_toy_0.3.0.rds were
## computed with 0.3.0 on the toy data of make_toy_mpa() (test-mpa_formula.R).

ref <- readRDS(test_path("fixtures", "reference_toy_0.3.0.rds"))
## The reference values were computed with lme4 1.1-35. The toy models are
## singular, and lme4 >= 2.0 ends at a slightly different point of the flat
## likelihood (with a logLik that is equal or higher in the third decimal). Under
## lme4 1.x the values are compared exactly; under lme4 >= 2.0 the fixed effects
## and the likelihood are compared with a tolerance and the boundary estimates
## (variances, correlations, singularity) are not compared.
strict <- utils::packageVersion("lme4") < "2.0"
expect_same_fit <- function(m, r, tol = 1e-10) {
  if (strict) {
    expect_equal(as.data.frame(m$fixed), r$fixed, tolerance = tol, ignore_attr = TRUE)
    expect_equal(as.data.frame(m$variance), r$variance, tolerance = tol, ignore_attr = TRUE)
    if (is.null(r$correlations)) expect_null(m$correlations)
    else expect_equal(as.data.frame(m$correlations), r$correlations, tolerance = tol, ignore_attr = TRUE)
    expect_equal(m$logLik, r$logLik, tolerance = tol); expect_equal(m$AIC, r$AIC, tolerance = tol)
    expect_equal(m$singular, r$singular); expect_equal(m$converged, r$converged)
  } else {
    expect_equal(m$fixed$estimate, r$fixed$estimate, tolerance = 5e-3)
    expect_equal(m$fixed$se, r$fixed$se, tolerance = 5e-2)
    expect_equal(m$logLik, r$logLik, tolerance = 1e-5)
    expect_gte(m$logLik, r$logLik - 1e-3)
  }
  expect_equal(m$n_obs, r$n_obs); expect_equal(m$n_units, r$n_units); expect_equal(mpa_rows(m), r$rows)
}
loose <- function(strict_tol, loose_tol) if (strict) strict_tol else loose_tol

test_that("engine objects: constructors, resolution from names and further arguments", {
  e <- engine_lme4()
  expect_s3_class(e, c("mpa_engine_lme4", "mpa_engine"))
  expect_equal(e$optimizer, "bobyqa"); expect_false(e$REML); expect_equal(e$df, "z")
  expect_output(print(e), "lme4")
  e2 <- profilemlm:::.mpa_resolve_engine("lme4", list(optimizer = "nloptwrap", REML = TRUE))
  expect_equal(e2$optimizer, "nloptwrap"); expect_true(e2$REML)
  e3 <- profilemlm:::.mpa_resolve_engine(engine_lme4(check_derivs = FALSE), list(REML = TRUE))   # dots override
  expect_true(e3$REML); expect_false(e3$check_derivs)
  expect_error(profilemlm:::.mpa_resolve_engine("glmmTMB"), "lme4")
  expect_error(engine_lme4(df = "kenward"))
  expect_s3_class(engine_brms(chains = 2), c("mpa_engine_brms", "mpa_engine"))
  expect_equal(engine_brms(chains = 2)$chains, 2)
})

test_that("the lme4 engine reproduces profilemlm 0.3.0 exactly (toy data)", {
  skip_if_not_installed("lme4")
  d <- make_toy_mpa()
  f0 <- mpa_formula(d, criterion = "self", random_units = c("perceiver", "target"))
  f1 <- mpa_formula(d, criterion = "self", predictors = "des.n", interactions = "criterion", random_units = c("perceiver", "target"))
  m0 <- mpa_fit(f0, d); m1 <- mpa_fit(f1, d)
  expect_same_fit(m0, ref$toy0); expect_same_fit(m1, ref$toy1)
  expect_same_fit(suppressWarnings(mpa_fit(f0, d, REML = TRUE)), ref$toy_reml)
  expect_same_fit(suppressWarnings(mpa_fit(f0, d, engine = engine_lme4(REML = TRUE))), ref$toy_reml)
  expect_same_fit(suppressWarnings(mpa_fit(f0, d, optimizer = "nloptwrap")), ref$toy_nlopt)
  cmp <- mpa_compare(m0, m1)
  rt <- suppressWarnings(mpa_random_test(mpa_formula(d, criterion = "self", random_units = "perceiver"), d))
  if (strict) {
    expect_equal(cmp$variance, ref$toy_cmp01$variance, tolerance = 1e-10)
    expect_equal(cmp$lrt, ref$toy_cmp01$lrt, tolerance = 1e-10); expect_equal(cmp$delta, ref$toy_cmp01$delta, tolerance = 1e-10)
    expect_equal(as.data.frame(rt), ref$toy_rt, tolerance = 1e-8, ignore_attr = TRUE)
  } else {
    expect_equal(cmp$lrt$chisq, ref$toy_cmp01$lrt$chisq, tolerance = 5e-2)
    expect_equal(dim(as.data.frame(rt)), dim(ref$toy_rt))
  }
  skip_if_not_installed("lmerTest")
  expect_same_fit(suppressWarnings(mpa_fit(f0, d, df = "satterthwaite")), ref$toy_satt, tol = 1e-6)   # lme4 2.x warns about convergence on the toy data
})

test_that("the canonical fit object: engine, settings, rows, labels, convenience fields", {
  skip_if_not_installed("lme4")
  d <- make_toy_mpa()
  m <- mpa_fit(mpa_formula(d, criterion = "self", random_units = "perceiver"), d, level = 0.9)
  expect_s3_class(m, c("mpa_fit_lme4", "mpa_fit"))
  expect_equal(m$engine, "lme4"); expect_s3_class(m$settings, "mpa_engine_lme4")
  expect_equal(m$rows, seq_len(nrow(as.data.frame(d))))
  expect_equal(m$level, 0.9)
  expect_equal(m$fixed$ci_upper - m$fixed$estimate, qnorm(0.95) * m$fixed$se)
  expect_equal(m$labels$interval, "90% Wald CI")
  expect_true(all(c("converged", "singular", "messages", "details") %in% names(m$diagnostics)))
  expect_equal(m$fit_stats$logLik, m$logLik); expect_true(is.numeric(m$fit_stats$df))
  expect_false(m$REML); expect_equal(m$df, "z")
  expect_output(print(m), "ML fit \\(lme4\\)")
  expect_output(print(m), "ML estimates, Wald z-tests, 90% Wald CI")
  expect_equal(mpa_fit(m$formula, d, engine = m)$logLik, m$logLik)       # "same engine as that fit"
})

test_that("missing values: rows are chosen by the package, refits keep the engine settings", {
  skip_if_not_installed("lme4")
  d <- make_toy_mpa()
  d2 <- as.data.frame(d); d2$know.n[d2$target %in% c("t1", "t2")] <- NA
  d2 <- structure(d2, mpa = attr(d, "mpa"), class = class(d))
  g0 <- mpa_formula(d2, criterion = "self", random_units = "perceiver")
  g1 <- mpa_formula(d2, criterion = "self", predictors = "know.n", interactions = "criterion", random_units = "perceiver")
  m0 <- mpa_fit(g0, d2); m1 <- suppressWarnings(mpa_fit(g1, d2))
  expect_same_fit(m0, ref$na_m0); expect_same_fit(m1, ref$na_m1)
  expect_equal(nrow(m1$data), nrow(d2)); expect_lt(length(m1$rows), nrow(d2))   # data kept whole, rows recorded
  expect_same_fit(mpa_fit(g0, d2, same_rows_as = m1), ref$na_m0b)
  expect_same_fit(mpa_fit(g1, d2, start_from = m1, check_derivs = FALSE), ref$na_m1s)
  cmp <- suppressMessages(mpa_compare(m0, m1))
  expect_true(cmp$refitted)
  if (strict) {
    expect_equal(cmp$variance, ref$na_cmp$variance, tolerance = 1e-8)
    expect_equal(cmp$lrt, ref$na_cmp$lrt, tolerance = 1e-8)
  } else expect_equal(cmp$lrt$chisq, ref$na_cmp$lrt$chisq, tolerance = 5e-2)
  m0n <- mpa_fit(g0, d2, optimizer = "nloptwrap", check_derivs = FALSE)
  cmpn <- suppressMessages(mpa_compare(m0n, m1))
  expect_equal(cmpn$fit0$settings$optimizer, "nloptwrap"); expect_false(cmpn$fit0$settings$check_derivs)   # refit keeps the settings
})

test_that("mpa_plot computes slopes from the fixed effects (agrees with emmeans' exact derivative)", {
  skip_if_not_installed("lme4"); skip_if_not_installed("ggplot2")
  d <- make_toy_mpa()
  d2 <- as.data.frame(d); d2$know.n[d2$target %in% c("t1", "t2")] <- NA
  d2 <- structure(d2, mpa = attr(d, "mpa"), class = class(d))
  g1 <- mpa_formula(d2, criterion = "self", predictors = "know.n", interactions = "criterion", random_units = "perceiver")
  m1 <- suppressWarnings(mpa_fit(g1, d2))
  pd <- function(p) as.data.frame(attr(p, "data"))[, intersect(c(".x", ".by", ".y", ".lcl", ".ucl"), names(attr(p, "data")))]
  expect_equal(pd(mpa_plot(m1, slope = "self.dist", x = "know.n")), ref$na_plot, tolerance = loose(1e-7, 1e-2), ignore_attr = TRUE)
  expect_equal(lapply(mpa_plots(m1), pd), ref$na_plots, tolerance = loose(1e-7, 1e-2), ignore_attr = TRUE)
  ## a slope with a quadratic term: the derivative at the mean of the variable
  f <- mpa_formula(d, criterion = "self", predictors = c("des.n", "know.n"), poly = "des.n", interactions = "all",
                   random_units = "perceiver", correlated = FALSE)
  m <- suppressWarnings(mpa_fit(f, d))
  p <- mpa_plot(m, slope = "des.n", x = "know.n", n = 5)
  skip_if_not_installed("emmeans")
  e <- suppressMessages(emmeans::emtrends(m$fit, specs = "know.n", var = "des.n", at = list(know.n = attr(p, "data")$.x),
                                          data = m$data, lmer.df = "asymp", delta.var = 1e-7))
  s <- as.data.frame(summary(e, infer = c(TRUE, FALSE)))
  expect_equal(attr(p, "data")$.y, s$des.n.trend, tolerance = 1e-6)
  expect_equal(attr(p, "data")$.lcl, s$asymp.LCL, tolerance = 1e-6)
})

test_that("mpa_random_test works with parallel workers and the engine layer", {
  skip_if_not_installed("lme4")
  skip_on_cran()
  d <- make_toy_mpa()
  rt <- suppressWarnings(mpa_random_test(mpa_formula(d, criterion = "self", random_units = "perceiver"), d, cores = 2))
  expect_equal(as.data.frame(rt), ref$toy_rt, tolerance = 1e-8, ignore_attr = TRUE)
  expect_error(mpa_random_test(mpa_formula(d, criterion = "self", random_units = "perceiver"), d, engine = "brms"), "lme4 engine only")
})

test_that("the default engine can be set by an option", {
  skip_if_not_installed("lme4")
  d <- make_toy_mpa()
  old <- options(profilemlm.engine = engine_lme4(REML = TRUE)); on.exit(options(old))
  m <- mpa_fit(mpa_formula(d, criterion = "self", random_units = "perceiver"), d)
  expect_true(m$REML); expect_equal(m$labels$estimate, "REML estimates")
})
