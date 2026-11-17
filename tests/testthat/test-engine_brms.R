## The brms engine. Needs brms and a Stan backend (cmdstanr with CmdStan, or
## rstan with a C++ toolchain); the fit is small (one chain) but compiling the
## Stan model takes about a minute, so this file is skipped on CRAN.

has_backend <- function() {
  (requireNamespace("cmdstanr", quietly = TRUE) &&
     !is.null(tryCatch(cmdstanr::cmdstan_version(error_on_NA = FALSE), error = function(e) NULL))) ||
    requireNamespace("rstan", quietly = TRUE)
}

test_that("brms coefficient renaming is reproduced", {
  skip_if_not_installed("brms")
  x <- c("(Intercept)", "I(des^2)", "self.dist:I(DES.n^2):knowing.n", "a-b", "x:z", "log(x)")
  expect_equal(profilemlm:::.brms_rename(x), brms:::rename(x))
})

test_that("the brms engine fits the tutorial-type model and feeds the same downstream functions", {
  skip_on_cran()
  skip_if_not_installed("brms"); skip_if_not_installed("posterior")
  skip_if_not(has_backend(), "no Stan backend")
  d <- mpa_simulate(n_perceiver = 16, n_target = 16, n_item = 6, group_size = 4, seed = 3, rating_scale = c(1, 101),
                    fixed = c(intercept = 51, norm = 0.7, dist = 0.25, des = 5, know_dist = 0.12))
  f <- mpa_formula(d, criterion = "self", random_slopes = "des", predictors = "know", poly = "des", interactions = "criterion",
                   random_units = c("perceiver", "target", "item"), correlated = FALSE)
  m_l <- suppressWarnings(mpa_fit(f, d))
  m_b <- suppressWarnings(mpa_fit(f, d, engine = "brms", chains = 1, iter = 600, seed = 1))
  expect_s3_class(m_b, c("mpa_fit_brms", "mpa_fit"))
  expect_equal(m_b$engine, "brms"); expect_s3_class(m_b$fit, "brmsfit")
  expect_equal(m_b$fixed$term, m_l$fixed$term)                       # canonical names
  expect_equal(m_b$variance[, c("unit", "effect")], m_l$variance[, c("unit", "effect")], ignore_attr = TRUE)
  expect_equal(m_b$n_units, m_l$n_units); expect_equal(m_b$n_obs, m_l$n_obs); expect_equal(mpa_rows(m_b), mpa_rows(m_l))
  expect_true(all(c("rhat", "ess_bulk", "ess_tail") %in% names(m_b$fixed)))
  expect_true(all(c("sd_lower", "sd_upper") %in% names(m_b$variance)))
  expect_equal(m_b$variance$variance, m_b$variance$sd^2)
  expect_false(m_b$singular); expect_true(is.na(m_b$logLik))
  expect_null(m_b$correlations)
  ## the two engines agree on the fixed effects (within sampling error of a short chain)
  z <- abs(m_b$fixed$estimate - m_l$fixed$estimate) / m_l$fixed$se
  expect_lt(max(z), 1.5)
  expect_output(print(m_b), "Bayesian fit \\(brms")
  expect_output(print(m_b), "posterior medians")
  expect_message(mpa_prune(m_b), "returned unchanged")
  rec <- mpa_recovery(m_b); expect_true(all(c("true", "estimate", "covered") %in% names(rec$fixed)))
  skip_if_not_installed("ggplot2")
  p <- mpa_plot(m_b, slope = "self.dist", x = "know", n = 5)
  pd <- attr(p, "data"); expect_equal(nrow(pd), 5); expect_true(all(pd$.lcl < pd$.y & pd$.y < pd$.ucl))
  pl <- mpa_plot(m_l, slope = "self.dist", x = "know", n = 5)
  expect_lt(max(abs(pd$.y - attr(pl, "data")$.y)), 0.1)
  ## comparison with a second brms fit; errors across engines
  f1 <- mpa_formula(d, criterion = "self", random_units = c("perceiver", "target", "item"), correlated = FALSE)
  m_b1 <- suppressWarnings(mpa_fit(f1, d, engine = engine_brms(chains = 1, iter = 600, seed = 1)))
  cmp <- mpa_compare(m_b1, m_b)
  expect_null(cmp$lrt); expect_true(is.na(cmp$delta[["AIC"]])); expect_equal(cmp$engine, "brms")
  expect_true("prop_reduced" %in% names(cmp$variance))
  expect_output(print(cmp), "no likelihood-ratio test")
  expect_error(mpa_compare(m_l, m_b), "different engines")
  expect_error(mpa_random_test(f1, d, engine = "brms"), "lme4 engine only")
  expect_message(suppressWarnings(mpa_fit(f1, d, engine = "brms", start_from = m_b1, chains = 1, iter = 300, seed = 1)), "start_from")
})
