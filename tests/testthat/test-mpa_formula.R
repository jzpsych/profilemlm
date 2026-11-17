test_that("base formula (step 1) has criterion components and four random units", {
  d <- make_toy_mpa()
  f <- mpa_formula(d, criterion = "self")
  s <- mpa_formula_string(f)
  expect_match(s, "^y ~ self.norm.c \\+ self.dist")
  expect_match(s, "\\(1 \\+ self.norm.c \\+ self.dist \\| perceiver\\)")
  expect_match(s, "\\| target\\)")
  expect_match(s, "\\| perceiver:target\\)")
  expect_match(s, "\\(1 \\+ self.dist \\| item\\)")
  expect_equal(f$dropped$item, "self.norm.c")
  expect_true(inherits(f$formula, "formula"))
  # the classical SAM: no item factor
  s3 <- mpa_formula_string(mpa_formula(d, criterion = "self", random_units = c("perceiver", "target", "dyad")))
  expect_false(grepl("item\\)", s3))
  # a group intercept (the toy data have no groups)
  expect_error(mpa_formula(d, criterion = "self", random_units = c("perceiver", "target", "dyad", "item", "group")), "no `group` column")
  dg <- mpa_simulate(n_perceiver = 8, n_item = 6, group_size = 4, seed = 3, rating_scale = c(1, 101))
  fg <- mpa_formula(dg, criterion = "self", random_units = c("perceiver", "target", "dyad", "item", "group"))
  expect_equal(fg$random$group, "1")
  expect_match(mpa_formula_string(fg), "\\(1 \\| group\\)")
  expect_false("group" %in% names(fg$dropped))
  # an intercept-only item factor via mpa_update()
  fi <- mpa_update(f, drop = list(item = "self.dist"))
  expect_equal(fi$random$item, "1")
  expect_match(mpa_formula_string(fi), "\\(1 \\| item\\)")
  # the old arguments are gone
  expect_error(mpa_formula(d, criterion = "self", random_item = "intercept"), "unused argument")
})

test_that("step 0 is intercept-only", {
  d <- make_toy_mpa()
  f <- mpa_formula(d, criterion = NULL)
  expect_match(mpa_formula_string(f), "^y ~ 1 \\+ \\(1 \\| perceiver\\) \\+ \\(1 \\| target\\) \\+ \\(1 \\| perceiver:target\\) \\+ \\(1 \\| item\\)$")
})

test_that("interactions, quadratic terms and identifiability filtering", {
  d <- make_toy_mpa()
  f <- mpa_formula(d, criterion = "self", random_slopes = c("des.n", "know.n"), poly = "des.n",
                   interactions = "all")
  # quadratic terms are fixed effects only, at every unit
  expect_false(any(grepl("^I\\(", unlist(f$random))))
  expect_true("I(des.n^2)" %in% f$fixed_terms)
  s <- mpa_formula_string(f)
  expect_match(s, "des.n:self.dist"); expect_match(s, "know.n:self.norm.c")
  expect_match(s, "I\\(des.n\\^2\\)"); expect_match(s, "des.n:know.n:self.dist")
  # dyad variable has no random slope at the dyad level; item variable none at the item level
  expect_false(grepl("know.n \\| perceiver:target", s))
  expect_true("know.n" %in% f$dropped$dyad)
  expect_true(all(c("self.norm.c", "des.n") %in% f$dropped$item))
  expect_match(s, "\\(1 \\+ self.dist \\+ know.n \\| item\\)")
  # uncorrelated random effects
  f2 <- mpa_formula(d, criterion = "self", correlated = FALSE)
  expect_match(mpa_formula_string(f2), "\\|\\| perceiver\\)")
})

test_that("dyadic random effects follow the design", {
  dy <- data.frame(perceiver = paste0("p", 1:6), target = rep(paste0("t", 1:3), each = 2),
                   item_1 = rnorm(6), item_2 = rnorm(6), stringsAsFactors = FALSE)
  d <- mpa_data(dy, item_cols = c("item_1", "item_2"), verbose = FALSE)
  expect_message(f <- mpa_formula(d, criterion = NULL), "Dyadic random effects omitted")
  expect_false(grepl("perceiver:target", mpa_formula_string(f)))
})

test_that("random_structure carries a pruned structure over and adds new slopes", {
  d <- make_toy_mpa()
  f1 <- mpa_update(mpa_formula(d, criterion = "self"), drop = list(item = "self.dist"))
  f1p <- mpa_update(f1, drop = list(perceiver = "self.dist", target = "self.norm.c"))
  f2 <- mpa_formula(d, criterion = "self", random_slopes = "des.n", interactions = "criterion", random_structure = f1p)
  expect_equal(f2$random$perceiver, c("1", "self.norm.c", "des.n"))
  expect_equal(f2$random$target, c("1", "self.dist", "des.n"))
  expect_equal(f2$random$item, "1")
  expect_equal(f2$random$dyad, c("1", "self.norm.c", "self.dist", "des.n"))
  f3 <- mpa_formula(d, criterion = "self", random_slopes = "know.n", random_structure = list(perceiver = "1", target = "1"))
  expect_equal(f3$random$perceiver, c("1", "know.n"))
  expect_false("dyad" %in% names(f3$random))
  # correlated is inherited from an mpa_formula template unless given explicitly
  f1u <- mpa_update(f1, correlated = FALSE)
  expect_false(mpa_formula(d, criterion = "self", random_slopes = "des.n", random_structure = f1u)$correlated)
  expect_true(mpa_formula(d, criterion = "self", random_slopes = "des.n", random_structure = f1u, correlated = TRUE)$correlated)
})

test_that("mpa_update edits random terms", {
  d <- make_toy_mpa()
  f <- mpa_formula(d, criterion = "self")
  f2 <- mpa_update(f, drop = list(dyad = "1", target = "self.norm.c"))
  s <- mpa_formula_string(f2)
  expect_match(s, "\\(0 \\+ self.norm.c \\+ self.dist \\| perceiver:target\\)")
  expect_match(s, "\\(1 \\+ self.dist \\| target\\)")
  f3 <- mpa_update(f, drop = list(dyad = c("1", "self.norm.c", "self.dist")))
  expect_false(grepl("perceiver:target", mpa_formula_string(f3)))
})

test_that("mpa_fit, mpa_compare and mpa_random_test run on toy data", {
  skip_if_not_installed("lme4")
  d <- make_toy_mpa()
  f0 <- mpa_formula(d, criterion = "self", random_units = c("perceiver", "target"))
  fit0 <- mpa_fit(f0, d)
  expect_s3_class(fit0, "mpa_fit")
  expect_true(all(c("term", "estimate", "se", "z", "p") %in% names(fit0$fixed)))
  expect_true(abs(fit0$fixed$estimate[fit0$fixed$term == "self.dist"] - 0.4) < 0.1)
  expect_true("Residual" %in% fit0$variance$unit)
  f1 <- mpa_formula(d, criterion = "self", predictors = "des.n", interactions = "criterion",
                    random_units = c("perceiver", "target"))
  fit1 <- mpa_fit(f1, d)
  cmp <- mpa_compare(fit0, fit1)
  expect_s3_class(cmp, "mpa_compare")
  expect_equal(cmp$lrt$df, 3)
  expect_output(print(cmp), "Likelihood-ratio test")
  rt <- suppressWarnings(mpa_random_test(mpa_formula(d, criterion = "self", random_units = "perceiver"), d))
  expect_s3_class(rt, "mpa_random_test")
  expect_equal(nrow(rt), 3)
  expect_true(all(rt$p_boundary <= rt$p, na.rm = TRUE))
  pruned <- mpa_prune(mpa_formula(d, criterion = "self", random_units = "perceiver"), rt, alpha = 0.05)
  expect_s3_class(pruned, "mpa_formula")
})

test_that("mpa_compare refits the reference model on common rows; mpa_plot works", {
  skip_if_not_installed("lme4")
  d <- make_toy_mpa()
  d2 <- as.data.frame(d); d2$know.n[d2$target %in% c("t1", "t2")] <- NA   # predictor with missing values
  d2 <- structure(d2, mpa = attr(d, "mpa"), class = class(d))
  f0 <- mpa_formula(d2, criterion = "self", random_units = "perceiver")
  f1 <- mpa_formula(d2, criterion = "self", predictors = "know.n", interactions = "criterion", random_units = "perceiver")
  m0 <- mpa_fit(f0, d2); m1 <- suppressWarnings(mpa_fit(f1, d2))
  expect_gt(m0$n_obs, m1$n_obs)
  expect_message(cmp <- mpa_compare(m0, m1), "refitted")
  expect_true(cmp$refitted)
  expect_equal(cmp$n_obs[1], m1$n_obs)   # after refit
  expect_false(is.null(cmp$lrt))
  m0b <- mpa_fit(f0, d2, same_rows_as = m1)
  expect_equal(m0b$n_obs, m1$n_obs)
  m1s <- mpa_fit(f1, d2, start_from = m1, check_derivs = FALSE)
  # variance table of a "||" fit uses the original grouping-factor names
  fu <- mpa_formula(d2, criterion = "self", correlated = FALSE)
  mu <- suppressWarnings(mpa_fit(fu, d2))
  expect_true(all(mu$variance$unit %in% c("perceiver", "target", "perceiver:target", "item", "Residual")))
  expect_equal(mu$variance$effect[mu$variance$unit == "perceiver"][1], "(Intercept)")
  expect_equal(m1s$logLik, m1$logLik, tolerance = 1e-4)
  m_other <- mpa_fit(mpa_formula(d2, criterion = "self", random_units = c("perceiver", "target")), d2)
  m_from_other <- suppressWarnings(mpa_fit(f1, d2, start_from = m_other))   # perceiver block matches, target block is dropped
  expect_s3_class(m_from_other, "mpa_fit")
  expect_true(inherits(stats::update(m1$fit, REML = TRUE), "merMod"))   # call is self-contained
  skip_if_not_installed("emmeans"); skip_if_not_installed("ggplot2")
  p <- mpa_plot(m1, slope = "self.dist", x = "know.n")
  expect_s3_class(p, "ggplot")
  pl <- mpa_plots(m1)
  expect_true(all(c("self.dist_by_know.n", "self.norm.c_by_know.n") %in% names(pl)))
  expect_equal(attr(pl[[1]], "slope"), "self.dist")
  al <- mpa_align_y(list(a = mpa_plot(m1, "self.dist", "know.n"), b = mpa_plot(m1, "self.dist", "know.n", n = 10)))
  expect_equal(al$a$coordinates$limits$y, al$b$coordinates$limits$y)
  expect_error(mpa_plot(m1, "self.dist", "know.n", by = "des.n"), "not a predictor")
})

test_that("the criterion source defaults to the one in the data", {
  d <- make_toy_mpa()
  expect_equal(mpa_formula(d)$criterion, "self")
  expect_equal(mpa_formula(d, criterion = "self")$formula, mpa_formula(d)$formula)
  expect_null(mpa_formula(d, criterion = NULL)$criterion)
})

test_that("random item slopes for the distinctive component", {
  d <- make_toy_mpa()
  f <- mpa_formula(d, criterion = "self")
  expect_equal(f$random$item, c("1", "self.dist"))
  expect_equal(f$dropped$item, "self.norm.c")
  f2 <- mpa_formula(d, criterion = "self", random_slopes = "des.n", random_structure = mpa_update(f, correlated = FALSE))
  expect_equal(f2$random$item, c("1", "self.dist"))        # an item property gets no item slope
  expect_equal(f2$dropped$item, "des.n")
})

test_that("start_from carries matching random terms over to a changed structure", {
  skip_if_not_installed("lme4")
  d <- make_toy_mpa()
  f1 <- mpa_formula(d, criterion = "self", random_units = c("perceiver", "target"), correlated = FALSE)
  m1 <- suppressWarnings(mpa_fit(f1, d))
  f2 <- mpa_formula(d, criterion = "self", random_slopes = "des.n", random_units = c("perceiver", "target"),
                    correlated = FALSE)
  expect_message(m2 <- suppressWarnings(mpa_fit(f2, d, start_from = m1)), "taken over for 6 of 8")
  expect_s3_class(m2, "mpa_fit")
  msgs <- capture_messages(suppressWarnings(mpa_fit(f1, d, start_from = m1)))   # identical structure: silent
  expect_false(any(grepl("start_from", msgs)))
})

test_that("mpa_compare warns when the models are not nested", {
  skip_if_not_installed("lme4")
  d <- make_toy_mpa()
  m_a <- mpa_fit(mpa_formula(d, criterion = "self", predictors = "des.n", random_units = "perceiver"), d)
  m_b <- mpa_fit(mpa_formula(d, criterion = "self", predictors = "know.n", random_units = "perceiver"), d)
  expect_warning(mpa_compare(m_a, m_b), "not nested")
})

test_that("perceiver x item random effects (random_units = \"perceiver:item\")", {
  d <- make_toy_mpa()
  f <- mpa_formula(d, criterion = "self", random_slopes = "know", random_units = c("perceiver", "target", "item", "perceiver:item"))
  expect_equal(f$random$perceiver_item, c("1", "self.dist"))        # intercept + distinctive slope, no other slopes
  expect_match(mpa_formula_string(f), "self.dist \\| perceiver:item\\)")
  expect_false("perceiver_item" %in% names(f$dropped))
  f2 <- mpa_update(mpa_formula(d, criterion = "self"), add = list(`perceiver:item` = "1"))
  expect_match(mpa_formula_string(f2), "\\(1 \\| perceiver:item\\)")
  expect_false(grepl("perceiver:item", mpa_formula_string(mpa_update(f2, drop = list(perceiver_item = "1")))))
  m <- mpa_fit(mpa_formula(d, criterion = "self", random_units = c("perceiver", "target", "perceiver:item"), correlated = FALSE), d)
  expect_true("perceiver:item" %in% m$variance$unit)
  expect_setequal(m$variance$effect[m$variance$unit == "perceiver:item"], c("(Intercept)", "self.dist"))
  expect_s3_class(mpa_prune(m), "mpa_formula")
  # omitted with a message when every perceiver judges a single target
  long <- expand.grid(perceiver = paste0("p", 1:8), item = paste0("i", 1:5), stringsAsFactors = FALSE)
  long$target <- sub("p", "t", long$perceiver)
  set.seed(8); long$rating <- round(runif(nrow(long), 1, 7))
  crit <- data.frame(target = rep(paste0("t", 1:8), each = 5), item = rep(paste0("i", 1:5), 8), value = round(runif(40, 1, 7)))
  d1 <- suppressMessages(mpa_data(judgments = long, format = "long", rating_col = "rating",
                                  criterion = list(self = crit), rating_scale = c(1, 7), verbose = FALSE))
  expect_message(f1 <- mpa_formula(d1, criterion = "self", random_units = c("perceiver", "target", "perceiver:item")),
                 "Perceiver x item random effects omitted")
  expect_false("perceiver_item" %in% names(f1$random))
})

test_that("target x item random intercepts (random_units = \"target:item\")", {
  d <- make_toy_mpa()
  f <- mpa_formula(d, criterion = "self", random_units = c("perceiver", "target", "dyad", "item", "target:item"))
  expect_equal(f$random$target_item, "1")
  expect_match(mpa_formula_string(f), "\\(1 \\| target:item\\)")
  expect_false("target_item" %in% names(f$dropped))
  # the internal name is accepted too, and the unit never gets random slopes
  f2 <- mpa_formula(d, criterion = "self", random_slopes = "know", random_units = c("perceiver", "target_item"))
  expect_equal(f2$random$target_item, "1")
  expect_equal(f2$random$perceiver, c("1", "self.norm.c", "self.dist", "know"))
  # mpa_update() with either spelling; a structure carried over keeps the unit
  f3 <- mpa_update(mpa_formula(d, criterion = "self"), add = list(`target:item` = "1"))
  expect_match(mpa_formula_string(f3), "\\(1 \\| target:item\\)")
  f4 <- mpa_formula(d, criterion = "self", random_slopes = "des", random_structure = f3)
  expect_equal(f4$random$target_item, "1")
  expect_match(mpa_formula_string(f4), "\\(1 \\| target:item\\)")
  f5 <- mpa_update(f3, drop = list(target_item = "1"))
  expect_false(grepl("target:item", mpa_formula_string(f5)))
  # fit, variance table and pruning work with the new unit
  m <- mpa_fit(mpa_formula(d, criterion = "self", random_units = c("perceiver", "target", "target:item"), correlated = FALSE), d)
  expect_true("target:item" %in% m$variance$unit)
  expect_equal(m$variance$effect[m$variance$unit == "target:item"], "(Intercept)")
  expect_s3_class(mpa_prune(m), "mpa_formula")
  # omitted with a message when every target is judged by a single perceiver
  long <- expand.grid(perceiver = paste0("p", 1:8), item = paste0("i", 1:5), stringsAsFactors = FALSE)
  long$target <- sub("p", "t", long$perceiver)
  set.seed(7); long$rating <- round(runif(nrow(long), 1, 7))
  crit <- data.frame(target = rep(paste0("t", 1:8), each = 5), item = rep(paste0("i", 1:5), 8), value = round(runif(40, 1, 7)))
  d1 <- suppressMessages(mpa_data(judgments = long, format = "long", rating_col = "rating",
                                  criterion = list(self = crit), rating_scale = c(1, 7), verbose = FALSE))
  expect_message(f1 <- mpa_formula(d1, criterion = "self", random_units = c("perceiver", "target", "target:item")),
                 "single perceiver")
  expect_false("target_item" %in% names(f1$random))
})
