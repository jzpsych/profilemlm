## further predictor profiles, pairs design, effect sizes

sim_similarity <- function(...) {
  mpa_simulate(n_perceiver = 24, n_item = 12, group_size = 4, seed = 5, rating_scale = c(1, 101),
               fixed = c(intercept = 51, norm = 0.7, dist = 0.2, similarity = 0.15, training_similarity = 0.1), ...)
}

make_pairs <- function(n = 60, k = 10, seed = 3) {
  set.seed(seed)
  items <- sprintf("item_%02d", seq_len(k))
  norm_i <- rnorm(k, 3.4, 0.6)
  prof <- function() matrix(norm_i, n, k, byrow = TRUE) + matrix(rnorm(n * k, 0, 0.9), n, k)
  anchor <- prof(); ex <- prof(); z <- rnorm(n)
  dev <- function(x) sweep(x, 2, colMeans(x))
  current <- 3.4 + rnorm(n, 0, 0.25) + 0.5 * matrix(norm_i - mean(norm_i), n, k, byrow = TRUE) +
    (0.25 + rnorm(n, 0, 0.1) - 0.1 * z) * dev(ex) + (0.2 + rnorm(n, 0, 0.1)) * dev(anchor) + matrix(rnorm(n * k, 0, 0.8), n, k)
  colnames(current) <- colnames(anchor) <- colnames(ex) <- items
  id <- sprintf("a%03d", seq_len(n))
  list(partners = data.frame(anchor = id, ex = paste0("x", id), length = runif(n), current),
       exes = data.frame(ex = paste0("x", id), ex),
       anchors = data.frame(anchor = id, extraversion = z, anchor))
}

test_that("mpa_data() builds a perceiver profile from the perceivers file", {
  st <- sim_similarity(output = "objects")
  expect_true("perceiver_self" %in% names(st))
  perc <- merge(st$perceivers, st$perceiver_self, by = "perceiver")
  d <- mpa_data(judgments = st$judgments, dyads = st$dyads, targets = st$targets, perceivers = perc, items = st$items,
                item_cols = "^item_", criterion = list(self = st$criterion), profiles = list(pself = "perceivers"),
                group = "group", rating_scale = c(1, 101), verbose = FALSE)
  m <- attr(d, "mpa")
  expect_equal(m$vars$profiles, "pself")
  expect_equal(unname(m$vars$profile_key), "perceiver")
  expect_true(all(c("pself", "pself.norm", "pself.dist") %in% names(d)))
  expect_equal(m$vars$perceiver, "training")                 # item columns are not perceiver variables
  ## one value per perceiver x item, the leave-one-out deviation, scaled like the criterion
  one <- d[d$perceiver == "p001" & d$item == "item_03", ]
  expect_equal(length(unique(one$pself.dist)), 1)
  others <- st$perceiver_self$item_03[st$perceiver_self$perceiver != "p001"]
  expect_equal(one$pself.dist[1], (st$perceiver_self$item_03[1] - mean(others)) / 50)
  ## the same object from a data frame keyed by the perceiver
  d2 <- mpa_data(judgments = st$judgments, dyads = st$dyads, targets = st$targets, perceivers = st$perceivers, items = st$items,
                 item_cols = "^item_", criterion = list(self = st$criterion), profiles = list(pself = st$perceiver_self),
                 group = "group", rating_scale = c(1, 101), verbose = FALSE)
  expect_equal(d2$pself.dist, d$pself.dist)
  expect_output(print(d), "Further predictor profiles: pself \\(perceiver x item")
  expect_error(mpa_data(judgments = st$judgments, item_cols = "^item_", criterion = list(self = st$criterion),
                        profiles = list(self = st$perceiver_self), group = "group", verbose = FALSE), "already used")
  expect_error(mpa_data(judgments = st$judgments, item_cols = "^item_", profiles = list(pself = "perceivers"),
                        group = "group", verbose = FALSE), "perceivers")
})

test_that("a target profile can be a criterion source or a profile", {
  st <- sim_similarity(output = "objects")
  second <- st$criterion; second[-1] <- second[-1] + rnorm(prod(dim(second[-1])), 0, 5)
  d <- mpa_data(judgments = st$judgments, item_cols = "^item_", criterion = list(self = st$criterion),
                profiles = list(informant = second), group = "group", rating_scale = c(1, 101), verbose = FALSE)
  expect_equal(unname(attr(d, "mpa")$vars$profile_key), "target")
  f <- mpa_formula(d, criterion = "self", profiles = "informant")
  expect_equal(unname(f$profile_level), "target_item")
  d2 <- mpa_data(judgments = st$judgments, item_cols = "^item_", criterion = list(self = st$criterion, informant = second),
                 criterion_combine = "none", group = "group", rating_scale = c(1, 101), verbose = FALSE)
  f2 <- mpa_formula(d2, criterion = "self", profiles = "informant")
  expect_equal(f2$profile_dist, "informant.dist")
  expect_equal(d2$informant.dist, d$informant.dist)
})

test_that("mpa_formula() treats a further profile like the distinctive component", {
  d <- sim_similarity()
  f <- mpa_formula(d, criterion = "self", profiles = "pself", predictors = "training", interactions = "criterion")
  expect_equal(f$fixed_terms, c("self.norm.c", "self.dist", "pself.dist", "training",
                                "training:self.norm.c", "training:self.dist", "training:pself.dist"))
  for (u in c("perceiver", "target", "dyad", "item")) expect_true("pself.dist" %in% f$random[[u]])
  expect_equal(unname(f$profile_level), "perceiver_item")
  expect_output(print(f), "Further predictor profiles: pself")
  expect_error(mpa_formula(d, criterion = "self", predictors = "pself.dist"), "profiles")
  expect_error(mpa_formula(d, criterion = NULL, profiles = "pself"), "criterion")
  expect_error(mpa_formula(d, criterion = "self", profiles = "nothere"), "not found")
  ## carried-over structure: the profile slope is added where it varies
  base <- mpa_formula(d, criterion = "self", random_units = c("perceiver", "item"), correlated = FALSE)
  f2 <- mpa_formula(d, criterion = "self", profiles = "pself", random_structure = base)
  expect_true("pself.dist" %in% f2$random$perceiver && "pself.dist" %in% f2$random$item)
  expect_false(f2$correlated)
})

test_that("the slope of a further profile is reported per scale point", {
  d <- sim_similarity()
  f <- mpa_formula(d, criterion = "self", profiles = "pself", predictors = "training", interactions = "criterion",
                   random_units = c("perceiver", "item"), correlated = FALSE)
  m <- mpa_fit(f, d, check_derivs = FALSE)
  ma <- mpa_fit(f, d, check_derivs = FALSE, report = "analysis")
  expect_true("pself.dist" %in% m$criterion_vars)
  pick <- function(x, t) x$fixed$estimate[x$fixed$term == t]
  expect_equal(pick(m, "pself.dist"), pick(ma, "pself.dist") / 50)
  expect_equal(pick(m, "pself.dist:training"), pick(ma, "pself.dist:training") / 50)
  sdp <- function(x) x$variance$sd[x$variance$unit == "perceiver" & x$variance$effect == "pself.dist"]
  expect_equal(sdp(m), sdp(ma) / 50)
  rec <- mpa_recovery(m)
  expect_equal(rec$fixed$true[rec$fixed$term == "pself.dist"], 0.15)
  expect_equal(rec$fixed$true[rec$fixed$term == "pself.dist:training"], 0.1)
  skip_if_not_installed("ggplot2")
  pl <- mpa_plots(m)
  expect_true("pself.dist_by_training" %in% names(pl))
})

test_that("mpa_simulate() without assumed similarity is unchanged", {
  a <- mpa_simulate(n_perceiver = 8, n_item = 6, group_size = 4, seed = 1)
  expect_false("pself.dist" %in% names(a))
  expect_false(any(grepl("pself", mpa_sim_parameters(a)$fixed$term)))
  expect_equal(nrow(mpa_sim_parameters(a)$variance), 15)
  b <- sim_similarity()
  expect_true("pself.dist" %in% names(b))
  rr <- mpa_simulate(n_target = 8, n_item = 6, group_size = 4, design = "round-robin", seed = 1, fixed = c(similarity = 0.1))
  one <- rr[rr$perceiver == "s001" & rr$item == "item_01", ]
  expect_equal(one$pself[1], rr$self[rr$target == "s001" & rr$item == "item_01"][1])
})

test_that("a pairs design keeps one person-level unit", {
  p <- make_pairs()
  expect_message(
    d <- mpa_data(judgments = p$partners, targets = p$exes, perceivers = p$anchors, id = c(perceiver = "anchor", target = "ex"),
                  item_cols = "^item_", criterion = list(ex = "targets"), profiles = list(self = "perceivers"),
                  norm = "sample"),
    "Pairs design")
  des <- attr(d, "mpa")$design
  expect_true(des$pairs)
  expect_match(des$type, "^pairs")
  expect_true("length" %in% attr(d, "mpa")$vars$perceiver)        # a variable of the pair
  expect_message(f <- mpa_formula(d, criterion = "ex", profiles = "self", random_slopes = "length",
                                  predictors = "extraversion", interactions = "criterion", correlated = FALSE),
                 "Person-level random effects are kept for `perceiver` only")
  expect_equal(names(f$random), c("perceiver", "item"))
  expect_equal(f$random$perceiver, c("1", "ex.norm.c", "ex.dist", "self.dist"))   # no slope for the pair's variable
  expect_true("length" %in% f$random$item)
  ft <- suppressMessages(mpa_formula(d, criterion = "ex", random_units = c("target", "item")))
  expect_equal(names(ft$random), c("target", "item"))
  m <- mpa_fit(mpa_formula(d, criterion = "ex", profiles = "self", predictors = "extraversion", interactions = "criterion",
                           random_units = "perceiver", correlated = FALSE), d, check_derivs = FALSE)
  est <- stats::setNames(m$fixed$estimate, m$fixed$term)
  expect_lt(abs(est[["ex.dist"]] - 0.25), 0.08)
  expect_lt(abs(est[["self.dist"]] - 0.20), 0.08)
  es <- mpa_effect_size(m)
  expect_true(all(es$reference %in% c("perceiver", NA)))
  sd_ex <- m$variance$sd[m$variance$unit == "perceiver" & m$variance$effect == "ex.dist"]
  z <- p$anchors$extraversion
  expect_equal(es$d[es$term == "ex.dist:extraversion"], est[["ex.dist:extraversion"]] * 2 * sd(z) / sd_ex)
})

test_that("a single target per perceiver gives target variables no perceiver slope", {
  d <- mpa_simulate(n_target = 10, n_item = 6, group_size = 3, design = "one-with-many", seed = 2)
  f <- suppressMessages(mpa_formula(d, criterion = "self", random_slopes = c("talkativeness", "training")))
  expect_false("talkativeness" %in% f$random$perceiver)
  expect_false("training" %in% f$random$perceiver)
  expect_true("training" %in% f$random$target)
})

test_that("mpa_effect_size() standardizes by the SD at the predictor's level", {
  d <- mpa_simulate(n_perceiver = 40, n_item = 12, group_size = 4, seed = 7, rating_scale = c(1, 101),
                    fixed = c(intercept = 51, norm = 0.7, dist = 0.25, training_dist = 0.15, des_dist = 0.1, know_dist = 0.1),
                    sd = list(perceiver = c(intercept = 10, dist = 0.15), target = c(intercept = 4, dist = 0.1),
                              dyad = c(intercept = 6, dist = 0.1), item = c(intercept = 3, dist = 0.25)),
                    derive = list(female = ~ as.numeric(training > 0)))
  f <- mpa_formula(d, criterion = "self", predictors = c("training", "des", "know", "female"), interactions = "criterion",
                   correlated = FALSE)
  expect_equal(f$predictors$level[f$predictors$term == "female"], "perceiver")
  m <- mpa_fit(f, d, check_derivs = FALSE)
  es <- mpa_effect_size(m)
  expect_s3_class(es, "mpa_effect_size")
  fx <- stats::setNames(m$fixed$estimate, m$fixed$term)
  v <- function(u, e) m$variance$variance[m$variance$unit == u & m$variance$effect == e]
  tr <- d$training[!duplicated(d$perceiver)]
  expect_equal(es$d[es$term == "self.dist:training"], fx[["self.dist:training"]] * 2 * sd(tr) / sqrt(v("perceiver", "self.dist")))
  expect_equal(es$d[es$term == "training"], fx[["training"]] * 2 * sd(tr) / sqrt(v("perceiver", "(Intercept)")))
  ## two values: the difference, not two standard deviations
  expect_equal(es$d[es$term == "self.dist:female"], fx[["self.dist:female"]] / sqrt(v("perceiver", "self.dist")))
  ## item property: SD across items
  ds <- d$des[!duplicated(d$item)]
  expect_equal(es$d[es$term == "self.dist:des"], fx[["self.dist:des"]] * 2 * sd(ds) / sqrt(v("item", "self.dist")))
  expect_true(is.na(es$d[es$term == "self.norm.c:des"]))            # items have no normative slope
  ## a dyadic variable that was not decomposed: perceiver + target + dyad
  kn <- d$know[!duplicated(paste(d$perceiver, d$target))]
  vsum <- sum(c(v("perceiver", "self.dist"), v("target", "self.dist"), v("perceiver:target", "self.dist")))
  expect_equal(es$d[es$term == "self.dist:know"], fx[["self.dist:know"]] * 2 * sd(kn) / sqrt(vsum), tolerance = 1e-6)
  ## interval and total
  i <- which(es$term == "self.dist:training")
  expect_equal(es$ci_lower[i] / es$d[i], m$fixed$ci_lower[m$fixed$term == "self.dist:training"] / fx[["self.dist:training"]])
  tot <- mpa_effect_size(m, denominator = "total")
  expect_true(all(abs(tot$d) <= abs(es$d) + 1e-12 | is.na(es$d)))
  expect_output(print(es), "mpa_effect_size")
  expect_error(mpa_effect_size(mpa_fit(mpa_formula(d, criterion = "self", random_units = "perceiver"), d)), "no predictors")
})

test_that("mpa_effect_size() sums perceiver and target variance in a one-with-many design", {
  d <- mpa_simulate(n_target = 40, n_item = 16, group_size = 5, design = "one-with-many", seed = 4, rating_scale = c(1, 101),
                    fixed = c(intercept = 51, norm = 0.7, dist = 0.25, training_dist = 0.15), resid_sd = 10,
                    sd = list(perceiver = c(intercept = 8, dist = 0.2), target = c(intercept = 5, dist = 0.25),
                              item = c(intercept = 3)))
  f <- suppressMessages(mpa_formula(d, criterion = "self", predictors = c("training", "talkativeness"), interactions = "criterion",
                                    random_units = c("perceiver", "target"), correlated = FALSE))
  m <- mpa_fit(f, d, check_derivs = FALSE)
  es <- mpa_effect_size(m)
  v <- function(u, e) m$variance$variance[m$variance$unit == u & m$variance$effect == e]
  fx <- stats::setNames(m$fixed$estimate, m$fixed$term)
  tr <- d$training[!duplicated(d$perceiver)]
  expect_equal(es$reference[es$term == "self.dist:training"], "perceiver + target")
  expect_equal(es$d[es$term == "self.dist:training"],
               fx[["self.dist:training"]] * 2 * sd(tr) / sqrt(v("perceiver", "self.dist") + v("target", "self.dist")))
  expect_equal(es$reference[es$term == "self.dist:talkativeness"], "target")
})
