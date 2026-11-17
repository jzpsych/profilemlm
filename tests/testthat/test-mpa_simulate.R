test_that("designs are generated and detected correctly", {
  d <- mpa_simulate(n_perceiver = 12, n_target = 12, n_item = 5, group_size = 3, seed = 1)
  ds <- attr(d, "mpa")$design
  expect_match(ds$type, "half-block")
  expect_equal(ds$n_group, 4); expect_equal(ds$n_dyad, 4 * 9)
  expect_equal(nrow(d), 36 * 5)
  expect_true("des" %in% attr(d, "mpa")$vars$item); expect_true("know" %in% attr(d, "mpa")$vars$dyad)

  d1 <- mpa_simulate(n_target = 6, n_item = 4, design = "one-with-many", group_size = 3, seed = 2)
  expect_match(attr(d1, "mpa")$design$type, "one-with-many")
  expect_false(attr(d1, "mpa")$design$dyad_re_identifiable)
  expect_equal(attr(d1, "mpa")$design$n_perceiver, 18)
  expect_true("know" %in% attr(d1, "mpa")$vars$perceiver)   # one-with-many: dyad-level == perceiver-level

  d2 <- mpa_simulate(n_target = 8, n_item = 4, design = "round-robin", group_size = 4, seed = 3)
  expect_match(attr(d2, "mpa")$design$type, "round-robin")
  expect_equal(attr(d2, "mpa")$design$n_dyad, 2 * 4 * 3)
  expect_false(attr(d2, "mpa")$design$self_ratings)
})

test_that("objects output, pairs argument and perceiver/target moderators", {
  ob <- mpa_simulate(n_perceiver = 12, n_target = 12, n_item = 5, group_size = 3, seed = 1, output = "objects",
                     fixed = c(intercept = 51, norm = 0.7, dist = 0.2, training_dist = 0.1, talkativeness = 3))
  expect_equal(names(ob), c("judgments", "criterion", "items", "perceivers", "targets", "dyads"))
  expect_equal(nrow(ob$judgments), 36); expect_equal(nrow(ob$dyads), 36); expect_equal(nrow(ob$criterion), 12)
  expect_false("know" %in% names(ob$judgments))
  d <- with(ob, mpa_data(judgments = judgments, dyads = dyads, targets = targets, perceivers = perceivers, items = items,
                         item_cols = "^item_", criterion = list(self = criterion), group = "group", verbose = FALSE))
  d <- mpa_sim_attach(d, ob)
  m <- attr(d, "mpa")
  expect_equal(m$vars$perceiver, "training"); expect_equal(m$vars$target, "talkativeness"); expect_equal(m$vars$dyad, "know")
  par <- mpa_sim_parameters(d)
  expect_equal(par$fixed$true[par$fixed$term == "self.dist:training"], 0.1)
  expect_equal(par$fixed$true[par$fixed$term == "talkativeness"], 3)
  d2 <- mpa_simulate(n_perceiver = 12, n_target = 12, n_item = 5, group_size = 3, seed = 1,
                     fixed = c(intercept = 51, norm = 0.7, dist = 0.2, training_dist = 0.1, talkativeness = 3))
  expect_equal(as.data.frame(d)$y, as.data.frame(d2)$y)      # same data either way
  tmp <- file.path(tempdir(), "mpa_sim_files")
  expect_message(mpa_simulate(n_perceiver = 12, n_target = 12, n_item = 5, group_size = 3, seed = 1, write_dir = tmp), "written")
  expect_true(all(file.exists(file.path(tmp, paste0(c("judgments", "criterion", "items", "perceivers", "targets", "dyads"), ".csv")))))
})

test_that("study layout (loop format) reproduces the same data set", {
  st <- mpa_simulate(n_perceiver = 12, n_target = 12, n_item = 5, group_size = 3, seed = 1, output = "study")
  expect_equal(names(st), c("perceiver_survey", "target_survey", "items"))
  expect_equal(nrow(st$perceiver_survey), 12)
  expect_true(all(c("t1_id", "t1_know", "t1_item_01", "t3_item_05", "training", "group") %in% names(st$perceiver_survey)))
  expect_false("t4_id" %in% names(st$perceiver_survey))
  d <- mpa_data(judgments = st$perceiver_survey, format = "loop", loop_id = "id",
                targets = st$target_survey, items = st$items, item_cols = "^item_",
                criterion = list(self = "targets"), group = "group", verbose = FALSE)
  m <- attr(d, "mpa")
  expect_equal(m$design$n_dyad, 36)
  expect_equal(m$vars$perceiver, "training")        # detected: constant within perceivers
  expect_equal(m$vars$target, "talkativeness")
  expect_equal(m$vars$dyad, "know")
  d2 <- mpa_simulate(n_perceiver = 12, n_target = 12, n_item = 5, group_size = 3, seed = 1)
  expect_equal(as.data.frame(d)$y, as.data.frame(d2)$y)
  expect_equal(as.data.frame(d)$know, as.data.frame(d2)$know)
  expect_equal(as.data.frame(d)$self.dist, as.data.frame(d2)$self.dist)
})

test_that("criterion_scale is reflected in the generating parameters", {
  d <- mpa_simulate(n_perceiver = 12, n_target = 12, n_item = 5, group_size = 3, seed = 1, criterion_scale = 10)
  par <- mpa_sim_parameters(d, scale = "analysis")
  expect_equal(par$fixed$true[par$fixed$term == "self.norm.c"], 7)
  expect_equal(par$variance$true_sd[par$variance$unit == "perceiver" & par$variance$effect == "self.norm.c"], 2.5)
  expect_equal(mpa_sim_parameters(d)$fixed$true[2], 0.7)
})

test_that("parameters are stored and criterion is recovered", {
  skip_if_not_installed("lme4")
  d <- mpa_simulate(n_perceiver = 40, n_target = 40, n_item = 20, group_size = 4, seed = 11,
                    fixed = c(intercept = 51, norm = 0.7, dist = 0.2, des = 5, know_dist = 0.1),
                    sd = list(perceiver = c(intercept = 8, norm = 0.15, dist = 0.05), target = c(intercept = 2),
                              dyad = c(intercept = 0), item = c(intercept = 2)),
                    resid_sd = 10)
  par <- mpa_sim_parameters(d)
  expect_equal(par$fixed$true[par$fixed$term == "self.dist:know"], 0.1)
  f <- mpa_formula(d, criterion = "self", predictors = c("des", "know"), interactions = "criterion",
                   random_units = c("perceiver", "target", "item"), correlated = FALSE)
  m <- mpa_fit(f, d)
  rec <- mpa_recovery(m)
  expect_s3_class(rec, "mpa_recovery")
  est <- function(term) rec$fixed$estimate[rec$fixed$term == term]
  expect_lt(abs(est("self.norm.c") - 0.7), 0.15)
  expect_lt(abs(est("self.dist") - 0.2), 0.08)
  expect_lt(abs(est("des") - 5), 4)
  expect_gt(mean(rec$fixed$covered), 0.6)
  expect_true(all(c("true_sd", "estimated_sd") %in% names(rec$variance)))
  expect_output(print(rec), "fixed effects")
})

test_that("scaled estimation is reported per criterion point; pruning by singularity", {
  skip_if_not_installed("lme4")
  sdl <- list(perceiver = c(intercept = 8, norm = 0.15), target = c(intercept = 0, dist = 0), item = c(intercept = 2))
  d1  <- mpa_simulate(n_perceiver = 40, n_target = 40, n_item = 20, group_size = 4, seed = 5, sd = sdl, resid_sd = 10, criterion_scale = 1)
  d10 <- mpa_simulate(n_perceiver = 40, n_target = 40, n_item = 20, group_size = 4, seed = 5, sd = sdl, resid_sd = 10, criterion_scale = 10)
  expect_equal(d1$self.dist, d10$self.dist * 10)
  f1  <- mpa_formula(d1,  criterion = "self", random_units = c("perceiver", "target", "item"), correlated = FALSE)
  f10 <- mpa_formula(d10, criterion = "self", random_units = c("perceiver", "target", "item"), correlated = FALSE)
  m1  <- suppressWarnings(mpa_fit(f1, d1)); m10 <- suppressWarnings(mpa_fit(f10, d10))
  expect_equal(m1$logLik, m10$logLik, tolerance = 1e-5)
  expect_equal(m1$fixed$estimate, m10$fixed$estimate, tolerance = 1e-2)       # back-transformed
  expect_equal(m1$variance$sd, m10$variance$sd, tolerance = 5e-2)
  m10a <- suppressWarnings(mpa_fit(f10, d10, report = "analysis"))
  expect_equal(m10a$fixed$estimate[m10a$fixed$term == "self.norm.c"],
               10 * m10$fixed$estimate[m10$fixed$term == "self.norm.c"], tolerance = 1e-6)
  ## target random effects were simulated as zero -> singular -> pruned away
  expect_true(m10$singular)
  expect_message(fp <- mpa_prune(m10), "near-\\)zero")
  expect_false(isTRUE("1" %in% fp$random$target))
  expect_true(all(c("1", "self.norm.c") %in% fp$random$perceiver))
  ## auto scaling uses half the rating range
  da <- mpa_simulate(n_perceiver = 12, n_target = 12, n_item = 5, group_size = 3, seed = 1, rating_scale = c(1, 101))
  expect_equal(attr(da, "mpa")$settings$criterion_scale, 50)
  expect_true(all(da$y >= 1 & da$y <= 101))
})

test_that("item random slopes for the distinctive component are simulated and reported", {
  d <- mpa_simulate(n_perceiver = 12, n_target = 12, n_item = 5, group_size = 3, seed = 1,
                    sd = list(item = c(intercept = 2, dist = 0.05, norm = 0.3)))
  par <- mpa_sim_parameters(d)
  expect_equal(par$variance$true_sd[par$variance$unit == "item" & par$variance$effect == "self.dist"], 0.05)
  expect_false(any(par$variance$unit == "item" & par$variance$effect == "self.norm.c"))   # ignored
  d0 <- mpa_simulate(n_perceiver = 12, n_target = 12, n_item = 5, group_size = 3, seed = 1,
                     sd = list(item = c(intercept = 2)))
  expect_false(isTRUE(all.equal(as.data.frame(d)$y, as.data.frame(d0)$y)))
})

test_that("a shared target x item component is simulated, reported and recovered", {
  skip_if_not_installed("lme4")
  args <- list(n_perceiver = 24, n_target = 24, n_item = 8, group_size = 4, seed = 2, resid_sd = 5)
  d <- do.call(mpa_simulate, c(args, list(sd = list(target_item = c(intercept = 4)))))
  par <- mpa_sim_parameters(d)
  expect_equal(par$variance$true_sd[par$variance$unit == "target:item"], 4)
  expect_equal(par$variance$effect[par$variance$unit == "target:item"], "(Intercept)")   # intercept only
  d2 <- do.call(mpa_simulate, c(args, list(sd = list("target:item" = c(intercept = 4, dist = 1)))))  # alias; slope ignored
  expect_equal(as.data.frame(d2)$y, as.data.frame(d)$y)
  d0 <- do.call(mpa_simulate, args)
  expect_false(any(mpa_sim_parameters(d0)$variance$unit == "target:item"))   # default 0: not listed
  expect_false(isTRUE(all.equal(as.data.frame(d0)$y, as.data.frame(d)$y)))
  expect_error(do.call(mpa_simulate, c(args, list(sd = list(group = c(intercept = 1))))), "Unknown unit")
  f <- mpa_formula(d, criterion = "self", random_units = c("perceiver", "target", "item", "target:item"),
                   correlated = FALSE)
  m <- mpa_fit(f, d)
  rec <- mpa_recovery(m)
  ti <- rec$variance[rec$variance$unit == "target:item", ]
  expect_equal(nrow(ti), 1)
  expect_equal(ti$true_sd, 4)
  expect_lt(abs(ti$estimated_sd - 4), 1.5)
  # perceiver x item: intercept and distinctive slope
  dp <- do.call(mpa_simulate, c(args, list(sd = list("perceiver:item" = c(intercept = 4, dist = 0.05, norm = 1)))))   # norm ignored
  pp <- mpa_sim_parameters(dp)$variance
  expect_equal(pp$effect[pp$unit == "perceiver:item"], c("(Intercept)", "self.dist"))
  expect_equal(pp$true_sd[pp$unit == "perceiver:item"], c(4, 0.05))
  expect_false(any(mpa_sim_parameters(d0)$variance$unit == "perceiver:item"))
  fp <- mpa_formula(dp, criterion = "self", random_units = c("perceiver", "target", "item", "perceiver:item"), correlated = FALSE)
  recp <- mpa_recovery(mpa_fit(fp, dp))
  pi_int <- recp$variance[recp$variance$unit == "perceiver:item" & recp$variance$effect == "(Intercept)", ]
  expect_equal(pi_int$true_sd, 4)
  expect_lt(abs(pi_int$estimated_sd - 4), 1.5)
})

test_that("the Wessels data sets are complete", {
  expect_equal(dim(wessels_judgments), c(1352 * 111, 8))
  expect_equal(nrow(wessels_items), 111)
  expect_false(anyNA(wessels_judgments)); expect_false(anyNA(wessels_items))
  expect_equal(as.integer(table(unique(wessels_judgments[, c("perceiver", "type")])$type)), c(943, 409))
  expect_equal(length(unique(wessels_judgments$target)), 189)
  expect_true(all(wessels_judgments$rating %in% 1:5))
})
