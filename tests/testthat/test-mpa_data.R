## Tests for mpa_data(). The reproduction test needs the raw study files and is
## skipped when they are not available; set the environment variable
## MPA_TESTDATA to the folder containing data_MPA_neu.xlsx, Itemsample.xlsx and
## data_restruct.csv to run it.

## ---- small synthetic example --------------------------------------------

make_toy <- function() {
  set.seed(1)
  perceivers <- paste0("p", 1:4); targets <- paste0("t", 1:3); items <- paste0("item_", 1:5)
  dy <- expand.grid(perceiver = perceivers, target = targets, stringsAsFactors = FALSE)
  for (i in items) dy[[i]] <- round(runif(nrow(dy), 0, 100))
  dy$know <- runif(nrow(dy), 1, 101)
  dy$role <- "perceiver"
  self <- data.frame(perceiver = targets, target = targets, stringsAsFactors = FALSE)
  for (i in items) self[[i]] <- round(runif(3, 0, 100))
  self$know <- NA; self$role <- "self"
  list(dyads = rbind(dy, self),
       items = data.frame(item = items, des = c(10, 30, 50, 70, 90)),
       targets = data.frame(target = targets, sex = c("f", "m", "f")))
}

test_that("basic structure and design detection", {
  toy <- make_toy()
  d <- mpa_data(toy$dyads, items = toy$items, targets = toy$targets,
                item_cols = "^item_", role = "role", judgment_role = "perceiver",
                criterion = list(self = "self"), norm = "sample", rating_scale = c(0, 100), criterion_scale = 1,
                rescale = list(des = c(0, 100), know = c(1, 101)), verbose = FALSE)
  expect_s3_class(d, "mpa_data")
  expect_equal(nrow(d), 4 * 3 * 5)
  expect_true(all(c("perceiver", "target", "item", "item_index", "y",
                    "self", "self.norm", "self.norm.c", "self.dist", "des", "des.n", "know.n", "sex") %in% names(d)))
  m <- attr(d, "mpa")
  expect_equal(m$design$type, "half-block (fully crossed)")
  expect_true(m$design$dyad_re_identifiable)
  expect_equal(m$vars$item, c("des", "self.sd", "des.n"))
  expect_equal(m$vars$dyad, c("know", "know.n"))
  expect_equal(m$vars$target, "sex")
  expect_equal(m$vars$sd, "self.sd")
  # normative profile = mean over targets per item, distinctive = criterion - norm
  self_wide <- toy$dyads[toy$dyads$role == "self", ]
  expect_equal(unique(d$self.norm[d$item == "item_1"]), mean(self_wide$item_1))
  expect_equal(d$self.dist, d$self - d$self.norm)
  expect_equal(mean(unique(d[, c("item", "self.norm.c")])$self.norm.c), 0)
  # item SD of the criterion = SD over targets per item, one value per item
  expect_equal(unique(d$self.sd[d$item == "item_2"]), sd(self_wide$item_2))
  expect_equal(nrow(unique(d[, c("item", "self.sd")])), 5)
  # rescaling maps scale range to -1..1
  expect_equal(range(d$des.n), c(-0.8, 0.8))
})

test_that("level detection registers perceiver- and target-level columns of the pair file", {
  toy <- make_toy()
  dy <- toy$dyads
  dy$age <- c(30, 40, 50, 60)[match(dy$perceiver, paste0("p", 1:4))]      # perceiver-level (NA for self rows)
  dy$age[dy$role == "self"] <- NA
  dy$sex <- c("f", "m", "f")[match(dy$target, paste0("t", 1:3))]           # target-level
  d <- mpa_data(dy, item_cols = "^item_", role = "role", judgment_role = "perceiver",
                criterion = list(self = "self"), verbose = FALSE)
  m <- attr(d, "mpa")
  expect_equal(m$vars$perceiver, "age"); expect_equal(m$vars$target, "sex"); expect_equal(m$vars$dyad, "know")
  d0 <- mpa_data(dy, item_cols = "^item_", role = "role", judgment_role = "perceiver",
                 criterion = list(self = "self"), detect_levels = FALSE, verbose = FALSE)
  expect_true(all(c("age", "sex", "know") %in% attr(d0, "mpa")$vars$dyad))
})

test_that("derived variables inherit the level", {
  toy <- make_toy()
  d <- mpa_data(toy$dyads, items = toy$items, item_cols = "^item_", role = "role", judgment_role = "perceiver",
                criterion = list(self = "self"), rescale = list(des = c(0, 100)),
                derive = list(eval = ~ abs(des.n), know2 = function(x) x$know^2), verbose = FALSE)
  expect_equal(d$eval, abs(d$des.n))
  expect_true("eval" %in% attr(d, "mpa")$vars$item)
  expect_equal(attr(d, "mpa")$vars$derived, c("eval", "know2"))
})

test_that("the item SD of the criterion is an item-level variable usable in formulas", {
  toy <- make_toy()
  d <- mpa_data(toy$dyads, items = toy$items, item_cols = "^item_", role = "role", judgment_role = "perceiver",
                criterion = list(self = "self"), rescale = list(des = c(0, 100)),
                derive = list(sd.c = ~ self.sd - mean(self.sd)), verbose = FALSE)
  expect_true("sd.c" %in% attr(d, "mpa")$vars$item)
  f <- mpa_formula(d, criterion = "self", predictors = c("des.n", "sd.c"), interactions = "criterion",
                   correlated = FALSE)
  s <- mpa_formula_string(f)
  expect_true(grepl("sd.c:self.dist", s, fixed = TRUE))
  item_term <- regmatches(s, regexpr("\\([^()]*\\|\\| item\\)", s))
  expect_false(grepl("sd.c", item_term, fixed = TRUE))        # one value per item: no item slope
  expect_output(summary(d), "item SD")
})

test_that("automatic criterion_scale is half the rating range", {
  toy <- make_toy()
  d <- mpa_data(toy$dyads, item_cols = "^item_", role = "role", judgment_role = "perceiver",
                criterion = list(self = "self"), rating_scale = c(0, 100), verbose = FALSE)
  expect_equal(attr(d, "mpa")$settings$criterion_scale, 50)
  expect_equal(d$self.dist * 50, d$self - d$self.norm)
  d0 <- mpa_data(toy$dyads, item_cols = "^item_", role = "role", judgment_role = "perceiver",
                 criterion = list(self = "self"), verbose = FALSE)              # no rating_scale -> 1
  expect_equal(attr(d0, "mpa")$settings$criterion_scale, 1)
})

test_that("dyadic decomposition and criterion_scale", {
  toy <- make_toy()
  d <- mpa_data(toy$dyads, item_cols = "^item_", role = "role", judgment_role = "perceiver",
                criterion = list(self = "self"), decompose = "know", criterion_scale = 10, verbose = FALSE)
  expect_true(all(c("know.p", "know.t", "know.r") %in% names(d)))
  expect_equal(d$know.p + d$know.t + d$know.r + mean(unique(d[, c("perceiver", "target", "know")])$know), d$know)
  m <- attr(d, "mpa")
  expect_true("know.p" %in% m$vars$perceiver); expect_true("know.t" %in% m$vars$target); expect_true("know.r" %in% m$vars$dyad)
  # component means are zero over dyads
  dy <- unique(d[, c("perceiver", "target", "know.p", "know.t", "know.r")])
  expect_equal(mean(dy$know.p), 0); expect_equal(mean(dy$know.t), 0); expect_equal(mean(dy$know.r), 0)
  d1 <- mpa_data(toy$dyads, item_cols = "^item_", role = "role", judgment_role = "perceiver",
                 criterion = list(self = "self"), verbose = FALSE)
  expect_equal(d$self.dist, d1$self.dist / 10)
  expect_equal(d$self.norm, d1$self.norm)          # descriptive variables unchanged
  # identifiability: the perceiver component gets no random slope at the perceiver level
  f <- mpa_formula(d, criterion = "self", random_slopes = c("know.p", "know.t", "know.r"))
  expect_false("know.p" %in% f$random$perceiver); expect_true("know.p" %in% f$random$target)
  expect_false("know.t" %in% f$random$target);    expect_true("know.t" %in% f$random$perceiver)
  expect_false("know.r" %in% f$random$dyad);      expect_true("know.r" %in% f$random$perceiver)
})

test_that("leave-one-out norm excludes the own target", {
  toy <- make_toy()
  d <- mpa_data(toy$dyads, item_cols = "^item_", role = "role", judgment_role = "perceiver",
                criterion = list(self = "self"), norm = "loo", verbose = FALSE)
  self_wide <- toy$dyads[toy$dyads$role == "self", ]
  loo_t1 <- mean(self_wide$item_1[self_wide$target != "t1"])
  expect_equal(unique(d$self.norm[d$item == "item_1" & d$target == "t1"]), loo_t1)
})

test_that("criterion from the targets object, external norm and midpoint centering", {
  toy <- make_toy()
  dy <- toy$dyads[toy$dyads$role == "perceiver", ]
  tg <- toy$dyads[toy$dyads$role == "self", c("target", paste0("item_", 1:5))]
  ext <- data.frame(item = paste0("item_", 1:5), norm = c(20, 40, 60, 80, 100))
  d <- mpa_data(dy, targets = tg, item_cols = "^item_", criterion = list(self = "targets"),
                norm = ext, center_norm = "midpoint", rating_scale = c(0, 100), criterion_scale = 1, verbose = FALSE)
  expect_equal(unique(d$self.norm[d$item == "item_3"]), 60)
  expect_equal(unique(d$self.norm.c[d$item == "item_3"]), 10)
})

test_that("one-with-many design is detected and dyadic REs flagged", {
  dy <- data.frame(perceiver = paste0("p", 1:6), target = rep(paste0("t", 1:3), each = 2),
                   item_1 = 1:6, item_2 = 6:1, stringsAsFactors = FALSE)
  d <- mpa_data(dy, item_cols = c("item_1", "item_2"), verbose = FALSE)
  expect_match(attr(d, "mpa")$design$type, "one-with-many")
  expect_false(attr(d, "mpa")$design$dyad_re_identifiable)
})

test_that("informative errors", {
  toy <- make_toy()
  expect_error(mpa_data(toy$dyads, item_cols = "^nothing"), "No columns")
  expect_error(mpa_data(toy$dyads, item_cols = "^item_", criterion = list(self = "self")), "`role` is not specified")
  expect_error(mpa_data(toy$dyads, item_cols = "^item_", role = "role", criterion = list(self = "nope")), "not found")
})

## ---- reproduction of prepare.txt / data_restruct.csv ---------------------

test_that("mpa_data reproduces data_restruct.csv", {
  dir <- Sys.getenv("MPA_TESTDATA", unset = "")
  skip_if(dir == "" || !file.exists(file.path(dir, "data_MPA_neu.xlsx")), "raw study data not available")
  skip_if_not_installed("readxl")

  raw <- as.data.frame(readxl::read_excel(file.path(dir, "data_MPA_neu.xlsx")))
  raw$Target <- as.character(raw$Target); raw$Role <- as.character(raw$Role)
  raw <- raw[grepl("^[0-9]+$", raw$Target), ]
  raw$role <- ifelse(grepl("^[0-9]+$", raw$Role), "self", ifelse(grepl("^TNI", raw$Role), "TNI", "perceiver"))
  items <- as.data.frame(readxl::read_excel(file.path(dir, "Itemsample.xlsx")))[-c(26, 67), ]

  ref <- read.csv2(file.path(dir, "data_restruct.csv")); ref$Target <- as.character(ref$Target)
  key_ref <- paste(ref$Role, ref$Target, ref$Item)

  build <- function(self_spec) {
    mpa_data(raw, items = items, id = c(perceiver = "Role", target = "Target"), item_cols = "^item_",
             role = "role", judgment_role = "perceiver",
             criterion = list(self = self_spec, other = "TNI"), norm = "sample", center_norm = "grand_mean", criterion_scale = 1,
             rename = c(DES = "mean_Pos", KNOW = "SD_1", LIKE = "SD_2"),
             rescale = list(DES = c(0, 100), KNOW = c(1, 101), LIKE = c(1, 101)),
             group = "Group", verbose = FALSE)
  }
  d1 <- build(list(role = "self", max_n = 1))
  idx <- match(key_ref, paste(as.character(d1$perceiver), as.character(d1$target), d1$item_index))
  expect_false(anyNA(idx)); expect_equal(nrow(d1), nrow(ref))
  expect_equal(d1$y[idx], ref$value)
  for (v in c("self", "self.norm", "self.dist", "self.norm.c", "other", "other.norm", "other.dist", "other.norm.c",
              "DES.n", "KNOW.n", "LIKE.n"))
    expect_equal(d1[[v]][idx], ref[[v]], tolerance = 1e-8, info = v)
  expect_equal(attr(d1, "mpa")$design$n_group, 32)
  expect_match(attr(d1, "mpa")$design$type, "half-block")

  d0 <- build("self")   # combined.* in the reference were computed before removing multiple self-ratings
  idx0 <- match(key_ref, paste(as.character(d0$perceiver), as.character(d0$target), d0$item_index))
  for (v in c("combined", "combined.norm", "combined.dist", "combined.norm.c"))
    expect_equal(d0[[v]][idx0], ref[[v]], tolerance = 1e-8, info = v)
})

test_that("dyadic decomposition is carried out within groups", {
  st <- mpa_simulate(n_perceiver = 12, n_target = 12, n_item = 5, group_size = 4, seed = 7, output = "study")
  ps <- st$perceiver_survey
  shift <- c(-1, 0, 1)[ps$group]                       # group-level differences in acquaintance
  kc <- grep("_know$", names(ps)); ps[kc] <- ps[kc] + shift
  d <- mpa_data(judgments = ps, format = "loop", loop_id = "id", targets = st$target_survey, items = st$items,
                item_cols = "^item_", criterion = list(self = "targets"), group = "group",
                decompose = "know", verbose = FALSE)
  dy <- unique(as.data.frame(d)[, c("perceiver", "target", "group", "know", "know.g", "know.p", "know.t", "know.r")])
  expect_equal(dy$know.g + dy$know.p + dy$know.t + dy$know.r + mean(dy$know), dy$know)
  expect_equal(as.numeric(tapply(dy$know.g, dy$group, mean)), as.numeric(tapply(dy$know, dy$group, mean)) - mean(dy$know))
  # within every group the perceiver, target and relationship components average to zero
  expect_true(all(abs(tapply(dy$know.p, dy$group, mean)) < 1e-10))
  expect_true(all(abs(tapply(dy$know.t, dy$group, mean)) < 1e-10))
  expect_true(all(abs(tapply(dy$know.r, dy$group, mean)) < 1e-10))
  # the group shift does not leak into the relationship component
  expect_lt(abs(cor(dy$know.r, shift[match(dy$perceiver, ps$perceiver)])), 1e-8)
  m <- attr(d, "mpa")
  expect_equal(m$vars$group, "know.g")
  f <- mpa_formula(d, criterion = "self", predictors = "know.g", random_slopes = c("know.p", "know.r"))
  expect_false(any(grepl("know.g", unlist(f$random))))
  # across the whole sample on request
  d0 <- mpa_data(judgments = ps, format = "loop", loop_id = "id", targets = st$target_survey, items = st$items,
                 item_cols = "^item_", criterion = list(self = "targets"), group = "group",
                 decompose = "know", decompose_within_group = FALSE, verbose = FALSE)
  expect_false("know.g" %in% names(d0))
  dy0 <- unique(as.data.frame(d0)[, c("perceiver", "target", "know", "know.p", "know.t", "know.r")])
  expect_equal(dy0$know.p + dy0$know.t + dy0$know.r + mean(dy0$know), dy0$know)
})
