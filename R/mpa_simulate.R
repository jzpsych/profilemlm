#' Simulate person judgment data from the multilevel profile analysis model
#'
#' Generates perceiver ratings of targets on items according to the cross-
#' classified profile model: a normative profile (item means), target-specific
#' distinctive deviations, an item property (`des`, e.g. desirability), a
#' dyadic (dyad-level) variable (`know`, e.g. acquaintance), a perceiver
#' variable (`training`) and a target variable (`talkativeness`), fixed effects
#' of the criterion components and their moderation by these variables, random
#' intercepts and slopes for perceivers, targets and dyads, random item
#' intercepts and random item slopes for the distinctive component,
#' optionally shared target x item and perceiver x item deviations (see
#' `sd`), and residual error. By default the generated objects are passed through [mpa_data()], so
#' the result is a regular `mpa_data` object whose criterion source is the targets'
#' (optionally noisy) self-ratings. With `output = "objects"` the six input
#' data sets are returned instead (judgments, criterion, items, perceivers,
#' targets, pairs) -- e.g., to write them to files and start a tutorial from
#' raw data. The generating parameters are stored in the metadata and can be
#' retrieved with [mpa_sim_parameters()] and compared with a fit via
#' [mpa_recovery()].
#'
#' @param n_perceiver,n_target,n_item Numbers of perceivers, targets and items.
#'   For `design = "one-with-many"` the number of perceivers is
#'   `n_target * group_size`; for `"round-robin"` perceivers and targets are
#'   the same `n_target` persons.
#' @param design `"half-block"` (groups of `group_size` perceivers rating the
#'   targets assigned to their group, fully crossed within groups),
#'   `"one-with-many"` (each perceiver rates one target, `group_size`
#'   perceivers per target) or `"round-robin"` (groups of `group_size`
#'   persons, everyone rates everyone else).
#' @param group_size Group size (see `design`).
#' @param targets_per_group For `"half-block"`: targets per group (default
#'   `group_size`). `n_target` is adjusted to `n_group * targets_per_group`.
#' @param midpoint Scale midpoint around which item norms are generated.
#' @param norm_sd SD of the item norms (normative profile) across items.
#' @param dist_sd SD of the target-specific deviations from the norm
#'   (distinctive profiles).
#' @param des_norm_cor Correlation between item desirability and item norm.
#' @param fixed Named numeric vector of fixed effects. Any subset of
#'   `intercept`, `norm`, `dist`; for each moderator `m` in `des`, `know`,
#'   `training`, `talkativeness` the main effect `m` and the moderations `m_norm`
#'   and `m_dist`. Unspecified effects are 0 (except the defaults shown in
#'   Usage). `similarity` is the slope of the perceiver's own distinctive
#'   profile (assumed similarity) and `m_similarity` its moderation by `m`. If
#'   any of them (or a `similarity` entry of `sd`) is nonzero, the perceivers
#'   get self-ratings on the items, which enter the data as the further
#'   predictor profile `pself` (see the `profiles` argument of [mpa_data()]).
#'   In a round-robin design these are the persons' criterion profiles.
#' @param sd Named list of random-effect standard deviations per unit
#'   (`perceiver`, `target`, `dyad`, `item`, `target_item`, `perceiver_item`), each a named vector
#'   with any of `intercept`, `norm`, `dist`, `des`, `similarity`. For items only `intercept`,
#'   `dist` and `similarity` are meaningful (an item has one norm value and one desirability
#'   value, so item slopes for them are not identifiable and are ignored). For
#'   `target_item` only `intercept` is used: a shared deviation of all
#'   perceivers of a target on an item from the criterion profile (consensus,
#'   the component that `random_units = "target:item"` estimates); the default
#'   is 0, and [mpa_sim_parameters()] lists the unit only when it is nonzero.
#'   For `perceiver_item` `intercept` and `dist` are used: a perceiver's
#'   item-specific view of people in general and item-specific differences
#'   between perceivers in distinctive agreement (the components that
#'   `random_units = "perceiver:item"` estimates); default 0, listed only
#'   when nonzero. Random effects are independent.
#' @param resid_sd Residual SD.
#' @param self_error_sd SD of measurement error added to the self-ratings that
#'   serve as criterion source (0 = the true profiles are observed).
#' @param rating_scale Optional `c(min, max)` of the rating scale. If given,
#'   simulated judgments and self-ratings are clipped to this range (as a
#'   bounded rating scale would do) and the value is passed on to [mpa_data()]
#'   (which then also scales the criterion components automatically).
#' @param seed Optional random seed.
#' @param output `"mpa_data"` (default): the analysis data set. `"objects"`: a
#'   list of six tidy input data sets `judgments` (one row per perceiver x
#'   target with the item ratings), `criterion` (one row per target with the
#'   item self-ratings), `items`, `perceivers`, `targets`, `pairs`. `"study"`:
#'   the three files a typical study produces -- `perceiver_survey` (one row
#'   per perceiver, with the judgments and the dyad-level variable looped over
#'   targets in columns `t1_id`, `t1_know`, `t1_item_01`, ..., plus
#'   perceiver-level variables), `target_survey` (one row per target with
#'   self-ratings and target-level variables) and `items`. Both lists carry
#'   the call needed to rebuild the `mpa_data` object in
#'   `attr(x, "mpa_data_call")`. With assumed similarity (see `fixed`) both
#'   lists have the additional data set `perceiver_self` (one row per
#'   perceiver with the item self-ratings).
#' @param write_dir Optional directory; if given, the objects of the chosen
#'   layout (`"objects"` or `"study"`; six or three files) are written as CSV
#'   files.
#' @param ... Passed to [mpa_data()] (e.g., `norm = "sample"`).
#'
#' @return An `mpa_data` object with criterion source `self`, item variable `des`,
#'   dyadic variable `know`, perceiver variable `training`, target variable
#'   `talkativeness` and (for grouped designs) `group`; or the list described
#'   under `output`.
#' @examples
#' dat <- mpa_simulate(n_perceiver = 24, n_item = 12, group_size = 4, seed = 1,
#'                   fixed = c(intercept = 51, norm = 0.7, dist = 0.25, des = 5, know_dist = 0.1),
#'                   rating_scale = c(1, 101))
#' dat
#' f <- mpa_formula(dat, criterion = "self", predictors = c("des", "know"), interactions = "criterion",
#'                  random_units = c("perceiver", "target", "item"),
#'                  correlated = FALSE)
#' m <- mpa_fit(f, dat)
#' mpa_recovery(m)
#' @export
mpa_simulate <- function(n_perceiver = 80, n_target = 80, n_item = 40,
                         design = c("half-block", "one-with-many", "round-robin"),
                         group_size = 4, targets_per_group = group_size,
                         midpoint = 51, norm_sd = 15, dist_sd = 12, des_norm_cor = 0.5,
                         fixed = c(intercept = 51, norm = 0.7, dist = 0.2),
                         sd = list(perceiver = c(intercept = 10, norm = 0.25, dist = 0.10),
                                   target    = c(intercept = 2,  norm = 0.05, dist = 0.08),
                                   dyad      = c(intercept = 6,  norm = 0.20, dist = 0.10),
                                   item      = c(intercept = 3)),
                         resid_sd = 18, self_error_sd = 0, rating_scale = NULL, seed = NULL,
                         output = c("mpa_data", "objects", "study"), write_dir = NULL, ...) {
  design <- match.arg(design)
  output <- match.arg(output)
  if (!is.null(seed)) set.seed(seed)

  mods <- c("des", "know", "training", "talkativeness")
  fx <- stats::setNames(rep(0, 4 + 4 * length(mods)),
                        c("intercept", "norm", "dist", as.vector(rbind(mods, paste0(mods, "_norm"), paste0(mods, "_dist"))),
                          "similarity", paste0(mods, "_similarity")))
  bad <- setdiff(names(fixed), names(fx)); if (length(bad)) stop("Unknown fixed effect(s): ", paste(bad, collapse = ", "))
  fx[names(fixed)] <- fixed
  bad <- setdiff(names(sd), c("perceiver", "target", "dyad", "item", "target_item", "target:item", "perceiver_item", "perceiver:item"))
  if (length(bad)) stop("Unknown unit(s) in `sd`: ", paste(bad, collapse = ", "))
  if (!is.null(sd[["target:item"]])) { sd$target_item <- sd[["target:item"]]; sd[["target:item"]] <- NULL }
  if (!is.null(sd[["perceiver:item"]])) { sd$perceiver_item <- sd[["perceiver:item"]]; sd[["perceiver:item"]] <- NULL }
  sd_all <- list()
  for (u in c("perceiver", "target", "dyad", "item", "target_item", "perceiver_item")) {
    v <- c(intercept = 0, norm = 0, dist = 0, des = 0, similarity = 0)
    if (!is.null(sd[[u]])) { bad <- setdiff(names(sd[[u]]), names(v)); if (length(bad)) stop("Unknown random effect(s) for ", u, ": ", paste(bad, collapse = ", ")); v[names(sd[[u]])] <- sd[[u]] }
    if (u == "item") v[c("norm", "des")] <- 0          # not identifiable within items; see @param sd
    if (u == "target_item") v[c("norm", "dist", "des", "similarity")] <- 0   # a shared target x item view: intercept only
    if (u == "perceiver_item") v[c("norm", "des", "similarity")] <- 0       # a perceiver's item-specific view: intercept and distinctive slope
    sd_all[[u]] <- v
  }
  ## assumed similarity: the perceivers' own profiles are generated only if asked for,
  ## so that data simulated without it do not depend on this option
  sim_names <- c("similarity", paste0(mods, "_similarity"))
  has_sim <- any(fx[sim_names] != 0) || any(vapply(sd_all, function(v) v[["similarity"]] > 0, logical(1)))

  ## ---- design: perceiver x target pairs -----------------------------------
  if (design == "half-block") {
    n_group <- ceiling(n_perceiver / group_size)
    n_target <- n_group * targets_per_group
    perceivers <- sprintf("p%03d", seq_len(n_perceiver)); targets <- sprintf("t%03d", seq_len(n_target))
    p_group <- rep(seq_len(n_group), each = group_size)[seq_len(n_perceiver)]
    t_group <- rep(seq_len(n_group), each = targets_per_group)
    pairs <- do.call(rbind, lapply(seq_len(n_group), function(g)
      expand.grid(perceiver = perceivers[p_group == g], target = targets[t_group == g], stringsAsFactors = FALSE)))
    pairs$group <- p_group[match(pairs$perceiver, perceivers)]
  } else if (design == "one-with-many") {
    n_perceiver <- n_target * group_size
    perceivers <- sprintf("p%03d", seq_len(n_perceiver)); targets <- sprintf("t%03d", seq_len(n_target))
    pairs <- data.frame(perceiver = perceivers, target = rep(targets, each = group_size), stringsAsFactors = FALSE)
    pairs$group <- as.integer(factor(pairs$target, levels = targets))
  } else {
    n_group <- ceiling(n_target / group_size)
    n_target <- n_group * group_size; n_perceiver <- n_target
    persons <- sprintf("s%03d", seq_len(n_target)); perceivers <- targets <- persons
    grp <- rep(seq_len(n_group), each = group_size)
    pairs <- do.call(rbind, lapply(seq_len(n_group), function(g)
      expand.grid(perceiver = persons[grp == g], target = persons[grp == g], stringsAsFactors = FALSE)))
    pairs <- pairs[pairs$perceiver != pairs$target, ]
    pairs$group <- grp[match(pairs$perceiver, persons)]
  }
  rownames(pairs) <- NULL
  items <- sprintf("item_%02d", seq_len(n_item))

  ## ---- items: norm and desirability ---------------------------------------
  z1 <- stats::rnorm(n_item); z2 <- stats::rnorm(n_item)
  norm_i <- midpoint + norm_sd * z1
  des_i  <- 2 * stats::pnorm(des_norm_cor * z1 + sqrt(1 - des_norm_cor^2) * z2) - 1     # in (-1, 1)
  normc_i <- norm_i - mean(norm_i)

  ## ---- targets: true profiles ---------------------------------------------
  dist_ti <- matrix(stats::rnorm(n_target * n_item, 0, dist_sd), n_target, n_item)
  criterion_ti <- sweep(dist_ti, 2, norm_i, "+")

  ## ---- perceivers: own profiles (assumed similarity) ----------------------
  sdist_pi <- NULL
  if (has_sim) {
    sdist_pi <- if (design == "round-robin") dist_ti else matrix(stats::rnorm(n_perceiver * n_item, 0, dist_sd), n_perceiver, n_item)
  }

  ## ---- moderators: acquaintance (pair), training (perceiver), talkativeness (target)
  pairs$know <- 2 * stats::pnorm(stats::rnorm(nrow(pairs))) - 1
  training_p  <- 2 * stats::pnorm(stats::rnorm(n_perceiver)) - 1
  talkativeness_t <- 2 * stats::pnorm(stats::rnorm(n_target)) - 1

  ## ---- random effects -------------------------------------------------------
  re <- function(n, sdv) sapply(sdv, function(s) if (s > 0) stats::rnorm(n, 0, s) else rep(0, n))
  re_p <- re(n_perceiver, sd_all$perceiver); rownames(re_p) <- perceivers
  re_t <- re(n_target, sd_all$target);       rownames(re_t) <- targets
  re_d <- re(nrow(pairs), sd_all$dyad)
  re_i <- re(n_item, sd_all$item)
  re_ti <- matrix(if (sd_all$target_item[["intercept"]] > 0) stats::rnorm(n_target * n_item, 0, sd_all$target_item[["intercept"]]) else 0,
                  n_target, n_item)                 # shared target x item deviation (consensus beyond the criterion)
  re_pi <- function(s) matrix(if (s > 0) stats::rnorm(n_perceiver * n_item, 0, s) else 0, n_perceiver, n_item)
  re_pi_int <- re_pi(sd_all$perceiver_item[["intercept"]])   # perceiver x item: item-specific view of people in general
  re_pi_dist <- re_pi(sd_all$perceiver_item[["dist"]])      # perceiver x item: item-specific distinctive agreement

  ## ---- ratings ---------------------------------------------------------------
  ip <- match(pairs$perceiver, perceivers); it <- match(pairs$target, targets)
  y <- matrix(NA_real_, nrow(pairs), n_item, dimnames = list(NULL, items))
  for (i in seq_len(n_item)) {
    b0 <- fx[["intercept"]] + re_p[ip, "intercept"] + re_t[it, "intercept"] + re_d[, "intercept"] + re_i[i, "intercept"] + re_ti[it, i] + re_pi_int[ip, i]
    b1 <- fx[["norm"]] + re_p[ip, "norm"] + re_t[it, "norm"] + re_d[, "norm"]
    b2 <- fx[["dist"]] + re_p[ip, "dist"] + re_t[it, "dist"] + re_d[, "dist"] + re_i[i, "dist"] + re_pi_dist[ip, i]
    b3 <- fx[["des"]]  + re_p[ip, "des"]  + re_t[it, "des"]  + re_d[, "des"]
    D  <- dist_ti[it, i]
    ex <- training_p[ip]; vi <- talkativeness_t[it]; kn <- pairs$know
    y[, i] <- b0 + b1 * normc_i[i] + b2 * D + b3 * des_i[i] +
      fx[["des_norm"]] * des_i[i] * normc_i[i] + fx[["des_dist"]] * des_i[i] * D +
      fx[["know"]] * kn + fx[["know_norm"]] * kn * normc_i[i] + fx[["know_dist"]] * kn * D +
      fx[["training"]] * ex + fx[["training_norm"]] * ex * normc_i[i] + fx[["training_dist"]] * ex * D +
      fx[["talkativeness"]] * vi + fx[["talkativeness_norm"]] * vi * normc_i[i] + fx[["talkativeness_dist"]] * vi * D +
      stats::rnorm(nrow(pairs), 0, resid_sd)
    if (has_sim) {
      b4 <- fx[["similarity"]] + re_p[ip, "similarity"] + re_t[it, "similarity"] + re_d[, "similarity"] + re_i[i, "similarity"] +
        fx[["des_similarity"]] * des_i[i] + fx[["know_similarity"]] * kn + fx[["training_similarity"]] * ex +
        fx[["talkativeness_similarity"]] * vi
      y[, i] <- y[, i] + b4 * sdist_pi[ip, i]
    }
  }

  ## ---- the six input data sets ----------------------------------------------
  self <- criterion_ti + if (self_error_sd > 0) stats::rnorm(length(criterion_ti), 0, self_error_sd) else 0
  if (!is.null(rating_scale)) {
    y    <- pmin(pmax(y, rating_scale[1]), rating_scale[2])
    self <- pmin(pmax(self, rating_scale[1]), rating_scale[2])
  }
  colnames(self) <- items
  pself_df <- NULL
  if (has_sim) {
    pself <- if (design == "round-robin") self else sweep(sdist_pi, 2, norm_i, "+") +
      if (self_error_sd > 0) stats::rnorm(length(sdist_pi), 0, self_error_sd) else 0
    if (!is.null(rating_scale)) pself <- pmin(pmax(pself, rating_scale[1]), rating_scale[2])
    colnames(pself) <- items
    pself_df <- cbind(data.frame(perceiver = perceivers, stringsAsFactors = FALSE), as.data.frame(pself))
  }
  judgments  <- cbind(pairs[, c("perceiver", "target", "group")], as.data.frame(y))
  pairs_df   <- pairs[, c("perceiver", "target", "know")]
  criterion_df   <- cbind(data.frame(target = targets, stringsAsFactors = FALSE), as.data.frame(self))
  items_df   <- data.frame(item = items, des = des_i, stringsAsFactors = FALSE)
  perc_df    <- data.frame(perceiver = perceivers, training = training_p, stringsAsFactors = FALSE)
  targ_df    <- data.frame(target = targets, talkativeness = talkativeness_t, stringsAsFactors = FALSE)
  sim <- list(design = design, fixed = fx, sd = sd_all, resid_sd = resid_sd,
              norm_sd = norm_sd, dist_sd = dist_sd, des_norm_cor = des_norm_cor,
              self_error_sd = self_error_sd, rating_scale = rating_scale, seed = seed, similarity = has_sim,
              item_norm = stats::setNames(norm_i, items), item_des = stats::setNames(des_i, items))

  ## "study" layout: one row per perceiver, targets looped in columns
  grp_p <- pairs$group[match(perceivers, pairs$perceiver)]
  idx_by_p <- split(seq_len(nrow(pairs)), factor(pairs$perceiver, levels = perceivers))
  k <- max(lengths(idx_by_p))
  perceiver_survey <- data.frame(perceiver = perceivers, group = grp_p, training = training_p, stringsAsFactors = FALSE)
  for (sl in seq_len(k)) {
    rows <- vapply(idx_by_p, function(ix) if (length(ix) >= sl) ix[sl] else NA_integer_, integer(1))
    perceiver_survey[[paste0("t", sl, "_id")]]   <- pairs$target[rows]
    perceiver_survey[[paste0("t", sl, "_know")]] <- pairs$know[rows]
    for (i in seq_len(n_item)) perceiver_survey[[paste0("t", sl, "_", items[i])]] <- y[rows, i]
  }
  target_survey <- cbind(data.frame(target = targets, talkativeness = talkativeness_t, stringsAsFactors = FALSE), as.data.frame(self))

  if (!is.null(write_dir)) {
    if (!dir.exists(write_dir)) dir.create(write_dir, recursive = TRUE)
    files <- if (output == "study") list(perceiver_survey = perceiver_survey, target_survey = target_survey, items = items_df)
             else list(judgments = judgments, criterion = criterion_df, items = items_df,
                       perceivers = perc_df, targets = targ_df, dyads = pairs_df)
    for (nm in names(files)) utils::write.csv(files[[nm]], file.path(write_dir, paste0(nm, ".csv")), row.names = FALSE)
    message(length(files), " data sets written to ", normalizePath(write_dir))
  }
  if (output == "study") {
    out <- list(perceiver_survey = perceiver_survey, target_survey = target_survey, items = items_df)
    if (has_sim) out$perceiver_self <- pself_df
    attr(out, "simulation") <- sim
    attr(out, "mpa_data_call") <- if (has_sim) bquote(
      mpa_data(judgments = perceiver_survey, format = "loop", loop_id = "id", targets = target_survey, items = items,
               item_cols = "^item_", criterion = list(self = "targets"), profiles = list(pself = perceiver_self),
               group = "group", rating_scale = .(rating_scale))) else bquote(
      mpa_data(judgments = perceiver_survey, format = "loop", loop_id = "id", targets = target_survey, items = items,
               item_cols = "^item_", criterion = list(self = "targets"), group = "group", rating_scale = .(rating_scale)))
    return(out)
  }
  if (output == "objects") {
    out <- list(judgments = judgments, criterion = criterion_df, items = items_df,
                perceivers = perc_df, targets = targ_df, dyads = pairs_df)
    if (has_sim) out$perceiver_self <- pself_df
    attr(out, "simulation") <- sim
    attr(out, "mpa_data_call") <- if (has_sim) bquote(
      mpa_data(judgments = judgments, dyads = dyads, targets = targets, perceivers = perceivers, items = items,
               item_cols = "^item_", criterion = list(self = criterion), profiles = list(pself = perceiver_self),
               group = "group", rating_scale = .(rating_scale))) else bquote(
      mpa_data(judgments = judgments, dyads = dyads, targets = targets, perceivers = perceivers, items = items,
               item_cols = "^item_", criterion = list(self = criterion), group = "group", rating_scale = .(rating_scale)))
    return(out)
  }

  d <- mpa_data(judgments = judgments, dyads = pairs_df, targets = targ_df, perceivers = perc_df, items = items_df,
                item_cols = items, criterion = list(self = criterion_df),
                profiles = if (has_sim) list(pself = pself_df) else NULL, group = "group",
                rating_scale = rating_scale, verbose = FALSE, ...)
  m <- attr(d, "mpa")
  sim$criterion_scale <- m$settings$criterion_scale
  m$simulation <- sim
  attr(d, "mpa") <- m
  d
}

#' Attach simulation parameters to an mpa_data object built from simulated files
#'
#' When the six data sets from `mpa_simulate(output = "objects")` (or the CSV
#' files written by `write_dir`) are read back and combined with [mpa_data()],
#' the generating parameters are no longer attached. `mpa_sim_attach()` adds
#' them so that [mpa_recovery()] works.
#'
#' @param data An `mpa_data` object.
#' @param objects The list returned by `mpa_simulate(output = "objects")`.
#' @return `data` with the simulation parameters in its metadata.
#' @export
mpa_sim_attach <- function(data, objects) {
  sim <- attr(objects, "simulation")
  if (is.null(sim)) stop("`objects` carries no simulation parameters.")
  m <- attr(data, "mpa")
  sim$criterion_scale <- m$settings$criterion_scale
  m$simulation <- sim
  attr(data, "mpa") <- m
  data
}

#' Generating parameters of a simulated data set
#'
#' @param data An `mpa_data` object created by [mpa_simulate()].
#' @param scale `"criterion"` (default): parameters per criterion point, matching
#'   `mpa_fit(..., report = "criterion")`; `"analysis"`: parameters on the scale of
#'   the (possibly divided) analysis variables.
#' @return A list with `fixed` (data frame: `term`, `true`, with `lme4`-style
#'   term names) and `variance` (data frame: `unit`, `effect`, `true_sd`).
#' @export
mpa_sim_parameters <- function(data, scale = c("criterion", "analysis")) {
  scale <- match.arg(scale)
  sim <- attr(data, "mpa")$simulation
  if (is.null(sim)) stop("`data` was not created by mpa_simulate().")
  fx <- sim$fixed
  ts <- if (scale == "analysis" && is.numeric(sim$criterion_scale)) sim$criterion_scale else 1
  mods <- c("des", "know", "training", "talkativeness")
  fixed <- data.frame(
    term = c("(Intercept)", "self.norm.c", "self.dist",
             as.vector(rbind(mods, paste0("self.norm.c:", mods), paste0("self.dist:", mods)))),
    true = unname(fx[c("intercept", "norm", "dist", as.vector(rbind(mods, paste0(mods, "_norm"), paste0(mods, "_dist"))))]) *
      c(1, ts, ts, rep(c(1, ts, ts), length(mods))),
    stringsAsFactors = FALSE)
  has_sim <- isTRUE(sim$similarity)
  if (has_sim) fixed <- rbind(fixed, data.frame(
    term = c("pself.dist", paste0("pself.dist:", mods)),
    true = unname(fx[c("similarity", paste0(mods, "_similarity"))]) * ts, stringsAsFactors = FALSE))
  eff_names <- c(intercept = "(Intercept)", norm = "self.norm.c", dist = "self.dist", des = "des", similarity = "pself.dist")
  unit_names <- c(perceiver = "perceiver", target = "target", dyad = "perceiver:target", item = "item", target_item = "target:item", perceiver_item = "perceiver:item")
  variance <- do.call(rbind, lapply(names(sim$sd), function(u) {
    v <- sim$sd[[u]]
    v <- v[setdiff(names(v), if (!has_sim) "similarity")]
    if (u == "item") v <- v[c("intercept", "dist", if (has_sim) "similarity")]
    if (u == "target_item") { v <- v["intercept"]; if (v[["intercept"]] == 0) return(NULL) }   # not a default unit: listed only when present
    if (u == "perceiver_item") { v <- v[c("intercept", "dist")]; if (all(v == 0)) return(NULL) }
    v[intersect(c("norm", "dist", "similarity"), names(v))] <- v[intersect(c("norm", "dist", "similarity"), names(v))] * ts
    data.frame(unit = unit_names[[u]], effect = eff_names[names(v)], true_sd = unname(v), stringsAsFactors = FALSE)
  }))
  variance <- rbind(variance, data.frame(unit = "Residual", effect = "", true_sd = sim$resid_sd))
  rownames(variance) <- NULL
  list(fixed = fixed, variance = variance)
}

#' Compare a fitted model with the generating parameters of simulated data
#'
#' @param fit An `mpa_fit` fitted to data from [mpa_simulate()].
#' @param data Optional: the simulated `mpa_data`. By default the metadata
#'   stored in `fit` are used.
#' @return A list of two data frames: `fixed` (true value, estimate, CI and
#'   whether the CI covers the true value; only terms present in the model) and
#'   `variance` (true and estimated SDs).
#' @export
mpa_recovery <- function(fit, data = NULL) {
  stopifnot(inherits(fit, "mpa_fit"))
  sim <- if (!is.null(data)) attr(data, "mpa")$simulation else fit$meta$simulation
  if (is.null(sim)) stop("Simulation parameters not found; pass the `mpa_data` object created by mpa_simulate().")
  par <- mpa_sim_parameters(structure(list(), mpa = list(simulation = sim)),
                            scale = if (isTRUE(fit$report == "analysis")) "analysis" else "criterion")
  key <- function(x) vapply(strsplit(x, ":", fixed = TRUE), function(z) paste(sort(z), collapse = ":"), character(1))

  fx <- fit$fixed
  fx$key <- key(fx$term); par$fixed$key <- key(par$fixed$term)
  fixed <- merge(fx, par$fixed[, c("key", "true")], by = "key", all.x = TRUE)
  fixed$true[is.na(fixed$true)] <- 0
  fixed$covered <- fixed$true >= fixed$ci_lower & fixed$true <= fixed$ci_upper
  fixed <- fixed[match(fx$key, fixed$key), c("term", "true", "estimate", "se", "ci_lower", "ci_upper", "covered")]
  rownames(fixed) <- NULL

  vt <- fit$variance
  variance <- merge(vt, par$variance, by = c("unit", "effect"), all.x = TRUE)
  variance$true_sd[is.na(variance$true_sd)] <- 0
  variance <- variance[order(match(variance$unit, unique(vt$unit)), variance$effect != "(Intercept)"),
                       c("unit", "effect", "true_sd", "sd", "variance")]
  names(variance)[4] <- "estimated_sd"
  rownames(variance) <- NULL
  structure(list(fixed = fixed, variance = variance), class = "mpa_recovery")
}

#' @export
print.mpa_recovery <- function(x, digits = 3, ...) {
  cat("<mpa_recovery>  fixed effects: true vs. estimated (95% CI)\n")
  fx <- x$fixed
  fx[c("true", "estimate", "se", "ci_lower", "ci_upper")] <- lapply(fx[c("true", "estimate", "se", "ci_lower", "ci_upper")], round, digits)
  print(fx, row.names = FALSE)
  cat("\nRandom-effect SDs: true vs. estimated\n")
  vt <- x$variance
  vt$true_sd <- round(vt$true_sd, digits); vt$estimated_sd <- round(vt$estimated_sd, digits); vt$variance <- round(vt$variance, digits)
  print(vt, row.names = FALSE)
  invisible(x)
}
