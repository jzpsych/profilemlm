#' Build the canonical long-format data set for multilevel profile analysis
#'
#' `mpa_data()` takes the judgments plus optional predictor files -- one row
#' per dyad, per target, per perceiver and per item -- and returns a single long-format
#' data frame (one row per perceiver x target x item) that all other functions
#' of the package operate on. It derives the judgment profile (the dependent variable), one or more
#' "criterion" (predictor) profiles, their normative and distinctive components,
#' merges item, perceiver, target and dyad variables, optionally rescales
#' variables, and detects the design type (one-with-many, half-block,
#' round-robin).
#'
#' @param judgments The judgment data: by default one row per rater x target
#'   (dyad) with the perceiver and target identifier columns named in `id` and
#'   the item rating columns given by `item_cols` (see `format` for other
#'   layouts). All remaining columns are treated as predictors; their level
#'   (perceiver, target or dyad) is detected (see `detect_levels`). If `role`
#'   is given, the data frame may also contain rows of other rater types (e.g.,
#'   self-ratings, target-nominated informants) that serve as criterion sources.
#' @param dyads Optional data frame with one row per perceiver x target dyad
#'   holding dyadic variables such as acquaintance or liking, if they are not
#'   already columns of `judgments`. Merged by the perceiver and target
#'   identifiers.
#' @param targets Optional data frame with one row per target: the target
#'   identifier plus target-level variables. May also contain item columns
#'   (same names as `item_cols`) holding self-ratings; refer to them with
#'   `criterion = list(self = "targets")`.
#' @param perceivers Optional data frame with one row per perceiver: the
#'   perceiver identifier plus perceiver-level variables. May also contain
#'   item columns (same names as `item_cols`) holding the perceivers'
#'   self-ratings; refer to them with `profiles = list(pself = "perceivers")`.
#' @param items Optional data frame with one row per item: item-level variables
#'   (item characteristics). If it contains the item identifier column named in
#'   `id`, its values must match `item_cols`; otherwise rows are matched to
#'   `item_cols` by position.
#' @param id Named character vector giving the identifier column names for
#'   `perceiver`, `target` and `item`. Defaults to `c(perceiver = "perceiver",
#'   target = "target", item = "item")`. Partial specification is allowed.
#' @param item_cols Character vector of the item rating columns in `dyads`
#'   (and, if used, `targets`), or a single regular expression such as
#'   `"^item_"` that selects them. The order defines the item order.
#' @param format Layout of `judgments`. `"wide"` (default): one row per rater x
#'   target with one column per item. `"long"`: one row per rater x target x
#'   item, with the item identifier in the column named in `id` and the rating
#'   in `rating_col`. `"loop"`: one row per perceiver, as exported by survey
#'   software that loops the same block of questions over several targets; the
#'   looped columns are identified by `loop_pattern` and the target of each loop
#'   by `loop_id`. Long and loop data are cast to wide internally.
#' @param rating_col Name of the rating column when `format = "long"`.
#' @param loop_pattern For `format = "loop"`: a regular expression with two
#'   capture groups, the first matching the loop index and the second the
#'   variable name, e.g. the default `"^t(\\d+)_(.+)$"` for columns
#'   `t1_item_01`, `t1_know`, `t2_item_01`, ... Columns not matching the pattern
#'   are perceiver-level (repeated for every target).
#' @param loop_id For `format = "loop"`: the looped variable that holds the
#'   target identifier (e.g. `"id"` for columns `t1_id`, `t2_id`, ...). If
#'   `NULL`, the loop index itself is used as target identifier (all
#'   perceivers rated the same targets in the same order).
#' @param detect_levels If `TRUE` (default), variables in `judgments` (or `dyads`)
#'   that are constant within perceivers are registered as perceiver-level
#'   variables and variables that are constant within targets as target-level
#'   variables; only variables that vary within both are treated as dyadic. This matters for the identifiability of random slopes in
#'   [mpa_formula()]. In a one-with-many design every dyad-level variable is
#'   constant within perceivers and is therefore registered as perceiver-level
#'   (the two are equivalent there).
#' @param role Optional name of a column in `dyads` that classifies rows into
#'   rater types (e.g., `"perceiver"`, `"self"`, `"TNI"`).
#' @param judgment_role Values of `role` identifying the rows that provide the
#'   judgments. If `NULL`, all rows whose role is not used by a criterion
#'   source are treated as judgment rows.
#' @param criterion Named list of criterion (predictor) sources. Each element is one of
#'   * a character string naming a value of `role` (rows in `dyads` with that
#'     role are aggregated per target and item),
#'   * the keyword `"targets"` (item columns in `targets`),
#'   * a list `list(role = "TNI", max_n = NULL, fun = mean)` where `max_n` sets
#'     the maximum number of raters per target (targets with more raters get
#'     `NA`, e.g., to remove ambiguous multiple self-ratings),
#'   * a data frame, either wide (target id + `item_cols`) or long (target id,
#'     item id, `value`).
#' @param criterion_combine If several criterion sources are given and
#'   `criterion_combine = "mean"`, an additional source `combined` (row mean over
#'   sources, `NA` ignored) is created.
#' @param profiles Named list of further *predictor profiles*: profiles on
#'   the same items as the judgments that enter a model next to the criterion
#'   (see the `profiles` argument of [mpa_formula()]), for example the
#'   perceivers' self-ratings, whose slope is assumed similarity beyond actual
#'   similarity. Each element
#'   is one of
#'   * the keyword `"perceivers"` (item columns in `perceivers`): a profile
#'     that belongs to the perceiver,
#'   * the keyword `"targets"` (item columns in `targets`) or a value of
#'     `role`: a profile that belongs to the target,
#'   * a data frame, wide (identifier + `item_cols`) or long (identifier, item
#'     id, `value`). It belongs to the perceiver if it has the perceiver
#'     identifier column and not the target identifier column, and to the
#'     target otherwise.
#'
#'   For each profile `<p>` the columns `<p>`, `<p>.norm` and `<p>.dist` are
#'   created. The normative profile is the mean of that profile across its
#'   perceivers (or targets), computed as `norm` says (an external `norm`
#'   applies to profiles too), and `<p>.dist` is the deviation from it, scaled
#'   like the distinctive component of the criterion. Profiles are not part of `criterion_combine`. A criterion
#'   source that is not the criterion of a model can serve as a further
#'   profile too, without being named here.
#' @param norm How the normative profile of each criterion source is computed:
#'   `"loo"` (leave-one-out: for each target the mean over all *other* targets;
#'   default), `"sample"` (mean over all targets), or a data frame with an
#'   `item` column and one column per criterion source (or a single column `norm`
#'   applied to all sources) holding an external normative profile.
#' @param norm_by Optional name of a target-level variable in `targets`; if
#'   given, normative profiles are computed within its levels (e.g., separate
#'   norms by gender).
#' @param center_norm How the normative profile is centered to obtain
#'   `<source>.norm.c`: `"grand_mean"` (the mean across items of the sample
#'   normative profile, per source and `norm_by` level; default -- with
#'   `norm = "loo"` the same constant is subtracted from every target's
#'   leave-one-out profile, so that the centering does not vary between
#'   targets), `"midpoint"` (midpoint of `rating_scale`), or `"none"`.
#' @param rating_scale Optional numeric vector `c(min, max)` of the rating
#'   scale of the item columns; required for `center_norm = "midpoint"` and
#'   stored as metadata.
#' @param rename Optional named character vector `c(new = "old")` renaming
#'   variables after merging (e.g., `c(DES = "mean_Pos")`).
#' @param rescale Optional named list `list(var = c(min, max))`. For each
#'   entry a column `<var>.n` is added that maps the theoretical scale range to
#'   -1..1 (scale midpoint = 0). Applied after `rename`.
#' @param derive Optional named list of derived variables, computed after
#'   `rename` and `rescale`. Each element is a one-sided formula or a function
#'   of the data, e.g. `list(EVAL.n = ~ abs(DES.n))` for evaluativeness as the
#'   absolute desirability, or `list(KNOW.c = function(d) d$KNOW.n - mean(d$KNOW.n))`.
#'   A derived variable inherits the level (item, dyad, perceiver, target) of
#'   the variables it is computed from, provided they share one level.
#' @param decompose Optional character vector of dyadic variables to
#'   decompose into a perceiver component (`<var>.p`: the perceiver's mean over
#'   its dyads, centered), a target component (`<var>.t`) and a relationship
#'   component (`<var>.r = value - perceiver mean - target mean + grand mean`;
#'   Social Relations Model terminology). Components are computed on the dyad
#'   level (one value per dyad) and registered as perceiver-, target- and
#'   dyad-level variables, so that [mpa_formula()] gives them random slopes
#'   only where identifiable. Applied after `rename` and `rescale`; the
#'   decomposition requires perceivers and targets to be crossed (half-block,
#'   round-robin). The components are simple means; in round-robin data the
#'   Social Relations Model's corrections for the missing self-dyad are not
#'   applied.
#' @param decompose_within_group If `TRUE` (default) and the data have a
#'   `group` column, the decomposition is carried out within groups: perceiver
#'   and target means are taken relative to their group's mean, and an
#'   additional group component `<var>.g` (group mean minus grand mean) is
#'   created and registered as a group-level variable (no random slopes).
#'   In a half-block design perceivers and targets are crossed only within
#'   groups, so without this a group difference would enter the perceiver and
#'   target components twice and the relationship component with the opposite
#'   sign. `FALSE` decomposes across the whole sample.
#' @param criterion_scale Scaling of the *analysis* variables `<src>.norm.c` and
#'   `<src>.dist` (the descriptive variables `<src>` and `<src>.norm` remain on
#'   the rating scale). `"auto"` (default): half the range of `rating_scale`,
#'   i.e., the same transformation that `rescale` applies to the other
#'   predictors (a rating scale 1-101 gives 50, so a criterion difference of one
#'   scale-half equals 1); falls back to 1 with a message if `rating_scale` is
#'   not given. `"sd"`: the pooled SD of the distinctive components (a
#'   data-dependent standardisation). A positive number divides by that
#'   number. Scaling puts all predictors on comparable ranges, which avoids
#'   optimiser warnings and speeds up estimation; it does not change the model.
#'   [mpa_fit()] reports results per criterion point by default (`report = "criterion"`),
#'   so the choice is invisible in the output.
#' @param group Optional name of a column in `dyads` or `perceivers` holding a
#'   group identifier (e.g., subgroups within which perceivers and targets are
#'   crossed).
#' @param na_action `"keep"` (default) or `"drop"` rows with a missing
#'   judgment.
#' @param verbose Print messages about data processing decisions.
#'
#' @return A data frame of class `mpa_data` with the reserved columns
#'   `perceiver`, `target`, `item`, `item_index`, `y` (judgment),
#'   optionally `group`, and for each criterion source `<src>`, `<src>.norm`,
#'   `<src>.norm.c`, `<src>.dist` and `<src>.sd`, followed by all item, dyad,
#'   perceiver and target variables. `<src>.sd` is the standard deviation of the
#'   criterion across targets, one value per item (and `norm_by` level): the
#'   *item variance* of the criterion profile, which de Vries et al. (2016) found
#'   to be the strongest item-level correlate of self-other agreement. It is
#'   registered as an item-level variable, so it can enter `mpa_formula()` as a
#'   predictor or moderator, for instance in a sensitivity analysis that asks
#'   whether an item property moderates agreement beyond the spread of the
#'   criterion on that item. Because it is computed from the same criterion
#'   values whose agreement the model estimates, it is a descriptive covariate,
#'   not a cause. Metadata (variable roles, settings, design, criterion rater
#'   counts) is stored in `attr(x, "mpa")` and shown by `print()` and
#'   `summary()`, which also reports the item-level correlations of the
#'   normative profile and of `<src>.sd` with the item variables.
#' @references de Vries, R. E., Realo, A., & Allik, J. (2016). Using personality
#'   item characteristics to predict single-item internal reliability, retest
#'   reliability, and self-other agreement. *European Journal of Personality,
#'   30*(6), 618-636. \doi{10.1002/per.2083}
#'
#' @examples
#' ## the three files of the tutorial study (perceiver survey in loop format,
#' ## target survey with self-ratings, item properties)
#' st <- mpa_simulate(n_perceiver = 16, n_item = 10, group_size = 4, seed = 1,
#'                    rating_scale = c(1, 101), output = "study")
#' dat <- mpa_data(judgments = st$perceiver_survey, format = "loop", loop_id = "id",
#'               targets = st$target_survey, items = st$items, item_cols = "^item_",
#'               criterion = list(self = "targets"), group = "group",
#'               rating_scale = c(1, 101))
#' dat
#' summary(dat)
#'
#' ## one row per informant x item, self-ratings repeated in every row
#' ## (Wessels et al., 2020, Study 2; see ?wessels_judgments)
#' self_long <- unique(wessels_judgments[, c("target", "item", "self")])
#' names(self_long)[3] <- "value"
#' dw <- mpa_data(judgments = wessels_judgments[, c("perceiver", "target", "item",
#'                                                   "rating", "knowing", "liking", "type")],
#'                format = "long", rating_col = "rating", items = wessels_items,
#'                criterion = list(self = self_long), rating_scale = c(1, 5),
#'                rescale = list(DES = c(1, 10), OBS = c(1, 10)))
#' dw
#' @export
mpa_data <- function(judgments,
                     dyads = NULL,
                     targets = NULL,
                     perceivers = NULL,
                     items = NULL,
                     id = c(perceiver = "perceiver", target = "target", item = "item"),
                     item_cols = NULL,
                     format = c("wide", "long", "loop"),
                     rating_col = "rating",
                     loop_pattern = "^t(\\d+)_(.+)$",
                     loop_id = NULL,
                     detect_levels = TRUE,
                     role = NULL,
                     judgment_role = NULL,
                     criterion = NULL,
                     criterion_combine = c("mean", "none"),
                     profiles = NULL,
                     norm = c("loo", "sample"),
                     norm_by = NULL,
                     center_norm = c("grand_mean", "midpoint", "none"),
                     rating_scale = NULL,
                     rename = NULL,
                     rescale = NULL,
                     derive = NULL,
                     decompose = NULL,
                     decompose_within_group = TRUE,
                     criterion_scale = "auto",
                     group = NULL,
                     na_action = c("keep", "drop"),
                     verbose = TRUE) {

  criterion_combine <- match.arg(criterion_combine)
  center_norm   <- match.arg(center_norm)
  na_action     <- match.arg(na_action)
  msg <- function(...) if (verbose) message(...)

  ## ---- identifiers -------------------------------------------------------
  id  <- .mpa_complete_id(id)
  pid <- id[["perceiver"]]; tid <- id[["target"]]; iid <- id[["item"]]

  if (!is.data.frame(judgments)) stop("`judgments` must be a data frame.")
  pairs <- dyads                                  # dyad-level variables (optional)
  dyads <- as.data.frame(judgments, stringsAsFactors = FALSE)   # internal name for the judgment rows
  if (!pid %in% names(dyads)) stop("Column `", pid, "` not found in `judgments`. Adjust the `id` argument.")
  if (!is.null(targets)) {
    targets <- as.data.frame(targets, stringsAsFactors = FALSE)
    if (!tid %in% names(targets)) stop("Column `", tid, "` not found in `targets`.")
    targets[[tid]] <- as.character(targets[[tid]])
    if (anyDuplicated(targets[[tid]])) stop("`targets` must have exactly one row per target.")
  }
  if (!is.null(perceivers)) {
    perceivers <- as.data.frame(perceivers, stringsAsFactors = FALSE)
    if (!pid %in% names(perceivers)) stop("Column `", pid, "` not found in `perceivers`.")
    perceivers[[pid]] <- as.character(perceivers[[pid]])
    if (anyDuplicated(perceivers[[pid]])) stop("`perceivers` must have exactly one row per perceiver.")
  }

  ## ---- long / loop -> wide ----------------------------------------------
  format <- match.arg(format)
  if (format == "long") {
    cast <- .mpa_long_to_wide(dyads, pid, tid, iid, rating_col, item_cols, items, msg)
    dyads <- cast$wide
    item_cols <- cast$item_cols
  } else if (format == "loop") {
    dyads <- .mpa_loop_to_wide(dyads, pid, tid, loop_pattern, loop_id, msg)
  }
  for (v in c(pid, tid)) {
    if (!v %in% names(dyads)) stop("Column `", v, "` not found in `judgments`. Adjust the `id` argument.")
    dyads[[v]] <- as.character(dyads[[v]])
  }

  ## ---- dyad-level variables from a separate object ------------------------
  if (!is.null(pairs)) {
    pairs <- as.data.frame(pairs, stringsAsFactors = FALSE)
    for (v in c(pid, tid)) {
      if (!v %in% names(pairs)) stop("Column `", v, "` not found in `dyads`.")
      pairs[[v]] <- as.character(pairs[[v]])
    }
    if (anyDuplicated(pairs[, c(pid, tid)])) stop("`dyads` must have exactly one row per perceiver x target dyad.")
    clash <- intersect(setdiff(names(pairs), c(pid, tid)), names(dyads))
    if (length(clash)) stop("Column(s) in both `judgments` and `dyads`: ", paste(clash, collapse = ", "))
    dyads <- merge(dyads, pairs, by = c(pid, tid), all.x = TRUE, sort = FALSE)
  }

  ## ---- item columns ------------------------------------------------------
  item_cols <- .mpa_resolve_item_cols(item_cols, dyads, items, iid)
  n_items <- length(item_cols)

  reserved <- c("perceiver", "target", "item", "item_index", "y", "group")
  clash <- setdiff(intersect(reserved, names(dyads)), c(pid, tid, group))
  if (length(clash))
    stop("`judgments` contains reserved column name(s): ", paste(clash, collapse = ", "),
         ". Please rename them before calling mpa_data().")

  ## ---- roles, criterion specification ---------------------------------------
  if (!is.null(role)) {
    if (!role %in% names(dyads)) stop("Column `", role, "` not found in `judgments`.")
    roles <- as.character(dyads[[role]])
  } else {
    roles <- rep("judgment", nrow(dyads))
  }
  criterion_spec  <- .mpa_parse_criterion(criterion)
  profile_spec <- .mpa_parse_profiles(profiles, pid, tid)
  criterion_roles <- unlist(lapply(c(criterion_spec, profile_spec), function(s) if (s$type == "role") s$role else NULL))
  if (length(criterion_roles) && is.null(role))
    stop("Criterion sources refer to rater roles, but `role` is not specified.")
  if (length(criterion_roles) && !all(criterion_roles %in% roles))
    stop("Role value(s) not found in `judgments$", role, "`: ",
         paste(setdiff(criterion_roles, roles), collapse = ", "))

  if (is.null(judgment_role)) {
    crit_idx <- !(roles %in% criterion_roles)
  } else {
    if (!all(judgment_role %in% roles))
      stop("`judgment_role` value(s) not found in `judgments$", role, "`.")
    crit_idx <- roles %in% judgment_role
  }
  if (!any(crit_idx)) stop("No judgment rows left in `judgments`.")

  ## ---- judgments: wide -> long ------------------------------------------
  crit <- dyads[crit_idx, , drop = FALSE]
  dup  <- duplicated(crit[, c(pid, tid)])
  if (any(dup))
    warning(sum(dup), " duplicated perceiver-target pair(s) among the judgment rows. ",
            "Each dyad should appear only once; check `role`/`judgment_role`.")

  keep_cols <- setdiff(names(crit), c(item_cols, role))
  long <- .mpa_melt_items(crit, item_cols, keep_cols)
  names(long)[names(long) == pid] <- "perceiver"
  names(long)[names(long) == tid] <- "target"
  long$item_index <- match(long$item, item_cols)

  ## group: from dyads (already in `long`) or from perceivers
  group_in_perceivers <- FALSE
  if (!is.null(group)) {
    if (group %in% names(long)) {
      names(long)[names(long) == group] <- "group"
    } else if (!is.null(perceivers) && group %in% names(perceivers)) {
      group_in_perceivers <- TRUE
    } else {
      stop("Group column `", group, "` not found in `judgments` or `perceivers`.")
    }
  }
  dyad_vars <- setdiff(names(long), c("perceiver", "target", "item", "item_index", "y", "group"))

  if (!is.null(rating_scale)) {
    if (length(rating_scale) != 2 || !is.numeric(rating_scale)) stop("`rating_scale` must be c(min, max).")
    rng <- range(long$y, na.rm = TRUE)
    if (rng[1] < rating_scale[1] || rng[2] > rating_scale[2])
      warning("Observed judgments (", rng[1], "-", rng[2], ") lie outside `rating_scale` (",
              rating_scale[1], "-", rating_scale[2], "). Check the scale definition.")
  }

  ## ---- criterion sources: target x item tables ------------------------------
  criterion_n  <- list()
  ti       <- NULL
  src_names <- character(0)
  for (nm in names(criterion_spec)) {
    res <- .mpa_criterion_table(criterion_spec[[nm]], nm, dyads, roles, targets, tid, iid, item_cols, msg)
    criterion_n[[nm]] <- res$n_raters
    tab <- res$table
    names(tab)[names(tab) == "value"] <- nm
    ti <- if (is.null(ti)) tab else merge(ti, tab, by = c("target", "item"), all = TRUE)
    src_names <- c(src_names, nm)
  }
  if (!is.null(ti)) {
    grid <- unique(data.frame(target = long$target, item = long$item, stringsAsFactors = FALSE))
    ti <- merge(ti, grid, by = c("target", "item"), all = TRUE)
  }
  if (length(src_names) >= 2 && criterion_combine == "mean") {
    cm <- rowMeans(ti[, src_names, drop = FALSE], na.rm = TRUE)
    cm[is.nan(cm)] <- NA_real_
    ti$combined <- cm
    src_names <- c(src_names, "combined")
  }

  ## ---- further predictor profiles: perceiver x item or target x item tables ----
  if (length(bad <- intersect(names(profile_spec), c(src_names, "combined"))))
    stop("Profile name(s) already used by a criterion source: ", paste(bad, collapse = ", "))
  profile_tabs <- list(); profile_n <- list()
  profile_key <- vapply(profile_spec, function(sp) sp$key, character(1))
  for (nm in names(profile_spec)) {
    sp <- profile_spec[[nm]]
    res <- if (sp$key == "perceiver")
      .mpa_criterion_table(sp, nm, dyads, roles, perceivers, pid, iid, item_cols, msg, file = "perceivers")
    else
      .mpa_criterion_table(sp, nm, dyads, roles, targets, tid, iid, item_cols, msg)
    profile_n[[nm]] <- res$n_raters
    tab <- res$table
    names(tab)[names(tab) == "value"] <- nm
    profile_tabs[[nm]] <- tab                      # columns: target (the owner's id), item, <nm>
  }

  ## ---- normative and distinctive profiles -------------------------------
  norm_method <- if (is.data.frame(norm)) "external" else match.arg(norm, c("loo", "sample"))
  by_vec <- NULL
  if (!is.null(norm_by) && !is.null(ti)) {
    if (is.null(targets) || !norm_by %in% names(targets))
      stop("`norm_by` must name a column in `targets`.")
    by_vec <- as.character(targets[[norm_by]])[match(ti$target, targets[[tid]])]
    if (anyNA(by_vec))
      warning("Some targets have missing `", norm_by, "`; their normative profiles are set to NA.")
  }
  if (center_norm == "midpoint" && is.null(rating_scale))
    stop("`center_norm = \"midpoint\"` requires `rating_scale = c(min, max)`.")
  midpoint <- if (!is.null(rating_scale)) mean(rating_scale) else NA_real_
  if (!is.null(ti)) {
    for (src in src_names) {
      ti <- .mpa_add_norm(ti, src, norm_method, by_vec, norm, item_cols, center_norm, midpoint)
    }
    long <- merge(long, ti, by = c("target", "item"), all.x = TRUE, sort = FALSE)
  }
  for (nm in names(profile_tabs)) {
    tab <- profile_tabs[[nm]]
    ## a profile's norm is the mean of that profile across its owners, or the external norm
    by_p <- if (profile_key[[nm]] == "target" && !is.null(by_vec)) as.character(targets[[norm_by]])[match(tab$target, targets[[tid]])] else NULL
    tab <- .mpa_add_norm(tab, nm, norm_method, by_p, norm, item_cols, "none", midpoint)
    tab <- tab[, c("target", "item", nm, paste0(nm, c(".norm", ".dist"))), drop = FALSE]
    names(tab)[1] <- profile_key[[nm]]
    long <- merge(long, tab, by = c(profile_key[[nm]], "item"), all.x = TRUE, sort = FALSE)
  }

  ## ---- merge item, target, perceiver variables --------------------------
  item_vars <- character(0)
  if (!is.null(items)) {
    items <- as.data.frame(items, stringsAsFactors = FALSE)
    if (iid %in% names(items)) {
      items[[iid]] <- as.character(items[[iid]])
      if (!all(item_cols %in% items[[iid]]))
        stop("Not all `item_cols` are found in `items$", iid, "`.")
      items$.item <- items[[iid]]
      items <- items[, setdiff(names(items), iid), drop = FALSE]
    } else {
      if (nrow(items) != n_items)
        stop("`items` has ", nrow(items), " rows but there are ", n_items,
             " item columns; add an item identifier column or match the rows.")
      msg("`items` has no `", iid, "` column: rows are matched to item columns by position.")
      items$.item <- item_cols
    }
    item_vars <- setdiff(names(items), ".item")
    clash <- intersect(item_vars, names(long))
    if (length(clash)) stop("Item variable name(s) clash with existing columns: ", paste(clash, collapse = ", "))
    names(items)[names(items) == ".item"] <- "item"
    long <- merge(long, items, by = "item", all.x = TRUE, sort = FALSE)
  }

  target_vars <- character(0)
  if (!is.null(targets)) {
    tv <- targets[, setdiff(names(targets), item_cols), drop = FALSE]
    target_vars <- setdiff(names(tv), tid)
    clash <- intersect(target_vars, names(long))
    if (length(clash)) stop("Target variable name(s) clash with existing columns: ", paste(clash, collapse = ", "))
    if (length(target_vars)) {
      names(tv)[names(tv) == tid] <- "target"
      long <- merge(long, tv, by = "target", all.x = TRUE, sort = FALSE)
    }
  }

  perceiver_vars <- character(0)
  if (!is.null(perceivers)) {
    pv <- perceivers[, setdiff(names(perceivers), item_cols), drop = FALSE]   # item columns are self-ratings (see `profiles`)
    perceiver_vars <- setdiff(names(pv), c(pid, if (group_in_perceivers) group))
    clash <- intersect(perceiver_vars, names(long))
    if (length(clash)) stop("Perceiver variable name(s) clash with existing columns: ", paste(clash, collapse = ", "))
    names(pv)[names(pv) == pid] <- "perceiver"
    if (group_in_perceivers) names(pv)[names(pv) == group] <- "group"
    if (ncol(pv) > 1) long <- merge(long, pv, by = "perceiver", all.x = TRUE, sort = FALSE)
  }

  ## ---- level detection for variables that came with the judgments --------
  if (detect_levels && length(dyad_vars)) {
    det <- .mpa_detect_levels(long, dyad_vars)
    if (length(det$perceiver)) {
      perceiver_vars <- c(perceiver_vars, det$perceiver)
      msg("Constant within perceivers, registered as perceiver-level: ", paste(det$perceiver, collapse = ", "))
    }
    if (length(det$target)) {
      target_vars <- c(target_vars, det$target)
      msg("Constant within targets, registered as target-level: ", paste(det$target, collapse = ", "))
    }
    if (length(det$constant)) msg("Constant overall (not informative): ", paste(det$constant, collapse = ", "))
    if (isTRUE(det$pairs) && length(det$perceiver))
      msg("Pairs design: every perceiver has one target and every target one perceiver, so these variables describe the pair.")
    dyad_vars <- det$dyad
  }

  ## ---- rename / rescale -------------------------------------------------
  ## the item SDs of the criterion sources are item-level variables like any other
  ## (no item random slope; usable as moderators, e.g. in a sensitivity analysis)
  sd_vars <- if (!is.null(ti)) paste0(src_names, ".sd") else character(0)
  var_roles <- list(item = c(item_vars, sd_vars), dyad = dyad_vars, target = target_vars, perceiver = perceiver_vars)
  if (!is.null(rename)) {
    if (is.null(names(rename)) || any(names(rename) == ""))
      stop("`rename` must be a named character vector c(new = \"old\").")
    miss <- setdiff(rename, names(long))
    if (length(miss)) stop("`rename`: column(s) not found: ", paste(miss, collapse = ", "))
    names(long)[match(rename, names(long))] <- names(rename)
    var_roles <- lapply(var_roles, function(v) {
      hit <- rename %in% v
      v[match(rename[hit], v)] <- names(rename)[hit]
      v
    })
  }
  rescaled <- character(0)
  if (!is.null(rescale)) {
    if (is.null(names(rescale))) stop("`rescale` must be a named list, e.g. list(DES = c(0, 100)).")
    for (nm in names(rescale)) {
      if (!nm %in% names(long)) stop("`rescale`: column `", nm, "` not found.")
      rng <- rescale[[nm]]
      if (length(rng) != 2 || !is.numeric(rng)) stop("`rescale$", nm, "` must be c(min, max).")
      new <- paste0(nm, ".n")
      long[[new]] <- (as.numeric(long[[nm]]) - rng[1]) / (rng[2] - rng[1]) * 2 - 1
      rescaled <- c(rescaled, new)
      for (r in names(var_roles)) if (nm %in% var_roles[[r]]) var_roles[[r]] <- c(var_roles[[r]], new)
    }
  }

  decomposed <- character(0)
  var_roles$group <- character(0)
  if (!is.null(decompose)) {
    key <- paste(long$perceiver, long$target, sep = "\r")
    ## Decompose within groups when the design has them: in a half-block design
    ## perceivers and targets are crossed only within their group, so perceiver
    ## and target means must be taken relative to the group mean. Otherwise a
    ## group difference would enter the .p and .t components twice and the
    ## relationship component with the opposite sign.
    by_group <- isTRUE(decompose_within_group) && "group" %in% names(long)
    grp <- if (by_group) as.character(long$group) else rep("all", nrow(long))
    for (v in decompose) {
      if (!v %in% names(long)) stop("`decompose`: column `", v, "` not found.")
      if (!v %in% var_roles$dyad) warning("`", v, "` is not registered as a dyad-level variable; decomposing anyway.")
      dy <- data.frame(perceiver = as.character(long$perceiver), target = as.character(long$target),
                       group = grp, x = long[[v]], stringsAsFactors = FALSE)[!duplicated(key), ]
      dy <- dy[!is.na(dy$x), ]
      gm <- mean(dy$x)
      mg <- tapply(dy$x, dy$group, mean)
      ## a perceiver (target) belongs to one group; its mean is taken relative to that group's mean
      mp <- tapply(dy$x, dy$perceiver, mean) - mg[tapply(dy$group, dy$perceiver, `[`, 1)]
      mt <- tapply(dy$x, dy$target, mean)    - mg[tapply(dy$group, dy$target, `[`, 1)]
      comp_g <- as.numeric(mg[grp]) - gm
      comp_p <- as.numeric(mp[as.character(long$perceiver)])
      comp_t <- as.numeric(mt[as.character(long$target)])
      comp_d <- long[[v]] - gm - comp_g - comp_p - comp_t
      if (all(abs(comp_d) < 1e-10, na.rm = TRUE) || all(abs(comp_p) < 1e-10, na.rm = TRUE))
        warning("`", v, "`: relationship component is degenerate - perceivers and targets are not crossed in this design.")
      long[[paste0(v, ".p")]] <- comp_p
      long[[paste0(v, ".t")]] <- comp_t
      long[[paste0(v, ".r")]] <- comp_d
      var_roles$perceiver <- c(var_roles$perceiver, paste0(v, ".p"))
      var_roles$target    <- c(var_roles$target, paste0(v, ".t"))
      var_roles$dyad      <- c(var_roles$dyad, paste0(v, ".r"))
      if (by_group) {
        long[[paste0(v, ".g")]] <- comp_g
        var_roles$group <- c(var_roles$group, paste0(v, ".g"))
      }
      decomposed <- c(decomposed, v)
    }
    if (by_group) msg("Dyadic variable(s) decomposed within groups: ", paste(decompose, collapse = ", "),
                      " (components .g, .p, .t, .r)")
  }

  if (is.character(criterion_scale)) {
    criterion_scale <- match.arg(criterion_scale, c("auto", "sd"))
    if (criterion_scale == "auto") {
      if (!is.null(rating_scale)) {
        criterion_scale <- diff(rating_scale) / 2
      } else {
        if (length(src_names)) msg("`criterion_scale = \"auto\"` needs `rating_scale`; criterion components are not rescaled.")
        criterion_scale <- 1
      }
    } else {
      dv <- unlist(lapply(src_names, function(src) long[[paste0(src, ".dist")]]))
      criterion_scale <- if (length(dv) && any(!is.na(dv))) stats::sd(dv, na.rm = TRUE) else 1
    }
  }
  if (!is.numeric(criterion_scale) || length(criterion_scale) != 1 || !is.finite(criterion_scale) || criterion_scale <= 0)
    stop("`criterion_scale` must be \"auto\", \"sd\" or a positive number.")

  ## one line per level, so that the messages say where every predictor was registered
  lv <- c(item = "item", dyad = "dyad", perceiver = "perceiver", target = "target", group = "group")
  reg <- vapply(names(lv), function(u) if (length(var_roles[[u]])) paste0(lv[[u]], ": ", paste(var_roles[[u]], collapse = ", ")) else "", character(1))
  if (any(nzchar(reg))) msg("Predictors registered by level: ", paste(reg[nzchar(reg)], collapse = "; "))
  if (criterion_scale != 1) for (src in src_names) {
    long[[paste0(src, ".norm.c")]] <- long[[paste0(src, ".norm.c")]] / criterion_scale
    long[[paste0(src, ".dist")]]   <- long[[paste0(src, ".dist")]] / criterion_scale
  }
  if (criterion_scale != 1) for (nm in names(profile_spec))
    long[[paste0(nm, ".dist")]] <- long[[paste0(nm, ".dist")]] / criterion_scale
  if (length(profile_spec))
    msg("Further predictor profiles: ", paste0(names(profile_spec), " (", profile_key, " x item)", collapse = ", "))

  derived <- character(0)
  if (!is.null(derive)) {
    if (is.null(names(derive))) stop("`derive` must be a named list.")
    for (nm in names(derive)) {
      ex <- derive[[nm]]
      if (inherits(ex, "formula")) {
        if (length(ex) != 2) stop("`derive$", nm, "` must be a one-sided formula, e.g. ~ abs(DES.n).")
        vals <- eval(ex[[2]], envir = long, enclos = environment(ex))
        used <- all.vars(ex[[2]])
      } else if (is.function(ex)) {
        vals <- ex(long)
        used <- character(0)
      } else stop("`derive$", nm, "` must be a formula or a function.")
      if (length(vals) != nrow(long)) stop("`derive$", nm, "` did not return one value per row.")
      long[[nm]] <- vals
      lv <- names(var_roles)[vapply(var_roles, function(v) any(used %in% v), logical(1))]
      if (length(lv) == 1) var_roles[[lv]] <- c(var_roles[[lv]], nm)
      else if (length(used)) warning("Level of derived variable `", nm, "` is ambiguous; not registered (no identifiability filtering).")
      derived <- c(derived, nm)
    }
  }

  ## ---- finalize ---------------------------------------------------------
  if (na_action == "drop") {
    n_before <- nrow(long)
    long <- long[!is.na(long$y), , drop = FALSE]
    msg("Dropped ", n_before - nrow(long), " rows with missing judgment.")
  }
  long <- long[order(long$perceiver, long$target, long$item_index), , drop = FALSE]
  rownames(long) <- NULL
  long$perceiver <- factor(long$perceiver)
  long$target    <- factor(long$target)
  long$item      <- factor(long$item, levels = item_cols)
  if ("group" %in% names(long)) long$group <- factor(long$group)

  front <- c("perceiver", "target", if ("group" %in% names(long)) "group", "item", "item_index", "y",
             unlist(lapply(src_names, function(s) paste0(s, c("", ".norm", ".norm.c", ".dist", ".sd")))),
             unlist(lapply(names(profile_spec), function(s) paste0(s, c("", ".norm", ".dist")))))
  long <- long[, c(front, setdiff(names(long), front)), drop = FALSE]

  design <- mpa_design(long)
  ipd <- table(paste(long$perceiver, long$target)[!is.na(long$y)])
  design$items_per_dyad <- if (length(ipd)) range(ipd) else c(NA, NA)

  attr(long, "mpa") <- list(
    id        = id,
    item_cols = item_cols,
    vars      = c(list(judgment = "y",
                       criterion  = src_names,
                       norm   = paste0(src_names, ".norm"),
                       norm_c = paste0(src_names, ".norm.c"),
                       dist   = paste0(src_names, ".dist"),
                       sd     = sd_vars,
                       profiles = names(profile_spec),
                       profile_key = profile_key),
                  var_roles,
                  list(rescaled = rescaled, derived = derived, decomposed = decomposed)),
    settings  = list(norm = norm_method, norm_by = norm_by, center_norm = center_norm,
                     rating_scale = rating_scale, criterion_combine = criterion_combine, criterion_scale = criterion_scale,
                     judgment_role = if (is.null(role)) NULL else if (is.null(judgment_role)) unique(roles[crit_idx]) else judgment_role,
                     rescale = rescale),
    criterion_n   = criterion_n,
    profile_n = profile_n,
    design    = design,
    created   = Sys.time()
  )
  class(long) <- c("mpa_data", "data.frame")
  long
}


## ---- internal helpers ----------------------------------------------------

.mpa_complete_id <- function(id) {
  default <- c(perceiver = "perceiver", target = "target", item = "item")
  if (is.null(names(id))) {
    if (length(id) != 3) stop("`id` must be a named vector or have exactly 3 elements (perceiver, target, item).")
    names(id) <- names(default)
  }
  bad <- setdiff(names(id), names(default))
  if (length(bad)) stop("Unknown `id` element(s): ", paste(bad, collapse = ", "))
  default[names(id)] <- id
  default
}

.mpa_resolve_item_cols <- function(item_cols, dyads, items, iid) {
  if (is.null(item_cols)) {
    if (!is.null(items) && iid %in% names(items)) {
      cand <- as.character(items[[iid]])
      if (all(cand %in% names(dyads))) return(cand)
    }
    stop("`item_cols` must be specified (column names or a regular expression such as \"^item_\").")
  }
  if (length(item_cols) == 1 && !item_cols %in% names(dyads)) {
    found <- grep(item_cols, names(dyads), value = TRUE)
    if (!length(found)) stop("No columns in `judgments` match the regular expression `", item_cols, "`.")
    return(found)
  }
  miss <- setdiff(item_cols, names(dyads))
  if (length(miss)) stop("Item column(s) not found in `judgments`: ", paste(miss, collapse = ", "))
  item_cols
}

.mpa_melt_items <- function(df, item_cols, keep) {
  n <- nrow(df)
  out <- df[rep(seq_len(n), times = length(item_cols)), keep, drop = FALSE]
  out$item <- rep(item_cols, each = n)
  vals <- lapply(df[item_cols], function(col) suppressWarnings(as.numeric(as.character(col))))
  out$y <- unlist(vals, use.names = FALSE)
  rownames(out) <- NULL
  out
}

.mpa_long_to_wide <- function(df, pid, tid, iid, rating_col, item_cols, items, msg) {
  for (v in c(iid, rating_col))
    if (!v %in% names(df)) stop("`format = \"long\"` requires column `", v, "` in `dyads`.")
  it <- as.character(df[[iid]])
  if (is.null(item_cols)) {
    item_cols <- if (!is.null(items) && iid %in% names(items)) as.character(items[[iid]]) else unique(it)
  }
  unknown <- setdiff(unique(it), item_cols)
  if (length(unknown)) stop("Item identifier(s) in `dyads` not covered by `item_cols`/`items`: ",
                            paste(head(unknown, 5), collapse = ", "), if (length(unknown) > 5) ", ...")
  key <- paste(df[[pid]], df[[tid]], sep = "\r")
  if (anyDuplicated(paste(key, it, sep = "\r")))
    stop("Long data contain duplicated perceiver x target x item rows.")
  ukey <- unique(key)
  mat <- matrix(NA_real_, length(ukey), length(item_cols), dimnames = list(NULL, item_cols))
  mat[cbind(match(key, ukey), match(it, item_cols))] <- suppressWarnings(as.numeric(as.character(df[[rating_col]])))
  other <- setdiff(names(df), c(iid, rating_col))
  first <- !duplicated(key)
  meta <- df[first, other, drop = FALSE]
  meta <- meta[match(ukey, key[first]), , drop = FALSE]
  ## dyad-level variables must be constant within dyad
  for (v in setdiff(other, c(pid, tid))) {
    nun <- tapply(df[[v]], key, function(z) length(unique(z[!is.na(z)])))
    if (any(nun > 1, na.rm = TRUE))
      warning("Column `", v, "` varies within perceiver x target pairs in long data; first value used.")
  }
  msg("Long data cast to wide: ", length(ukey), " rater x target rows, ", length(item_cols), " items.")
  wide <- cbind(meta, as.data.frame(mat, stringsAsFactors = FALSE))
  rownames(wide) <- NULL
  list(wide = wide, item_cols = item_cols)
}

.mpa_loop_to_wide <- function(df, pid, tid, pattern, loop_id, msg) {
  if (!pid %in% names(df)) stop("Column `", pid, "` not found in `judgments`.")
  if (anyDuplicated(df[[pid]])) stop("`format = \"loop\"` expects one row per perceiver.")
  cols <- names(df)
  m <- regmatches(cols, regexec(pattern, cols))
  hit <- lengths(m) == 3
  if (!any(hit)) stop("No columns match `loop_pattern` (", pattern, ").")
  slots <- vapply(m[hit], `[`, "", 2); vars <- vapply(m[hit], `[`, "", 3)
  loop_cols <- cols[hit]; other <- cols[!hit]
  if (tid %in% other) stop("`judgments` has a column `", tid, "` although `format = \"loop\"`; the target is taken from the loop.")
  if (!is.null(loop_id) && !loop_id %in% vars) stop("`loop_id` = `", loop_id, "` is not among the looped variables.")
  all_vars <- unique(vars)
  blocks <- lapply(unique(slots), function(sl) {
    sel <- loop_cols[slots == sl]; vn <- vars[slots == sl]
    block <- df[, sel, drop = FALSE]; names(block) <- vn
    for (v in setdiff(all_vars, vn)) block[[v]] <- NA
    block <- block[, all_vars, drop = FALSE]
    block[[tid]] <- if (is.null(loop_id)) sl else as.character(block[[loop_id]])
    cbind(df[, other, drop = FALSE], block, stringsAsFactors = FALSE)
  })
  wide <- do.call(rbind, blocks)
  if (!is.null(loop_id)) wide[[loop_id]] <- NULL
  keep <- !is.na(wide[[tid]]) & wide[[tid]] != ""
  wide <- wide[keep, , drop = FALSE]
  rownames(wide) <- NULL
  msg("Loop data cast to wide: ", nrow(df), " perceivers x ", length(unique(slots)), " loops -> ",
      nrow(wide), " perceiver x target rows; looped variables: ", paste(all_vars, collapse = ", "))
  wide
}

## classify dyad-level variables by whether they vary within perceivers / targets
.mpa_detect_levels <- function(long, vars) {
  key <- paste(long$perceiver, long$target, sep = "\r")
  pt <- long[!duplicated(key), c("perceiver", "target", vars), drop = FALSE]
  n_per_p <- table(pt$perceiver); n_per_t <- table(pt$target)
  can_p <- any(n_per_p > 1); can_t <- any(n_per_t > 1)     # detection needs replication
  nun <- function(z) length(unique(z[!is.na(z)]))
  out <- list(perceiver = character(0), target = character(0), constant = character(0), dyad = character(0))
  for (v in vars) {
    if (!is.numeric(pt[[v]]) && !is.character(pt[[v]]) && !is.factor(pt[[v]]) && !is.logical(pt[[v]])) { out$dyad <- c(out$dyad, v); next }
    ## In a one-with-many design every dyad-level variable is constant within
    ## perceivers (each perceiver has one pair) and therefore equivalent to a
    ## perceiver-level variable for the purpose of identifiability; the same
    ## holds for targets in the reverse design.
    const_p <- (!can_p) || all(tapply(pt[[v]], pt$perceiver, nun) <= 1)
    const_t <- (!can_t) || all(tapply(pt[[v]], pt$target, nun) <= 1)
    ## pairs design (no replication on either side): a variable that differs
    ## between pairs describes the pair, which is the person-level unit there
    if (!can_p && !can_t && nun(pt[[v]]) > 1) { out$perceiver <- c(out$perceiver, v); out$pairs <- TRUE; next }
    if (const_p && const_t)       out$constant  <- c(out$constant, v)
    else if (const_p)             out$perceiver <- c(out$perceiver, v)
    else if (const_t)             out$target    <- c(out$target, v)
    else                          out$dyad      <- c(out$dyad, v)
  }
  out$dyad <- c(out$dyad, out$constant)      # keep constants in the data, but as uninformative dyad variables
  out
}

.mpa_parse_criterion <- function(criterion) {
  if (is.null(criterion)) return(list())
  if (!is.list(criterion) || is.null(names(criterion)) || any(names(criterion) == ""))
    stop("`criterion` must be a named list, e.g. list(self = \"targets\", other = \"TNI\").")
  lapply(criterion, function(x) {
    if (is.data.frame(x)) return(list(type = "data", data = x))
    if (is.character(x) && length(x) == 1) {
      if (x == "targets") return(list(type = "targets"))
      return(list(type = "role", role = x, max_n = NULL, fun = mean))
    }
    if (is.list(x) && !is.null(x$role)) {
      return(list(type = "role", role = x$role,
                  max_n = x$max_n,
                  fun = if (is.null(x$fun)) mean else match.fun(x$fun)))
    }
    stop("Invalid criterion specification. See ?mpa_data.")
  })
}

## further predictor profiles: like criterion sources, plus the owner (`key`) of each profile
.mpa_parse_profiles <- function(profiles, pid, tid) {
  if (is.null(profiles)) return(list())
  if (!is.list(profiles) || is.data.frame(profiles) || is.null(names(profiles)) || any(names(profiles) == ""))
    stop("`profiles` must be a named list, e.g. list(pself = \"perceivers\").")
  lapply(profiles, function(x) {
    if (is.data.frame(x)) {
      key <- if (pid %in% names(x) && !tid %in% names(x)) "perceiver" else "target"
      return(list(type = "data", data = x, key = key))
    }
    if (is.character(x) && length(x) == 1) {
      if (x == "perceivers") return(list(type = "targets", key = "perceiver"))
      if (x == "targets") return(list(type = "targets", key = "target"))
      return(list(type = "role", role = x, max_n = NULL, fun = mean, key = "target"))
    }
    stop("Invalid profile specification. See ?mpa_data.")
  })
}

.mpa_aggregate_ti <- function(lg, fun = mean) {
  f <- function(v) { v <- v[!is.na(v)]; if (!length(v)) NA_real_ else fun(v) }
  aggregate(list(value = lg$y), by = list(target = lg$target, item = lg$item), FUN = f)
}

## `targets`/`tid` are the file and identifier of the profile's owners: the
## targets for criterion sources, the perceivers for perceiver profiles. The
## owner's id is returned in the column `target` either way.
.mpa_criterion_table <- function(spec, name, dyads, roles, targets, tid, iid, item_cols, msg, file = "targets") {
  if (spec$type == "role") {
    sub <- dyads[roles %in% spec$role, c(tid, item_cols), drop = FALSE]
    n_raters <- table(sub[[tid]])
    if (!is.null(spec$max_n)) {
      bad <- names(n_raters)[n_raters > spec$max_n]
      if (length(bad)) {
        msg("Criterion source `", name, "`: ", length(bad), " target(s) have more than ",
            spec$max_n, " rater(s) and are set to NA.")
        sub <- sub[!(sub[[tid]] %in% bad), , drop = FALSE]
      }
    }
    lg <- .mpa_melt_items(sub, item_cols, tid)
    names(lg)[names(lg) == tid] <- "target"
  } else if (spec$type == "targets") {
    if (is.null(targets)) stop("Source `", name, "` refers to `", file, "`, but `", file, "` is NULL.")
    if (!all(item_cols %in% names(targets)))
      stop("Source `", name, "`: `", file, "` does not contain all item columns.")
    sub <- targets[, c(tid, item_cols), drop = FALSE]
    n_raters <- table(sub[[tid]])
    lg <- .mpa_melt_items(sub, item_cols, tid)
    names(lg)[names(lg) == tid] <- "target"
  } else {  # data frame
    d <- as.data.frame(spec$data, stringsAsFactors = FALSE)
    if (!tid %in% names(d)) stop("Source `", name, "`: data frame lacks the id column `", tid, "`.")
    d[[tid]] <- as.character(d[[tid]])
    if (all(item_cols %in% names(d))) {
      sub <- d[, c(tid, item_cols), drop = FALSE]
      n_raters <- table(sub[[tid]])
      lg <- .mpa_melt_items(sub, item_cols, tid)
      names(lg)[names(lg) == tid] <- "target"
    } else {
      if (!all(c(iid, "value") %in% names(d)))
        stop("Criterion source `", name, "`: long data frame needs columns `", tid, "`, `", iid, "` and `value`.")
      it <- as.character(d[[iid]])
      if (!all(it %in% item_cols)) {
        if (is.numeric(d[[iid]]) && all(d[[iid]] %in% seq_along(item_cols))) it <- item_cols[d[[iid]]]
        else stop("Criterion source `", name, "`: item identifiers do not match `item_cols`.")
      }
      lg <- data.frame(target = d[[tid]], item = it, y = as.numeric(d$value), stringsAsFactors = FALSE)
      n_raters <- apply(table(lg$target, lg$item), 1, max)
    }
  }
  list(table = .mpa_aggregate_ti(lg, spec$fun %||% mean), n_raters = n_raters)
}

`%||%` <- function(a, b) if (is.null(a)) b else a

.mpa_add_norm <- function(ti, src, method, by_vec, norm_ext, item_cols, center_norm, midpoint) {
  v  <- ti[[src]]
  ok <- !is.na(v)
  by <- if (is.null(by_vec)) rep("all", nrow(ti)) else by_vec
  key <- data.frame(item = as.character(ti$item), by = by, stringsAsFactors = FALSE)

  if (method == "external") {
    if (!"item" %in% names(norm_ext)) stop("External `norm` data frame needs an `item` column.")
    col <- if (src %in% names(norm_ext)) src else if ("norm" %in% names(norm_ext)) "norm" else
      stop("External `norm` has no column `", src, "` or `norm`.")
    ext_item <- as.character(norm_ext$item)
    if (!all(ext_item %in% item_cols)) {
      if (is.numeric(norm_ext$item) && all(norm_ext$item %in% seq_along(item_cols))) ext_item <- item_cols[norm_ext$item]
      else stop("External `norm`: item identifiers do not match `item_cols`.")
    }
    nrm <- as.numeric(norm_ext[[col]])[match(key$item, ext_item)]
    gm_row <- rep(mean(as.numeric(norm_ext[[col]]), na.rm = TRUE), nrow(ti))
  } else {
    st <- aggregate(list(S = ifelse(ok, v, 0), N = as.numeric(ok)), by = key, FUN = sum)
    st$sn <- st$S / st$N
    gm <- aggregate(list(gm = st$sn), by = list(by = st$by), FUN = function(z) mean(z[is.finite(z)]))
    idx <- match(paste(key$item, key$by), paste(st$item, st$by))
    S <- st$S[idx]; N <- st$N[idx]
    nrm <- if (method == "loo") ifelse(ok, (S - v) / (N - 1), S / N) else S / N
    nrm[!is.finite(nrm)] <- NA_real_
    gm_row <- gm$gm[match(key$by, gm$by)]
  }
  if (!is.null(by_vec)) nrm[is.na(by_vec)] <- NA_real_

  ## item variance of the criterion: SD across targets per item (and `by` level),
  ## an item-level descriptor of how much targets differ on the item
  sdt <- aggregate(list(sd = v), by = key, FUN = function(z) { z <- z[!is.na(z)]; if (length(z) < 2) NA_real_ else stats::sd(z) })
  sd_item <- sdt$sd[match(paste(key$item, key$by), paste(sdt$item, sdt$by))]
  if (!is.null(by_vec)) sd_item[is.na(by_vec)] <- NA_real_

  ti[[paste0(src, ".norm")]]   <- nrm
  ti[[paste0(src, ".norm.c")]] <- switch(center_norm,
                                         grand_mean = nrm - gm_row,
                                         midpoint   = nrm - midpoint,
                                         none       = nrm)
  ti[[paste0(src, ".dist")]]   <- v - nrm
  ti[[paste0(src, ".sd")]]     <- sd_item
  ti
}
