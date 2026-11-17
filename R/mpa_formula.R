#' Build a multilevel profile analysis model formula
#'
#' Constructs the `lme4` formula for a multilevel profile analysis from an
#' `mpa_data` object. The judgments `y` are predicted by the normative and
#' distinctive components of a criterion profile -- note that the criterion is
#' a *predictor* here, not the dependent variable -- and by further predictors
#' at the item, perceiver, target or dyad level. A predictor can have a
#' *main effect* on the judgments (`predictors`, `random_slopes`) and/or a
#' *moderating effect* on profile agreement (`interactions`). Random effects are
#' specified per unit (perceivers, targets, dyads, items, groups, and target x
#' item cells); random slopes that are not identifiable at a given unit (because
#' the variable does not vary within that unit) are dropped automatically and
#' reported.
#'
#' The tutorial "steps" map onto the arguments as follows (variable names from
#' the tutorial data, see `vignette("mpa-tutorial")`):
#' * Step 0 (variance decomposition): `criterion = NULL` -- random intercepts
#'   for perceivers, targets, dyads and items
#' * Step 1 (base model): `criterion = "self"` -- the criterion components with
#'   random slopes for perceivers, targets and dyads, plus random item
#'   intercepts and random item slopes for the distinctive component
#' * Step 2 (main effects of an item property): `random_slopes = "des"`
#' * Step 3 (moderation of agreement by an item property):
#'   `random_slopes = "des", poly = "des", interactions = "criterion"`
#' * Step 4 (main and moderating effects of perceiver, target and dyad-level
#'   variables): `predictors = c("know", "training", "talkativeness"), interactions = "criterion"`
#'
#' One rule governs the random part: every unit named in `random_units` gets a
#' random intercept and a random slope for every term (the criterion
#' components and `random_slopes`) that varies within that unit. Items
#' therefore get an intercept and a slope for the distinctive component (the
#' normative profile and item properties take one value per item and get
#' none), and a group gets an intercept only. Item random slopes matter
#' whenever item properties are entered as moderators of agreement: the
#' distinctive-agreement slope varies between items over and above what the
#' item properties explain, and if that variation is not modelled, the
#' standard errors of item property x criterion terms are computed as if the
#' judgments rather than the items were the units of the item-level inference
#' (Judd, Westfall & Kenny, 2012). To fit the classical SAM without item
#' effects, leave `"item"` out of `random_units`; to keep item intercepts but
#' drop the item slope, use [mpa_update()], e.g.
#' `mpa_update(f, drop = list(item = "self.dist"))`.
#'
#' A model can contain more than one *predictor profile*. The criterion is
#' one: a profile on the same items as the judgments, whose slope is
#' agreement. `profiles` names further ones, for example the perceivers'
#' self-ratings (Human & Biesanz, 2011). Their slope shows whether perceivers
#' ascribe to targets what is distinctive about themselves. Because the
#' normative profile and the target's distinctive profile are in the same
#' model, this is assumed similarity beyond actual similarity, not the
#' classical index that compares judgments and self-ratings without these
#' controls (see `vignette("predictor-profiles")`, Section 2.4). A
#' further profile enters with its distinctive component only, because the
#' normative profiles of two profiles on the same items are nearly collinear.
#' The component is treated like the distinctive component of the criterion.
#' It gets a fixed slope, random slopes at every unit within which it varies,
#' the same scaling per scale point, and, with `interactions = "criterion"`,
#' an interaction with every predictor. Interactions between two profiles are
#' not formed; add them with `extra_terms` if needed. A profile is entered as
#' `profiles` rather than as `predictors` because a predictor takes one value
#' per perceiver, target, dyad or item, whereas a profile takes one value per
#' person and item.
#'
#' In a pairs design every perceiver has one target and every target one
#' perceiver (e.g., Park & MacDonald, 2019). Perceivers, targets and dyads are
#' then the same units. One set of person-level random effects is kept, under
#' the name of the first of `"perceiver"` and `"target"` in `random_units`,
#' and variables that describe the pair get no random slope there.
#'
#' `"target:item"` in `random_units` adds a random intercept for the target x
#' item cells: the part of the judgments of a target on an item that its
#' perceivers share beyond the criterion profile (their consensus, or the
#' target's reputation on that item). The unit gets an intercept only. It is
#' identifiable when targets are judged by several perceivers and is omitted,
#' with a message, when every target has a single perceiver. It is not part of
#' the default structure because the Social Accuracy Model does not contain it;
#' the accompanying article discusses it under *Extensions*. To
#' add it to a structure that is carried over with `random_structure`, use
#' `mpa_update(f, add = list(target_item = "1"))`.
#'
#' `"perceiver:item"` is its counterpart at the perceiver level: a random intercept
#' for the perceiver x item cells, the part of a perceiver's judgments of an item that is shared
#' across all their targets beyond the normative profile (the perceiver's
#' idiosyncratic view of people in general on that item, the item-specific
#' part of the Social Relations Model's perceiver effect). Because the
#' distinctive profile varies across the targets of one perceiver, the unit
#' also gets a random slope for the distinctive component (item-specific
#' differences between perceivers in distinctive agreement, Funder's judge x
#' trait interaction); other random slopes are not formed for it. It is
#' identifiable when perceivers judge several targets and is omitted, with a
#' message, when every perceiver judges a single target (one-with-many
#' designs). Inside the object the unit is called `perceiver_item`.
#'
#' @param data An `mpa_data` object (recommended) or a data frame in canonical
#'   format with columns `y`, `perceiver`, `target`, `item` and, optionally,
#'   `group`.
#' @param criterion Name of the criterion source whose normative and distinctive
#'   components are used (e.g., `"self"`); `NULL` for an intercept-only
#'   variance decomposition. If not given, the criterion source of `data` is
#'   used (`"combined"` if the data have several sources). The variables
#'   `<criterion>.norm.c` (or `<criterion>.norm`) and `<criterion>.dist` must
#'   exist in `data`; override with `norm` and `dist`.
#' @param norm,dist Optional explicit names of the normative and distinctive
#'   predictor columns.
#' @param y Name of the judgment (dependent variable) column.
#' @param profiles Character vector naming further predictor profiles (see
#'   Details): profiles created with the `profiles` argument of [mpa_data()]
#'   (e.g., `"pself"` for the perceivers' self-ratings) or criterion sources
#'   other than `criterion`. The column `<profile>.dist` must exist in `data`.
#' @param predictors Character vector of predictors (item, perceiver, target or
#'   dyad-level variables) whose *main effects* on the judgments are estimated
#'   as fixed effects only (random slopes for them are requested with
#'   `random_slopes`, moderating effects with `interactions`).
#' @param random_slopes Character vector of predictors whose main effects
#'   receive random slopes at every unit where they are identifiable (in
#'   addition to the fixed effect).
#' @param interactions Which *moderating effects* to include. `"none"`: main
#'   effects only. `"criterion"`: each predictor (and quadratic term) interacts
#'   with the normative and the distinctive component, i.e., moderates
#'   normative and distinctive agreement. `"all"`: additionally, all pairs of
#'   predictors from different levels interact with each other (moderation of
#'   main effects) and, as three-way terms, with the criterion components.
#' @param extra_terms Character vector of additional fixed-effect terms in
#'   `lme4` syntax (e.g., `"DES.n:KNOW.n:combined.dist"`), appended as given.
#' @param poly Character vector of predictors for which quadratic terms
#'   `I(x^2)` are added (and treated like predictors in interactions). Quadratic
#'   terms always enter as fixed effects only: being non-negative, they compete
#'   with the random intercepts, which readily produces degenerate solutions.
#' @param random_units Units that receive random effects: any of
#'   `"perceiver"`, `"target"`, `"dyad"`, `"item"`, `"group"`,
#'   `"target:item"` and `"perceiver:item"`. Each unit gets a random intercept and random slopes for
#'   the criterion components and for `random_slopes` wherever they are
#'   identifiable (i.e., vary within the unit). The default is the
#'   design-implied maximal structure for the four units of the model; add
#'   `"group"` for a random intercept of the groups of a half-block or
#'   round-robin design (a `group` column must exist), and `"target:item"` for
#'   a random intercept of the target x item cells and `"perceiver:item"` for
#'   a random intercept and distinctive slope of the perceiver x item cells
#'   (see Details). In the returned object and in [mpa_update()] these two
#'   units are called `target_item` and `perceiver_item`.
#' @param random_dyad Include dyadic (`perceiver:target`) random effects? If
#'   `NULL` (default), taken from the design detected by [mpa_data()]: they are
#'   omitted in one-with-many designs where they are not separable from
#'   perceiver or target effects.
#' @param correlated Estimate correlations among the random effects within a
#'   unit (`TRUE`, `|`) or not (`FALSE`, `||`)? Inherited from `random_structure`
#'   when that is an `mpa_formula` and `correlated` is not given explicitly.
#' @param random_structure Optional `mpa_formula` (or its `$random` list) whose
#'   random-effects structure is carried over: units and terms are taken from
#'   it instead of being derived from `random_units` and `random_dyad`, and
#'   only random slopes for new
#'   `random_slopes` are added (where identifiable). A predictor whose slopes
#'   the template already contains, or from which they were pruned, is not
#'   added again. Use this to fix the random structure once (e.g., after
#'   pruning the base model) and keep it constant across the tutorial steps,
#'   updating it only when a step adds a new random slope.
#' @param fixed Deprecated name of `predictors` (accepted until the next
#'   minor version).
#' @references Human, L. J., & Biesanz, J. C. (2011). Through the looking
#'   glass clearly: Accuracy and assumed similarity in well-adjusted
#'   individuals' first impressions. *Journal of Personality and Social
#'   Psychology, 100*(2), 349-364. \doi{10.1037/a0021850}
#'
#'   Park, Y., & MacDonald, G. (2019). Consistency between individuals' past
#'   and current romantic partners' own reports of their personalities.
#'   *Proceedings of the National Academy of Sciences, 116*(26), 12793-12797.
#'   \doi{10.1073/pnas.1902937116}
#' @param var_roles Optional named list giving the level of each predictor
#'   (`item`, `dyad`, `perceiver`, `target`, `group`); taken from `data` if it
#'   is an `mpa_data` object. Needed to filter non-identifiable random slopes
#'   (a group-level variable, such as the `.g` component of a decomposed
#'   dyadic variable, varies within no unit and gets no random slope).
#'
#' @return An object of class `mpa_formula`: a list with the `formula`, its
#'   components (`fixed_terms`, `random` per unit), the predictor table, the
#'   options used, and `dropped` (random slopes removed for identifiability).
#'   Use `print()` to inspect and [mpa_update()] to edit the random part.
#' @examples
#' dat <- mpa_simulate(n_perceiver = 16, n_item = 10, group_size = 4, seed = 1,
#'                   rating_scale = c(1, 101))
#' mpa_formula(dat, criterion = NULL)               # Step 0: intercepts at all four units
#' f1 <- mpa_formula(dat, criterion = "self")       # Step 1: the SAM plus item effects
#' f1
#' mpa_formula(dat, criterion = "self", random_units = c("perceiver", "target", "dyad"))  # the SAM
#' mpa_formula(dat, criterion = "self", random_units = c("perceiver", "target", "dyad", "item", "group"))
#' mpa_formula(dat, criterion = "self", random_units = c("perceiver", "target", "dyad", "item", "target:item"))
#' f_base <- mpa_update(f1, correlated = FALSE)
#' mpa_formula(dat, criterion = "self", random_slopes = "des", poly = "des",
#'             predictors = "know", interactions = "criterion", random_structure = f_base)
#'
#' ## a further predictor profile: the perceivers' self-ratings (assumed similarity)
#' ds <- mpa_simulate(n_perceiver = 16, n_item = 10, group_size = 4, seed = 1,
#'                    rating_scale = c(1, 101), fixed = c(intercept = 51, norm = 0.7,
#'                    dist = 0.2, similarity = 0.1))
#' mpa_formula(ds, criterion = "self", profiles = "pself")
#' @export
mpa_formula <- function(data = NULL,
                        criterion = "combined",
                        norm = NULL, dist = NULL, y = "y",
                        profiles = NULL,
                        predictors = NULL,
                        random_slopes = NULL,
                        interactions = c("none", "criterion", "all"),
                        extra_terms = NULL,
                        poly = NULL,
                        random_units = c("perceiver", "target", "dyad", "item"),
                        random_dyad = NULL,
                        correlated = TRUE,
                        random_structure = NULL,
                        var_roles = NULL,
                        fixed = NULL) {
  if (!is.null(fixed)) {                      # `fixed` was the name of `predictors` up to 0.5.7
    if (!is.null(predictors) && !identical(predictors, fixed))
      stop("`fixed` is the old name of `predictors`; give the predictors once, as `predictors`.")
    predictors <- fixed
  }

  interactions <- match.arg(interactions)
  random_units[random_units == "target_item"] <- "target:item"
  random_units[random_units == "perceiver_item"] <- "perceiver:item"
  random_units <- match.arg(random_units, c("perceiver", "target", "dyad", "item", "group", "target:item", "perceiver:item"), several.ok = TRUE)
  random_units <- .mpa_unit_key(random_units)

  meta <- if (!is.null(data)) attr(data, "mpa") else NULL
  if (is.null(var_roles) && !is.null(meta)) var_roles <- meta$vars[intersect(c("item", "dyad", "perceiver", "target", "group"), names(meta$vars))]
  cols <- if (!is.null(data)) names(data) else NULL
  if ("group" %in% random_units && !is.null(cols) && !"group" %in% cols)
    stop("`random_units` names \"group\" but there is no `group` column.")

  ## default criterion source: the one in the data (or `combined` if there are several)
  if (missing(criterion) && !is.null(meta)) {
    src <- meta$vars$criterion
    if (!length(src)) stop("The data have no criterion source; use `criterion = NULL` for an intercept-only model.")
    criterion <- if ("combined" %in% src) "combined" else src[1]
  }

  ## ---- criterion components -------------------------------------------------
  if (!is.null(criterion)) {
    if (is.null(norm)) {
      cand <- paste0(criterion, c(".norm.c", ".norm"))
      norm <- if (!is.null(cols)) cand[cand %in% cols][1] else cand[1]
      if (is.na(norm)) stop("No normative profile found for criterion source `", criterion, "`.")
    }
    if (is.null(dist)) dist <- paste0(criterion, ".dist")
    if (!is.null(cols) && !all(c(norm, dist) %in% cols))
      stop("Column(s) not found in data: ", paste(setdiff(c(norm, dist), cols), collapse = ", "))
  } else {
    norm <- dist <- NULL
  }

  ## ---- further predictor profiles: distinctive component only -------------
  prof_dist <- character(0); prof_level <- character(0)
  if (length(profiles)) {
    if (is.null(criterion)) stop("`profiles` are further predictor profiles next to a criterion; give a `criterion`.")
    if (criterion %in% profiles) stop("`", criterion, "` is the criterion and cannot also be a further profile.")
    prof_dist <- paste0(profiles, ".dist")
    if (!is.null(cols) && !all(prof_dist %in% cols))
      stop("Column(s) not found in data: ", paste(setdiff(prof_dist, cols), collapse = ", "),
           ". Create profiles with the `profiles` (or `criterion`) argument of mpa_data().")
    key <- if (!is.null(meta$vars$profile_key)) meta$vars$profile_key[profiles] else rep(NA_character_, length(profiles))
    if (!is.null(var_roles$perceiver_item)) key[prof_dist %in% var_roles$perceiver_item] <- "perceiver"
    key[is.na(key)] <- "target"                    # criterion sources belong to the target
    prof_level <- stats::setNames(paste0(key, "_item"), prof_dist)
  }

  ## ---- predictors ---------------------------------------------------------
  preds <- unique(c(predictors, random_slopes))
  profile_cols <- c(norm, dist, prof_dist, paste0(c(meta$vars$criterion, meta$vars$profiles), ".dist"))
  if (length(bad <- intersect(preds, profile_cols)))
    stop("Profile component(s) given as predictors: ", paste(bad, collapse = ", "),
         ". Name the profile in `profiles` (or `criterion`) instead.")
  if (!is.null(cols)) {
    miss <- setdiff(c(preds, poly), cols)
    if (length(miss)) stop("Predictor(s) not found in data: ", paste(miss, collapse = ", "))
  }
  if (length(setdiff(poly, preds)))
    stop("`poly` variables must also be listed in `predictors` or `random_slopes`.")

  level_of <- function(v) {
    base <- sub("^I\\((.*)\\^2\\)$", "\\1", v)
    if (!is.null(norm) && base == norm) return("item")
    if (!is.null(dist) && base == dist) return("target_item")
    if (base %in% prof_dist) return(prof_level[[base]])
    for (lv in names(var_roles)) if (base %in% var_roles[[lv]]) return(lv)
    "unknown"
  }
  polys <- if (length(poly)) paste0("I(", poly, "^2)") else character(0)
  all_preds <- c(preds, polys)
  pred_tab <- if (length(all_preds)) data.frame(
    term  = all_preds,
    level = vapply(all_preds, level_of, character(1)),
    role  = ifelse(sub("^I\\((.*)\\^2\\)$", "\\1", all_preds) %in% random_slopes, "random slope", "fixed only"),
    stringsAsFactors = FALSE, row.names = NULL) else NULL
  if (!is.null(pred_tab) && any(pred_tab$level == "unknown"))
    warning("Level of predictor(s) unknown (not registered in the mpa_data metadata): ",
            paste(pred_tab$term[pred_tab$level == "unknown"], collapse = ", "),
            ". Random slopes for them are not filtered for identifiability.")

  ## ---- fixed part ---------------------------------------------------------
  fixed_terms <- c(norm, dist, prof_dist, all_preds)
  if (interactions %in% c("criterion", "all") && length(all_preds) && !is.null(criterion)) {
    fixed_terms <- c(fixed_terms, paste0(all_preds, ":", norm), paste0(all_preds, ":", dist),
                     unlist(lapply(prof_dist, function(pd) paste0(all_preds, ":", pd))))
  }
  if (interactions == "all" && length(all_preds) > 1) {
    pairs <- utils::combn(all_preds, 2, simplify = FALSE)
    lv <- vapply(all_preds, level_of, character(1))
    pairs <- Filter(function(p) lv[[p[1]]] != lv[[p[2]]] || lv[[p[1]]] == "unknown", pairs)
    for (p in pairs) {
      fixed_terms <- c(fixed_terms, paste0(p[1], ":", p[2]))
      if (!is.null(criterion)) fixed_terms <- c(fixed_terms, paste0(p[1], ":", p[2], ":", norm), paste0(p[1], ":", p[2], ":", dist),
                                                 unlist(lapply(prof_dist, function(pd) paste0(p[1], ":", p[2], ":", pd))))
    }
  }
  fixed_terms <- unique(c(fixed_terms, extra_terms))

  ## ---- random part --------------------------------------------------------
  varies_within <- list(
    perceiver = c("item", "dyad", "target", "target_item", "perceiver_item", "unknown"),
    target    = c("item", "dyad", "perceiver", "target_item", "perceiver_item", "unknown"),
    dyad      = c("item", "target_item", "perceiver_item", "unknown"),
    item      = c("dyad", "perceiver", "target", "target_item", "perceiver_item", "unknown"),
    group     = character(0),
    target_item = character(0),                     # a shared target x item view: intercept only
    perceiver_item = "target_item"                  # a perceiver's item-specific view: intercept + distinctive slope
  )
  ## a perceiver with a single target has one value of every target and dyad
  ## variable (and a target with a single perceiver one value of every
  ## perceiver and dyad variable), so these get no random slope at that unit
  tpp <- if (!is.null(meta)) meta$design$targets_per_perceiver else NULL
  ppt <- if (!is.null(meta)) meta$design$perceivers_per_target else NULL
  if (!is.null(tpp) && max(tpp, na.rm = TRUE) < 2) varies_within$perceiver <- setdiff(varies_within$perceiver, c("target", "dyad"))
  if (!is.null(ppt) && max(ppt, na.rm = TRUE) < 2) varies_within$target <- setdiff(varies_within$target, c("perceiver", "dyad"))
  ## pairs design: perceiver, target and dyad are the same unit; keep one
  is_pairs <- !is.null(tpp) && !is.null(ppt) && max(tpp, na.rm = TRUE) < 2 && max(ppt, na.rm = TRUE) < 2
  if (is_pairs && all(c("perceiver", "target") %in% random_units) && is.null(random_structure)) {
    person <- intersect(random_units, c("perceiver", "target"))
    message("Pairs design: perceivers, targets and dyads are the same units. Person-level random effects are kept for `",
            person[1], "` only.")
    random_units <- setdiff(random_units, person[-1])
  }
  unit_var <- c(perceiver = "perceiver", target = "target", dyad = "perceiver:target", item = "item", group = "group",
                target_item = "target:item", perceiver_item = "perceiver:item")

  if (is.null(random_dyad)) {
    random_dyad <- if (!is.null(meta)) isTRUE(meta$design$dyad_re_identifiable) else TRUE
    if (!random_dyad && "dyad" %in% random_units && !is_pairs)
      message("Dyadic random effects omitted: not separable from perceiver/target effects in this design (",
              meta$design$type, ").")
  }
  if (!random_dyad) random_units <- setdiff(random_units, "dyad")
  ## target x item intercepts need several perceivers per target (otherwise they are the residual)
  if ("target_item" %in% random_units) {
    ppt <- if (!is.null(meta)) meta$design$perceivers_per_target else NULL
    if (is.null(ppt) && !is.null(data) && all(c("perceiver", "target") %in% cols))
      ppt <- tapply(data$perceiver, data$target, function(x) length(unique(x)))
    if (!is.null(ppt) && max(ppt, na.rm = TRUE) < 2) {
      message("Target x item random intercepts omitted: every target is judged by a single perceiver, so they are not separable from the residual.")
      random_units <- setdiff(random_units, "target_item")
    }
  }
  ## perceiver x item effects need several targets per perceiver (otherwise they are the residual)
  if ("perceiver_item" %in% random_units) {
    tpp <- if (!is.null(meta)) meta$design$targets_per_perceiver else NULL
    if (is.null(tpp) && !is.null(data) && all(c("perceiver", "target") %in% cols))
      tpp <- tapply(data$target, data$perceiver, function(x) length(unique(x)))
    if (!is.null(tpp) && max(tpp, na.rm = TRUE) < 2) {
      message("Perceiver x item random effects omitted: every perceiver judges a single target, so they are not separable from the residual.")
      random_units <- setdiff(random_units, "perceiver_item")
    }
  }

  slope_cands <- c(norm, dist, prof_dist, random_slopes)
  random  <- list()
  dropped <- list()
  if (!is.null(random_structure)) {
    if (inherits(random_structure, "mpa_formula") && missing(correlated)) correlated <- random_structure$correlated
    tmpl <- if (inherits(random_structure, "mpa_formula")) random_structure$random else random_structure
    if (!is.list(tmpl) || is.null(names(tmpl))) stop("`random_structure` must be an mpa_formula or a named list.")
    ## predictors whose slopes the template already decided on (present in it,
    ## or pruned from it) are not re-added: the structure is carried over as is
    decided <- if (inherits(random_structure, "mpa_formula")) random_structure$decided else character(0)
    new_slopes <- setdiff(c(prof_dist, random_slopes), c(unlist(tmpl), decided))
    for (u in names(tmpl)) {
      if (!u %in% names(unit_var)) stop("Unknown unit `", u, "` in `random_structure`.")
      add <- new_slopes[vapply(new_slopes, function(v) level_of(v) %in% varies_within[[u]], logical(1))]
      random[[u]] <- unique(c(tmpl[[u]], add))
      dropped[[u]] <- setdiff(new_slopes, add)
    }
    random_units <- names(random)
    random_dyad <- "dyad" %in% names(random)
  }
  ## one rule for every unit: an intercept, plus a slope for every term that varies within the unit
  for (u in if (is.null(random_structure)) random_units else character(0)) {
    keep <- slope_cands[vapply(slope_cands, function(v) level_of(v) %in% varies_within[[u]], logical(1))]
    if (!u %in% c("group", "target_item", "perceiver_item")) dropped[[u]] <- setdiff(slope_cands, keep)
    random[[u]] <- c("1", keep)
  }
  dropped <- dropped[vapply(dropped, length, integer(1)) > 0]

  out <- list(y = y, criterion = criterion, norm = norm, dist = dist,
              profiles = profiles, profile_dist = prof_dist, profile_level = prof_level,
              predictors = pred_tab, fixed_terms = fixed_terms, random = random,
              unit_var = unit_var, correlated = correlated, interactions = interactions,
              random_dyad = random_dyad, dropped = dropped,
              decided = unique(c(if (!is.null(random_structure) && inherits(random_structure, "mpa_formula")) random_structure$decided,
                                 prof_dist, random_slopes)))
  out$formula <- .mpa_build_formula(out)
  class(out) <- "mpa_formula"
  out
}

## internal unit names: "target:item" / "perceiver:item" (user-facing) are `target_item` / `perceiver_item` inside the object
.mpa_unit_key <- function(u) { u[u == "target:item"] <- "target_item"; u[u == "perceiver:item"] <- "perceiver_item"; u }

.mpa_build_formula <- function(x) {
  re <- vapply(names(x$random), function(u) {
    terms <- x$random[[u]]
    bar <- if (isTRUE(x$correlated) || length(terms) == 1) "|" else "||"
    if (!"1" %in% terms) terms <- c("0", terms)
    paste0("(", paste(terms, collapse = " + "), " ", bar, " ", x$unit_var[[u]], ")")
  }, character(1))
  rhs <- c(if (length(x$fixed_terms)) x$fixed_terms else "1", re)
  stats::as.formula(paste(x$y, "~", paste(rhs, collapse = " + ")), env = globalenv())
}

#' Deparse an mpa_formula to a single string
#' @param x An `mpa_formula`.
#' @export
mpa_formula_string <- function(x) paste(deparse(x$formula, width.cutoff = 500L), collapse = " ")

#' @export
print.mpa_formula <- function(x, ...) {
  cat("<mpa_formula>\n")
  cat(strwrap(mpa_formula_string(x), width = 90, exdent = 4), sep = "\n")
  cat("\nJudgments: ", x$y,
      if (!is.null(x$criterion)) paste0("   Criterion: ", x$criterion, " (normative: ", x$norm, ", distinctive: ", x$dist, ")")
      else "   (intercept-only variance decomposition)", "\n", sep = "")
  if (length(x$profiles))
    cat("Further predictor profiles: ",
        paste0(x$profiles, " (distinctive: ", x$profile_dist, ", ", sub("_item$", " x item", x$profile_level), ")", collapse = "; "),
        "\n", sep = "")
  if (!is.null(x$predictors)) {
    cat("Predictors:\n")
    print(x$predictors, row.names = FALSE)
  }
  cat("Interactions: ", x$interactions, "\n", sep = "")
  cat("Random effects", if (!x$correlated) " (uncorrelated)", ":\n", sep = "")
  for (u in names(x$random))
    cat("  ", format(x$unit_var[[u]], width = 18), paste(x$random[[u]], collapse = " + "), "\n", sep = "")
  if (length(x$dropped)) {
    cat("Random slopes dropped (not identifiable at that unit):\n")
    for (u in names(x$dropped)) cat("  ", format(x$unit_var[[u]], width = 18), paste(x$dropped[[u]], collapse = ", "), "\n", sep = "")
  }
  invisible(x)
}

#' Modify the random-effects part of an mpa_formula
#'
#' @param x An `mpa_formula`.
#' @param drop Named list: for each unit (`perceiver`, `target`, `dyad`, `item`,
#'   `group`, `target_item`, `perceiver_item`) a character vector of random terms to remove
#'   (`"1"` removes the intercept; removing all terms removes the unit).
#' @param add Named list of random terms to add per unit, e.g.
#'   `list(target_item = "1")` for a random intercept of the target x item
#'   cells (see [mpa_formula()]).
#' @param correlated Optionally change the correlation setting.
#' @return The modified `mpa_formula`.
#' @export
mpa_update <- function(x, drop = NULL, add = NULL, correlated = NULL) {
  stopifnot(inherits(x, "mpa_formula"))
  if (!is.null(drop)) names(drop) <- .mpa_unit_key(names(drop))
  if (!is.null(add)) names(add) <- .mpa_unit_key(names(add))
  if (is.null(x$unit_var["target_item"]) || is.na(x$unit_var["target_item"])) x$unit_var <- c(x$unit_var, target_item = "target:item")  # objects from earlier versions
  if (is.null(x$unit_var["perceiver_item"]) || is.na(x$unit_var["perceiver_item"])) x$unit_var <- c(x$unit_var, perceiver_item = "perceiver:item")
  for (u in names(drop)) {
    if (!u %in% names(x$random)) next
    x$random[[u]] <- setdiff(x$random[[u]], drop[[u]])
    if (!length(x$random[[u]])) x$random[[u]] <- NULL
  }
  for (u in names(add)) {
    if (!u %in% names(x$unit_var)) stop("Unknown unit `", u, "`.")
    x$random[[u]] <- unique(c(x$random[[u]], add[[u]]))
  }
  if (!is.null(correlated)) x$correlated <- correlated
  x$formula <- .mpa_build_formula(x)
  x
}
