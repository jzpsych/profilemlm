#' Standardized effect sizes for main and moderating effects
#'
#' Expresses the fixed effects of the predictors of a model as standardized
#' effect sizes *d*: the change in the judgments (main effects) or in a slope
#' such as distinctive agreement (moderating effects) across two standard
#' deviations of the predictor, divided by the standard deviation of that
#' intercept or slope across the units on which the predictor is measured.
#'
#' The numerator follows Gelman (2008). A numeric predictor is scaled by two
#' of its standard deviations, computed across the units of its level (one
#' value per perceiver, target, dyad or item). A predictor with two values is
#' left as it is, so that its coefficient is the difference between the two
#' groups. Both then compare a typical "high" with a typical "low" group of
#' equal size.
#'
#' The denominator is the standard deviation of the intercept or slope that
#' the predictor explains. With `denominator = "level"` (default) it is the
#' random-effect standard deviation at the level of the predictor: across
#' perceivers for a perceiver variable, across targets for a target variable,
#' across dyads for the relationship component of a decomposed dyadic
#' variable, and across items for an item property. For a dyadic variable
#' that was not decomposed it is the square root of the summed perceiver,
#' target and dyad variances, because such a variable varies at all three
#' levels. In a one-with-many design, in which every perceiver judges one
#' target, perceivers are nested in targets. A perceiver variable then varies
#' between and within targets, and the denominator is the square root of the
#' summed perceiver and target variances (as in Wessels et al., 2020). In a
#' pairs design there is one person-level unit. The model contains the predictor, so this is the *residual*
#' standard deviation: the variation among units that share a value of the
#' predictor. *d* then reads like Cohen's *d* for two equally large groups of
#' units, and it is tied to the variance that the predictor explains at its
#' level (see `prop_reduced` in [mpa_compare()]) by
#' \eqn{d = 2r / \sqrt{1 - r^2}}, where \eqn{r^2} is the proportion of
#' explained variance. This is the convention of studies with the Social
#' Accuracy Model (e.g., Human & Biesanz, 2011; Wessels et al., 2020;
#' Zimmermann et al., 2018, who define the denominator for a dyadic
#' predictor as the square root of the summed perceiver, target and dyad
#' variances of the slope).
#'
#' With `denominator = "total"` the denominator is the standard deviation of
#' the slope across all judgments, that is, the square root of the summed
#' slope variances of all units (for main effects: all intercept variances
#' plus the residual variance; cf. Westfall et al., 2014). This *d* is smaller
#' and does not depend on the level of the predictor.
#'
#' *d* is `NA` when the variance component in the denominator is absent from
#' the model or estimated as (near) zero, as after [mpa_prune()]. A variance
#' component that is estimated too small makes *d* too large, so *d* should be
#' read together with the precision of the variance estimate. The interval is
#' the interval of the fixed effect on the scale of *d*; it treats the
#' denominator as known.
#'
#' Effect sizes are given for terms with one predictor and at most one profile
#' component. Quadratic terms and interactions among predictors are left out.
#'
#' @param fit An `mpa_fit` whose formula is an `mpa_formula` with predictors.
#' @param denominator `"level"` (default) or `"total"`, see Details.
#' @param tol A variance component counts as zero if its standard deviation
#'   is below `tol` times the residual standard deviation (on the scale used
#'   for estimation).
#' @return A data frame of class `mpa_effect_size` with one row per term:
#'   `term`, `level` (of the predictor), `change` (the fixed effect times the
#'   scaling of the predictor: two standard deviations, or the difference
#'   between the two values), `sd_ref` (the denominator), `reference` (the
#'   unit or units whose variance forms the denominator), `d`, `ci_lower` and
#'   `ci_upper`.
#' @references Gelman, A. (2008). Scaling regression inputs by dividing by two
#'   standard deviations. *Statistics in Medicine, 27*(15), 2865-2873.
#'   \doi{10.1002/sim.3107}
#'
#'   Westfall, J., Kenny, D. A., & Judd, C. M. (2014). Statistical power and
#'   optimal design in experiments in which samples of participants respond
#'   to samples of stimuli. *Journal of Experimental Psychology: General,
#'   143*(5), 2020-2045. \doi{10.1037/xge0000014}
#'
#'   Zimmermann, J., Schindler, S., Klaus, G., & Leising, D. (2018). The
#'   effect of dislike on accuracy and bias in person perception. *Social
#'   Psychological and Personality Science, 9*(1), 80-88.
#'   \doi{10.1177/1948550617703167}
#' @examples
#' dat <- mpa_simulate(n_perceiver = 24, n_item = 12, group_size = 4, seed = 1,
#'                   fixed = c(intercept = 51, norm = 0.7, dist = 0.25, training_dist = 0.1),
#'                   rating_scale = c(1, 101))
#' f <- mpa_formula(dat, criterion = "self", predictors = "training", interactions = "criterion",
#'                  random_units = c("perceiver", "target", "item"), correlated = FALSE)
#' m <- mpa_fit(f, dat)
#' mpa_effect_size(m)
#' mpa_effect_size(m, denominator = "total")
#' @export
mpa_effect_size <- function(fit, denominator = c("level", "total"), tol = 1e-3) {
  stopifnot(inherits(fit, "mpa_fit"))
  if (!inherits(fit$formula, "mpa_formula")) stop("The fit must be based on an mpa_formula.")
  denominator <- match.arg(denominator)
  f <- fit$formula
  preds <- f$predictors
  if (is.null(preds)) stop("The model has no predictors; there are no main or moderating effects to standardize.")
  dat <- fit$data[mpa_rows(fit), , drop = FALSE]
  crit <- c(f$norm, f$dist, f$profile_dist)
  vt <- as.data.frame(fit$variance)
  ts <- if (is.numeric(fit$criterion_scale)) fit$criterion_scale else 1
  sigma <- vt$sd[vt$unit == "Residual"]
  vt_re <- vt[vt$unit != "Residual", , drop = FALSE]

  ## a variance component, NA if absent or (near) zero
  comp <- function(unit, effect) {
    v <- vt_re$variance[vt_re$unit == unit & vt_re$effect == effect]
    if (!length(v)) return(NA_real_)
    k <- if (effect == "(Intercept)") 0 else 1
    if (sqrt(v) * ts^k / sigma < tol) NA_real_ else v
  }
  unit_of <- c(perceiver = "perceiver", target = "target", dyad = "perceiver:target", item = "item", group = "group")
  ## pairs design: perceiver, target and dyad are one unit, kept under one name
  des <- fit$meta$design
  is_pairs <- !is.null(des) && max(des$targets_per_perceiver) < 2 && max(des$perceivers_per_target) < 2
  if (is_pairs) {
    person <- intersect(c("perceiver", "target"), vt_re$unit)[1]
    unit_of[c("perceiver", "target", "dyad")] <- person
  }
  ## one-with-many: the units of the many side are nested in those of the one side
  nested <- character(0)
  if (!is.null(des) && !is_pairs) {
    if (max(des$targets_per_perceiver) < 2) nested <- "perceiver"
    if (max(des$perceivers_per_target) < 2) nested <- "target"
  }
  decomposed_r <- paste0(fit$meta$vars$decomposed, ".r")
  unit_id <- list(perceiver = dat$perceiver, target = dat$target, dyad = paste(dat$perceiver, dat$target),
                  item = dat$item, group = if ("group" %in% names(dat)) dat$group else NULL)

  rows <- list()
  for (j in seq_len(nrow(fit$fixed))) {
    term <- fit$fixed$term[j]
    parts <- strsplit(term, ":", fixed = TRUE)[[1]]
    pr <- parts[!parts %in% crit]
    cr <- parts[parts %in% crit]
    if (length(pr) != 1 || length(cr) > 1 || !pr %in% preds$term || grepl("^I\\(", pr)) next
    level <- preds$level[match(pr, preds$term)]
    x <- dat[[pr]]
    if (!is.numeric(x)) next
    ## one value per unit of the predictor's level
    xu <- if (level %in% names(unit_id) && !is.null(unit_id[[level]])) x[!duplicated(unit_id[[level]])] else x
    xu <- xu[!is.na(xu)]
    scale <- if (length(unique(xu)) == 2) diff(range(xu)) else 2 * stats::sd(xu)
    effect <- if (length(cr)) cr else "(Intercept)"

    if (denominator == "level") {
      units <- if ((level == "dyad" && !pr %in% decomposed_r && !is_pairs) || level %in% nested)
        intersect(c("perceiver", "target", "perceiver:target"), vt_re$unit)
      else if (level %in% names(unit_of)) unit_of[[level]] else character(0)
      v <- vapply(units, comp, numeric(1), effect = effect)
      ## a single unit needs its component; a variable that varies at several levels sums what is there
      v_ref <- if (!length(v) || all(is.na(v))) NA_real_ else if (length(v) == 1) v else sum(v, na.rm = TRUE)
      reference <- paste(units[!is.na(v)], collapse = " + ")
    } else {
      units <- unique(vt_re$unit[vt_re$effect == effect])
      v <- vapply(units, comp, numeric(1), effect = effect)
      v_ref <- if (!length(v) || all(is.na(v))) NA_real_ else sum(v, na.rm = TRUE)
      if (effect == "(Intercept)" && !is.na(v_ref)) v_ref <- v_ref + sigma^2
      reference <- paste(c(units[!is.na(v)], if (effect == "(Intercept)") "Residual"), collapse = " + ")
    }
    sd_ref <- sqrt(v_ref)
    fx <- fit$fixed[j, ]
    rows[[length(rows) + 1]] <- data.frame(
      term = term, level = level, change = fx$estimate * scale, sd_ref = sd_ref,
      reference = if (is.na(sd_ref)) NA_character_ else reference,
      d = fx$estimate * scale / sd_ref, ci_lower = fx$ci_lower * scale / sd_ref, ci_upper = fx$ci_upper * scale / sd_ref,
      stringsAsFactors = FALSE)
  }
  if (!length(rows)) stop("No term with a single predictor found; see ?mpa_effect_size.")
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  attr(out, "denominator") <- denominator
  attr(out, "level") <- if (is.null(fit$level)) 0.95 else fit$level
  class(out) <- c("mpa_effect_size", "data.frame")
  out
}

#' @export
print.mpa_effect_size <- function(x, digits = 2, ...) {
  cat("<mpa_effect_size>  d = change across 2 SD of the predictor (two values: their difference) / SD of the intercept or slope\n")
  cat("Denominator: ", if (identical(attr(x, "denominator"), "total")) "SD across all judgments (all units)"
                       else "random-effect SD at the level of the predictor", "; ",
      round(100 * (attr(x, "level") %||% 0.95)), "% interval of the fixed effect, denominator treated as known\n", sep = "")
  y <- as.data.frame(x)
  y$change <- signif(x$change, 3); y$sd_ref <- signif(x$sd_ref, 3)
  y$d <- round(x$d, digits); y$ci_lower <- round(x$ci_lower, digits); y$ci_upper <- round(x$ci_upper, digits)
  print(y, row.names = FALSE, ...)
  if (anyNA(x$d)) cat("NA: the variance component in the denominator is not in the model or was estimated as zero.\n")
  invisible(x)
}
