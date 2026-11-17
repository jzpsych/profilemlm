#' Plot how a slope changes with a moderator
#'
#' Shows the estimated effect (slope) of one predictor -- typically the
#' distinctive or normative component -- as a function of a moderator, with a
#' confidence (or credible) band, optionally for three levels (minimum,
#' median, maximum) of a second moderator. Requires the package `ggplot2`.
#'
#' @param fit An `mpa_fit`.
#' @param slope Name of the predictor whose slope is plotted on the y-axis
#'   (e.g., `"combined.dist"`).
#' @param x Name of the moderator on the x-axis.
#' @param by Optional second moderator shown as separate lines at its minimum,
#'   median and maximum (or at `by_values`).
#' @param by_values Optional numeric values of `by` at which lines are drawn.
#' @param n Number of points along `x`.
#' @param level Level of the band (default: the level of the fit's intervals).
#' @param labels Optional named character vector mapping variable names to
#'   axis labels, e.g. `c(LIKE.n = "Liking", combined.dist = "Distinctive agreement")`.
#' @param hold Value at which all *other* predictors of the model (those in
#'   `fit$formula$predictors`, other than `x`, `by` and `slope`) are held.
#'   Default `0`, i.e., the centering point (scale midpoint / grand mean), so
#'   that the plotted slopes correspond to the fixed-effect coefficients.
#'   `NULL` holds them at their sample means. Note that for a
#'   0/1 dummy such as an informant type, `0` is the reference category, not
#'   the average; use `hold = NULL` or plot `by` that variable if the average
#'   is wanted.
#' @details The slope at each point of the grid is a linear combination of
#'   the fixed effects, obtained as the derivative of the model's fixed part
#'   with respect to `slope` (by central differences, which are exact for the
#'   linear and quadratic terms the package builds). For `lme4` fits the band
#'   is a Wald interval (asymptotic, as `emmeans::emtrends(..., lmer.df =
#'   "asymp")` would give), for `brms` fits the equal-tailed credible interval
#'   of the posterior of the linear combination. Variables of the fixed part
#'   that are neither `x`, `by`, `slope` nor held at `hold` (e.g. the other
#'   criterion component) are set to their means over the rows used by the fit;
#'   they do not affect the plotted slope unless they interact with `slope`.
#' @return A `ggplot` object. The underlying estimates are in `attr(p, "data")`.
#' @export
mpa_plot <- function(fit, slope, x, by = NULL, by_values = NULL, n = 50, level = NULL, labels = NULL, hold = 0) {
  stopifnot(inherits(fit, "mpa_fit"))
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Package `ggplot2` is required for mpa_plot().")
  if (is.null(level)) level <- if (is.null(fit$level)) 0.95 else fit$level
  dat <- fit$data[mpa_rows(fit), , drop = FALSE]
  in_model <- unique(unlist(strsplit(fit$fixed$term, ":", fixed = TRUE)))
  in_model <- unique(c(in_model, sub("^I\\((.*)\\^2\\)$", "\\1", in_model)))
  for (v in c(slope, x, by)) {
    if (!v %in% names(dat)) stop("Variable `", v, "` not found in the data of the fit.")
    if (!v %in% in_model) stop("Variable `", v, "` is not a predictor in this model.")
  }
  lab <- function(v) if (!is.null(labels) && v %in% names(labels)) labels[[v]] else v

  at <- list(seq(min(dat[[x]], na.rm = TRUE), max(dat[[x]], na.rm = TRUE), length.out = n))
  names(at) <- x
  if (!is.null(by))
    at[[by]] <- if (is.null(by_values)) unname(stats::quantile(dat[[by]], c(0, 0.5, 1), na.rm = TRUE)) else by_values
  if (!is.null(hold) && inherits(fit$formula, "mpa_formula") && !is.null(fit$formula$predictors)) {
    others <- fit$formula$predictors$term
    others <- others[!grepl("^I\\(", others)]
    others <- setdiff(others, c(x, by, slope))
    for (o in others) at[[o]] <- hold
  }
  emm <- .mpa_slopes(fit, slope, at, level)
  ## slopes of a scaled criterion component are reported per criterion point
  ts <- if (isTRUE(fit$criterion_scale != 1) && slope %in% fit$criterion_vars) fit$criterion_scale else 1
  names(emm)[names(emm) == x] <- ".x"
  if (!is.null(by)) names(emm)[names(emm) == by] <- ".by"
  emm$.y <- emm$estimate / ts; emm$.lcl <- emm$lower / ts; emm$.ucl <- emm$upper / ts

  p <- ggplot2::ggplot(emm, ggplot2::aes(x = .x, y = .y)) +
    ggplot2::geom_hline(yintercept = 0, linetype = 2) +
    ggplot2::theme_bw() +
    ggplot2::labs(x = lab(x), y = if (!is.null(labels) && slope %in% names(labels)) lab(slope) else paste("Slope of", slope))
  if (is.null(by)) {
    p <- p + ggplot2::geom_ribbon(ggplot2::aes(ymin = .lcl, ymax = .ucl), alpha = 0.1, colour = NA) +
      ggplot2::geom_line()
  } else {
    p <- p + ggplot2::geom_ribbon(ggplot2::aes(ymin = .lcl, ymax = .ucl, group = .by, fill = .by),
                                  alpha = 0.1, colour = NA) +
      ggplot2::geom_line(ggplot2::aes(group = .by, colour = .by)) +
      ggplot2::labs(colour = lab(by), fill = lab(by))
  }
  attr(p, "data") <- emm
  attr(p, "slope") <- slope
  p
}

## The slope of `slope` on a grid of the variables in `at` (a named list of
## values, expanded to all combinations): estimate and interval of the linear
## combination L %*% beta, where L is the derivative of the fixed-effects
## design matrix with respect to `slope`, by central differences. All other
## variables of the fixed part are set to their means over the rows used.
.mpa_slopes <- function(fit, slope, at, level) {
  dat <- fit$data[mpa_rows(fit), , drop = FALSE]
  f <- if (inherits(fit$formula, "mpa_formula")) fit$formula$formula else stats::as.formula(fit$formula)
  tt <- stats::delete.response(stats::terms(.mpa_nobars(f)))
  vars <- all.vars(tt)
  grid <- expand.grid(at, KEEP.OUT.ATTRS = FALSE)
  for (v in setdiff(vars, names(grid))) {
    if (!is.numeric(dat[[v]])) stop("mpa_plot() supports numeric predictors only; `", v, "` is not numeric.")
    grid[[v]] <- mean(dat[[v]])
  }
  h <- 1e-3 * max(stats::sd(dat[[slope]]), 1e-8)
  up <- grid; up[[slope]] <- grid[[slope]] + h
  lo <- grid; lo[[slope]] <- grid[[slope]] - h
  L <- (stats::model.matrix(tt, up) - stats::model.matrix(tt, lo)) / (2 * h)
  miss <- setdiff(fit$fixed$term, colnames(L))
  if (length(miss)) stop("Internal error: fixed effects not found in the design matrix: ", paste(miss, collapse = ", "))
  res <- .engine_lincomb(.mpa_fit_engine(fit), fit$fit, L[, fit$fixed$term, drop = FALSE], level)
  cbind(grid[names(at)], res)
}

#' Give plots of the same slope a common y-axis
#'
#' Plots produced by [mpa_plot()] carry their estimates as attributes; this
#' function sets identical y-limits for all plots that show the same slope
#' (e.g., all panels of distinctive agreement), so that panels in a row can be
#' compared. Plots of different slopes keep separate limits.
#'
#' @param plots Named list of plots from [mpa_plot()] / [mpa_plots()].
#' @param expand Proportional padding added to the common range.
#' @return The list with adjusted plots.
#' @export
mpa_align_y <- function(plots, expand = 0.02) {
  slopes <- vapply(plots, function(p) { s <- attr(p, "slope"); if (is.null(s)) NA_character_ else s }, character(1))
  for (s in unique(stats::na.omit(slopes))) {
    idx <- which(slopes == s)
    rng <- range(unlist(lapply(plots[idx], function(p) { d <- attr(p, "data"); c(d$.lcl, d$.ucl, 0) })), na.rm = TRUE)
    pad <- diff(rng) * expand
    for (i in idx) {
      a1 <- attr(plots[[i]], "data"); a2 <- attr(plots[[i]], "slope")
      plots[[i]] <- plots[[i]] + ggplot2::coord_cartesian(ylim = rng + c(-pad, pad))
      attr(plots[[i]], "data") <- a1; attr(plots[[i]], "slope") <- a2
    }
  }
  plots
}

#' All marginal-effect plots implied by a fitted MPA model
#'
#' Convenience wrapper around [mpa_plot()] that creates, for a model built with
#' [mpa_formula()], (a) plots of the distinctive and normative slopes, and of
#' the slopes of further predictor profiles, as functions of every additional
#' predictor (if the formula contains criterion interactions), (b) plots of each predictor's own effect as a function of
#' every predictor from another level, and (c) for `interactions = "all"`,
#' plots of the criterion slopes against one predictor at three levels of another.
#'
#' @param fit An `mpa_fit` whose `formula` is an `mpa_formula`.
#' @param labels Optional named character vector of axis labels, see [mpa_plot()].
#' @param common_y Give all panels showing the same slope a common y-axis
#'   (default `TRUE`; see [mpa_align_y()]).
#' @param ... Passed to [mpa_plot()].
#' @return A named list of `ggplot` objects. Names have the form
#'   `"<slope>_by_<x>"` or `"<slope>_by_<x>_at_<by>"`. Combine them with, e.g.,
#'   `patchwork::wrap_plots()`.
#' @export
mpa_plots <- function(fit, labels = NULL, common_y = TRUE, ...) {
  stopifnot(inherits(fit, "mpa_fit"), inherits(fit$formula, "mpa_formula"))
  f <- fit$formula
  preds <- f$predictors
  if (is.null(preds)) stop("The model has no additional predictors; nothing to plot.")
  preds <- preds[!grepl("^I\\(", preds$term), ]     # quadratic terms are absorbed by emtrends via the base variable
  out <- list()
  slopes <- c(f$dist, f$profile_dist, f$norm)
  slope_lab <- c(if (!is.null(f$dist)) stats::setNames("Distinctive agreement", f$dist),
                 if (length(f$profile_dist)) stats::setNames(paste0("Slope of the distinctive ", f$profiles, " profile"), f$profile_dist),
                 if (!is.null(f$norm)) stats::setNames("Normative agreement", f$norm))
  labels <- c(labels, slope_lab[setdiff(names(slope_lab), names(labels))])

  if (f$interactions %in% c("criterion", "all") && length(slopes)) {
    for (m in preds$term) for (s in slopes)
      out[[paste0(s, "_by_", m)]] <- mpa_plot(fit, slope = s, x = m, labels = labels, ...)
  }
  if (nrow(preds) > 1 && f$interactions == "all") {
    pairs <- utils::combn(seq_len(nrow(preds)), 2, simplify = FALSE)
    for (pr in pairs) {
      a <- preds$term[pr[1]]; b <- preds$term[pr[2]]
      if (preds$level[pr[1]] == preds$level[pr[2]]) next
      out[[paste0(a, "_by_", b)]] <- mpa_plot(fit, slope = a, x = b, labels = labels, ...)
      out[[paste0(b, "_by_", a)]] <- mpa_plot(fit, slope = b, x = a, labels = labels, ...)
      for (s in slopes) {
        out[[paste0(s, "_by_", b, "_at_", a)]] <- mpa_plot(fit, slope = s, x = b, by = a, labels = labels, ...)
        out[[paste0(s, "_by_", a, "_at_", b)]] <- mpa_plot(fit, slope = s, x = a, by = b, labels = labels, ...)
      }
    }
  }
  if (!length(out)) stop("No interaction terms in the model; use mpa_plot() directly for single effects.")
  if (common_y) out <- mpa_align_y(out)
  out
}

#' Save a list of plots to files
#'
#' @param plots A named list of `ggplot` objects, e.g. from [mpa_plots()], or a
#'   single plot (also a `patchwork` composition).
#' @param dir Directory (created if necessary).
#' @param prefix Optional file-name prefix.
#' @param width,height,dpi Passed to [ggplot2::ggsave()]; in inches.
#' @param device File type, e.g. `"png"` or `"pdf"`.
#' @return Invisibly, the file paths written.
#' @export
mpa_save_plots <- function(plots, dir = "figures", prefix = "", width = 6, height = 4.5, dpi = 300, device = "png") {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Package `ggplot2` is required.")
  if (inherits(plots, "ggplot")) plots <- list(plot = plots)
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
  nm <- names(plots); if (is.null(nm)) nm <- paste0("plot", seq_along(plots))
  files <- file.path(dir, paste0(prefix, gsub("[^A-Za-z0-9_.-]", "_", nm), ".", device))
  for (i in seq_along(plots)) ggplot2::ggsave(files[i], plots[[i]], width = width, height = height, dpi = dpi)
  message(length(files), " file(s) written to ", normalizePath(dir))
  invisible(files)
}

utils::globalVariables(c(".x", ".y", ".lcl", ".ucl", ".by"))
