#' Fit a multilevel profile analysis model
#'
#' Fits an [mpa_formula()] to an [mpa_data()] object with one of two
#' estimation engines, `lme4` (maximum likelihood, the default) or `brms`
#' (Bayesian, Stan), and returns tidy tables of fixed effects and variance
#' components whose layout does not depend on the engine.
#'
#' @section Engines:
#' The engine is chosen with `engine` (`"lme4"`, `"brms"`, or an engine
#' object from [engine_lme4()] / [engine_brms()]) or, for a whole session,
#' with `options(profilemlm.engine = "brms")`. Settings of the engine
#' (`optimizer`, `REML`, `check_derivs`, `chains`, `prior`, ...) can be given
#' directly as further arguments of `mpa_fit()`; see [mpa_engine] for the
#' complete list. Both engines fit the same model to the same rows. The tables
#' of the fit have the same columns `term`, `estimate`, `se`, `ci_lower`,
#' `ci_upper` (fixed effects) and `unit`, `effect`, `variance`, `sd` (variance
#' components), and differ in what these contain and in the extra columns:
#'
#' * `lme4`: ML or REML estimates, standard errors, Wald intervals, `z` and
#'   `p` (or `df`, `t`, `p` with `df = "satterthwaite"`); `logLik`, `AIC`,
#'   `BIC`; convergence messages and a singularity flag.
#' * `brms`: posterior medians, posterior standard deviations (in `se`),
#'   equal-tailed credible intervals, `rhat`, `ess_bulk`, `ess_tail`; the
#'   variance table additionally has `sd_lower`, `sd_upper`; no likelihood,
#'   no singularity (a posterior has no boundary solution; `singular` is
#'   always `FALSE` and [mpa_prune()] leaves the structure unchanged).
#'
#' @param formula An `mpa_formula` object (or a plain `lme4`-style formula).
#' @param data An `mpa_data` object or data frame.
#' @param engine `"lme4"` (default), `"brms"`, or an engine object; see
#'   [mpa_engine]. An `mpa_fit` is accepted too and means "the engine and
#'   settings of that fit". The default is `getOption("profilemlm.engine", "lme4")`.
#' @param start_from Optional `mpa_fit` whose variance-parameter estimates are
#'   used as starting values, typically the previous step's model (`lme4`
#'   only; ignored by `brms` with a message). Random terms
#'   that exist identically in both models (same unit, same effects) take over
#'   the previous estimates; terms that are new (e.g., the random slope of a
#'   predictor added in this step) start from lme4's defaults. A message says
#'   how many terms were matched.
#' @param report `"criterion"` (default): coefficients, standard errors, confidence
#'   intervals and variance components that involve the criterion components are
#'   reported per *criterion point*, i.e., back-transformed if the data were built
#'   with a `criterion_scale` other than 1 (each term is divided by `criterion_scale^k`,
#'   where `k` is the number of criterion variables in the term, counting the
#'   distinctive components of further predictor profiles; slope SDs by
#'   `criterion_scale`). `"analysis"`: report on the scale used for estimation.
#' @param same_rows_as Optional `mpa_fit`; the model is then fitted to exactly
#'   the observations used by that fit. Use this to fit a reference model
#'   (e.g., the Step 0 variance decomposition) on the same rows as a model
#'   whose predictors have missing values, so that [mpa_compare()] compares
#'   like with like.
#' @param level Level of the confidence / credible intervals (default 0.95).
#' @param ... Engine settings, e.g. `REML = TRUE`, `optimizer = "nloptwrap"`,
#'   `check_derivs = FALSE`, `df = "satterthwaite"` for `lme4`, or
#'   `chains = 4`, `prior = ...` for `brms`; see [mpa_engine].
#'
#' @return An object of class `mpa_fit` (and `mpa_fit_lme4` or
#'   `mpa_fit_brms`) with elements `engine` (`"lme4"` or `"brms"`),
#'   `settings` (the engine object), `fit` (the native fit: `merMod` or
#'   `brmsfit`), `formula`, `data` (the data frame passed to `mpa_fit()`),
#'   `rows` (the rows of `data` that were used), `fixed`, `variance`
#'   (incl. residual) and `correlations` (data frames, see Engines),
#'   `diagnostics` (list: `converged`, `singular`, `messages`, `details`),
#'   `fit_stats` (`logLik`, `df`, `AIC`, `BIC` for `lme4`; sampling
#'   information for `brms`), `n_obs`, `n_units`, `level`, and, for
#'   convenience, `converged`, `singular`, `conv_messages`, `logLik`, `AIC`,
#'   `BIC`, `REML`, `df`.
#' @examples
#' dat <- mpa_simulate(n_perceiver = 16, n_item = 10, group_size = 4, seed = 1, rating_scale = c(1, 101))
#' f <- mpa_formula(dat, criterion = "self", random_units = c("perceiver", "target"),
#'                  correlated = FALSE)
#' m <- mpa_fit(f, dat)
#' m
#' m$fixed
#' \dontrun{
#' mpa_fit(f, dat, engine = "brms", chains = 2)      # the same model, Bayesian
#' }
#' @export
mpa_fit <- function(formula, data, engine = getOption("profilemlm.engine", "lme4"),
                    start_from = NULL, report = c("criterion", "analysis"), same_rows_as = NULL,
                    level = 0.95, ...) {
  report <- match.arg(report)
  engine <- .mpa_resolve_engine(engine, list(...))
  .engine_check(engine)
  stopifnot(is.numeric(level), length(level) == 1, level > 0, level < 1)
  f <- if (inherits(formula, "mpa_formula")) formula$formula else stats::as.formula(formula)
  dat <- as.data.frame(data)
  rownames(dat) <- NULL

  ## the rows used: complete cases on the model variables, or those of another fit
  vars <- intersect(all.vars(f), names(dat))
  rows <- which(stats::complete.cases(dat[vars]))
  if (!is.null(same_rows_as)) {
    stopifnot(inherits(same_rows_as, "mpa_fit"))
    rows <- intersect(mpa_rows(same_rows_as), rows)
  }
  dat_used <- dat[rows, , drop = FALSE]

  start <- NULL
  if (!is.null(start_from)) {
    stopifnot(inherits(start_from, "mpa_fit"))
    start <- .engine_start(engine, f, dat_used, start_from)
  }
  native <- .engine_fit(engine, f, dat_used, start)

  meta <- attr(data, "mpa")
  ts <- if (report == "criterion" && is.numeric(meta$settings$criterion_scale)) meta$settings$criterion_scale else 1
  criterion_vars <- if (inherits(formula, "mpa_formula")) c(formula$norm, formula$dist, formula$profile_dist) else character(0)
  tabs <- .engine_tables(engine, native, level)
  diag <- .engine_diagnostics(engine, native)
  stats <- .engine_fit_stats(engine, native)
  out <- list(engine = engine$name, settings = engine,
              fit = native,
              formula = formula,
              data = dat, rows = rows,
              meta = meta,
              report = report, criterion_scale = ts, criterion_vars = criterion_vars,
              fixed = .mpa_table(.mpa_rescale_fixed(tabs$fixed, criterion_vars, ts)),
              variance = .mpa_table(.mpa_rescale_variance(tabs$variance, criterion_vars, ts)),
              correlations = .mpa_table(tabs$correlations),
              labels = tabs$labels,
              diagnostics = diag, fit_stats = stats,
              n_obs = length(rows), n_units = tabs$n_units, level = level,
              ## convenience fields (also kept for backward compatibility)
              converged = diag$converged, conv_messages = diag$messages, singular = diag$singular,
              REML = engine$REML, df = engine$df, check_derivs = engine$check_derivs,
              logLik = if (is.null(stats$logLik)) NA_real_ else stats$logLik,
              AIC = if (is.null(stats$AIC)) NA_real_ else stats$AIC,
              BIC = if (is.null(stats$BIC)) NA_real_ else stats$BIC)
  class(out) <- c(paste0("mpa_fit_", engine$name), "mpa_fit")
  out
}

#' Row indices (in the data passed to mpa_fit) of the observations used by a fit
#' @param x An `mpa_fit`.
#' @return Integer vector of row indices into `x$data`.
#' @export
mpa_rows <- function(x) {
  stopifnot(inherits(x, "mpa_fit"))
  if (!is.null(x$rows)) return(x$rows)
  as.integer(rownames(stats::model.frame(x$fit)))     # fits saved by profilemlm <= 0.3.0
}

## number of criterion variables in each term (for back-transformation of scaled criterion components)
.mpa_criterion_power <- function(terms, criterion_vars) {
  if (!length(criterion_vars)) return(rep(0L, length(terms)))
  vapply(strsplit(terms, ":", fixed = TRUE), function(z) sum(z %in% criterion_vars), integer(1))
}

.mpa_rescale_fixed <- function(tab, criterion_vars, ts) {
  if (ts == 1) return(tab)
  f <- ts^.mpa_criterion_power(tab$term, criterion_vars)
  for (v in intersect(c("estimate", "se", "ci_lower", "ci_upper"), names(tab))) tab[[v]] <- tab[[v]] / f
  tab
}

.mpa_rescale_variance <- function(tab, criterion_vars, ts) {
  if (ts == 1) return(tab)
  k <- .mpa_criterion_power(tab$effect, criterion_vars)
  for (v in intersect(c("sd", "sd_lower", "sd_upper"), names(tab))) tab[[v]] <- tab[[v]] / ts^k
  for (v in intersect(c("variance", "variance_lower", "variance_upper"), names(tab))) tab[[v]] <- tab[[v]] / ts^(2 * k)
  tab
}

## The tables of an mpa_fit are plain data frames with an extra class whose
## print method rounds, so that `m$fixed` reads as cleanly as `print(m)`.
.mpa_table <- function(x) {
  if (is.null(x)) return(NULL)
  class(x) <- c("mpa_table", "data.frame")
  x
}

#' @export
print.mpa_table <- function(x, digits = 3, row.names = FALSE, ...) {
  y <- as.data.frame(x)
  num <- vapply(y, is.numeric, logical(1))
  y[num] <- lapply(y[num], function(v) round(v, digits))
  if ("p" %in% names(y)) y$p <- format.pval(x$p, digits = 3, eps = 1e-4)
  for (v in intersect(c("ess_bulk", "ess_tail"), names(y))) y[[v]] <- round(x[[v]])
  print(y, row.names = row.names, ...)
  invisible(x)
}

## engine settings of a fit, also for fits saved by profilemlm <= 0.3.0 (no `settings`)
.mpa_fit_engine <- function(x) {
  if (!is.null(x$settings)) return(x$settings)
  engine_lme4(REML = isTRUE(x$REML), df = if (is.null(x$df)) "z" else x$df,
              check_derivs = !isFALSE(x$check_derivs))
}

#' @export
print.mpa_fit <- function(x, digits = 3, ...) {
  eng <- .mpa_fit_engine(x)
  cat("<mpa_fit>  ", .engine_label(eng, x$fit), ", ", format(x$n_obs, big.mark = ","), " observations; ",
      paste(paste0(names(x$n_units), ": ", x$n_units), collapse = ", "), "\n", sep = "")
  fstr <- if (inherits(x$formula, "mpa_formula")) mpa_formula_string(x$formula)
          else paste(deparse(stats::formula(x$fit), width.cutoff = 500L), collapse = " ")
  cat(strwrap(fstr, width = 90, exdent = 4), sep = "\n")
  .engine_print_diagnostics(eng, x)
  cat("\n")

  if (isTRUE(x$criterion_scale != 1))
    cat("Criterion components were divided by ", signif(x$criterion_scale, 4),
        " for estimation; results are reported per criterion point.\n", sep = "")
  lab <- x$labels
  if (is.null(lab)) lab <- list(estimate = if (isTRUE(x$REML)) "REML estimates" else "ML estimates",
                                test = if (identical(x$df, "satterthwaite")) "Satterthwaite df" else "Wald z-tests",
                                interval = "95% Wald CI")
  cat("Fixed effects (", paste(unlist(lab[c("estimate", "test", "interval")]), collapse = ", "), "):\n", sep = "")
  print(.mpa_table(x$fixed), digits = digits)

  cat("\nVariance components", if (!is.null(lab$variance)) paste0(" (", lab$variance, ")"), ":\n", sep = "")
  print(.mpa_table(x$variance), digits = digits)
  if (!is.null(x$correlations)) {
    cat("\nRandom-effect correlations:\n")
    print(.mpa_table(x$correlations), digits = 2)
  }
  invisible(x)
}

#' @export
summary.mpa_fit <- function(object, ...) {
  print(object, ...)
  invisible(object)
}

#' Compare two multilevel profile analysis models
#'
#' Reports, for each variance component, the proportion by which it is reduced
#' from a reference model to a comparison model (`prop_reduced`, a
#' "pseudo R-squared" per component), the change in information criteria,
#' and -- if both models were fitted with ML -- a likelihood-ratio test.
#' For `brms` fits the variance table is computed from the posterior medians
#' and the likelihood-based parts are omitted; see Details.
#'
#' @param fit0 Reference `mpa_fit` (e.g., the base model).
#' @param fit1 Comparison `mpa_fit` (e.g., with additional predictors).
#' @param refit If the two fits use different observations (typically because
#'   predictors added in `fit1` have missing values), `refit = TRUE` (default)
#'   refits the *reference* model on the observations used by `fit1` before
#'   comparing (with the engine and settings of `fit0`); `FALSE` compares as
#'   is (variance components are then not
#'   strictly comparable and no likelihood-ratio test is computed).
#' @param loo For `brms` fits: compare the models by leave-one-out cross
#'   validation ([loo::loo_compare()]; slow for large data and requires fits
#'   estimated with `save_group_draws = TRUE`). Default `FALSE`.
#' @return An object of class `mpa_compare` with `variance` (data frame with
#'   the columns `unit`, `effect`, `variance_0`, `variance_1` and
#'   `prop_reduced`), `lrt` (data frame or `NULL`), `delta` (AIC/BIC
#'   differences, `NA` for `brms`) and, with `loo = TRUE`, `loo` (the
#'   `loo_compare` table).
#' @details `prop_reduced` is `1 - variance_1 / variance_0`: the share of the
#'   reference model's variance component that the predictors added in the
#'   comparison model account for. For example, a `prop_reduced` of 0.68 for
#'   the item slope of the distinctive component means that 68 percent of the
#'   item-level differences in distinctive agreement are explained by the
#'   added item properties. Its upper bound is 1; it is `NA` for random terms
#'   that exist in only one of the two models.
#'
#'   Variance components can also *increase* when predictors are added (a
#'   negative `prop_reduced`). This is not an error: with ML estimation and
#'   cross-classified random effects, adding a predictor redistributes variance
#'   among components, and variance explained at one level may reappear as
#'   variance at another. Interpret `prop_reduced` as descriptive.
#'
#'   Bayesian fits have no likelihood-ratio test and no AIC/BIC. The
#'   variance table and `prop_reduced` are computed from the posterior medians
#'   of the two fits; `loo = TRUE` adds the expected log predictive density
#'   difference with its standard error (negative `elpd_diff` of the
#'   comparison model favours the reference model). Two fits must come from
#'   the same engine to be compared.
#' @export
mpa_compare <- function(fit0, fit1, refit = TRUE, loo = FALSE) {
  stopifnot(inherits(fit0, "mpa_fit"), inherits(fit1, "mpa_fit"))
  e0 <- .mpa_fit_engine(fit0); e1 <- .mpa_fit_engine(fit1)
  if (e0$name != e1$name)
    stop("The two fits were estimated with different engines (", e0$name, ", ", e1$name,
         "); refit one of them with the other engine to compare them.")
  refitted <- FALSE
  if (fit0$n_obs != fit1$n_obs && refit) {
    r0 <- mpa_rows(fit0); r1 <- mpa_rows(fit1)
    if (!all(r1 %in% r0))
      stop("The comparison model uses observations not used by the reference model; cannot refit.")
    message("Reference model refitted on the ", format(fit1$n_obs, big.mark = ","),
            " observations used by the comparison model (", format(fit0$n_obs - fit1$n_obs, big.mark = ","),
            " dropped).")
    fit0 <- mpa_fit(fit0$formula, fit0$data, engine = e0, report = if (is.null(fit0$report)) "criterion" else fit0$report,
                    same_rows_as = fit1, level = if (is.null(fit0$level)) 0.95 else fit0$level)
    refitted <- TRUE
  }
  v0 <- fit0$variance; v1 <- fit1$variance
  names(v0)[3:4] <- c("variance_0", "sd_0"); names(v1)[3:4] <- c("variance_1", "sd_1")
  vt <- merge(v0[, c("unit", "effect", "variance_0")], v1[, c("unit", "effect", "variance_1")],
              by = c("unit", "effect"), all = TRUE, sort = FALSE)
  vt$prop_reduced <- 1 - vt$variance_1 / vt$variance_0
  ## order: units as in fit1, residual last
  ord <- order(match(vt$unit, unique(c(v1$unit, v0$unit))), vt$effect != "(Intercept)")
  vt <- vt[ord, ]; rownames(vt) <- NULL

  ## the LRT presupposes nested models: the fixed effects of the reference
  ## model must be a subset of those of the comparison model (random terms may
  ## differ only by added terms, which is not checked here)
  key <- function(x) vapply(strsplit(x, ":", fixed = TRUE), function(z) paste(sort(z), collapse = ":"), character(1))
  not_nested <- setdiff(key(fit0$fixed$term), key(fit1$fixed$term))
  if (length(not_nested))
    warning("The reference model (0) has fixed effects that the comparison model (1) lacks (",
            paste(not_nested, collapse = ", "), "); the models are not nested and the likelihood-ratio test ",
            "and `prop_reduced` is not interpretable as such.")

  lrt <- NULL
  s0 <- fit0$fit_stats; s1 <- fit1$fit_stats
  has_ll <- is.numeric(s0$logLik) && is.numeric(s1$logLik)
  if (has_ll && isFALSE(e0$REML) && isFALSE(e1$REML) && fit0$n_obs == fit1$n_obs) {
    df0 <- s0$df; df1 <- s1$df
    chisq <- 2 * (s1$logLik - s0$logLik)
    ddf <- df1 - df0
    if (ddf > 0 && chisq < -1e-3)
      warning("The larger model (1) has a LOWER log-likelihood than the smaller model (0) although it has more ",
              "parameters. If the models are nested this is impossible at the optimum: the larger model has ",
              "most likely not converged (note that check_derivs = FALSE suppresses lme4's convergence warnings). ",
              "Try another optimizer (e.g. optimizer = \"nloptwrap\"), correlated = FALSE, or a simpler random structure.")
    if (ddf < 0 && chisq > 1e-3)
      warning("The smaller model (0) has a HIGHER log-likelihood than the larger model (1); check convergence of model 0.")
    if (ddf < 0) { chisq <- -chisq; ddf <- -ddf }
    lrt <- data.frame(chisq = chisq, df = ddf, p = stats::pchisq(chisq, ddf, lower.tail = FALSE))
  } else if (has_ll && fit0$n_obs != fit1$n_obs) {
    warning("Models were fitted to different numbers of observations; no likelihood-ratio test.")
  }
  delta <- if (has_ll) c(AIC = s1$AIC - s0$AIC, BIC = s1$BIC - s0$BIC) else c(AIC = NA_real_, BIC = NA_real_)
  loo_tab <- if (isTRUE(loo)) .mpa_loo_compare(fit0, fit1) else NULL
  structure(list(variance = vt, lrt = lrt, delta = delta, loo = loo_tab, engine = e0$name,
                 n_obs = c(fit0$n_obs, fit1$n_obs), refitted = refitted, fit0 = fit0),
            class = "mpa_compare")
}

#' @export
print.mpa_compare <- function(x, digits = 3, ...) {
  cat("<mpa_compare>  variance components: reference model (0) vs. comparison model (1); n = ",
      format(x$n_obs[2], big.mark = ","), if (x$refitted) " (reference model refitted on these rows)", "\n", sep = "")
  vt <- x$variance
  vt$variance_0 <- round(vt$variance_0, digits); vt$variance_1 <- round(vt$variance_1, digits)
  vt$prop_reduced <- round(vt$prop_reduced, 3)
  print(vt, row.names = FALSE)
  cat("prop_reduced = 1 - variance_1 / variance_0: the share of each component explained by the added predictors",
      if (any(vt$prop_reduced < 0, na.rm = TRUE)) ";\n  negative values mark components that increased (see ?mpa_compare)" else "",
      ".\n", sep = "")
  if (!is.na(x$delta[["AIC"]]))
    cat("\nDelta AIC = ", round(x$delta[["AIC"]], 1), ", Delta BIC = ", round(x$delta[["BIC"]], 1),
        " (negative favours the comparison model)\n", sep = "")
  if (!is.null(x$lrt))
    cat("Likelihood-ratio test: chi-square(", x$lrt$df, ") = ", round(x$lrt$chisq, 2),
        ", p ", format.pval(x$lrt$p, digits = 3, eps = 1e-4), "\n", sep = "")
  if (!is.null(x$loo)) {
    cat("\nLeave-one-out cross validation (elpd difference to the better model, with SE):\n")
    print(round(x$loo, 1))
  } else if (identical(x$engine, "brms")) {
    cat("(Bayesian fits: no likelihood-ratio test and no AIC/BIC; mpa_compare(..., loo = TRUE) compares by LOO.)\n")
  }
  invisible(x)
}

#' Test random-effect terms by single-term deletion
#'
#' Refits the model once for every random term (intercept or slope) of every
#' unit with that term removed and compares it with the full model by a
#' likelihood-ratio test, AIC and BIC. Because variance components are tested
#' on the boundary of the parameter space, the naive LRT p-value is
#' conservative; `p_boundary` halves it (exact for a single variance component
#' without correlations; an approximation otherwise; Self & Liang, 1987; Stram
#' & Lee, 1994).
#'
#' This is a *diagnostic* tool. Our default recommendation is to keep the
#' design-implied random-effects structure and simplify only in case of
#' convergence or singularity problems (see the package vignette). It is
#' available for the `lme4` engine only: a Bayesian fit reports the posterior
#' interval of every random-effect SD instead.
#'
#' @param formula An `mpa_formula`.
#' @param data An `mpa_data` object or data frame.
#' @param units Units whose terms are tested (default: all in the formula).
#' @param cores Number of parallel workers for the reduced fits (default 1).
#'   Uses a socket cluster (works on Windows); each worker refits one reduced
#'   model at a time, so memory use scales with `cores`.
#' @param ... Passed to [mpa_fit()], e.g. `REML = TRUE` (the default is ML,
#'   as required for comparing models by LRT when random structures differ
#'   only in variance components; REML is also admissible here; ML keeps
#'   results comparable with [mpa_compare()]).
#' @return An object of class `mpa_random_test`: a data frame with one row per
#'   tested term (`unit`, `term`, `n_par`, `chisq`, `df`, `p`, `p_boundary`,
#'   `delta_AIC`, `delta_BIC`, `singular_reduced`), with the full fit stored
#'   in `attr(x, "fit")`.
#' @export
mpa_random_test <- function(formula, data, units = NULL, cores = 1, ...) {
  stopifnot(inherits(formula, "mpa_formula"))
  dots <- list(...)
  engine <- .mpa_resolve_engine(if ("engine" %in% names(dots)) dots$engine else getOption("profilemlm.engine", "lme4"),
                                dots[setdiff(names(dots), "engine")])
  if (engine$name != "lme4")
    stop("mpa_random_test() needs a likelihood and is available for the lme4 engine only. ",
         "For a Bayesian fit, judge the random-effect SDs by their posterior intervals (`$variance`).")
  dots <- list(engine = engine)
  full <- do.call(mpa_fit, c(list(formula = formula, data = data), dots))
  units <- if (is.null(units)) names(formula$random) else intersect(units, names(formula$random))
  jobs <- list()
  for (u in units) for (term in formula$random[[u]]) jobs[[paste(u, term, sep = ":")]] <- list(unit = u, term = term)

  dat <- full$data
  fit_reduced <- function(job) {
    red_f <- mpa_update(formula, drop = stats::setNames(list(job$term), job$unit))
    tryCatch(do.call(mpa_fit, c(list(formula = red_f, data = dat), dots)),
             error = function(e) NULL)
  }
  if (cores > 1 && length(jobs) > 1) {
    cl <- parallel::makeCluster(min(cores, length(jobs)))
    on.exit(parallel::stopCluster(cl), add = TRUE)
    ## Ship the package functions to the workers inside a private environment,
    ## so that this also works with devtools::load_all() (no installed package)
    ## and nothing is written to the workers' global environment.
    ns <- environment(mpa_random_test)
    objs <- Filter(is.function, mget(ls(ns, all.names = TRUE), envir = ns))
    objs <- lapply(objs, function(o) { environment(o) <- emptyenv(); o })   # detach from the namespace for serialisation
    objs$fit_reduced <- fit_reduced
    environment(objs$fit_reduced) <- emptyenv()
    payload <- list(objs = objs, formula = formula, dat = dat, dots = dots, pkgs = .engine_packages(engine))
    reduced <- parallel::clusterApplyLB(cl, jobs, function(job, payload) {
      for (p in payload$pkgs) loadNamespace(p)
      e <- new.env(parent = globalenv())
      for (n in names(payload$objs)) { o <- payload$objs[[n]]; environment(o) <- e; e[[n]] <- o }
      e$formula <- payload$formula; e$dat <- payload$dat; e$dots <- payload$dots
      e$fit_reduced(job)
    }, payload = payload)
    names(reduced) <- names(jobs)
  } else {
    reduced <- lapply(jobs, fit_reduced)
  }

  rows <- list(); fits <- list()
  for (key in names(jobs)) {
      u <- jobs[[key]]$unit; term <- jobs[[key]]$term
      red <- reduced[[key]]
      if (is.null(red)) {
        rows[[key]] <- data.frame(unit = u, term = term, n_par = NA, chisq = NA, df = NA, p = NA,
                                  p_boundary = NA, delta_AIC = NA, delta_BIC = NA, singular_reduced = NA)
        next
      }
      ddf <- full$fit_stats$df - red$fit_stats$df
      chisq <- max(0, 2 * (full$logLik - red$logLik))
      p <- stats::pchisq(chisq, ddf, lower.tail = FALSE)
      rows[[key]] <- data.frame(unit = u, term = term, n_par = ddf, chisq = chisq, df = ddf, p = p,
                                p_boundary = p / 2,
                                delta_AIC = red$AIC - full$AIC, delta_BIC = red$BIC - full$BIC,
                                singular_reduced = red$singular, stringsAsFactors = FALSE)
      fits[[key]] <- red
  }
  out <- do.call(rbind, rows); rownames(out) <- NULL
  attr(out, "fit") <- full
  attr(out, "reduced_fits") <- fits
  class(out) <- c("mpa_random_test", "data.frame")
  out
}

#' @export
print.mpa_random_test <- function(x, digits = 2, ...) {
  cat("<mpa_random_test>  single-term deletions from the full random-effects structure\n")
  cat("(positive delta AIC/BIC: removing the term worsens the criterion, i.e. keep it)\n\n")
  df <- as.data.frame(x)
  df$chisq <- round(df$chisq, digits)
  df$p <- format.pval(df$p, digits = 3, eps = 1e-4); df$p_boundary <- format.pval(df$p_boundary, digits = 3, eps = 1e-4)
  df$delta_AIC <- round(df$delta_AIC, 1); df$delta_BIC <- round(df$delta_BIC, 1)
  print(df, row.names = FALSE)
  invisible(x)
}

#' Prune random terms
#'
#' Two modes. (1) After [mpa_random_test()]: remove terms that fail a
#' likelihood-ratio, AIC or BIC criterion. (2) After a singular fit: pass the
#' `mpa_fit` (as first argument) and terms whose estimated standard deviation is
#' essentially zero (`sd / residual sd < tol`) are removed -- no refits needed.
#' Mode (2) is the natural companion of the workflow in which the random
#' structure is fixed and carried over with `random_structure`: whenever a step
#' adds a new random slope and the fit is singular, prune the degenerate terms
#' and pass the pruned formula on. A Bayesian (`brms`) fit has no boundary
#' solutions; in mode (2) its structure is returned unchanged, with a note on
#' the terms whose posterior median SD is near zero.
#'
#' @param formula The `mpa_formula` that was tested, or an `mpa_fit` (singular
#'   mode).
#' @param test Result of [mpa_random_test()] (test mode); ignored in singular
#'   mode.
#' @param criterion `"lrt"` (keep terms with `p_boundary < alpha`), `"aic"` or
#'   `"bic"` (keep terms whose removal increases the criterion). Singular mode
#'   is chosen automatically when `formula` is an `mpa_fit`.
#' @param alpha Significance level for `criterion = "lrt"` (e.g., `0.05`,
#'   `0.001`).
#' @param keep_intercepts Never remove random intercepts. Default `TRUE` in test
#'   mode and `FALSE` in singular mode (a zero intercept variance is exactly the
#'   degenerate case worth removing).
#' @param tol Relative tolerance in singular mode.
#' @return The pruned `mpa_formula`.
#' @export
mpa_prune <- function(formula, test = NULL, criterion = c("lrt", "aic", "bic"),
                      alpha = 0.05, keep_intercepts = NULL, tol = 1e-3) {
  if (inherits(formula, "mpa_fit")) {
    fit <- formula; formula <- fit$formula
    if (!inherits(formula, "mpa_formula")) stop("The fit must be based on an mpa_formula.")
    if (is.null(keep_intercepts)) keep_intercepts <- FALSE
    vt <- fit$variance
    sigma <- vt$sd[vt$unit == "Residual"]
    vt <- vt[vt$unit != "Residual", ]
    unit_key <- names(formula$unit_var)[match(vt$unit, unlist(formula$unit_var))]
    term <- ifelse(vt$effect == "(Intercept)", "1", vt$effect)
    k <- .mpa_criterion_power(vt$effect, fit$criterion_vars)           # undo the reporting back-transformation
    ratio <- vt$sd * fit$criterion_scale^k / sigma
    remove <- ratio < tol
    if (keep_intercepts) remove <- remove & term != "1"
    drop <- split(term[remove], unit_key[remove])
    if (identical(fit$engine, "brms")) {
      near <- split(term[ratio < 0.05], unit_key[ratio < 0.05])
      message("Bayesian fit: no boundary solutions, the random structure is returned unchanged",
              if (length(near)) paste0(" (posterior median SD below 5% of the residual SD: ",
                                       paste(paste0(names(near), ": ", vapply(near, paste, character(1), collapse = ", ")), collapse = "; "),
                                       "; see the posterior intervals in $variance)") else "", ".")
      return(formula)
    }
    if (length(drop)) message("Removing random terms with (near-)zero variance: ",
                              paste(paste0(names(drop), ": ", vapply(drop, paste, character(1), collapse = ", ")), collapse = "; "))
    else message("No degenerate random terms found.")
    return(mpa_update(formula, drop = drop))
  }
  criterion <- match.arg(criterion)
  if (is.null(test)) stop("Provide an mpa_random_test result as `test`, or an mpa_fit as first argument.")
  if (is.null(keep_intercepts)) keep_intercepts <- TRUE
  t <- as.data.frame(test)
  remove <- switch(criterion,
                   lrt = !(t$p_boundary < alpha),
                   aic = !(t$delta_AIC > 0),
                   bic = !(t$delta_BIC > 0))
  remove[is.na(remove)] <- FALSE
  if (keep_intercepts) remove <- remove & t$term != "1"
  drop <- split(t$term[remove], t$unit[remove])
  if (length(drop)) message("Removing random terms: ",
                            paste(paste0(names(drop), ": ", vapply(drop, paste, character(1), collapse = ", ")), collapse = "; "))
  else message("No random terms removed.")
  mpa_update(formula, drop = drop)
}
