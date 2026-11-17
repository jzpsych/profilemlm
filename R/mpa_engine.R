#' Estimation engines
#'
#' [mpa_fit()] estimates a model with one of two engines: `lme4` (maximum
#' likelihood, the default) or `brms` (Bayesian, Stan via cmdstanr or rstan).
#' An engine object holds the settings of the estimation; it is created by
#' `engine_lme4()` or `engine_brms()` and passed to `mpa_fit(engine = )`, or
#' the settings are given directly to `mpa_fit()` as further arguments, which
#' is the same thing (`mpa_fit(f, dat, optimizer = "nloptwrap")` is
#' `mpa_fit(f, dat, engine = engine_lme4(optimizer = "nloptwrap"))`). The option
#' `profilemlm.engine` sets the default engine of all `mpa_fit()` calls of a
#' session: `options(profilemlm.engine = "brms")`.
#'
#' Both engines fit the same model: the same fixed effects, the same random
#' terms per unit, the same correlation assumptions (`|` vs. `||`) and the
#' same data rows. What differs is the estimation and therefore what the
#' tables of the fit contain; see [mpa_fit()], section "Engines".
#'
#' @param REML Use REML (`TRUE`) or ML (`FALSE`, default; required for
#'   likelihood-ratio comparisons of models with different fixed effects).
#' @param df How to obtain p-values for fixed effects: `"z"` (Wald z-tests;
#'   default, adequate for the typically large number of observations) or
#'   `"satterthwaite"` (requires the `lmerTest` package).
#' @param optimizer Optimizer passed to [lme4::lmerControl()]. The default,
#'   `"bobyqa"`, converged without gradient warnings in all tutorial models and
#'   reproduced published results where lme4's default (`"nloptwrap"`) landed
#'   on a singular solution; `"nloptwrap"` is somewhat faster on large data
#'   sets. If a fit reports a gradient slightly above lme4's tolerance, refit
#'   with the other optimizer; the estimates should agree.
#' @param check_derivs Run lme4's post-fit gradient/Hessian check
#'   (`calc.derivs`)? `FALSE` skips it, which can save a substantial share of
#'   the fitting time for large models, at the price of losing the convergence
#'   warnings that are based on it (singularity is still detected).
#' @param control An optional `lmerControl` object (overrides `optimizer` and
#'   `check_derivs`).
#' @param ... `engine_lme4()`: passed on to [lme4::lmer()]. `engine_brms()`:
#'   passed on to [brms::brm()].
#' @return An object of class `mpa_engine` (and `mpa_engine_lme4` or
#'   `mpa_engine_brms`): a list of the settings.
#' @examples
#' engine_lme4()
#' engine_lme4(optimizer = "nloptwrap", check_derivs = FALSE)
#' @name mpa_engine
#' @export
engine_lme4 <- function(REML = FALSE, df = c("z", "satterthwaite"), optimizer = "bobyqa",
                        check_derivs = TRUE, control = NULL, ...) {
  df <- match.arg(df)
  stopifnot(is.logical(REML), length(REML) == 1, is.character(optimizer), length(optimizer) == 1,
            is.logical(check_derivs), length(check_derivs) == 1)
  structure(list(name = "lme4", REML = REML, df = df, optimizer = optimizer, check_derivs = check_derivs,
                 control = control, extra = list(...)),
            class = c("mpa_engine_lme4", "mpa_engine"))
}

#' @export
print.mpa_engine <- function(x, ...) {
  cat("<mpa_engine: ", x$name, ">\n", sep = "")
  s <- x[setdiff(names(x), c("name", "extra", "control"))]
  for (nm in names(s)) {
    v <- s[[nm]]
    cat("  ", format(nm, width = 16), if (is.null(v)) "NULL" else if (is.function(v)) "<function>" else if (length(v) > 1 || is.list(v)) paste(deparse(v), collapse = " ") else format(v), "\n", sep = "")
  }
  if (!is.null(x$control)) cat("  ", format("control", width = 16), "<", class(x$control)[1], ">\n", sep = "")
  if (length(x$extra)) cat("  ", format("...", width = 16), paste(names(x$extra), collapse = ", "), "\n", sep = "")
  invisible(x)
}

## Turn `engine` (a name or an engine object) plus further settings into an
## engine object. Settings given as `...` to mpa_fit() override those of the
## object, so that `mpa_fit(f, dat, engine = e, REML = TRUE)` works as expected.
.mpa_resolve_engine <- function(engine, dots = list()) {
  if (inherits(engine, "mpa_fit")) engine <- engine$settings      # "same engine as that fit"
  if (is.character(engine)) {
    if (length(engine) != 1 || !engine %in% c("lme4", "brms"))
      stop("`engine` must be \"lme4\" or \"brms\" (or an engine object, see ?mpa_engine).")
    if (engine == "lme4") return(do.call(engine_lme4, dots))
    return(do.call(engine_brms, dots))
  }
  if (!inherits(engine, "mpa_engine"))
    stop("`engine` must be \"lme4\", \"brms\" or an engine object created by engine_lme4() / engine_brms().")
  if (!length(dots)) return(engine)
  ctor <- if (engine$name == "lme4") engine_lme4 else engine_brms
  args <- unclass(engine); args$name <- NULL
  extra <- args$extra; args$extra <- NULL
  do.call(ctor, c(utils::modifyList(args, dots[names(dots) %in% names(args)]), extra, dots[!names(dots) %in% names(args)]))
}

## ---- the adapter contract --------------------------------------------------
## Every engine implements these internal generics. All other code of the
## package works on the canonical `mpa_fit` object they produce and never
## touches the native fit (merMod, brmsfit) directly.

## check that the engine's packages are installed; error with instructions otherwise
.engine_check <- function(engine, ...) UseMethod(".engine_check")
## packages a parallel worker must load to refit with this engine
.engine_packages <- function(engine, ...) UseMethod(".engine_packages")
## short label for the print method, e.g. "ML fit (lme4)"
.engine_label <- function(engine, native, ...) UseMethod(".engine_label")
## starting values derived from a previous fit (NULL if not supported)
.engine_start <- function(engine, formula, data, old_fit, ...) UseMethod(".engine_start")
## fit the model: lme4-syntax formula, data (used rows only) -> native object
.engine_fit <- function(engine, formula, data, start = NULL, ...) UseMethod(".engine_fit")
## canonical tables: list(fixed, variance, correlations, n_units, labels)
.engine_tables <- function(engine, native, level, ...) UseMethod(".engine_tables")
## list(converged, singular, messages, details)
.engine_diagnostics <- function(engine, native, ...) UseMethod(".engine_diagnostics")
## list(logLik, df, AIC, BIC) for likelihood engines; sampling information otherwise
.engine_fit_stats <- function(engine, native, ...) UseMethod(".engine_fit_stats")
## estimates and intervals of linear combinations L %*% beta of the fixed effects
.engine_lincomb <- function(engine, native, L, level, ...) UseMethod(".engine_lincomb")
## print the convergence / sampling block of print.mpa_fit()
.engine_print_diagnostics <- function(engine, x, ...) UseMethod(".engine_print_diagnostics")

## lme4 >= 2.0 moved the formula helpers to the reformulas package and warns when
## its own copies are called; use reformulas where it is installed
.mpa_nobars <- function(f) {
  if (requireNamespace("reformulas", quietly = TRUE)) reformulas::nobars(f) else lme4::nobars(f)
}
.mpa_findbars <- function(f) {
  if (requireNamespace("reformulas", quietly = TRUE)) reformulas::findbars(f) else lme4::findbars(f)
}
