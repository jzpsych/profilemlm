## The brms engine: Bayesian estimation via brms::brm() (Stan, through
## cmdstanr or rstan). This file holds all code that touches brmsfit objects.

#' @rdname mpa_engine
#' @param chains,iter,warmup Number of Markov chains, iterations per chain and
#'   warmup iterations (defaults 4, 2000, 1000: 4,000 post-warmup draws).
#' @param cores Number of cores for parallel chains; default: as many as
#'   chains, but at most the number of available cores minus one.
#' @param threads Optional number of threads per chain for within-chain
#'   parallelisation (`brms::threading()`; needs the cmdstanr backend or a
#'   recent rstan).
#' @param prior Optional `brmsprior` object ([brms::set_prior()]); `NULL`
#'   uses brms' default priors: flat on the fixed effects, Student-t(3, 0, s)
#'   on the standard deviations (s derived from the spread of the judgments)
#'   and LKJ(1) on correlations. With these priors the posterior medians of
#'   the fixed effects are close to the ML estimates when the data are
#'   informative; the random-effect SDs are kept away from zero.
#' @param adapt_delta,max_treedepth Sampler settings passed as `control` to
#'   [brms::brm()].
#' @param seed Random seed for the sampler.
#' @param backend `"cmdstanr"` or `"rstan"`; `NULL` (default) uses cmdstanr
#'   if it is installed with a working CmdStan, otherwise rstan.
#' @param file Optional file name (without extension) under which brms saves
#'   the fit and from which it is reloaded on the next call with the same
#'   file, so that a long fit is not repeated (see [brms::brm()]).
#' @param save_group_draws Save the draws of the group-specific effects
#'   (random effects of every perceiver, target, dyad and item)? `FALSE`
#'   (default) keeps the fit small; `TRUE` is required for the leave-one-out
#'   comparison in [mpa_compare()] (`loo = TRUE`).
#' @param verbose Show Stan's sampling progress?
#' @export
engine_brms <- function(chains = 4, iter = 2000, warmup = floor(iter / 2), cores = NULL, threads = NULL,
                        prior = NULL, adapt_delta = 0.9, max_treedepth = 10, seed = NULL,
                        backend = NULL, file = NULL, save_group_draws = FALSE, verbose = FALSE, ...) {
  stopifnot(is.numeric(chains), is.numeric(iter), is.numeric(warmup), warmup < iter)
  if (!is.null(backend)) backend <- match.arg(backend, c("cmdstanr", "rstan"))
  structure(list(name = "brms", chains = as.integer(chains), iter = as.integer(iter), warmup = as.integer(warmup),
                 cores = cores, threads = threads, prior = prior, adapt_delta = adapt_delta,
                 max_treedepth = as.integer(max_treedepth), seed = seed, backend = backend, file = file,
                 save_group_draws = save_group_draws, verbose = verbose, extra = list(...)),
            class = c("mpa_engine_brms", "mpa_engine"))
}

## which Stan backend to use: the one asked for, else cmdstanr if usable, else rstan
.brms_backend <- function(engine) {
  has_cmdstanr <- requireNamespace("cmdstanr", quietly = TRUE) &&
    !is.null(tryCatch(cmdstanr::cmdstan_version(error_on_NA = FALSE), error = function(e) NULL))
  has_rstan <- requireNamespace("rstan", quietly = TRUE)
  if (is.null(engine$backend)) {
    if (has_cmdstanr) return("cmdstanr")
    if (has_rstan) return("rstan")
    stop("The brms engine needs a Stan backend: install cmdstanr (install.packages(\"cmdstanr\", ",
         "repos = c(\"https://stan-dev.r-universe.dev\", getOption(\"repos\"))); cmdstanr::install_cmdstan()) ",
         "or rstan.")
  }
  if (engine$backend == "cmdstanr" && !has_cmdstanr)
    stop("backend = \"cmdstanr\" needs the cmdstanr package and an installed CmdStan (cmdstanr::install_cmdstan()).")
  if (engine$backend == "rstan" && !has_rstan) stop("backend = \"rstan\" needs the rstan package.")
  engine$backend
}

.engine_check.mpa_engine_brms <- function(engine, ...) {
  if (!requireNamespace("brms", quietly = TRUE))
    stop("The brms engine needs the brms package (install.packages(\"brms\")).")
  if (!requireNamespace("posterior", quietly = TRUE)) stop("The brms engine needs the posterior package.")
  .brms_backend(engine)
  invisible(TRUE)
}

.engine_packages.mpa_engine_brms <- function(engine, ...) c("brms", "posterior", .brms_backend(engine))

.engine_label.mpa_engine_brms <- function(engine, native, ...)
  paste0("Bayesian fit (brms, ", if (!is.null(native$backend)) native$backend else .brms_backend(engine), ")")

.engine_start.mpa_engine_brms <- function(engine, formula, data, old_fit, ...) {
  message("`start_from` is not used by the brms engine (the warmup finds its own starting values).")
  NULL
}

.engine_fit.mpa_engine_brms <- function(engine, formula, data, start = NULL, ...) {
  backend <- .brms_backend(engine)
  cores <- engine$cores
  if (is.null(cores)) {
    nc <- tryCatch(parallel::detectCores(), error = function(e) NA_integer_)
    cores <- if (is.na(nc)) 1L else max(1L, min(engine$chains, nc - 1L))
  }
  environment(formula) <- globalenv()
  args <- list(formula = formula, data = data, family = stats::gaussian(),
               chains = engine$chains, iter = engine$iter, warmup = engine$warmup, cores = cores,
               control = list(adapt_delta = engine$adapt_delta, max_treedepth = engine$max_treedepth),
               save_pars = brms::save_pars(group = engine$save_group_draws),
               backend = backend, refresh = if (engine$verbose) NULL else 0, silent = if (engine$verbose) 1 else 2)
  if (!is.null(engine$prior)) args$prior <- engine$prior
  if (!is.null(engine$seed)) args$seed <- engine$seed
  if (!is.null(engine$file)) args$file <- engine$file
  if (!is.null(engine$threads)) args$threads <- brms::threading(engine$threads)
  args$refresh <- args$refresh   # NULL entries are dropped by do.call
  args <- c(args[!vapply(args, is.null, logical(1))], engine$extra)
  fit <- do.call(brms::brm, args)
  if (stats::nobs(fit) != nrow(data))
    stop("Internal error: brms used ", stats::nobs(fit), " of ", nrow(data), " rows.")
  fit
}

## brms' renaming of coefficient names (brms:::rename with its defaults)
.brms_rename <- function(x) {
  pattern <- c(" ", "(", ")", "[", "]", ",", "\"", "'", "?", "+", "-", "*", "/", "^", "=")
  replacement <- c(rep("", 9), "P", "M", "MU", "D", "E", "EQ")
  for (i in seq_along(pattern)) x <- gsub(pattern[i], replacement[i], x, fixed = TRUE)
  x
}

## the canonical (model.matrix) names of the fixed effects of a brmsfit, in
## the order of its b_ parameters, checked against brms' own renaming
.brms_fixed_names <- function(native) {
  f <- native$formula$formula
  tt <- stats::delete.response(stats::terms(.mpa_nobars(f)))
  canon <- colnames(stats::model.matrix(tt, native$data))
  bn <- sub("^b_", "", grep("^b_", brms::variables(native), value = TRUE))
  if (!identical(.brms_rename(canon), bn))
    stop("Internal error: the fixed-effect names of the brms fit (", paste(bn, collapse = ", "),
         ") do not match the model matrix (", paste(canon, collapse = ", "), ").")
  stats::setNames(canon, bn)
}

## lme4 orders the random-effect groups by decreasing number of levels (ties:
## reversed formula order); the brms tables follow the same order so that both
## engines print the units identically
.brms_unit_order <- function(native) {
  re <- native$ranef
  groups <- unique(re$group)
  bars <- .mpa_findbars(native$formula$formula)
  term_groups <- vapply(bars, function(b) paste(deparse(b[[3]]), collapse = ""), character(1))
  nl <- vapply(term_groups, function(g) length(attr(re, "levels")[[g]]), integer(1))
  if (any(diff(nl) > 0)) term_groups <- term_groups[rev(order(nl))]
  c(unique(term_groups), setdiff(groups, term_groups))
}

.brms_summary <- function(native, level, pattern) {
  draws <- brms::as_draws_array(native, variable = pattern, regex = TRUE)
  probs <- c((1 - level) / 2, 1 - (1 - level) / 2)
  s <- posterior::summarise_draws(draws, median = stats::median, sd = stats::sd,
                                  ~posterior::quantile2(.x, probs = probs),
                                  rhat = posterior::rhat, ess_bulk = posterior::ess_bulk, ess_tail = posterior::ess_tail)
  s <- as.data.frame(s)
  names(s)[4:5] <- c("lower", "upper")
  s
}

.engine_tables.mpa_engine_brms <- function(engine, native, level, ...) {
  canon <- .brms_fixed_names(native)
  s <- .brms_summary(native, level, "^(b_|sd_|cor_|sigma$)")
  ## fixed effects
  fx <- s[grepl("^b_", s$variable), ]
  fixed <- data.frame(term = unname(canon[sub("^b_", "", fx$variable)]), estimate = fx$median, se = fx$sd,
                      ci_lower = fx$lower, ci_upper = fx$upper, rhat = fx$rhat, ess_bulk = fx$ess_bulk,
                      ess_tail = fx$ess_tail, stringsAsFactors = FALSE)
  ## variance components
  re <- native$ranef
  back <- function(coef) { m <- canon[match(coef, names(canon))]; ifelse(is.na(m), coef, unname(m)) }
  sds <- s[grepl("^sd_", s$variable), ]
  key <- sub("^sd_", "", sds$variable)
  unit <- sub("__.*$", "", key); coef <- sub("^.*?__", "", key)
  variance <- data.frame(unit = unit, effect = back(coef), variance = sds$median^2, sd = sds$median,
                         sd_lower = sds$lower, sd_upper = sds$upper, rhat = sds$rhat, ess_bulk = sds$ess_bulk,
                         stringsAsFactors = FALSE)
  sg <- s[s$variable == "sigma", ]
  variance <- rbind(variance, data.frame(unit = "Residual", effect = "", variance = sg$median^2, sd = sg$median,
                                         sd_lower = sg$lower, sd_upper = sg$upper, rhat = sg$rhat, ess_bulk = sg$ess_bulk))
  units <- .brms_unit_order(native)
  variance <- variance[order(match(variance$unit, c(units, "Residual")), variance$effect != "(Intercept)",
                             match(variance$effect, fixed$term)), ]
  rownames(variance) <- NULL
  ## correlations
  cors <- s[grepl("^cor_", s$variable), ]
  correlations <- NULL
  if (nrow(cors)) {
    key <- sub("^cor_", "", cors$variable)
    parts <- strsplit(key, "__", fixed = TRUE)
    correlations <- data.frame(unit = vapply(parts, `[`, character(1), 1), effect1 = back(vapply(parts, `[`, character(1), 2)),
                               effect2 = back(vapply(parts, `[`, character(1), 3)), correlation = cors$median,
                               lower = cors$lower, upper = cors$upper, rhat = cors$rhat, stringsAsFactors = FALSE)
    correlations <- correlations[order(match(correlations$unit, units)), ]
    rownames(correlations) <- NULL
  }
  n_units <- vapply(units, function(g) length(attr(re, "levels")[[g]]), integer(1))
  pct <- round(100 * level)
  list(fixed = fixed, variance = variance, correlations = correlations, n_units = n_units,
       labels = list(estimate = "posterior medians", test = "posterior SD in `se`", interval = paste0(pct, "% credible interval"),
                     variance = paste0("posterior medians; ", pct, "% credible interval of the SD")))
}

.engine_diagnostics.mpa_engine_brms <- function(engine, native, ...) {
  s <- .brms_summary(native, 0.95, "^(b_|sd_|cor_|sigma$)")
  np <- tryCatch(brms::nuts_params(native), error = function(e) NULL)
  div <- if (!is.null(np)) sum(np$Value[np$Parameter == "divergent__"]) else NA_integer_
  td <- if (!is.null(np)) sum(np$Value[np$Parameter == "treedepth__"] >= engine$max_treedepth) else NA_integer_
  max_rhat <- max(s$rhat, na.rm = TRUE); min_ess <- min(s$ess_bulk, s$ess_tail, na.rm = TRUE)
  msgs <- character(0)
  if (is.finite(max_rhat) && max_rhat > 1.05) msgs <- c(msgs, sprintf("largest R-hat %.3f (> 1.05): chains have not mixed", max_rhat))
  if (is.finite(min_ess) && min_ess < 100 * dim(posterior::as_draws_array(native, variable = "sigma"))[2])
    msgs <- c(msgs, sprintf("smallest effective sample size %.0f: increase `iter`", min_ess))
  if (isTRUE(div > 0)) msgs <- c(msgs, sprintf("%d divergent transitions: increase `adapt_delta`", div))
  if (isTRUE(td > 0)) msgs <- c(msgs, sprintf("%d transitions hit the maximum tree depth: increase `max_treedepth`", td))
  list(converged = length(msgs) == 0, singular = FALSE, messages = msgs,
       details = list(max_rhat = max_rhat, min_ess = min_ess, divergent = div, max_treedepth_hits = td))
}

.engine_fit_stats.mpa_engine_brms <- function(engine, native, ...) {
  dr <- dim(posterior::as_draws_array(native, variable = "sigma"))
  elapsed <- tryCatch({
    if (identical(native$backend, "rstan")) sum(rstan::get_elapsed_time(native$fit)) / dr[2]   # per chain (parallel)
    else native$fit$time()$total
  }, error = function(e) NA_real_)
  list(chains = dr[2], draws = dr[1] * dr[2], iter = engine$iter, warmup = engine$warmup,
       backend = native$backend, elapsed_seconds = elapsed)
}

## draws of L %*% beta: median, SD and equal-tailed interval
.engine_lincomb.mpa_engine_brms <- function(engine, native, L, level, ...) {
  canon <- .brms_fixed_names(native)
  D <- posterior::as_draws_matrix(brms::as_draws_array(native, variable = "^b_", regex = TRUE))
  colnames(D) <- unname(canon[sub("^b_", "", colnames(D))])
  D <- D[, colnames(L), drop = FALSE]
  S <- D %*% t(L)
  probs <- c((1 - level) / 2, 1 - (1 - level) / 2)
  q <- apply(S, 2, stats::quantile, probs = probs)
  data.frame(estimate = apply(S, 2, stats::median), se = apply(S, 2, stats::sd), lower = q[1, ], upper = q[2, ])
}

.engine_print_diagnostics.mpa_engine_brms <- function(engine, x, ...) {
  d <- x$diagnostics; s <- x$fit_stats
  cat("Sampling: ", s$chains, " chains x ", format(s$draws / s$chains, big.mark = ","), " post-warmup draws (", s$backend,
      if (is.finite(s$elapsed_seconds)) paste0(", ", round(s$elapsed_seconds / 60, 1), " min") else "", "); ",
      "max R-hat ", formatC(d$details$max_rhat, digits = 3, format = "f"), ", min ESS ", round(d$details$min_ess),
      if (!is.na(d$details$divergent)) paste0(", ", d$details$divergent, " divergent transitions") else "", ".\n", sep = "")
  if (!d$converged) cat("WARNING: ", paste(d$messages, collapse = "; "), "\n", sep = "")
  else cat("No sampling problems detected.\n")
  invisible(NULL)
}

## leave-one-out comparison of two brms fits (loo::loo_compare)
.mpa_loo_compare <- function(fit0, fit1) {
  if (!identical(fit0$engine, "brms") || !identical(fit1$engine, "brms"))
    stop("`loo = TRUE` is for brms fits; lme4 fits are compared by the likelihood-ratio test and AIC/BIC.")
  if (!requireNamespace("loo", quietly = TRUE)) stop("`loo = TRUE` needs the loo package.")
  get_loo <- function(f, nm) {
    if (!isTRUE(f$settings$save_group_draws))
      stop("`loo = TRUE` needs fits estimated with `save_group_draws = TRUE` (", nm, " was not).")
    brms::loo(f$fit, pointwise = f$n_obs > 20000)
  }
  l0 <- get_loo(fit0, "the reference model"); l1 <- get_loo(fit1, "the comparison model")
  tab <- loo::loo_compare(list(model_0 = l0, model_1 = l1))
  out <- as.data.frame(unclass(tab))
  out <- cbind(model = rownames(out), out); rownames(out) <- NULL
  out
}
