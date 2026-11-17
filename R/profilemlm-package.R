#' profilemlm: Multilevel Profile Analysis of Person Judgment Data
#'
#' Cross-classified linear mixed-effects models for person judgment data in the
#' tradition of the Social Accuracy Model. The workflow is
#' [mpa_data()] (data preparation) -> [mpa_formula()] (model specification) ->
#' [mpa_fit()] (estimation) -> [mpa_compare()] (model comparison) ->
#' [mpa_plots()] (figures), with [mpa_random_test()] / [mpa_prune()] for the
#' random-effects structure and [mpa_simulate()] / [mpa_recovery()] for
#' simulation. See `vignette("mpa-tutorial", package = "profilemlm")`.
#'
#' @keywords internal
#' @importFrom stats aggregate
#' @importFrom utils head
"_PACKAGE"
