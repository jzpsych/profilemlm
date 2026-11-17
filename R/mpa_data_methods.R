#' Detect the design type of a perceiver x target data set
#'
#' Classifies the crossing structure of perceivers and targets as pairs,
#' one-with-many, half-block or round-robin and reports whether dyadic
#' (perceiver:target) random effects are separable from perceiver and target
#' random effects. Dyadic effects can be separated from the residual only if
#' several items were rated per dyad; they can be separated from perceiver and
#' target effects only if perceivers rated several targets *and* targets were
#' rated by several perceivers.
#'
#' @param x An `mpa_data` object, or any data frame with columns `perceiver`,
#'   `target` and optionally `group`.
#' @return A list with elements `type`, `n_perceiver`, `n_target`, `n_dyad`,
#'   `n_group`, `targets_per_perceiver`, `perceivers_per_target`,
#'   `overlap` (number of persons who are both perceiver and target),
#'   `self_ratings` (any perceiver == target rows), `crossed_within_groups`,
#'   `pairs` (every perceiver has one target and every target one perceiver,
#'   so that perceiver, target and dyad are the same unit), and
#'   `dyad_re_identifiable`.
#' @export
mpa_design <- function(x) {
  cols <- c("perceiver", "target", if ("group" %in% names(x)) "group")
  pt <- unique(as.data.frame(x)[, cols, drop = FALSE])
  pt$perceiver <- as.character(pt$perceiver)
  pt$target    <- as.character(pt$target)

  n_p <- length(unique(pt$perceiver))
  n_t <- length(unique(pt$target))
  n_d <- nrow(unique(pt[, c("perceiver", "target")]))
  tpp <- table(pt$perceiver)
  ppt <- table(pt$target)
  overlap <- length(intersect(unique(pt$perceiver), unique(pt$target)))
  self_ratings <- any(pt$perceiver == pt$target)

  n_group <- NA_integer_
  crossed_within <- NA
  if ("group" %in% names(pt)) {
    n_group <- length(unique(pt$group))
    crossed_within <- all(vapply(split(pt, pt$group), function(d) {
      nrow(unique(d[, c("perceiver", "target")])) ==
        length(unique(d$perceiver)) * length(unique(d$target))
    }, logical(1)))
  }

  one_target_each    <- all(tpp == 1)
  one_perceiver_each <- all(ppt == 1)
  type <- if (one_target_each && one_perceiver_each) {
    "pairs (each perceiver rated one target, each target was rated once)"
  } else if (one_target_each) {
    "one-with-many (each perceiver rated exactly one target)"
  } else if (one_perceiver_each) {
    "one-with-many (each target was rated by exactly one perceiver)"
  } else if (overlap / min(n_p, n_t) > 0.5) {
    "round-robin (perceivers and targets are largely the same persons)"
  } else if (n_d == n_p * n_t) {
    "half-block (fully crossed)"
  } else if (isTRUE(crossed_within)) {
    "half-block (fully crossed within groups)"
  } else {
    "half-block (partially crossed)"
  }

  list(type = type,
       n_perceiver = n_p, n_target = n_t, n_dyad = n_d, n_group = n_group,
       targets_per_perceiver = as.integer(range(tpp)),
       perceivers_per_target = as.integer(range(ppt)),
       overlap = overlap,
       self_ratings = self_ratings,
       crossed_within_groups = crossed_within,
       pairs = one_target_each && one_perceiver_each,
       dyad_re_identifiable = !(one_target_each || one_perceiver_each))
}


#' @export
print.mpa_data <- function(x, hint = TRUE, ...) {
  m <- attr(x, "mpa")
  if (is.null(m)) { print(as.data.frame(x)); return(invisible(x)) }
  d <- m$design
  fmt <- function(n) format(n, big.mark = ",")
  rng <- function(r) if (length(r) == 2 && !anyNA(r)) { if (r[1] == r[2]) r[1] else paste(r[1], r[2], sep = "-") } else "NA"

  cat("<mpa_data> ", fmt(nrow(x)), " rows: ",
      d$n_perceiver, " perceivers x ", d$n_target, " targets x ", length(m$item_cols), " items",
      " (", fmt(d$n_dyad), " dyads, ", rng(d$items_per_dyad), " rated items per dyad)\n", sep = "")
  cat("Design:  ", d$type, "\n", sep = "")
  if (!is.na(d$n_group)) cat("         ", d$n_group, " groups\n", sep = "")
  cat("         ", rng(d$targets_per_perceiver), " targets per perceiver; ",
      rng(d$perceivers_per_target), " perceivers per target\n", sep = "")
  cat("         dyadic random effects (perceiver:target) ",
      if (d$dyad_re_identifiable) "are" else "are NOT", " separable from perceiver/target effects\n", sep = "")
  if (isTRUE(d$self_ratings)) cat("         NOTE: some perceivers rated themselves (perceiver == target)\n")

  cat("Judgments: y",
      if (!is.null(m$settings$rating_scale)) paste0(" (scale ", m$settings$rating_scale[1], "-", m$settings$rating_scale[2], ")"),
      if (length(m$settings$judgment_role)) paste0(", rater role: ", paste(m$settings$judgment_role, collapse = ", ")),
      "\n", sep = "")

  if (length(m$vars$criterion)) {
    cat("Criterion sources: ")
    src <- vapply(m$vars$criterion, function(s) {
      n <- m$criterion_n[[s]]
      if (is.null(n)) return(paste0(s, " (mean of ", paste(setdiff(m$vars$criterion, s), collapse = ", "), ")"))
      nr <- range(as.integer(n))
      paste0(s, " (", length(n), " targets, ", rng(nr), if (nr[2] == 1) " rater)" else " raters)")
    }, character(1))
    cat(paste(src, collapse = "; "), "\n", sep = "")
    cat("Normative profile: ", switch(m$settings$norm,
                                      loo = "leave-one-out sample mean", sample = "sample mean", external = "external"),
        if (!is.null(m$settings$norm_by)) paste0(", computed within `", m$settings$norm_by, "`"),
        "; centered at ", switch(m$settings$center_norm, grand_mean = "grand mean", midpoint = "scale midpoint", none = "(not centered)"),
        "\n", sep = "")
    cat("Distinctive profile: criterion - normative profile",
        if (isTRUE(m$settings$criterion_scale != 1)) paste0("; analysis variables .norm.c and .dist divided by ", signif(m$settings$criterion_scale, 4), " (results are reported per criterion point)"),
        "\n", sep = "")
  } else {
    cat("Criterion sources: none\n")
  }
  if (length(m$vars$profiles))
    cat("Further predictor profiles: ",
        paste0(m$vars$profiles, " (", m$vars$profile_key[m$vars$profiles], " x item, ",
               vapply(m$vars$profiles, function(p) length(m$profile_n[[p]]), integer(1)), " ",
               m$vars$profile_key[m$vars$profiles], "s)", collapse = "; "),
        "; distinctive component: <name>.dist\n", sep = "")
  vr <- function(v) if (length(v)) paste(v, collapse = ", ") else "-"
  cat("Variables: item: ", vr(m$vars$item), "\n",
      "           dyad: ", vr(m$vars$dyad), "\n",
      "           perceiver: ", vr(m$vars$perceiver), "\n",
      "           target: ", vr(m$vars$target), "\n", sep = "")
  if (length(m$vars$group)) cat("           group: ", vr(m$vars$group), "\n", sep = "")
  if (hint) cat("Use summary() for descriptive statistics; as.data.frame() for the plain data.\n")
  invisible(x)
}


#' @export
summary.mpa_data <- function(object, digits = 2, ...) {
  m <- attr(object, "mpa")
  print(object, hint = FALSE)
  cat("\n")

  if (length(m$criterion_n)) {
    cat("Raters per target and criterion source:\n")
    for (s in names(m$criterion_n)) {
      n <- as.integer(m$criterion_n[[s]])
      tab <- table(factor(n))
      cat("  ", s, ": ", length(n), " targets; raters per target: ",
          paste0(names(tab), ifelse(names(tab) == "1", " rater (", " raters ("), as.integer(tab), " targets)", collapse = ", "), "\n", sep = "")
    }
    no_criterion <- setdiff(levels(object$target), unlist(lapply(m$criterion_n, names)))
    if (length(no_criterion)) cat("  ", length(no_criterion), " target(s) without any criterion data\n", sep = "")
    cat("\n")
  }

  ## item-level correlations of the normative profiles and the item SDs of the
  ## criterion with item variables (e.g., desirability): how much of "normative"
  ## is "desirable", and do item properties go along with the spread of the criterion?
  sd_vars <- m$vars$sd %||% character(0)
  item_vars <- setdiff(m$vars$item, sd_vars)
  item_vars <- item_vars[vapply(item_vars, function(v) is.numeric(object[[v]]), logical(1))]
  rows <- c(m$vars$norm, sd_vars)
  if (length(rows) && length(item_vars)) {
    it <- object[!duplicated(object$item), c("item", rows, item_vars), drop = FALSE]
    ok <- vapply(it[item_vars], function(v) stats::sd(v, na.rm = TRUE) > 0, logical(1))
    if (any(ok)) {
      cr <- sapply(item_vars[ok], function(v) sapply(rows, function(nm) stats::cor(it[[nm]], it[[v]], use = "complete.obs")))
      cr <- matrix(cr, nrow = length(rows), dimnames = list(rows, item_vars[ok]))
      cat("Item-level correlations of the normative profile(s) and the item SD(s) of the criterion with item variables (",
          nrow(it), " items):\n", sep = "")
      print(round(cr, 2))
      cat("\n")
    }
  }

  desc <- function(v) {
    v <- as.numeric(v)
    c(n = sum(!is.na(v)), missing = sum(is.na(v)),
      mean = mean(v, na.rm = TRUE), sd = stats::sd(v, na.rm = TRUE),
      min = suppressWarnings(min(v, na.rm = TRUE)), max = suppressWarnings(max(v, na.rm = TRUE)))
  }
  num_vars <- c(m$vars$judgment, m$vars$criterion, m$vars$norm, m$vars$norm_c, m$vars$dist,
                if (length(m$vars$profiles)) paste0(m$vars$profiles, ".dist"),
                m$vars$item, m$vars$dyad, m$vars$perceiver, m$vars$target, m$vars$group)
  num_vars <- num_vars[num_vars %in% names(object)]
  num_vars <- num_vars[vapply(object[num_vars], is.numeric, logical(1))]
  if (length(num_vars)) {
    tab <- t(vapply(object[num_vars], desc, numeric(6)))
    tab <- round(tab, digits)
    cat("Descriptive statistics (row level):\n")
    print(tab)
  }
  invisible(list(design = m$design, criterion_n = m$criterion_n,
                 descriptives = if (length(num_vars)) tab else NULL))
}


#' @export
as.data.frame.mpa_data <- function(x, ...) {
  attr(x, "mpa") <- NULL
  class(x) <- "data.frame"
  x
}

## Subsetting returns a plain data frame: the design metadata would no longer
## describe the subset. Rebuild with mpa_data() if an mpa_data object is needed.
#' @export
`[.mpa_data` <- function(x, ...) {
  out <- NextMethod()
  if (is.data.frame(out)) {
    attr(out, "mpa") <- NULL
    class(out) <- "data.frame"
  }
  out
}
