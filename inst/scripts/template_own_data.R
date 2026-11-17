###########################################################
### profilemlm: template for running the five steps on your own data
###
### Copy this file, replace the parts marked <...> and run it top to bottom.
### Every step prints its model; the comments say what to look at. The
### tutorial vignette (vignette("mpa-tutorial", package = "profilemlm"))
### walks through the same steps on simulated data with full output.
###########################################################

library(profilemlm)
# options(profilemlm.engine = "brms")   # uncomment to estimate every model below with brms (Bayesian) instead of lme4

## ---- 1. Read the files as your study produced them --------------------------
## Judgments: one row per perceiver with the target blocks ("loops") in columns
## (format = "loop"), one row per perceiver-target pair (format = "wide"), or
## one row per rating (format = "long", with rating_col).
perceiver_survey <- read.csv("<perceiver_survey.csv>")
target_survey    <- read.csv("<target_survey.csv>")        # one row per target; the criterion (e.g., self-ratings) in item columns
items            <- read.csv("<items.csv>")                # one row per item; item properties (e.g., desirability)

## ---- 2. Build the analysis data set -----------------------------------------
dat <- mpa_data(judgments = perceiver_survey, format = "loop",
              loop_id   = "<id>",                           # the column inside each loop that holds the target's ID
              targets   = target_survey, items = items,
              item_cols = "<^item_>",                       # regular expression that picks the item columns
              criterion = list(self = "targets"),           # the criterion source: the target file ...
              group     = "<group>",                        # ... or NULL if there are no groups
              rating_scale = c(<1>, <101>),                 # minimum and maximum of the rating scale
              rescale = list(<predictor> = c(<min>, <max>)))   # rescales predictors to -1...+1; omit if not wanted
dat                     # check: design, criterion source, level of every predictor
summary(dat)            # descriptive statistics; item-level correlations of the normative profile with item properties

## ---- 3. Step 0: where is the variance? --------------------------------------
f0 <- mpa_formula(dat, criterion = NULL)
m0 <- mpa_fit(f0, dat)
m0                    # shares of perceiver, target, dyad, item and residual variance

## ---- 4. Step 1: the base model ----------------------------------------------
f1 <- mpa_formula(dat, criterion = "self")
m1 <- mpa_fit(f1, dat)
m1                    # normative and distinctive agreement; their SDs across units
mpa_compare(m0, m1)   # how much item and residual variance does the criterion profile explain?

## Deciding on the random-effects structure: if m1 is singular, remove the
## zero-variance terms and/or drop the correlations; keep the result (f_base)
## fixed for all later steps.
f_base <- mpa_update(f1, correlated = FALSE)
m_base <- mpa_fit(f_base, dat)
if (m_base$singular) { f_base <- mpa_prune(m_base); m_base <- mpa_fit(f_base, dat) }

## ---- 5. Step 2: main effects of item properties ---------------------------
f2 <- mpa_formula(dat, criterion = "self", random_slopes = "<des>", random_structure = f_base)
m2 <- mpa_fit(f2, dat)
if (m2$singular) { f_base <- mpa_prune(m2); m2 <- mpa_fit(f_base, dat) }   # a new random slope may be degenerate
m2
mpa_compare(m_base, m2)

## ---- 6. Step 3: item properties as moderators of agreement -----------------
f3 <- mpa_formula(dat, criterion = "self", random_slopes = "<des>", poly = "<des>",
                  interactions = "criterion", random_structure = f_base)
m3 <- mpa_fit(f3, dat, start_from = m2)
m3                    # <property> x self.dist: change in distinctive agreement per unit of the property
mpa_compare(m2, m3)

## ---- 7. Step 4: perceiver, target and dyad-level predictors -----------------
f4 <- mpa_formula(dat, criterion = "self", random_slopes = "<des>", poly = "<des>",
                  predictors = c("<know>", "<training>", "<talkativeness>"),
                  interactions = "criterion", random_structure = f_base)
m4 <- mpa_fit(f4, dat, start_from = m3)
m4
mpa_compare(m3, m4)      # what Step 4 adds
mpa_compare(m_base, m4)  # how much of the slope variance of Step 1 do all predictors explain?
mpa_effect_size(m4)      # standardized effect sizes d (see ?mpa_effect_size)

## Optional: decompose a dyadic variable into group, perceiver, target and
## relationship components (needs perceivers who judge several targets and
## targets judged by several perceivers, i.e., a half-block or round-robin design)
## dd  <- mpa_data(..., decompose = "<know>")
## f4d <- mpa_formula(dd, ..., predictors = c("<know>.g", "<know>.p", "<know>.t", "<know>.r", ...), ...)

## Optional: a shared target x item component (what the perceivers of a target
## agree on beyond the criterion profile; see the paper, Extensions)
## f4ti <- mpa_update(f4, add = list(target_item = "1"))
## m4ti <- mpa_fit(f4ti, dat, start_from = m4)
## ... and its perceiver-level counterpart, a perceiver x item component (a perceiver's item-specific
## view of people in general, plus item-specific distinctive agreement); needs several
## targets per perceiver
## f4pi <- mpa_update(f4ti, add = list(perceiver_item = c("1", "self.dist")))
## m4pi <- mpa_fit(f4pi, dat, start_from = m4ti)

## ---- 8. Figures -------------------------------------------------------------
if (requireNamespace("ggplot2", quietly = TRUE)) {
  p <- mpa_plots(m4, labels = c(<des> = "Desirability", <know> = "Acquaintance"))
  print(p[["self.dist_by_<know>"]])
  mpa_save_plots(p, dir = "<figures>")
}
