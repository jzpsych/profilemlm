# profilemlm

Multilevel profile analysis (MPA) of person judgment data in R.

Several people describe someone they know on a list of items, and their descriptions agree with a criterion, for example that person's self-ratings, to some degree. MPA asks who agrees with whom about what. It regresses each judgment profile on the normative and the distinctive component of the criterion profile, in the tradition of the Social Accuracy Model (Biesanz, 2010), and lets the two agreement slopes vary across perceivers, targets, dyads, and items. Items are a crossed random factor, so conclusions generalize beyond the items of a study, and item properties such as desirability can be tested as predictors of agreement.

The package accompanies the tutorial

> Zimmermann, J., Wessels, N. M., Rau, R., Freyer, N., & Leising, D. (in preparation). *Who agrees with whom about what? A tutorial on multilevel profile analysis of person judgments with the R package profilemlm.*

## Installation

```r
install.packages("remotes")
remotes::install_github("jzpsych/profilemlm", build_vignettes = TRUE)
```

The package requires lme4. ggplot2 and patchwork are needed for the figures, and brms (with cmdstanr or rstan) for Bayesian estimation. Building the vignettes needs knitr and rmarkdown and takes a few minutes. Without `build_vignettes = TRUE` the installation is faster, and the vignette sources can be read in the folder `vignettes/`.

## Getting started

```r
library(profilemlm)
vignette("mpa-tutorial", package = "profilemlm")
```

The tutorial vignette runs a complete analysis of a simulated study in five steps, with every call and its output. The example below is its beginning.

```r
path <- system.file("extdata", package = "profilemlm")
perceiver_survey <- read.csv(file.path(path, "perceiver_survey.csv"))
target_survey <- read.csv(file.path(path, "target_survey.csv"))
items <- read.csv(file.path(path, "items.csv"))

dat <- mpa_data(judgments = perceiver_survey, format = "loop", loop_id = "id",
                targets = target_survey, items = items, item_cols = "^item_",
                criterion = list(self = "targets"), group = "group",
                rating_scale = c(1, 101))

f1 <- mpa_formula(dat, criterion = "self")   # base model: normative and distinctive agreement
m1 <- mpa_fit(f1, dat)
m1
```

## Functions

| Function | What it does |
|---|---|
| `mpa_data()` | Reads judgments, criterion, and predictors as a study produces them, builds the normative and distinctive profiles, registers the level of each predictor, and detects the design |
| `mpa_formula()` | Builds the model formula from arguments |
| `mpa_update()` | Changes the random part of a formula |
| `mpa_fit()` | Estimates a model with lme4 or brms and returns tables of fixed effects and variance components |
| `mpa_prune()`, `mpa_random_test()` | Remove or test random terms |
| `mpa_compare()` | Compares two models (likelihood-ratio test, AIC, BIC, reduction of every variance component) |
| `mpa_effect_size()` | Standardized effect sizes for main and moderating effects |
| `mpa_plots()`, `mpa_plot()` | Agreement slopes as functions of the moderators |
| `mpa_simulate()`, `mpa_recovery()` | Simulate data from the model, for examples and for planning a study, and compare estimates with the generating values |

## Vignettes, scripts, and data

- `mpa-tutorial`: the five steps on simulated data, planning a study, further model components, Bayesian estimation.
- `predictor-profiles`: more than one predictor profile (assumed similarity), designs with unconnected pairs, effect sizes.
- `wessels-application`: replication of Model 2.2c of Wessels et al. (2020) with the data sets `wessels_judgments` and `wessels_items`.
- `system.file("scripts", package = "profilemlm")` contains a commented template for one's own data (`template_own_data.R`), the simulations reported in the tutorial (`sim_identifiability.R`, `sim_item_slope.R`), the figure script (`make_figures.R`), and the replication script (`wessels_study2.R`).

## Citation

`citation("profilemlm")` gives the reference. Please cite the tutorial when you use the package.

## License

GPL (>= 3). The data sets `wessels_judgments` and `wessels_items` were collected by Nele M. Wessels and colleagues (Wessels et al., 2020, 2025; https://osf.io/s5e6m).

## References

Biesanz, J. C. (2010). The social accuracy model of interpersonal perception: Assessing individual differences in perceptive and expressive accuracy. *Multivariate Behavioral Research, 45*(5), 853–885. https://doi.org/10.1080/00273171.2010.519262

Wessels, N. M., Zimmermann, J., Biesanz, J. C., & Leising, D. (2020). Differential associations of knowing and liking with accuracy and positivity bias in person perception. *Journal of Personality and Social Psychology, 118*(1), 149–171. https://doi.org/10.1037/pspp0000218
