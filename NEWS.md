# profilemlm 1.0.1

* The tutorial is now available as a preprint (PsyArXiv, https://doi.org/10.31234/osf.io/5czt9_v1). `citation("profilemlm")`, the README, and the vignettes cite it.
* Removed the script `sim_item_slope.R` and its mention in the tutorial vignette.
* Reference lists of the vignettes and the README aligned with the article (issue numbers and DOIs added).

# profilemlm 1.0.0

First public release, accompanying the tutorial by Zimmermann, Wessels, Rau, Freyer, and Leising, *Who agrees with whom about what? A tutorial on multilevel profile analysis of person judgments with the R package profilemlm*.

* Data preparation from survey exports (`mpa_data()`), with normative and distinctive profiles, design detection, and the decomposition of dyadic variables.
* Model formulas from arguments (`mpa_formula()`, `mpa_update()`), with items as a crossed random factor and random slopes wherever a term varies within a unit.
* Estimation with lme4 or brms (`mpa_fit()`), with the same output tables for both engines.
* Random-effects structure (`mpa_prune()`, `mpa_random_test()`), model comparison (`mpa_compare()`), standardized effect sizes (`mpa_effect_size()`), and figures (`mpa_plots()`, `mpa_plot()`).
* Simulation for examples and for planning a study (`mpa_simulate()`, `mpa_recovery()`).
* Further predictor profiles, designs with unconnected pairs, and target × item and perceiver × item components.
* Vignettes `mpa-tutorial`, `predictor-profiles`, and `wessels-application`, and the data of Wessels et al. (2020, 2025), Study 2.
