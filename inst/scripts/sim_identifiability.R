###########################################################
### Simulation study: can perceiver-level variance in distinctive agreement be
### separated from dyad-level variance? Design factor: targets per
### perceiver (half-block, group_size perceivers x targets_per_group targets).
###
### Requires the package to be INSTALLED (workers call library(profilemlm)):
###   devtools::install(build_vignettes = FALSE)
### Runtime: roughly n_rep x 12 cells x (5-40 s) / cores.
### Note: the data are generated with the random-effect SDs of the tutorial
### example (Table 4 of the paper), including item-level variance in
### distinctive agreement, and the fitted model is the base model with the
### item slope and uncorrelated random effects.
###########################################################

library(profilemlm)
library(parallel)

cores  <- 4
n_rep  <- 30                                    # run 1 (Sept 2026): 50 reps, group_size 4 only
grid <- expand.grid(targets_per_perceiver = c(4, 8, 16), perceivers_per_target = c(4, 8),
                    n_item = c(30, 60), rep = seq_len(n_rep))

## true SDs (per criterion point), as in the tutorial example:
## distinctive slopes perceiver 0.08, dyad 0.08, target 0.06, item 0.04
criterion_sd <- list(perceiver = c(intercept = 8, norm = 0.25, dist = 0.08),
                 target    = c(intercept = 2, norm = 0.05, dist = 0.06),
                 dyad      = c(intercept = 5, norm = 0.15, dist = 0.08),
                 item      = c(intercept = 3, dist = 0.04))

one_cell <- function(i) {
  g <- grid[i, ]
  dat <- suppressWarnings(mpa_simulate(n_perceiver = 64, group_size = g$perceivers_per_target,
                    targets_per_group = g$targets_per_perceiver,
                    n_item = g$n_item, seed = 1000 + i,
                    norm_sd = 12, dist_sd = 10,
                    fixed = c(intercept = 51, norm = 0.6, dist = 0.15),
                    sd = criterion_sd, resid_sd = 12, rating_scale = c(1, 101)))
  f <- mpa_formula(dat, criterion = "self", correlated = FALSE)   # base model with the item slope
  m <- tryCatch(suppressWarnings(mpa_fit(f, dat, check_derivs = FALSE)), error = function(e) NULL)
  if (is.null(m)) return(NULL)
  v <- m$variance
  get <- function(u, e) { x <- v$sd[v$unit == u & v$effect == e]; if (length(x)) x else NA }
  data.frame(g, n_target = attr(dat, "mpa")$design$n_target,
             sd_perc_dist = get("perceiver", "self.dist"), sd_rel_dist = get("perceiver:target", "self.dist"),
             sd_targ_dist = get("target", "self.dist"), sd_item_dist = get("item", "self.dist"),
             sd_perc_norm = get("perceiver", "self.norm.c"), sd_rel_norm = get("perceiver:target", "self.norm.c"),
             b_dist = m$fixed$estimate[m$fixed$term == "self.dist"],
             singular = m$singular, seconds = NA)
}

cl <- makeCluster(cores)
clusterEvalQ(cl, library(profilemlm))
clusterExport(cl, c("grid", "criterion_sd", "one_cell"))
t0 <- Sys.time()
res <- do.call(rbind, parLapply(cl, seq_len(nrow(grid)), one_cell))
stopCluster(cl)
cat("Elapsed:", format(Sys.time() - t0), "\n")
saveRDS(res, "sim_identifiability.rds")

## ---- summary: mean, SD and share of boundary (zero) estimates per cell ------
agg <- aggregate(cbind(sd_perc_dist, sd_rel_dist, sd_targ_dist, sd_item_dist) ~ targets_per_perceiver + perceivers_per_target + n_item, data = res,
                 FUN = function(x) c(mean = mean(x), sd = sd(x), zero = mean(x < 1e-3)))
print(agg, digits = 3)
cat("\nTrue values: perceiver 0.08, dyad 0.08, target 0.06, item 0.04\n")

## ---- figure -----------------------------------------------------------------
if (requireNamespace("ggplot2", quietly = TRUE)) {
  library(ggplot2)
  cols <- c("targets_per_perceiver", "perceivers_per_target", "n_item")
  long <- rbind(
    data.frame(res[, cols], component = "perceiver",    sd = res$sd_perc_dist, true = 0.08),
    data.frame(res[, cols], component = "dyad",         sd = res$sd_rel_dist,  true = 0.08),
    data.frame(res[, cols], component = "target",       sd = res$sd_targ_dist, true = 0.06),
    data.frame(res[, cols], component = "item",         sd = res$sd_item_dist, true = 0.04))
  p <- ggplot(long, aes(factor(targets_per_perceiver), sd, fill = factor(perceivers_per_target))) +
    geom_boxplot(outlier.size = 0.5, position = position_dodge(0.8)) +
    geom_hline(aes(yintercept = true), linetype = 2) +
    facet_grid(component ~ n_item, labeller = label_both) +
    scale_fill_grey(start = 0.6, end = 0.95, name = "Perceivers per target") +
    labs(x = "Targets per perceiver", y = "Estimated SD of distinctive-agreement slopes (per criterion point)") +
    theme_bw()
  print(p)
  ggsave("sim_identifiability.png", p, width = 8, height = 9, dpi = 300)
}
