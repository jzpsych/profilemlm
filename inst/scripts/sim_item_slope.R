###########################################################
### What happens when the item slope is left out? (article, Step 1)
### False-positive rate of an item-property moderator of distinctive agreement
### (true effect 0) in the design of the tutorial example, with and without the
### item slope, for an item-level SD of the distinctive slopes of .04 and .08.
### 400 data sets, two fits each; about 45 minutes on four cores.
###########################################################
library(profilemlm); library(parallel)
sdx <- function(item_dist) list(perceiver = c(intercept = 8, norm = 0.25, dist = 0.08),
                   target = c(intercept = 2, norm = 0.05, dist = 0.06),
                   dyad = c(intercept = 5, norm = 0.15, dist = 0.08),
                   item = c(intercept = 3, dist = item_dist))
one <- function(seed, item_dist) {
  dat <- mpa_simulate(n_perceiver = 48, n_item = 30, group_size = 8, norm_sd = 12, dist_sd = 10,
                      fixed = c(intercept = 51, norm = 0.6, dist = 0.15, des_dist = 0),
                      sd = sdx(item_dist), resid_sd = 12, rating_scale = c(1, 101), seed = seed)
  f <- mpa_formula(dat, criterion = "self", predictors = "des", interactions = "criterion", correlated = FALSE)
  f0 <- mpa_update(f, drop = list(item = "self.dist"))
  g <- function(ff) { m <- suppressWarnings(suppressMessages(mpa_fit(ff, dat))); i <- m$fixed$term == "self.dist:des"
    c(m$fixed$p[i], m$fixed$se[i]) }
  a <- g(f); b <- g(f0)
  c(seed = seed, item_dist = item_dist, p_with = a[1], se_with = a[2], p_without = b[1], se_without = b[2])
}
g <- expand.grid(s = 1:200, sd = c(.04, .08))
res <- mclapply(seq_len(nrow(g)), function(i) tryCatch(one(g$s[i], g$sd[i]), error = function(e) rep(NA, 6)), mc.cores = 4)
res <- as.data.frame(do.call(rbind, res))
names(res) <- c("seed", "item_sd", "p_with", "se_with", "p_without", "se_without")
saveRDS(res, "sim_item_slope.rds")
aggregate(cbind(false_pos_with = p_with < .05, false_pos_without = p_without < .05, se_with, se_without) ~ item_sd, res, mean)
