###########################################################
### Wessels et al. (2020, JPSP; 2025, Collabra), Study 2
### Replication of Model 2.2c (Table 5, Wessels et al., 2020) with profilemlm
###
### Data: the package data sets `wessels_judgments` (one row per informant x
### item; 1352 informants, 189 targets, 111 items; judgment `rating`,
### self-ratings `self`, knowing, liking, type (0 = group informant, 1 =
### target-nominated)) and `wessels_items` (item features DES OBS IMP ST BR,
### 1-10). See ?wessels_judgments.
###
### Runtime: five to ten minutes, depending on the computer (one model on
### 150,072 rows). The fit and the
### comparison table are written to `out_dir`. The results are also shipped
### with the package and shown, with commentary, in the vignette
### "wessels-application" (vignette("wessels-application", package = "profilemlm")).
###########################################################

library(profilemlm)
out_dir <- "wessels_output"                      # adapt
dir.create(out_dir, showWarnings = FALSE)

## ---- data --------------------------------------------------------------------
self_long <- unique(wessels_judgments[, c("target", "item", "self")])   # the criterion: one row per target x item
names(self_long)[3] <- "value"
wj <- wessels_judgments[, c("perceiver", "target", "item", "rating", "knowing", "liking", "type")]   # the judgments

## ---- replication of Model 2.2c ----------------------------------------------
## = sample mean; desirability linearly rescaled to the 1-5 range of the ratings and,
## like knowing and liking, grand-mean centered; random intercepts and slopes (dist,
## norm, DES) for perceivers and targets, correlated; no item random effects; REML.

dA <- mpa_data(judgments = wj, format = "long", rating_col = "rating",
               items = wessels_items, criterion = list(self = self_long),
               norm = "sample", rating_scale = c(1, 5),
               derive = list(DES.c = ~ (DES - mean(DES)) * 4/9,    # desirability rescaled from 1-10 to the 1-5 range of the ratings, as in the paper
                             knowing.c = ~ knowing - mean(knowing),
                             liking.c = ~ liking - mean(liking)))
dA        # one-with-many: dyadic random effects off; type, liking, knowing registered as perceiver-level

## Their model has no DES x criterion terms, only the cross-level products of the
## three perceiver variables with intercept, dist, norm and DES -> interactions = "none"
## plus the nine products as extra_terms.
mods <- c("knowing.c", "liking.c", "type")
fA <- mpa_formula(dA, criterion = "self", random_slopes = "DES.c", predictors = mods, interactions = "none",
                  extra_terms = as.vector(outer(c("self.dist", "self.norm.c", "DES.c"), mods, paste, sep = ":")),
                  random_units = c("perceiver", "target"), correlated = TRUE)   # the classical SAM: no item random effects
fA
## The optimizer matters here: with lme4's default optimizer (nloptwrap) the fit
## is singular under some lme4 versions (the target intercept variance goes to
## zero) and misses the published values in the third decimal; bobyqa, the
## profilemlm default, reproduces them exactly
mA <- mpa_fit(fA, dA, REML = TRUE)
mA

## Reference values, Wessels et al. (2020) Table 5, Model 2.2c (checked against the printed table):
##  (Intercept) 2.567  self.dist 0.134  self.norm.c 0.572  DES 0.294
##  knowing 0.015  knowing:dist 0.026  knowing:norm 0.122  knowing:DES -0.098
##  liking -0.115  liking:dist -0.019  liking:norm -0.090  liking:DES  0.301
##  type    0.052  type:dist   0.074   type:norm    0.136  type:DES   -0.131
##  SDs perceiver: int .190 dist .117 norm .274 DES .282; target: int .059 dist .112 norm .219 DES .153; residual .853
ref <- c("(Intercept)" = 2.567, "self.dist" = 0.134, "self.norm.c" = 0.572, "DES.c" = 0.294,
         "knowing.c" = 0.015, "self.dist:knowing.c" = 0.026, "self.norm.c:knowing.c" = 0.122, "DES.c:knowing.c" = -0.098,
         "liking.c" = -0.115, "self.dist:liking.c" = -0.019, "self.norm.c:liking.c" = -0.090, "DES.c:liking.c" = 0.301,
         "type" = 0.052, "self.dist:type" = 0.074, "self.norm.c:type" = 0.136, "DES.c:type" = -0.131)
key <- function(x) vapply(strsplit(x, ":", fixed = TRUE), function(z) paste(sort(z), collapse = ":"), character(1))
cmp <- data.frame(term = mA$fixed$term, profilemlm = round(mA$fixed$estimate, 3),
                  wessels2020 = ref[match(key(mA$fixed$term), key(names(ref)))])
print(cmp)
print(mA$variance)
saveRDS(list(fixed = mA$fixed, variance = mA$variance, cmp = cmp), file.path(out_dir, "mA_replication.rds"))
saveRDS(mA, file.path(out_dir, "mA_fit.rds"))

## ---- reduced results as shipped with the package ------------------------------
## (inst/extdata/wessels_results.rds, shown in the vignette "wessels-application")
saveRDS(list(session = list(R = R.version.string,
                            lme4 = as.character(utils::packageVersion("lme4")),
                            profilemlm = as.character(utils::packageVersion("profilemlm")),
                            date = as.character(Sys.Date())),
             mA = list(fixed = mA$fixed, variance = mA$variance, correlations = mA$correlations,
                       singular = mA$singular, cmp = cmp),
             formulas = list(fA = utils::capture.output(print(fA)))),
        file.path(out_dir, "wessels_results.rds"))
