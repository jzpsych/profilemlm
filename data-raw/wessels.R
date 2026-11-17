## Build the package data sets `wessels_judgments` and `wessels_items` from the
## files of Wessels et al. (2020, JPSP; 2025, Collabra), Study 2.
##   datatwo_selfasval.csv            informant-level file (OSF s5e6m; reduced to the columns used)
##   S2_data_raw_agreement_itemwise.csv  item features (OSF fxk7c)
## Run from the package root: source("data-raw/wessels.R")
path <- "data-raw/wessels_src/"   # adapt
w <- read.csv2(file.path(path, "datatwo_selfasval.csv"), fileEncoding = "UTF-8-BOM")
names(w)[1] <- "item"
wessels_judgments <- data.frame(
  perceiver = as.integer(w$perceiver), target = as.integer(w$target), item = as.integer(w$item),
  rating = as.numeric(w$info), self = as.numeric(w$self),
  type = as.integer(w$info_type), knowing = as.numeric(w$knowing), liking = as.numeric(w$liking))
wessels_judgments <- wessels_judgments[order(wessels_judgments$perceiver, wessels_judgments$item), ]
rownames(wessels_judgments) <- NULL
stopifnot(nrow(wessels_judgments) == 1352 * 111, !anyDuplicated(wessels_judgments[, c("perceiver", "item")]))

feat <- read.csv2(file.path(path, "S2_data_raw_agreement_itemwise.csv"), fileEncoding = "UTF-8-BOM")
names(feat)[1] <- "item"
wessels_items <- unique(feat[, c("item", "DES", "OBS", "IMP", "ST", "BR")])
wessels_items[] <- lapply(wessels_items, as.numeric); wessels_items$item <- as.integer(wessels_items$item)
wessels_items <- wessels_items[order(wessels_items$item), ]; rownames(wessels_items) <- NULL
stopifnot(nrow(wessels_items) == 111)
## the self-ratings and judgments in the item-feature file equal those in the informant file (checked)
save(wessels_judgments, file = "data/wessels_judgments.rda", compress = "xz")
save(wessels_items, file = "data/wessels_items.rda", compress = "xz")
