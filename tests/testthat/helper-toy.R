## Toy data set used by several test files
make_toy_mpa <- function(n_p = 20, n_t = 10, n_i = 8, seed = 42) {
  set.seed(seed)
  perceivers <- paste0("p", seq_len(n_p)); targets <- paste0("t", seq_len(n_t)); items <- paste0("item_", seq_len(n_i))
  norm_true <- runif(n_i, 30, 70)
  criterion <- matrix(norm_true, n_t, n_i, byrow = TRUE) + matrix(rnorm(n_t * n_i, 0, 10), n_t, n_i)
  dy <- expand.grid(perceiver = perceivers, target = targets, stringsAsFactors = FALSE)
  dy$know <- runif(nrow(dy), 1, 101)
  b_dist <- 0.4 + rnorm(n_p, 0, 0.1)[match(dy$perceiver, perceivers)]
  for (i in seq_len(n_i)) {
    dy[[items[i]]] <- 20 + 0.6 * (norm_true[i] - mean(norm_true)) +
      b_dist * (criterion[match(dy$target, targets), i] - norm_true[i]) + rnorm(nrow(dy), 0, 8)
  }
  dy$role <- "perceiver"
  self <- data.frame(perceiver = targets, target = targets, know = NA, stringsAsFactors = FALSE)
  for (i in seq_len(n_i)) self[[items[i]]] <- criterion[, i]
  self$role <- "self"
  mpa_data(rbind(dy, self), items = data.frame(item = items, des = seq(0, 100, length.out = n_i)),
           item_cols = "^item_", role = "role", judgment_role = "perceiver",
           criterion = list(self = "self"), norm = "loo",
           rescale = list(des = c(0, 100), know = c(1, 101)), verbose = FALSE)
}
