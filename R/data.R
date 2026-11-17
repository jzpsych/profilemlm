#' Wessels et al. (2020, 2025), Study 2: informant judgments
#'
#' Informant descriptions of 189 targets on 111 person-descriptive items,
#' with the targets' self-ratings, from Study 2 of Wessels, Zimmermann,
#' Biesanz and Leising (2020), reanalysed by Wessels, Leising and Zimmermann
#' (2025). Each of 1,352 informants described one target (a one-with-many
#' design) on a 1--5 scale and reported how well they knew and how much they
#' liked the target. 943 informants were recruited from the targets' university
#' groups (`type = 0`), 409 were nominated by the targets (`type = 1`). The
#' 111 items combine established instruments (30 Big Five adjectives, the
#' short form of the Interpersonal Adjective List, the Rosenberg Self-Esteem
#' Scale, single items on self-esteem and narcissism, the PHQ-2, one item on
#' self-criticism) with 50 self-developed items on maladaptive traits; their
#' rated features are in [wessels_items].
#'
#' One row per informant x item (150,072 rows). The self-rating of the target
#' is repeated in every row of that target; take `unique(wessels_judgments[,
#' c("target", "item", "self")])` to obtain the criterion profiles (see the
#' example). The data are not used by the tutorial on simulated data; the
#' vignette `wessels-application` (`vignette("wessels-application", package =
#' "profilemlm")`) uses them to reproduce Model 2.2c of Wessels et al. (2020,
#' Table 5) exactly.
#'
#' @format A data frame with 150,072 rows and 8 variables:
#' \describe{
#'   \item{perceiver}{informant identifier (1--1352)}
#'   \item{target}{target identifier (1--189)}
#'   \item{item}{item identifier (1--111), matching `wessels_items$item`}
#'   \item{rating}{the informant's judgment of the target on the item (1--5)}
#'   \item{self}{the target's self-rating on the item (1--5)}
#'   \item{type}{informant type: 0 = recruited from the target's group, 1 = nominated by the target}
#'   \item{knowing}{how well the informant knows the target (1--5; mean of four items)}
#'   \item{liking}{how much the informant likes the target (1--5; mean of four items)}
#' }
#' @source Wessels, N. M., Zimmermann, J., Biesanz, J. C., & Leising, D. (2020).
#'   Differential associations of knowing and liking with accuracy and
#'   positivity bias in person perception. *Journal of Personality and Social
#'   Psychology, 118*(1), 149--171. \doi{10.1037/pspp0000218}. Data:
#'   <https://osf.io/s5e6m>. Reduced to the variables used in the replication and
#'   included with the consent of the first author, who is listed as data contributor of the package.
#' @seealso [wessels_items]
#' @examples
#' str(wessels_judgments)
#' self_long <- unique(wessels_judgments[, c("target", "item", "self")])
#' names(self_long)[3] <- "value"
#' dat <- mpa_data(judgments = wessels_judgments[, c("perceiver", "target", "item", "rating",
#'                                                  "knowing", "liking", "type")],
#'               format = "long", rating_col = "rating", items = wessels_items,
#'               criterion = list(self = self_long), rating_scale = c(1, 5),
#'               rescale = list(DES = c(1, 10), OBS = c(1, 10), IMP = c(1, 10),
#'                              ST = c(1, 10), BR = c(1, 10), knowing = c(1, 5), liking = c(1, 5)))
#' dat
"wessels_judgments"

#' Wessels et al. (2025), Study 2: rated item features
#'
#' Mean ratings of five features of the 111 items of [wessels_judgments] by 30
#' independent raters on a 1--10 scale (Wessels, Leising, & Zimmermann, 2025):
#' social desirability, observability, importance, stability and base rate.
#'
#' @format A data frame with 111 rows and 6 variables:
#' \describe{
#'   \item{item}{item identifier, matching `wessels_judgments$item`}
#'   \item{DES}{social desirability (1 = very undesirable, 10 = very desirable)}
#'   \item{OBS}{observability}
#'   \item{IMP}{importance}
#'   \item{ST}{stability}
#'   \item{BR}{base rate (how many people the item applies to)}
#' }
#' @source Wessels, N. M., Leising, D., & Zimmermann, J. (2025). A nuanced
#'   perspective on how item features are associated with different forms of
#'   agreement. *Collabra: Psychology, 11*(1), Article 127423.
#'   \doi{10.1525/collabra.127423}. Data: <https://osf.io/fxk7c>.
#' @seealso [wessels_judgments]
#' @examples
#' summary(wessels_items)
#' cor(wessels_items[, -1])
"wessels_items"
