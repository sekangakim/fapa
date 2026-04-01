#' Synthetic EDI-2 profile data for FAPA examples
#'
#' A simulated dataset approximating the structure of the calibration sample
#' used in Kim (in preparation).  It contains no real clinical records.
#' The data comprise 500 synthetic cases on 22 variables: 11 pre-treatment
#' and 11 post-treatment administrations of the Eating Disorder Inventory-2
#' (EDI-2) subscales.
#'
#' @format A data frame with 500 rows and 22 columns.  The first 11 columns
#'   contain pre-treatment EDI-2 subscale scores (Drive for Thinness through
#'   Social Insecurity) and the remaining 11 columns contain the corresponding
#'   post-treatment scores.  Column names follow the convention
#'   \code{Before_<n>_<abbr>} and \code{After_<n>_<abbr>}, where \code{n} is
#'   the subscale index and \code{abbr} is a two-letter abbreviation.
#'   Scores are integers in the range 0--40.
#'
#' @details
#' The 11 EDI-2 subscales are:
#' \itemize{
#'   \item \strong{Dt} Drive for Thinness
#'   \item \strong{Bu} Bulimia
#'   \item \strong{Bd} Body Dissatisfaction
#'   \item \strong{In} Ineffectiveness
#'   \item \strong{Pf} Perfectionism
#'   \item \strong{Id} Interpersonal Distrust
#'   \item \strong{Ia} Interoceptive Awareness
#'   \item \strong{Mf} Maturity Fears
#'   \item \strong{As} Asceticism
#'   \item \strong{Ir} Impulse Regulation
#'   \item \strong{Si} Social Insecurity
#' }
#'
#' The latent structure was constructed to approximate two components: a
#' normative symptom gradient (CP1) and a pre-/post-treatment change contrast
#' (CP2).
#'
#' @source Simulated. See \code{data-raw/simulate_fapa_data.R}.
#'
#' @references
#' Garner, D. M. (1991). \emph{Eating Disorder Inventory-2 Manual}.
#' Psychological Assessment Resources.
#'
#' @examples
#' data(fapa_simdata)
#' dim(fapa_simdata)
#'
#' ## Quick ipsatization check
#' Xtilde <- as.matrix(fapa_simdata) - rowMeans(as.matrix(fapa_simdata))
#' range(rowSums(Xtilde))   # should be ~0
"fapa_simdata"
