# =============================================================================
# FAPA-package.R
# Package-level documentation, imports, and global declarations.
# =============================================================================

#' FAPA: Factor Analytic Profile Analysis of Ipsatized Data
#'
#' @description
#' Implements Factor Analytic Profile Analysis of Ipsatized Data (FAPA), a
#' metric inferential framework for pattern detection and person-level
#' reconstruction in multivariate profile data.
#'
#' After row-centering (ipsatization) to remove profile elevation, FAPA
#' applies singular value decomposition (SVD) to recover shared \emph{core
#' profiles} and individual pattern weights, supporting the following
#' workflow:
#'
#' \enumerate{
#'   \item \strong{Ipsatization} --- \code{\link{load_and_ipsatize}} removes
#'     person-level elevation, isolating within-person pattern structure.
#'   \item \strong{Core estimation} --- \code{\link{fapa_core}} performs SVD
#'     and returns the core-profile matrix, person weights, and
#'     variance decomposition.
#'   \item \strong{Stage 1: Dimensionality} --- \code{\link{fapa_pa}}
#'     applies variance-matched Horn's parallel analysis to determine
#'     how many components to retain.
#'   \item \strong{Stage 2: Subspace stability} ---
#'     \code{\link{fapa_procrustes}} assesses dimensional stability via
#'     bootstrap principal angles.
#'   \item \strong{Stage 3: Profile replicability} ---
#'     \code{\link{fapa_tucker}} computes Tucker's congruence coefficients
#'     across bootstrap resamples.
#'   \item \strong{Inference} --- \code{\link{fapa_bca}} provides BCa
#'     bootstrap confidence intervals for core-profile coordinates using
#'     the canonical \pkg{boot} implementation.
#'   \item \strong{Reconstruction} --- \code{\link{fapa_person}} projects
#'     each person onto the retained core profiles and reports
#'     reconstruction \eqn{R^2} and optional bootstrap CIs for selected
#'     participants.
#' }
#'
#' @section Key references:
#' \itemize{
#'   \item Horn, J. L. (1965). A rationale and test for the number of
#'     factors in factor analysis. \emph{Psychometrika, 30}(2), 179--185.
#'   \item Davison, A. C., & Hinkley, D. V. (1997).
#'     \emph{Bootstrap Methods and Their Application}.
#'     Cambridge University Press.
#'   \item Lorenzo-Seva, U., & ten Berge, J. M. F. (2006). Tucker's
#'     congruence coefficient as a meaningful index of factor similarity.
#'     \emph{Methodology, 2}(2), 57--64.
#'   \item Kim, S.-K. (2024). Factorization of person response profiles
#'     to identify summative profiles carrying central response patterns.
#'     \emph{Psychological Methods, 29}(4), 723--730.
#'     \doi{10.1037/met0000568}
#' }
#'
#' @author Se-Kang Kim \email{se-kang.kim@@bcm.edu}
#'
#' @importFrom graphics abline axis hist legend lines par
#' @importFrom stats quantile rnorm sd
#' @importFrom utils read.csv write.csv
#'
#' @docType package
#' @name FAPA-package
#' @aliases FAPA
"_PACKAGE"
