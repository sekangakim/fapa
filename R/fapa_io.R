# =============================================================================
# fapa_io.R
# CSV output helpers for core profiles and verification results.
# =============================================================================


#' Write core-profile BCa CI tables to CSV
#'
#' Writes one CSV file per retained core profile, containing the original
#' coordinates together with bootstrap mean, SE, percentile, and BCa
#' confidence bounds.
#'
#' @param bca    A list returned by \code{\link{fapa_bca}}.
#' @param prefix Character. Base name for output files (e.g., \code{"fapa"}).
#'   Files are named \code{<prefix>_CoreProfile1.csv}, etc.
#'
#' @return Invisibly returns a character vector of file paths written.
#' @export
write_fapa_results <- function(bca, prefix) {
  paths <- character(bca$K)
  for (k in seq_len(bca$K)) {
    fname <- paste0(prefix, "_CoreProfile", k, ".csv")
    df    <- bca$ci[[k]]
    rownames(df) <- bca$varnames
    write.csv(df, file = fname)
    message("Wrote: ", fname)
    paths[k] <- fname
  }
  invisible(paths)
}


#' Write three-stage verification results to CSV
#'
#' Writes one CSV file for each of the three verification stages.
#'
#' @param pa     A list returned by \code{\link{fapa_pa}}.
#' @param pr     A list returned by \code{\link{fapa_procrustes}}.
#' @param tc     A list returned by \code{\link{fapa_tucker}}.
#' @param prefix Character. Base name for output files.
#' @param K_pa   Integer or \code{NULL}. If supplied, a \code{PA_Retained}
#'   column is appended to the Stage 2 and Stage 3 CSVs.
#'
#' @return Invisibly returns a named character vector of the three file paths.
#' @export
write_verification <- function(pa, pr, tc, prefix, K_pa = NULL) {

  ## Stage 1
  f1 <- paste0(prefix, "_Stage1_PA.csv")
  write.csv(
    data.frame(Component    = seq_along(pa$obs_sv2),
               Obs_sv2      = pa$obs_sv2,
               Obs_PropVar  = pa$prop_obs,
               Rand_95pct   = pa$thresh,
               Rand_PropVar = pa$prop_rand,
               Retained     = pa$obs_sv2 > pa$thresh),
    file = f1, row.names = FALSE)

  ## Stage 2
  f2  <- paste0(prefix, "_Stage2_Procrustes.csv")
  df2 <- data.frame(Dimension   = seq_len(pr$K),
                    Mean_deg    = pr$angle_mean,
                    SD_deg      = pr$angle_sd,
                    CI_lo       = pr$angle_q025,
                    CI_hi       = pr$angle_q975,
                    Angle_OK    = pr$angle_mean < pr$angle_thresh,
                    Prop_stable = pr$prop_stable)
  if (!is.null(K_pa))
    df2$PA_Retained <- seq_len(pr$K) <= K_pa
  write.csv(df2, file = f2, row.names = FALSE)

  ## Stage 3
  f3  <- paste0(prefix, "_Stage3_TuckerCC.csv")
  df3 <- data.frame(Core_Profile = seq_len(tc$K),
                    Mean_CC      = tc$cc_mean,
                    SD_CC        = tc$cc_sd,
                    CI_lo        = tc$cc_q025,
                    CI_hi        = tc$cc_q975,
                    CC_OK        = tc$cc_mean >= tc$cc_thresh)
  if (!is.null(K_pa))
    df3$PA_Retained <- seq_len(tc$K) <= K_pa
  write.csv(df3, file = f3, row.names = FALSE)

  message("Verification results written to: ", f1, ", ", f2, ", ", f3)
  invisible(c(Stage1 = f1, Stage2 = f2, Stage3 = f3))
}
