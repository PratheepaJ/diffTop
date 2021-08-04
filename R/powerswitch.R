#' A vector
#'
#' Contais the number of pvalues (out of 1000) less than 0.05 and 0.01 for the PERMANOVA in the simulated data with strain switching.
#'
#'
#' @format A vector of integers.
#' \describe{
#'   \item{first element}{number of pvalues (out of 1000) less than 0.05.}
#'   \item{Second element}{number of pvalues (out of 1000) less than 0.01.}
#'   }
#' @source Simulated data based on the data from Fitzpatrick et al. (2018)
#' @examples
#' \dontrun{
#' print(round(powerswitch[1]/1000,2))
#' }
"powerswitch"
