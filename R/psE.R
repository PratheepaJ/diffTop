#' \code{'phyloseq'} object for the data published in Fitzpatrick et al. (2018).
#'
#' A \code{'phyloseq'} class object with 6929 taxa in 111 samples that are from seven rounds of six-fold dilutions (1:1 up to 1:279,936) from standard ZymoBIOMICS microbial community and 10 negative extraction controls.
#'
#' @format An object of class \code{'phyloseq'}.
#' \describe{
#'   \item{SampleType}{Whether Asteraceae or non-Asteraceae}
#'   \item{pna}{Type of pPNA used universal or Asteraceae-modified}
#'   \item{species_names}{Sample ID}
#'   }
#' @source Fitzpatrick et al. (2018)
#' @examples
#' \dontrun{
#' sample_data(psE)
#' otu_table(psE)
#' }
"psE"
