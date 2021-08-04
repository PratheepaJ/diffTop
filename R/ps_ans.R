#' \code{'phyloseq'} object for the Anscombe transformed data published in Fitzpatrick et al. (2018) after removing DNA contamination.
#'
#' A \code{'phyloseq'} class object with 1418 taxa in 86 samples .
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
#' sample_data(ps_ans)
#' otu_table(ps_ans)
#' }
"ps_ans"
