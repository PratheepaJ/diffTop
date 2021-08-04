#' Differential topic analysis
#'
#' @param subsetSample Character vector. Sample names to be considered.
#' @inheritParams plotTopicProportion
#' @param ... Arguments to DESeq2::DeSeq.
#'
#' @return DESeq object.
#' @importFrom DESeq2 DESeq DESeqDataSetFromMatrix
#' @import phyloseq
#' @export
#'
#'

diffTopAnalysis <- function(
  design = ~ pna,
  ps,
  theta_aligned,
  subsetSample = NULL,
  ...
){

  median <- NULL
  if(is.null(subsetSample)){
    subsetSample = sample_names(ps)
  }
  # posterior samples
  theta_aligned_abun <- theta_aligned
  # median Bayesian posterior of topic proportions in each specimen
  median_theta_aligned_abun <- apply(
    theta_aligned_abun,
    c(2,3),
    median)

  # Multiply the proportions by the library size and round to an integer
  lib_size <- colSums(
    otu_table(ps) %>%
      data.frame()
  )

  median_theta_aligned_abun <- median_theta_aligned_abun * lib_size

  median_theta_aligned_abun <- ceiling(
    median_theta_aligned_abun
  )

  sam_data <- phyloseq::sample_data(ps) %>%
    data.frame()

  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = t(median_theta_aligned_abun),
    colData = sam_data,
    design = design
  )
  dds <- dds[, as.character(dds$X) %in% subsetSample]

  dds <- DESeq2::DESeq(
    dds,
    ...
  )

  return(dds)
}
