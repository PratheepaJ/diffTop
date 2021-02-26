#' Differential topic analysis
#'
#' @param subsetSample Character vector. Sample names to be considered.
#' @inheritParams plotTopicProportion
#'
#' @return DESeq object.
#' @export
#'
diffTop <- function(
  design = ~ pna,
  ps,
  theta_aligned,
  subsetSample = sample_names(psE_BARBI)
){
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

  sam_data <- sample_data(ps) %>%
    data.frame()

  dds <- DESeqDataSetFromMatrix(
    countData = t(median_theta_aligned_abun),
    colData = sam_data,
    design= ~ pna
  )
  dds <- dds[, as.character(dds$X) %in% subsetSample]

  dds <- DESeq(
    dds,
    fitType = "mean"
  )

  return(dds)
}
