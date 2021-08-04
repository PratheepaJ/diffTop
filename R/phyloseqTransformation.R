#' Anscombe (1948) transformation
#'
#' @param cc Double. a constant used in the Anscombe (1948) transformation if dispersion is <= 2
#' @inheritParams alignmentMatrix
#' @return phyloseq object. otu table has transformed values of counts.
#' @importFrom DESeq2 estimateSizeFactors estimateDispersions sizeFactors dispersions
#' @importFrom phyloseq phyloseq phyloseq_to_deseq2
#' @importFrom SummarizedExperiment assay
#' @export
#'
phyloseqTransformation <- function(
  ps,
  cc = .4
  ){

  dds <- phyloseq::phyloseq_to_deseq2(
    ps,
    design = ~ 1
    )
  dds <-  DESeq2::estimateSizeFactors(
    dds,
    type = "poscounts"
    )
  dds <- DESeq2::estimateDispersions(
    dds,
    fitType = "local"
    )
  abund <- SummarizedExperiment::assay(dds)
  abund_temp <- matrix(
    nrow = nrow(abund),
    ncol = ncol(abund)
    )

  # Column co-factor

  for(i in 1:dim(abund)[1]){
    for(j in 1:dim(abund)[2]){
      abund_temp[i,j] <- abund[i,j]/sizeFactors(dds)[j]
    }
  }


  # Row co-factor

  # k > 2
  ind <- which(dispersions(dds) > 2)
  # K <= 2
  ind2 <- which(dispersions(dds) <= 2)

  for(i in ind){
    for(j in 1:dim(abund_temp)[2]){
      if(abund_temp[i,j] > 0){
        abund_temp[i,j] <- sqrt(dispersions(dds)[i]- .5)*asinh(sqrt((abund_temp[i, j]+ cc)/(dispersions(dds)[i]-2*cc)))
      }else{
        abund_temp[i,j] <- sqrt(dispersions(dds)[i]- .5)*asinh(sqrt(abund_temp[i, j]))
      }

    }
  }


  if(length(ind2) > 0){
    for(i in ind2){
      for(j in 1:dim(abund_temp)[2]){
        if(abund_temp[i,j] > 0){
          abund_temp[i,j] <- log2(abund_temp[i, j] + .5*dispersions(dds)[i])
        }

      }
    }
  }


  rownames(abund_temp) <- rownames(abund)
  colnames(abund_temp) <- colnames(abund)


  ps <- phyloseq::phyloseq(
    otu_table(abund_temp, taxa_are_rows = TRUE),
    sample_data(ps),
    tax_table(ps),
    phy_tree(ps)
    )

  return(ps)

}






