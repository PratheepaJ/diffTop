#' \code{'DESeqDataSet'} object for the differential topic analysis on all samples in psE_BARBI.
#'
#' A By applying DEseq2::results() on \code{'DESeqDataSet'}, we get a data frame with topics on rows and test results on the columns (baseMean, log2FoldChange, lfcSE, stat, pvalue, padj)
#'
#' @format An object of class \code{'DESeqDataSet'}.
#' \describe{
#'   \item{baseMean}{mean of the topic abundance}
#'   \item{log2FoldChange}{log fold change of topic abundance in experimental condions tested.}
#'   \item{lfcSE}{Standard error of log fold change.}
#'   \item{stat}{Test statistic value under the null hypothesis.}
#'   \item{pvalue}{P-value for the test.}
#'   \item{padj}{Adjusted P-value.}
#'   }
#' @source Data is from Fitzpatrick et al. (2018)
#' @examples
#' \dontrun{
#' results(dds_all_K_11)
#' }
"dds_all_K_11"
