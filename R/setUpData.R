#' Set-up data for latent Dirichlet allocation
#'
#' @param ps phyloseq class.
#' @param K  Integer. Number of topics.
#' @param alpha Non-negative numeric. Hyper-parameter alpha is less than one to generate sparse mixtures that are different from each other and avoid generating unrealistic topics.
#' @param gamma Non-negative numeric. Hyper-parameter gamma is less than one to generate sparse mixtures of ASVs in each topic.
#'
#' @return List. List of object with number of topics, number of ASVs,number of specimens, count table, hyper parameters for topic and ASVs distribution.
#' @export
#'
setUpData <- function(
  ps,
  K,
  alpha = 0.9,
  gamma = 0.5
){
  x <- t(get_taxa(ps))
  dimnames(x) <- NULL
  stan.data <- list(K = K,
                    V = ncol(x),
                    D = nrow(x),
                    n = x,
                    alpha = rep(alpha, K),
                    gamma = rep(gamma, ncol(x))
  )
  return(stan.data)
}
