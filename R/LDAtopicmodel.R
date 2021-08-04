#' Bayesian latent Dirichlet allocation for 16S and metagenomics data.
#'
#' @param stan_data A list of stan data.
#' @param ... Arguments passed to `rstan::sampling` (e.g. iter, chains).
#'
#' @return An object of class `stanfit` returned by `rstan::sampling`.
#' @importFrom rstan sampling
#' @export
LDAtopicmodel <- function(stan_data, ...){
        stan.fit <- rstan::sampling(
          object = stanmodels$lda,
          data = stan_data,
          ...
        )
        return(stan.fit)
}


