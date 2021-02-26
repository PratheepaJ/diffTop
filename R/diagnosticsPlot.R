#' Diagnostics plots - ESS and Rhat.
#'
#' @param theta_aligned array. Three dimensional array returned by thetaAligned().
#' @param beta_aligned array. Three dimensional array returned by betaAligned().
#'
#' @return list. ggplots of effective sample size and rhat plots.
#' @export
#'
diagnosticsPlot <- function(
  theta_aligned,
  beta_aligned
){
  Rhat_theta <- matrix(
    nrow = dim(theta_aligned)[2],
    ncol = dim(theta_aligned)[3]
  )

  ESS_bulk_theta <- matrix(
    nrow = dim(theta_aligned)[2],
    ncol = dim(theta_aligned)[3]
  )


  for(sam in 1:dim(theta_aligned)[2]){
    for(top in 1:dim(theta_aligned)[3]){
      sims_theta <- matrix(
        theta_aligned[ ,sam , top],
        nrow = (iter/2),
        ncol = 4,
        byrow = FALSE
      )
      Rhat_theta[sam, top] <- Rhat(sims_theta)
      ESS_bulk_theta[sam, top] <- ess_bulk(sims_theta)
    }

  }

  Rhat_theta <- as.vector(Rhat_theta)
  ESS_bulk_theta <- as.vector(ESS_bulk_theta)


  Rhat_beta <- matrix(
    nrow = dim(beta_aligned)[2],
    ncol = dim(beta_aligned)[3]
  )
  ESS_bulk_beta <- matrix(
    nrow = dim(beta_aligned)[2],
    ncol = dim(beta_aligned)[3]
  )


  for(top in 1:dim(beta_aligned)[2]){
    for(fea in 1:dim(beta_aligned)[3]){
      sims_beta <- matrix(
        beta_aligned[ , top, fea],
        nrow = (iter/2),
        ncol = 4,
        byrow = FALSE)
      Rhat_beta[top, fea] <- Rhat(sims_beta)
      ESS_bulk_beta[top, fea] <- ess_bulk(sims_beta)

    }

  }

  Rhat_beta <- as.vector(Rhat_beta)
  ESS_bulk_beta <- as.vector(ESS_bulk_beta)


  Rhat <- c(Rhat_theta, Rhat_beta)


  ESS_bulk <- c(ESS_bulk_theta, ESS_bulk_beta)




  # R hat ~ 1.05
  p_rhat <- ggplot(
    data.frame(Rhat = Rhat)
  ) +
    geom_histogram(
      aes(x = Rhat),
      fill = "lavender",
      colour = "black",
      bins = 100
    ) +
    theme(
      plot.title = element_text(hjust = 0.5)
    )  +
    theme_minimal(base_size = 20) +
    xlab("")




  # ESS bulk and ESS tail at least 100 per Markov Chain in order to be reliable and indicate that estimates of respective posterior quantiles are reliable

  p_ess_bulk <- ggplot(
    data.frame(ESS_bulk = ESS_bulk)
  ) +
    geom_histogram(
      aes(x = ESS_bulk),
      fill = "lavender",
      colour = "black",
      bins = 100
    ) +
    theme(
      plot.title = element_text(hjust = 0.5)
    )   +
    theme_minimal(
      base_size = 20
    ) +
    xlab("")

  return(list(p_ess_bulk, p_rhat))
}

