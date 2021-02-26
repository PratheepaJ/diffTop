#' Model assessement for LDA
#'
#' @param stan.fit An instance of stanfit.
#' @param statistic Function name. Use this statistic to assess the goodness of fit.
#' @param ASVsIndexToPlot An integer vector. Use to select ASVs to show the goodness of fit in histograms.
#' @inheritParams alignmentMatrix
#' @return
#' @export
#'
modelAssessment <- function(
  ps,
  stan.fit,
  iter,
  statistic = "max",
  ASVsIndexToPlot = c(1, 3, 10:14, 19:26, 36, 51:53, 148)
){
  samples <- rstan::extract(
    stan.fit,
    permuted = TRUE,
    inc_warmup = FALSE,
    include = TRUE)# samples is a list

  x <- get_taxa(ps) %>%
    t()

  dimnames(x) <- NULL
  x_max_asv <- apply(
    x %>% data.frame ,
    2,
    max)

  x_max_asv <- data.frame(
    x_max_asv = x_max_asv
  )

  # draws from posterior predictive distribution
  x_sim <- samples$x_sim # iteration * samples * ASVs
  # Choose only the first chain
  x_sim <- x_sim[1:(iter/2), ,] # For each iteration, simulated data is x_sim[i, ,]
  sim_ite <- dim(x_sim)[1]

  #Find maximum of each asv for each replication
  max_all <- data.frame()
  x_sim_i <- x_sim[1, ,]
  max_x_sim_i <- apply(
    x_sim_i %>% data.frame ,
    2,
    max
  )
  max_all <- data.frame(max_x_sim_i)

  for(i in 1:sim_ite){
    x_sim_i <- x_sim[i, ,]
    max_x_sim_i <- apply(
      x_sim_i %>% data.frame ,
      2,
      max)
    max_all <- cbind(max_all, max_x_sim_i)
  }

  colnames(max_all) <- c(
    paste0(
      "x_max_rep",
      seq(1, sim_ite)
    )
  )

  rownames(max_all) <- paste0(
    "ASV_",
    seq(1, ntaxa(ps)),
    "_",
    tax_table(ps)[,"Class"]
  )

  max_all <- max_all %>% t()
  ## choose some of asvs for Posterior predictive check (these indexes are considered to avoid if ASV's Class is missing)
  max_all <- max_all[, ASVsIndexToPlot]
  max_all_long <- melt(max_all)

  # observed stat data frame
  x_max_asv <- data.frame(
    Var1 = rep(
      "x_max_obs",
      dim(x_max_asv)[1]
    ),
    Var2 = paste0(
      "ASV_",
      seq(1, ntaxa(ps)),
      "_",
      tax_table(ps)[,"Class"]
    ),
    value = x_max_asv$x_max_asv
  )

  x_max_asv <- x_max_asv[ASVsIndexToPlot,]


  p_hist <- ggplot(
    data = max_all_long
  ) +
    geom_histogram(
      aes(
        x = value,
        group = Var2
      ),
      color = "blue",
      fill = "blue",
      bins = 50) +
    xlab("maximum")+
    facet_wrap(~Var2, nrow = 4) +
    geom_vline(
      data = x_max_asv,
      aes(xintercept = value),
      color = "purple"
    ) +
    theme_update(
      text = element_text(size = 8)
    )
  return(p_hist)

}
