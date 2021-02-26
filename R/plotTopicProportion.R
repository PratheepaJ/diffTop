#' Plot the topic distribution in all specimens.
#'
#' @param theta_aligned List. returned object of thetaAligned().
#' @param col_names_theta_all Character vector. Column names for long formatted theta. An example is c("iteration", "Sample", "Topic", "topic.dis")
#' @param warmup Integer. Number of iterations in each chain used for warmup samples.
#' @param SampleIDVarInphyloseq Character. SampleID variable name in phyloseq object.
#' @param design Factors of interests.
#' @inheritParams alignmentMatrix
#'
#' @return
#' @export
#'

plotTopicProportion <- function(
  ps,
  theta_aligned,
  K,
  col_names_theta_all = c("iteration", "Sample", "Topic", "topic.dis"),
  chain = 4,
  warmup = iter/2,
  SampleIDVarInphyloseq = "unique_names",
  design = ~pna
){
  dimnames(theta_aligned)[[2]] <- sample_names(ps)
  dimnames(theta_aligned)[[3]] <- c(paste0("Topic_", seq(1,K)))

  # array to a dataframe
  theta_all <- reshape2::melt(theta_aligned)

  colnames(theta_all) <- col_names_theta_all

  theta_all$Chain <- paste0(
    "Chain ",
    rep(seq(1, chain),
        each = warmup
    )
  )

  sam <- sample_data(ps) %>%
    data.frame()

  theta_all$Sample <- theta_all$Sample %>%
    as.character()

  theta_all <- left_join(
    theta_all,
    sam,
    by =c("Sample" = "unique_names")
  )

  theta_all$Chain <- factor(theta_all$Chain)
  theta_all$Topic <- factor(theta_all$Topic)
  theta_all$Sample <- factor(theta_all$Sample)
  theta_all$pna <- factor(theta_all$pna)


  theta_summary <- theta_all %>%
    group_by(
      Sample,
      Topic,
      pna
    ) %>%
    summarize(
      median.topic.dis = median(topic.dis)
    ) %>%
    ungroup() %>%
    mutate(
      Topic = factor(
        Topic,
        levels = rev(str_c("Topic_",1:K)
        )
      )
    )

  p <- ggplot(
    theta_summary,
    aes(x = pna,
        y = Topic,
        fill = pna)
  )
  p <- p+
    geom_tile(
      aes(alpha = median.topic.dis)
    ) +
    facet_grid(
      .~Sample,
      scale = "free"
    )+
    xlab("pna") +
    scale_fill_manual(
      name = "pna",
      values = c("steelblue1","green3")
    ) +
    scale_alpha(
      name = "median topic distribution"
    ) +
    theme_minimal(
      base_size = 20
    ) +
    theme(
      plot.title = element_text(hjust = 0.5),
      strip.text.x = element_text(angle = 90)
    )
  p
}
