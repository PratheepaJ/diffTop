#' Plot the ASVs distribution in each topic.
#'
#' @param beta_aligned List. returned object of betaAligned().
#' @param varnames Character vector. Use to name column of long formatted data frame of ASVs distribution. An example is c("iterations", "topic", "rsv_ix")
#' @param value_name Character. Use to name the summary of posterior samples of ASVs proportions.
#' @param taxonomylevel Character. Use to label ASVs.
#' @param thresholdASVprop Numberic. Use to select ASVs for circular plot.
#' @inheritParams alignmentMatrix
#' @return
#' @export
#'

plotASVCircular <- function(
  ps,
  K,
  iter = 2000,
  chain = 4,
  beta_aligned,
  varnames = c("iterations", "topic", "rsv_ix"),
  value_name = "beta_h",
  taxonomylevel = "Class",
  thresholdASVprop = 0.008
){
  # array to data frame
  beta_hat <- beta_aligned %>%
    reshape2::melt(
      varnames = c("iterations",
                   "topic",
                   "rsv_ix"),
      value.name = "beta_h") %>%
    as_tibble()

  beta_hat$rsv <- rownames(
    tax_table(ps)
  )[beta_hat$rsv_ix]

  # join taxa_table with beta_hat

  ## If we use a taxonomy level with NA,
  ## we can replace the taxonomy level with
  ## one level before this level
  #taxa$Class[which(is.na(taxa$Class))] = taxa$Phylum[which(is.na(taxa$Class))]
  taxa <- tax_table(ps) %>%
    as.data.frame() %>%
    as_tibble()

  taxa$rsv <- rownames(
    tax_table(ps)
  )

  beta_hat <- beta_hat %>%
    left_join(
      taxa,
      by = "rsv"
    ) %>%
    mutate(
      topic = paste("Topic", topic)
    )

  # Filtering ASVs for visualization:
  ## We can filter by ASVs proportion in all topics.

  beta_hat$Class <- factor(beta_hat$Class)
  beta_hat$rsv <- factor(beta_hat$rsv)
  beta_hat$rsv_ix <- factor(beta_hat$rsv_ix)
  beta_hat$topic <- factor(beta_hat$topic)

  # We plot the median of posterior sample and keep all taxonomy for each ASV.
  beta_summary <- beta_hat %>%
    dplyr::group_by(
      rsv_ix,
      topic
    ) %>%
    dplyr::summarise(
      beta_median = median(beta_h),
      rsv = rsv[1],
      Phylum = Phylum[1],
      Class = Class[1],
      Order = Order[1],
      Family = Family[1],
      Genus = Genus[1]
    )

  usethis::use_data(beta_summary)

  #######
  beta_subset <- beta_summary

  beta_subset$rsv_ix <- rep(
    seq_len(
      nrow(beta_subset) / K
    ),
    each = K
  )


  beta_subset <- beta_subset %>%
    arrange(
      rsv_ix,
      topic
    )

  beta_subset <- beta_subset %>%
    mutate(
      Class = factor(
        Class,
        levels = unique(beta_subset$Class)
      ),
      Topic = str_remove(topic, "Topic ")
    )


  beta_subset$Class <- beta_subset$Class %>%
    as.character()

  beta_subset$Class <- ifelse(
    is.na(beta_subset$Class),
    "other",
    beta_subset$Class
  )

  beta_subset$Class <- factor(
    beta_subset$Class
  )


  beta_subset$Topic <- factor(
    beta_subset$Topic,
    levels = seq(1,K) %>%
      as.character()
  )

  ####

  beta_subset_wide <- dplyr::select(
    beta_subset,
    rsv_ix,
    rsv,
    topic,
    Class,
    beta_median
  ) %>%
    as.data.frame()

  beta_subset_wide <- dplyr::select(
    beta_subset_wide,
    -rsv_ix
  )

  beta_subset_wide <- tidyr::spread(
    beta_subset_wide,
    key = topic,
    value = beta_median
  )

  rownames(beta_subset_wide) <- beta_subset_wide$rsv

  beta_subset_wide$Class <- beta_subset_wide$Class %>%
    as.character()

  beta_subset_wide$Class <- ifelse(
    is.na(beta_subset_wide$Class),
    "other",
    beta_subset_wide$Class)

  beta_subset_wide$Class <- factor(
    beta_subset_wide$Class
  )

  # Choose rsv if at least in one topic the feature distribution is greater than the cutoff = 0.008
  choose_rsv_logical <- apply(
    beta_subset_wide[, 3:dim(beta_subset_wide)[2]],
    1,
    function(x){
      sum(x >= 0.008) >= 1
    }
  )

  choose_rsv <- beta_subset_wide$rsv[which(
    choose_rsv_logical
  )] %>%
    as.character()


  beta_subset_wide <- beta_subset_wide[choose_rsv, ]

  tree <- phyloseq::phy_tree(
    prune_taxa(
      choose_rsv,
      ps
    )
  )

  circ <- ggtree(
    tree,
    layout = "circular"
  )

  beta_subset_wide <- beta_subset_wide[tree$tip.label, ]

  beta_subset_wide <- beta_subset_wide[ , c("rsv", "Class", paste0("Topic ", seq(1,K)))]

  color <- distinctColorPalette(
    length(
      unique(
        beta_subset_wide$Class
      )
    )
  )

  color_topic <- rainbow(K)

  ## in circular plot, ASVs will be labeled by Class
  df <- data.frame(
    Class = beta_subset_wide$Class
  )

  rownames(df) <- tree$tip.label


  df2 <- beta_subset_wide[, 3:dim(beta_subset_wide)[2]]

  rownames(df2) <- tree$tip.label


  ## append Class to phylogenetic tree
  p <- gheatmap(
    p = circ,
    data = df,
    offset = -.1,
    width = .1,
    colnames_angle = 95,
    colnames_offset_y = .5,
    font.size = 5) +
    scale_fill_manual(
      values = color,
      name = "Class"
    )


  for(i in 1:K){
    if(i == 1){
      p <- p +
        new_scale_fill()# to make once fill scale
    }
    df2_i <- dplyr::select(
      beta_subset_wide, (i+2)
    )
    rownames(df2_i) <- tree$tip.label
    p <- gheatmap(
      p,
      df2_i,
      offset = i*.08,
      width = .1,
      colnames_angle = 90,
      colnames_offset_y = .25,
      font.size = 6,
      high = "dodgerblue",
      low = "gray98",
      legend_title = expression(beta[k]))
  }

  p <- p +
    theme_minimal(
      base_size = 20
    ) +
    theme(
      plot.title = element_text(hjust = 0.5)
    )

  return(p)
}
