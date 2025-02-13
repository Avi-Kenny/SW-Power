# Archived 2025-02-13: Test dataset generating function
if (F) {
  
  # Generate data
  dat <- generate_dataset(
    # data_type = "normal",
    # mu = 5,
    # sigma = 0.5,
    # tau = 0.2,
    # beta_j = c(0,2,4,6),
    # delta_s = c(5,5,5),
    # gamma_j = c(0,0,0,0),
    # n_sequences = 6,
    # n_clust_per_seq = 4,
    # n_ind_per_cluster = 10,
    # n_extra_time_points = 0
  )
  
  # Visualize cluster means
  dat_plot <- dat %>%
    dplyr::group_by(i,j,seq) %>%
    dplyr::summarize(y_bar=mean(y))
  ggplot(
    dat_plot,
    aes(x=j, y=y_bar, color=factor(seq), group=factor(i))
  ) +
    geom_line(alpha=0.5) +
    labs(color="Sequence")
  
}
