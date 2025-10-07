# Archived 2025-02-25: Old versions of "extra time" plots
if (F) {
  
  df_plot_1 <- create_df(
    iccs = c(0, 0.005, 0.01, 0.05, 0.1),
    cacs = 1,
    effect_sizes = 0.15
  )
  plot_1 <- ggplot(df_plot_1, aes(x=time, y=power, color=which)) +
    geom_line() +
    geom_point() +
    scale_color_manual(values=c("#009E73", "#56B4E9")) +
    # scale_color_manual(values=c("#009E73", "#56B4E9", "#CC79A7")) +
    scale_y_continuous(
      labels = scales::percent,
      limits = c(0.4,1),
      breaks = seq(0.4,1,0.1)
    ) +
    facet_grid(cols=dplyr::vars(icc)) +
    labs(
      x = "# Extra time points",
      y = "Power",
      color = NULL
    ) +
    theme(legend.position="bottom")
  
  ggsave(
    filename = paste0("../Figures + Tables/", cfg2$d,
                      " fig_power_extra_time_tate.pdf"),
    plot=plot_1, device="pdf", width=9, height=4
  )
  
  df_plot_1 <- create_df(
    iccs = c(0.01, 0.05, 0.1),
    cacs = c(0.25,0.5,0.75,1),
    effect_sizes = 0.15
  )
  plot_1 <- ggplot(df_plot_1, aes(x=time, y=power, color=which)) +
    geom_line() +
    geom_point() +
    scale_color_manual(values=c("#009E73", "#56B4E9")) +
    scale_y_continuous(
      labels = scales::percent,
      limits = c(0.2,1),
      breaks = seq(0.2,1,0.1)
    ) +
    facet_grid(rows=dplyr::vars(icc), cols=dplyr::vars(cac)) +
    labs(
      x = "# Extra time points",
      y = "Power",
      color = NULL
    ) +
    theme(legend.position="bottom")
  
  ggsave(
    filename = paste0("../Figures + Tables/", cfg2$d,
                      " fig_power_extra_time_new.pdf"),
    plot=plot_1, device="pdf", width=9, height=6
  )
  
}

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
