########################################.
##### Summarize simulation results #####
########################################.

# sim <- readRDS("SimEngine.out/power_1_20240613.rds")

summ <- SimEngine::summarize(sim,
  list(stat="mean", x="reject", name="power"),
  list(stat="bias_pct", estimate="est", truth="true_tate"),
  list(stat="sd", x="est"),
  list(stat="coverage", estimate="est", se="se", truth="true_tate", name="cov")
)
summ %>% arrange(n_sequences) %>% print()


# > summ %>% arrange(n_sequences)
# level_id data_type sigma  tau n_sequences n_clust_per_seq n_ind_per_cell      re                model n_reps power     sd_est bias_pct_est   cov
# 1         1    normal     1 0.25           6               2             10 cluster                   IT    200 0.435 0.11596684 -0.011056058 0.935
# 2         4    normal     1 0.25           6               2             10 cluster                  ETI    200 0.220 0.16959785 -0.067438266 0.925
# 3         7    normal     1 0.25           6               2             10 cluster ETI, linear cal time    200 0.260 0.16319505 -0.038650225 0.925
# 4        10    normal     1 0.25           6               2             10 cluster     ETI, omit last 1    200 0.265 0.15187016 -0.046303965 0.935
# 5        13    normal     1 0.25           6               2             10 cluster     ETI, omit last 2    200 0.310 0.13923631 -0.056804375 0.935
# 6         2    normal     1 0.25          10               2             10 cluster                   IT    200 0.830 0.06628194  0.011475086 0.955
# 7         5    normal     1 0.25          10               2             10 cluster                  ETI    200 0.435 0.11327592 -0.025863126 0.930
# 8         8    normal     1 0.25          10               2             10 cluster ETI, linear cal time    200 0.465 0.10641800 -0.016315271 0.915
# 9        11    normal     1 0.25          10               2             10 cluster     ETI, omit last 1    200 0.500 0.10216643 -0.017485257 0.935
# 10       14    normal     1 0.25          10               2             10 cluster     ETI, omit last 2    200 0.560 0.09545552 -0.013669123 0.945
# 11        3    normal     1 0.25          14               2             10 cluster                   IT    200 0.975 0.05220920  0.005975401 0.930
# 12        6    normal     1 0.25          14               2             10 cluster                  ETI    200 0.700 0.08389943  0.019414842 0.950
# 13        9    normal     1 0.25          14               2             10 cluster ETI, linear cal time    200 0.750 0.07863489  0.019534617 0.940
# 14       12    normal     1 0.25          14               2             10 cluster     ETI, omit last 1    200 0.770 0.08027104  0.020513546 0.925
# 15       15    normal     1 0.25          14               2             10 cluster     ETI, omit last 2    200 0.785 0.07644749  0.019354199 0.920



############################################.
##### Test dataset generating function #####
############################################.

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



############################.
##### VIZ: Scatterplot #####
############################.

# Figures produced: fig_1, fig_2

if (F) {
  
  # sim <- readRDS("SimEngine.out/estimation_1_20231112.rds")
  
  summ <- sim %>% SimEngine::summarize(
    list(stat="mean", x="est"),
    list(stat="bias", estimate="est", truth="true_tate")
  )
  
  print(summ)
  
}
